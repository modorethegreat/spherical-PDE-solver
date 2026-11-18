import gzip
import json
from pathlib import Path
import unittest

import numpy as np

from exoplanet_experiment import MEAN_DEPTH, RADIUS, TIMESCALES, WAVE_TIME, TidallyLockedWater, run_case
from rotating_swe import GaussSphere


class ExoplanetTests(unittest.TestCase):
    def test_day_night_source_and_drag(self):
        sphere = GaussSphere(8, RADIUS)
        model = TidallyLockedWater(sphere, WAVE_TIME, 2 * WAVE_TIME, 0.2)
        rest = np.zeros(sphere.shape)
        depth_source, vorticity_source, divergence_source = model._forcing_tendency(
            MEAN_DEPTH + rest, rest, rest,
        )
        heating = sphere.synthesis(depth_source)
        np.testing.assert_allclose(heating, (model.equilibrium - MEAN_DEPTH) / WAVE_TIME, atol=1e-15)
        self.assertAlmostEqual(sphere.mean(heating), 0, places=15)
        np.testing.assert_array_equal(vorticity_source, 0)
        np.testing.assert_array_equal(divergence_source, 0)

        east = np.broadcast_to(10 * sphere.cos_lat[:, None], sphere.shape)
        _, vorticity_source, divergence_source = model._forcing_tendency(model.equilibrium, east, rest)
        np.testing.assert_allclose(vorticity_source, -sphere.vorticity(east, rest) / model.tau_drag, atol=1e-24)
        np.testing.assert_allclose(divergence_source, 0, atol=1e-24)

        for offset in (-20, 20):
            depth = model.equilibrium + offset
            _, vorticity_source, divergence_source = model._forcing_tendency(depth, east, rest)
            damping = 1 / model.tau_drag + max(-offset, 0) / (model.tau_rad * depth)
            np.testing.assert_allclose(vorticity_source, sphere.vorticity(-damping * east, rest), atol=1e-24)
            np.testing.assert_allclose(divergence_source, sphere.divergence(-damping * east, rest), atol=1e-24)

    def test_zero_heating_control_and_forced_response(self):
        sphere = GaussSphere(8, RADIUS)
        control, _ = run_case(sphere, 0.5, 2, 4, 0.05, 0)
        self.assertLess(control["wind_rms_m_s"], 1e-10)
        self.assertIsNone(control["equatorial_thickness_offset_deg_east"])
        forced, fields = run_case(sphere, 0.5, 2, 4, 0.05, 0.2)
        self.assertGreater(forced["wind_rms_m_s"], 1)
        self.assertGreater(forced["equatorial_thickness_amplitude_m"], 10)
        self.assertAlmostEqual(forced["mean_depth_m"], MEAN_DEPTH, places=9)
        self.assertGreater(forced["minimum_depth_m"], 0)
        np.testing.assert_allclose(fields["depth_m"], fields["depth_m"][::-1], atol=1e-9)
        np.testing.assert_allclose(fields["north_m_s"], -fields["north_m_s"][::-1], atol=1e-9)

    def test_unstable_timestep_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Timestep"):
            run_case(GaussSphere(8, RADIUS), 0.5, 2, 4, 2, 0.2)

    def test_published_globe_fields_match_each_runs_diagnostics(self):
        dataset = Path(__file__).resolve().parents[1] / "docs" / "exoplanet-sweep.json"
        data = json.loads(dataset.read_text(encoding="utf-8"))
        cases = data["cases"]
        self.assertEqual(data["schema"], 2)
        self.assertEqual([(row["tau_rad_wave"], row["tau_drag_wave"]) for row in cases],
                         [(rad, drag) for rad in TIMESCALES for drag in TIMESCALES])
        sphere = GaussSphere(data["degree"], RADIUS)
        np.testing.assert_allclose(data["grid"]["latitude_deg"], np.rad2deg(sphere.lat))
        self.assertEqual(data["grid"]["longitude_count"], sphere.nlon)
        for maps, summary in zip(data["fields"], cases, strict=True):
            with self.subTest(radiative=summary["tau_rad_wave"], drag=summary["tau_drag_wave"]):
                anomaly = np.array(maps["thickness_m"]).reshape(sphere.shape)
                east = np.array(maps["zonal_wind_m_s"]).reshape(sphere.shape)
                self.assertAlmostEqual(anomaly.min() + MEAN_DEPTH, summary["minimum_depth_m"], delta=0.0005)
                equatorial_wind = np.interp(0, sphere.lat, east.mean(axis=1))
                self.assertAlmostEqual(equatorial_wind, summary["equatorial_zonal_wind_m_s"], delta=0.0005)
                for layer, values in maps.items():
                    self.assertLessEqual(max(map(abs, values)), data["scales"][layer])
        for values in data["control_fields"].values():
            np.testing.assert_array_equal(values, np.zeros(np.prod(sphere.shape)))

    def test_published_playback_starts_at_rest_and_matches_solver(self):
        docs = Path(__file__).resolve().parents[1] / "docs"
        data = json.loads((docs / "exoplanet-sweep.json").read_text(encoding="utf-8"))
        timeline = data["playback"]
        times = timeline["times_wave"]
        self.assertEqual(len(times), 73)
        self.assertEqual(times[0], 0)
        self.assertEqual(times[-1], data["duration_wave"])
        self.assertTrue(all(np.diff(times) > 0))
        sphere = GaussSphere(data["degree"], RADIUS)
        quantum = np.array(timeline["quantum"])[:, None, None]
        for summary, maps, run in zip(data["cases"], data["fields"], timeline["runs"], strict=True):
            with self.subTest(radiative=summary["tau_rad_wave"], drag=summary["tau_drag_wave"]):
                raw = gzip.decompress((docs / run["file"]).read_bytes())
                frames = np.frombuffer(raw, dtype="<i2").reshape(len(times), 2, *sphere.shape) * quantum
                np.testing.assert_array_equal(frames[0], np.zeros_like(frames[0]))
                self.assertEqual(run["equatorial_wind_m_s"][0], 0)
                self.assertIsNone(run["thickness_offset_deg_east"][0])
                for index, layer in enumerate(timeline["layers"]):
                    self.assertLessEqual(np.max(np.abs(frames[:, index])), data["scales"][layer])
                    limit = quantum[index, 0, 0] / 2
                    np.testing.assert_allclose(frames[-1, index].ravel(), maps[layer], atol=limit + 0.0005, rtol=0)
                frame = times.index(1.0)
                expected, fields = run_case(sphere, summary["tau_rad_wave"], summary["tau_drag_wave"],
                                            times[frame], data["max_dt_wave"], summary["forcing_amplitude"])
                np.testing.assert_allclose(frames[frame, 0], fields["depth_m"] - MEAN_DEPTH, atol=0.005, rtol=0)
                np.testing.assert_allclose(frames[frame, 1], fields["east_m_s"], atol=0.005, rtol=0)
                self.assertAlmostEqual(run["equatorial_wind_m_s"][frame], expected["equatorial_zonal_wind_m_s"])
                self.assertAlmostEqual(run["thickness_offset_deg_east"][frame], expected["equatorial_thickness_offset_deg_east"])
