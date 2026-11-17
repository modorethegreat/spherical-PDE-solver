import unittest

import numpy as np

from rotating_benchmarks import run
from rotating_swe import GaussSphere


class RotatingTests(unittest.TestCase):
    def test_vector_operators_and_helmholtz_inversion(self):
        sphere = GaussSphere(10, 2.5)
        potential = np.zeros((11, 11), dtype=complex)
        potential[3, 1] = 0.5 - 0.25j
        potential[5, 0] = 0.15
        streamfunction = np.zeros_like(potential)
        streamfunction[4, 2] = 0.19 + 0.1j
        zeta = sphere.laplacian * streamfunction
        divergence = sphere.laplacian * potential
        east, north = sphere.wind(zeta, divergence)
        np.testing.assert_allclose(sphere.divergence(east, north), divergence, atol=3e-13)
        np.testing.assert_allclose(sphere.vorticity(east, north), zeta, atol=3e-13)
        pe, pn = sphere.gradient(potential)
        se, sn = sphere.gradient(streamfunction)
        np.testing.assert_allclose(east, pe-sn, atol=3e-13)
        np.testing.assert_allclose(north, pn+se, atol=3e-13)

    def test_padded_grid_does_not_alias_highest_zonal_product(self):
        sphere = GaussSphere(12, 1.0)
        coefficients = np.zeros((13, 13), dtype=complex)
        coefficients[-1, -1] = 0.5
        field = sphere.synthesis(coefficients)
        projected = sphere.analysis(field**2)
        # Its two azimuthal frequencies are zero and 24. Order 24 is discarded.
        np.testing.assert_allclose(projected[:, 1:], 0, atol=3e-13)
        self.assertGreater(projected[0, 0].real, 0)

    def test_area_weighted_rms_matches_parseval(self):
        sphere = GaussSphere(8, 1.0)
        coefficients = np.zeros((9, 9), dtype=complex)
        coefficients[0, 0] = 2 * np.sqrt(4 * np.pi)
        coefficients[8, 8] = 1 - 2j
        expected = np.sqrt(4 + 10 / (4 * np.pi))
        self.assertAlmostEqual(sphere.spectral_rms(coefficients), expected, places=13)
        self.assertAlmostEqual(sphere.rms(sphere.synthesis(coefficients)), expected, places=13)

    def test_steady_tilted_williamson_case(self):
        metrics, _ = run("williamson2", degree=10, hours=120, max_dt=600, diffusion=0)
        for name in ("height", "vorticity", "wind"):
            self.assertLess(metrics["relative_l2_" + name], 5e-10)
        self.assertLess(abs(metrics["relative_mean_depth_change"]), 1e-13)

    def test_unsteady_jet_preserves_mass_and_generates_divergence(self):
        metrics, fields = run("galewsky", degree=20, hours=4, max_dt=300, diffusion=0)
        self.assertLess(abs(metrics["relative_mean_depth_change"]), 1e-13)
        self.assertGreater(metrics["final"]["divergence_rms_s1"], 1e-8)
        self.assertGreater(np.min(fields["depth_m"]), 0)

    def test_four_hour_galewsky_agrees_with_independent_rms(self):
        metrics, _ = run("galewsky", degree=85, hours=4, max_dt=60, diffusion=0)
        final = metrics["final"]
        # Independent T85 AB3 reference from validation/check_galewsky_torch.py.
        self.assertAlmostEqual(final["divergence_rms_s1"] / 2.07557257e-7, 1, delta=0.01)
        self.assertAlmostEqual(final["divergence_spectral_rms_s1"] / final["divergence_rms_s1"], 1, places=12)
        self.assertAlmostEqual(final["height_min_m"], 9052, delta=2)
        self.assertAlmostEqual(final["height_max_m"], 10182, delta=2)
