"""Sweep radiative relaxation and drag on an idealized tidally locked planet."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from rotating_swe import GaussSphere, ShallowWater


RADIUS = 1.0e7
REDUCED_GRAVITY = 10.0
MEAN_DEPTH = 1_000.0
ROTATION_PERIOD = 3 * 86_400
OMEGA = 2 * np.pi / ROTATION_PERIOD
WAVE_SPEED = math.sqrt(REDUCED_GRAVITY * MEAN_DEPTH)
WAVE_TIME = RADIUS / WAVE_SPEED
DIFFUSION = 1.0e6
TIMESCALES = (0.5, 2.0, 8.0)
REPRESENTATIVE = (0.5, 8.0)


class TidallyLockedWater(ShallowWater):
    """An active atmospheric layer above a quiescent deep layer."""

    def __init__(self, sphere, tau_rad, tau_drag, amplitude):
        super().__init__(sphere, REDUCED_GRAVITY, OMEGA, DIFFUSION)
        if not all(math.isfinite(tau) and tau > 0 for tau in (tau_rad, tau_drag)):
            raise ValueError("Radiative and drag times must be finite and positive.")
        if not 0 <= amplitude < 1:
            raise ValueError("Forcing amplitude must be in [0, 1).")
        self.tau_rad = tau_rad
        self.tau_drag = tau_drag
        illumination = sphere.cos_lat[:, None] * np.cos(sphere.lon)
        self.equilibrium = MEAN_DEPTH * (1 + amplitude * illumination)

    def _forcing_tendency(self, depth, east, north):
        mass_source = (self.equilibrium - depth) / self.tau_rad
        # Inflow enters with zero relative wind. Outflow carries its local wind.
        damping = 1 / self.tau_drag + np.maximum(mass_source, 0) / depth
        east_accel = -damping * east
        north_accel = -damping * north
        return np.stack((
            self.sphere.analysis(mass_source),
            self.sphere.vorticity(east_accel, north_accel),
            self.sphere.divergence(east_accel, north_accel),
        ))


def diagnostics(sphere, depth, east, north):
    zonal_east = np.mean(east, axis=1)
    equator_wind = np.interp(0, sphere.lat, zonal_east)
    zonal_momentum = sphere.radius * sphere.cos_lat * (OMEGA * sphere.radius * sphere.cos_lat + zonal_east)
    harmonic = np.mean((depth - MEAN_DEPTH) * np.exp(-1j * sphere.lon), axis=1)
    equator_mode = np.interp(0, sphere.lat, harmonic)
    return {
        "mean_depth_m": sphere.mean(depth),
        "minimum_depth_m": float(depth.min()),
        "wind_rms_m_s": math.sqrt(sphere.mean(east**2 + north**2)),
        "equatorial_zonal_wind_m_s": float(equator_wind),
        "equatorial_momentum_excess_m2_s": float(sphere.radius * equator_wind),
        "peak_zonal_momentum_latitude_deg": float(np.rad2deg(sphere.lat[np.argmax(zonal_momentum)])),
        "equatorial_thickness_amplitude_m": float(2 * abs(equator_mode)),
        "equatorial_thickness_offset_deg_east": float(-np.rad2deg(np.angle(equator_mode)))
        if abs(equator_mode) > 1e-8 else None,
    }


def run_case(sphere, tau_rad, tau_drag, duration, max_dt, amplitude):
    """Use units of R/c for all time arguments and return SI fields."""
    if not all(math.isfinite(value) and value > 0 for value in (duration, max_dt)):
        raise ValueError("Duration and timestep must be finite and positive.")
    model = TidallyLockedWater(sphere, tau_rad * WAVE_TIME, tau_drag * WAVE_TIME, amplitude)
    step_count = duration / max_dt
    if not math.isfinite(step_count):
        raise ValueError("The requested run exceeds the representable step count.")
    steps = max(1, math.ceil(step_count))
    dt = duration * WAVE_TIME / steps
    wavenumber = math.sqrt(-sphere.laplacian[-1, 0])
    rate = 2 * WAVE_SPEED * wavenumber + 2 * OMEGA + DIFFUSION * wavenumber**2
    rate += 1 / model.tau_rad + 1 / model.tau_drag
    if dt * rate > 2.5:
        raise ValueError("Timestep exceeds the estimated wave, diffusion, or relaxation limit.")

    modes = sphere.laplacian.shape[0]
    state = np.zeros((3, modes, modes), dtype=complex)
    state[0, 0, 0] = MEAN_DEPTH * math.sqrt(4 * np.pi)
    # Track the final four wave times at five evenly spaced samples, including the endpoint.
    checkpoints = set()
    if duration >= 4:
        checkpoints = {round(steps * (duration - offset) / duration) for offset in range(1, 5)}
    snapshots = []
    for index in range(steps + 1):
        if index in checkpoints:
            depth, _, _, east, north = model.fields(state)
            snapshots.append((depth, east, north))
        if index < steps:
            state = model.advance(state, dt)

    depth, zeta, divergence, east, north = model.fields(state)
    if not np.isfinite(state).all() or depth.min() <= 0:
        raise FloatingPointError("The final state has invalid depth or wind. Reduce the timestep.")
    summary = diagnostics(sphere, depth, east, north)
    summary.update({"tau_rad_wave": tau_rad, "tau_drag_wave": tau_drag,
                    "forcing_amplitude": amplitude, "steps": steps, "dt_wave": dt / WAVE_TIME})
    if snapshots:
        wind_drift = max(math.sqrt(sphere.mean((east-u)**2 + (north-v)**2)) for _, u, v in snapshots)
        depth_drift = max(sphere.rms(depth-h) for h, _, _ in snapshots)
        relative_wind_drift = wind_drift / max(summary["wind_rms_m_s"], 1e-8 * WAVE_SPEED)
        relative_depth_drift = depth_drift / max(sphere.rms(depth-MEAN_DEPTH), 1e-8 * MEAN_DEPTH)
        summary["late_window"] = {
            "duration_wave": 4,
            "relative_wind_drift": relative_wind_drift,
            "relative_thickness_anomaly_drift": relative_depth_drift,
            "below_two_percent": max(relative_wind_drift, relative_depth_drift) <= 0.02,
        }
    else:
        summary["late_window"] = None
    return summary, {"latitude_rad": sphere.lat, "longitude_rad": sphere.lon,
                     "depth_m": depth, "east_m_s": east, "north_m_s": north,
                     "divergence_s1": divergence, "vorticity_s1": zeta}


def sweep(degree, duration, max_dt, amplitude, refine=False, web_data=None):
    sphere = GaussSphere(degree, RADIUS)
    if refine and degree + 16 > 128:
        raise ValueError("Use T112 or lower when refining by 16 degrees.")
    cases = []
    grids = []
    for tau_rad in TIMESCALES:
        for tau_drag in TIMESCALES:
            summary, fields = run_case(sphere, tau_rad, tau_drag, duration, max_dt, amplitude)
            cases.append(summary)
            if web_data is not None:
                grids.append(viewer_grid(fields))
            if (tau_rad, tau_drag) == REPRESENTATIVE:
                representative_fields = fields
    control, control_fields = run_case(sphere, *REPRESENTATIVE, duration, max_dt, 0.0)
    report = {
        "planet": {"radius_m": RADIUS, "reduced_gravity_m_s2": REDUCED_GRAVITY, "mean_depth_m": MEAN_DEPTH,
                   "rotation_period_s": ROTATION_PERIOD, "wave_time_s": WAVE_TIME, "diffusion_m2_s": DIFFUSION},
        "degree": degree, "duration_wave": duration, "max_dt_wave": max_dt,
        "representative_timescales_wave": REPRESENTATIVE,
        "timescales_wave": TIMESCALES, "cases": cases, "zero_heating_control": control,
    }
    if refine:
        halved, _ = run_case(sphere, *REPRESENTATIVE, duration, max_dt / 2, amplitude)
        comparisons = [{"degree": degree, **halved}]
        for finer_degree in (degree + 8, degree + 16):
            finer = GaussSphere(finer_degree, RADIUS)
            refined, _ = run_case(finer, *REPRESENTATIVE, duration, max_dt / 2, amplitude)
            comparisons.append({"degree": finer_degree, **refined})
        report["representative_refinement"] = comparisons
    if web_data is not None:
        grids.append(viewer_grid(control_fields))
        save_web_data(report, sphere, grids, web_data)
    return report, representative_fields


def viewer_grid(fields):
    return {
        "thickness_m": (fields["depth_m"] - MEAN_DEPTH).round(3).ravel().tolist(),
        "zonal_wind_m_s": fields["east_m_s"].round(3).ravel().tolist(),
    }


def save_web_data(report, sphere, grids, path):
    """Export rounded latitude-major fields with a common scale for each layer."""
    scales = {
        name: max(1, math.ceil(max(abs(value) for grid in grids for value in grid[name])))
        for name in grids[0]
    }
    data = {
        "schema": 1,
        **report,
        "grid": {"latitude_deg": np.rad2deg(sphere.lat).tolist(), "longitude_count": sphere.nlon},
        "scales": scales,
        "fields": grids[:-1],
        "control_fields": grids[-1],
    }
    path.write_text(json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n", encoding="utf-8")


def save_figure(report, fields, path):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.ticker import MultipleLocator

    figure = Figure(figsize=(10, 8), layout="constrained", facecolor="#f7f9fc")
    FigureCanvasAgg(figure)
    panels = figure.add_gridspec(2, 2, height_ratios=(1, 1.15))
    cases = report["cases"]
    timescales = report["timescales_wave"]
    shape = (len(timescales), len(timescales))
    for column, (name, title, unit) in enumerate((
        ("equatorial_zonal_wind_m_s", "Equatorial zonal-mean wind", "m/s"),
        ("equatorial_thickness_offset_deg_east", "Equatorial thickness offset", "° east"),
    )):
        axis = figure.add_subplot(panels[0, column])
        values = np.array([np.nan if case[name] is None else case[name] for case in cases]).reshape(shape)
        limit = max((abs(value) for value in values.flat if np.isfinite(value)), default=1.0) or 1.0
        image = axis.imshow(values, cmap="RdBu_r", vmin=-limit, vmax=limit, interpolation="nearest")
        for (row, col), value in np.ndenumerate(values):
            color = image.cmap(image.norm(value))
            ink = "white" if np.dot(color[:3], (0.213, 0.715, 0.072)) < 0.45 else "#182838"
            label = f"{value:+.3f}" if column == 0 else f"{value:+.1f}"
            if not np.isfinite(value):
                label = "—"
            late_window = cases[row * shape[1] + col]["late_window"]
            if not late_window or not late_window["below_two_percent"]:
                label += "*"
            axis.text(col, row, label, ha="center", va="center", color=ink, fontsize=11)
        axis.set(xticks=range(shape[1]), yticks=range(shape[0]),
                 xticklabels=[f"{t:g}" for t in timescales], yticklabels=[f"{t:g}" for t in timescales],
                 xlabel="Drag time (R/c)", ylabel="Radiative time (R/c)", title=title)
        figure.colorbar(image, ax=axis, label=unit, shrink=0.83)

    axis = figure.add_subplot(panels[1, :])
    longitude = (np.rad2deg(fields["longitude_rad"]) + 180) % 360 - 180
    order = np.argsort(longitude)
    lon = longitude[order]
    lat = np.rad2deg(fields["latitude_rad"])
    anomaly = fields["depth_m"][:, order] - MEAN_DEPTH
    limit = max(abs(anomaly).max(), 1.0)
    levels = np.linspace(-limit, limit, 25)
    periodic_lon = np.r_[lon[-1] - 360, lon, lon[0] + 360]
    periodic_anomaly = np.concatenate((anomaly[:, -1:], anomaly, anomaly[:, :1]), axis=1)
    image = axis.contourf(periodic_lon, lat, periodic_anomaly, cmap="RdBu_r", levels=levels)
    figure.colorbar(image, ax=axis, label="Thickness anomaly (m)", shrink=0.92)
    row_stride, col_stride = max(1, len(lat) // 12), max(1, len(lon) // 24)
    rows, cols = slice(None, None, row_stride), order[::col_stride]
    arrows = axis.quiver(lon[::col_stride], lat[rows], fields["east_m_s"][rows][:, cols],
                         fields["north_m_s"][rows][:, cols], color="#162739", alpha=0.8,
                         angles="uv", scale=620, width=0.0024)
    axis.quiverkey(arrows, 0.98, 1.05, 20, "20 m/s", labelpos="W")
    axis.axvline(-90, linestyle="--", linewidth=0.8, color="#48586e")
    axis.axvline(90, linestyle="--", linewidth=0.8, color="#48586e")
    axis.scatter(0, 0, marker="*", s=90, facecolor="#ffdf66", edgecolor="#182838",
                 zorder=5, label="Substellar point")
    axis.legend(loc="lower right", fontsize=8)
    radiative, drag = report["representative_timescales_wave"]
    axis.set(xlabel="Longitude (° east)", ylabel="Latitude (° north)",
             title=f"Flow and thickness at radiative time {radiative:g} R/c and drag time {drag:g} R/c",
             xlim=(-180, 180), ylim=(-90, 90))
    axis.xaxis.set_major_locator(MultipleLocator(60))
    axis.yaxis.set_major_locator(MultipleLocator(30))
    figure.suptitle(
        f"Circulation on a tidally locked sphere\n"
        f"T{report['degree']} after {report['duration_wave']:g} R/c  ·  "
        "* Late-window drift above 2% or insufficient duration",
        fontsize=15,
    )
    figure.savefig(path, dpi=180)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, default=16)
    parser.add_argument("--wave-times", type=float, default=32, help="Duration in wave crossing times R/c.")
    parser.add_argument("--dt-wave", type=float, default=0.025, help="Maximum step in wave crossing times.")
    parser.add_argument("--amplitude", type=float, default=0.2, help="Fractional equilibrium thickness anomaly.")
    parser.add_argument("--refine", action="store_true", help="Also refine the representative run in time and space.")
    parser.add_argument("--fields", type=Path, help="Save representative Gaussian-grid SI fields as NPZ.")
    parser.add_argument("--web-data", type=Path, help="Save all sweep fields and diagnostics as viewer JSON.")
    parser.add_argument("--figure", type=Path, help="Save the sweep and representative circulation as a figure.")
    args = parser.parse_args()
    try:
        report, fields = sweep(args.degree, args.wave_times, args.dt_wave, args.amplitude, args.refine, args.web_data)
    except (ValueError, FloatingPointError) as error:
        parser.error(str(error))
    if args.fields:
        np.savez_compressed(args.fields, **fields)
    if args.figure:
        save_figure(report, fields, args.figure)
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
