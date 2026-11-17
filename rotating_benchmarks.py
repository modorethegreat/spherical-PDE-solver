"""Run Williamson case 2 or the Galewsky jet with reproducible SI parameters."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from rotating_swe import GaussSphere, ShallowWater


RADIUS = 6.37122e6
GRAVITY = 9.80616
OMEGA = 7.292e-5
G04_DEPTH = 10_000.0
G04_SOUTH = np.pi / 7
G04_NORTH = np.pi / 2 - G04_SOUTH


def williamson_case2(sphere, tilt=np.pi / 4):
    """Exact steady height and wind, with the axis tilted by alpha."""
    speed = 2 * np.pi * sphere.radius / (12 * 86400)
    sin_tilt, cos_tilt = np.sin(tilt), np.cos(tilt)
    sin_lat = sphere.sin_lat[:, None]
    cos_lat = sphere.cos_lat[:, None]
    cos_lon, sin_lon = np.cos(sphere.lon), np.sin(sphere.lon)
    axis_dot = cos_tilt * sin_lat - sin_tilt * cos_lat * cos_lon
    depth = (29_400 - (sphere.radius * OMEGA * speed + speed**2 / 2) * axis_dot**2) / GRAVITY
    east = speed * (cos_tilt * cos_lat + sin_tilt * sin_lat * cos_lon)
    north = np.broadcast_to(-speed * sin_tilt * sin_lon, sphere.shape)
    vorticity = 2 * speed / sphere.radius * axis_dot
    return depth, vorticity, east, north


def jet_speed(lat):
    if lat <= G04_SOUTH or lat >= G04_NORTH:
        return 0.0
    return 80 * np.exp(1 / ((lat - G04_SOUTH) * (lat - G04_NORTH))
                        + 4 / (G04_NORTH - G04_SOUTH)**2)


def galewsky_initial(sphere, amplitude=120.0):
    """Integrate gradient-wind balance, then apply the published height bump."""
    def height_gradient(lat):
        east = jet_speed(lat)
        return -east * (2 * OMEGA * sphere.radius * np.sin(lat) + east * np.tan(lat)) / GRAVITY

    profile = np.array([
        quad(height_gradient, G04_SOUTH, min(lat, G04_NORTH), epsabs=1e-9)[0]
        if lat > G04_SOUTH else 0 for lat in sphere.lat
    ])
    profile += G04_DEPTH - float(np.dot(sphere.weights, profile) / 2)
    lon = (sphere.lon + np.pi) % (2 * np.pi) - np.pi
    depth = profile[:, None] + amplitude * sphere.cos_lat[:, None] * np.exp(
        -(lon / (1 / 3))**2 - ((sphere.lat[:, None] - np.pi / 4) / (1 / 15))**2
    )
    east = np.broadcast_to(np.array([jet_speed(lat) for lat in sphere.lat])[:, None], sphere.shape)
    return depth, east


def initial_state(case, sphere, diffusion):
    if case == "williamson2":
        tilt = np.pi / 4
        model = ShallowWater(sphere, GRAVITY, OMEGA, diffusion, axis_tilt=tilt)
        depth, zeta, east, north = williamson_case2(sphere, tilt)
        vorticity = sphere.analysis(zeta)
        expected = (depth, zeta, east, north)
        return model, np.stack((sphere.analysis(depth), vorticity, np.zeros_like(vorticity))), expected
    model = ShallowWater(sphere, GRAVITY, OMEGA, diffusion)
    depth, east = galewsky_initial(sphere)
    zeta = sphere.vorticity(east, np.zeros_like(east))
    return model, np.stack((sphere.analysis(depth), zeta, np.zeros_like(zeta))), None


def run(case, degree, hours, max_dt, diffusion):
    if case not in ("williamson2", "galewsky"):
        raise ValueError("Unknown benchmark case.")
    if not math.isfinite(hours) or hours <= 0 or not math.isfinite(max_dt) or max_dt <= 0:
        raise ValueError("Hours and timestep must be finite and positive.")
    if case == "williamson2" and diffusion != 0:
        raise ValueError("Williamson's steady reference assumes zero diffusion.")
    if case == "williamson2" and degree < 2:
        raise ValueError("Williamson case 2 contains degree-2 height terms.")
    duration = 3600 * hours
    step_count = duration / max_dt
    if not math.isfinite(step_count):
        raise ValueError("The requested run exceeds the representable step count.")
    steps = max(1, math.ceil(step_count))
    dt = duration / steps
    sphere = GaussSphere(degree, RADIUS)
    model, state, expected = initial_state(case, sphere, diffusion)
    height, _, _, east, north = model.fields(state)
    speed = np.max(np.hypot(east, north))
    # Conservative initial explicit-step check. Later nonphysical depth fails fast.
    max_rate = (np.sqrt(GRAVITY * np.max(height)) + speed) * np.sqrt(degree * (degree+1)) / RADIUS
    max_rate += 2 * OMEGA + diffusion * degree * (degree+1) / RADIUS**2
    if dt * max_rate > 2.5:
        raise ValueError("Timestep is too large for the initial wave/advection estimate.")
    start = model.diagnostics(state)
    for _ in range(steps):
        state = model.advance(state, dt)
    end = model.diagnostics(state)
    depth, vorticity, divergence, east, north = model.fields(state)
    result = {
        "case": case, "degree": degree, "hours": hours, "steps": steps, "dt_s": dt,
        "diffusion_m2_s": diffusion, "radius_m": RADIUS, "gravity_m_s2": GRAVITY,
        "rotation_rate_s1": OMEGA, "grid_shape": sphere.shape, "initial": start, "final": end,
        "relative_mean_depth_change": end["mean_depth_m"] / start["mean_depth_m"] - 1,
    }
    if expected is not None:
        h, z, u, v = expected
        result["relative_l2_height"] = sphere.rms(depth - h) / sphere.rms(h)
        result["relative_l2_vorticity"] = sphere.rms(vorticity - z) / sphere.rms(z)
        result["relative_l2_wind"] = math.sqrt(sphere.mean((east-u)**2 + (north-v)**2)
                                               / sphere.mean(u**2 + v**2))
    return result, {"latitude_rad": sphere.lat, "longitude_rad": sphere.lon,
                    "depth_m": depth, "vorticity_s1": vorticity, "divergence_s1": divergence,
                    "east_m_s": east, "north_m_s": north}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=["williamson2", "galewsky"])
    parser.add_argument("--degree", default=42, type=int)
    parser.add_argument("--hours", type=float)
    parser.add_argument("--dt", type=float, default=300, help="Maximum timestep in seconds.")
    parser.add_argument("--diffusion", type=float, default=0, help="Laplacian coefficient in m²/s.")
    parser.add_argument("--fields", type=Path, help="Optional NPZ output of the final Gaussian-grid fields.")
    args = parser.parse_args()
    hours = args.hours if args.hours is not None else (120 if args.case == "williamson2" else 4)
    try:
        result, fields = run(args.case, args.degree, hours, args.dt, args.diffusion)
    except (ValueError, FloatingPointError) as error:
        parser.error(str(error))
    if args.fields:
        np.savez_compressed(args.fields, **fields)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
