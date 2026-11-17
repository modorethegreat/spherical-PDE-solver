"""Compare the four-hour inviscid Galewsky case with torch-harmonics 0.8.0."""

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from scipy.integrate import quad
from torch_harmonics.examples import ShallowWaterSolver


def reference(slice_path):
    with np.load(slice_path) as data:
        latitude = data["latitude_rad"]
        longitude = data["longitude_rad"]
        fork = data["divergence_s1"]
    degree = 85
    duration, dt = 14_400, 30
    torch.set_num_threads(1)
    solver = ShallowWaterSolver(
        latitude.size, longitude.size, dt,
        lmax=degree + 1, mmax=degree + 1, grid="legendre-gauss",
    ).double()
    if fork.shape != (latitude.size, longitude.size) or not np.isfinite(fork).all():
        raise ValueError("Expected a finite Gaussian-grid divergence field.")
    np.testing.assert_allclose(solver.lats.numpy()[::-1], latitude, atol=1e-12, rtol=0)
    np.testing.assert_allclose(solver.lons.numpy(), longitude, atol=1e-12, rtol=0)
    solver.hyperdiff.fill_(1)
    gravity, radius, omega = 9.80616, 6_371_220, 7.292e-5
    south, north = np.pi/7, np.pi/2 - np.pi/7
    width = north - south

    def wind(phi):
        if not south < phi < north:
            return 0.0
        position = (2 * phi - south - north) / width
        return 80 * np.exp((4 / width**2) * (1 - 1 / (1 - position**2)))

    def gradient(phi):
        u = wind(phi)
        return -u * (2 * radius * omega * np.sin(phi) + u * np.tan(phi)) / gravity

    lat = solver.lats.numpy()
    lon = solver.lons.numpy()
    zonal_depth = np.array([
        quad(gradient, south, min(phi, north), epsabs=1e-9)[0] if phi > south else 0
        for phi in lat
    ])
    weights = solver.quad_weights.numpy().ravel() / 2
    zonal_depth += 10000 - np.dot(zonal_depth, weights)
    centered_lon = (lon + np.pi) % (2*np.pi) - np.pi
    # Eq. 4 squares the entire latitude offset divided by beta.
    lat_offset = lat[:, None] - np.pi / 4
    bump = 120 * np.cos(lat[:, None]) * np.exp(-(3 * centered_lon)**2 - (15 * lat_offset)**2)
    east = np.broadcast_to(np.array([wind(phi) for phi in lat])[:, None], fork.shape)
    uv = torch.from_numpy(np.stack((east, np.zeros_like(east))))
    height = solver.grid2spec(torch.from_numpy(gravity * (zonal_depth[:, None] + bump)))
    vort, div = solver.vrtdivspec(uv)
    state = torch.stack((height, vort, div))
    state = solver.timestep(state, duration // dt)
    independent = solver.spec2grid(state[2]).numpy()[::-1]
    height = solver.spec2grid(state[0]).numpy()[::-1] / gravity
    weights = weights[::-1]

    def rms(values):
        return float(np.sqrt(np.dot(weights, np.mean(values**2, axis=1))))

    difference = rms(independent - fork) / rms(independent)
    return {
        "degree": degree, "hours": 4, "independent_dt_s": dt,
        "independent_rms_s1": rms(independent), "fork_rms_s1": rms(fork),
        "relative_field_l2": difference, "relative_field_l2_limit": 0.01,
        "independent_divergence_extrema_s1": [float(independent.min()), float(independent.max())],
        "independent_height_extrema_m": [float(height.min()), float(height.max())],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fork_fields", type=Path, help="T85, four hours, no diffusion")
    args = parser.parse_args()
    result = reference(args.fork_fields)
    print(json.dumps(result, indent=2, allow_nan=False))
    if result["relative_field_l2"] > result["relative_field_l2_limit"]:
        raise SystemExit(1)
