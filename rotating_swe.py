"""Pseudospectral rotating shallow water on a triangularly truncated sphere."""

import math

import numpy as np
from pyshtools.legendre import PlmON_d1
from scipy.fft import irfft, next_fast_len, rfft
from scipy.special import roots_legendre


class GaussSphere:
    """Complex-normalized harmonics, packed as [degree, nonnegative order]."""

    def __init__(self, degree, radius):
        if not isinstance(degree, int) or not 1 <= degree <= 128:
            raise ValueError("Choose a truncation between T1 and T128.")
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError("Radius must be finite and positive.")
        count = degree + 1
        nlat = (3 * degree + 1) // 2 + 1
        self.nlon = next_fast_len(3 * degree + 1, real=True)
        self.radius = radius
        self.sin_lat, self.weights = roots_legendre(nlat)
        self.lat = np.arcsin(self.sin_lat)
        self.cos_lat = np.sqrt(1 - self.sin_lat**2)
        self.lon = 2 * np.pi * np.arange(self.nlon) / self.nlon
        self.shape = (nlat, self.nlon)
        self.laplacian = -(np.arange(count, dtype=float) * np.arange(1, count + 1))[:, None] / radius**2
        self.inverse_laplacian = np.zeros_like(self.laplacian)
        self.inverse_laplacian[1:] = 1 / self.laplacian[1:]
        self.zonal_factor = 1j * np.arange(count) / radius
        indexes = np.tril_indices(count)
        self.basis = np.zeros((count, count, nlat))
        self.north = np.zeros_like(self.basis)
        for i, node in enumerate(self.sin_lat):
            p, dp = PlmON_d1(degree, node, csphase=1, cnorm=1)
            self.basis[indexes[0], indexes[1], i] = p
            self.north[indexes[0], indexes[1], i] = dp * self.cos_lat[i] / radius
        area_weights = 2 * np.pi * self.weights
        self.project = self.basis * area_weights
        self.project_east = self.project / self.cos_lat
        self.project_north = self.north * area_weights

    def _fft(self, field):
        return rfft(field, axis=-1)[..., :self.basis.shape[1]] / self.nlon

    def _ifft(self, coefficients):
        return irfft(coefficients, n=self.nlon, axis=-1) * self.nlon

    def analysis(self, field):
        return np.einsum("lmi,...im->...lm", self.project, self._fft(field))

    def synthesis(self, coefficients):
        return self._ifft(np.einsum("lmi,...lm->...im", self.basis, coefficients))

    def gradient(self, coefficients):
        amplitudes = np.einsum("lmi,...lm->...im", self.basis, coefficients)
        east = self._ifft(amplitudes * self.zonal_factor / self.cos_lat[:, None])
        north = self._ifft(np.einsum("lmi,...lm->...im", self.north, coefficients))
        return east, north

    def divergence(self, east, north):
        """Project divergence by integrating by parts on the closed sphere."""
        zonal = np.einsum("lmi,...im->...lm", self.project_east, self._fft(east))
        meridional = np.einsum("lmi,...im->...lm", self.project_north, self._fft(north))
        return self.zonal_factor * zonal - meridional

    def vorticity(self, east, north):
        return self.divergence(north, -east)

    def wind(self, vorticity, divergence):
        east, north = self.gradient(np.stack((divergence, vorticity)) * self.inverse_laplacian)
        return east[0] - north[1], north[0] + east[1]

    def mean(self, field):
        """Area-normalized integral for fields on this Gaussian grid."""
        return float(np.dot(self.weights, np.mean(field, axis=-1)) / 2)

    def rms(self, field):
        return math.sqrt(self.mean(field**2))

    def spectral_rms(self, coefficients):
        """Area-normalized Parseval norm for a real field packed with m >= 0."""
        power = np.abs(coefficients)**2
        return math.sqrt((power[:, 0].sum() + 2 * power[:, 1:].sum()) / (4 * np.pi))


class ShallowWater:
    """State order is layer depth, relative vorticity, horizontal divergence."""

    def __init__(self, sphere, gravity, rotation_rate, diffusion=0.0, axis_tilt=0.0):
        if not all(map(math.isfinite, (gravity, rotation_rate, diffusion, axis_tilt))) or gravity <= 0 or diffusion < 0:
            raise ValueError("Use positive gravity and nonnegative diffusion with finite coefficients.")
        self.sphere = sphere
        self.gravity = gravity
        axis_dot = (np.cos(axis_tilt) * sphere.sin_lat[:, None]
                    - np.sin(axis_tilt) * sphere.cos_lat[:, None] * np.cos(sphere.lon))
        self.coriolis = 2 * rotation_rate * axis_dot
        self.decay = diffusion * sphere.laplacian

    def fields(self, state):
        depth, vorticity, divergence = self.sphere.synthesis(state)
        east, north = self.sphere.wind(state[1], state[2])
        return depth, vorticity, divergence, east, north

    def _forcing_tendency(self, depth, east, north):
        """Optional spectral sources in depth, vorticity, divergence order."""
        return 0.0

    def tendency(self, state):
        depth, vorticity, _, east, north = self.fields(state)
        if not np.isfinite(depth).all() or np.min(depth) <= 0:
            raise FloatingPointError("Depth became non-finite or nonpositive. Reduce the timestep.")
        absolute = vorticity + self.coriolis
        # Vector-invariant momentum uses acceleration (absolute*north, -absolute*east).
        east_fluxes = np.stack((depth * east, absolute * east, absolute * north))
        north_fluxes = np.stack((depth * north, absolute * north, -absolute * east))
        continuity, curl_flux, momentum = self.sphere.divergence(east_fluxes, north_fluxes)
        energy_head = self.gravity * state[0] + self.sphere.analysis((east**2 + north**2) / 2)
        transport = np.stack((-continuity, -curl_flux, momentum - self.sphere.laplacian * energy_head))
        return transport + state * self.decay + self._forcing_tendency(depth, east, north)

    def advance(self, state, dt):
        """Classical fourth-order Runge–Kutta with fixed positive timestep."""
        if not math.isfinite(dt) or dt <= 0:
            raise ValueError("Timestep must be finite and positive.")
        a = self.tendency(state)
        b = self.tendency(state + dt / 2 * a)
        c = self.tendency(state + dt / 2 * b)
        d = self.tendency(state + dt * c)
        return state + dt / 6 * (a + 2*b + 2*c + d)

    def diagnostics(self, state):
        depth, zeta, divergence, east, north = self.fields(state)
        mean = self.sphere.mean
        mean_depth = mean(depth)
        if not np.isfinite(depth).all() or np.min(depth) <= 0:
            raise FloatingPointError("Invalid diagnosed layer depth.")
        return {
            "mean_depth_m": mean_depth,
            "height_rms_m": self.sphere.rms(depth),
            "height_max_m": float(np.max(depth)),
            "height_min_m": float(np.min(depth)),
            "energy_per_area_m3_s2": mean(depth * (east**2 + north**2) / 2
                                            + self.gravity * (depth - mean_depth)**2 / 2),
            "potential_enstrophy": mean((zeta + self.coriolis)**2 / (2 * depth)),
            "divergence_rms_s1": self.sphere.rms(divergence),
            "divergence_spectral_rms_s1": self.sphere.spectral_rms(state[2]),
            "divergence_max_s1": float(np.max(divergence)),
            "divergence_min_s1": float(np.min(divergence)),
            "vorticity_rms_s1": self.sphere.rms(zeta),
            "vorticity_max_s1": float(np.max(zeta)),
            "vorticity_min_s1": float(np.min(zeta)),
        }
