import numpy as np
import pyshtools as pysh 
from scipy.special import roots_legendre

def squarize(array_1d):
    matrix = np.zeros((lmax + 1, lmax + 1), dtype=array_1d.dtype)
    l = np.arange(lmax + 1)  # All possible l values
    indices = l * (l + 1) // 2  # Start indices for each l in the 1D array
    triangle_indices = np.concatenate([np.arange(idx, idx + ell + 1) for idx, ell in zip(indices, l)])
    lower_tri_rows = np.concatenate([np.full(ell + 1, ell) for ell in l])  # Row indices (l)
    lower_tri_cols = np.concatenate([np.arange(ell + 1) for ell in l])  # Column indices (m)
    matrix[lower_tri_rows, lower_tri_cols] = array_1d[triangle_indices]

    return matrix

lmax = 100
#naxo
ls = np.arange(lmax + 1)
cos_theta, weights = roots_legendre(180)
theta = np.arccos(cos_theta)
phi = np.linspace(0, 2*np.pi, 361)[:-1]
theta_mesh, phi_mesh = np.meshgrid(theta, phi, indexing='ij')
Plm = np.zeros((lmax+1, lmax+1, len(theta)))
Plm_1d = np.zeros((lmax+1, lmax+1, len(theta)))
for i, cos_theta_val in enumerate(cos_theta):
    Plm_val = pysh.legendre.PlmON(lmax, cos_theta_val)
    Plm[:, :, i] = squarize(Plm_val)

def Y(l, m):
	return np.tile(Plm[l, np.abs(m)][:, np.newaxis], (1, len(phi))) * np.exp(1j * m * phi_mesh)

def SHT(f):
    f_fourier = np.fft.fft(f, axis=1) / len(phi)
    m_factor = np.ones(lmax+1)
    m_factor[1:] = 0.5
    coeff = np.zeros((lmax+1, lmax+1), dtype=complex)
    for l in range(0, lmax+1):
        coeff[l, :] = np.sum(weights[:, None] * f_fourier[:, :lmax+1] * Plm[l].T, axis=0) * 2 * np.pi * m_factor
    return coeff

def ISHT(F):
    f_theta_m = np.zeros((len(theta), lmax+1), dtype=complex)
    for i, theta_val in enumerate(theta):
        f_theta_m[i] += np.sum(F * Plm[:, :, i], axis=0)
    f_recon = np.zeros(theta_mesh.shape, dtype=complex)
    f_recon[:, :lmax+1] = f_theta_m
    f_recon[:, -lmax:][:, ::-1] = np.conj(f_theta_m[:, 1:])
    f_recon = np.fft.ifft(f_recon, axis=1) * len(phi)
    f_recon = f_recon.real
    return f_recon