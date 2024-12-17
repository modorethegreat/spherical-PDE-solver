# spherical-PDE-solver
Numerical program to solve PDEs on a spherical surface by a spectral methods using spherical harmonics.

## The shallow water equations as an example

The shallaw water equations are a special case of the Navier-Stokes equations, when the depth of water is much smaller than horizontal scales. The simplest linearized shallow water equations read

$$
\frac{∂h}{∂t} + H(\nabla \cdot \mathbf v) = 0
$$

$$
\frac{∂\mathbf v}{∂t} = -g\nabla h - b\mathbf v
$$

Solution fields: $h$ is the water height, $\mathbf v$ is the velocity field. Parameters: $H$ is the mean height, $g$ is the acceleration due to gravity, and $b$ is the rate of viscous dissipation.

## Spherical harmonics expansion

To solve these PDEs on a spherical surface, this program used a spectral method. In this spherical surface, we transform the solution fields into **spherical harmonics** $Y_l^m(\theta,\phi)$, which are analogs to Fourier modes. Spherical harmonics are the eigenfunctions of the Laplace operator on a spherical surface:

$$
\nabla^2 Y_l^m(\theta,\phi) = -l(l+1)Y_l^m
$$

The spherical harmonics expansions yield

$$
h(\theta,\phi,t) = \sum_{l,m} h_l^m(t) Y_l^m(\theta,\phi)
$$

$$
\mathbf v(\theta,\phi,t) = \sum_{l,m} v_l^m(t) \pmb \Psi_l^m(\theta,\phi)
$$

The summation $\sum_{l,m}$ is an abbreviation of a double summation of $l$ from $0$ to $l_\max$ and $m$ from $-l$ to $l$.

## PDEs in spherical harmonics space

Now write the PDEs in the spherical harmonics space (the RHS del operators become linear)

$$
\frac{∂h_l^m}{∂t} = Hl(l+1) v_l^m
$$

$$
\frac{∂v_l^m}{∂t} = gh_l^m - bv_l^m
$$

