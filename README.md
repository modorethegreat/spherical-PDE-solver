# spherical-PDE-solver
Numerical program to solve PDEs on a spherical surface by a spectral methods using spherical harmonics.

## Shallow water equations as an example

The shallaw water equations are a special case of the Navier-Stokes equations, when the depth of water is much smaller than horizontal scales. The simplest linearized shallow water equations read

$$
\frac{∂h}{∂t} + H(\nabla \cdot \mathbf v) = 0
$$

$$
\frac{∂\mathbf v}{∂t} = -g\nabla h - b\mathbf v
$$

Solution fields: $h$ is the water height, $\mathbf v$ is the velocity field. Parameters: $H$ is the mean height, $g$ is the acceleration due to gravity, and $b$ is the rate of viscous dissipation.
