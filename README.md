# spherical-PDE-solver
Numerical program to solve PDEs on a spherical surface by a spectral methods using spherical harmonics.



https://github.com/user-attachments/assets/771b6360-3eda-46e1-8ebd-580be7bab89e



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
\nabla^2 Y_l^m(\theta,\phi) = -\frac{l(l+1)}{R^2}Y_l^m
$$

The spherical harmonics expansions yield

$$
h(\theta,\phi,t) = \sum_{l,m} h_l^m(t) Y_l^m(\theta,\phi)
$$

$$
\mathbf v(\theta,\phi,t) = \sum_{l,m} v_l^m(t) \Psi_l^m(\theta,\phi)
$$

The summation $\sum_{l,m}$ is an abbreviation of a double summation of $l$ from $0$ to $l_\max$ and $m$ from $-l$ to $l$.

## PDEs in spherical harmonics space

Now write the PDEs in the spherical harmonics space (the RHS del operators become linear)

$$
\frac{∂h_l^m}{∂t} = \frac H R l(l+1) v_l^m
$$

$$
\frac{∂v_l^m}{∂t} = \frac g R h_l^m - bv_l^m
$$

## Technical notes about spherical harmonics

> Spherical harmonics are defined as
> 
> $$
> Y_l^m(\theta,\phi) = \tilde P_l^m(\cos\theta)e^{im\phi} \tag{scalar}
> $$
> 
> where $\tilde P_l^m(x)$ is the orthonormalized associated Legendre polynomial. The vector spherical harmonics are defined using the gradient of scalar spherical harmonics (on a unit sphere):
> 
> $$
> \Psi_l^m(\theta,\phi) = \nabla Y_l^m(\theta,\phi) \tag{curl-free vector}
> $$
>
> $$
> \Phi_l^m(\theta,\phi) = \mathbf {R}\times \nabla Y_l^m(\theta,\phi) \tag{divergence-free vector}
> $$
>
> The first is the curl-free component, and the second is the divergence-free component. It is easy to show that in shallow water equations, the velocity field is always curl-free, so only $\Psi_l^m$ terms left.
>
> In the spherical harmonics space, del operators are easier to deal with. For a scalar field,
> 
> $$
> \nabla f \to f_l^m \nabla Y_l^m = \frac 1 R f_l^m\Psi_l^m \tag{scalar gradient}
> $$
>
> $$
> \nabla^2 f \to \nabla^2 Y_l^m = -\frac{l(l+1)}{R^2}f_l^mY_l^m \tag{scalar Laplacian}
> $$
>
> For a vector field,
> 
> $$
> \nabla\cdot \mathbf v \to v_l^m(R\nabla\cdot \Psi_l^m) = v_l^m(R\nabla^2Y_l^m) = -\frac{l(l+1)}{R}v_l^m Y_l^m \tag{vector divergence}
> $$
