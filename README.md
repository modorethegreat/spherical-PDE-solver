# spherical-PDE-solver
A spherical-harmonic solver for linear, nonrotating shallow-water waves on a sphere, forked from [Carl Zhu’s original solver](https://github.com/hongkangcarl/spherical-PDE-solver).

## Interactive exact-mode benchmark

Open the [interactive globe on GitHub Pages](https://modorethegreat.github.io/spherical-PDE-solver/). It computes each run in the browser.

From the repository root, start the local viewer with Python 3.10 or later. The viewer and single-mode command-line benchmark use the standard library.

```bash
python3 serve_viewer.py
```

Open [http://127.0.0.1:4210/](http://127.0.0.1:4210/). Choose a harmonic degree and order, drag either globe (or use its arrow keys), and scrub time. The globes compare the height anomaly on the same scale and with the same camera. Increasing the number of timesteps reduces the timestep at fixed final time. **Export run** saves the settings, initial conditions, synchronized coefficients, and refinement table.

Run the same benchmark headlessly:

```bash
python3 mode_benchmark.py --degree 4 --order 2 --steps 100 --drag 0.5 --output mode.json
```

The initial height coefficient is 1 and the initial velocity coefficient is 0. Both globe fields multiply the corresponding height coefficient by the same orthonormal real cosine harmonic with no Condon–Shortley phase. The final-time error combines height and velocity coefficients using $e=\sqrt{(h-h_*)^2+(a/c)(v-v_*)^2}$, where $a=H\ell(\ell+1)/R$ and $c=g/R$. The refinement table halves $\Delta t$ three times and reports $\log_2(e_{\rm coarse}/e_{\rm fine})$. The scheme is first-order in time. Coarse steps may be outside the asymptotic range. The case uses illustrative model units.

## Tests and original movie

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

The tests cover exact underdamped, critically damped, and overdamped modes, constant mean height, timestep refinement, spectral packing, and grid-to-grid convergence. To run the standard-library benchmark tests alone, use `python3 -m unittest discover -s tests -p test_benchmark.py`.

The original forced example remains available with `python main.py` after installing the dependencies and putting FFmpeg on your PATH. It writes `test.mp4`.

https://github.com/user-attachments/assets/7e2edfef-69dc-49f8-b979-d54618c5b571

## Rotating nonlinear benchmarks

The rotating solver evolves layer depth, relative vorticity, and horizontal divergence. It reconstructs wind from a velocity potential and streamfunction, evaluates nonlinear products on a padded Gaussian grid, and advances the retained harmonics with fourth-order Runge–Kutta. The conserved depth mode is separate from the optional Laplacian diffusion of the other modes.

```bash
python rotating_benchmarks.py williamson2 --degree 10 --hours 120 --dt 600
python rotating_benchmarks.py galewsky --degree 85 --hours 4 --dt 60 --fields galewsky-4h.npz
python rotating_benchmarks.py galewsky --degree 85 --hours 144 --dt 360 --diffusion 100000 --fields galewsky-t85.npz
```

The commands print diagnostics as JSON. The optional NPZ contains latitude and longitude in radians plus the final depth, velocity, divergence, and relative vorticity, all in SI units. Fields are on the padded Gauss–Legendre grid. The solver omits topography, moisture, and external forcing.

- **Williamson case 2.** The axis is tilted by 45° so the benchmark also exercises nonzonal modes. At T10 over five days the relative area-weighted height and wind errors were about `1.5e-15` and `1.2e-14`. This exact solution contains low-degree modes, so this check primarily exercises operator signs, orientation, and steady balance. [Williamson et al.](https://www.sciencedirect.com/science/article/pii/S0021999105800166)
- **Galewsky adjustment, four hours.** Gaussian area quadrature and spectral Parseval both give a T85 divergence RMS of `2.07781e-7 s⁻¹`. An independent T85 torch-harmonics run with AB3 stepping at 30 s gives `2.07557e-7 s⁻¹` and differs from this solver by 0.19% in area-weighted relative divergence L2 on the common grid. This solver gives divergence extrema of `−2.14e-6` and `+3.76e-6 s⁻¹` and height extrema of 9052.25 and 10182.40 m, close to the published extrema. The source of the published divergence RMS of `4.0e-7 s⁻¹` remains unresolved. The table also reports a height RMS of 9778 m, below the 10,000 m unperturbed area mean. Equations (5)–(6) require the RMS to be at least the magnitude of the area mean, and the inviscid test adds a positive bump and conserves mass. [Galewsky et al., Table 1](https://www.columbia.edu/~lmp/paps/galewsky%2Betal-TELLUS-2004.pdf)
- **Galewsky roll-up, six days.** With `ν = 100,000 m²/s` and a 360 s timestep, the relative vorticity-field differences from the T85 run fell from 25.8% at T42 to 8.8% at T63 on the common T85 grid. The T85 vorticity extrema were `−7.36e-5` and `+9.28e-5 s⁻¹`, close to the published `−7.3e-5` and `+9.3e-5 s⁻¹`. These are comparisons with published summary values and an internal refinement run. Independent pointwise agreement remains unverified. [Galewsky et al., Table 1](https://www.columbia.edu/~lmp/paps/galewsky%2Betal-TELLUS-2004.pdf) · [Corrected reference contours](https://www.columbia.edu/~lmp/paps/galewsky%2Betal-TELLUS-2004-CORRIGENDUM.pdf)

To repeat the four-hour comparison, generate `galewsky-4h.npz` with the command above and run this validation script in a separate environment with the benchmark dependencies and `torch-harmonics==0.8.0`. It independently initializes the jet from the published equations and checks that the two divergence fields differ by at most 1% in area-weighted relative L2.

```bash
python validation/check_galewsky_torch.py galewsky-4h.npz
```

## The shallow water equations as an example

For a thin fluid layer linearized about rest and constant depth $H$ on a nonrotating sphere, the equations are

$$
\frac{∂h}{∂t} + H(\nabla \cdot \mathbf v) = 0
$$

$$
\frac{∂\mathbf v}{∂t} = -g\nabla h - b\mathbf v
$$

Here $h$ is the height anomaly relative to mean depth $H$, $\mathbf v$ is the tangential velocity, $g$ is gravity, and $b$ is the linear drag rate. The benchmark evolves the anomaly. Adding the constant $H$ to the height in the original movie does not change the pressure gradient.

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
\frac{∂v_l^m}{∂t} = -\frac g R h_l^m - bv_l^m
$$

Eliminating velocity gives a damped oscillator for each height mode:

$$
\ddot h_l^m + b\dot h_l^m + \frac{gH l(l+1)}{R^2}h_l^m=0.
$$

The degree-zero coefficient is the constant mean in unforced runs. The gradient of the degree-zero harmonic vanishes, so there is no corresponding physical velocity mode.

## Technical notes about spherical harmonics

Spherical harmonics are defined as
 
$$
Y_l^m(\theta,\phi) = \tilde P_l^m(\cos\theta)e^{im\phi} \tag{scalar}
$$
 
where $\tilde P_l^m(x)$ is the orthonormalized associated Legendre polynomial. The vector spherical harmonics are defined using the tangential gradient on the physical sphere of radius $R$:
 
$$
\Psi_l^m(\theta,\phi) = R\nabla Y_l^m(\theta,\phi) \tag{curl-free vector}
$$

$$
\Phi_l^m(\theta,\phi) = \hat{\mathbf r}\times R\nabla Y_l^m(\theta,\phi) \tag{divergence-free vector}
$$

The first basis represents the curl-free component and the second represents the divergence-free component. In this linear, nonrotating model, gradient forcing produces no vorticity. With curl-free initial velocity, only the first basis is needed.

In the spherical harmonics space, del operators are easier to deal with. For a scalar field,

$$
\nabla f \to f_l^m \nabla Y_l^m = \frac 1 R f_l^m\Psi_l^m \tag{scalar gradient}
$$

$$
\nabla^2 f \to \nabla^2 Y_l^m = -\frac{l(l+1)}{R^2}f_l^mY_l^m \tag{scalar Laplacian}
$$

For a vector field,

$$
\nabla\cdot \mathbf v \to v_l^m(\nabla\cdot \Psi_l^m) = v_l^m(R\nabla^2Y_l^m) = -\frac{l(l+1)}{R}v_l^m Y_l^m \tag{vector divergence}
$$
