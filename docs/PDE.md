# Partial differential equations: two bounded 1D finite-difference examples

This section implements **two specific explicit numerical schemes**, not a
general-purpose PDE solver. Both work on a uniform one-dimensional mesh,
with fixed Dirichlet values at the endpoints. The independent analytic and
discrete witnesses are in `tests/test_pde.py`.

## Heat/diffusion equation

For a constant, nonnegative diffusivity \(\kappa\), consider

\[
\frac{\partial u}{\partial t}=\kappa\frac{\partial^2 u}{\partial x^2},
\qquad 0<x<L.
\]

At nodes \(x_j=j\Delta x\) and times \(t_n=n\Delta t\), the forward-time,
centered-space (FTCS) update is

\[
u^{n+1}_j=u^n_j+r(u^n_{j-1}-2u^n_j+u^n_{j+1}),\qquad
r=\frac{\kappa\Delta t}{(\Delta x)^2}.
\]

For this linear constant-coefficient scheme with these boundary conditions,
the module **rejects** \(r>\tfrac12\). At \(0\le r\le\tfrac12\), each interior
update is a convex combination of the old node and its neighbours, so an
interior value cannot exceed their old range. The truncation error for smooth
solutions is first-order in time and second-order in space.

```python
import math
from calc.pde import heat_1d_dirichlet

grid = [math.sin(math.pi*j/100) for j in range(101)]
history = heat_1d_dirichlet(grid, diffusivity=1.0,
                            dx=0.01, dt=0.00002, steps=1000)
# Each history[k] is the profile at t = k*dt.
print(history[-1][50])  # about exp(-pi**2 * 0.02), not exactly equal
```

Independent continuum witness: with zero endpoint values and initial
\(\sin(\pi x/L)\), the exact heat solution is
\(u(x,t)=e^{-\kappa(\pi/L)^2t}\sin(\pi x/L)\). The suite compares the numerical
center node against this result. It also checks the separate discrete-mode
amplification factor \(1-4r\sin^2(\pi/[2(N-1)])\).

## Wave equation

For nonnegative constant propagation speed \(c\), consider

\[
\frac{\partial^2u}{\partial t^2}=c^2\frac{\partial^2u}{\partial x^2}.
\]

With Courant number \(\nu=c\Delta t/\Delta x\), the interior update is

\[
u^{n+1}_j=2u^n_j-u^{n-1}_j+
\nu^2(u^n_{j-1}-2u^n_j+u^n_{j+1}).
\]

The first time slice incorporates the initial velocity \(v_j\):

\[
u^1_j=u^0_j+\Delta t\,v_j+
\frac{\nu^2}{2}(u^0_{j-1}-2u^0_j+u^0_{j+1}).
\]

The code enforces the usual one-dimensional Courant stability bound
\(\nu\le1\). For smooth fields its truncation error is second order in time
and space. Fixed endpoint displacement requires zero endpoint velocity.

```python
import math
from calc.pde import wave_1d_dirichlet

initial = [math.sin(math.pi*j/100) for j in range(101)]
history = wave_1d_dirichlet(initial, [0.0]*101,
                            wave_speed=2.0, dx=0.01, dt=0.003, steps=100)
print(history[-1][50])  # about cos(2*pi*0.3), not exact
```

Independent witness: for \(u(x,0)=\sin(\pi x/L)\) and zero initial velocity,
the exact solution is
\(u(x,t)=\sin(\pi x/L)\cos(c\pi t/L)\). We also verify a
zero-speed case, steady affine profiles, and the Courant equality case.

## Contracts and limitations

- Both APIs return immutable tuples of **every** time slice, including the
  initial profile. They reject a request exceeding five million grid values;
  this deliberately caps accidental memory exhaustion.
- The input displacement includes both boundary nodes and must agree with the
  specified fixed boundary values; a tiny floating-point discrepancy at a
  mathematically zero endpoint is tolerated and normalized.
- Inputs must be finite real numbers. `steps` is a nonnegative integer; grid
  spacing and time step are strictly positive. Nonfinite intermediate states
  raise `ArithmeticError`.
- These are fixed-coefficient, 1D, linear models on equal spacing. There is
  no adaptive mesh, Neumann/Robin boundary, multidimensional domain,
  forcing, nonlinear term, inverse problem, spectral/FEM discretization,
  uncertainty quantification, or stiffness handling.
- A CFL bound is a **stability condition**, not an assurance of accuracy.
  Coarse stable meshes still yield poor approximations; smoothness and
  convergence studies matter. No rigorous global error certification is
  provided.

## Research references

The equations and explanations above were independently authored; no
third-party source text or code was copied.

1. [MIT OpenCourseWare 18.303: *Linear Partial Differential Equations — Analysis and Numerics*](https://ocw.mit.edu/courses/18-303-linear-partial-differential-equations-analysis-and-numerics-fall-2014/), covering heat, wave, finite differences and stability.
2. [MIT OpenCourseWare 18.086, Lecture 7: *Finite Differences for the Heat Equation*](https://ocw.mit.edu/courses/18-086-mathematical-methods-for-engineers-ii-spring-2006/resources/lecture-7-finite-differences-for-the-heat-equation/).
3. [MIT OpenCourseWare 18.336: *Numerical Methods for Partial Differential Equations*, lecture notes](https://ocw.mit.edu/courses/18-336-numerical-methods-for-partial-differential-equations-spring-2009/pages/lecture-notes/), including heat, wave, stability and consistency.
