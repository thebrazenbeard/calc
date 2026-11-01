# Finite-volume advection: a conservative periodic example

`calc.finite_volumes` implements one narrowly scoped finite-volume method:
first-order explicit **upwind advection** of scalar cell averages on a uniform
**one-dimensional periodic** grid with **constant velocity**. It is an
executable introduction to fluxes, conservation, and Courant–Friedrichs–Lewy
(CFL) stability, not a general PDE solver.

## Model and derivation

The conservation law is

\[
\partial_t q + \partial_x(aq)=0, \qquad a\in\mathbb{R}\text{ constant}.
\]

For cells \(i=0,\ldots,N-1\) of common width \(\Delta x>0\), let \(Q_i^n\)
be the *average over cell* \(i\) at time step \(n\), not a nodal sample. The
periodic convention makes every cell index modulo \(N\).

The donor-cell interface flux at the right face of cell \(i\) is

\[
F_{i+1/2}^n=\begin{cases}
a Q_i^n, & a\geq 0,\\
a Q_{i+1}^n, & a<0.
\end{cases}
\]

Conservative finite-volume updates satisfy

\[
Q_i^{n+1}=Q_i^n-\frac{\Delta t}{\Delta x}
\big(F_{i+1/2}^n-F_{i-1/2}^n\big).
\]

Let \(\nu=|a|\Delta t/\Delta x\). For the constant-velocity equation this
becomes the *convex combination*

\[
Q_i^{n+1}=(1-\nu)Q_i^n+
\begin{cases}
\nu Q_{i-1}^n, & a\geq 0,\\
\nu Q_{i+1}^n, & a<0.
\end{cases}
\]

when \(0\leq\nu\leq1\). The implementation checks this CFL bound before
advancing, explicitly rejecting \(\nu>1\). It computes the update in the
convex-combination form so that cell averages obey the discrete maximum
principle, apart from ordinary floating-point roundoff. Because face fluxes
cancel in pairs around the periodic grid, total discrete mass
\(M=\Delta x\sum_iQ_i\) is conserved up to rounding.

At \(\nu=1\) the method transports the array by exactly one cell per step;
for \(0<\nu<1\), first-order upwinding introduces numerical diffusion even
for smooth fields.

## API and examples

```python
from calc.finite_volumes import (
    upwind_fluxes, upwind_step, advance_upwind,
    cell_mass, periodic_total_variation,
)

cells = (1.0, 2.0, 3.0, 4.0)

upwind_fluxes(cells, velocity=2)           # (2.0, 4.0, 6.0, 8.0)
upwind_step(cells, dx=2, velocity=2, dt=1) # (4.0, 1.0, 2.0, 3.0)
upwind_step(cells, dx=2, velocity=-2, dt=1)  # (2.0, 3.0, 4.0, 1.0)

start = (0.0, 1.0, 0.0, 0.0)
finish = advance_upwind(start, dx=0.25, velocity=1.0,
                        dt=0.125, steps=10)
assert abs(cell_mass(finish, dx=0.25) - cell_mass(start, dx=0.25)) < 1e-12
assert periodic_total_variation(finish) <= periodic_total_variation(start) + 1e-12
```

| Function | Meaning | Caveats |
| --- | --- | --- |
| `upwind_fluxes(cells, velocity)` | Flux at each cell's right face | Every flux must be representable as a finite float. |
| `upwind_step(cells, *, dx, velocity, dt)` | One CFL-checked explicit step | Requires periodic uniform cells and `dt >= 0`. |
| `advance_upwind(cells, *, dx, velocity, dt, steps)` | Repeated steps | Nonnegative integer `steps`; bounded work. |
| `cell_mass(cells, *, dx)` | Uniform-grid discrete integral | Returns a finite float or raises; exact binary64 fallback handles intermediate overflow/cancellation. |
| `periodic_total_variation(cells)` | Sum of adjacent absolute jumps, including periodic seam | Returns a finite float or raises. |

The grid must contain 2–100,000 finite real cell averages. The grid
spacing is strictly positive; its total period must be finite and
representable. The number of steps is bounded by 100,000 and the total
cell-update work by 2,000,000. Booleans, nonfinite numbers, and unsized
iterables are rejected. Rational inputs (including `Fraction`) are accepted
but converted to finite binary floating-point values before the update;
this is **not an exact-rational time evolution method**. CFL is checked
using exact fractions reconstructed from those finite binary64 inputs, so
intermediate overflow cannot falsely bypass stability validation.

## Verified properties

`tests/test_finite_volumes.py` exercises:

1. Both upwind directions, conservative flux differences, and \(\nu=1\) periodic translations.
2. Zero-velocity and zero-time-step identities.
3. Constant-state preservation and near-conservation of mass over hundreds of steps.
4. A discrete maximum principle and nonincreasing periodic total variation.
5. Mesh refinement of a transported sinusoid: the finer grid has a smaller
   error than the coarse grid, consistent with first-order convergence.
6. Invalid inputs, out-of-CFL requests, resource limits, and recoverable
   finite sums of individually overflowing float products.

These are numerical and algebraic **tests**, not machine-checked proofs.

## Not implemented

Variable-coefficient or nonlinear conservation laws; nonuniform/adaptive or
multidimensional meshes; arbitrary boundary conditions; Riemann solvers;
shock-capturing limiters; higher-order reconstructions; source terms; and
general finite-volume systems are outside this module's scope. For other
schemes, the admissible CFL condition and stability theory may differ.

## Mathematical references

- Randall J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*,
  Cambridge University Press, Ch. 4 on first-order finite-volume and upwind
  methods: <https://www.cambridge.org/core/books/abs/finite-volume-methods-for-hyperbolic-problems/finite-volume-methods/CB7B0A27A6D37AE3B906D4AE7C60A87E>
- [Clawpack project documentation](https://www.clawpack.org/), a widely used
  implementation and teaching resource for finite-volume conservation-law
  algorithms; our example is an independent small implementation, not a copy
  of its production solver.
