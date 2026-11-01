# Chebyshev spectral collocation and interpolation

The dependency-free `calc.collocation` module adds an actual (small) **spectral
collocation building block** to the atlas. In contrast, `calc.spectral`
computes symmetric-matrix eigenpairs. These are related uses of the word
"spectral" but solve different problems.

This module constructs Chebyshev--Lobatto nodes, interpolates values with the
second barycentric formula, and differentiates the unique degree-`N`
interpolating polynomial at those nodes. It does **not** solve a general
boundary-value problem, discretize arbitrary PDEs, or provide certified error
bounds.

## Mathematical conventions

The grid is descending, with `N+1` values:

\[
x_j=\cos(\pi j/N),\quad j=0,\ldots,N,\quad
x_0=1,\quad x_N=-1.
\]

Data `values[j]` belongs to the point `x_j`. With barycentric weights

\[
w_j=(-1)^j c_j,\qquad c_0=c_N=\tfrac12,\quad c_j=1
\quad(0<j<N),
\]

the interpolating polynomial is evaluated away from the nodes by

\[
p(x)=\frac{\displaystyle\sum_{j=0}^{N}\frac{w_j f_j}{x-x_j}}
{\displaystyle\sum_{j=0}^{N}\frac{w_j}{x-x_j}}.
\]

At a grid node, its sample value is returned unchanged. The code scales the
two sums to avoid some intermediate floating-point overflow. The resulting
values remain approximate binary64 floats.

The nodal differentiation matrix uses

\[
D_{ij}=\frac{w_j}{w_i(x_i-x_j)}\quad(i\ne j),\qquad
D_{ii}=-\sum_{j\ne i}D_{ij},\qquad p'(x_i)=\sum_jD_{ij}f_j.
\]

The high-level `differentiate_at_nodes` instead accumulates
\(\sum_{j\ne i}D_{ij}(f_j-f_i)\), so data constant to binary64 precision has
an exactly zero derivative and near-constant samples avoid one source of
diagonal cancellation. These formulas are based on polynomial interpolation,
not automatic differentiation or differentiation of the original function
outside the grid.

## Usage

```python
import math
from calc.collocation import (
    chebyshev_lobatto_nodes,
    barycentric_interpolate,
    chebyshev_differentiation_matrix,
    differentiate_at_nodes,
)

nodes = chebyshev_lobatto_nodes(12)
samples = tuple(math.exp(x) for x in nodes)
at_third = barycentric_interpolate(samples, 1/3)
derivatives = differentiate_at_nodes(samples)
D = chebyshev_differentiation_matrix(12)

assert abs(at_third - math.exp(1/3)) < 1e-10
assert max(abs(estimate - math.exp(x))
           for x, estimate in zip(nodes, derivatives)) < 1e-8
assert len(D) == len(nodes) == 13
```

All functions are deterministic, pure, and leave input iterables untouched
apart from consuming up to 130 values from a sample iterator. Inputs are
finite real numbers, with booleans and textual numerals rejected. Grid
degree is an integer from 1 through 128 inclusive; sample arrays must have
2 through 129 entries. Evaluation points are restricted to `[-1, 1]` to
avoid unstable extrapolation. A numerical overflow raises `ValueError`.

## Accuracy, cost, and limitations

- For any polynomial of degree at most `N`, interpolation and differentiation
  reproduce that polynomial mathematically. In floating-point arithmetic
  there is rounding error, especially with larger `N` and near endpoints.
- For sufficiently smooth functions, Chebyshev interpolation often converges
  rapidly as `N` grows. It need not do so for discontinuous or nonsmooth
  functions; Gibbs behavior and the interpolation error must be evaluated.
- Constructing the dense differentiation matrix costs **O(N²)** memory and
  time; evaluating the barycentric formula at one point costs **O(N)**.
  `differentiate_at_nodes` also uses a dense matrix and costs **O(N²)**.
- The coordinate domain is fixed at `[-1, 1]`; a different interval requires
  an affine coordinate change and derivative scale factor, which this API
  does not silently apply.
- The analytic source function is not stored. The derivative produced here is
  **of its interpolating polynomial**, not the exact original derivative.
  This distinction is crucial for sparse data, noisy data, or nonanalytic
  functions.
- There is no PDE boundary-condition enforcement, matrix conditioning
  estimate, adaptive mesh, or spectral error certificate.

## Independent tests

`tests/test_collocation.py` validates the known three-node differentiation
matrix, degree-five polynomial evaluation/derivatives against direct
formulas, monotone grids, endpoint reproduction, derivative of a constant,
fast convergence for `exp`, invalid inputs, and resource bounds.

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_collocation.py -v
```

## Scholarly primary references

1. Jean-Paul Berrut and Lloyd N. Trefethen, **Barycentric Lagrange
   Interpolation**, *SIAM Review* 46 (2004), 501–517,
   [DOI:10.1137/S0036144502417715](https://doi.org/10.1137/S0036144502417715).
   The paper explains barycentric evaluation and stable forms.
2. Lloyd N. Trefethen, **Spectral Methods in MATLAB**, SIAM (2000),
   [Oxford author page](https://people.maths.ox.ac.uk/trefethen/spectral.html).
   Its Chebyshev differentiation examples motivate the nodal matrix, while
   the author explicitly points readers to barycentric interpolation work.
3. Lloyd N. Trefethen, **Chebyshev Differentiation Matrices**, Chapter 6 of
   *Spectral Methods in MATLAB*,
   [SIAM DOI:10.1137/1.9780898719598.ch6](https://doi.org/10.1137/1.9780898719598.ch6).

The implementation and tests here are newly written, not copied from those
sources or third-party repositories. External materials retain their own
copyright and license terms.
