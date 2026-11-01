# Advanced methods: small executable footholds

The code in [`src/calc/advanced.py`](../src/calc/advanced.py) deliberately
implements **five narrowly scoped operations**. Advanced calculus is not one
algorithm: fractional integration, discrete difference operators, exterior
algebra, and calculus of variations have different domains and mathematical
objects. This module offers inspectable examples, not a universal solver.

Use Python 3.11+; no runtime dependencies beyond the standard library.

## Fractional calculus: Riemann–Liouville power integral

For `alpha > 0`, `beta > -1`, and `x > 0`, the left-sided fractional integral
with lower terminal zero is

\[
I_{0+}^{\alpha}[t^{\beta}](x)
= \frac{1}{\Gamma(\alpha)}\int_0^x (x-t)^{\alpha-1}t^\beta\,dt
= \frac{\Gamma(\beta+1)}{\Gamma(\beta+\alpha+1)}x^{\beta+\alpha}.
\]

`riemann_liouville_power_integral(alpha, beta, x)` evaluates **only that
power-function closed form** with the logarithm of gamma to reduce intermediate
overflow. At `x=0`, only a finite zero limit is accepted (`alpha + beta > 0`).
It does not integrate arbitrary functions, implement Caputo derivatives, or
numerically resolve singular integral kernels.

Example: a half-integral of the constant function at `x=4` equals
`4 / sqrt(pi)`.

## Discrete and q-calculus

`finite_difference(f, x, step=h, order=n)` implements the normalized forward
difference

\[
h^{-n}\Delta_h^n f(x)
=h^{-n}\sum_{k=0}^{n}(-1)^{n-k}{n\choose k}f(x+kh).
\]

It is exact for the chosen difference grid (apart from numerical evaluation);
its interpretation as the corresponding continuous derivative is an
**approximation** as `h` tends to zero for suitably smooth functions. High
orders and very small steps can suffer severe floating-point cancellation.

`q_derivative(f, x, q=q)` uses the Jackson q-derivative

\[
D_q f(x)=\frac{f(qx)-f(x)}{(q-1)x}.
\]

This code restricts `q > 0`, `q != 1`, and `x != 0`. Mathematical limits at
`q=1` or `x=0` exist for appropriate functions, but require separate
analytic/numerical treatment and are **not** silently approximated here.
For `f(x) = x²`, `D_q f(x) = (q+1)x`.

## Exterior algebra: wedge product

An exterior k-form is represented by a dictionary from ascending tuples of
basis indices to coefficients:

```python
from calc.advanced import wedge

dx = {(0,): 1}
dy = {(1,): 1}
assert wedge(dx, dy) == {(0, 1): 1}
assert wedge(dy, dx) == {(0, 1): -1}
assert wedge(dx, dx) == {}
```

The implementation sorts the concatenated basis indices and uses permutation
parity for the sign. It enforces homogeneous form degree and kills repeated
basis vectors, so the anti-symmetric and graded-commutative wedge identities
hold. Scalar forms have the empty basis tuple `()`.

This is **exterior algebra over a fixed formal basis**; it does not include
coordinate charts, manifold geometry, exterior differentiation, pullbacks,
Hodge stars, or Stokes' theorem.

## Variational calculus: Euler–Lagrange residual

For a scalar path `q(t)` and smooth Lagrangian `L(t,q,v)`, a stationary path
for the action (under appropriate boundary conditions) obeys

\[
R(t)=\frac{d}{dt}\left(\frac{\partial L}{\partial v}\right)
-\frac{\partial L}{\partial q}=0,
\quad v=\dot{q}.
\]

`euler_lagrange_residual(L, q, t, step=1e-3)` computes a **local central-
difference estimate** of `R(t)`. It cannot prove global optimality, solve an
Euler–Lagrange boundary-value problem, or account for nonsmooth constraints.
Numerical cancellation becomes more severe for extreme step choices.

```python
import math
from calc.advanced import euler_lagrange_residual

L = lambda t, q, v: (v * v - q * q) / 2
assert abs(euler_lagrange_residual(L, math.cos, 0.3)) < 1e-4
```

## Verification

The standard-library suite [`tests/test_advanced.py`](../tests/test_advanced.py)
checks independent closed forms, algebraic identities, invalid domains,
stationary/nonstationary paths, and normal finite differences.

```sh
PYTHONPATH=src python -m unittest discover -s tests -p test_advanced.py -v
```

The test suite gives evidence only for these operations and cases, not a claim
that any entire field of advanced calculus has been implemented.

## Mathematical references (linked, not copied)

- [NIST DLMF §2.6, Distributional Methods](https://dlmf.nist.gov/2.6) — Riemann–Liouville fractional integral definition.
- [NIST DLMF, Gamma Function chapter](https://dlmf.nist.gov/5) — gamma identities and evaluation context.
- [NIST DLMF §17.2(iv), q-Derivatives](https://dlmf.nist.gov/17.2) — q-calculus terminology and derivative operators.
- [MIT OCW 18.101 Analysis II, lecture 19](https://ocw.mit.edu/courses/18-101-analysis-ii-fall-2005/resources/lecture19/) — wedge product and exterior algebra.
- [MIT OCW 18.086, §7.2 Calculus of Variations](https://ocw.mit.edu/courses/18-086-mathematical-methods-for-engineers-ii-spring-2006/e94c05947ed036cd6ad0150102087062_am72.pdf) — Euler–Lagrange variational theory.

These references are educational/technical sources. Their content has not
been vendored or reproduced in this repository.
