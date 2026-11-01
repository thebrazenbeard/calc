# Bounded nonsmooth and convex calculus

`calc.nonsmooth` provides **exact scalar convex piecewise-linear** operations,
not a general-purpose nonsmooth calculus engine. It uses only Python's standard
library and returns rational numbers (`fractions.Fraction`).

## Scope and formulas

For a nonempty list of rational affine functions,

\[
f(x)=\max_{1\leq i\leq n}\{m_i x+b_i\},\qquad
I(x)=\{i:m_i x+b_i=f(x)\}.
\]

The maximum is always finite and convex. At a query point, the **convex
subdifferential** (all global supporting-line slopes) is exactly

\[
\partial f(x)=\operatorname{conv}\{m_i:i\in I(x)\}
=\left[\min_{i\in I(x)}m_i,\max_{i\in I(x)}m_i\right].
\]

The directional derivative is
`f'(x;d) = max(m_i * d for i in I(x))`. At a differentiable point the
subdifferential interval collapses to a single slope; at a kink it need not.
For instance, `abs(x)` has interval `[-1, 1]` at `x=0`, whereas
`max(0, x)` (the hinge) has `[0, 1]` there.

The **Clarke generalized gradient** is a separate, more general notion that
exists for locally Lipschitz nonconvex functions as well. For **these finite
convex functions alone**, the Clarke generalized gradient agrees with the
convex subdifferential. The method `clarke_subdifferential` exposes this
special-case equality; it does not implement nonconvex Clarke analysis.

## Example

```python
from fractions import Fraction
from calc.nonsmooth import (
    ConvexPiecewiseLinear,
    absolute_subgradient,
    soft_threshold,
    prox_hinge,
    moreau_envelope_absolute,
)

f = ConvexPiecewiseLinear(((1, 0), (-2, 1)))
assert f.value(Fraction(1, 3)) == Fraction(1, 3)
assert f.subdifferential(Fraction(1, 3)) == (Fraction(-2), Fraction(1))
assert absolute_subgradient(0) == (Fraction(-1), Fraction(1))

assert soft_threshold(-3, 2) == -1
assert prox_hinge(1, 3) == 0
assert moreau_envelope_absolute(Fraction(1, 2), 1) == Fraction(1, 8)
```

### Proximal operations

For nonnegative `weight = t`, `soft_threshold(v,t)` minimizes
`t*abs(x) + (x-v)**2/2`, yielding `sign(v)*max(abs(v)-t, 0)`.
`prox_hinge(v,t)` similarly minimizes `t*max(0,x) + (x-v)**2/2`:

| Input | Proximal output |
|---|---|
| `v <= 0` | `v` |
| `0 < v <= t` | `0` |
| `v > t` | `v-t` |

For strictly positive `t`, `moreau_envelope_absolute(v,t)` computes
`min_x abs(x) + (x-v)**2/(2*t)`. It is `v**2/(2*t)` when `abs(v)<=t`, and
`abs(v)-t/2` otherwise. These formulas have exact rational output.

## Boundaries and failure modes

- **Only one real scalar variable.** No multivariate, stochastic,
  set-valued, weak, distributional, fractional, or nonconvex nonsmooth calculus.
- **Convexity by construction.** A finite maximum of affine lines is convex;
  there is no API to pass arbitrary nonconvex segments or function callbacks.
- **Exact rational inputs only.** `int` and `Fraction` are accepted; booleans,
  floats, decimals, infinities, and complex numbers raise `TypeError`. In
  particular, floating-point near-kink heuristics are not applied.
- **Explicit resource bounds.** At most 1,024 affine pieces, with each input
  rational numerator and denominator at most 4,096 bits; oversized inputs
  raise `ValueError`. Query cost is linear in the number of pieces.
- **Parameter domains.** Proximal weights must be nonnegative. The Moreau
  envelope parameter must be strictly positive.

## Mathematical sources

- R. T. Rockafellar, *Characterization of the subdifferentials of convex
  functions*, Pacific Journal of Mathematics **17** (1966), pp. 497–510:
  https://doi.org/10.2140/pjm.1966.17.497
- F. H. Clarke, *A new approach to Lagrange multipliers*, Mathematics of
  Operations Research **1** (1976), pp. 165–174:
  https://doi.org/10.1287/moor.1.2.165
- N. Parikh and S. Boyd, *Proximal Algorithms*, Foundations and Trends in
  Optimization **1** (2014), pp. 127–239:
  https://stanford.edu/~boyd/papers/prox_algs.html

The implementation and independent analytic tests are original and do not
copy code or text from these references.
