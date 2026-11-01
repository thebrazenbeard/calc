# Exact asymptotic comparisons for rational functions

The module [`calc.asymptotics`](../src/calc/asymptotics.py) makes **symbolic, exact** conclusions for a deliberately restricted class: ratios of bounded-degree polynomials with rational coefficients. It neither estimates asymptotic relations from samples nor decides them for arbitrary real functions.

## Derivation

For nonzero polynomials (P,Q\in\mathbb{Q}[x]), write

\[
P(x)=a_m x^m+\text{lower powers},\qquad Q(x)=b_n x^n+\text{lower powers},
\]

with (a_m,b_n\neq 0). Then, as (x\to+\infty),

\[
\frac{P(x)}{Q(x)}\sim\frac{a_m}{b_n}x^{m-n}.
\]

If instead the lowest nonzero terms are (a_i x^i) and (b_j x^j), then near (x\to 0^+\), the leading term is ((a_i/b_j)x^{i-j}\). The denominator need only be nonzero in a sufficiently small punctured neighborhood or sufficiently far out; finite poles elsewhere are permitted.

This elementary leading-coefficient fact gives exact comparisons for two nonzero rational functions (f,g): inspect the leading power and rational coefficient of (f/g\). At (+\infty\), a negative power gives (f=o(g)\), zero power gives (f=\Theta(g)\), and positive power gives neither. At (0^+\), the signs reverse. Asymptotic equivalence (f\sim g\) requires leading power zero **and** coefficient exactly one.

Definitions follow the [NIST Digital Library of Mathematical Functions, §2.1(i)](https://dlmf.nist.gov/2.1.i): (f\sim g\) means the quotient tends to one; (f=o(g)\) means it tends to zero; and (f=O(g)\) means its absolute value stays bounded. This module uses `theta` only for a finite **nonzero** quotient limit. Its conclusions result from finite exact polynomial algebra, not extrapolation from plotted data.

## Executable examples

```python
from fractions import Fraction
from calc.asymptotics import RationalAsymptotic, compare_growth

x = RationalAsymptotic((0, 1))
x_squared = RationalAsymptotic((0, 0, 1))

large_x = compare_growth(x, x_squared)
assert large_x.big_o and large_x.little_o and not large_x.theta
assert large_x.ratio_limit == Fraction(0)

near_zero = compare_growth(x, x_squared, at="zero_plus")
assert not near_zero.big_o and near_zero.ratio_limit is None

f = RationalAsymptotic((1, -1), (1, 1))
assert f.leading_term() == (0, Fraction(-1))
assert compare_growth(f, RationalAsymptotic((-1,))).asymptotically_equivalent
```

Polynomial coefficients follow the repository's convention: **ascending** powers. `RationalAsymptotic((3, 2, 1), (1,))` represents (3+2x+x^2\). Its `leading_term()` returns `(power, coefficient)` for the asymptotic monomial.

## Scope and safety contracts

- Only `int`, `Fraction`, and existing `calc.algebra.Polynomial` inputs are accepted; floats and booleans are rejected. Inputs are exact rational values, never rounded floating approximations.
- Degree is bounded at 64; coefficient numerator/denominator bit lengths are limited to 4,096. Iterators are read only through the first 66 coefficients, rejecting overlong or infinite iterables without hanging.
- The denominator polynomial cannot be identically zero. A zero numerator is valid but has no nonzero leading monomial.
- The **reference** of `compare_growth(f, g)` must be nonzero, since the comparison uses the quotient (f/g\). A zero `f` is both `O(g)` and `o(g)` for any nonzero rational `g`. This avoids convention-dependent statements about `o(0)`.
- The `ratio_limit` field is an exact `Fraction` for finite limits. `None` denotes a divergent-power case—not an indeterminate or complex limit.
- Only (x\to+\infty\) and (x\to0^+\) are implemented. Negative infinity, oscillatory behavior, logarithms, exponentials, multivariable limits, and general analytic functions remain outside scope.
- The comparison object is a machine-checkable application of a known mathematical theorem to an exact restricted class, **not** a formal proof artifact for arbitrary analysis.

## Verification

```bash
PYTHONPATH=src python -m unittest tests.test_asymptotics -v
```

Tests include analytically derived finite examples, zero/reference conventions, float and resource-bound failures, and 200 reproducible seeded rational-function pairs checked at **both** asymptotic approaches against independently enumerated nonzero leading exponents and coefficients.

Background: [NIST DLMF Chapter 2, Asymptotic Approximations](https://dlmf.nist.gov/2) and [§2.1, Definitions and Elementary Properties](https://dlmf.nist.gov/2.1).
