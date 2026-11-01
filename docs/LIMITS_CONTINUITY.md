# Limits and continuity: exact witnesses versus numerical evidence

The module [`calc.limits`](../src/calc/limits.py) deliberately has **two different evidence levels**:

1. An exact rational **ε–δ certificate** for polynomial limits and polynomial continuity at rational centers.
2. A finite one-sided sampling diagnostic for arbitrary finite-real callbacks. Sampling is **not a mathematical proof of a limit**, regardless of how convincing the table looks.

This split matters: a graph or a million nearby samples cannot rule out a spike closer to the point than any sample.

## Exact polynomial ε–δ witness

Write a polynomial with rational coefficients in ascending powers,

\[
P(x)=\sum_{k=0}^n a_k x^k.
\]

For a rational center \(c\), positive rational radius \(r\), and positive rational \(\varepsilon\), set

\[
B=|c|+r,\qquad
M=\sum_{k=1}^n k|a_k|B^{k-1}.
\]

For \(|x-c|\le r\), both \(|x|\) and \(|c|\) are at most \(B\). The finite factorization

\[
x^k-c^k=(x-c)\sum_{j=0}^{k-1}x^{k-1-j}c^j
\]

gives \(|P(x)-P(c)|\le M|x-c|\). Thus choose

\[
\delta=\begin{cases}
r,&M=0,\\
\min\{r,\varepsilon/(2M)\},&M>0.
\end{cases}
\]

For every **real** \(x\) with \(0<|x-c|<\delta\), the algebraic bound yields \(|P(x)-P(c)|<\varepsilon\). Although the code works with exact rational representations, this inequality establishes the real-variable limit, since the factorization is valid over the reals as well.

```python
from fractions import Fraction
from calc.limits import polynomial_limit_certificate

# P(x) = x^2, center c = 1, radius 1, epsilon = 1/8.
certificate = polynomial_limit_certificate((0, 0, 1), 1, Fraction(1, 8))
assert certificate.value == 1
assert certificate.lipschitz_bound == 4
assert certificate.delta == Fraction(1, 64)

# A single rational witness can be checked against the certified inequality.
assert certificate.error_bound_at(Fraction(1) + certificate.delta / 2) < certificate.epsilon
```

This is **conservative**, not the mathematically largest possible \(\delta\). The returned object records normalized polynomial coefficients, center, value, ε, δ, local radius and proven bound. It is immutable. The `error_bound_at(x)` helper is only defined for rational points inside the open δ-ball; it is an illustration of the already-proven inequality, not a method for discovering the proof by sampling.

### Important restrictions

- Coefficients, center, ε and radius must be `int` or `fractions.Fraction`, never implicit floats.
- Polynomial degree is at most 32. Rational intermediates have an explicit 8,192-bit numerator/denominator budget; exceeding it raises `ValueError` instead of hanging or silently rounding.
- This certifies a specific polynomial at a specific rational center. It does **not** certify arbitrary functions, rational functions with poles, non-polynomial compositions, discontinuities or limits at infinity.
- The zero polynomial and constant polynomials use \(M=0\) and may take the full supplied radius.

## One-sided finite sampling (nonproof)

```python
from calc.limits import one_sided_samples

# A jump is apparent in the finite observations.
samples = one_sided_samples(lambda x: -1 if x < 0 else 2,
                            0, initial_step=0.5, count=6)
assert samples.left == (-1.0,) * 6
assert samples.right == (2.0,) * 6

# A removable point need not be evaluated at all.
samples = one_sided_samples(lambda x: (x*x-1)/(x-1),
                            1, initial_step=0.25, count=8)
```

The diagnostic returns only the sampled offsets and left/right values. It intentionally has **no `converges`, `limit`, or `continuous` verdict**. No callback invocation occurs at the center. If an offset rounds back to the center or a callback result is nonfinite, the diagnostic rejects the sample rather than suggesting that the machine reached a distinct mathematical point.

This is useful for education, plotting and hypothesis generation—but cannot certify or refute a general limit from finite data.

## Verification

```bash
PYTHONPATH=src python -m unittest tests.test_limits -v
```

The tests independently calculate linear, quadratic and rational cubic bounds, evaluate distinct nearby rational witnesses, exercise constant and zero polynomials, diagnose jump/removable examples, and reject invalid or nonrepresentable inputs.

## Sources and attribution

- OpenStax, *Calculus Volume 1*, [§2.5 The Precise Definition of a Limit](https://openstax.org/books/calculus-volume-1/pages/2-5-the-precise-definition-of-a-limit), epsilon–delta definition and one-sided limits.
- OpenStax, *Calculus Volume 1*, [§2.4 Continuity](https://openstax.org/books/calculus-volume-1/pages/2-4-continuity), polynomial continuity and relationship between the limit and function value.

These are reference links only. The coefficient-bound implementation and explanatory argument were written for this repository, without copying textbook source code.
