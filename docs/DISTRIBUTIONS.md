# Distributional derivatives of finite rational step functions

The module [`calc.distributions`](../src/calc/distributions.py) computes a small but exact part of distribution theory: the **weak derivative** of a compactly supported, piecewise-constant function on the real line. It deliberately does not claim to implement arbitrary distributions, Sobolev spaces, or weak PDE solvers.

## Why jumps become Dirac atoms

For any smooth, compactly supported test function \(\varphi\), the distributional derivative of a locally integrable function \(f\) is defined by

\[
  \langle Df,\varphi\rangle = -\int_{\mathbb R} f(x)\varphi'(x)\,dx.
\]

When \(f=h\mathbf 1_{[a,b)}\), the fundamental theorem of calculus gives

\[
  \langle Df,\varphi\rangle
     = -h(\varphi(b)-\varphi(a))
     = \langle h\delta_a-h\delta_b,\varphi\rangle.
\]

For a sum of disjoint steps, adjacent boundary atoms are added. For example, heights \(3\) on \([0,2)\) and \(-1\) on \([2,3)\) create a jump of \(-4\) at \(2\), not two unrelated discontinuities.

The returned object is a finite *signed atomic measure*. This is a special kind of distribution, not an ordinary differentiable function.

## Executable example

```python
from fractions import Fraction
from calc.distributions import PiecewiseConstant

f = PiecewiseConstant([
    (0, 2, Fraction(3, 2)),
    (2, 3, -1),
])

assert f.integral() == 2
assert f.weak_derivative().atoms == (
    (Fraction(0), Fraction(3, 2)),
    (Fraction(2), Fraction(-5, 2)),
    (Fraction(3), Fraction(1)),
)

# p(x) = 1 + 2x + x², ascending coefficient order.
p = (1, 2, 1)
assert f.integration_by_parts_residual(p) == 0
assert f.weak_derivative().total_mass() == 0

# Translation and scalar multiplication commute with D.
assert f.shift(Fraction(1, 3)).weak_derivative().atoms == (
    f.weak_derivative().shift(Fraction(1, 3)).atoms
)
```

## Exact API and contracts

| Operation | Meaning |
| --- | --- |
| `PiecewiseConstant(intervals)` | Rational half-open intervals `(left, right, height)` in increasing order; zero outside intervals. |
| `f.value_at(x)` | Selected half-open pointwise representative; distributional derivatives do not depend on endpoint point values. |
| `f.integral()` | Exact rational total integral. |
| `f.pair_polynomial(coefficients)` | Exact \(\int f(x)p(x)\,dx\) for coefficients in ascending powers. |
| `f.weak_derivative()` | Exact finite signed Dirac measure of endpoint jumps. |
| `f.integration_by_parts_residual(coefficients)` | Exact \(\langle Df,p\rangle+\langle f,p'\rangle\), always zero by the definition. |
| `DiracMeasure(atoms)` | Sorted, merged signed masses `(point, weight)`; zero atoms removed. |
| `mu.pair_polynomial(coefficients)` | Exact \(\sum_i w_i p(x_i)\). |
| `mu.total_mass()` / `mu.total_variation()` | Sum of signed weights / sum of absolute weights. |
| `.shift(d)` / `.scale(c)` | Exact translation and coefficient scaling. |

**Polynomial test-function nuance.** A nonzero global polynomial is not compactly supported, so is *not* literally a member of the usual test space \(C_c^\infty\). Because all supported functions and atomic derivatives here have compact support, the polynomial pairings equal pairings with a smooth compactly supported test function that agrees with that polynomial on a neighborhood of the support. This cutoff extension is an interpretation, not an algorithm for constructing a smooth cutoff.

**Scope and safeguards.** Inputs must be integer or `fractions.Fraction`, never float or bool. Rational numerators and denominators are limited to 4096 bits on input; at most 2048 ordered nonoverlapping input intervals, 4096 input atoms, and 65 polynomial coefficients are accepted. Adjacent equal-height intervals are merged. Zero intervals are dropped after validation. Overlapping, reversed, zero-width, unsorted, malformed or non-exact intervals are rejected. Results can grow in bit length from repeated rational arithmetic, so users should still avoid pathological symbolic inputs. We do not approximate a Dirac delta using floating-point samples.

## Mathematical sources

- MIT OpenCourseWare, *Step and Delta Functions: Integrals and Generalized Derivatives* ([course resource](https://ocw.mit.edu/courses/18-03sc-differential-equations-fall-2011/resources/step-and-delta-functions-integration-and-generalized-derivatives/)). Introduces step functions, impulses and generalized derivatives.
- MIT OpenCourseWare 18.303, *When functions have no value(s)* ([lecture notes PDF](https://ocw.mit.edu/courses/18-303-linear-partial-differential-equations-analysis-and-numerics-fall-2014/df601974d067ae5b900c1ae01c5ba416_MIT18_303F14_delta_notes.pdf)). Derivatives of distributions are defined through integration by parts against compactly supported tests.
- Stanford Mathematics, *18.03 Class 31* ([notes](https://math.stanford.edu/~vakil/18.03final/class31.html)). The difference of two Heaviside steps produces a finite-width indicator.

Sources are mathematical references, not evidence of copied implementation or a formal verification of this module. The identities here are verified by deterministic, exact-rational unit tests.
