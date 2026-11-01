# Discrete calculus of finite rational sequences

This module provides **exact operators on an already sampled, uniformly spaced finite sequence**. It covers forward/backward/centered differences, discrete antiderivatives, signed finite sums, summation by parts, and polynomial power sums. It is not a continuous differentiator, symbolic summation system, general time-scales calculus, or an infinite-series convergence test.

```python
from fractions import Fraction
from calc.discrete import (
    backward_differences,
    central_differences,
    cumulative_sum,
    discrete_integral,
    forward_differences,
    power_sum,
    summation_by_parts,
)

f = (0, 1, 4, 9, 16)  # x**2 at x = 0, 1, 2, 3, 4
assert forward_differences(f) == (1, 3, 5, 7)
assert forward_differences(f, order=2) == (2, 2, 2)
assert backward_differences(f) == (1, 3, 5, 7)
assert central_differences(f) == (2, 4, 6)

increments = (Fraction(1, 2), Fraction(-1, 3), Fraction(2, 3))
F = cumulative_sum(increments, initial=7)
assert F == (7, Fraction(15, 2), Fraction(43, 6), Fraction(47, 6))
assert discrete_integral(increments, 1, 3) == Fraction(1, 3)
assert discrete_integral(increments, 3, 1) == -Fraction(1, 3)
assert power_sum(5, 3) == 100  # 0**3 + 1**3 + ... + 4**3

lhs, boundary, correction = summation_by_parts((1, 2, 4), (0, 3, 5))
assert (lhs, boundary, correction) == (7, 20, 13)
assert lhs == boundary - correction
```

## Difference operators and alignment

Given values \(f_k=f(x_0+kh)\) and nonzero rational \(h\), first differences are

\[
\frac{\Delta_h f_k}{h}=\frac{f_{k+1}-f_k}{h},\qquad
\frac{\nabla_h f_k}{h}=\frac{f_k-f_{k-1}}{h},\qquad
D_0f_k=\frac{f_{k+1}-f_{k-1}}{2h}.
\]

`forward_differences(values, order=m, step=h)` applies the normalized forward stencil \(m\) times and returns values at *starting* indices \(0,\ldots,N-m-1\). `backward_differences` produces the same numerical values but aligns them at *ending* indices \(m,\ldots,N-1\). `central_differences` returns only the first centered stencil at the interior indices \(1,\ldots,N-2\). For example, the same tuple `(1, 3, 5, 7)` represents \(\Delta f\) at \(x=0,1,2,3\) or \(\nabla f\) at \(x=1,2,3,4\).

For `order=0`, samples pass through unchanged. Missing stencils give empty tuples. A negative nonzero step is allowed **only when the samples really follow a descending regular grid**; changing the sign of `step` does not reorder data. An exact difference quotient is not automatically the exact derivative of a continuous interpolating function.

## Finite sums and a discrete fundamental theorem

For a finite input \((a_0,\ldots,a_{N-1})\), `cumulative_sum(values, initial=C)` returns \((F_0,\ldots,F_N)\) with

\[
F_0=C,\qquad F_{k+1}-F_k=a_k,\qquad
F_j=C+\sum_{k=0}^{j-1}a_k.
\]

Hence `F[k+1] - F[k] == values[k]` exactly. `discrete_integral(values, start, stop)` uses half-open **index bounds**, not physical coordinates:

\[
\operatorname{Sum}(a;i,j)=
\begin{cases}
\sum_{k=i}^{j-1} a_k,&i\le j,\\
-\sum_{k=j}^{i-1} a_k,&i>j.
\end{cases}
\]

This makes reversed bounds antisymmetric and empty intervals zero. For an empty sequence `cumulative_sum(()) == (0,)` and `discrete_integral(()) == 0`. The function sums *sample values* directly: no multiplication by \(h\) is implicit, so it is not a numerical quadrature routine.

## Summation by parts

For nonempty equally sized endpoint sequences \(u_0,\ldots,u_n\) and \(v_0,\ldots,v_n\),

\[
\underbrace{\sum_{k=0}^{n-1}u_k(v_{k+1}-v_k)}_{\mathrm{lhs}}
=\underbrace{u_nv_n-u_0v_0}_{\mathrm{boundary}}
-\underbrace{\sum_{k=0}^{n-1}v_{k+1}(u_{k+1}-u_k)}_{\mathrm{correction}}.
\]

The `summation_by_parts(u, v)` return value is `(lhs, boundary, correction)` to make all three terms independently inspectable. This follows by adding and subtracting \(u_kv_{k+1}\) in \(u_{k+1}v_{k+1}-u_kv_k\), then telescoping. The single-endpoint case is the empty sum. Unlike classical integration by parts, the correction uses the **forward-shifted** \(v_{k+1}\); dropping that shift changes the identity.

## Polynomial sums (Faulhaber)

`power_sum(n, p)` exactly calculates \(\sum_{k=0}^{n-1}k^p\) for \(0\le n\le 10^{12}\) and integer \(0\le p\le 10\). This is evaluated as a polynomial, so no loop over a potentially huge upper bound is required:

\[
\sum_{k=0}^{n-1}k^p
=\frac{1}{p+1}\sum_{j=0}^{p}\binom{p+1}{j}B_j n^{p+1-j},
\]

where the Bernoulli-number convention is \(B_0=1\), \(B_1=-\tfrac12\), \(B_2=\tfrac16\). In particular,

\[
S_0(n)=n,\quad S_1(n)=\frac{n(n-1)}2,\quad
S_2(n)=\frac{n(n-1)(2n-1)}6,\quad
S_3(n)=\left(\frac{n(n-1)}2\right)^2.
\]

`p=0` counts every included index once; the conventional combinatorial \(0^0=1\) is used. Independent test oracles compare against direct finite sums for small \(n\) and against exact polynomial formulas/telescoping for large \(n\). This routine does not perform indefinite symbolic summation, nor does it handle general summands.

## Domain, resource limits, and nonuniform grids

- Inputs must be Python `int` or `fractions.Fraction`, excluding `bool` and floating-point values. Results are rational `Fraction`s (power sums return integer `int`). No silent floating-point approximation is introduced.
- Sequences may contain up to **10,000** samples; the numerator and denominator of any supplied coefficient may use up to **1,024 bits** each. Difference orders are capped at **32**. There is no length-32 prerequisite: an insufficient stencil yields `()`.
- Sum bounds are integer indices between zero and the sample count (inclusive); invalid indices are rejected rather than silently clipped. A zero sampling step is invalid.
- All difference operators assume a **single nonzero constant** \(h\). On a *nonuniform* grid, substituting adjacent distances into the denominator of an already-calculated uniform finite-difference stencil is generally incorrect for higher orders and centered stencils. Use a nonuniform-grid-specific method or the separate [`time_scales`](../src/calc/time_scales.py) module for the API and mathematical setting it actually implements; do not assume that module is a complete arbitrary-time-scale solver.
- Discrete finite sums are not the time-scale delta integral in general. Here `discrete_integral` sums index values without graininess factors; arbitrary closed subsets of real numbers, jump operators and dynamic equations are **outside the scope** of this module.

## Independent mathematical references

- [NIST Digital Library of Mathematical Functions, §24.4(iii), Sums of Powers](https://dlmf.nist.gov/24.4.iii) — standard Bernoulli-polynomial power sum relations.
- [MIT OpenCourseWare, 18.303, notes on finite differences](https://ocw.mit.edu/courses/18-303-linear-partial-differential-equations-analysis-and-numerics-fall-2014/edea18830013f8c73172ae6bc0eeaafb_MIT18_303F14_difference.pdf) — forward, backward and centered derivative approximations and accuracy caveats.

The descriptions, formulas, tests and implementation here are original; no source implementation or lesson text has been copied.
