# Classical umbral calculus over exact rational polynomials

The module `calc.umbral` gives a small, test-backed implementation of **finite operator calculus**: shifts, forward differences, falling-factorial polynomial bases, signed/unsigned Stirling transforms, and Newton forward interpolation. Every coefficient is a Python `int` or `fractions.Fraction`, and every result is exact rational arithmetic. This module does **not** purport to implement the entire modern theory of umbral operators or Sheffer sequences.

## Basic examples

```python
from fractions import Fraction
from calc.umbral import (
    falling_factorial, rising_factorial, forward_difference,
    stirling_first, stirling_second, monomial_to_falling,
    falling_to_monomial, newton_forward_interpolate,
)

assert falling_factorial(5, 3) == 60            # 5 * 4 * 3
assert rising_factorial(3, 4) == 360            # 3 * 4 * 5 * 6
assert stirling_first(4, 2) == 11               # signed first kind
assert stirling_second(4, 2) == 7

# Coefficients are ascending: (0, 0, 0, 1) is x**3.
assert forward_difference((0, 0, 0, 1)) == (1, 3, 3)
assert forward_difference((0, 0, 0, 1), order=3) == (6,)

# x**4 = x_(1) + 7*x_(2) + 6*x_(3) + x_(4).
falling = monomial_to_falling((0, 0, 0, 0, 1))
assert falling == (0, 1, 7, 6, 1)
assert falling_to_monomial(falling) == (0, 0, 0, 0, 1)

# Values at x = -2, 1, 4, 7 for p(x) = x**2 + 4*x - 2.
samples = (-6, 3, 30, 75)
assert newton_forward_interpolate(samples, start=-2, step=3) == (-2, 4, 1)

# Rational-grid interpolation remains rational.
assert newton_forward_interpolate((1, Fraction(3, 2)),
                                  start=0, step=Fraction(1, 2)) == (1, 1)
```

## Falling factorials and Stirling basis change

For nonnegative integers \(n\), define

\[
x_{\underline n}=\prod_{j=0}^{n-1}(x-j),\qquad
x^{\overline n}=\prod_{j=0}^{n-1}(x+j),\qquad
x_{\underline 0}=x^{\overline 0}=1.
\]

The **signed** first-kind Stirling numbers \(s(n,k)\) and second-kind numbers \(S(n,k)\) satisfy

\[
x_{\underline n}=\sum_{k=0}^n s(n,k)x^k,
\qquad
x^n=\sum_{k=0}^n S(n,k)x_{\underline k}.
\]

The two coefficient transformations are mutual inverses over \(\mathbb Q[x]\). In particular, \(s(4,k)=(0,-6,11,-6,1)\) while \(S(4,k)=(0,1,7,6,1)\). Tests compare these rows to independently known values and check inverses on many nontrivial rational polynomials.

Both `monomial_to_falling` and `falling_to_monomial` receive *ascending* coefficient tuples. Only the **basis** changes; tuple position \(k\) continues to identify degree \(k\). Trailing zero coefficients are removed, and a zero polynomial is represented by `(Fraction(0),)`.

## Forward difference as an operator

Let \(E_h p(x)=p(x+h)\), and let \(\Delta_h=E_h-I\). The public function

```python
forward_difference(coefficients, step=h, order=m)
```

returns the monomial coefficients of \(\Delta_h^m p\). Unlike a numerical derivative, there is **no division by \(h\)**. For unit shift, the fundamental factorial-basis lowering identity is

\[
\Delta_1 x_{\underline n}=n\,x_{\underline{n-1}}.
\]

For a general nonzero rational step, \(\Delta_h x^3=3h x^2+3h^2x+h^3\). Higher differences eventually vanish identically on a polynomial, and `order=0` leaves coefficients unchanged. A negative rational step is valid. This is exact manipulation of a polynomial, not taking differences of already rounded samples.

## Newton forward interpolation

For samples \(y_j=p(x_0+jh)\) on a nonzero rational regular grid, the unique polynomial of degree at most \(N-1\) through \(N\) points is

\[
p(x)=\sum_{k=0}^{N-1}\frac{\Delta^k y_0}{k!}
\prod_{j=0}^{k-1}\frac{x-x_0-jh}{h}.
\]

`newton_forward_interpolate(values, start=x0, step=h)` computes this polynomial in ascending *monomial* coefficients. For example, four exact samples of a cubic reconstruct the cubic exactly. Interpolation is a finite algebraic identity, not a convergence guarantee outside the supplied nodes or a fit to noisy observations.

## Explicit bounds and what is not implemented

- Scalar inputs must be `int` or `Fraction`. Floats, booleans, complex numbers and automatically rounded conversions are rejected.
- Degree, factorial order and Stirling row index are at most **32**. Interpolation accepts from **1 to 33** sample values. The iterable readers stop after at most 34 values, so an infinite iterator is safely rejected rather than fully consumed.
- Each supplied rational numerator and denominator is limited to **4,096 bits**. The implementation does not certify an output-memory ceiling under every composition of large rationals; callers should treat this as bounded educational computation, not an untrusted high-throughput algebra service.
- Interpolation nodes must be an arithmetic progression with nonzero rational spacing. Irregular nodes, Hermite multiplicities, least-squares fitting, interpolation-error bounds and floating-point conditioning are outside scope.
- No claims are made about arbitrary Sheffer sequences, general delta operators, q-umbral variants or infinite formal power series.

## Primary mathematical references

- [NIST Digital Library of Mathematical Functions, §26.8 — Stirling numbers](https://dlmf.nist.gov/26.8) gives the standard first- and second-kind Stirling numbers and their factorial-basis identities.
- [NIST DLMF, §3.3 — Interpolation](https://dlmf.nist.gov/3.3) presents Newton polynomial interpolation and divided differences, of which equispaced Newton forward interpolation is a special case.

These descriptions and examples were independently written for the repository. Sources are linked for mathematical provenance, not copied as implementation code.
