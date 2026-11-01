# Exact polynomial integration over Q

This module adds a bounded, exact integration API for rational univariate
polynomials. It uses the repository's existing `calc.algebra.Polynomial`
representation: coefficients are ascending powers of `x`, and all coefficients,
endpoints, and integration constants are integers or `fractions.Fraction`.

```python
from fractions import Fraction
from calc.algebra import Polynomial
from calc.integrals import (
    polynomial_antiderivative, definite_integral,
    substitution_integral, integration_by_parts,
)

p = Polynomial((1, 0, 3))      # 1 + 3*x^2
F = polynomial_antiderivative(p)
assert F.differentiate() == p
assert definite_integral(p, 0, 1) == 2

# An oriented integral negates on reversal.
assert definite_integral(p, 1, 0) == -2

# Substitution: integral of p(q(x)) q'(x), q=x^2+1.
outer = Polynomial((0, 0, 1))  # p(y)=y^2
inner = Polynomial((1, 0, 1))  # q(x)=1+x^2
assert substitution_integral(outer, inner, 0, 1) == Fraction(7, 3)

# ∫ u'(x) v(x) dx, computed from [u*v] - ∫u*v'.
u = Polynomial((0, 0, 1))      # x^2
v = Polynomial((0, 0, 0, 1))   # x^3
assert integration_by_parts(u, v, 0, 2) == Fraction(64, 5)
```

## Mathematical contract

For a polynomial \(p(x)=\sum_{k=0}^n a_kx^k\),

\[
P(x)=C+\sum_{k=0}^n\frac{a_k}{k+1}x^{k+1},\quad P'(x)=p(x).
\]

The definite integral is exactly \(P(b)-P(a)\), including rational and
reversed bounds. The polynomial is continuous everywhere, so this use of the
Fundamental Theorem of Calculus has no exceptional discontinuity domain.

The change-of-variable utility accepts **two separate polynomials** `outer=p`
and `inner=q` and computes

\[
\int_a^b p(q(x))q'(x)\,dx=P(q(b))-P(q(a)).
\]

It does not accept a pre-composed integrand or infer a substitution from an
arbitrary expression. It avoids expanding \(p(q(x))q'(x)\).

The integration-by-parts API uses an explicit order:

\[
\texttt{integration_by_parts}(u,v,a,b)
=\int_a^b u'(x)v(x)\,dx
=[u(x)v(x)]_a^b-\int_a^b u(x)v'(x)\,dx.
\]

In this restricted setting the correction polynomial is integrated exactly,
not approximated by quadrature.

## Resource and precision boundaries

- Maximum input polynomial degree: 64.
- Maximum input rational numerator/denominator: 4,096 bits.
- Maximum computed rational numerator/denominator: 65,536 bits, with a
  conservative pre-evaluation estimate to reject explosively large endpoint
  exponentiation.
- Compositions and polynomial products also have a degree budget.
- Inexact floats, booleans, expressions, nonpolynomial functions, nonfinite
  bounds, and general symbolic integration are deliberately unsupported.
- Results are exact Fractions for supported arguments. There is no numerical
  truncation error, but resource-bound rejection is part of the API contract.

For approximate integration of callable functions, see `calc.numerical`,
which has different mathematical guarantees and error considerations.

## Verification

`tests/test_integrals.py` uses independently known elementary antiderivatives,
80 seeded monomial closed-form witnesses, reversed-bound cancellation,
linearity, 60 seeded integration-by-parts comparisons, substitution examples,
and invalid/resource cases. Tests demonstrate bounded operations, not a formal
proof of an unrestricted integration algorithm.

## Primary mathematical references

- OpenStax, *Calculus Volume 1*, [§5.3: Fundamental Theorem of Calculus](https://openstax.org/books/calculus-volume-1/pages/5-3-the-fundamental-theorem-of-calculus).
- OpenStax, *Calculus Volume 2*, [§3.1: Integration by Parts](https://openstax.org/books/calculus-volume-2/pages/3-1-integration-by-parts).
- OpenStax, *Calculus Volume 1*, [§5.5: Substitution](https://openstax.org/books/calculus-volume-1/pages/5-5-substitution).

These references explain the mathematical identities. The implementation and
tests are original and do not copy upstream code. No reuse license is implied.
