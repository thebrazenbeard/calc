# Exact computational algebra: Gröbner bases

This module adds a **small, testable demonstration** of computational commutative algebra. It does not implement a general-purpose computer algebra system. Calculus appears only indirectly: multivariate polynomial ideals, elimination, and symbolic polynomial relations are tools that later connect to differential algebra and algebraic geometry.

## Mathematical contract

Work in the polynomial ring \(\mathbb Q[x_0,\ldots,x_{n-1}]\), with \(1\le n\le8\). Coefficients are Python integers or `fractions.Fraction`; floats and booleans are rejected. Sparse terms use exponent tuples:

\[
3x_0^2x_1-\frac12x_1^3 \quad\leftrightarrow\quad
\{(2,1):3,\ (0,3):-\tfrac12\}.
\]

A monomial order determines which term is *leading*. The implementation supports:

- `lex`: ordinary lexicographic comparison of exponents, with \(x_0>x_1>\cdots\).
- `grlex`: first compare total degree, then lexicographically.

For monomials \(x^a\), \(x^b\), divisibility means \(a_i\le b_i\) in every coordinate. The multivariate division algorithm repeatedly cancels a leading term if an allowed divisor's leading monomial divides it; otherwise it moves that term to the remainder. For an arbitrary list of divisors the remainder may depend on their order. **When the divisors form a Gröbner basis**, the normal form is unique for that monomial order.

For nonzero polynomials \(f,g\), their S-polynomial cancels the leading terms using the least common multiple of the leading monomials. Buchberger's criterion says a finite set is a Gröbner basis precisely when every pair's S-polynomial reduces to zero with respect to that set. The module applies this criterion to construct a basis, then minimizes leading terms and reduces its tails. No alternate F4/F5, FGLM conversion, modular methods, or numerical conditioning are provided.

## Original worked example

Consider the ideal

\[
I=\langle x^2-y,\ xy-1\rangle \subseteq \mathbb Q[x,y],
\quad x>y \text{ (lex)}.
\]

The first equation gives \(x^2=y\); multiplying the second by \(x\) gives \(x^2y=x\), hence \(x=y^2\). Substitution into \(xy=1\) gives \(y^3=1\). The reduced lex Gröbner basis is therefore

\[
G=\{x-y^2,\ y^3-1\}.
\]

This is not a floating-point guess; all coefficients and reductions remain in \(\mathbb Q\).

```python
from calc.grobner import variables, grobner_basis, normal_form

x, y = variables(2)
generators = (x*x - y, x*y - 1)
G = grobner_basis(generators, order="lex")

assert set(G) == {x-y*y, y*y*y-1}
assert not normal_form(x*x-y, G)
assert normal_form(x+y, G) == y*y+y
```

The tests independently check generator membership, S-polynomial reductions, lex versus graded lex, rational coefficients, polynomial identities, and a separate three-variable system.

## Public API

| Name | Behavior |
| --- | --- |
| `MultiPolynomial(mapping, nvars=...)` | Sparse immutable exact multivariate polynomial |
| `variables(nvars)` | Canonical ordered indeterminates |
| `leading_term(poly, order)` | Leading monomial/coefficient for `lex` or `grlex` |
| `s_polynomial(f, g, order)` | Leading-term cancellation polynomial |
| `normal_form(f, divisors, order, ...)` | Exact multivariate long-division remainder |
| `grobner_basis(generators, order, ...)` | Monic reduced basis via Buchberger, within budgets |

`MultiPolynomial` supports addition, subtraction, multiplication, nonnegative powers, and exact evaluation. Scalar integer/fraction coefficients are accepted in arithmetic and as generators once another polynomial establishes `nvars`.

## Failure and resource boundaries

- **Not a theorem prover:** Tests establish selected algebraic contracts, not the universal correctness of arbitrary ideal problems.
- **No inexact coefficients:** Floating point is rejected because Buchberger reduction depends on exact zero tests.
- **Not a generic polynomial parser:** Callers construct term mappings or arithmetic on `MultiPolynomial` objects; no `eval`, `exec`, or third-party notebook code is used.
- **Memory/time constraints:** Maximum eight variables; an exponent in one monomial can be at most 512; at most 4,096 sparse terms per polynomial; at most 100,000 multiplication-term pairs; coefficient numerator/denominator bit length at most 4,096. The caller can further restrict Buchberger pairs, division steps, basis size, and intermediate term count. Limits raise `ResourceLimitError`, not partial success.
- **No proof of solution geometry:** An ideal membership or elimination result does not enumerate all complex or real solutions, multiplicities, or singularities.
- **Fixed monomial order:** A basis computed for `lex` is not automatically a Gröbner basis for `grlex`. Recompute for the desired order.
- **Potential extreme complexity:** Gröbner computations can grow dramatically in degree and size. This educational implementation deliberately fails closed under its explicit budgets.

## Source comparisons

- [SymPy polynomial tools reference](https://docs.sympy.org/latest/modules/polys/reference.html) documents its `groebner` interface and monomial-order options. This implementation is independently written and does not import SymPy.
- [SageMath multivariate polynomial ideals](https://doc.sagemath.org/html/en/reference/polynomial_rings/sage/rings/polynomial/multi_polynomial_ideal.html) is a mature comparison target for larger problems and explicitly discusses reduced bases and ideal operations.
- [ETH Zürich course notes: *Effective Methods in Polynomial Rings*](https://metaphor.ethz.ch/x/2025/hs/401-3132-00L/ex/Groebnerbases.pdf) provides Buchberger's criterion and the S-polynomial construction at proof level.
- [University of Illinois at Chicago lecture notes](https://homepages.math.uic.edu/~jan/mcs563s14/lec06.pdf) separately states the S-polynomial remainder criterion.

These links are mathematical references, not code imports or copied educational content. License review is necessary before incorporating source text or implementation fragments from them.
