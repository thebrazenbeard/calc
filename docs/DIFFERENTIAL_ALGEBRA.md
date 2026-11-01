# Differential algebra: bounded rational jet polynomials

This module gives one **exact, finite example** of differential algebra: a
commutative rational polynomial ring in formal jet symbols

\[
  \mathbb{Q}[u_0,u_1,\ldots,u_m], \qquad u_j\text{ denotes a formal }j\text{th derivative}.
\]

The total derivation is the formal operator

\[
  D=\sum_{j=0}^{m-1} u_{j+1}\frac{\partial}{\partial u_j},
  \qquad D(fg)=D(f)g+fD(g).
\]

For example, `D(u_0^2 u_1) = 2 u_0 u_1^2 + u_0^2 u_2`.
The coefficients stay in `fractions.Fraction`; there are no floating-point
approximations hidden in the polynomial operations.

## Runnable example

From the repository root:

```python
from fractions import Fraction
from calc.differential_algebra import JetPolynomial

u = JetPolynomial.jet(0, max_jet=4)  # u
ux = JetPolynomial.jet(1, max_jet=4) # u'
uxx = JetPolynomial.jet(2, max_jet=4) # u''

f = Fraction(3, 2) * u**2 * ux
df = f.total_derivative()
assert df == 3 * u * ux**2 + Fraction(3, 2) * u**2 * uxx
assert (u * ux).total_derivative() == ux**2 + u * uxx

# Exact rational substitution supplies numerical values for *jet symbols*,
# not proof that they come from a particular differentiable function.
assert df.evaluate((2, 3, 4, 0, 0)) == 3*2*9 + Fraction(3, 2)*4*4
```

Polynomials are stored in an immutable sparse canonical representation. A
monomial's exponent tuple is ordered as `(power(u_0), power(u_1), ...)`.
For instance:

```python
f = JetPolynomial(2, {(2, 1, 0): Fraction(3, 4)})
assert f == Fraction(3, 4) * JetPolynomial.jet(0, max_jet=2)**2 * JetPolynomial.jet(1, max_jet=2)
```

The `partial_derivative(j)` operation is differentiation with respect to one
formal jet variable. It differs from the total derivation `total_derivative()`:

\[
  \frac{\partial}{\partial u_0}(u_0 u_1)=u_1,
  \qquad D(u_0 u_1)=u_1^2+u_0 u_2.
\]

## Finite-jet boundary is explicit

The full differential polynomial ring ordinarily includes arbitrarily high
jet orders. Our ring represents only `u_0` through `u_m`. Its total
derivation is therefore **partial as an operation on the truncated ring**:
if a polynomial depends on `u_m`, the true derivative needs `u_(m+1)` and
raises `ValueError`. No highest-jet contribution is silently removed.

```python
top = JetPolynomial.jet(2, max_jet=2)
try:
    top.total_derivative()
except ValueError:
    pass  # u_3 is unavailable in this finite ring
else:
    raise AssertionError("truncation must not discard the missing derivative")
```

Additional finite limits protect the implementation:

| Constraint | Bound |
| --- | --- |
| Highest jet index | 8 |
| Total monomial degree | 16 |
| Nonzero monomials | 1,024 |
| Pair products in one multiplication | 20,000 |
| Numerator/denominator bit length | 4,096 each |

Only integer and `Fraction` coefficients are accepted. Zero terms are
normalized away. No infinite generators are consumed by `evaluate`, and
out-of-ring operations fail explicitly. The code does **not** implement
differential elimination, Ritt–Kolchin algorithms, solving ODEs/PDEs, a
differential field, analytic derivatives, or full arbitrary-jet closure.

## Independent mathematical checks

`tests/test_differential_algebra.py` checks the Leibniz and linearity
identities, polynomial chain and power examples, 70 seeded rational Leibniz
witnesses, repeated jet advancement, exact substitution, commutativity,
truncation errors, input rejection and resource limits.

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_differential_algebra.py' -v
```

## Research sources and attribution

- [Reid Dale, *Notes on Differential Algebra* (UC Berkeley)](https://math.berkeley.edu/~reiddale/differential_algebra_notes.pdf) — differential polynomial rings and the extension of derivations by the product rule.
- [Gleb Pogudin, *Differential Algebra* lecture notes (École Polytechnique)](https://www.lix.polytechnique.fr/Labo/Gleb.POGUDIN/files/da_notes.pdf) — differential rings and derivations.
- [Joseph Fels Ritt, *Differential Algebra* (American Mathematical Society)](https://doi.org/10.1090/coll/033) — historical monograph on differential algebra.

These are references, not copied implementation sources. Their individual
copyright and reuse terms remain separate from this repository's license.
