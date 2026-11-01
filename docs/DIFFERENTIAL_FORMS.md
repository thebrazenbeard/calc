# Exact exterior forms and the Euclidean Hodge star

`calc.forms` implements the finite-dimensional *exterior algebra* of an
oriented, Euclidean vector space with a fixed **orthonormal basis**. Its
coefficients are exact rational numbers (`fractions.Fraction`). The code uses
zero-based basis indices and a sparse representation, so it can manipulate
forms without constructing the entire `2**n`-dimensional algebra.

**Scope:** This is linear algebra on forms with constant coefficients, **not**
an implementation of the exterior derivative on fields, chart transitions,
smooth manifolds, curved metrics, differential geometry, de Rham cohomology,
connections, or Hodge theory for arbitrary Riemannian manifolds. In particular,
the Hodge star formulas here require an *oriented orthonormal* basis.

## Quick example

```python
from fractions import Fraction
from calc.forms import Form, basis_form, hodge_star, wedge, inner_product, norm_squared

dx = basis_form(3, (0,))
dy = basis_form(3, (1,))
dz = basis_form(3, (2,))

assert wedge(dx, dy) == basis_form(3, (0, 1))
assert wedge(dy, dx) == -wedge(dx, dy)
assert wedge(dx, dx) == Form(3, 2)  # Zero 2-form

assert hodge_star(dx) == wedge(dy, dz)
assert hodge_star(dy) == -wedge(dx, dz)
assert hodge_star(wedge(dx, dy)) == dz

alpha = Form(3, 1, {(0,): Fraction(1, 2), (2,): -3})
assert norm_squared(alpha) == Fraction(37, 4)
assert inner_product(alpha, dz) == -3

# Reversing orientation flips the Hodge dual, without changing the metric.
assert hodge_star(dx, orientation=-1) == -wedge(dy, dz)
```

A `Form(n, k, terms)` represents a homogeneous `k`-form in `n` dimensions;
`terms` is a dictionary from **strictly increasing tuples** of `k` indices to
integer or `Fraction` coefficients. For example:

```python
Form(4, 2, {(0, 2): Fraction(2, 3), (1, 3): -5})
```

is \(\frac23\,dx_0\wedge dx_2-5\,dx_1\wedge dx_3\). The sorted sparse
representation is immutable; `.as_dict()` returns a separate mutable copy.
Missing basis elements are zero, and explicitly zero coefficients disappear.
`Form(n, k)` creates the zero form. The empty tuple is the 0-form basis
(scalar unit): `Form(3, 0, {(): 1})`.

## Algebraic rules

For an ordered exterior basis element
\(e_I=dx_{i_1}\wedge\cdots\wedge dx_{i_k}\), with
\(i_1<\cdots<i_k\):

- `wedge(a, b)` implements the exterior product, expanding bilinearly. It is
  associative and graded-commutative:
  \(\alpha\wedge\beta=(-1)^{pq}\beta\wedge\alpha\) for degrees \(p,q\).
- If the two basis monomials share an index, their wedge is zero. Otherwise,
  the coefficient's sign is the parity of the swaps necessary to put the
  concatenated indices in ascending order.
- The wedge product of degrees summing beyond `n` is the zero form in that
  resulting degree. A zero form may have degree greater than `n`; `hodge_star`
  rejects it because its Hodge degree \(n-k\) would be negative.
- Forms of equal degree can be added, subtracted, negated, multiplied by an
  integer or `Fraction`, and divided by a nonzero integer or `Fraction`.
- `inner_product(a, b)` uses the basis wedge monomials as an orthonormal basis:
  \(\langle\alpha,\beta\rangle=\sum_I \alpha_I\beta_I\). Forms of
  different degree are orthogonal; different ambient dimensions raise an
  error. `norm_squared(a)` returns \(\langle a,a\rangle\), remaining exact.

## Hodge dual, sign, and orientation

Write `n` for the ambient dimension and `k` for the input degree.
`hodge_star(alpha)` returns the `(n-k)`-form \(\star\alpha\). On a basis
monomial, let \(I^c\) be its **ordered complement** in `(0, ..., n-1)`.
Then

\[
\star e_I = \operatorname{sgn}(I,I^c)\,e_{I^c},
\]

where `sgn` is the sign of the permutation obtained by concatenating `I` and
`I^c`. This fixes the standard orientation. Setting `orientation=-1` negates
the entire Hodge dual, corresponding to the opposite orientation. Only `+1`
and `-1` are accepted. For the chosen orientation and Euclidean metric,

\[
\alpha\wedge\star\beta
  =\langle\alpha,\beta\rangle\,\operatorname{vol},
\qquad
\star\star\alpha=(-1)^{k(n-k)}\alpha,
\]

where the first identity holds when `alpha` and `beta` have equal degree and
\(\operatorname{vol}=dx_0\wedge\cdots\wedge dx_{n-1}\) in standard
orientation (its negative in the reversed orientation). Both calls in the
second identity must use the **same** orientation.

For example, in a 2D orthonormal frame, \(\star dx_0=dx_1\) and
\(\star dx_1=-dx_0\), so \(\star^2=-1\) on 1-forms. In 3D, it is `+1`
on all degrees. These are properties of Euclidean signature: Lorentzian or
other indefinite metrics require modified signs and are not supported.

## Input validation and limits

- Dimensions must be positive integers and degrees nonnegative integers.
- Indices must be tuples of the stated degree, strictly increasing, unique,
  integer-valued (not Boolean), and in the range `0..n-1`.
- Coefficients are only `int` or `Fraction`. Floats, even if they happen to
  represent an exact integer, are deliberately rejected to avoid silent
  approximation.
- The constructor copies and normalizes its input; callers cannot mutate
  stored terms by editing the original mapping or `.as_dict()` results.
- The implementation is sparse but **not** complexity-free: an exterior
  product may involve every pair of nonzero terms and the result can still
  have up to `binomial(n, k)` basis coefficients. It is not intended for huge
  dense algebras.

## Tests

```sh
PYTHONPATH=src python -m unittest discover -s tests -p test_forms.py -v
```

The suite checks basis signs, wedge associativity/anticommutativity and
nilpotence, rational coefficients, orientation, star-square signs across
dimensions/degrees, metric pairing, grade overflow, immutable normalization,
and malformed inputs.

## Mathematical references

- [MIT 18.155 differential forms / Hodge star lecture notes](https://math.mit.edu/~dyatlov/18.155/155-notes.pdf)
- [Stanford: The Hodge star operator](https://math.stanford.edu/~conrad/diffgeomPage/handouts/star.pdf)
- [UC San Diego: Discrete Differential Geometry notes](https://cseweb.ucsd.edu/~alchern/teaching/DDG.pdf)

These references motivate the finite-dimensional algebraic conventions.
Their broader manifold-level mathematics is outside this module's scope.
