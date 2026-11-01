# Free associative noncommutative algebra over ℚ

`calc.free_algebra` supplies a deliberately small and exact model of a free
associative algebra. It represents elements of

\[
\mathbb Q\langle x_1,\ldots,x_m\rangle
=\bigoplus_{n\ge0} V^{\otimes n},
\]

where each ordered word in the generators is an independent monomial. The
empty word is the multiplicative identity. Multiplication concatenates words
and distributes across sums. In particular, **`xy` and `yx` do not coincide**
without an explicitly supplied relation; there is no automatic simplification
to commutative polynomials.

## Executable example

```python
from fractions import Fraction
from calc.free_algebra import FreeAlgebra

A = FreeAlgebra(("x", "y"))
x, y = A.variable("x"), A.variable("y")

assert x*y != y*x
assert (x+y)**2 == x*x + x*y + y*x + y*y
assert (x*y*x).derivative("x") == y*x + x*y
assert (Fraction(2, 3)*x + Fraction(1, 3)*x) == x

assert x.commutator(y) == x*y - y*x
assert (x.commutator(y)).commutator(x) + \
       (y.commutator(x)).commutator(x) + \
       (x.commutator(x)).commutator(y) == A.zero()
```

This last specially chosen Jacobi sum is zero, but `tests/test_free_algebra.py`
also uses three distinct generators for an independent Jacobi check.

## Formal derivations

For a selected generator \(x\), the method `P.derivative("x")` is the unique
derivation \(D_x\) taking \(D_x(x)=1\) and all other generators to zero.
The definition on words is

\[
D_x(a_1\cdots a_n)
=\sum_{i:a_i=x} a_1\cdots a_{i-1}a_{i+1}\cdots a_n.
\]

Every term preserves its surviving word order. The Leibniz rule holds:

\[
D_x(PQ)=D_x(P)Q+P D_x(Q).
\]

`P.derivation({"x": image_x, "y": image_y})` instead extends an arbitrary
specified map from generators to algebra elements to a derivation. Omitted
generators have zero image. This is an algebraic derivation, not a pointwise
numerical derivative or a noncommutative difference quotient. It is also not
the related *free difference quotient* mapping into a tensor product.

The commutator \([P,Q]=PQ-QP\) forms a Lie bracket because the underlying
concatenation product is associative. Antisymmetry and Jacobi hold exactly.

## Encoding and resource limits

- Generator names are short ASCII identifiers; 1–8 distinct generators.
- Words are tuples of generator names, at most 64 symbols in each word.
- Coefficients must be `int` or `fractions.Fraction`; floats and booleans
  are rejected. Numerator and denominator each have at most 4,096 bits.
- A polynomial contains at most 4,096 nonzero sparse terms. Multiplication
  and generalized derivation limit pair operations to 100,000. Nonnegative
  integer exponentiation is capped at 32.
- Arithmetic is exact within those limits. The resource ceilings protect
  examples from combinatorial explosion; they are not mathematical degree
  theorems. A deliberately empty zero polynomial is distinct from the unit.
- Polynomials are immutable; `as_dict()` returns a fresh copy. Cross-algebra
  operations are rejected unless their ordered generator lists agree.

This module does **not** provide noncommutative Gröbner bases, quotienting by
relations, universal enveloping algebras, completed tensor algebras, general
representations, differential operator algebras, or a classification of
noncommutative rings. The separate `calc.grobner` module concerns *commutative*
multivariate rational polynomials and must not be silently reused as an
algorithm for this free algebra.

## References and provenance

- Massachusetts Institute of Technology, *Noncommutative Algebra (18.706)*,
  lecture notes: https://ocw.mit.edu/courses/18-706-noncommutative-algebra-spring-2023/lists/lecture-notes/
- MIT Mathematics, *Free associative algebras*:
  https://math.mit.edu/~dav/freeassoc.pdf

These are background references, not copied implementation sources. This
repository does not inherit permission to redistribute third-party notes or
code merely by linking to them; attribution is separate from licensing.
