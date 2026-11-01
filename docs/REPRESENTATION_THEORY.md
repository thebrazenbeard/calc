# Finite-group representation theory over the rationals

The `calc.representations` module provides a **bounded exact example**, not a
general representation-theory engine. It constructs small groups from complete
multiplication tables and checks that assigned rational square matrices form a
homomorphism into `GL_d(ℚ)`.

## Mathematics

For a finite group \(G\), a representation \(\rho\) satisfies

\[
\rho(e)=I,\qquad \rho(gh)=\rho(g)\rho(h).
\]

The implementation verifies these conditions for **every pair** of group
elements with exact `fractions.Fraction` arithmetic. The group constructor
first checks closure, both identity laws, two-sided inverses, and all ordered
triples for associativity.

The character is \(\chi_\rho(g)=\operatorname{tr}\rho(g)\). For two
rational representations of the same group, the finite sum

\[
\langle\chi,\psi\rangle_G=
  \frac{1}{|G|}\sum_{g\in G}\chi(g)\psi(g^{-1})
\]

is exact. It is useful for examples of character orthogonality; **do not
interpret arbitrary values as an irreducibility certificate over ℚ**.
Standard irreducible character orthogonality is normally stated over a
splitting field such as \(\mathbb C\), and rational irreducibility has
additional subtleties.

The left regular representation has a basis indexed by \(G\) and
\(\rho(g)e_h=e_{gh}\). Its character is \(|G|\) at the identity and zero
elsewhere. The tests also construct the nonabelian symmetric group \(S_3\)
independently from permutation composition and verify its regular action.

## Executable examples

```python
from fractions import Fraction
from calc.representations import (
    cyclic_group, RationalRepresentation, trivial_representation,
    regular_representation, character_pairing, direct_sum, intertwines,
)

g = cyclic_group(2)
trivial = trivial_representation(g)
sign = RationalRepresentation(g, {"e": ((1,),), "g": ((-1,),)})
regular = regular_representation(g)

assert regular.matrix("g") == ((0, 1), (1, 0))
assert regular.character() == {"e": Fraction(2), "g": Fraction(0)}
assert character_pairing(trivial, sign) == 0
assert character_pairing(regular, trivial) == 1

# The C2 regular action is isomorphic to trivial ⊕ sign.
sum_rep = direct_sum(trivial, sign)
assert intertwines(sum_rep, regular, ((1, 1), (1, -1)))
```

`FiniteGroup(elements, multiplication, identity)` also accepts a caller's
fully specified multiplication table. Invalid tables, wrong dimensions,
nonhomomorphic assignments, and inexact floating-point coefficients fail
with explicit exceptions.

## Exact limits and performance

- Up to **16 group elements**, **16-dimensional** rational matrices, and
  **4096-bit** individual numerator/denominator budgets.
- Group validation is cubic in group order; checking all pairs of matrices
  is approximately \(O(|G|^2d^3)\).
- No complex coefficients, general field extensions, arbitrary character
  tables, representation induction, symbolic decomposition, or infinite groups.
- `intertwines` checks an explicitly supplied linear map. It **does not**
  search for all intertwiners or guarantee an isomorphism merely because a
  noninvertible equivariant map exists.
- Source code and API examples here are original. The external material below
  is linked for mathematical reference, not copied into the repository.

## Mathematical reference

- P. Etingof et al., *Introduction to Representation Theory*, MIT OpenCourseWare
  18.712 (Fall 2010), [Chapter 3: Representations of Finite Groups — Basic
  Results](https://ocw.mit.edu/courses/18-712-introduction-to-representation-theory-fall-2010/resources/mit18_712f10_ch3/).
  Sections 3.2–3.5 cover characters and orthogonality.

Run `PYTHONPATH=src python -m unittest tests.test_representations -v` for
the independent exact C₂/C₃/S₃ examples, 19 tests, and group-law
rejection cases.
