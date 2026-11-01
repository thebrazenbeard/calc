# Exact matrix Lie algebra

This module handles **finite-dimensional rational matrix representations** of
Lie algebras. A Lie algebra is a vector space with a bilinear alternating
bracket satisfying the Jacobi identity. Square matrices give a concrete example
using the commutator

\[
[A,B]=AB-BA.
\]

The ambient matrix algebra has an *associative* product, but its Lie bracket is
not associative. Rather, the commutator obeys

\[
[[A,B],C]+[[B,C],A]+[[C,A],B]=0.
\]

## Quick start

```python
from fractions import Fraction
from calc.lie import (
    sl2_basis, so3_basis, commutator, structure_constants,
    adjoint_matrix, basis_coordinates, linear_combination,
)

H, E, F = sl2_basis()
assert commutator(H, E) == linear_combination((2,), (E,))
assert commutator(H, F) == linear_combination((-2,), (F,))
assert commutator(E, F) == H

constants = structure_constants((H, E, F))
assert constants[0][1] == (0, 2, 0)
assert constants[1][2] == (1, 0, 0)

ad_H = adjoint_matrix(H, (H, E, F))
assert ad_H == ((0, 0, 0), (0, 2, 0), (0, 0, -2))

X = linear_combination((Fraction(1, 3), 2, -5), (H, E, F))
assert basis_coordinates(X, (H, E, F)) == (Fraction(1, 3), 2, -5)

Lx, Ly, Lz = so3_basis()
assert commutator(Lx, Ly) == Lz
```

All returned matrices and structure constants are immutable tuples containing
`fractions.Fraction` scalars. Integers are normalized exactly. Floats, complex
numbers, booleans, and non-square matrices are rejected. The fixed-size limits
are **8×8 matrices**, **16 basis vectors**, and **2048-bit input rational
numerators/denominators**.

## Basis and structure constant semantics

For a chosen ordered basis \(e_0,\ldots,e_{r-1}\), the structure constants
obey

\[
[e_i,e_j] = \sum_{k=0}^{r-1} c_{ij}^{\;k} e_k.
\]

`structure_constants(basis)` performs exact Gaussian elimination over rational
numbers. It rejects linearly dependent bases and any basis span that is **not
closed under commutators**. Closure matters: arbitrary subspaces of a matrix
space do not necessarily form Lie subalgebras.

`adjoint_matrix(X, basis)` requires `X` to belong to the same closed subalgebra.
Its **column `j`** is the coordinate vector of \([X,e_j]\). The routine does
not infer a canonical basis or choose a numerical tolerance.

The built-in \(\mathfrak{sl}_2\) basis is ordered `(H, E, F)` with
traceless 2×2 matrices. The \(\mathfrak{so}_3\) basis is ordered `(Lx, Ly, Lz)`
with real skew-symmetric 3×3 matrices. These provide computational examples of
the relations; the library does not classify other Lie algebras or construct
continuous Lie groups.

## Limits and verification

- Rational input means exact algebra for every supported operation, but output
  numerator/denominator sizes can grow; dimensions are bounded for this reason.
- Matrix products are square and equal-dimension. `basis_coordinates` raises
  on targets outside the rational span; it is not a least-squares solver.
- No matrix exponential/logarithm, Baker–Campbell–Hausdorff expansion, Killing
  form, integration to Lie groups, or representation-classification claims.
- `tests/test_lie.py` verifies the Jacobi identity, bracket antisymmetry,
  \(\mathfrak{sl}_2\) and \(\mathfrak{so}_3\) relations, structure constants,
  adjoint convention, rational coefficients, and invalid inputs.

## References

1. Pavel Etingof, *Lie Groups and Lie Algebras I & II*, [MIT OpenCourseWare
   18.755, full lecture notes](https://ocw.mit.edu/courses/18-755-lie-groups-and-lie-algebras-ii-spring-2024/mit18_755_s24_lec_full.pdf).
2. [Encyclopedia of Mathematics: Lie algebra](https://encyclopediaofmath.org/wiki/Lie_algebra),
   definition and the Jacobi identity.

The examples and implementation here are original; the sources are references
to definitions and standard identities, not copied code.
