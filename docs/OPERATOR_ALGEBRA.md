# Finite rational operator algebras

This example is the exact, finite-matrix corner of **operator algebra**. It
works over rational entries in a real orthonormal coordinate basis. Its adjoint
is the ordinary transpose, not a conjugate transpose over complex numbers.
It does **not** implement infinite-dimensional bounded operators, Hilbert-space
functional calculus, spectral measures, C\*-algebra completions, or arbitrary
von Neumann algebras.

## Algebraic rules

For square rational matrices of the same size,

\[
A^* = A^T,\qquad [A,B]=AB-BA,\qquad \{A,B\}=AB+BA.
\]

The involution satisfies \((AB)^*=B^*A^*\), and the trace obeys
\(\operatorname{tr}(AB)=\operatorname{tr}(BA)\). For a real matrix, a standard
orthogonal projection is exactly a self-adjoint idempotent:
\(P^T=P=P^2\). A matrix is normal when \(AA^T=A^TA\).

The `is_positive_semidefinite_2x2` test uses the exact rational criterion for
the **symmetric** two-by-two matrix

\[
M=\begin{pmatrix}a&b\\b&d\end{pmatrix}:
\quad a\ge0,\ d\ge0,\ ad-b^2\ge0.
\]

Larger matrices and nonsymmetric matrices are rejected by that predicate.
`quadratic_form(A, x)` calculates \(x^TAx\) exactly, but evaluates a particular
vector; it is *not* a general positive-definiteness certificate.

## Executable examples

```python
from fractions import Fraction
from calc.operator_algebra import (
    adjoint, product, trace, is_orthogonal_projection,
    is_positive_semidefinite_2x2, star_algebra,
)

P = ((Fraction(1, 2), Fraction(1, 2)),
     (Fraction(1, 2), Fraction(1, 2)))
assert adjoint(P) == P
assert product(P, P) == P
assert is_orthogonal_projection(P)
assert is_positive_semidefinite_2x2(P)
assert trace(P) == 1

E11 = ((1, 0), (0, 0))
E22 = ((0, 0), (0, 1))
diagonal = star_algebra((E11, E22))
assert diagonal.unit == (1, 1)
assert diagonal.multiplication[0][1] == (0, 0)
assert diagonal.involution[0] == (1, 0)
```

The `star_algebra` validator takes an **independent rational basis** of
matrices. It verifies that every basis product \(E_iE_j\) and every adjoint
\(E_i^T\) belongs to its rational span. It then returns exact structure
coefficients in the *input basis*, with the convention

\[
E_iE_j=\sum_k m_{ijk}E_k,\qquad E_i^T=\sum_k s_{ik}E_k.
\]

By default, the identity matrix must lie in the span (`require_unit=True`).
Use `require_unit=False` to validate a nonunital matrix \(*\)-algebra; its
`unit` is `None` when the ambient identity is absent. A span failing closure
is rejected, never silently projected onto the basis. A finite rational
\(*\)-subalgebra is not thereby asserted to be norm-complete as an arbitrary
operator algebra; the scope here is purely algebraic matrix operations.

## Resource and exactness boundaries

- Matrices must be square, at most 8 by 8, with `int` or `Fraction` entries.
- The supplied independent algebra basis contains at most 16 matrices.
- Rational coefficient numerators/denominators are bounded to 2,048 bits;
  outputs exceeding that bound are rejected.
- Float, complex, bool, ragged matrices, and mixed-size operators are refused.
- The library has no approximate equality, eigenvalue certification, or
  inference of continuum operator properties.

Run focused tests with:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_operator_algebra.py -v
```

Tests use exact rational witnesses for the involution anti-automorphism,
cyclic trace, projections, positive-semidefinite boundaries, diagonal and
full matrix algebras, invalid nonclosed spans, and resource limits.

## Sources and reuse

- [MIT 18.712, Introduction to Representation Theory — lecture notes](https://www.ocw.mit.edu/courses/18-712-introduction-to-representation-theory-fall-2010/pages/lecture-notes/): matrix-algebra representations and finite-dimensional algebra structure.
- [Gilbert Strang, MIT Linear Algebra lecture notes](https://math.mit.edu/~gs/LectureNotes/): matrices, transpose, orthogonality, and projections.
- [Encyclopedia of Mathematics — Involution algebra](https://encyclopediaofmath.org/wiki/Involution_algebra): algebraic involution concepts.

These are **linked mathematical references**, not code dependencies. No
third-party source code or document text was copied into this implementation.
