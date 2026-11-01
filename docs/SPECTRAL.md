# Spectral methods: real symmetric matrices

The `calc.spectral` module supplies an educational, pure-Python **Jacobi
eigensolver** for small real symmetric matrices, plus tools to inspect
candidate eigenpairs. It requires no numerical computing dependencies.

It is intentionally scoped: no general nonsymmetric eigensolver, singular-value
decomposition, generalized eigenproblem, or certified interval bound is
implemented. For larger or performance-sensitive numerical tasks, use a
well-tested numerical library such as LAPACK.

## Symmetric spectral theorem

Every real symmetric matrix \(A\in\mathbb R^{n\times n}\) admits a decomposition

\[
A = V\,\mathrm{diag}(\lambda_1,\ldots,\lambda_n)\,V^T,
\qquad V^TV=I.
\]

The eigenvalues are real and the eigenvectors can be selected orthonormal.
For repeated eigenvalues, any orthonormal basis of the relevant eigenspace is
valid; signs and individual columns are not uniquely determined.

The **Jacobi algorithm** repeatedly applies orthogonal plane rotations to
eliminate off-diagonal entries. Each full sweep visits every pair of column
indices. Its output approaches a diagonal matrix whose diagonal entries are
the eigenvalues; the product of rotations supplies the eigenvectors.

## Usage

```python
from calc.spectral import (
    jacobi_eigen,
    rayleigh_quotient,
    eigen_residuals,
    orthogonality_error,
)

a = ((2.0, 1.0), (1.0, 2.0))
values, vectors = jacobi_eigen(a, tolerance=1e-12, max_sweeps=50)

assert all(abs(found - expected) < 1e-10
           for found, expected in zip(values, (1.0, 3.0)))
assert max(eigen_residuals(a, values, vectors)) < 1e-10
assert orthogonality_error(vectors) < 1e-10
assert abs(rayleigh_quotient(a, (1.0, 1.0)) - 3.0) < 1e-10
```

The matrix `vectors` is **row-major**, with eigenvectors in its **columns**:
`vectors[i][j]` is coordinate `i` of eigenvector `j`. Eigenvalues are returned
in ascending order with corresponding eigenvector columns. Inputs are never
modified. Return values are immutable tuples.

### API contracts

| Function | Output | Important preconditions |
| --- | --- | --- |
| `jacobi_eigen(A, *, tolerance=1e-12, max_sweeps=50)` | `(eigenvalues, eigenvectors)` | Nonempty, square, **exactly symmetric**, finite real matrix; `0 < tolerance < 1`; positive integer sweep limit |
| `rayleigh_quotient(A, v)` | `float` | Symmetric finite `A`; finite nonzero `v` of matching dimension |
| `eigen_residuals(A, values, vectors)` | Tuple of nonnegative floats | Same symmetric matrix; one value per nonzero eigenvector column; finite entries |
| `orthogonality_error(vectors)` | Nonnegative float | Nonempty, square, finite real matrix; candidate vectors are columns |

If the solver cannot reduce the maximum off-diagonal entry to at most
`tolerance` times the maximum absolute input matrix entry within `max_sweeps`
complete passes, it raises `JacobiConvergenceError` (a subclass of
`ArithmeticError`). In a `2 x 2` matrix, one rotation is often enough. Invalid
shapes, inputs, or parameters raise `TypeError` or `ValueError`. An eigenvalue
or quotient too large for finite Python floating point raises
`ArithmeticError`.

The input symmetry check is exact: `A[i][j] == A[j][i]`. If you constructed a
theoretically symmetric matrix through numerically different expressions,
symmetrize it before calling these functions.

## Interpreting diagnostics

The `j`th residual is

\[
r_j = \frac{\lVert Av_j-\lambda_j v_j\rVert_\infty}
{\max(1,\max_{i,k}|A_{ik}|,|\lambda_j|)\,\lVert v_j\rVert_\infty}.
\]

This scaling gives a useful *backward-error-style diagnostic*, not a
certification of exact eigenvalues or a guarantee on an individual eigenvector
inside a nearly repeated eigenspace. The orthogonality diagnostic is
\(\lVert V^TV-I\rVert_{\max}\), the largest entrywise absolute deviation.
For a credible decomposition, check both and optionally reconstruct
\(V\operatorname{diag}(\lambda)V^T\).

The **Rayleigh quotient** for nonzero \(v\) is

\[
R_A(v) = \frac{v^TAv}{v^Tv}, \qquad
\lambda_{\min} \le R_A(v) \le \lambda_{\max}
\]

for symmetric \(A\). The extrema occur at the corresponding eigenvectors.
The code scales the matrix and vector before floating-point arithmetic to
reduce avoidable overflow or underflow. It cannot guarantee accuracy when
entries undergo severe cancellation or results lie outside floating-point
range.

## Independently checkable examples

- **Two by two:** \(\begin{pmatrix}2&1\\1&2\end{pmatrix}\) has eigenvalues
  \(1,3\), with eigenvectors proportional to \((1,-1)\) and \((1,1)\).
- **Three by three with a repeated eigenvalue:** \(3I+J\), where \(J\) is the
  \(3\times3\) all-ones matrix, has eigenvalues \(3,3,6\).
- **Three by three with distinct eigenvalues:** the symmetric tridiagonal
  matrix with diagonal entries 2 and adjacent off-diagonal entries 1 has
  eigenvalues \(2-\sqrt 2,2,2+\sqrt 2\).
- **General cross-check:** every symmetric eigenvalue lies in a Gershgorin
  interval defined by diagonal entries and sums of off-diagonal magnitudes;
  also \(\sum_j\lambda_j=\operatorname{trace}(A)\).

These examples, reconstruction, Rayleigh bounds, and malformed-input behavior
are exercised in `tests/test_spectral.py`:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_spectral.py' -v
```

## Further reading

- [MIT OpenCourseWare 18.06: Lecture 25, Symmetric Matrices and Positive
  Definiteness](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/video_galleries/video-lectures/)
- [Netlib LAPACK `DSYEV`: eigenvalues and eigenvectors of real symmetric
  matrices](https://www.netlib.org/lapack/explore-html/d8/d1c/group__heev_ga8995c47a7578fef733189df3490258ff.html)
- [Reference LAPACK implementation on GitHub](https://github.com/Reference-LAPACK/lapack/blob/master/SRC/dsyev.f)

These are background references, not source-code dependencies or a claim of
equivalence with production LAPACK routines.
