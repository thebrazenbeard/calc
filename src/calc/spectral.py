"""Numerical spectral tools for finite, real, symmetric matrices.

The pure-Python Jacobi rotations here are educational and appropriate for
small dense matrices.  This module does not solve non-normal eigenproblems,
perform a general singular-value decomposition, or certify exact eigenpairs.
"""

from __future__ import annotations

import math
from fractions import Fraction
from numbers import Integral, Real


class JacobiConvergenceError(ArithmeticError):
    """A Jacobi eigensolver exhausted its prescribed sweep budget."""


def _finite_real(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real number")
    try:
        number = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be finite and representable") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _matrix(matrix: object) -> list[list[float]]:
    """Copy and validate an exactly symmetric, nonempty square matrix."""
    try:
        rows = [list(row) for row in matrix]
    except (TypeError, ValueError) as exc:
        raise TypeError("matrix must be a sequence of rows") from exc
    n = len(rows)
    if n == 0 or any(len(row) != n for row in rows):
        raise ValueError("matrix must be nonempty and square")
    data = [[_finite_real(entry, "matrix entry") for entry in row]
            for row in rows]
    for i in range(n):
        for j in range(i + 1, n):
            if data[i][j] != data[j][i]:
                raise ValueError("matrix must be symmetric")
    return data


def _vector(vector: object, n: int, name: str) -> list[float]:
    try:
        values = list(vector)
    except TypeError as exc:
        raise TypeError(f"{name} must be a sequence") from exc
    if len(values) != n:
        raise ValueError(f"{name} must have length {n}")
    return [_finite_real(item, f"{name} entry") for item in values]


def _eigenvectors(vectors: object, n: int) -> list[list[float]]:
    try:
        rows = [list(row) for row in vectors]
    except TypeError as exc:
        raise TypeError("eigenvectors must be a square matrix") from exc
    if len(rows) != n or any(len(row) != n for row in rows):
        raise ValueError(f"eigenvectors must be a {n} by {n} matrix")
    return [[_finite_real(item, "eigenvector entry") for item in row]
            for row in rows]


def _identity(n: int) -> list[list[float]]:
    return [[float(i == j) for j in range(n)] for i in range(n)]


def jacobi_eigen(
    matrix: object, *, tolerance: float = 1e-12, max_sweeps: int = 50,
) -> tuple[tuple[float, ...], tuple[tuple[float, ...], ...]]:
    """Return ascending eigenvalues and corresponding eigenvector columns.

    The input must be exactly real symmetric and finite.  Jacobi sweeps
    diagonalize a copy of the matrix.  ``tolerance`` controls the largest
    remaining off-diagonal entry *relative to the largest input entry*;
    ``max_sweeps`` limits complete passes over all (p, q) index pairs.

    Returns ``(values, vectors)`` where ``vectors[i][j]`` is component *i*
    of eigenvector *j*.  Repeated eigenspaces have non-unique bases.
    """
    tol = _finite_real(tolerance, "tolerance")
    if not 0.0 < tol < 1.0:
        raise ValueError("tolerance must lie strictly between 0 and 1")
    if isinstance(max_sweeps, bool) or not isinstance(max_sweeps, Integral):
        raise TypeError("max_sweeps must be a positive integer")
    if max_sweeps <= 0:
        raise ValueError("max_sweeps must be positive")

    a = _matrix(matrix)
    n = len(a)
    # A matrix which is already diagonal needs no rotations or global
    # scaling: scaling 1e-300 by a 1e300 maximum would erase a real
    # nonzero eigenvalue before the method even starts.
    if all(a[i][j] == 0.0 for i in range(n) for j in range(n) if i != j):
        order = sorted(range(n), key=lambda j: a[j][j])
        return (tuple(a[j][j] for j in order),
                tuple(tuple(float(i == j) for j in order) for i in range(n)))
    scale = max(abs(entry) for row in a for entry in row)
    vectors = _identity(n)
    if scale == 0.0:
        return (0.0,) * n, tuple(tuple(row) for row in vectors)

    a = [[entry / scale for entry in row] for row in a]

    def off_diagonal_max() -> float:
        return max((abs(a[i][j]) for i in range(n)
                    for j in range(i + 1, n)), default=0.0)

    for sweep in range(max_sweeps + 1):
        if off_diagonal_max() <= tol:
            break
        if sweep == max_sweeps:
            raise JacobiConvergenceError(
                f"Jacobi method did not converge within {max_sweeps} sweeps"
            )
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if apq == 0.0:
                    continue
                tau = (a[q][q] - a[p][p]) / (2.0 * apq)
                t = math.copysign(1.0, tau) / (abs(tau) + math.hypot(1.0, tau))
                c = 1.0 / math.hypot(1.0, t)
                s = t * c
                a[p][p] -= t * apq
                a[q][q] += t * apq
                a[p][q] = a[q][p] = 0.0
                for k in range(n):
                    if k != p and k != q:
                        akp, akq = a[k][p], a[k][q]
                        a[k][p] = a[p][k] = c * akp - s * akq
                        a[k][q] = a[q][k] = s * akp + c * akq
                    vkp, vkq = vectors[k][p], vectors[k][q]
                    vectors[k][p] = c * vkp - s * vkq
                    vectors[k][q] = s * vkp + c * vkq

    order = sorted(range(n), key=lambda j: a[j][j])
    eigenvalues = tuple(a[j][j] * scale for j in order)
    if not all(math.isfinite(value) for value in eigenvalues):
        raise ArithmeticError("eigenvalue cannot be represented as a finite float")
    eigenvectors = tuple(tuple(vectors[i][j] for j in order) for i in range(n))
    return eigenvalues, eigenvectors


def rayleigh_quotient(matrix: object, vector: object) -> float:
    """Return ``v.T @ A @ v / (v.T @ v)`` for a nonzero real vector.

    Both the matrix and vector are scaled before arithmetic to reduce
    avoidable overflow and underflow.  The result is a floating estimate.
    """
    a = _matrix(matrix)
    original_v = _vector(vector, len(a), "vector")
    vscale = max(abs(x) for x in original_v)
    if vscale == 0.0:
        raise ValueError("vector must not be zero")
    v = [x / vscale for x in original_v]
    ascale = max(abs(x) for row in a for x in row)
    if ascale == 0.0:
        return 0.0
    numerator = sum(v[i] * sum((a[i][j] / ascale) * v[j]
                               for j in range(len(a))) for i in range(len(a)))
    denominator = sum(x * x for x in v)
    result = (numerator / denominator) * ascale
    if result == 0.0:
        # Relative scaling can lose a tiny, but representable, quotient.
        # Reconstruct exact binary64 inputs on this exceptional path.
        rational_v = [Fraction.from_float(x) for x in original_v]
        exact_numerator = sum(
            (Fraction.from_float(a[i][j]) * rational_v[i] * rational_v[j]
             for i in range(len(a)) for j in range(len(a))), Fraction(0))
        exact_denominator = sum((value * value for value in rational_v), Fraction(0))
        try:
            result = float(exact_numerator / exact_denominator)
        except OverflowError as exc:
            raise ArithmeticError("Rayleigh quotient is not finite") from exc
    if not math.isfinite(result):
        raise ArithmeticError("Rayleigh quotient is not finite")
    return result


def eigen_residuals(
    matrix: object, eigenvalues: object, eigenvectors: object,
) -> tuple[float, ...]:
    """Return scaled max-entry residuals of proposed eigenvector columns.

    Residual j is ``||A v_j - lambda_j v_j||_inf`` divided by
    ``max(1, ||A||_max, |lambda_j|) * ||v_j||_inf``.  It is a numerical
    diagnostic, not proof that the eigenvalues are exact or exhaustive.
    """
    a = _matrix(matrix)
    n = len(a)
    lambdas = _vector(eigenvalues, n, "eigenvalues")
    v = _eigenvectors(eigenvectors, n)
    matrix_scale = max(abs(x) for row in a for x in row)
    residuals = []
    for j in range(n):
        vec_scale = max(abs(v[i][j]) for i in range(n))
        if vec_scale == 0.0:
            raise ValueError("eigenvector columns must be nonzero")
        denom = max(1.0, matrix_scale, abs(lambdas[j]))
        vec = [v[i][j] / vec_scale for i in range(n)]
        error = max(abs(sum((a[i][k] / denom) * vec[k]
                            for k in range(n)) - (lambdas[j] / denom) * vec[i])
                    for i in range(n))
        residuals.append(error)
    return tuple(residuals)


def orthogonality_error(eigenvectors: object) -> float:
    """Return the max-entry deviation of a square matrix's ``V.T @ V`` from I.

    Eigenvectors are columns.  A well-conditioned orthonormal basis gives
    values close to zero; this is a diagnostic, not exact certification.
    """
    try:
        rows = list(eigenvectors)
    except TypeError as exc:
        raise TypeError("eigenvectors must be a matrix") from exc
    n = len(rows)
    if n == 0:
        raise ValueError("eigenvectors must be nonempty")
    v = _eigenvectors(rows, n)
    error = max(abs(sum(v[k][i] * v[k][j] for k in range(n))
                    - float(i == j)) for i in range(n) for j in range(n))
    return error


__all__ = [
    "JacobiConvergenceError", "jacobi_eigen", "rayleigh_quotient",
    "eigen_residuals", "orthogonality_error",
]
