"""Exact finite-dimensional real-rational matrix *-algebra examples.

The involution is transpose, because the coefficient field is Q embedded in
the real numbers.  These functions manipulate bounded matrix *representations*;
they do not construct general complex C*-algebras, bounded operators on
infinite-dimensional Hilbert spaces, or an operator norm completion.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

from .lie import Matrix, as_matrix, basis_coordinates, matrix_product

Rational = int | Fraction
_MAX_BASIS = 16
_MAX_BITS = 2048


def _scalar(value: Rational) -> Fraction:
    if type(value) is not int and not isinstance(value, Fraction):
        raise TypeError("scalar must be int or Fraction, not float or bool")
    fraction = Fraction(value)
    if (fraction.numerator.bit_length() > _MAX_BITS or
            fraction.denominator.bit_length() > _MAX_BITS):
        raise ValueError("rational scalar exceeds bit-size bound")
    return fraction


def _bounded_matrix(rows: Sequence[Sequence[Rational]]) -> Matrix:
    return as_matrix(rows)


def identity(dimension: int) -> Matrix:
    """Construct the rational identity matrix, dimensions 1 through 8."""
    if type(dimension) is not int or not 1 <= dimension <= 8:
        raise ValueError("dimension must be an integer from 1 to 8")
    return _bounded_matrix(tuple(tuple(int(i == j) for j in range(dimension))
                                 for i in range(dimension)))


def adjoint(matrix: Sequence[Sequence[Rational]]) -> Matrix:
    """The real adjoint A* = A^T in the standard Euclidean inner product."""
    a = _bounded_matrix(matrix)
    return tuple(tuple(a[j][i] for j in range(len(a))) for i in range(len(a)))


def product(left: Sequence[Sequence[Rational]],
            right: Sequence[Sequence[Rational]]) -> Matrix:
    """Exact rational operator composition AB, with size/precision bounds."""
    return _bounded_matrix(matrix_product(left, right))


def _pair(left: Sequence[Sequence[Rational]],
          right: Sequence[Sequence[Rational]]) -> tuple[Matrix, Matrix]:
    a, b = _bounded_matrix(left), _bounded_matrix(right)
    if len(a) != len(b):
        raise ValueError("operators must act on the same finite-dimensional space")
    return a, b


def commutator(left: Sequence[Sequence[Rational]],
               right: Sequence[Sequence[Rational]]) -> Matrix:
    """[A,B] = AB - BA; an exact rational Lie bracket."""
    a, b = _pair(left, right)
    ab, ba = product(a, b), product(b, a)
    return _bounded_matrix(tuple(tuple(ab[i][j] - ba[i][j] for j in range(len(a)))
                                 for i in range(len(a))))


def anticommutator(left: Sequence[Sequence[Rational]],
                   right: Sequence[Sequence[Rational]]) -> Matrix:
    """{A,B} = AB + BA, unscaled (no implicit factor of 1/2)."""
    a, b = _pair(left, right)
    ab, ba = product(a, b), product(b, a)
    return _bounded_matrix(tuple(tuple(ab[i][j] + ba[i][j] for j in range(len(a)))
                                 for i in range(len(a))))


def trace(matrix: Sequence[Sequence[Rational]]) -> Fraction:
    """Exact matrix trace, without normalizing by dimension."""
    a = _bounded_matrix(matrix)
    return _scalar(sum((a[i][i] for i in range(len(a))), Fraction(0)))


def quadratic_form(matrix: Sequence[Sequence[Rational]],
                   vector: Sequence[Rational]) -> Fraction:
    """Compute x^T A x exactly; A need not be symmetric."""
    a = _bounded_matrix(matrix)
    if isinstance(vector, (str, bytes)) or not isinstance(vector, Sequence):
        raise TypeError("vector must be a finite sequence of rationals")
    if len(vector) != len(a):
        raise ValueError("vector length must equal operator dimension")
    x = tuple(_scalar(v) for v in vector)
    value = sum((x[i] * a[i][j] * x[j]
                 for i in range(len(a)) for j in range(len(a))), Fraction(0))
    return _scalar(value)


def is_self_adjoint(matrix: Sequence[Sequence[Rational]]) -> bool:
    """Determine whether A=A^T (real standard inner product)."""
    a = _bounded_matrix(matrix)
    return a == adjoint(a)


def is_orthogonal_projection(matrix: Sequence[Sequence[Rational]]) -> bool:
    """Check exact self-adjoint idempotence A=A*=A²."""
    a = _bounded_matrix(matrix)
    return a == adjoint(a) and a == product(a, a)


def is_normal(matrix: Sequence[Sequence[Rational]]) -> bool:
    """Check exact normality AA*=A*A over the real rational field."""
    a = _bounded_matrix(matrix)
    star = adjoint(a)
    return product(a, star) == product(star, a)


def is_positive_semidefinite_2x2(matrix: Sequence[Sequence[Rational]]) -> bool:
    """Exact PSD test for a real symmetric 2x2 matrix.

    By Sylvester's criterion for nonnegative principal minors, a,d >= 0 and
    ad-b² >= 0 suffice. Nonsymmetric/other-sized inputs are invalid.
    """
    a = _bounded_matrix(matrix)
    if len(a) != 2 or a[0][1] != a[1][0]:
        raise ValueError("PSD predicate requires a symmetric 2x2 matrix")
    return a[0][0] >= 0 and a[1][1] >= 0 and (a[0][0] * a[1][1] - a[0][1] ** 2) >= 0


@dataclass(frozen=True, slots=True)
class RationalStarAlgebra:
    """Validated structure constants for a small rational matrix *-algebra.

    `multiplication[i][j][k]` is coefficient k of E_i E_j.
    `involution[i][k]` is coefficient k of E_i^T.
    `unit` is the coordinate vector of I, or None if require_unit=False and
    the identity is outside the span.
    """

    basis: tuple[Matrix, ...]
    multiplication: tuple[tuple[tuple[Fraction, ...], ...], ...]
    involution: tuple[tuple[Fraction, ...], ...]
    unit: tuple[Fraction, ...] | None


def star_algebra(
    basis: Sequence[Sequence[Sequence[Rational]]], *, require_unit: bool = True,
) -> RationalStarAlgebra:
    """Validate a rational matrix span closed under product and transpose.

    Input basis elements must be linearly independent, same-sized square
    matrices; the rational vector span must contain every pairwise product
    and transpose. By default the identity is also required, yielding a
    finite-dimensional unital real-rational matrix *-algebra example.

    Failure to close is an explicit ValueError rather than projection into
    the closest span. This checks finitely many basis relations, not arbitrary
    or operator-norm closure in an infinite-dimensional topology.
    """
    if type(require_unit) is not bool:
        raise TypeError("require_unit must be bool")
    if isinstance(basis, (str, bytes)) or not isinstance(basis, Sequence):
        raise TypeError("basis must be a finite sequence of matrices")
    if not 1 <= len(basis) <= _MAX_BASIS:
        raise ValueError("basis must contain between 1 and 16 matrices")
    e = tuple(_bounded_matrix(entry) for entry in basis)
    n = len(e[0])
    if any(len(item) != n for item in e):
        raise ValueError("basis matrices must be equally dimensioned")
    if len(e) > n * n:
        raise ValueError("too many basis elements for this matrix space")
    basis_coordinates(e[0], e)  # Forces linear independence even for abelian spans.
    multiplication = tuple(tuple(basis_coordinates(product(a, b), e) for b in e)
                           for a in e)
    involution = tuple(basis_coordinates(adjoint(a), e) for a in e)
    try:
        unit: tuple[Fraction, ...] | None = basis_coordinates(identity(n), e)
    except ValueError:
        if require_unit:
            raise ValueError("unital star algebra basis span must include identity") from None
        unit = None
    return RationalStarAlgebra(e, multiplication, involution, unit)


__all__ = (
    "RationalStarAlgebra", "identity", "adjoint", "product", "commutator",
    "anticommutator", "trace", "quadratic_form", "is_self_adjoint",
    "is_orthogonal_projection", "is_normal", "is_positive_semidefinite_2x2",
    "star_algebra",
)
