"""Exact, bounded *matrix* Lie algebra computations over the rationals.

The Lie bracket is ``[A, B] = AB - BA``. Every value is a ``Fraction``;
floating inputs are rejected instead of being silently approximated. This
module deals only with matrix representations and their finite rational spans,
not arbitrary Lie algebras, exponentials or classification problems.

Structure constants use the convention ``[e_i, e_j] = sum_k c[i][j][k] e_k``.
Adjoint matrices are column-oriented: column j represents ``[X, e_j]``.
"""

from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction

Rational = int | Fraction
Matrix = tuple[tuple[Fraction, ...], ...]
Basis = tuple[Matrix, ...]
StructureConstants = tuple[tuple[tuple[Fraction, ...], ...], ...]

_MAX_DIMENSION = 8
_MAX_BASIS = 16
_MAX_SCALAR_BITS = 2048


def _rational(value: Rational) -> Fraction:
    if type(value) is int or isinstance(value, Fraction):
        result = Fraction(value)
        if (result.numerator.bit_length() > _MAX_SCALAR_BITS or
                result.denominator.bit_length() > _MAX_SCALAR_BITS):
            raise ValueError("rational coefficient exceeds bit-size limit")
        return result
    raise TypeError("matrix entries and coefficients must be int or Fraction")


def as_matrix(values: Sequence[Sequence[Rational]]) -> Matrix:
    """Normalize a nonempty square rational matrix into immutable fractions.

    Dimensions larger than eight are intentionally refused to bound the exact
    linear-algebra work required by basis and structure-constant operations.
    """
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise TypeError("matrix must be a sequence of rows")
    n = len(values)
    if not 1 <= n <= _MAX_DIMENSION:
        raise ValueError("matrix dimension must be between 1 and 8")
    rows: list[tuple[Fraction, ...]] = []
    for row in values:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence):
            raise TypeError("each matrix row must be a sequence")
        if len(row) != n:
            raise ValueError("matrix must be square")
        rows.append(tuple(_rational(value) for value in row))
    return tuple(rows)


def _same_size(left: Matrix, right: Matrix) -> None:
    if len(left) != len(right):
        raise ValueError("matrices must have equal dimensions")


def _add(left: Matrix, right: Matrix, factor: Fraction = Fraction(1)) -> Matrix:
    _same_size(left, right)
    return tuple(tuple(a + factor * b for a, b in zip(row_a, row_b))
                 for row_a, row_b in zip(left, right))


def matrix_product(left: Sequence[Sequence[Rational]],
                   right: Sequence[Sequence[Rational]]) -> Matrix:
    """Exact matrix multiplication; square, equally dimensioned matrices."""
    a, b = as_matrix(left), as_matrix(right)
    _same_size(a, b)
    n = len(a)
    return tuple(tuple(sum((a[i][k] * b[k][j] for k in range(n)), Fraction(0))
                       for j in range(n)) for i in range(n))


def commutator(left: Sequence[Sequence[Rational]],
               right: Sequence[Sequence[Rational]]) -> Matrix:
    """Lie bracket ``[A, B] = AB - BA`` for rational square matrices."""
    a, b = as_matrix(left), as_matrix(right)
    _same_size(a, b)
    return _add(matrix_product(a, b), matrix_product(b, a), Fraction(-1))


def jacobiator(a: Sequence[Sequence[Rational]],
              b: Sequence[Sequence[Rational]],
              c: Sequence[Sequence[Rational]]) -> Matrix:
    """Compute ``[[A,B],C] + [[B,C],A] + [[C,A],B]`` exactly.

    This is identically zero for the commutator bracket, and is useful as a
    computational witness; it is not a theorem prover.
    """
    aa, bb, cc = as_matrix(a), as_matrix(b), as_matrix(c)
    _same_size(aa, bb)
    _same_size(aa, cc)
    return _add(_add(commutator(commutator(aa, bb), cc),
                     commutator(commutator(bb, cc), aa)),
                commutator(commutator(cc, aa), bb))


def _as_basis(values: Sequence[Sequence[Sequence[Rational]]]) -> Basis:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise TypeError("basis must be a sequence of matrices")
    if not 1 <= len(values) <= _MAX_BASIS:
        raise ValueError("basis must contain between 1 and 16 elements")
    basis = tuple(as_matrix(value) for value in values)
    n = len(basis[0])
    if any(len(matrix) != n for matrix in basis):
        raise ValueError("basis matrices must have a common dimension")
    if len(basis) > n * n:
        raise ValueError("basis cannot be independent in this matrix space")
    return basis


def _coordinates(target: Matrix, basis: Basis) -> tuple[Fraction, ...]:
    """Solve an overdetermined rational system with independent columns.

    Pivots are scanned by variable (basis column); all unused equations are
    checked afterwards, so a target outside the span is never projected in.
    """
    n = len(target)
    _same_size(target, basis[0])
    k = len(basis)
    table = [
        [*(basis[col][i][j] for col in range(k)), target[i][j]]
        for i in range(n) for j in range(n)
    ]
    pivots: list[int] = []
    pivot_row = 0
    for col in range(k):
        pivot = next((row for row in range(pivot_row, len(table))
                      if table[row][col]), None)
        if pivot is None:
            raise ValueError("basis matrices are linearly dependent")
        table[pivot_row], table[pivot] = table[pivot], table[pivot_row]
        scale = table[pivot_row][col]
        table[pivot_row] = [value / scale for value in table[pivot_row]]
        for row in range(len(table)):
            if row == pivot_row:
                continue
            factor = table[row][col]
            if factor:
                table[row] = [value - factor * basis_value
                              for value, basis_value in zip(table[row], table[pivot_row])]
        pivots.append(pivot_row)
        pivot_row += 1
    if any(row[-1] for row in table[pivot_row:]):
        raise ValueError("matrix is not in the rational span of the basis")
    return tuple(table[row][-1] for row in pivots)


def basis_coordinates(target: Sequence[Sequence[Rational]],
                      basis: Sequence[Sequence[Sequence[Rational]]]
                      ) -> tuple[Fraction, ...]:
    """Return unique rational coordinates of a matrix in an independent basis.

    Raises ``ValueError`` if the basis is dependent or the target is outside
    the span; no least-squares approximation or floating solve is used.
    """
    e = _as_basis(basis)
    return _coordinates(as_matrix(target), e)


def linear_combination(coefficients: Sequence[Rational],
                       basis: Sequence[Sequence[Sequence[Rational]]]) -> Matrix:
    """Construct ``sum_i coefficients[i] * basis[i]`` exactly."""
    e = _as_basis(basis)
    if isinstance(coefficients, (str, bytes)) or not isinstance(coefficients, Sequence):
        raise TypeError("coefficients must be a sequence")
    if len(coefficients) != len(e):
        raise ValueError("coefficient count must match the basis")
    cs = tuple(_rational(value) for value in coefficients)
    n = len(e[0])
    return tuple(tuple(sum((cs[k] * e[k][i][j] for k in range(len(e))),
                           Fraction(0)) for j in range(n)) for i in range(n))


def structure_constants(basis: Sequence[Sequence[Sequence[Rational]]]
                        ) -> StructureConstants:
    """Return exact ``c[i][j][k]`` for a rational matrix Lie subalgebra.

    Besides independence, every pairwise bracket must remain in the supplied
    span; a nonclosed span is rejected rather than misrepresented as a Lie
    subalgebra. Complexity is bounded by basis and matrix dimension limits.
    """
    e = _as_basis(basis)
    _coordinates(e[0], e)  # Validate linear independence even for abelian cases.
    return tuple(tuple(_coordinates(commutator(a, b), e) for b in e) for a in e)


def adjoint_matrix(element: Sequence[Sequence[Rational]],
                   basis: Sequence[Sequence[Sequence[Rational]]]) -> Matrix:
    """Matrix of ``ad_X(Y) = [X,Y]`` in a rational Lie subalgebra basis.

    Column ``j`` contains the basis coordinates of ``[X,e_j]``. Requires X
    itself in the subalgebra, and verifies closure of the whole basis.
    """
    e = _as_basis(basis)
    x = as_matrix(element)
    _coordinates(x, e)
    structure_constants(e)
    columns = tuple(_coordinates(commutator(x, element_j), e) for element_j in e)
    return tuple(tuple(columns[j][i] for j in range(len(e))) for i in range(len(e)))


def sl2_basis() -> Basis:
    """Return the conventional ``(H,E,F)`` traceless 2x2 matrix basis.

    ``[H,E]=2E``, ``[H,F]=-2F``, and ``[E,F]=H``.
    """
    return tuple(as_matrix(matrix) for matrix in (
        ((1, 0), (0, -1)),
        ((0, 1), (0, 0)),
        ((0, 0), (1, 0)),
    ))


def so3_basis() -> Basis:
    """Return standard ``(Lx,Ly,Lz)`` 3x3 real skew matrices.

    With this orientation ``[Lx,Ly]=Lz`` and cyclic permutations hold.
    """
    return tuple(as_matrix(matrix) for matrix in (
        ((0, 0, 0), (0, 0, -1), (0, 1, 0)),
        ((0, 0, 1), (0, 0, 0), (-1, 0, 0)),
        ((0, -1, 0), (1, 0, 0), (0, 0, 0)),
    ))


__all__ = (
    "Matrix", "Basis", "StructureConstants", "as_matrix", "matrix_product",
    "commutator", "jacobiator", "basis_coordinates", "linear_combination",
    "structure_constants", "adjoint_matrix", "sl2_basis", "so3_basis",
)
