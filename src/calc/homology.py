"""Small, exact chain complexes of rational vector spaces.

For each n, ``d_n: C_n -> C_(n-1)`` is a matrix with columns indexed by
the basis of C_n.  The constructor verifies ``d_(n-1) @ d_n = 0``.
Only finite, supplied chain groups are supported; no topology or cell
complex is inferred from geometric input.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction


Matrix = tuple[tuple[Fraction, ...], ...]
_MAX_DEGREES = 8
_MAX_DIMENSION = 32
_MAX_RATIONAL_BITS = 2048


def _coefficient(value: int | Fraction) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("chain coefficients must be integers or rational Fractions")
    result = Fraction(value)
    if (result.numerator.bit_length() > _MAX_RATIONAL_BITS or
            result.denominator.bit_length() > _MAX_RATIONAL_BITS):
        raise ValueError("chain coefficient exceeds rational bit limit")
    return result


def _normalize_matrix(
    rows: Sequence[Sequence[int | Fraction]], count_rows: int,
    count_columns: int,
) -> Matrix:
    if isinstance(rows, (str, bytes)) or not isinstance(rows, Sequence):
        raise TypeError("boundary must be a sequence of numeric rows")
    if any(isinstance(row, (str, bytes)) or not isinstance(row, Sequence) for row in rows):
        raise TypeError("boundary rows must be numeric sequences")
    if len(rows) != count_rows or any(len(row) != count_columns for row in rows):
        raise ValueError(f"boundary must have shape {count_rows} x {count_columns}")
    supplied = tuple(tuple(row) for row in rows)
    return tuple(tuple(_coefficient(value) for value in row) for row in supplied)


def rational_rank(rows: Sequence[Sequence[int | Fraction]]) -> int:
    """Rank of a non-ragged rational matrix by exact row elimination.

    An empty matrix has rank zero. Each dimension is bounded to avoid
    accidental unbounded rational arithmetic.
    """
    if isinstance(rows, (str, bytes)) or not isinstance(rows, Sequence):
        raise TypeError("matrix must be a sequence of numeric rows")
    if len(rows) > _MAX_DIMENSION:
        raise ValueError("matrix row limit exceeded")
    if any(isinstance(row, (str, bytes)) or not isinstance(row, Sequence)
           for row in rows):
        raise TypeError("matrix rows must be numeric sequences")
    supplied = tuple(tuple(row) for row in rows)
    if not supplied:
        return 0
    columns = len(supplied[0])
    if columns > _MAX_DIMENSION or any(len(row) != columns for row in supplied):
        raise ValueError("matrix must be rectangular and within dimension limit")
    matrix = [[_coefficient(value) for value in row] for row in supplied]
    pivot_row = 0
    for column in range(columns):
        pivot = next((i for i in range(pivot_row, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        divisor = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / divisor for entry in matrix[pivot_row]]
        for i in range(pivot_row + 1, len(matrix)):
            factor = matrix[i][column]
            if factor:
                matrix[i] = [entry - factor * p for entry, p in zip(matrix[i], matrix[pivot_row])]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


@dataclass(frozen=True, slots=True, init=False)
class ChainComplex:
    """A validated bounded finite chain complex over the field Q.

    ``dimensions[n]`` is dim(C_n), with C_0 at position zero.
    ``boundaries[n]`` supplies the matrix of d_n for n >= 1; omitted maps
    are the correctly shaped zero map. All coefficients must be exact.
    
    Example: an oriented interval has dimensions (2, 1) and d_1=(-1, 1)^T.
    """

    dimensions: tuple[int, ...]
    boundary_matrices: tuple[Matrix, ...]

    def __init__(
        self,
        dimensions: Sequence[int],
        boundaries: Mapping[int, Sequence[Sequence[int | Fraction]]] | None = None,
    ) -> None:
        if isinstance(dimensions, (str, bytes)) or not isinstance(dimensions, Sequence):
            raise TypeError("dimensions must be a sequence of nonnegative integers")
        if not 1 <= len(dimensions) <= _MAX_DEGREES:
            raise ValueError(f"chain degrees must number between 1 and {_MAX_DEGREES}")
        sizes = tuple(dimensions)
        for size in sizes:
            if isinstance(size, bool) or not isinstance(size, int):
                raise TypeError("chain dimensions must be integers, not bool")
            if not 0 <= size <= _MAX_DIMENSION:
                raise ValueError(f"chain dimension must lie in 0..{_MAX_DIMENSION}")
        if boundaries is None:
            boundaries = {}
        if not isinstance(boundaries, Mapping):
            raise TypeError("boundaries must map degree to a matrix")
        if len(boundaries) >= _MAX_DEGREES:
            raise ValueError("too many supplied boundary matrices")
        for degree in boundaries:
            if isinstance(degree, bool) or not isinstance(degree, int):
                raise TypeError("boundary degree must be an integer, not bool")
            if not 1 <= degree < len(sizes):
                raise ValueError("boundary degree is outside the chain complex")
        matrices: list[Matrix] = [()]
        for degree in range(1, len(sizes)):
            rows, cols = sizes[degree - 1], sizes[degree]
            if degree in boundaries:
                matrix = _normalize_matrix(boundaries[degree], rows, cols)
            else:
                matrix = tuple(tuple(Fraction(0) for _ in range(cols)) for _ in range(rows))
            matrices.append(matrix)
        for degree in range(2, len(sizes)):
            earlier = matrices[degree - 1]
            later = matrices[degree]
            # The composition has shape dim(C_(degree-2)) x dim(C_degree).
            for i in range(sizes[degree - 2]):
                for j in range(sizes[degree]):
                    if sum((earlier[i][k] * later[k][j] for k in range(sizes[degree - 1])),
                           Fraction(0)):
                        raise ValueError(f"chain identity d_{degree-1} d_{degree} = 0 violated")
        object.__setattr__(self, "dimensions", sizes)
        object.__setattr__(self, "boundary_matrices", tuple(matrices))

    def boundary(self, degree: int) -> Matrix:
        """Return immutable d_n matrix (n>=1); no d_0 matrix is stored."""
        if isinstance(degree, bool) or not isinstance(degree, int):
            raise TypeError("degree must be an integer")
        if not 1 <= degree < len(self.dimensions):
            raise ValueError("boundary degree out of range")
        return self.boundary_matrices[degree]

    def boundary_rank(self, degree: int) -> int:
        """Rank of d_n; d_0 and d_(top+1) are zero by convention."""
        if isinstance(degree, bool) or not isinstance(degree, int):
            raise TypeError("degree must be an integer")
        if degree == 0 or degree == len(self.dimensions):
            return 0
        return rational_rank(self.boundary(degree))

    def betti_numbers(self) -> tuple[int, ...]:
        """Dimensions of H_n(C;Q) = ker d_n / im d_(n+1)."""
        ranks = [0] + [rational_rank(mat) for mat in self.boundary_matrices[1:]] + [0]
        return tuple(size - ranks[n] - ranks[n+1] for n, size in enumerate(self.dimensions))

    def euler_characteristic(self) -> int:
        """Alternating chain rank, equal to alternating rational Betti rank."""
        return sum(((-1) ** n) * size for n, size in enumerate(self.dimensions))

    def homology_euler_characteristic(self) -> int:
        """Alternating sum of Betti numbers, as an independent identity."""
        return sum(((-1) ** n) * betti for n, betti in enumerate(self.betti_numbers()))


__all__ = ["rational_rank", "ChainComplex"]
