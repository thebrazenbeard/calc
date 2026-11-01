"""Bounded exact-rational tensor *components* over finite index sets.

This module manipulates dense arrays that represent tensors in a fixed basis.
Contractions identify one covariant and one contravariant index *only when the
caller explicitly selects the two axes*. Variance is not tracked here, so
the caller is responsible for mathematically meaningful index placement.
No metric, connection, coordinate transformation, or differential geometry is
implicitly supplied by an array of numbers.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import islice, product
from collections.abc import Iterable

_MAX_RANK = 6
_MAX_AXIS = 32
_MAX_ELEMENTS = 4096
_MAX_COEFFICIENT_BITS = 4096


def _exact(value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("tensor coefficients must be exact ints or Fractions")
    result = Fraction(value)
    if (result.numerator.bit_length() > _MAX_COEFFICIENT_BITS or
            result.denominator.bit_length() > _MAX_COEFFICIENT_BITS):
        raise ValueError("tensor coefficient exceeds the rational size limit")
    return result


def _shape(shape: Iterable[int]) -> tuple[tuple[int, ...], int]:
    if isinstance(shape, (str, bytes)):
        raise TypeError("shape must be an iterable of positive integer dimensions")
    try:
        dims = tuple(islice(shape, _MAX_RANK + 1))
    except TypeError as exc:
        raise TypeError("shape must be iterable") from exc
    if len(dims) > _MAX_RANK:
        raise ValueError(f"tensor rank cannot exceed {_MAX_RANK}")
    count = 1
    for axis in dims:
        if isinstance(axis, bool) or not isinstance(axis, int):
            raise TypeError("shape entries must be integers (not bool)")
        if not 1 <= axis <= _MAX_AXIS:
            raise ValueError(f"axis dimension must be between 1 and {_MAX_AXIS}")
        count *= axis
        if count > _MAX_ELEMENTS:
            raise ValueError(f"tensor cannot exceed {_MAX_ELEMENTS} elements")
    return dims, count


def _indices(shape: tuple[int, ...]):
    return product(*(range(dimension) for dimension in shape))


@dataclass(frozen=True, slots=True, init=False)
class Tensor:
    """Immutable row-major dense tensor with exact rational components.

    ``Tensor((2, 2), [1, 2, 3, 4])`` has entries ``T[0, 0] == 1`` and
    ``T[1, 0] == 3``. Shape ``()`` is a rank-zero scalar with one component.
    """

    shape: tuple[int, ...]
    values: tuple[Fraction, ...]

    def __init__(self, shape: Iterable[int], values: Iterable[int | Fraction]):
        dims, count = _shape(shape)
        if isinstance(values, (str, bytes)):
            raise TypeError("values must be a finite iterable of exact rationals")
        try:
            raw = tuple(islice(values, count + 1))
        except TypeError as exc:
            raise TypeError("values must be iterable") from exc
        if len(raw) != count:
            raise ValueError(f"shape {dims} requires exactly {count} values")
        object.__setattr__(self, "shape", dims)
        object.__setattr__(self, "values", tuple(_exact(value) for value in raw))

    @property
    def rank(self) -> int:
        return len(self.shape)

    @staticmethod
    def scalar(value: int | Fraction) -> Tensor:
        return Tensor((), (value,))

    @staticmethod
    def zeros(shape: Iterable[int]) -> Tensor:
        dims, count = _shape(shape)
        return Tensor(dims, (0 for _ in range(count)))

    def __getitem__(self, indices: int | tuple[int, ...]) -> Fraction:
        coordinates = (indices,) if isinstance(indices, int) and not isinstance(indices, bool) else indices
        if not isinstance(coordinates, tuple) or len(coordinates) != self.rank:
            raise IndexError("index must be a tuple containing one index per axis")
        offset = 0
        for index, dim in zip(coordinates, self.shape):
            if isinstance(index, bool) or not isinstance(index, int):
                raise TypeError("tensor indices must be integers, not bool")
            if index < 0 or index >= dim:
                raise IndexError("tensor index is out of range")
            offset = offset * dim + index
        return self.values[offset]

    def __add__(self, other: Tensor) -> Tensor:
        if not isinstance(other, Tensor):
            return NotImplemented
        if self.shape != other.shape:
            raise ValueError("tensor addition requires identical shape")
        return Tensor(self.shape, (left + right for left, right in zip(self.values, other.values)))

    def __neg__(self) -> Tensor:
        return Tensor(self.shape, (-value for value in self.values))

    def __sub__(self, other: Tensor) -> Tensor:
        if not isinstance(other, Tensor):
            return NotImplemented
        return self + (-other)

    def __mul__(self, scalar: int | Fraction) -> Tensor:
        coefficient = _exact(scalar)
        return Tensor(self.shape, (coefficient * value for value in self.values))

    __rmul__ = __mul__

    def __truediv__(self, scalar: int | Fraction) -> Tensor:
        coefficient = _exact(scalar)
        if coefficient == 0:
            raise ZeroDivisionError("cannot divide tensor by zero")
        return Tensor(self.shape, (value / coefficient for value in self.values))

    def reshape(self, new_shape: Iterable[int]) -> Tensor:
        dims, count = _shape(new_shape)
        if count != len(self.values):
            raise ValueError("reshape cannot change the component count")
        return Tensor(dims, self.values)

    def permute_axes(self, ordering: Iterable[int]) -> Tensor:
        try:
            order = tuple(islice(ordering, _MAX_RANK + 1))
        except TypeError as exc:
            raise TypeError("axis ordering must be iterable") from exc
        if (len(order) != self.rank or
                any(isinstance(axis, bool) or not isinstance(axis, int) for axis in order) or
                set(order) != set(range(self.rank))):
            raise ValueError("axis ordering must be a permutation of the tensor axes")
        out_shape = tuple(self.shape[axis] for axis in order)
        result: list[Fraction] = []
        for output_indices in _indices(out_shape):
            original_indices = [0] * self.rank
            for out_axis, original_axis in enumerate(order):
                original_indices[original_axis] = output_indices[out_axis]
            result.append(self[tuple(original_indices)])
        return Tensor(out_shape, result)

    def transpose(self) -> Tensor:
        """Reverse all axes; for rank-two arrays this is matrix transpose."""
        return self.permute_axes(reversed(range(self.rank)))

    def tensor_product(self, other: Tensor) -> Tensor:
        """Outer (Kronecker) tensor product ``A ⊗ B`` in this fixed basis."""
        if not isinstance(other, Tensor):
            raise TypeError("tensor product requires a Tensor")
        shape, _ = _shape(self.shape + other.shape)
        return Tensor(shape, (left * right for left in self.values for right in other.values))

    def contract(self, first_axis: int, second_axis: int) -> Tensor:
        """Explicitly sum over two same-sized axes, lowering rank by two.

        For a matrix ``A``, ``A.contract(0, 1)`` returns its scalar trace.
        This is a *component array operation* and makes no variance assertion.
        """
        if (isinstance(first_axis, bool) or isinstance(second_axis, bool) or
                not isinstance(first_axis, int) or not isinstance(second_axis, int)):
            raise TypeError("contraction axes must be integer positions")
        if first_axis == second_axis or not (0 <= first_axis < self.rank) or not (0 <= second_axis < self.rank):
            raise ValueError("contraction needs two distinct existing axes")
        if self.shape[first_axis] != self.shape[second_axis]:
            raise ValueError("contracted axes must have matching dimensions")
        remaining = tuple(axis for axis in range(self.rank) if axis not in (first_axis, second_axis))
        out_shape = tuple(self.shape[axis] for axis in remaining)
        result: list[Fraction] = []
        for out_indices in _indices(out_shape):
            original = [0] * self.rank
            for position, axis in enumerate(remaining):
                original[axis] = out_indices[position]
            total = Fraction(0)
            for diagonal in range(self.shape[first_axis]):
                original[first_axis] = diagonal
                original[second_axis] = diagonal
                total += self[tuple(original)]
            result.append(total)
        return Tensor(out_shape, result)

    def symmetrize(self, first_axis: int, second_axis: int) -> Tensor:
        """Average a tensor with the swap of two same-sized axes."""
        return self._symmetry_part(first_axis, second_axis, sign=1)

    def antisymmetrize(self, first_axis: int, second_axis: int) -> Tensor:
        """Alternating component part under one transposition."""
        return self._symmetry_part(first_axis, second_axis, sign=-1)

    def _symmetry_part(self, first_axis: int, second_axis: int, sign: int) -> Tensor:
        if (isinstance(first_axis, bool) or isinstance(second_axis, bool) or
                not isinstance(first_axis, int) or not isinstance(second_axis, int)):
            raise TypeError("symmetry axes must be integer positions")
        if first_axis == second_axis or not (0 <= first_axis < self.rank) or not (0 <= second_axis < self.rank):
            raise ValueError("symmetry requires two distinct existing axes")
        if self.shape[first_axis] != self.shape[second_axis]:
            raise ValueError("symmetry axes must have matching dimensions")
        order = list(range(self.rank))
        order[first_axis], order[second_axis] = order[second_axis], order[first_axis]
        swapped = self.permute_axes(order)
        return (self + (sign * swapped)) / 2

    def matrix_multiply(self, other: Tensor) -> Tensor:
        """Ordinary exact matrix multiplication, with no implicit tensor metric."""
        if not isinstance(other, Tensor) or self.rank != 2 or other.rank != 2:
            raise TypeError("matrix multiplication requires two rank-two tensors")
        rows, middle = self.shape
        if middle != other.shape[0]:
            raise ValueError("matrix inner dimensions must match")
        columns = other.shape[1]
        return Tensor((rows, columns),
                      (sum((self[i, k] * other[k, j] for k in range(middle)), Fraction(0))
                       for i in range(rows) for j in range(columns)))

    __matmul__ = matrix_multiply


def identity_tensor(dimension: int) -> Tensor:
    """Rank-two Kronecker delta in a fixed finite basis."""
    _shape((dimension, dimension))
    return Tensor((dimension, dimension),
                  (int(i == j) for i in range(dimension) for j in range(dimension)))


__all__ = ["Tensor", "identity_tensor"]
