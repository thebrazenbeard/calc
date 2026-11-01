"""Exact rational representations of small, explicitly enumerated finite groups.

This module checks algebraic identities from supplied finite tables. It does
not construct general irreducible decompositions, complex representations,
or representations of infinite groups.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from itertools import islice

Rational = int | Fraction
Matrix = tuple[tuple[Fraction, ...], ...]

_MAX_GROUP = 16
_MAX_DIMENSION = 16
_MAX_BITS = 4096


def _labels(values: Iterable[str]) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)):
        raise TypeError("group elements must be an iterable of labels")
    try:
        labels = tuple(islice(iter(values), _MAX_GROUP + 1))
    except TypeError as exc:
        raise TypeError("group elements must be iterable") from exc
    if not 1 <= len(labels) <= _MAX_GROUP:
        raise ValueError("group must have between 1 and 16 elements")
    if any(type(item) is not str or not item for item in labels):
        raise TypeError("group labels must be nonempty strings")
    if len(set(labels)) != len(labels):
        raise ValueError("group elements must be distinct")
    return labels


def _rational(value: Rational) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("matrix coefficients must be exact integers or Fractions")
    exact = Fraction(value)
    if max(exact.numerator.bit_length(), exact.denominator.bit_length()) > _MAX_BITS:
        raise ValueError("rational coefficient exceeds bit budget")
    return exact


def _matrix(values: Sequence[Sequence[Rational]], rows: int | None = None,
            cols: int | None = None) -> Matrix:
    if isinstance(values, (str, bytes)):
        raise TypeError("matrix must be a rectangular rational array")
    try:
        row_iter = iter(values)
        raw_rows = tuple(islice(row_iter, _MAX_DIMENSION + 1))
        raw = tuple(tuple(islice(iter(row), _MAX_DIMENSION + 1)) for row in raw_rows)
    except TypeError as exc:
        raise TypeError("matrix must be a rectangular rational array") from exc
    if not raw or not raw[0] or any(len(row) != len(raw[0]) for row in raw):
        raise ValueError("matrix must be nonempty and rectangular")
    if len(raw) > _MAX_DIMENSION or len(raw[0]) > _MAX_DIMENSION:
        raise ValueError("matrix exceeds dimension budget")
    if rows is not None and len(raw) != rows:
        raise ValueError("matrix has an incorrect row count")
    if cols is not None and len(raw[0]) != cols:
        raise ValueError("matrix has an incorrect column count")
    return tuple(tuple(_rational(value) for value in row) for row in raw)


def _identity(dimension: int) -> Matrix:
    return tuple(tuple(Fraction(i == j) for j in range(dimension))
                 for i in range(dimension))


def _multiply(left: Matrix, right: Matrix) -> Matrix:
    if len(left[0]) != len(right):
        raise ValueError("matrix product dimensions do not match")
    output = tuple(tuple(sum((left[i][k] * right[k][j] for k in range(len(right))),
                             Fraction(0))
                         for j in range(len(right[0]))) for i in range(len(left)))
    for row in output:
        for value in row:
            if max(value.numerator.bit_length(), value.denominator.bit_length()) > _MAX_BITS:
                raise ValueError("matrix product exceeds rational bit budget")
    return output


@dataclass(frozen=True, slots=True, init=False)
class FiniteGroup:
    """A bounded group with an explicit, fully validated multiplication table."""

    elements: tuple[str, ...]
    identity: str
    table: tuple[tuple[str, ...], ...]

    def __init__(self, elements: Iterable[str],
                 multiplication: Mapping[tuple[str, str], str], identity: str):
        names = _labels(elements)
        if type(identity) is not str or identity not in names:
            raise ValueError("identity must name a group element")
        if not isinstance(multiplication, Mapping):
            raise TypeError("multiplication must be a complete mapping")
        if len(multiplication) != len(names) ** 2:
            raise ValueError("multiplication table has the wrong number of entries")
        expected = {(left, right) for left in names for right in names}
        if set(multiplication) != expected:
            raise ValueError("multiplication table must contain exactly all ordered pairs")
        table = tuple(tuple(multiplication[left, right] for right in names)
                      for left in names)
        if any(type(value) is not str or value not in names
               for row in table for value in row):
            raise ValueError("multiplication must close on group elements")
        index = {name: i for i, name in enumerate(names)}
        e = index[identity]
        for a in range(len(names)):
            if table[e][a] != names[a] or table[a][e] != names[a]:
                raise ValueError("identity laws are violated")
            if not any(table[a][b] == identity and table[b][a] == identity
                       for b in range(len(names))):
                raise ValueError("every group element needs a two-sided inverse")
        for a in range(len(names)):
            for b in range(len(names)):
                for c in range(len(names)):
                    lhs = table[index[table[a][b]]][c]
                    rhs = table[a][index[table[b][c]]]
                    if lhs != rhs:
                        raise ValueError("group multiplication is not associative")
        object.__setattr__(self, "elements", names)
        object.__setattr__(self, "identity", identity)
        object.__setattr__(self, "table", table)

    def multiply(self, left: str, right: str) -> str:
        """Return ``left * right``; reject unknown element labels."""
        try:
            return self.table[self.elements.index(left)][self.elements.index(right)]
        except ValueError as exc:
            raise KeyError("unknown group element") from exc

    def inverse(self, element: str) -> str:
        """Return the unique inverse of an element in the validated group."""
        if element not in self.elements:
            raise KeyError(element)
        return next(candidate for candidate in self.elements
                    if self.multiply(element, candidate) == self.identity)


@dataclass(frozen=True, slots=True, init=False)
class RationalRepresentation:
    """Finite group homomorphism into invertible rational square matrices.

    All group laws are checked in exact arithmetic at construction time.
    The matrix data are stored in immutable group-element order.
    """

    group: FiniteGroup
    matrices: tuple[Matrix, ...]
    dimension: int

    def __init__(self, group: FiniteGroup,
                 matrices: Mapping[str, Sequence[Sequence[Rational]]]):
        if not isinstance(group, FiniteGroup):
            raise TypeError("group must be a FiniteGroup")
        if (not isinstance(matrices, Mapping) or len(matrices) != len(group.elements)
                or set(matrices) != set(group.elements)):
            raise ValueError("one matrix is required for every group element")
        parsed = tuple(_matrix(matrices[label]) for label in group.elements)
        dimension = len(parsed[0])
        if any(len(matrix) != dimension or len(matrix[0]) != dimension
               for matrix in parsed):
            raise ValueError("representation matrices must be equally-sized square matrices")
        if parsed[group.elements.index(group.identity)] != _identity(dimension):
            raise ValueError("the identity must map to the identity matrix")
        index = {label: i for i, label in enumerate(group.elements)}
        for left in group.elements:
            for right in group.elements:
                if _multiply(parsed[index[left]], parsed[index[right]]) != parsed[index[group.multiply(left, right)]]:
                    raise ValueError("matrices violate the representation homomorphism")
        object.__setattr__(self, "group", group)
        object.__setattr__(self, "matrices", parsed)
        object.__setattr__(self, "dimension", dimension)

    def matrix(self, element: str) -> Matrix:
        try:
            return self.matrices[self.group.elements.index(element)]
        except ValueError as exc:
            raise KeyError(element) from exc

    def character(self) -> dict[str, Fraction]:
        """Return exact traces χ(g) as a fresh mapping."""
        return {name: sum((matrix[i][i] for i in range(self.dimension)), Fraction(0))
                for name, matrix in zip(self.group.elements, self.matrices)}


def cyclic_group(order: int) -> FiniteGroup:
    """Construct C_n with generator g and identity e, for 1 <= n <= 16."""
    if type(order) is not int or not 1 <= order <= _MAX_GROUP:
        raise ValueError("cyclic group order must be in 1..16")
    names = tuple("e" if i == 0 else "g" if i == 1 else f"g{i}"
                  for i in range(order))
    table = {(a, b): names[(i + j) % order]
             for i, a in enumerate(names) for j, b in enumerate(names)}
    return FiniteGroup(names, table, "e")


def trivial_representation(group: FiniteGroup) -> RationalRepresentation:
    """The dimension-one trivial representation for any finite group."""
    if not isinstance(group, FiniteGroup):
        raise TypeError("group must be a FiniteGroup")
    return RationalRepresentation(group, {name: ((1,),) for name in group.elements})


def regular_representation(group: FiniteGroup) -> RationalRepresentation:
    """Left regular action ρ(g)e_h=e_(gh) as permutation matrices."""
    if not isinstance(group, FiniteGroup):
        raise TypeError("group must be a FiniteGroup")
    names = group.elements
    matrices = {}
    for left in names:
        matrix = [[0 for _ in names] for _ in names]
        for column, right in enumerate(names):
            row = names.index(group.multiply(left, right))
            matrix[row][column] = 1
        matrices[left] = matrix
    return RationalRepresentation(group, matrices)


def character_pairing(left: RationalRepresentation,
                      right: RationalRepresentation) -> Fraction:
    """Exact finite-group character pairing sum_g χ1(g)χ2(g^-1)/|G|.

    This is a bilinear character identity, not an irreducibility detector
    for arbitrary rational representations.
    """
    if not isinstance(left, RationalRepresentation) or not isinstance(right, RationalRepresentation):
        raise TypeError("both inputs must be rational representations")
    if left.group != right.group:
        raise ValueError("character pairing needs the same finite group")
    a, b = left.character(), right.character()
    return sum((a[element] * b[left.group.inverse(element)]
                for element in left.group.elements), Fraction(0)) / len(left.group.elements)


def direct_sum(left: RationalRepresentation,
               right: RationalRepresentation) -> RationalRepresentation:
    """Block-diagonal sum; character equals the sum of the characters."""
    if not isinstance(left, RationalRepresentation) or not isinstance(right, RationalRepresentation):
        raise TypeError("both inputs must be rational representations")
    if left.group != right.group:
        raise ValueError("direct sum needs representations of the same group")
    if left.dimension + right.dimension > _MAX_DIMENSION:
        raise ValueError("direct sum exceeds dimension budget")
    matrices = {}
    for element in left.group.elements:
        a, b = left.matrix(element), right.matrix(element)
        matrices[element] = (tuple(tuple(row) + (Fraction(0),) * right.dimension for row in a)
                             + tuple((Fraction(0),) * left.dimension + tuple(row) for row in b))
    return RationalRepresentation(left.group, matrices)


def intertwines(source: RationalRepresentation, target: RationalRepresentation,
                map_matrix: Sequence[Sequence[Rational]]) -> bool:
    """Test exact equivariance target(g)*S = S*source(g) for every g."""
    if not isinstance(source, RationalRepresentation) or not isinstance(target, RationalRepresentation):
        raise TypeError("both inputs must be rational representations")
    if source.group != target.group:
        raise ValueError("intertwining map requires the same group")
    transform = _matrix(map_matrix, rows=target.dimension, cols=source.dimension)
    return all(_multiply(target.matrix(g), transform) == _multiply(transform, source.matrix(g))
               for g in source.group.elements)


__all__ = ["FiniteGroup", "RationalRepresentation", "cyclic_group",
           "trivial_representation", "regular_representation",
           "character_pairing", "direct_sum", "intertwines"]
