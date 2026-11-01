"""Finite free modules over the residue ring ``Z/nZ``.

This is a *finite*, exact, deliberately small model of module theory:
``(Z/nZ)^d`` with componentwise operations.  Submodules are constructed by
additive closure of generators, and quotient modules are actual finite coset
partitions.  No arbitrary modules, infinite-dimensional spaces, Smith normal
form, or PID structure theorem is claimed.

Vectors are canonical tuples of ordinary Python integers in ``range(n)``.
Floats and booleans are rejected instead of being silently rounded/coerced.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable
from dataclasses import dataclass
from itertools import islice, product


Vector = tuple[int, ...]
_MAX_CARDINALITY = 4096
_MAX_RANK = 8
_MAX_GENERATORS = 16
_FACTORY_TOKEN = object()


def _integer(value: object, name: str) -> int:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer (not bool or float)")
    return value


def _finite_coordinates(values: Iterable[int], rank: int) -> tuple[int, ...]:
    if isinstance(values, (str, bytes)):
        raise TypeError("vector coordinates must be integers")
    try:
        result = tuple(islice(iter(values), rank + 1))
    except TypeError as exc:
        raise TypeError("vector must be iterable") from exc
    if len(result) != rank:
        raise ValueError(f"vector requires exactly {rank} coordinates")
    for value in result:
        _integer(value, "coordinate")
    return result


@dataclass(frozen=True, slots=True)
class ResidueFreeModule:
    """The finite free module ``(Z/nZ)^rank`` for ``2 <= n``.

    Constructions are capped at 4096 ambient elements, so enumeration and
    coset arithmetic stay bounded.  The underlying ring acts by multiplication
    of each coordinate modulo ``n``.
    """

    modulus: int
    rank: int

    def __post_init__(self) -> None:
        _integer(self.modulus, "modulus")
        _integer(self.rank, "rank")
        if not 2 <= self.modulus <= _MAX_CARDINALITY:
            raise ValueError("modulus must lie in [2, 4096]")
        if not 1 <= self.rank <= _MAX_RANK:
            raise ValueError(f"rank must lie in [1, {_MAX_RANK}]")
        if self.modulus ** self.rank > _MAX_CARDINALITY:
            raise ValueError("module exceeds the finite ambient cardinality budget")

    @property
    def order(self) -> int:
        """Total number of module elements, ``n**rank``."""
        return self.modulus ** self.rank

    @property
    def zero(self) -> Vector:
        return (0,) * self.rank

    def normalize(self, vector: Iterable[int]) -> Vector:
        return tuple(value % self.modulus for value in _finite_coordinates(vector, self.rank))

    def add(self, left: Iterable[int], right: Iterable[int]) -> Vector:
        a, b = self.normalize(left), self.normalize(right)
        return tuple((x + y) % self.modulus for x, y in zip(a, b))

    def neg(self, vector: Iterable[int]) -> Vector:
        return tuple((-value) % self.modulus for value in self.normalize(vector))

    def subtract(self, left: Iterable[int], right: Iterable[int]) -> Vector:
        return self.add(left, self.neg(right))

    def scalar_multiply(self, scalar: int, vector: Iterable[int]) -> Vector:
        coefficient = _integer(scalar, "scalar") % self.modulus
        return tuple((coefficient * value) % self.modulus for value in self.normalize(vector))

    def elements(self) -> tuple[Vector, ...]:
        """Enumerate all canonical elements in lexicographic order."""
        return tuple(product(range(self.modulus), repeat=self.rank))

    def generated_submodule(self, *generators: Iterable[int]) -> FiniteSubmodule:
        """Smallest submodule containing the supplied generators.

        Repeated addition of each generator in the finite additive group
        produces its inverses as well, so a breadth-first search from zero
        computes the exact generated submodule over ``Z/nZ``.
        """
        if len(generators) > _MAX_GENERATORS:
            raise ValueError(f"at most {_MAX_GENERATORS} generators are supported")
        canonical = tuple(dict.fromkeys(self.normalize(g) for g in generators))
        seen = {self.zero}
        queue = deque([self.zero])
        while queue:
            vector = queue.popleft()
            for generator in canonical:
                candidate = self.add(vector, generator)
                if candidate not in seen:
                    seen.add(candidate)
                    queue.append(candidate)
        return FiniteSubmodule(self, tuple(sorted(seen)), _FACTORY_TOKEN)


@dataclass(frozen=True, slots=True, init=False)
class FiniteSubmodule:
    """A generated finite submodule, constructed by ``generated_submodule``."""

    ambient: ResidueFreeModule
    elements: tuple[Vector, ...]

    def __init__(self, ambient: ResidueFreeModule, elements: tuple[Vector, ...], token: object = None) -> None:
        if token is not _FACTORY_TOKEN:
            raise TypeError("construct submodules using ResidueFreeModule.generated_submodule")
        object.__setattr__(self, "ambient", ambient)
        object.__setattr__(self, "elements", elements)

    @property
    def order(self) -> int:
        return len(self.elements)

    @property
    def index(self) -> int:
        """The additive quotient cardinality ``|M/N|``."""
        return self.ambient.order // self.order

    def contains(self, vector: Iterable[int]) -> bool:
        return self.ambient.normalize(vector) in self.elements

    def quotient(self) -> QuotientModule:
        """Enumerate the disjoint additive cosets in ``M/N``."""
        seen: set[Vector] = set()
        classes: list[tuple[Vector, ...]] = []
        for representative in self.ambient.elements():
            if representative in seen:
                continue
            coset = tuple(sorted(self.ambient.add(representative, x)
                                 for x in self.elements))
            seen.update(coset)
            classes.append(coset)
        return QuotientModule(self, tuple(classes), _FACTORY_TOKEN)


@dataclass(frozen=True, slots=True, init=False)
class QuotientModule:
    """Canonical additive cosets with an induced ``Z/nZ`` module action.

    Classes are identified by deterministic integer indices ``0..order-1``;
    the index is a presentation choice, *not* a canonical module isomorphism.
    """

    submodule: FiniteSubmodule
    cosets: tuple[tuple[Vector, ...], ...]

    def __init__(self, submodule: FiniteSubmodule,
                 cosets: tuple[tuple[Vector, ...], ...], token: object = None) -> None:
        if token is not _FACTORY_TOKEN:
            raise TypeError("construct quotients using FiniteSubmodule.quotient")
        object.__setattr__(self, "submodule", submodule)
        object.__setattr__(self, "cosets", cosets)

    @property
    def order(self) -> int:
        return len(self.cosets)

    @property
    def zero_class(self) -> int:
        return self.class_of(self.submodule.ambient.zero)

    def class_of(self, vector: Iterable[int]) -> int:
        target = self.submodule.ambient.normalize(vector)
        for index, coset in enumerate(self.cosets):
            if target in coset:
                return index
        raise ArithmeticError("coset partition is incomplete")

    def _class(self, index: int) -> tuple[Vector, ...]:
        _integer(index, "coset index")
        if not 0 <= index < self.order:
            raise ValueError("coset index is out of range")
        return self.cosets[index]

    def add(self, left: int, right: int) -> int:
        a, b = self._class(left)[0], self._class(right)[0]
        return self.class_of(self.submodule.ambient.add(a, b))

    def neg(self, index: int) -> int:
        return self.class_of(self.submodule.ambient.neg(self._class(index)[0]))

    def scalar_multiply(self, scalar: int, index: int) -> int:
        return self.class_of(self.submodule.ambient.scalar_multiply(
            scalar, self._class(index)[0]))


__all__ = ["ResidueFreeModule", "FiniteSubmodule", "QuotientModule"]
