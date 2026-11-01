"""Small exact rational Clifford (geometric) algebras of real signatures.

``CliffordAlgebra((1, 1, -1))`` has orthogonal generators e_0, e_1, e_2
with e_i**2 equal to the specified signs.  The basis mask ``1 << i``
represents e_i, and ordered mask bits represent e_{i_1} ... e_{i_k}.

This finite dimensional *rational coefficient* model covers nondegenerate
orthogonal real signatures but not arbitrary metrics or smooth Clifford
bundles.  It is exact at the coefficient level, subject to resource bounds.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction

_MAX_DIMENSION = 8
_MAX_COEFFICIENT_BITS = 4096


def _rational(value: int | Fraction) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("Clifford coefficients must be exact int or Fraction")
    exact = Fraction(value)
    if (exact.numerator.bit_length() > _MAX_COEFFICIENT_BITS or
            exact.denominator.bit_length() > _MAX_COEFFICIENT_BITS):
        raise ValueError("Clifford coefficient exceeds bit-length resource limit")
    return exact


@dataclass(frozen=True, slots=True, init=False)
class CliffordAlgebra:
    """Cl(p,q) with exact rational multivector coefficients.

    ``signature`` lists +1 and -1 for the squares of the ordered generators.
    Dimension is bounded at eight (256 basis blades).  Zero metric directions
    and arbitrary nonorthogonal metrics are deliberately unsupported.
    """

    signature: tuple[int, ...]

    def __init__(self, signature: Sequence[int]) -> None:
        try:
            components = tuple(signature)
        except TypeError as exc:
            raise TypeError("signature must be a sequence of +1/-1 integers") from exc
        if not 1 <= len(components) <= _MAX_DIMENSION:
            raise ValueError("Clifford dimension must be between 1 and 8")
        if any(type(sign) is not int for sign in components):
            raise TypeError("signature entries must be +1/-1 integers")
        if any(sign not in (-1, 1) for sign in components):
            raise ValueError("signature entries must be nonzero unit signs")
        object.__setattr__(self, "signature", components)

    @property
    def dimension(self) -> int:
        return len(self.signature)

    @property
    def positive(self) -> int:
        return self.signature.count(1)

    @property
    def negative(self) -> int:
        return self.signature.count(-1)

    def multivector(self, coefficients: Mapping[int, int | Fraction]) -> Multivector:
        return Multivector(self, coefficients)

    def scalar(self, coefficient: int | Fraction) -> Multivector:
        return self.multivector({0: coefficient})

    def zero(self) -> Multivector:
        return self.multivector({})

    def basis(self, mask: int) -> Multivector:
        """Return the positively ordered basis blade identified by ``mask``."""
        if isinstance(mask, bool) or not isinstance(mask, int):
            raise TypeError("blade masks must be nonnegative integers")
        if not 0 <= mask < (1 << self.dimension):
            raise ValueError("blade mask is outside the Clifford dimension")
        return self.multivector({mask: 1})

    def generator(self, index: int) -> Multivector:
        if isinstance(index, bool) or not isinstance(index, int):
            raise TypeError("generator index must be integer")
        if not 0 <= index < self.dimension:
            raise ValueError("generator index is outside the Clifford dimension")
        return self.basis(1 << index)


def _blade_geometric_product(left: int, right: int,
                             signature: tuple[int, ...]) -> tuple[int, int]:
    """Return mask and sign for two orthogonal ordered basis blade products."""
    sign = 1
    remaining = left
    while remaining:
        lowest = remaining & -remaining
        i = lowest.bit_length() - 1
        if (right & (lowest - 1)).bit_count() & 1:
            sign = -sign
        remaining ^= lowest
    overlap = left & right
    while overlap:
        lowest = overlap & -overlap
        sign *= signature[lowest.bit_length() - 1]
        overlap ^= lowest
    return left ^ right, sign


@dataclass(frozen=True, slots=True, init=False)
class Multivector:
    """An immutable sparse sum of ordered basis blades over the rationals.

    The stored ``coefficients`` are canonical, sorted ``(mask, Fraction)``
    pairs with zero coefficients removed.  All operations preserve the owning
    signature and reject mixing distinct algebras.
    """

    algebra: CliffordAlgebra
    coefficients: tuple[tuple[int, Fraction], ...] = field(repr=True)

    def __init__(self, algebra: CliffordAlgebra,
                 coefficients: Mapping[int, int | Fraction]) -> None:
        if not isinstance(algebra, CliffordAlgebra):
            raise TypeError("algebra must be a CliffordAlgebra")
        if not isinstance(coefficients, Mapping):
            raise TypeError("coefficients must map basis masks to exact rationals")
        terms: list[tuple[int, Fraction]] = []
        for mask, value in coefficients.items():
            if isinstance(mask, bool) or not isinstance(mask, int):
                raise TypeError("basis mask must be a nonnegative integer")
            if not 0 <= mask < 1 << algebra.dimension:
                raise ValueError("basis mask outside Clifford dimension")
            exact = _rational(value)
            if exact:
                terms.append((mask, exact))
        object.__setattr__(self, "algebra", algebra)
        object.__setattr__(self, "coefficients", tuple(sorted(terms)))

    def as_dict(self) -> dict[int, Fraction]:
        """Return a fresh mapping of masks to exact coefficients."""
        return dict(self.coefficients)

    def _other(self, other: object) -> Multivector | object:
        if isinstance(other, Multivector):
            if other.algebra != self.algebra:
                raise ValueError("cannot combine different Clifford signatures")
            return other
        if isinstance(other, (int, Fraction)) and not isinstance(other, bool):
            return self.algebra.scalar(other)
        return NotImplemented

    def __add__(self, other: object) -> Multivector:
        rhs = self._other(other)
        if rhs is NotImplemented:
            return NotImplemented
        result = self.as_dict()
        for mask, coefficient in rhs.coefficients:
            result[mask] = result.get(mask, Fraction(0)) + coefficient
        return self.algebra.multivector(result)

    __radd__ = __add__

    def __neg__(self) -> Multivector:
        return self.algebra.multivector({mask: -c for mask, c in self.coefficients})

    def __sub__(self, other: object) -> Multivector:
        rhs = self._other(other)
        if rhs is NotImplemented:
            return NotImplemented
        return self + (-rhs)

    def __rsub__(self, other: object) -> Multivector:
        rhs = self._other(other)
        if rhs is NotImplemented:
            return NotImplemented
        return rhs - self

    def _product(self, rhs: Multivector, *, product: str) -> Multivector:
        result: dict[int, Fraction] = {}
        for left_mask, left_coeff in self.coefficients:
            left_grade = left_mask.bit_count()
            for right_mask, right_coeff in rhs.coefficients:
                right_grade = right_mask.bit_count()
                mask, sign = _blade_geometric_product(left_mask, right_mask,
                                                       self.algebra.signature)
                if product == "wedge" and mask.bit_count() != left_grade + right_grade:
                    continue
                if product == "left" and (left_grade > right_grade or
                                          mask.bit_count() != right_grade - left_grade):
                    continue
                result[mask] = result.get(mask, Fraction(0)) + sign * left_coeff * right_coeff
        return self.algebra.multivector(result)

    def __mul__(self, other: object) -> Multivector:
        rhs = self._other(other)
        if rhs is NotImplemented:
            return NotImplemented
        return self._product(rhs, product="geometric")

    def __rmul__(self, other: object) -> Multivector:
        lhs = self._other(other)
        if lhs is NotImplemented:
            return NotImplemented
        return lhs._product(self, product="geometric")

    def wedge(self, other: Multivector) -> Multivector:
        """Exterior product A∧B: the grade-r+s component for homogeneous blades."""
        rhs = self._other(other)
        if rhs is NotImplemented:
            return NotImplemented
        return self._product(rhs, product="wedge")

    def left_contract(self, other: Multivector) -> Multivector:
        """Hestenes-style left contraction: ⟨AB⟩_(s-r) for grades r≤s.

        This is the *left* contraction, not a symmetric full inner product.
        In particular, a bivector contracting a vector contributes zero.
        """
        rhs = self._other(other)
        if rhs is NotImplemented:
            return NotImplemented
        return self._product(rhs, product="left")

    def grade(self, degree: int) -> Multivector:
        if isinstance(degree, bool) or not isinstance(degree, int):
            raise TypeError("grade must be a nonnegative integer")
        if degree < 0:
            raise ValueError("grade must be nonnegative")
        return self.algebra.multivector({mask: c for mask, c in self.coefficients
                                         if mask.bit_count() == degree})

    def scalar_part(self) -> Fraction:
        return self.as_dict().get(0, Fraction(0))

    def reverse(self) -> Multivector:
        """Reversion anti-automorphism: grade k receives (-1)^(k(k-1)/2)."""
        return self.algebra.multivector({mask: (-1 if (mask.bit_count() *
                                  (mask.bit_count()-1)//2) % 2 else 1) * c
                                        for mask, c in self.coefficients})

    def grade_involution(self) -> Multivector:
        """Grade automorphism: negate all odd-grade terms."""
        return self.algebra.multivector({mask: -c if mask.bit_count() & 1 else c
                                        for mask, c in self.coefficients})

    def clifford_conjugate(self) -> Multivector:
        """Composition of reversion and grade involution."""
        return self.reverse().grade_involution()

    def reversion_quadratic(self) -> Fraction:
        """Scalar part of reverse(A) A; NOT a positive norm for mixed signatures.

        For an all-positive signature, this equals the sum of the squared
        coefficients.  Indefinite signatures can yield negative or zero.
        """
        return (self.reverse() * self).scalar_part()
