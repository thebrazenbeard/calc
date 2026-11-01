"""Exact, bounded, *truncated* differential polynomial algebra over Q.

Jet symbols ``u_0, ..., u_m`` stand formally for a function and its
successive ordinary derivatives.  The total derivation maps ``u_j`` to
``u_(j+1)`` and satisfies Leibniz's rule.  Since the jet ring is truncated,
the derivation rejects terms depending on the highest represented jet.

No numerical functions, ODE solutions, or differential-field algorithms are
provided.  All coefficients are exact rational values.
"""

from __future__ import annotations

from collections.abc import Mapping, Iterable
from dataclasses import dataclass, field
from fractions import Fraction
from types import NotImplementedType

_MAX_JET = 8
_MAX_DEGREE = 16
_MAX_TERMS = 1024
_MAX_PRODUCT_PAIRS = 20000
_MAX_BITS = 4096

Exponent = tuple[int, ...]
Rational = int | Fraction


def _rational(value: Rational, name: str = "coefficient") -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be int or Fraction, not float/bool")
    number = Fraction(value)
    if max(number.numerator.bit_length(), number.denominator.bit_length()) > _MAX_BITS:
        raise ValueError(f"{name} exceeds exact-rational bit limit")
    return number


def _jet_count(max_jet: int) -> int:
    if isinstance(max_jet, bool) or not isinstance(max_jet, int):
        raise TypeError("max_jet must be a nonnegative integer")
    if not 0 <= max_jet <= _MAX_JET:
        raise ValueError(f"max_jet must be between 0 and {_MAX_JET}")
    return max_jet + 1


@dataclass(frozen=True, slots=True, init=False)
class JetPolynomial:
    """A finite rational polynomial in jet symbols ``u_0, ..., u_m``.

    The ``terms`` mapping uses exponent tuples in increasing jet order:
    ``JetPolynomial(2, {(2,1,0): 3})`` denotes ``3 u_0² u_1``.
    Zero terms are discarded.  The canonical internal data is immutable.
    """

    max_jet: int
    coefficients: tuple[tuple[Exponent, Fraction], ...] = field(repr=True)

    def __init__(self, max_jet: int,
                 terms: Mapping[Exponent, Rational] | None = None) -> None:
        size = _jet_count(max_jet)
        if terms is None:
            terms = {}
        if not isinstance(terms, Mapping):
            raise TypeError("terms must map exponent tuples to exact rationals")
        if len(terms) > _MAX_TERMS:
            raise ValueError("polynomial exceeds monomial limit")
        canonical: list[tuple[Exponent, Fraction]] = []
        for powers, coefficient in terms.items():
            if not isinstance(powers, tuple) or len(powers) != size:
                raise ValueError(f"each exponent tuple must have length {size}")
            if any(type(exponent) is not int or exponent < 0 for exponent in powers):
                raise TypeError("jet exponents must be nonnegative integers, not booleans")
            if sum(powers) > _MAX_DEGREE:
                raise ValueError("polynomial total degree exceeds limit")
            value = _rational(coefficient)
            if value:
                canonical.append((powers, value))
        object.__setattr__(self, "max_jet", max_jet)
        object.__setattr__(self, "coefficients", tuple(sorted(canonical)))

    @classmethod
    def constant(cls, value: Rational, *, max_jet: int = 4) -> JetPolynomial:
        size = _jet_count(max_jet)
        return cls(max_jet, {(0,) * size: _rational(value)})

    @classmethod
    def jet(cls, index: int, *, max_jet: int = 4) -> JetPolynomial:
        size = _jet_count(max_jet)
        if type(index) is not int or not 0 <= index < size:
            raise ValueError(f"jet index must be between 0 and {size - 1}")
        powers = [0] * size
        powers[index] = 1
        return cls(max_jet, {tuple(powers): Fraction(1)})

    def as_dict(self) -> dict[Exponent, Fraction]:
        """Return a copy of the canonical term mapping."""
        return dict(self.coefficients)

    def _coerce(self, other: object) -> JetPolynomial | NotImplementedType:
        if isinstance(other, JetPolynomial):
            if self.max_jet != other.max_jet:
                raise ValueError("jet polynomials use different ambient jet orders")
            return other
        if isinstance(other, (int, Fraction)) and not isinstance(other, bool):
            return JetPolynomial.constant(other, max_jet=self.max_jet)
        return NotImplemented

    def __add__(self, other: object) -> JetPolynomial:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        result = self.as_dict()
        for powers, coefficient in value.coefficients:
            result[powers] = result.get(powers, Fraction(0)) + coefficient
        return JetPolynomial(self.max_jet, result)

    __radd__ = __add__

    def __neg__(self) -> JetPolynomial:
        return JetPolynomial(self.max_jet, {e: -c for e, c in self.coefficients})

    def __sub__(self, other: object) -> JetPolynomial:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return self + (-value)

    def __rsub__(self, other: object) -> JetPolynomial:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return value - self

    def __mul__(self, other: object) -> JetPolynomial:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        if len(self.coefficients) * len(value.coefficients) > _MAX_PRODUCT_PAIRS:
            raise ValueError("polynomial multiplication exceeds work limit")
        result: dict[Exponent, Fraction] = {}
        for left, c in self.coefficients:
            for right, d in value.coefficients:
                powers = tuple(x + y for x, y in zip(left, right))
                if sum(powers) > _MAX_DEGREE:
                    raise ValueError("polynomial product exceeds degree limit")
                result[powers] = result.get(powers, Fraction(0)) + c * d
        return JetPolynomial(self.max_jet, result)

    __rmul__ = __mul__

    def __truediv__(self, scalar: Rational) -> JetPolynomial:
        denominator = _rational(scalar, "divisor")
        if not denominator:
            raise ZeroDivisionError("cannot divide by zero")
        return JetPolynomial(self.max_jet,
                             {powers: coefficient / denominator
                              for powers, coefficient in self.coefficients})

    def __pow__(self, exponent: int) -> JetPolynomial:
        if type(exponent) is not int:
            raise TypeError("polynomial exponent must be an integer")
        if not 0 <= exponent <= _MAX_DEGREE:
            raise ValueError("polynomial exponent exceeds nonnegative degree bound")
        output = JetPolynomial.constant(1, max_jet=self.max_jet)
        base = self
        power = exponent
        while power:
            if power & 1:
                output = output * base
            power //= 2
            if power:
                base = base * base
        return output

    def partial_derivative(self, jet_index: int) -> JetPolynomial:
        """Ordinary formal partial derivative w.r.t. one jet indeterminate."""
        size = self.max_jet + 1
        if type(jet_index) is not int or not 0 <= jet_index < size:
            raise ValueError("partial jet index outside represented jets")
        result: dict[Exponent, Fraction] = {}
        for exponents, coefficient in self.coefficients:
            power = exponents[jet_index]
            if not power:
                continue
            powers = list(exponents)
            powers[jet_index] -= 1
            key = tuple(powers)
            result[key] = result.get(key, Fraction(0)) + power * coefficient
        return JetPolynomial(self.max_jet, result)

    def total_derivative(self) -> JetPolynomial:
        """Derivation ``D=Σ u_(j+1)∂/∂u_j`` on representable terms.

        The finite jet ring is not closed under D.  If the polynomial depends
        on the top jet ``u_m``, D needs the unrepresented ``u_(m+1)`` and
        raises ValueError rather than silently discarding that contribution.
        """
        result: dict[Exponent, Fraction] = {}
        for exponents, coefficient in self.coefficients:
            if exponents[-1]:
                raise ValueError("total derivative needs a higher unrepresented jet")
            for index in range(self.max_jet):
                power = exponents[index]
                if not power:
                    continue
                powers = list(exponents)
                powers[index] -= 1
                powers[index + 1] += 1
                key = tuple(powers)
                result[key] = result.get(key, Fraction(0)) + power * coefficient
        return JetPolynomial(self.max_jet, result)

    def evaluate(self, jet_values: Iterable[Rational]) -> Fraction:
        """Substitute exact rational values for every jet symbol."""
        iterator = iter(jet_values)
        values = []
        for _ in range(self.max_jet + 2):
            try:
                values.append(next(iterator))
            except StopIteration:
                break
        if len(values) != self.max_jet + 1:
            raise ValueError("expected exactly max_jet + 1 jet values")
        numbers = tuple(_rational(value, "jet value") for value in values)
        total = sum((coefficient *
                     (math_prod(numbers[i] ** exponent
                                for i, exponent in enumerate(exponents)))
                     for exponents, coefficient in self.coefficients), Fraction(0))
        return _rational(total, "evaluation")


def math_prod(values: Iterable[Fraction]) -> Fraction:
    """Exact rational product (internal; independent of binary floating math)."""
    result = Fraction(1)
    for value in values:
        result *= value
    return result


__all__ = ["JetPolynomial"]
