"""Outward-rounded intervals for exact operations on finite binary64 inputs.

Finite IEEE-754 binary64 endpoints are converted into their *exact* rational
values before each basic operation. Rational extrema are then rounded down and
up to surrounding binary64 numbers using ``math.nextafter`` as needed. This
avoids relying on a process-global floating-point rounding mode.

Bounds can become infinite when a mathematical result exceeds binary64 range.
Such intervals support inspection and intersection, but further arithmetic is
rejected because indeterminate extended-real operations need separate rules.
There are no certified transcendental functions in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
import sys


_MAX_FLOAT = sys.float_info.max
_MAX_EXACT = Fraction.from_float(_MAX_FLOAT)


def _endpoint(value: int | float, name: str) -> float:
    """Accept binary64 endpoints or exactly representable integers only."""
    if type(value) not in (int, float):
        raise TypeError(f"{name} must be a float or exactly representable integer")
    try:
        answer = float(value)
    except OverflowError as exc:
        raise ValueError(f"{name} is outside the binary64 range") from exc
    if math.isnan(answer):
        raise ValueError(f"{name} cannot be NaN")
    if type(value) is int and (not math.isfinite(answer) or int(answer) != value):
        raise ValueError(f"{name} integer is not exactly representable as binary64")
    return answer


def _finite_fraction(value: float) -> Fraction:
    if not math.isfinite(value):
        raise ValueError("interval arithmetic requires finite input bounds")
    return Fraction.from_float(value)


def _down(exact: Fraction) -> float:
    """Greatest binary64 value not exceeding an exact rational."""
    if exact < -_MAX_EXACT:
        return -math.inf
    if exact > _MAX_EXACT:
        return _MAX_FLOAT
    rounded = float(exact)
    if Fraction.from_float(rounded) > exact:
        rounded = math.nextafter(rounded, -math.inf)
    return -0.0 if rounded == 0.0 else rounded


def _up(exact: Fraction) -> float:
    """Least binary64 value not less than an exact rational."""
    if exact > _MAX_EXACT:
        return math.inf
    if exact < -_MAX_EXACT:
        return -_MAX_FLOAT
    rounded = float(exact)
    if Fraction.from_float(rounded) < exact:
        rounded = math.nextafter(rounded, math.inf)
    return 0.0 if rounded == 0.0 else rounded


def _enclose(lower: Fraction, upper: Fraction) -> Interval:
    assert lower <= upper
    return Interval(_down(lower), _up(upper))


@dataclass(frozen=True, slots=True)
class Interval:
    """Closed interval ``[lower, upper]`` of binary64 or infinite endpoints.

    Infinite endpoints denote an unbounded *enclosure*, not a finite input
    suitable for arithmetic. Use ``Interval.point(x)`` for a singleton.
    """

    lower: float
    upper: float

    def __post_init__(self) -> None:
        low = _endpoint(self.lower, "lower")
        high = _endpoint(self.upper, "upper")
        if low > high:
            raise ValueError("interval lower bound exceeds upper bound")
        object.__setattr__(self, "lower", low)
        object.__setattr__(self, "upper", high)

    @classmethod
    def point(cls, value: int | float) -> Interval:
        """Construct the interval containing exactly one binary64 number."""
        return cls(value, value)

    def _exact(self) -> tuple[Fraction, Fraction]:
        return _finite_fraction(self.lower), _finite_fraction(self.upper)

    def __add__(self, other: Interval) -> Interval:
        if not isinstance(other, Interval):
            return NotImplemented
        lo, hi = self._exact()
        other_lo, other_hi = other._exact()
        return _enclose(lo + other_lo, hi + other_hi)

    def __sub__(self, other: Interval) -> Interval:
        if not isinstance(other, Interval):
            return NotImplemented
        lo, hi = self._exact()
        other_lo, other_hi = other._exact()
        return _enclose(lo - other_hi, hi - other_lo)

    def __neg__(self) -> Interval:
        # Negating infinities is exact and keeps the interval order.
        return Interval(-self.upper, -self.lower)

    def __mul__(self, other: Interval) -> Interval:
        if not isinstance(other, Interval):
            return NotImplemented
        lo, hi = self._exact()
        other_lo, other_hi = other._exact()
        extrema = (lo * other_lo, lo * other_hi, hi * other_lo, hi * other_hi)
        return _enclose(min(extrema), max(extrema))

    def __truediv__(self, other: Interval) -> Interval:
        if not isinstance(other, Interval):
            return NotImplemented
        lo, hi = self._exact()
        other_lo, other_hi = other._exact()
        if other_lo <= 0 <= other_hi:
            raise ZeroDivisionError("cannot divide by an interval containing zero")
        extrema = (lo / other_lo, lo / other_hi, hi / other_lo, hi / other_hi)
        return _enclose(min(extrema), max(extrema))

    def contains(self, value: int | float) -> bool:
        """Whether this interval contains a representable binary64 value."""
        point = _endpoint(value, "value")
        return self.lower <= point <= self.upper

    def intersection(self, other: Interval) -> Interval | None:
        """Return the common closed interval, or None if disjoint."""
        if not isinstance(other, Interval):
            raise TypeError("intersection expects another Interval")
        low, high = max(self.lower, other.lower), min(self.upper, other.upper)
        return None if low > high else Interval(low, high)

    def hull(self, other: Interval) -> Interval:
        """Smallest closed interval containing both intervals."""
        if not isinstance(other, Interval):
            raise TypeError("hull expects another Interval")
        return Interval(min(self.lower, other.lower), max(self.upper, other.upper))

