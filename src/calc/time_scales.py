"""Delta/nabla calculus on **finite, isolated** time scales.

A time scale is generally any nonempty closed subset of the real numbers.
This module implements only an explicitly enumerated finite subset. Every
nonterminal point is right-scattered, and every noninitial point is
left-scattered; no limits at right-dense/left-dense accumulation points or
continuous derivatives are modeled.

For t in T except max(T), sigma(t) is its successor, mu(t)=sigma(t)-t, and
f^Delta(t)=(f(sigma(t))-f(t))/mu(t). Nabla operators use predecessors.
Finite delta/nabla integrals are the corresponding left/right weighted sums.

The stored grid and callback values are represented as finite binary floats.
This is a finite-grid demonstration, not exact symbolic time-scale calculus.
"""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Real
from collections.abc import Iterable, Callable

RealFunction = Callable[[float], Real]
_MAX_GRID_POINTS = 100_000


def _real(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number, not bool or text")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} is outside finite float range") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _evaluate(function: RealFunction, point: float) -> float:
    if not callable(function):
        raise TypeError("function must be callable")
    return _real(function(point), "function value")


def _quotient(numerator: float, denominator: float) -> float:
    value = numerator / denominator
    if not math.isfinite(value):
        raise ValueError("finite-difference quotient is not representable")
    return value


def _difference_quotient(right_value: float, left_value: float, gap: float) -> float:
    numerator = right_value - left_value
    if math.isfinite(numerator):
        return _quotient(numerator, gap)
    # Finite endpoint values can have an overflowing difference even when
    # their quotient is finite. Use exact binary64 rational values only for
    # this exceptional path, avoiding an unnecessary false overflow.
    exact = ((Fraction.from_float(right_value) - Fraction.from_float(left_value))
             / Fraction.from_float(gap))
    try:
        value = float(exact)
    except OverflowError as exc:
        raise ValueError("finite-difference quotient is not representable") from exc
    if not math.isfinite(value):
        raise ValueError("finite-difference quotient is not representable")
    return value


def _finite_weighted_sum(pairs: list[tuple[float, float]], sign: float,
                         name: str) -> float:
    """Weighted float sum, recovering finite cancellation after term overflow."""
    terms = [value * gap for value, gap in pairs]
    if all(math.isfinite(term) for term in terms):
        try:
            result = sign * math.fsum(terms)
        except (OverflowError, ValueError):
            pass
        else:
            if math.isfinite(result):
                return result
    exact = sum((Fraction.from_float(value) * Fraction.from_float(gap)
                 for value, gap in pairs), Fraction(0))
    try:
        result = float(exact) * sign
    except OverflowError as exc:
        raise ValueError(f"{name} is outside finite float range") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} is outside finite float range")
    return result


@dataclass(frozen=True, slots=True, init=False)
class FiniteTimeScale:
    """Immutable ordered nonempty finite time scale, with float coordinates.

    Constructors reject duplicate, unsorted, unrepresentable or too-large
    gaps. The singleton time scale is valid but supports no derivatives.
    """

    points: tuple[float, ...]

    def __init__(self, points: Iterable[Real]):
        if isinstance(points, (str, bytes)):
            raise TypeError("points must be an iterable of finite real numbers")
        try:
            raw = tuple(points)
        except TypeError as exc:
            raise TypeError("points must be an iterable of finite real numbers") from exc
        if not raw:
            raise ValueError("time scale cannot be empty")
        if len(raw) > _MAX_GRID_POINTS:
            raise ValueError("finite time scale exceeds grid point limit")
        values = tuple(_real(point, "time-scale point") for point in raw)
        for left, right in zip(values, values[1:]):
            delta = right - left
            if not 0.0 < delta < math.inf:
                raise ValueError("points must remain strictly increasing with finite representable gaps")
        object.__setattr__(self, "points", values)

    def _index(self, point: Real) -> int:
        t = _real(point, "point")
        index = bisect_left(self.points, t)
        if index >= len(self.points) or self.points[index] != t:
            raise ValueError("point must belong exactly to the finite time scale")
        return index

    def sigma(self, point: Real) -> float:
        """Forward jump σ(t); by convention σ(max T)=max T."""
        i = self._index(point)
        return self.points[min(i + 1, len(self.points) - 1)]

    def rho(self, point: Real) -> float:
        """Backward jump ρ(t); by convention ρ(min T)=min T."""
        i = self._index(point)
        return self.points[max(i - 1, 0)]

    def mu(self, point: Real) -> float:
        """Forward graininess μ(t)=σ(t)-t, zero at terminal point."""
        i = self._index(point)
        return 0.0 if i == len(self.points) - 1 else self.points[i + 1] - self.points[i]

    def nu(self, point: Real) -> float:
        """Backward graininess ν(t)=t-ρ(t), zero at initial point."""
        i = self._index(point)
        return 0.0 if i == 0 else self.points[i] - self.points[i - 1]

    def is_right_scattered(self, point: Real) -> bool:
        return self._index(point) < len(self.points) - 1

    def is_left_scattered(self, point: Real) -> bool:
        return self._index(point) > 0

    def delta_derivative(self, function: RealFunction, point: Real) -> float:
        """Discrete delta derivative at nonterminal right-scattered t.

        At max(T), no delta derivative is defined on a finite time scale.
        """
        i = self._index(point)
        if i == len(self.points) - 1:
            raise ValueError("delta derivative undefined at finite time-scale maximum")
        left, right = self.points[i:i + 2]
        return _difference_quotient(_evaluate(function, right), _evaluate(function, left), right - left)

    def nabla_derivative(self, function: RealFunction, point: Real) -> float:
        """Discrete nabla derivative at noninitial left-scattered t."""
        i = self._index(point)
        if i == 0:
            raise ValueError("nabla derivative undefined at finite time-scale minimum")
        left, right = self.points[i - 1:i + 1]
        return _difference_quotient(_evaluate(function, right), _evaluate(function, left), right - left)

    def delta_integral(self, function: RealFunction, start: Real, stop: Real) -> float:
        """Finite delta integral Σ f(t_i)·(t_{i+1}-t_i), left-sampled.

        Endpoints must lie in T. Reversed bounds negate the integral.
        """
        if not callable(function):
            raise TypeError("function must be callable")
        i, j = self._index(start), self._index(stop)
        if i == j:
            return 0.0
        sign = 1.0
        if i > j:
            i, j, sign = j, i, -1.0
        pairs: list[tuple[float, float]] = []
        for k in range(i, j):
            pairs.append((_evaluate(function, self.points[k]),
                          self.points[k + 1] - self.points[k]))
        return _finite_weighted_sum(pairs, sign, "delta integral")

    def nabla_integral(self, function: RealFunction, start: Real, stop: Real) -> float:
        """Finite nabla integral Σ f(t_i)·(t_i-t_{i-1}), right-sampled."""
        if not callable(function):
            raise TypeError("function must be callable")
        i, j = self._index(start), self._index(stop)
        if i == j:
            return 0.0
        sign = 1.0
        if i > j:
            i, j, sign = j, i, -1.0
        pairs: list[tuple[float, float]] = []
        for k in range(i + 1, j + 1):
            pairs.append((_evaluate(function, self.points[k]),
                          self.points[k] - self.points[k - 1]))
        return _finite_weighted_sum(pairs, sign, "nabla integral")
