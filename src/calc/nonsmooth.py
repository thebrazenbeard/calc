"""Exact, one-dimensional convex piecewise-linear and proximal calculus.

This module is deliberately narrower than general nonsmooth analysis. It
supports finite maxima of affine functions with rational coefficients and a
few closed-form proximal operations; it does not calculate generalized
gradients for arbitrary nonconvex functions.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import islice
from typing import Iterable


_MAX_AFFINES = 1024
_MAX_BITS = 4096


def _rational(value: int | Fraction, *, name: str) -> Fraction:
    """Validate an exact rational scalar and bound its input size."""

    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an integer or Fraction (not a float)")
    result = Fraction(value)
    if (result.numerator.bit_length() > _MAX_BITS
            or result.denominator.bit_length() > _MAX_BITS):
        raise ValueError(f"{name} exceeds the {_MAX_BITS}-bit input limit")
    return result


class ConvexPiecewiseLinear:
    """A finite scalar convex function ``f(x) = max_i(m_i*x+b_i)``.

    ``affines`` is a nonempty iterable of ``(slope, intercept)`` pairs, all
    exact integers or :class:`fractions.Fraction` objects. The function is
    convex by construction, so arbitrary nonconvex piecewise functions are
    not silently accepted as if they were convex.

    Only the one-dimensional convex subdifferential is exposed. In this
    setting, and only because ``f`` is convex, its Clarke generalized
    gradient equals the convex subdifferential.
    """

    __slots__ = ("affines",)

    def __init__(self, affines: Iterable[tuple[int | Fraction, int | Fraction]]):
        if not isinstance(affines, Iterable):
            raise TypeError("affines must be an iterable of (slope, intercept) pairs")
        limited = tuple(islice(affines, _MAX_AFFINES + 1))
        if not limited:
            raise ValueError("at least one affine function is required")
        if len(limited) > _MAX_AFFINES:
            raise ValueError(f"at most {_MAX_AFFINES} affine functions are supported")
        normalized: list[tuple[Fraction, Fraction]] = []
        for i, pair in enumerate(limited):
            try:
                slope, intercept = pair
            except (TypeError, ValueError) as exc:
                raise TypeError(f"affine {i} must be a pair (slope, intercept)") from exc
            normalized.append((
                _rational(slope, name=f"slope {i}"),
                _rational(intercept, name=f"intercept {i}"),
            ))
        self.affines = tuple(normalized)

    def value(self, x: int | Fraction) -> Fraction:
        """Evaluate ``f(x)`` exactly."""

        x = _rational(x, name="x")
        return max(slope * x + intercept for slope, intercept in self.affines)

    def _active_slopes(self, x: Fraction) -> tuple[Fraction, ...]:
        values = tuple(slope * x + intercept for slope, intercept in self.affines)
        maximum = max(values)
        return tuple(
            slope for (slope, _), value in zip(self.affines, values)
            if value == maximum
        )

    def subdifferential(self, x: int | Fraction) -> tuple[Fraction, Fraction]:
        """Return the closed interval of convex subgradients at ``x``.

        It is the convex hull of slopes of the active affine functions,
        represented in one dimension as ``(minimum, maximum)``.
        """

        slopes = self._active_slopes(_rational(x, name="x"))
        return (min(slopes), max(slopes))

    def clarke_subdifferential(self, x: int | Fraction) -> tuple[Fraction, Fraction]:
        """Return the Clarke generalized gradient for *this convex class*.

        For finite convex functions, the Clarke gradient agrees with the
        ordinary convex subdifferential. No nonconvex Clarke rules are
        implemented here.
        """

        return self.subdifferential(x)

    def directional_derivative(self, x: int | Fraction, direction: int | Fraction) -> Fraction:
        """Return the one-sided directional derivative ``f'(x; direction)``."""

        x = _rational(x, name="x")
        direction = _rational(direction, name="direction")
        return max(slope * direction for slope in self._active_slopes(x))


_ABSOLUTE = ConvexPiecewiseLinear(((1, 0), (-1, 0)))
_HINGE = ConvexPiecewiseLinear(((0, 0), (1, 0)))


def absolute_subgradient(x: int | Fraction) -> tuple[Fraction, Fraction]:
    """Convex subgradient interval of ``|x|`` (``[-1,1]`` at zero)."""

    return _ABSOLUTE.subdifferential(x)


def hinge_subgradient(x: int | Fraction) -> tuple[Fraction, Fraction]:
    """Convex subgradient interval of ``max(0,x)`` (``[0,1]`` at zero)."""

    return _HINGE.subdifferential(x)


def soft_threshold(value: int | Fraction, weight: int | Fraction) -> Fraction:
    """Compute ``prox_{weight*|.|}(value)`` for ``weight >= 0``.

    This is the minimizer of ``weight*|x| + (x-value)**2/2``.
    """

    value = _rational(value, name="value")
    weight = _rational(weight, name="weight")
    if weight < 0:
        raise ValueError("weight must be nonnegative")
    if value > weight:
        return value - weight
    if value < -weight:
        return value + weight
    return Fraction(0)


def prox_hinge(value: int | Fraction, weight: int | Fraction) -> Fraction:
    """Compute ``prox_{weight*max(0,.)}(value)`` for ``weight >= 0``."""

    value = _rational(value, name="value")
    weight = _rational(weight, name="weight")
    if weight < 0:
        raise ValueError("weight must be nonnegative")
    if value <= 0:
        return value
    if value <= weight:
        return Fraction(0)
    return value - weight


def moreau_envelope_absolute(value: int | Fraction, weight: int | Fraction) -> Fraction:
    """Compute ``min_x |x| + (x-value)**2/(2*weight)`` for ``weight > 0``.

    This is the Huber-function form of the Moreau envelope of ``|.|``.
    """

    value = _rational(value, name="value")
    weight = _rational(weight, name="weight")
    if weight <= 0:
        raise ValueError("weight must be strictly positive")
    if abs(value) <= weight:
        return value * value / (2 * weight)
    return abs(value) - weight / 2


__all__ = [
    "ConvexPiecewiseLinear",
    "absolute_subgradient",
    "hinge_subgradient",
    "soft_threshold",
    "prox_hinge",
    "moreau_envelope_absolute",
]
