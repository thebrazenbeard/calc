"""Exact finite-jump examples of one-dimensional distributional calculus.

A compactly supported piecewise-constant rational function has a weak
derivative consisting entirely of signed Dirac atoms at its jumps.  This
module computes those atoms and pairings with polynomial germs exactly.

This is deliberately *not* a general distribution, Schwartz-test-function,
or weak PDE framework.  All coordinates and coefficients are rational.
"""

from __future__ import annotations

from bisect import bisect_right
from fractions import Fraction
from itertools import islice
from typing import Iterable


Rational = int | Fraction
_MAX_INTERVALS = 2048
_MAX_ATOMS = 4096
_MAX_POLYNOMIAL_TERMS = 65
_MAX_INPUT_BITS = 4096


def _rational(value: Rational, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an int or Fraction (not a float)")
    result = Fraction(value)
    if (result.numerator.bit_length() > _MAX_INPUT_BITS
            or result.denominator.bit_length() > _MAX_INPUT_BITS):
        raise ValueError(f"{name} exceeds the {_MAX_INPUT_BITS}-bit limit")
    return result


def _bounded_items(items: Iterable[object], maximum: int, name: str) -> tuple[object, ...]:
    if isinstance(items, (str, bytes)):
        raise TypeError(f"{name} must be an iterable, not text")
    try:
        result = tuple(islice(iter(items), maximum + 1))
    except TypeError as exc:
        raise TypeError(f"{name} must be an iterable") from exc
    if len(result) > maximum:
        raise ValueError(f"{name} exceeds the {maximum}-item limit")
    return result


def _polynomial(coefficients: Iterable[Rational]) -> tuple[Fraction, ...]:
    """Coefficients in ascending powers; empty means the zero polynomial."""
    values = _bounded_items(coefficients, _MAX_POLYNOMIAL_TERMS, "polynomial")
    return tuple(_rational(value, f"polynomial coefficient {i}")
                 for i, value in enumerate(values))


def _evaluate_polynomial(coefficients: tuple[Fraction, ...], x: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(coefficients):
        result = result * x + coefficient
    return result


class DiracMeasure:
    """A finite rational signed atomic measure ``Σ weight_i δ_{point_i}``.

    Repeated locations are combined, zero weights removed and atoms sorted.
    The resulting object is immutable in its public representation.
    """

    __slots__ = ("atoms",)

    def __init__(self, atoms: Iterable[tuple[Rational, Rational]] = ()):
        raw = _bounded_items(atoms, _MAX_ATOMS, "atoms")
        contributions: dict[Fraction, Fraction] = {}
        for i, atom in enumerate(raw):
            try:
                point, weight = atom
            except (TypeError, ValueError) as exc:
                raise TypeError(f"atom {i} must contain (point, weight)") from exc
            point = _rational(point, f"atom {i} location")
            weight = _rational(weight, f"atom {i} weight")
            contributions[point] = contributions.get(point, Fraction(0)) + weight
        self.atoms: tuple[tuple[Fraction, Fraction], ...] = tuple(
            (point, weight) for point, weight in sorted(contributions.items())
            if weight
        )

    def total_mass(self) -> Fraction:
        """Signed mass of the entire real line."""
        return sum((weight for _, weight in self.atoms), Fraction(0))

    def total_variation(self) -> Fraction:
        """Total variation of this finite signed atomic measure."""
        return sum((abs(weight) for _, weight in self.atoms), Fraction(0))

    def pair_polynomial(self, coefficients: Iterable[Rational]) -> Fraction:
        """Return ``Σ_i weight_i p(point_i)`` for rational polynomial ``p``.

        A global nonzero polynomial is not compactly supported.  For the
        compact atomic support here, this pairing is equivalently against a
        smooth compactly supported cutoff equal to ``p`` near the atoms.
        """
        p = _polynomial(coefficients)
        return sum((weight * _evaluate_polynomial(p, point)
                    for point, weight in self.atoms), Fraction(0))

    def shift(self, displacement: Rational) -> "DiracMeasure":
        """Push the measure forward under ``x ↦ x + displacement``."""
        d = _rational(displacement, "displacement")
        return DiracMeasure((point + d, weight) for point, weight in self.atoms)

    def scale(self, factor: Rational) -> "DiracMeasure":
        factor = _rational(factor, "factor")
        return DiracMeasure((point, factor * weight) for point, weight in self.atoms)


class PiecewiseConstant:
    """A compactly supported exact rational step function.

    Input triples ``(left, right, height)`` define half-open intervals
    ``[left, right)``; outside their union, ``f=0``.  Intervals must be in
    nondecreasing left-endpoint order, have strictly positive widths, and
    must not overlap.  Adjacent intervals are allowed. Zero-height intervals
    are ignored; touching intervals of equal height are merged.
    """

    __slots__ = ("intervals", "_starts")

    def __init__(self, intervals: Iterable[tuple[Rational, Rational, Rational]] = ()):
        raw = _bounded_items(intervals, _MAX_INTERVALS, "intervals")
        normalized: list[tuple[Fraction, Fraction, Fraction]] = []
        previous_right: Fraction | None = None
        for i, interval in enumerate(raw):
            try:
                left, right, height = interval
            except (TypeError, ValueError) as exc:
                raise TypeError(f"interval {i} must contain (left, right, height)") from exc
            left = _rational(left, f"interval {i} left")
            right = _rational(right, f"interval {i} right")
            height = _rational(height, f"interval {i} height")
            if left >= right:
                raise ValueError(f"interval {i} must have left < right")
            if previous_right is not None and left < previous_right:
                raise ValueError("intervals must be ordered and nonoverlapping")
            previous_right = right
            if height == 0:
                continue
            if normalized and normalized[-1][1] == left and normalized[-1][2] == height:
                previous_left, _, _ = normalized[-1]
                normalized[-1] = (previous_left, right, height)
            else:
                normalized.append((left, right, height))
        self.intervals: tuple[tuple[Fraction, Fraction, Fraction], ...] = tuple(normalized)
        self._starts: tuple[Fraction, ...] = tuple(left for left, _, _ in normalized)

    def value_at(self, point: Rational) -> Fraction:
        """Return the selected half-open representative at rational ``point``."""
        point = _rational(point, "point")
        index = bisect_right(self._starts, point) - 1
        if index >= 0:
            _, right, height = self.intervals[index]
            if point < right:
                return height
        return Fraction(0)

    def integral(self) -> Fraction:
        """Exact integral over the real line."""
        return sum((height * (right - left) for left, right, height in self.intervals),
                   Fraction(0))

    def pair_polynomial(self, coefficients: Iterable[Rational]) -> Fraction:
        """Compute ``∫_R f(x)p(x) dx`` using rational antiderivatives.

        Although a polynomial need not be globally compactly supported, the
        result is the same for any compactly supported smooth test function
        that agrees with ``p`` near this function's compact support.
        """
        p = _polynomial(coefficients)
        total = Fraction(0)
        for left, right, height in self.intervals:
            for degree, coefficient in enumerate(p):
                total += (height * coefficient * (right**(degree + 1) - left**(degree + 1))
                          / (degree + 1))
        return total

    def weak_derivative(self) -> DiracMeasure:
        """Return the distributional derivative ``Df = Σ jumps · δ``.

        A positive jump at an interval's left endpoint contributes +height;
        the drop to zero at its right endpoint contributes -height. Adjacent
        intervals aggregate correctly. No classical boundary derivative is
        asserted at discontinuities.
        """
        return DiracMeasure(
            atom for left, right, height in self.intervals
            for atom in ((left, height), (right, -height))
        )

    def integration_by_parts_residual(self, coefficients: Iterable[Rational]) -> Fraction:
        """Return ``⟨Df,p⟩ + ⟨f,p'⟩``; mathematically it is exactly zero."""
        p = _polynomial(coefficients)
        derivative = tuple((k + 1) * p[k + 1] for k in range(len(p) - 1))
        return self.weak_derivative().pair_polynomial(p) + self.pair_polynomial(derivative)

    def shift(self, displacement: Rational) -> "PiecewiseConstant":
        """Translate the step function by ``f_new(x)=f(x-displacement)``."""
        d = _rational(displacement, "displacement")
        return PiecewiseConstant((left + d, right + d, height)
                                 for left, right, height in self.intervals)

    def scale(self, factor: Rational) -> "PiecewiseConstant":
        factor = _rational(factor, "factor")
        return PiecewiseConstant((left, right, height * factor)
                                 for left, right, height in self.intervals)


__all__ = ["DiracMeasure", "PiecewiseConstant"]
