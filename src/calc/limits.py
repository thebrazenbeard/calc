"""Bounded real-limit examples with a strict proof/diagnostic distinction.

``polynomial_limit_certificate`` returns a rational epsilon-delta witness for
a *given polynomial*, justified by a finite local Lipschitz bound. This is
an exact algebraic proof method for this class, not sampling-based evidence.

``one_sided_samples`` returns only finite binary64 values from approaching
either side. No finite collection of samples proves convergence, divergence,
continuity, or the existence of a limit for an arbitrary function.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from fractions import Fraction
from itertools import islice
import math
from numbers import Real

from .algebra import Polynomial


Exact = int | Fraction
_MAX_POLYNOMIAL_DEGREE = 32
_MAX_RESULT_BITS = 8192
_MAX_SAMPLES = 64


def _checked(value: Fraction, name: str) -> Fraction:
    if max(value.numerator.bit_length(), value.denominator.bit_length()) > _MAX_RESULT_BITS:
        raise ValueError(f"{name} exceeds exact rational resource limit")
    return value


def _exact(value: Exact, name: str) -> Fraction:
    if type(value) is int or isinstance(value, Fraction):
        return _checked(Fraction(value), name)
    raise TypeError(f"{name} must be an exact int or Fraction, not float/bool")


def _polynomial(coefficients: Polynomial | Iterable[Exact]) -> Polynomial:
    if isinstance(coefficients, Polynomial):
        raw = coefficients.coefficients
    else:
        if isinstance(coefficients, (str, bytes)):
            raise TypeError("polynomial coefficients must be numeric")
        try:
            raw = tuple(islice(iter(coefficients), _MAX_POLYNOMIAL_DEGREE + 2))
        except TypeError as exc:
            raise TypeError("polynomial coefficients must be iterable") from exc
    if len(raw) > _MAX_POLYNOMIAL_DEGREE + 1:
        raise ValueError("polynomial exceeds maximum degree 32")
    return Polynomial(_exact(value, "polynomial coefficient") for value in raw)


@dataclass(frozen=True, slots=True)
class EpsilonDeltaCertificate:
    """Exact coefficient-derived witness for ``lim_(x→point) P(x)=P(point)``.

    The included ``lipschitz_bound`` proves
    ``|P(x)-P(point)| <= lipschitz_bound * |x-point|`` whenever
    ``|x-point| <= radius``. Therefore every ``0 < |x-point| < delta``
    implies ``|P(x)-value| < epsilon``. The bound is conservative.
    """

    coefficients: tuple[Fraction, ...]
    point: Fraction
    value: Fraction
    epsilon: Fraction
    delta: Fraction
    radius: Fraction
    lipschitz_bound: Fraction

    def error_bound_at(self, x: Exact) -> Fraction:
        """Upper bound valid only within this certificate's delta-ball.

        Evaluating this method at a rational witness does *not* create a
        second proof. The certificate is justified by its polynomial formula.
        """
        distance = _checked(abs(_exact(x, "x") - self.point), "distance")
        if distance >= self.delta:
            raise ValueError("point lies outside the open epsilon-delta witness")
        return _checked(self.lipschitz_bound * distance, "error upper bound")


def polynomial_limit_certificate(
    coefficients: Polynomial | Iterable[Exact],
    point: Exact,
    epsilon: Exact,
    *,
    radius: Exact = 1,
) -> EpsilonDeltaCertificate:
    """Construct an exact rational epsilon-delta witness for a polynomial.

    If ``P(x)=sum(a_k*x**k)`` and ``|x-c|<=r``, then ``|x|<=B=|c|+r``.
    The factorization of ``x**k-c**k`` proves

        ``|P(x)-P(c)| <= M*|x-c|``

    with ``M=sum_{k>=1} k*|a_k|*B**(k-1)``. Choose
    ``delta=min(r, epsilon/(2*M))`` for ``M>0`` or ``delta=r`` otherwise.
    This exact argument establishes continuity of rational polynomials at
    rational centers; it is not a general arbitrary-function limit solver.
    """
    polynomial = _polynomial(coefficients)
    center = _exact(point, "point")
    eps = _exact(epsilon, "epsilon")
    neighborhood = _exact(radius, "radius")
    if eps <= 0 or neighborhood <= 0:
        raise ValueError("epsilon and radius must be strictly positive")

    value = Fraction(0)
    for coefficient in reversed(polynomial.coefficients):
        value = _checked(value * center + coefficient, "polynomial value")

    magnitude = _checked(abs(center) + neighborhood, "local interval bound")
    monomial_power = Fraction(1)
    bound = Fraction(0)
    for degree, coefficient in enumerate(polynomial.coefficients[1:], start=1):
        term = _checked(degree * abs(coefficient) * monomial_power,
                        "Lipschitz coefficient")
        bound = _checked(bound + term, "local Lipschitz bound")
        monomial_power = _checked(monomial_power * magnitude,
                                  "local polynomial power")
    delta = neighborhood if bound == 0 else min(
        neighborhood, _checked(eps / _checked(2 * bound, "scaled Lipschitz bound"),
                               "epsilon-delta value")
    )
    return EpsilonDeltaCertificate(
        polynomial.coefficients, center, value, eps, delta, neighborhood, bound
    )


def _finite_real(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real number")
    try:
        number = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} is not representable in binary64") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


@dataclass(frozen=True, slots=True)
class OneSidedSamples:
    """Finite samples at ``point +/- offsets``; not a convergence verdict."""

    point: float
    offsets: tuple[float, ...]
    left: tuple[float, ...]
    right: tuple[float, ...]


def one_sided_samples(
    function: Callable[[float], Real],
    point: Real,
    *,
    initial_step: Real = 1.0,
    count: int = 8,
) -> OneSidedSamples:
    """Return finite geometric one-sided numerical samples, without inference.

    The callback is never evaluated at ``point``. If a perturbation rounds
    back to ``point`` or leaves the finite float range, sampling rejects it
    instead of silently claiming to have approached a new point.
    """
    if not callable(function):
        raise TypeError("function must be callable")
    if type(count) is not int or not 1 <= count <= _MAX_SAMPLES:
        raise ValueError(f"count must be an integer between 1 and {_MAX_SAMPLES}")
    center = _finite_real(point, "point")
    step = _finite_real(initial_step, "initial_step")
    if step <= 0:
        raise ValueError("initial_step must be positive")

    offsets: list[float] = []
    left: list[float] = []
    right: list[float] = []
    for i in range(count):
        h = step / (2 ** i)
        left_point, right_point = center - h, center + h
        if (h == 0.0 or left_point == center or right_point == center
                or not math.isfinite(left_point) or not math.isfinite(right_point)):
            raise ValueError("numerical perturbation is not representable")
        offsets.append(h)
        left.append(_finite_real(function(left_point), "left function output"))
        right.append(_finite_real(function(right_point), "right function output"))
    return OneSidedSamples(center, tuple(offsets), tuple(left), tuple(right))


__all__ = [
    "EpsilonDeltaCertificate", "OneSidedSamples", "polynomial_limit_certificate",
    "one_sided_samples",
]
