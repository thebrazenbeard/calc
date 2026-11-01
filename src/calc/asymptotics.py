"""Exact leading-term asymptotics for bounded rational functions over Q.

This is a *restricted symbolic* implementation.  A nonzero rational function
P(x)/Q(x) has the leading term ``c*x**k`` at positive infinity, where
``k = degree(P)-degree(Q)`` and ``c = lead(P)/lead(Q)``.  At ``0+`` the
lowest nonzero powers play the same role.  This yields exact O, o and Theta
comparisons without pretending finite samples prove an asymptotic relation.

No claims are made about arbitrary functions, essential singularities,
piecewise expressions, infinite series, or numerical asymptotic fitting.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction
from itertools import islice
from typing import Literal

from .algebra import Polynomial


Exact = int | Fraction
Approach = Literal["infinity", "zero_plus"]
_MAX_DEGREE = 64
_MAX_BITS = 4096


def _as_polynomial(value: Polynomial | Iterable[Exact], name: str) -> Polynomial:
    if isinstance(value, Polynomial):
        coefficients = value.coefficients
    else:
        if isinstance(value, (str, bytes)):
            raise TypeError(f"{name} must contain exact rational coefficients")
        try:
            coefficients = tuple(islice(iter(value), _MAX_DEGREE + 2))
        except TypeError as exc:
            raise TypeError(f"{name} must be a coefficient iterable") from exc
    if len(coefficients) > _MAX_DEGREE + 1:
        raise ValueError(f"{name} degree exceeds {_MAX_DEGREE}")
    checked: list[Fraction] = []
    for coefficient in coefficients:
        if type(coefficient) is not int and not isinstance(coefficient, Fraction):
            raise TypeError(f"{name} coefficients must be int or Fraction")
        rational = Fraction(coefficient)
        if (rational.numerator.bit_length() > _MAX_BITS or
                rational.denominator.bit_length() > _MAX_BITS):
            raise ValueError(f"{name} coefficient exceeds {_MAX_BITS}-bit limit")
        checked.append(rational)
    return Polynomial(checked)


@dataclass(frozen=True, slots=True, init=False)
class RationalAsymptotic:
    """A rational function ``numerator / denominator`` in the variable x.

    Inputs are exact rational polynomials with low-to-high coefficients.
    The denominator cannot be identically zero.  There may be finite poles;
    the asymptotic statements apply on a punctured eventual domain.
    """

    numerator: Polynomial
    denominator: Polynomial

    def __init__(self, numerator: Polynomial | Iterable[Exact],
                 denominator: Polynomial | Iterable[Exact] = (1,)) -> None:
        p = _as_polynomial(numerator, "numerator")
        q = _as_polynomial(denominator, "denominator")
        if q.degree < 0:
            raise ZeroDivisionError("rational-function denominator cannot be zero")
        object.__setattr__(self, "numerator", p)
        object.__setattr__(self, "denominator", q)

    @property
    def is_zero(self) -> bool:
        return self.numerator.degree < 0

    def leading_term(self, *, at: Approach = "infinity") -> tuple[int, Fraction]:
        """Return (power, coefficient) for ``f(x) ~ coefficient*x**power``.

        At ``infinity``, x tends to positive infinity; at ``zero_plus`` it
        tends to zero from the positive side. Zero has no nonzero leading term.
        """
        if at not in ("infinity", "zero_plus"):
            raise ValueError("at must be 'infinity' or 'zero_plus'")
        if self.is_zero:
            raise ValueError("zero rational function has no nonzero leading term")
        p, q = self.numerator.coefficients, self.denominator.coefficients
        if at == "infinity":
            return (self.numerator.degree - self.denominator.degree,
                    p[-1] / q[-1])
        p_order = next(i for i, c in enumerate(p) if c)
        q_order = next(i for i, c in enumerate(q) if c)
        return (p_order - q_order, p[p_order] / q[q_order])


@dataclass(frozen=True, slots=True)
class GrowthComparison:
    """Exact order relations for f relative to an eventually nonzero g.

    ``ratio_limit=None`` means f/g has no finite limit in Q because its
    leading power diverges; otherwise the limit is an exact Fraction.
    """

    approach: Approach
    big_o: bool
    little_o: bool
    theta: bool
    asymptotically_equivalent: bool
    ratio_limit: Fraction | None
    ratio_power: int | None
    ratio_leading_coefficient: Fraction | None


def compare_growth(
    first: RationalAsymptotic,
    second: RationalAsymptotic,
    *,
    at: Approach = "infinity",
) -> GrowthComparison:
    """Decide restricted O, o, Theta and asymptotic equality exactly.

    The comparator ``second`` must be nonzero, so f/g is defined on a
    sufficiently distant (or punctured near-zero) neighborhood.  ``Theta``
    means the ratio tends to a finite *nonzero* constant for this class.
    No numerical evaluation or heuristic sampling is used.
    """
    if not isinstance(first, RationalAsymptotic) or not isinstance(second, RationalAsymptotic):
        raise TypeError("compare_growth requires two RationalAsymptotic objects")
    if at not in ("infinity", "zero_plus"):
        raise ValueError("at must be 'infinity' or 'zero_plus'")
    if second.is_zero:
        raise ValueError("comparison reference must not be identically zero")
    if first.is_zero:
        return GrowthComparison(at, True, True, False, False, Fraction(0), None, None)

    first_power, first_coefficient = first.leading_term(at=at)
    second_power, second_coefficient = second.leading_term(at=at)
    power = first_power - second_power
    coefficient = first_coefficient / second_coefficient
    # At 0+, positive powers vanish; at +infinity, negative powers vanish.
    decays = power < 0 if at == "infinity" else power > 0
    diverges = power > 0 if at == "infinity" else power < 0
    limit = None if diverges else (Fraction(0) if decays else coefficient)
    return GrowthComparison(
        at, not diverges, decays, power == 0,
        power == 0 and coefficient == 1, limit, power, coefficient,
    )


__all__ = ["RationalAsymptotic", "GrowthComparison", "compare_growth"]
