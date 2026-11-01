"""Small exact affine-algebraic sets over Q, with explicit finite boundaries.

The zero locus of a polynomial system is represented by its equations, not by
enumerating every solution. ``points_on_grid`` checks only caller-supplied
rational coordinates and cannot certify that a variety has no other points.

Ideal membership delegates to the package's bounded Buchberger implementation.
The result is exact if it returns; exhaustion raises ResourceLimitError rather
than silently returning an incomplete answer. Neither a full algebraic closure
nor general algebraic-geometry theorem proving is implemented here.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from fractions import Fraction
from itertools import islice, product
from math import prod

from .grobner import (
    MultiPolynomial,
    ResourceLimitError,
    grobner_basis,
    normal_form,
)

Exact = int | Fraction
Point = tuple[Fraction, ...]

_MAX_VARS = 4
_MAX_EQUATIONS = 16
_MAX_DEGREE = 32
_MAX_TOTAL_TERMS = 1024
_MAX_AXIS_VALUES = 64
_MAX_GRID_POINTS = 10_000
_MAX_COORDINATE_BITS = 512
_MAX_RESULT_BITS = 32_768


def _fraction(value: Exact, name: str) -> Fraction:
    if type(value) is not int and not isinstance(value, Fraction):
        raise TypeError(f"{name} must be int or Fraction, never float/bool")
    result = Fraction(value)
    if max(result.numerator.bit_length(), result.denominator.bit_length()) > _MAX_COORDINATE_BITS:
        raise ResourceLimitError(f"{name} exceeds rational coordinate budget")
    return result


def _bounded_result(value: Fraction) -> Fraction:
    if max(value.numerator.bit_length(), value.denominator.bit_length()) > _MAX_RESULT_BITS:
        raise ResourceLimitError("evaluation exceeds exact arithmetic result budget")
    return value


def _point(values: Sequence[Exact], nvars: int) -> Point:
    if not isinstance(values, (tuple, list)) or len(values) != nvars:
        raise ValueError(f"point needs exactly {nvars} rational coordinates")
    return tuple(_fraction(value, "point coordinate") for value in values)


def _valid_polynomial(polynomial: object, *, nvars: int) -> MultiPolynomial:
    if not isinstance(polynomial, MultiPolynomial) or polynomial.nvars != nvars:
        raise ValueError("equations must be MultiPolynomial values in the same variables")
    if any(sum(exponents) > _MAX_DEGREE for exponents, _ in polynomial.terms):
        raise ResourceLimitError("total polynomial degree exceeds affine-geometry budget")
    return polynomial


def partial_derivative(polynomial: MultiPolynomial, variable: int) -> MultiPolynomial:
    """Exact formal partial derivative of a Q[x0,...]-polynomial."""
    if not isinstance(polynomial, MultiPolynomial):
        raise TypeError("polynomial must be a MultiPolynomial")
    if type(variable) is not int or not 0 <= variable < polynomial.nvars:
        raise ValueError("variable must be a zero-based coordinate index")
    terms: dict[tuple[int, ...], Fraction] = {}
    for exponents, coefficient in polynomial.terms:
        power = exponents[variable]
        if power:
            lower = list(exponents)
            lower[variable] -= 1
            terms[tuple(lower)] = coefficient * power
    return MultiPolynomial(terms, nvars=polynomial.nvars)


@dataclass(frozen=True, slots=True, init=False)
class AffineSystem:
    """A bounded finite list of rational polynomial equations ``f_i=0``.

    The equations define a set theoretically, over whatever extension of Q
    one later chooses. Only exact rational point membership, finite rational
    grid searches, Jacobians at rational points, and ideal membership over Q
    are calculated by this class.
    """

    equations: tuple[MultiPolynomial, ...]
    nvars: int

    def __init__(self, equations: Iterable[MultiPolynomial]):
        if isinstance(equations, (str, bytes)):
            raise TypeError("equations must be an iterable of polynomials")
        try:
            items = tuple(islice(iter(equations), _MAX_EQUATIONS + 1))
        except TypeError as exc:
            raise TypeError("equations must be an iterable of polynomials") from exc
        if not items:
            raise ValueError("at least one equation is needed to define variable dimension")
        if len(items) > _MAX_EQUATIONS:
            raise ResourceLimitError("too many defining equations")
        first = items[0]
        if not isinstance(first, MultiPolynomial):
            raise TypeError("equations must be exact MultiPolynomial objects")
        if first.nvars > _MAX_VARS:
            raise ResourceLimitError("ambient affine dimension exceeds budget")
        parsed = tuple(_valid_polynomial(poly, nvars=first.nvars) for poly in items)
        if sum(len(poly.terms) for poly in parsed) > _MAX_TOTAL_TERMS:
            raise ResourceLimitError("system term budget exceeded")
        object.__setattr__(self, "equations", parsed)
        object.__setattr__(self, "nvars", first.nvars)

    def contains(self, values: Sequence[Exact]) -> bool:
        """Whether all equations vanish at one exact rational point."""
        point = _point(values, self.nvars)
        for poly in self.equations:
            if _bounded_result(poly.evaluate(point)):
                return False
        return True

    def jacobian_at(self, values: Sequence[Exact]) -> tuple[tuple[Fraction, ...], ...]:
        """Exact Jacobian of the defining equations at a rational point.

        This is a derivative matrix, not by itself a general smoothness or
        dimension certification for the associated variety.
        """
        point = _point(values, self.nvars)
        return tuple(tuple(_bounded_result(partial_derivative(poly, index).evaluate(point))
                           for index in range(self.nvars)) for poly in self.equations)

    def points_on_grid(self, axes: Sequence[Sequence[Exact]]) -> tuple[Point, ...]:
        """Enumerate *only* rational grid points satisfying all equations.

        Axes are discrete candidate coordinates.  An empty result does not
        imply the true affine zero set is empty, even over Q.
        """
        if not isinstance(axes, (list, tuple)) or len(axes) != self.nvars:
            raise ValueError(f"grid needs exactly {self.nvars} coordinate axes")
        grid: list[tuple[Fraction, ...]] = []
        for raw_axis in axes:
            if not isinstance(raw_axis, (list, tuple)):
                raise TypeError("grid axes must be finite list/tuple sequences")
            if len(raw_axis) > _MAX_AXIS_VALUES:
                raise ResourceLimitError("grid axis exceeds candidate limit")
            axis = tuple(dict.fromkeys(_fraction(value, "grid coordinate") for value in raw_axis))
            grid.append(axis)
        if prod(len(axis) for axis in grid) > _MAX_GRID_POINTS:
            raise ResourceLimitError("rational grid search exceeds point budget")
        return tuple(point for point in product(*grid) if self.contains(point))

    def ideal_contains(self, polynomial: MultiPolynomial, *, order: str = "lex") -> bool:
        """Exact ideal-membership decision using a completed Gröbner basis.

        This answers ``polynomial in (f_1,...,f_k)`` over Q[x], not whether the
        polynomial vanishes on the same points: a radical is a distinct ideal.
        """
        candidate = _valid_polynomial(polynomial, nvars=self.nvars)
        if order not in ("lex", "grlex"):
            raise ValueError("monomial order must be 'lex' or 'grlex'")
        if not candidate:
            return True
        basis_generators = tuple(poly for poly in self.equations if poly)
        if not basis_generators:
            return False
        basis = grobner_basis(basis_generators, order, max_pairs=500,
                              max_steps=5000, max_terms=1024, max_basis=32)
        remainder = normal_form(candidate, basis, order,
                                max_steps=5000, max_terms=1024)
        return not remainder


__all__ = ["AffineSystem", "partial_derivative"]
