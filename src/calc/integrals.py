"""Exact, bounded integration identities for rational univariate polynomials.

This is a restricted symbolic integrator over ``Q[x]``. It deliberately does
not integrate arbitrary elementary functions, discontinuous integrands, or
improper integrals; numerical quadrature lives in :mod:`calc.numerical`.

Polynomials use :class:`calc.algebra.Polynomial`, whose coefficients are in
ascending powers. Endpoints and constants must be integers or Fractions.
"""

from __future__ import annotations

from fractions import Fraction

from calc.algebra import Polynomial

Exact = int | Fraction
_MAX_DEGREE = 64
_MAX_INPUT_BITS = 4096
_MAX_OUTPUT_BITS = 65536


def _checked_number(value: Exact, name: str) -> Fraction:
    if type(value) is not int and not isinstance(value, Fraction):
        raise TypeError(f"{name} must be an exact integer or Fraction")
    exact = Fraction(value)
    if max(exact.numerator.bit_length(), exact.denominator.bit_length()) > _MAX_INPUT_BITS:
        raise ValueError(f"{name} exceeds the exact input bit limit")
    return exact


def _result(value: Fraction, name: str) -> Fraction:
    if max(value.numerator.bit_length(), value.denominator.bit_length()) > _MAX_OUTPUT_BITS:
        raise ValueError(f"{name} exceeds the exact output bit limit")
    return value


def _checked_polynomial(polynomial: Polynomial, name: str) -> Polynomial:
    if not isinstance(polynomial, Polynomial):
        raise TypeError(f"{name} must be a calc.algebra.Polynomial")
    if polynomial.degree > _MAX_DEGREE:
        raise ValueError(f"{name} exceeds the {_MAX_DEGREE}-degree limit")
    for coefficient in polynomial.coefficients:
        if max(coefficient.numerator.bit_length(), coefficient.denominator.bit_length()) > _MAX_INPUT_BITS:
            raise ValueError(f"{name} exceeds the exact coefficient bit limit")
    return polynomial


def _evaluate(polynomial: Polynomial, argument: Fraction, name: str) -> Fraction:
    degree = polynomial.degree
    if degree >= 0:
        coefficient_bits = max(
            max(c.numerator.bit_length(), c.denominator.bit_length())
            for c in polynomial.coefficients
        )
        argument_bits = max(1, argument.numerator.bit_length(), argument.denominator.bit_length())
        # A conservative cap before raising a large rational to a high power.
        if coefficient_bits + degree * argument_bits > _MAX_OUTPUT_BITS - 64:
            raise ValueError(f"{name} would exceed the rational evaluation budget")
    return _result(polynomial(argument), name)


def polynomial_antiderivative(
    polynomial: Polynomial, *, constant: Exact = 0,
) -> Polynomial:
    """Return the exact polynomial ``F`` with ``F'=polynomial`` and ``F(0)=constant``."""
    p = _checked_polynomial(polynomial, "polynomial")
    c = _checked_number(constant, "constant")
    primitive = p.integrate(c)
    for coefficient in primitive.coefficients:
        _result(coefficient, "antiderivative coefficient")
    return primitive


def definite_integral(polynomial: Polynomial, start: Exact, stop: Exact) -> Fraction:
    """Return exactly ``∫_start^stop p(x) dx`` for ``p∈Q[x]``.

    Reversed bounds negate the result. Equal bounds give zero. Because the
    integrand is polynomial, the Fundamental Theorem of Calculus applies on
    every finite rational interval, without numerical approximation.
    """
    p = _checked_polynomial(polynomial, "polynomial")
    a, b = _checked_number(start, "start"), _checked_number(stop, "stop")
    if a == b:
        return Fraction(0)
    primitive = polynomial_antiderivative(p)
    return _result(_evaluate(primitive, b, "integral endpoint") -
                   _evaluate(primitive, a, "integral endpoint"), "definite integral")


def substitution_integral(
    outer: Polynomial, inner: Polynomial, start: Exact, stop: Exact,
) -> Fraction:
    """Return ``∫ p(q(x)) q'(x) dx`` using exact polynomial substitution.

    ``outer`` supplies ``p`` and ``inner`` supplies ``q``. The change-of-
    variables formula gives ``P(q(stop))-P(q(start))`` for ``P'=p``. This
    avoids symbolically expanding the composite integrand.
    """
    p = _checked_polynomial(outer, "outer")
    q = _checked_polynomial(inner, "inner")
    a, b = _checked_number(start, "start"), _checked_number(stop, "stop")
    if a == b or p.degree < 0 or q.degree <= 0:
        return Fraction(0)
    if max(1, p.degree + 1) * q.degree > _MAX_DEGREE:
        raise ValueError("composite polynomial exceeds supported degree budget")
    qa, qb = _evaluate(q, a, "inner endpoint"), _evaluate(q, b, "inner endpoint")
    primitive = polynomial_antiderivative(p)
    return _result(_evaluate(primitive, qb, "outer endpoint") -
                   _evaluate(primitive, qa, "outer endpoint"), "substitution integral")


def integration_by_parts(
    u: Polynomial, v: Polynomial, start: Exact, stop: Exact,
) -> Fraction:
    """Return ``∫ u'(x)v(x) dx`` via ``[uv] - ∫uv'`` exactly.

    This fixed convention is explicit: ``u`` is *not* differentiated by the
    caller, whereas ``v`` is differentiated in the correction term.
    """
    left = _checked_polynomial(u, "u")
    right = _checked_polynomial(v, "v")
    if max(0, left.degree) + max(0, right.degree) > _MAX_DEGREE:
        raise ValueError("polynomial product exceeds degree budget")
    a, b = _checked_number(start, "start"), _checked_number(stop, "stop")
    if a == b or left.degree < 0 or right.degree < 0:
        return Fraction(0)
    boundary = _result(_evaluate(left, b, "u(b)") * _evaluate(right, b, "v(b)") -
                       _evaluate(left, a, "u(a)") * _evaluate(right, a, "v(a)"),
                       "integration-by-parts boundary")
    correction = definite_integral(left * right.differentiate(), a, b)
    return _result(boundary - correction, "integration by parts")


__all__ = [
    "polynomial_antiderivative", "definite_integral",
    "substitution_integral", "integration_by_parts",
]
