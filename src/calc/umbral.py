"""Exact rational examples of classical umbral / finite operator calculus.

Polynomials use ascending monomial coefficients, e.g. (1, 2, 3) means
``1 + 2*x + 3*x**2``. Falling-factorial coefficient tuples instead mean
``sum(c[k] * x_(k))`` where ``x_(k)=x(x-1)...(x-k+1)``. Both domains are
bounded to degree at most 32 and use Fraction, never implicit float coercion.

The operators here use the ordinary shift E_h and forward difference
``Delta_h = E_h - I``. This is *classical* umbral calculus, not a generic
Sheffer, q-difference, or symbolic expression engine.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import islice
from math import comb
from collections.abc import Iterable


Exact = int | Fraction
_MAX_DEGREE = 32
_MAX_BITS = 4096


def _fraction(value: Exact, name: str) -> Fraction:
    if type(value) is int or isinstance(value, Fraction):
        exact = Fraction(value)
        if max(exact.numerator.bit_length(), exact.denominator.bit_length()) > _MAX_BITS:
            raise ValueError(f"{name} exceeds rational bit limit")
        return exact
    raise TypeError(f"{name} must be int or Fraction, never float/bool")


def _order(order: int, name: str = "order") -> int:
    if type(order) is not int:
        raise TypeError(f"{name} must be an integer, not bool or Fraction")
    if not 0 <= order <= _MAX_DEGREE:
        raise ValueError(f"{name} must be between 0 and {_MAX_DEGREE}")
    return order


def _coefficients(coefficients: Iterable[Exact], *, empty_ok: bool = True) -> tuple[Fraction, ...]:
    if isinstance(coefficients, (str, bytes)):
        raise TypeError("polynomial coefficients must be iterable rationals")
    try:
        iterator = iter(coefficients)
    except TypeError as exc:
        raise TypeError("polynomial coefficients must be iterable rationals") from exc
    raw = tuple(islice(iterator, _MAX_DEGREE + 2))
    if len(raw) > _MAX_DEGREE + 1:
        raise ValueError("polynomial degree exceeds 32")
    if not raw and not empty_ok:
        raise ValueError("Newton interpolation requires at least one value")
    exact = [_fraction(value, "coefficient") for value in raw]
    while len(exact) > 1 and not exact[-1]:
        exact.pop()
    return tuple(exact) if exact else (Fraction(0),)


def _trim(coefficients: list[Fraction]) -> tuple[Fraction, ...]:
    while len(coefficients) > 1 and not coefficients[-1]:
        coefficients.pop()
    return tuple(coefficients) if coefficients else (Fraction(0),)


def _multiply_linear(coefficients: tuple[Fraction, ...], root: Fraction) -> tuple[Fraction, ...]:
    result = [Fraction(0)] * (len(coefficients) + 1)
    for j, coefficient in enumerate(coefficients):
        result[j] -= root * coefficient
        result[j + 1] += coefficient
    return tuple(result)


def falling_factorial(x: Exact, order: int) -> Fraction:
    """Compute ``x_(n) = x(x-1)...(x-n+1)`` exactly, with ``x_(0)=1``."""
    argument, n = _fraction(x, "x"), _order(order)
    result = Fraction(1)
    for k in range(n):
        result *= argument - k
    return result


def rising_factorial(x: Exact, order: int) -> Fraction:
    """Compute ``x^(n) = x(x+1)...(x+n-1)`` exactly, with ``x^(0)=1``."""
    argument, n = _fraction(x, "x"), _order(order)
    result = Fraction(1)
    for k in range(n):
        result *= argument + k
    return result


@lru_cache(maxsize=64)
def _stirling_first_row(n: int) -> tuple[int, ...]:
    """Signed Stirling first-kind row s(n,k) via the factorial recurrence."""
    if n == 0:
        return (1,)
    previous = _stirling_first_row(n - 1)
    return tuple((previous[k - 1] if k else 0) -
                 (n - 1) * (previous[k] if k < n else 0)
                 for k in range(n + 1))


@lru_cache(maxsize=64)
def _stirling_second_row(n: int) -> tuple[int, ...]:
    """Stirling second-kind row S(n,k) via the partition recurrence."""
    if n == 0:
        return (1,)
    previous = _stirling_second_row(n - 1)
    return tuple((previous[k - 1] if k else 0) +
                 k * (previous[k] if k < n else 0)
                 for k in range(n + 1))


def stirling_first(n: int, k: int) -> int:
    """Return signed ``s(n,k)`` satisfying ``x_(n)=Σ s(n,k)*x**k``."""
    size = _order(n, "n")
    column = _order(k, "k")
    return _stirling_first_row(size)[column] if column <= size else 0


def stirling_second(n: int, k: int) -> int:
    """Return ``S(n,k)`` satisfying ``x**n=Σ S(n,k)*x_(k)``."""
    size = _order(n, "n")
    column = _order(k, "k")
    return _stirling_second_row(size)[column] if column <= size else 0


def monomial_to_falling(coefficients: Iterable[Exact]) -> tuple[Fraction, ...]:
    """Convert ascending power basis coefficients to falling-factorial basis."""
    monomial = _coefficients(coefficients)
    result = [Fraction(0)] * len(monomial)
    for n, coefficient in enumerate(monomial):
        if coefficient:
            row = _stirling_second_row(n)
            for k in range(n + 1):
                result[k] += coefficient * row[k]
    return _trim(result)


def falling_to_monomial(coefficients: Iterable[Exact]) -> tuple[Fraction, ...]:
    """Convert ascending falling-factorial basis coefficients to powers."""
    falling = _coefficients(coefficients)
    result = [Fraction(0)] * len(falling)
    for n, coefficient in enumerate(falling):
        if coefficient:
            row = _stirling_first_row(n)
            for k in range(n + 1):
                result[k] += coefficient * row[k]
    return _trim(result)


def forward_difference(coefficients: Iterable[Exact], *, step: Exact = 1,
                       order: int = 1) -> tuple[Fraction, ...]:
    """Return coefficients of ``Delta_h**order p = (E_h-I)**order p``.

    This operator is a raw difference, **not** the divided-difference
    ``(p(x+h)-p(x))/h``. The nonzero rational step may be negative.
    """
    polynomial = _coefficients(coefficients)
    h, count = _fraction(step, "step"), _order(order)
    if not h:
        raise ValueError("forward-difference step must be nonzero")
    for _ in range(count):
        if len(polynomial) == 1:
            polynomial = (Fraction(0),)
            continue
        output = [Fraction(0)] * (len(polynomial) - 1)
        for degree in range(1, len(polynomial)):
            coefficient = polynomial[degree]
            for power in range(degree):
                output[power] += coefficient * comb(degree, power) * h ** (degree - power)
        polynomial = _trim(output)
    return polynomial


def newton_forward_interpolate(
    values: Iterable[Exact], *, start: Exact = 0, step: Exact = 1,
) -> tuple[Fraction, ...]:
    """Interpolate data at ``start+j*step`` into ascending power coefficients.

    The polynomial of degree at most ``len(values)-1`` is constructed from
    Newton's exact forward-difference formula. Values at distinct nodes are
    necessary: ``step`` must be nonzero. At most 33 samples are accepted.
    """
    if isinstance(values, (str, bytes)):
        raise TypeError("values must be iterable rationals")
    try:
        raw = tuple(islice(iter(values), _MAX_DEGREE + 2))
    except TypeError as exc:
        raise TypeError("values must be iterable rationals") from exc
    if not raw:
        raise ValueError("Newton interpolation requires at least one value")
    if len(raw) > _MAX_DEGREE + 1:
        raise ValueError("Newton interpolation exceeds 33 nodes")
    data = tuple(_fraction(v, "sample value") for v in raw)
    first, h = _fraction(start, "start"), _fraction(step, "step")
    if not h:
        raise ValueError("Newton interpolation step must be nonzero")
    result = [Fraction(0)] * len(data)
    difference_row = data
    basis = (Fraction(1),)
    denominator = Fraction(1)
    for degree in range(len(data)):
        coefficient = difference_row[0] / denominator
        for j, term in enumerate(basis):
            result[j] += coefficient * term
        if degree == len(data) - 1:
            break
        difference_row = tuple(right - left for left, right in
                               zip(difference_row, difference_row[1:]))
        basis = _multiply_linear(basis, first + degree * h)
        denominator *= (degree + 1) * h
    return _trim(result)


__all__ = [
    "falling_factorial", "rising_factorial", "stirling_first", "stirling_second",
    "forward_difference", "monomial_to_falling", "falling_to_monomial",
    "newton_forward_interpolate",
]
