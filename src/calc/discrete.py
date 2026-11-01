"""Exact finite-sequence calculus over rational arithmetic.

These are difference and summation operators on an *already sampled uniform*
grid. They do not interpolate continuous functions, represent accumulation
points, or implement general dynamic equations on time scales. For callable
finite-difference approximations and finite nonuniform time scales, see the
separate ``advanced`` and ``time_scales`` modules.
"""

from __future__ import annotations

from collections.abc import Iterable
from fractions import Fraction
from math import comb


Exact = int | Fraction
_MAX_SAMPLES = 10_000
_MAX_ORDER = 32
_MAX_BITS = 1024
_MAX_POWER = 10
_MAX_INDEX = 10**12


def _fraction(value: Exact, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(f"{name} must be an integer or Fraction; floats are inexact")
    answer = Fraction(value)
    if (abs(answer.numerator).bit_length() > _MAX_BITS
            or answer.denominator.bit_length() > _MAX_BITS):
        raise ValueError(f"{name} exceeds exact coefficient size limit")
    return answer


def _sequence(values: Iterable[Exact]) -> tuple[Fraction, ...]:
    if isinstance(values, (str, bytes)):
        raise TypeError("samples must be an iterable of exact rational values")
    try:
        iterator = iter(values)
    except TypeError as exc:
        raise TypeError("samples must be iterable") from exc
    samples: list[Fraction] = []
    for raw in iterator:
        if len(samples) >= _MAX_SAMPLES:
            raise ValueError("sample sequence exceeds 10000-value resource cap")
        samples.append(_fraction(raw, "sample"))
    return tuple(samples)


def _order(order: int) -> int:
    if type(order) is not int:
        raise TypeError("difference order must be an integer")
    if order < 0 or order > _MAX_ORDER:
        raise ValueError("difference order must be between 0 and 32")
    return order


def _step(step: Exact) -> Fraction:
    h = _fraction(step, "step")
    if not h:
        raise ValueError("sample step must be nonzero")
    return h


def forward_differences(
    values: Iterable[Exact], *, order: int = 1, step: Exact = 1,
) -> tuple[Fraction, ...]:
    """Normalized forward differences on a uniform sequence.

    For samples ``f(x_0 + k*h)``, returns ``Delta_h^order f / h^order``
    at the nodes ``k=0,...,N-order-1``. ``order=0`` returns the samples;
    orders beyond the sample count have no valid stencil and return ``()``.
    This is an exact discrete operator, **not** a derivative limit.
    """
    current = _sequence(values)
    count, h = _order(order), _step(step)
    for _ in range(count):
        current = tuple((right-left)/h for left, right in zip(current, current[1:]))
    return current


def backward_differences(
    values: Iterable[Exact], *, order: int = 1, step: Exact = 1,
) -> tuple[Fraction, ...]:
    """Normalized backward differences aligned to ending nodes.

    The numerical stencil values match the forward differences, but the
    locations are ``k=order,...,N-1`` rather than their starting nodes.
    """
    return forward_differences(values, order=order, step=step)


def central_differences(
    values: Iterable[Exact], *, step: Exact = 1,
) -> tuple[Fraction, ...]:
    """First centered difference ``(f_{k+1}-f_{k-1})/(2h)``.

    Returns values at the interior sample positions. Fewer than three
    samples yield an empty result. A centered difference is generally an
    approximation to the continuous derivative, not an equality.
    """
    samples, h = _sequence(values), _step(step)
    return tuple((samples[i + 2]-samples[i])/(2*h)
                 for i in range(max(0, len(samples)-2)))


def cumulative_sum(values: Iterable[Exact], *, initial: Exact = 0) -> tuple[Fraction, ...]:
    """Discrete antiderivative with an explicit integration constant.

    For ``F_0=initial``, returns ``F_0,...,F_N`` where
    ``F_{k+1}-F_k = values[k]``. There are ``N+1`` output values.
    """
    samples, total = _sequence(values), _fraction(initial, "initial")
    result = [total]
    for value in samples:
        total += value
        result.append(total)
    return tuple(result)


def discrete_integral(
    values: Iterable[Exact], start: int = 0, stop: int | None = None,
) -> Fraction:
    """Exact signed sum over a half-open index interval.

    For ``start <= stop``, computes ``sum(values[k], k=start..stop-1)``.
    Reversed bounds negate the same sum. Both endpoints must be integers
    between ``0`` and the sequence length, inclusive. A missing ``stop``
    means the end of the sequence.
    """
    samples = _sequence(values)
    end = len(samples) if stop is None else stop
    if type(start) is not int or type(end) is not int:
        raise TypeError("summation bounds must be integers")
    if not 0 <= start <= len(samples) or not 0 <= end <= len(samples):
        raise IndexError("summation bounds must lie within the sampled grid")
    if start <= end:
        return sum(samples[start:end], Fraction(0))
    return -sum(samples[end:start], Fraction(0))


def summation_by_parts(
    u: Iterable[Exact], v: Iterable[Exact],
) -> tuple[Fraction, Fraction, Fraction]:
    """Return the three terms of discrete summation by parts.

    For equal-length endpoint sequences ``u_0..u_n``, ``v_0..v_n``:

    ``lhs = sum_{k=0}^{n-1} u_k (v_{k+1}-v_k)``
    ``boundary = u_n v_n - u_0 v_0``
    ``correction = sum_{k=0}^{n-1} v_{k+1}(u_{k+1}-u_k)``

    Exactly ``lhs == boundary - correction`` for rational inputs. A
    one-point pair represents the empty sum with zero boundary difference.
    """
    left, right = _sequence(u), _sequence(v)
    if not left or len(left) != len(right):
        raise ValueError("summation by parts needs equally sized nonempty endpoint arrays")
    lhs = sum((a*(b1-b0) for a,b0,b1 in zip(left, right, right[1:])), Fraction(0))
    boundary = left[-1]*right[-1] - left[0]*right[0]
    correction = sum((b*(a1-a0) for a0,a1,b in zip(left, left[1:], right[1:])), Fraction(0))
    return lhs, boundary, correction


_BERNOULLI: tuple[Fraction, ...] = (
    Fraction(1), Fraction(-1, 2), Fraction(1, 6), Fraction(0),
    Fraction(-1, 30), Fraction(0), Fraction(1, 42), Fraction(0),
    Fraction(-1, 30), Fraction(0), Fraction(5, 66),
)


def power_sum(n: int, power: int) -> int:
    """Exact Faulhaber sum ``0**power+...+(n-1)**power``.

    Supports nonnegative ``n<=10**12`` and integer ``0<=power<=10``.
    The empty sum is zero; for ``power=0`` each included term counts as one,
    including the conventional ``0**0=1`` combinatorial term. The result is
    evaluated from Bernoulli-number coefficients (with ``B_1=-1/2``), not
    from a potentially enormous loop over samples.
    """
    if type(n) is not int or type(power) is not int:
        raise TypeError("n and power must be integers")
    if n < 0 or n > _MAX_INDEX or power < 0 or power > _MAX_POWER:
        raise ValueError("requires 0<=n<=10**12 and 0<=power<=10")
    total = sum(
        (comb(power+1, j) * _BERNOULLI[j] * n**(power+1-j)
         for j in range(power+1)),
        Fraction(0),
    ) / (power+1)
    if total.denominator != 1:
        raise ArithmeticError("Faulhaber rational coefficients did not yield an integer")
    return total.numerator
