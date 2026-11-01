"""Small, explicit computational examples of advanced calculus.

These routines cover narrowly defined operations, **not** general symbolic
fractional, exterior, discrete, or variational calculus systems. Mathematical
domains and approximation limits are explained in ``docs/ADVANCED_METHODS.md``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from math import comb, exp, isfinite, lgamma, log
from numbers import Real
from typing import TypeVar


Number = TypeVar("Number")


def _finite_real(value: Real, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    return answer


def riemann_liouville_power_integral(alpha: Real, beta: Real, x: Real) -> float:
    """Evaluate ``I_0^alpha t^beta`` at ``x`` using a known closed form.

    The left Riemann–Liouville integral is valid for ``alpha > 0``,
    ``beta > -1`` and ``x > 0``. At ``x = 0`` only finite zero limits are
    returned (requiring ``alpha + beta > 0``). This is *not* a general
    fractional integration routine.

    Formula: ``Gamma(beta + 1) / Gamma(alpha + beta + 1) * x**(alpha+beta)``.
    """
    a = _finite_real(alpha, "alpha")
    b = _finite_real(beta, "beta")
    point = _finite_real(x, "x")
    if a <= 0 or b <= -1 or point < 0:
        raise ValueError("requires alpha > 0, beta > -1, x >= 0")
    exponent = a + b
    if point == 0:
        if exponent <= 0:
            raise ValueError("fractional integral has no zero limit at x=0")
        return 0.0
    return exp(lgamma(b + 1) - lgamma(a + b + 1) + exponent * log(point))


def finite_difference(
    function: Callable[[float], Number],
    x: Real,
    *,
    step: Real = 1,
    order: int = 1,
) -> Number:
    """Return the ``order``-th forward finite difference ``Delta_h^n f(x)/h^n``.

    ``order=0`` returns ``f(x)``; ``order >= 1`` approximates the corresponding
    derivative for smooth functions and a small nonzero step. It is an *exact*
    finite-difference operator for the chosen discrete step.
    """
    point = _finite_real(x, "x")
    h = _finite_real(step, "step")
    if h == 0:
        raise ValueError("step must be nonzero")
    if isinstance(order, bool) or not isinstance(order, int):
        raise TypeError("order must be an integer")
    if order < 0:
        raise ValueError("order must be nonnegative")
    if order == 0:
        return function(point)
    numerator = sum(
        (-1) ** (order - k) * comb(order, k) * function(point + k * h)
        for k in range(order + 1)
    )
    return numerator / h**order


def q_derivative(function: Callable[[float], Number], x: Real, *, q: Real) -> Number:
    """Jackson q-derivative ``(f(q*x)-f(x))/((q-1)*x)``.

    This implementation requires positive ``q != 1`` and nonzero ``x``;
    limiting definitions at ``x=0`` or ``q=1`` are intentionally not inferred.
    """
    point = _finite_real(x, "x")
    ratio = _finite_real(q, "q")
    if point == 0 or ratio <= 0 or ratio == 1:
        raise ValueError("requires x != 0 and q > 0, q != 1")
    return (function(ratio * point) - function(point)) / ((ratio - 1) * point)


def _validate_form(form: Mapping[tuple[int, ...], Number]) -> None:
    if not isinstance(form, Mapping):
        raise TypeError("a form must map basis tuples to coefficients")
    degrees: set[int] = set()
    for basis, coefficient in form.items():
        if not isinstance(basis, tuple) or any(
            isinstance(index, bool) or not isinstance(index, int) or index < 0
            for index in basis
        ):
            raise ValueError("basis must be a tuple of nonnegative integer indices")
        if tuple(sorted(set(basis))) != basis:
            raise ValueError("basis indices must be distinct and ascending")
        if coefficient:
            degrees.add(len(basis))
    if len(degrees) > 1:
        raise ValueError("each nonzero form must be homogeneous in degree")


def wedge(
    left: Mapping[tuple[int, ...], Number],
    right: Mapping[tuple[int, ...], Number],
) -> dict[tuple[int, ...], Number]:
    """Exterior wedge product of two homogeneous forms in an indexed basis.

    Inputs are maps ``{(i, j, ...): coefficient}`` with distinct ascending
    nonnegative basis indices. The empty tuple denotes scalar (degree zero).
    Zero coefficients are removed. Use exact numbers (e.g. ``Fraction``) for
    exact arithmetic. No differential or manifold structure is assumed.
    """
    _validate_form(left)
    _validate_form(right)
    result: dict[tuple[int, ...], Number] = {}
    for a, avalue in left.items():
        if not avalue:
            continue
        for b, bvalue in right.items():
            if not bvalue or set(a).intersection(b):
                continue
            inversions = sum(ai > bj for ai in a for bj in b)
            basis = tuple(sorted(a + b))
            contribution = avalue * bvalue * (-1 if inversions % 2 else 1)
            result[basis] = result.get(basis, 0) + contribution
            if not result[basis]:
                del result[basis]
    return result


def euler_lagrange_residual(
    lagrangian: Callable[[float, float, float], float],
    trajectory: Callable[[float], float],
    t: Real,
    *,
    step: Real = 1e-3,
) -> float:
    """Numerically estimate ``d/dt(dL/dqdot) - dL/dq`` for one scalar path.

    ``L(t, q, qdot)`` and ``q(t)`` must be smooth near ``t``. This is a local
    stationarity diagnostic, not a solver or proof of an extremal path.
    It uses central differences, including an inner numerical ``qdot``.
    """
    point = _finite_real(t, "t")
    h = _finite_real(step, "step")
    if h <= 0:
        raise ValueError("step must be positive")

    def velocity(time: float) -> float:
        return (trajectory(time + h) - trajectory(time - h)) / (2 * h)

    def momentum(time: float) -> float:
        q = trajectory(time)
        v = velocity(time)
        return (
            lagrangian(time, q, v + h) - lagrangian(time, q, v - h)
        ) / (2 * h)

    dp_dt = (momentum(point + h) - momentum(point - h)) / (2 * h)
    q = trajectory(point)
    v = velocity(point)
    dl_dq = (
        lagrangian(point, q + h, v) - lagrangian(point, q - h, v)
    ) / (2 * h)
    return dp_dt - dl_dq
