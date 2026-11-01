"""Bracketing-based scalar roots and unimodal minimization.

These routines require real, finite function values. Bracketing algorithms
assume continuity: a change of sign at a discontinuity does not imply a root.
For that reason, root solvers require a small function residual before they
report success, in addition to a sufficiently small horizontal bracket.

All tolerances are absolute except ``rel_tol``, which scales the root's
horizontal bracket. The algorithms are educational, not globally optimizing
or formally certified root-isolation systems.
"""

from __future__ import annotations

import math
from typing import Callable

Scalar = Callable[[float], float]


class ConvergenceError(ArithmeticError):
    """The specified accuracy was not demonstrated within the budget."""


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite real number, not bool")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be a finite real number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _count(value: object) -> int:
    if type(value) is not int or value <= 0:
        raise ValueError("max_iter must be a positive integer")
    return value


def _parameters(abs_tol: float, rel_tol: float, f_tol: float, max_iter: int) -> tuple[float, float, float, int]:
    abs_tol = _finite(abs_tol, "abs_tol")
    rel_tol = _finite(rel_tol, "rel_tol")
    f_tol = _finite(f_tol, "f_tol")
    if abs_tol <= 0 or rel_tol < 0 or f_tol < 0:
        raise ValueError("require abs_tol > 0, rel_tol >= 0, and f_tol >= 0")
    return abs_tol, rel_tol, f_tol, _count(max_iter)


def _evaluate(function: Scalar, x: float) -> float:
    return _finite(function(x), "function result")


def _midpoint(a: float, b: float) -> float:
    # Avoid a+b overflow and b-a overflow on opposite-sign huge endpoints.
    return a / 2 + b / 2 if a < 0 < b else a + (b - a) / 2


def _opposite(a: float, b: float) -> bool:
    # Multiplying signs can underflow or overflow; sign comparison cannot.
    return (a < 0) != (b < 0)


def _bracket(f: Scalar, a: float, b: float) -> tuple[float, float, float, float]:
    a, b = _finite(a, "a"), _finite(b, "b")
    if a == b:
        raise ValueError("root bracket must have nonzero width")
    if a > b:
        a, b = b, a
    fa, fb = _evaluate(f, a), _evaluate(f, b)
    if fa != 0 and fb != 0 and not _opposite(fa, fb):
        raise ValueError("root bracket endpoints must have opposite signs")
    return a, b, fa, fb


def _width_small(a: float, b: float, abs_tol: float, rel_tol: float) -> bool:
    return b - a <= abs_tol + rel_tol * max(abs(a), abs(b))


def bisection(
    f: Scalar, a: float, b: float, *, abs_tol: float = 1e-12,
    rel_tol: float = 1e-12, f_tol: float = 1e-12, max_iter: int = 200,
) -> float:
    """Find a bracketed continuous root by bisection.

    Raises ``ValueError`` for an invalid bracket or parameters and
    ``ConvergenceError`` when shrinking the bracket does not produce a small
    residual. This latter check helps reject sign-jump discontinuities.
    ``f_tol`` must be scaled to the problem's function range.
    """
    abs_tol, rel_tol, f_tol, max_iter = _parameters(abs_tol, rel_tol, f_tol, max_iter)
    a, b, fa, fb = _bracket(f, a, b)
    if fa == 0:
        return a
    if fb == 0:
        return b
    for _ in range(max_iter):
        mid = _midpoint(a, b)
        fm = _evaluate(f, mid)
        if fm == 0:
            return mid
        if mid == a or mid == b or _width_small(a, b, abs_tol, rel_tol):
            if abs(fm) <= f_tol:
                return mid
            raise ConvergenceError("bracket shrank without satisfying residual tolerance")
        if _opposite(fa, fm):
            b, fb = mid, fm
        else:
            a, fa = mid, fm
    raise ConvergenceError("bisection iteration budget exhausted")


def safeguarded_newton(
    f: Scalar, derivative: Scalar, a: float, b: float, *, x0: float | None = None,
    abs_tol: float = 1e-12, rel_tol: float = 1e-12,
    f_tol: float = 1e-12, max_iter: int = 200,
) -> float:
    """Bracket-preserving Newton iterations with bisection fallback.

    The derivative need not be accurate or nonzero everywhere: unusable or
    out-of-bracket steps fall back to bisection. The function must be continuous
    over a valid sign-changing bracket. An initial point, if supplied, must
    lie inside it.
    """
    abs_tol, rel_tol, f_tol, max_iter = _parameters(abs_tol, rel_tol, f_tol, max_iter)
    a, b, fa, fb = _bracket(f, a, b)
    x = _midpoint(a, b) if x0 is None else _finite(x0, "x0")
    if not a <= x <= b:
        raise ValueError("initial Newton point must be inside root bracket")
    if fa == 0:
        return a
    if fb == 0:
        return b
    for _ in range(max_iter):
        fx = _evaluate(f, x)
        if fx == 0:
            return x
        if x == a or x == b:
            # An initial endpoint is already part of the bracket; it cannot
            # replace the opposite endpoint without losing sign containment.
            pass
        elif _opposite(fa, fx):
            b, fb = x, fx
        else:
            a, fa = x, fx
        if _width_small(a, b, abs_tol, rel_tol):
            if abs(fx) <= f_tol:
                return x
            raise ConvergenceError("bracket shrank without satisfying residual tolerance")
        # Nonfinite or zero derivatives trigger safe bisection. The candidate
        # must stay away from the bracket's edges to ensure shrinkage.
        try:
            slope = float(derivative(x))
        except (ArithmeticError, TypeError, ValueError, OverflowError):
            slope = math.nan
        candidate = math.nan if not math.isfinite(slope) or slope == 0 else x - fx / slope
        margin = (b - a) * 0.1
        if math.isfinite(candidate) and a + margin < candidate < b - margin and candidate != x:
            x = candidate
        else:
            x = _midpoint(a, b)
    raise ConvergenceError("safeguarded Newton iteration budget exhausted")


def golden_section_minimize(
    f: Scalar, a: float, b: float, *, abs_tol: float = 1e-10,
    max_iter: int = 300,
) -> float:
    """Locate a minimum of an assumed unimodal scalar function on ``[a,b]``.

    Returns the approximate *argmin*, not ``f(argmin)``. Unimodality is an
    assumption, not something finite evaluations can establish. Endpoints
    are considered so monotone objectives may minimize at the boundary.
    """
    a, b = _finite(a, "a"), _finite(b, "b")
    abs_tol = _finite(abs_tol, "abs_tol")
    max_iter = _count(max_iter)
    if a >= b or abs_tol <= 0:
        raise ValueError("require a < b and abs_tol > 0")
    if not math.isfinite(b - a):
        raise ValueError("minimization interval width exceeds finite float range")
    ratio = (math.sqrt(5) - 1) / 2
    c = b - ratio * (b - a)
    d = a + ratio * (b - a)
    if not (a < c < d < b):
        raise ConvergenceError("minimization interval has insufficient float resolution")
    fc, fd = _evaluate(f, c), _evaluate(f, d)
    for _ in range(max_iter):
        if b - a <= abs_tol:
            break
        if c == d or c == a or d == b:
            raise ConvergenceError("floating-point resolution exhausted")
        if fc <= fd:
            b, d, fd = d, c, fc
            c = b - ratio * (b - a)
            fc = _evaluate(f, c)
        else:
            a, c, fc = c, d, fd
            d = a + ratio * (b - a)
            fd = _evaluate(f, d)
    else:
        raise ConvergenceError("golden-section iteration budget exhausted")
    candidates = ((a, _evaluate(f, a)), (b, _evaluate(f, b)), (c, fc), (d, fd))
    return min(candidates, key=lambda item: item[1])[0]
