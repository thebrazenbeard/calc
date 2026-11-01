"""Dependency-free numerical integration and initial-value solvers.

All methods are approximations, not symbolic proofs. Quadrature accepts a real
function on a finite interval; IVP solvers accept ``f(t, y)`` and a scalar or
fixed-size vector initial state. Values must be finite real numbers. ODE solvers
return the entire trajectory, including both endpoints.

The fixed-grid algorithms use O(n) evaluations and O(1) extra memory for
quadrature. Adaptive Simpson raises when its requested tolerance cannot be
confirmed within ``max_depth`` rather than silently returning an estimate.
"""

from __future__ import annotations

import math
from numbers import Real
from typing import Callable, TypeAlias


State: TypeAlias = float | tuple[float, ...]
Trajectory: TypeAlias = list[tuple[float, State]]

__all__ = ["trapezoid", "simpson", "adaptive_simpson", "euler", "rk4"]


def _real(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _count(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _value(f: Callable[[float], Real], x: float) -> float:
    return _real(f(x), "f(x)")


def _bounds(a: Real, b: Real) -> tuple[float, float]:
    return _real(a, "a"), _real(b, "b")


def trapezoid(f: Callable[[float], Real], a: Real, b: Real, n: int = 100) -> float:
    """Composite trapezoidal quadrature over ``[a, b]`` with ``n`` panels.

    Reversing the bounds reverses the sign. Exact for affine functions in exact
    arithmetic; generally second-order accurate for smooth integrands.
    """
    a, b = _bounds(a, b)
    n = _count(n, "n")
    if a == b:
        return 0.0
    h = (b - a) / n
    endpoints = 0.5 * (_value(f, a) + _value(f, b))
    interior = math.fsum(_value(f, a + i * h) for i in range(1, n))
    return _real(h * (endpoints + interior), "integral estimate")


def simpson(f: Callable[[float], Real], a: Real, b: Real, n: int = 100) -> float:
    """Composite Simpson 1/3 quadrature, requiring a positive even ``n``.

    Exact for cubic polynomials in exact arithmetic, with fourth-order global
    truncation error for sufficiently differentiable integrands.
    """
    a, b = _bounds(a, b)
    n = _count(n, "n")
    if n % 2:
        raise ValueError("Simpson's rule requires an even panel count")
    if a == b:
        return 0.0
    h = (b - a) / n
    odd = math.fsum(_value(f, a + i * h) for i in range(1, n, 2))
    even = math.fsum(_value(f, a + i * h) for i in range(2, n, 2))
    total = (_value(f, a) + _value(f, b) + 4 * odd + 2 * even) * h / 3
    return _real(total, "integral estimate")


def adaptive_simpson(
    f: Callable[[float], Real], a: Real, b: Real,
    tol: float = 1e-10, max_depth: int = 20,
) -> float:
    """Adaptive Simpson quadrature with absolute local error estimate.

    ``tol`` is a requested absolute error threshold estimated by Richardson
    extrapolation, not a rigorous certified error bound for arbitrary functions
    (e.g., highly oscillatory, discontinuous, or singular integrands). Stops
    with ``ValueError`` if the interval cannot meet this criterion by
    ``max_depth`` subdivisions.
    """
    a, b = _bounds(a, b)
    tol = _real(tol, "tol")
    if tol <= 0:
        raise ValueError("tol must be positive")
    if isinstance(max_depth, bool) or not isinstance(max_depth, int):
        raise TypeError("max_depth must be an integer")
    if max_depth < 0:
        raise ValueError("max_depth cannot be negative")
    if a == b:
        return 0.0

    def subdivide(
        left: float, right: float, center: float, fleft: float,
        fcenter: float, fright: float, whole: float,
        budget: float, depth: int,
    ) -> float:
        qleft = (left + center) / 2
        qright = (center + right) / 2
        fl, fr = _value(f, qleft), _value(f, qright)
        left_area = (center - left) * (fleft + 4 * fl + fcenter) / 6
        right_area = (right - center) * (fcenter + 4 * fr + fright) / 6
        difference = left_area + right_area - whole
        if abs(difference) <= 15 * budget:
            return left_area + right_area + difference / 15
        if depth == 0:
            raise ValueError("adaptive Simpson tolerance not met before max_depth")
        return subdivide(
            left, center, qleft, fleft, fl, fcenter,
            left_area, budget / 2, depth - 1,
        ) + subdivide(
            center, right, qright, fcenter, fr, fright,
            right_area, budget / 2, depth - 1,
        )

    # A single Simpson panel can alias an oscillation: cos(8*x) looks exactly
    # constant at the usual first three samples on [0, pi], and its quarter
    # points are also constant. Seed seven non-dyadic panels to disrupt common
    # power-of-two aliases. No finite grid can rule out all adversarial aliases.
    panels = 7
    h = (b - a) / panels
    segment_estimates = []
    for i in range(panels):
        left = a + i * h
        right = b if i == panels - 1 else a + (i + 1) * h
        center = (left + right) / 2
        fa, fm, fb = _value(f, left), _value(f, center), _value(f, right)
        whole = (right - left) * (fa + 4 * fm + fb) / 6
        segment_estimates.append(
            subdivide(left, right, center, fa, fm, fb, whole, tol / panels, max_depth)
        )
    return _real(math.fsum(segment_estimates), "integral estimate")


def _state(value: object, name: str = "state") -> State:
    if isinstance(value, (list, tuple)):
        if not value:
            raise ValueError(f"{name} vector cannot be empty")
        return tuple(_real(part, f"{name}[{i}]") for i, part in enumerate(value))
    return _real(value, name)


def _derivative(
    f: Callable[[float, State], object], t: float, y: State,
) -> State:
    derivative = _state(f(t, y), "derivative")
    if isinstance(y, tuple):
        if not isinstance(derivative, tuple) or len(derivative) != len(y):
            raise ValueError("derivative must have the same vector dimension as state")
    elif isinstance(derivative, tuple):
        raise ValueError("scalar state requires a scalar derivative")
    return derivative


def _add_scaled(y: State, derivative: State, h: float) -> State:
    if isinstance(y, tuple) and isinstance(derivative, tuple):
        return tuple(_real(a + h*b, "computed state") for a, b in zip(y, derivative))
    assert isinstance(y, float) and isinstance(derivative, float)
    return _real(y + h * derivative, "computed state")


def _rk4_step(f: Callable[[float, State], object], t: float, y: State, h: float) -> State:
    k1 = _derivative(f, t, y)
    k2 = _derivative(f, t + h / 2, _add_scaled(y, k1, h / 2))
    k3 = _derivative(f, t + h / 2, _add_scaled(y, k2, h / 2))
    k4 = _derivative(f, t + h, _add_scaled(y, k3, h))
    if isinstance(y, tuple):
        assert all(isinstance(k, tuple) for k in (k1, k2, k3, k4))
        components = tuple(
            _real((a + 2*b + 2*c + d) / 6, "weighted derivative")
            for a, b, c, d in zip(k1, k2, k3, k4)
        )
        return _add_scaled(y, components, h)
    assert all(isinstance(k, float) for k in (k1, k2, k3, k4))
    return _add_scaled(y, (k1 + 2*k2 + 2*k3 + k4) / 6, h)


def _solve_ivp(
    f: Callable[[float, State], object], t0: Real, y0: object,
    t1: Real, steps: int, method: str,
) -> Trajectory:
    t0 = _real(t0, "t0")
    t1 = _real(t1, "t1")
    steps = _count(steps, "steps")
    y = _state(y0, "y0")
    h = (t1 - t0) / steps
    trajectory: Trajectory = [(t0, y)]
    for i in range(steps):
        t = t0 + i * h
        if h != 0:
            if method == "euler":
                y = _add_scaled(y, _derivative(f, t, y), h)
            else:
                y = _rk4_step(f, t, y, h)
        next_t = t1 if i == steps - 1 else t0 + (i + 1) * h
        trajectory.append((_real(next_t, "time"), y))
    return trajectory


def euler(
    f: Callable[[float, State], object], t0: Real, y0: object,
    t1: Real, steps: int = 100,
) -> Trajectory:
    """Explicit Euler IVP integration (first-order global accuracy).

    ``f(t, y)`` computes dy/dt. Returns ``[(t0, y0), ..., (t1, y1)]``;
    vector states are immutable tuples to avoid accidental aliasing.
    """
    return _solve_ivp(f, t0, y0, t1, steps, "euler")


def rk4(
    f: Callable[[float, State], object], t0: Real, y0: object,
    t1: Real, steps: int = 100,
) -> Trajectory:
    """Classical fixed-step RK4 IVP integration (fourth-order global accuracy).

    ``f(t, y)`` computes dy/dt. Returns ``[(t0, y0), ..., (t1, y1)]``.
    This is a fixed-step nonstiff solver, not a stability guarantee.
    """
    return _solve_ivp(f, t0, y0, t1, steps, "rk4")
