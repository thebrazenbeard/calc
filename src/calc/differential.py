"""Finite-difference differential, multivariable and vector calculus.

These routines *approximate* derivatives of callable real-valued functions.
Functions of several variables accept one coordinate sequence, e.g.
``gradient(lambda p: p[0]**2 + p[1]**2, (1, 2))``. They are not a
symbolic differentiator; near discontinuities, singularities and badly
scaled inputs numerical derivatives can be unreliable.

Only Python's standard library is required.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Callable, Sequence

ScalarFunction = Callable[[float], float]
MultiFunction = Callable[[tuple[float, ...]], float]


def _real(value: float, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _point(point: Sequence[float]) -> tuple[float, ...]:
    try:
        values = tuple(_real(value, "coordinate") for value in point)
    except TypeError as exc:
        raise ValueError("point must be a nonempty sequence of finite reals") from exc
    if not values:
        raise ValueError("point must be nonempty")
    return values


def _step(x: float, step: float | None, *, second_order: bool = False,
          method: str = "central") -> float:
    if step is not None:
        h = _real(step, "step")
        if h <= 0:
            raise ValueError("step must be strictly positive")
    else:
        # The truncation/roundoff balance differs by the finite-difference
        # stencil: O(h), O(h²), and second derivatives, respectively.
        exponent = 0.25 if second_order else (0.5 if method != "central" else 1 / 3)
        h = sys.float_info.epsilon**exponent * max(1.0, abs(x))
    if not math.isfinite(x + h) or x + h == x or x - h == x:
        raise ValueError("step is not representable around coordinate")
    return h


def derivative(func: ScalarFunction, x: float, *, step: float | None = None,
               method: str = "central") -> float:
    """Approximate ``d func / dx`` at ``x``.

    ``method`` is ``central`` (default, second-order accurate), ``forward``
    or ``backward`` (first-order accurate). ``step`` is an absolute positive
    finite displacement; omitted steps are chosen by machine precision and
    coordinate scale. Nonfinite input/output values raise ``ValueError``.
    """
    x = _real(x, "x")
    if method not in ("central", "forward", "backward"):
        raise ValueError("method must be 'central', 'forward', or 'backward'")
    if not callable(func):
        raise TypeError("func must be callable")
    h = _step(x, step, method=method)
    if method == "central":
        return (_real(func(x + h), "function output") -
                _real(func(x - h), "function output")) / (2 * h)
    if method == "forward":
        return (_real(func(x + h), "function output") -
                _real(func(x), "function output")) / h
    return (_real(func(x), "function output") -
            _real(func(x - h), "function output")) / h


def partial_derivative(func: MultiFunction, point: Sequence[float], axis: int,
                       *, step: float | None = None,
                       method: str = "central") -> float:
    """Approximate the coordinate ``axis`` partial derivative at ``point``.

    ``func`` is called with a *tuple* of coordinates. It may evaluate a scalar
    function at every sampled point. The caller's ``point`` is never mutated.
    """
    p = _point(point)
    if isinstance(axis, bool) or not isinstance(axis, int):
        raise TypeError("axis must be an integer")
    if not 0 <= axis < len(p):
        raise IndexError("axis out of range")

    def slice_func(value: float) -> float:
        coordinates = list(p)
        coordinates[axis] = value
        return func(tuple(coordinates))

    return derivative(slice_func, p[axis], step=step, method=method)


def gradient(func: MultiFunction, point: Sequence[float], *,
             step: float | None = None,
             method: str = "central") -> tuple[float, ...]:
    """Approximate the gradient of a scalar-valued multivariable function."""
    p = _point(point)
    return tuple(partial_derivative(func, p, axis, step=step, method=method)
                 for axis in range(len(p)))


def jacobian(functions: Sequence[MultiFunction], point: Sequence[float], *,
             step: float | None = None,
             method: str = "central") -> tuple[tuple[float, ...], ...]:
    """Jacobian matrix with one row per scalar output function.

    For a vector-valued mapping ``F``, pass its component functions as a
    sequence. Rectangular Jacobians are supported. At least one component is
    required.
    """
    p = _point(point)
    components = tuple(functions)
    if not components:
        raise ValueError("functions must have at least one component")
    return tuple(gradient(f, p, step=step, method=method) for f in components)


def hessian(func: MultiFunction, point: Sequence[float], *,
            step: float | None = None) -> tuple[tuple[float, ...], ...]:
    """Approximate a scalar function's symmetric Hessian by central stencils.

    Diagonal entries use the three-point second derivative; mixed entries
    use four-point central differences. Uses the *same* mixed value for both
    symmetric positions, up to numerical error in the actual derivative.
    """
    p = _point(point)
    n = len(p)
    h = tuple(_step(x, step, second_order=True) for x in p)

    def evaluate(deltas: tuple[tuple[int, float], ...]) -> float:
        q = list(p)
        for axis, delta in deltas:
            q[axis] += delta
        return _real(func(tuple(q)), "function output")

    baseline = evaluate(())
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        matrix[i][i] = (evaluate(((i, h[i]),)) - 2 * baseline +
                        evaluate(((i, -h[i]),))) / (h[i] ** 2)
        for j in range(i + 1, n):
            mixed = (
                evaluate(((i, h[i]), (j, h[j])))
                - evaluate(((i, h[i]), (j, -h[j])))
                - evaluate(((i, -h[i]), (j, h[j])))
                + evaluate(((i, -h[i]), (j, -h[j])))
            ) / (4 * h[i] * h[j])
            matrix[i][j] = matrix[j][i] = mixed
    return tuple(tuple(row) for row in matrix)


def directional_derivative(func: MultiFunction, point: Sequence[float],
                           direction: Sequence[float], *,
                           normalize: bool = True,
                           step: float | None = None) -> float:
    """Gradient dotted with ``direction``; unit-direction by default.

    Set ``normalize=False`` for the coordinate directional derivative in the
    raw supplied vector. The zero vector is rejected in either mode.
    """
    p = _point(point)
    vector = _point(direction)
    if len(vector) != len(p):
        raise ValueError("direction and point dimensions differ")
    norm = math.hypot(*vector)
    if norm == 0.0:
        raise ValueError("direction must be nonzero")
    factor = norm if normalize else 1.0
    return sum(g * d / factor for g, d in zip(gradient(func, p, step=step), vector))


def divergence(components: Sequence[MultiFunction], point: Sequence[float], *,
               step: float | None = None) -> float:
    """Approximate the divergence of an n-dimensional vector field."""
    p = _point(point)
    fields = tuple(components)
    if len(fields) != len(p):
        raise ValueError("vector field must have one component per dimension")
    return sum(partial_derivative(f, p, i, step=step)
               for i, f in enumerate(fields))


def curl(components: Sequence[MultiFunction], point: Sequence[float], *,
         step: float | None = None) -> tuple[float, float, float]:
    """Approximate the Cartesian curl of a three-dimensional vector field."""
    p = _point(point)
    fields = tuple(components)
    if len(p) != 3 or len(fields) != 3:
        raise ValueError("curl requires a three-dimensional field and point")
    return (
        partial_derivative(fields[2], p, 1, step=step)
        - partial_derivative(fields[1], p, 2, step=step),
        partial_derivative(fields[0], p, 2, step=step)
        - partial_derivative(fields[2], p, 0, step=step),
        partial_derivative(fields[1], p, 0, step=step)
        - partial_derivative(fields[0], p, 1, step=step),
    )


def laplacian(func: MultiFunction, point: Sequence[float], *,
              step: float | None = None) -> float:
    """Approximate a scalar function's Cartesian Laplacian (Hessian trace)."""
    hh = hessian(func, point, step=step)
    return sum(row[i] for i, row in enumerate(hh))
