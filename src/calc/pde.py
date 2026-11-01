"""Bounded educational solvers for one-dimensional PDEs.

The solvers use explicit finite differences with fixed Dirichlet boundaries.
They produce approximations on a uniform grid; they do not perform general
PDE analysis, certify global error, or support arbitrary boundary conditions.
"""

from __future__ import annotations

import math
from numbers import Real
from typing import Sequence


History = tuple[tuple[float, ...], ...]
_MAX_VALUES = 5_000_000


def _finite(value: Real, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real number")
    try:
        number = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} is outside finite float range") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _grid(initial: Sequence[Real], left: float, right: float) -> tuple[float, ...]:
    if not isinstance(initial, Sequence) or isinstance(initial, (str, bytes)):
        raise TypeError("initial grid must be a sequence of real numbers")
    if len(initial) < 3:
        raise ValueError("grid requires at least two boundary nodes and one interior node")
    if len(initial) > _MAX_VALUES:
        raise ValueError("grid exceeds the five-million-value memory limit")
    try:
        values = tuple(_finite(value, f"initial[{j}]") for j, value in enumerate(initial))
    except TypeError as exc:
        raise TypeError("initial grid must be an iterable of real numbers") from exc
    if not math.isclose(values[0], left, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("initial left boundary disagrees with fixed Dirichlet value")
    if not math.isclose(values[-1], right, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("initial right boundary disagrees with fixed Dirichlet value")
    return (left, *values[1:-1], right)


def _parameters(
    *, scale: Real, dx: Real, dt: Real, steps: int, n: int, name: str,
    limit: float, spatial_power: int,
) -> tuple[float, float, float, float]:
    coefficient = _finite(scale, name)
    spacing = _finite(dx, "dx")
    duration = _finite(dt, "dt")
    if coefficient < 0 or spacing <= 0 or duration <= 0:
        raise ValueError(f"{name} must be nonnegative and dx, dt must be positive")
    if type(steps) is not int:
        raise TypeError("steps must be an integer, not a boolean or float")
    if steps < 0:
        raise ValueError("steps must be nonnegative")
    if n * (steps + 1) > _MAX_VALUES:
        raise ValueError("requested full history exceeds 5,000,000 grid values")
    try:
        ratio = coefficient * duration / spacing**spatial_power
    except ZeroDivisionError as exc:
        raise ValueError("dx is too small for floating-point discretization") from exc
    if not math.isfinite(ratio):
        raise ValueError("CFL ratio is not finite; rescale the grid and step")
    if ratio > limit:
        raise ValueError(f"explicit scheme unstable: CFL ratio {ratio:g} exceeds {limit:g}")
    return coefficient, spacing, duration, ratio


def _next_state(values: tuple[float, ...], r: float, left: float, right: float) -> tuple[float, ...]:
    if r == 0.0:
        return values
    result = [left]
    for j in range(1, len(values) - 1):
        # Under CFL <= 1/2, this is a convex combination. Regrouping
        # avoids premature 2*center overflow before exact cancellation.
        value = math.fsum(((1.0 - 2.0 * r) * values[j],
                           r * values[j - 1], r * values[j + 1]))
        if not math.isfinite(value):
            raise ArithmeticError("nonfinite PDE state; scale initial data")
        result.append(value)
    result.append(right)
    return tuple(result)


def heat_1d_dirichlet(
    initial: Sequence[Real], *, diffusivity: Real, dx: Real, dt: Real,
    steps: int, left: Real = 0, right: Real = 0,
) -> History:
    """Solve ``u_t = diffusivity * u_xx`` by explicit FTCS.

    ``initial[j]`` is ``u(j*dx, 0)`` including boundary nodes. The returned
    tuple contains states at ``t=k*dt`` for ``k=0,...,steps``. ``left`` and
    ``right`` are fixed Dirichlet values and must match the initial endpoints.
    For positive diffusivity the enforced stability bound is
    ``diffusivity*dt/dx**2 <= 1/2``. For zero diffusivity the state is fixed.
    This is an approximation with space/time truncation error.
    """
    l = _finite(left, "left")
    r = _finite(right, "right")
    state = _grid(initial, l, r)
    _, _, _, cfl = _parameters(
        scale=diffusivity, dx=dx, dt=dt, steps=steps,
        n=len(state), name="diffusivity", limit=0.5, spatial_power=2,
    )
    states = [state]
    for _ in range(steps):
        state = _next_state(state, cfl, l, r)
        states.append(state)
    return tuple(states)


def wave_1d_dirichlet(
    initial: Sequence[Real], initial_velocity: Sequence[Real], *,
    wave_speed: Real, dx: Real, dt: Real, steps: int,
    left: Real = 0, right: Real = 0,
) -> History:
    """Solve ``u_tt = wave_speed**2 * u_xx`` by centered differences.

    Initial displacement and velocity are supplied at every grid node.
    ``wave_speed * dt / dx <= 1`` is required for the explicit scheme.
    The first time slice uses the zero-time Taylor expansion with the spatial
    second derivative; subsequent slices use the centered three-time formula.
    Boundary displacement is fixed, so boundary velocity must be zero.
    The result contains states at each ``k*dt`` including time zero.
    """
    l = _finite(left, "left")
    r = _finite(right, "right")
    first = _grid(initial, l, r)
    if not isinstance(initial_velocity, Sequence) or isinstance(initial_velocity, (str, bytes)):
        raise TypeError("initial_velocity must be real sequence")
    try:
        velocities = tuple(_finite(v, f"initial_velocity[{j}]")
                           for j, v in enumerate(initial_velocity))
    except TypeError as exc:
        raise TypeError("initial_velocity must be a real sequence") from exc
    if len(first) != len(velocities):
        raise ValueError("velocity and displacement grids must have the same length")
    if not math.isclose(velocities[0], 0, abs_tol=1e-12, rel_tol=0) or not math.isclose(velocities[-1], 0, abs_tol=1e-12, rel_tol=0):
        raise ValueError("fixed Dirichlet boundaries require zero endpoint velocity")
    _, _, duration, courant = _parameters(
        scale=wave_speed, dx=dx, dt=dt, steps=steps,
        n=len(first), name="wave_speed", limit=1.0, spatial_power=1,
    )
    states = [first]
    if steps == 0:
        return tuple(states)
    if courant == 0.0:
        for step in range(1, steps + 1):
            elapsed = step * duration
            row = [l]
            for j in range(1, len(first) - 1):
                value = first[j] + elapsed * velocities[j]
                if not math.isfinite(value):
                    raise ArithmeticError("nonfinite wave state")
                row.append(value)
            row.append(r)
            states.append(tuple(row))
        return tuple(states)
    c2 = courant * courant
    second = [l]
    for j in range(1, len(first) - 1):
        value = math.fsum(((1.0 - c2) * first[j],
                           0.5 * c2 * first[j - 1],
                           0.5 * c2 * first[j + 1],
                           duration * velocities[j]))
        if not math.isfinite(value):
            raise ArithmeticError("nonfinite wave state")
        second.append(value)
    second.append(r)
    current = tuple(second)
    states.append(current)
    previous = first
    for _ in range(1, steps):
        following = [l]
        for j in range(1, len(first) - 1):
            value = math.fsum((2.0 * (1.0 - c2) * current[j],
                               -previous[j], c2 * current[j - 1],
                               c2 * current[j + 1]))
            if not math.isfinite(value):
                raise ArithmeticError("nonfinite wave state")
            following.append(value)
        following.append(r)
        previous, current = current, tuple(following)
        states.append(current)
    return tuple(states)
