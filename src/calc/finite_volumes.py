"""Conservative first-order finite volumes for 1D periodic linear advection.

This is a deliberately narrow solver for ``q_t + a q_x = 0`` with constant
signed velocity ``a``, uniformly spaced cells, and periodic boundary data.
Its state is a tuple of finite cell *averages*, not nodal point samples.

The upwind method uses the donor-cell interface flux and an explicit step
with Courant number ``nu = |a| dt / dx <= 1``. Under this restriction each
updated average is a convex combination of its old value and its upwind
neighbor; no claim of high-order accuracy or support for other PDEs is made.
"""

from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction
import math
from numbers import Real


_MAX_CELLS = 100_000
_MAX_STEPS = 100_000
_MAX_CELL_UPDATES = 2_000_000


def _finite_real(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{label} must be a real number (not bool)")
    try:
        result = float(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"{label} must be representable as a finite float") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be a finite real number")
    return result


def _cells(values: Sequence[Real]) -> tuple[float, ...]:
    # Requiring a sized sequence avoids consuming an unbounded generator.
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise TypeError("cells must be a finite sequence of real cell averages")
    count = len(values)
    if not 2 <= count <= _MAX_CELLS:
        raise ValueError("periodic grid requires 2 to 100000 cells")
    return tuple(_finite_real(value, "cell average") for value in values)


def _spacing(dx: Real, count: int) -> float:
    spacing = _finite_real(dx, "dx")
    if spacing <= 0:
        raise ValueError("dx must be strictly positive")
    if not math.isfinite(spacing * count):
        raise ValueError("period length is not representable as a finite float")
    return spacing


def _courant(velocity: float, dt: float, dx: float) -> tuple[float, float]:
    # Float multiplication can overflow even for a representable quotient.
    # Rational reconstruction computes CFL for the *input binary64 values*.
    nu_exact = (Fraction.from_float(abs(velocity)) * Fraction.from_float(dt)
                / Fraction.from_float(dx))
    if nu_exact > 1:
        raise ValueError("explicit upwind step violates the CFL limit |velocity|*dt/dx <= 1")
    return float(nu_exact), float(1 - nu_exact)


def _step_validated(state: tuple[float, ...], velocity: float,
                    move: float, stay: float) -> tuple[float, ...]:
    if move == 0.0:
        return state
    size = len(state)
    result: list[float] = []
    for i, old in enumerate(state):
        upwind = state[(i - 1) % size] if velocity >= 0 else state[(i + 1) % size]
        # Both weights are nonnegative and add to one in exact arithmetic.
        # Accumulate these two finite products without avoidable cancellation.
        value = math.fsum((stay * old, move * upwind))
        if not math.isfinite(value):
            raise ValueError("updated cell average exceeded finite float range")
        result.append(value)
    return tuple(result)


def upwind_fluxes(cells: Sequence[Real], velocity: Real) -> tuple[float, ...]:
    """Return periodic donor-cell fluxes ``F[i]`` at face ``i+1/2``.

    For nonnegative velocity, ``F[i]=a*q[i]``. For negative velocity,
    ``F[i]=a*q[(i+1) mod n]``. An unrepresentable face flux is rejected.
    """
    state = _cells(cells)
    a = _finite_real(velocity, "velocity")
    result = []
    for i in range(len(state)):
        donor = state[i] if a >= 0 else state[(i + 1) % len(state)]
        flux = a * donor
        if not math.isfinite(flux):
            raise ValueError("upwind flux exceeded finite float range")
        result.append(flux)
    return tuple(result)


def upwind_step(cells: Sequence[Real], *, dx: Real,
                velocity: Real, dt: Real) -> tuple[float, ...]:
    """Advance uniform periodic cell averages by one first-order CFL step.

    ``dt=0`` or ``velocity=0`` returns the same values. The computed update
    is algebraically equivalent to ``q_i - dt/dx*(F_i-F_{i-1})`` while the
    convex-combination form preserves the discrete maximum principle.
    """
    state = _cells(cells)
    spacing = _spacing(dx, len(state))
    a = _finite_real(velocity, "velocity")
    timestep = _finite_real(dt, "dt")
    if timestep < 0:
        raise ValueError("dt must be nonnegative")
    move, stay = _courant(a, timestep, spacing)
    return _step_validated(state, a, move, stay)


def advance_upwind(cells: Sequence[Real], *, dx: Real, velocity: Real,
                   dt: Real, steps: int) -> tuple[float, ...]:
    """Apply a bounded number of identical explicit periodic upwind steps.

    ``steps=0`` returns a validated copy of the original cell averages.
    Restricts total cell updates to avoid accidentally enormous workloads.
    """
    state = _cells(cells)
    spacing = _spacing(dx, len(state))
    a = _finite_real(velocity, "velocity")
    timestep = _finite_real(dt, "dt")
    if timestep < 0:
        raise ValueError("dt must be nonnegative")
    if type(steps) is not int:
        raise TypeError("steps must be an integer (not bool)")
    if steps < 0 or steps > _MAX_STEPS or len(state) * steps > _MAX_CELL_UPDATES:
        raise ValueError("steps exceed bounded finite-volume workload")
    move, stay = _courant(a, timestep, spacing)
    for _ in range(steps):
        state = _step_validated(state, a, move, stay)
    return state


def cell_mass(cells: Sequence[Real], *, dx: Real) -> float:
    """Return ``sum(q_i*dx)`` with exceptional exact-binary64 fallback.

    A finite mass is recoverable even if individual float products overflow
    before cancellation. The returned number remains a float approximation.
    """
    state = _cells(cells)
    spacing = _spacing(dx, len(state))
    terms = [value * spacing for value in state]
    if all(math.isfinite(term) for term in terms):
        try:
            result = math.fsum(terms)
        except (OverflowError, ValueError):
            pass
        else:
            if math.isfinite(result):
                return result
    exact = sum((Fraction.from_float(value) * Fraction.from_float(spacing)
                 for value in state), Fraction(0))
    try:
        result = float(exact)
    except OverflowError as exc:
        raise ValueError("cell mass is not representable as a finite float") from exc
    if not math.isfinite(result):
        raise ValueError("cell mass is not representable as a finite float")
    return result


def periodic_total_variation(cells: Sequence[Real]) -> float:
    """Return ``Σ_i |q_{i+1}-q_i|`` around the periodic grid."""
    state = _cells(cells)
    try:
        result = math.fsum(abs(state[(i + 1) % len(state)] - state[i])
                           for i in range(len(state)))
    except (OverflowError, ValueError) as exc:
        raise ValueError("periodic total variation is not representable") from exc
    if not math.isfinite(result):
        raise ValueError("periodic total variation is not representable")
    return result


__all__ = ["upwind_fluxes", "upwind_step", "advance_upwind",
           "cell_mass", "periodic_total_variation"]
