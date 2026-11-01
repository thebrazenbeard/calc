"""Bounded scalar Itô/Stratonovich conversion and finite-grid experiments.

The formulas assume a single standard Brownian driver and sufficiently smooth
state-dependent coefficients. All stochastic increments come from the caller;
this module makes no claim to construct a continuous-time stochastic process.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from fractions import Fraction
from itertools import islice


ScalarCoefficient = Callable[[float, float], float]
_MAX_STEPS = 1_000_000


def _real(value: float, name: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite real number, not bool")
    try:
        finite = float(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise TypeError(f"{name} must be a finite real number") from exc
    if not math.isfinite(finite):
        raise ValueError(f"{name} must be finite")
    return finite


def _finite_sum(values: Iterable[float], name: str) -> float:
    try:
        total = math.fsum(values)
    except OverflowError as exc:
        raise ValueError(f"{name} overflowed") from exc
    return _real(total, name)


def _coefficients(*values: tuple[str, ScalarCoefficient]) -> None:
    for name, coefficient in values:
        if not callable(coefficient):
            raise TypeError(f"{name} must be callable as {name}(t, y)")


def _converted_drift(drift: float, diffusion: ScalarCoefficient,
                     diffusion_dx: ScalarCoefficient, t: float, y: float,
                     sign: int, name: str) -> float:
    b = _real(diffusion(t, y), "diffusion value")
    bx = _real(diffusion_dx(t, y), "diffusion state derivative")
    correction = 0.5 * b * bx
    converted = drift + sign * correction
    if math.isfinite(converted):
        return converted
    # The correction, or its addition to the drift, may overflow even when
    # the full scalar conversion is finite due to cancellation. Reconstruct
    # exact binary64 inputs only for this exceptional path.
    exact = (Fraction.from_float(drift) +
             sign * Fraction.from_float(b) * Fraction.from_float(bx) / 2)
    try:
        converted = float(exact)
    except OverflowError as exc:
        raise ValueError(f"{name} is outside finite float range") from exc
    return _real(converted, name)


def stratonovich_to_ito_drift(
    drift_stratonovich: ScalarCoefficient,
    diffusion: ScalarCoefficient,
    diffusion_dx: ScalarCoefficient,
    t: float,
    y: float,
) -> float:
    """Evaluate ``a_Ito(t,y) = a_Strat(t,y) + 0.5*b(t,y)*d_y b(t,y)``.

    ``diffusion_dx`` is the **user-supplied analytic state derivative** of
    ``diffusion`` at fixed ``t``. No numerical differentiation is performed.
    The conversion applies to scalar, single-driver Brownian SDEs with the
    regularity required for classical Itô/Stratonovich equivalence.
    """
    _coefficients(("drift_stratonovich", drift_stratonovich),
                  ("diffusion", diffusion), ("diffusion_dx", diffusion_dx))
    time, state = _real(t, "time"), _real(y, "state")
    drift = _real(drift_stratonovich(time, state), "Stratonovich drift")
    return _converted_drift(drift, diffusion, diffusion_dx, time, state, 1, "Itô drift")


def ito_to_stratonovich_drift(
    drift_ito: ScalarCoefficient,
    diffusion: ScalarCoefficient,
    diffusion_dx: ScalarCoefficient,
    t: float,
    y: float,
) -> float:
    """Evaluate the inverse conversion ``a_Strat = a_Ito - 0.5*b*d_y b``."""
    _coefficients(("drift_ito", drift_ito), ("diffusion", diffusion),
                  ("diffusion_dx", diffusion_dx))
    time, state = _real(t, "time"), _real(y, "state")
    drift = _real(drift_ito(time, state), "Itô drift")
    return _converted_drift(drift, diffusion, diffusion_dx, time, state,
                            -1, "Stratonovich drift")


def _noise(increments: Iterable[float], steps: int) -> tuple[float, ...]:
    """Read at most steps+1 values, so infinite generators cannot hang us."""
    try:
        iterator = iter(increments)
    except TypeError as exc:
        raise TypeError("increments must be an iterable of real numbers") from exc
    noise = tuple(_real(dw, "Brownian increment") for dw in islice(iterator, steps + 1))
    if len(noise) != steps:
        raise ValueError("number of Brownian increments must equal steps")
    return noise


def stratonovich_heun(
    drift: ScalarCoefficient,
    diffusion: ScalarCoefficient,
    t0: float,
    y0: float,
    dt: float,
    steps: int,
    *,
    increments: Iterable[float],
) -> tuple[tuple[float, float], ...]:
    """Explicit stochastic Heun update for a scalar Stratonovich SDE.

    For ``dY = a_S(t,Y) dt + b(t,Y) ◦ dW``, first predict with left-endpoint
    Euler, then average both coefficients at the left endpoint and the
    *predicted* right endpoint. Uses precisely the supplied finite increments.
    For smooth coefficients and a single Brownian driver this is a standard
    Stratonovich-consistent predictor-corrector; it is not an exact solver or
    a convergence guarantee for arbitrary/drift-discontinuous coefficients.
    """
    _coefficients(("drift", drift), ("diffusion", diffusion))
    if not isinstance(steps, int) or isinstance(steps, bool):
        raise TypeError("steps must be a positive integer")
    if steps <= 0 or steps > _MAX_STEPS:
        raise ValueError(f"steps must be between 1 and {_MAX_STEPS}")
    delta = _real(dt, "dt")
    if delta <= 0:
        raise ValueError("dt must be positive")
    start_t = _real(t0, "t0")
    value = _real(y0, "y0")
    noise = _noise(increments, steps)

    result = [(start_t, value)]
    for i, dw in enumerate(noise):
        t = _real(start_t + i * delta, "time")
        next_t = _real(start_t + (i + 1) * delta, "time")
        if next_t <= t:
            raise ValueError("dt cannot advance time at this magnitude")

        a0 = _real(drift(t, value), "left drift")
        b0 = _real(diffusion(t, value), "left diffusion")
        predicted = _real(value + a0 * delta + b0 * dw, "predicted state")
        a1 = _real(drift(next_t, predicted), "predicted drift")
        b1 = _real(diffusion(next_t, predicted), "predicted diffusion")
        value = _real(value + 0.5 * (a0 + a1) * delta
                      + 0.5 * (b0 + b1) * dw, "updated state")
        result.append((next_t, value))
    return tuple(result)


@dataclass(frozen=True)
class BrownianSquareIntegrals:
    """Finite-grid Riemann sums for ``d(W²)`` from a path starting at W=0.

    The algebraic correction is observed pathwise at finite resolution, but
    finite sampled quadratic variation need *not* equal physical elapsed time.
    """

    terminal_displacement: float
    terminal_square: float
    ito_left: float
    stratonovich_midpoint: float
    quadratic_variation: float


def brownian_square_integrals(increments: Iterable[float]) -> BrownianSquareIntegrals:
    """Compare left-endpoint Itô and midpoint Stratonovich sums for ``2W dW``.

    From ``W_0=0`` and sampled increments ``dW_k`` the finite-grid identities
    are ``W_n² = 2*sum(W_k*dW_k) + sum(dW_k²)`` and
    ``W_n² = 2*sum((W_k+W_{k+1})/2 * dW_k)`` (roundoff aside).
    For actual Brownian motion under partition refinement, the quadratic
    variation converges toward elapsed time, giving the Itô ``+dt`` term.
    """
    left_terms: list[float] = []
    midpoint_terms: list[float] = []
    squares: list[float] = []
    position = 0.0
    count = 0
    for increment in increments:
        if count == _MAX_STEPS:
            raise ValueError(f"at most {_MAX_STEPS} increments are supported")
        dw = _real(increment, "Brownian increment")
        following = _real(position + dw, "Brownian position")
        left_terms.append(_real(2 * position * dw, "left-endpoint integral term"))
        midpoint_terms.append(_real((position + following) * dw,
                                    "midpoint integral term"))
        squares.append(_real(dw * dw, "quadratic variation term"))
        position = following
        count += 1
    if count == 0:
        raise ValueError("at least one Brownian increment is required")
    return BrownianSquareIntegrals(
        terminal_displacement=position,
        terminal_square=_real(position * position, "terminal Brownian square"),
        ito_left=_finite_sum(left_terms, "Itô left sum"),
        stratonovich_midpoint=_finite_sum(midpoint_terms, "Stratonovich midpoint sum"),
        quadratic_variation=_finite_sum(squares, "quadratic variation"),
    )
