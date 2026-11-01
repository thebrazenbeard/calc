"""Small, explicit Itô stochastic-calculus simulation primitives.

The Brownian generator uses Python's pseudorandom normal sampler. This module
models an *idealized discrete approximation* of Wiener noise; it does not
construct a continuous Brownian sample, certify statistics, or implement a
general stochastic-calculus theorem prover.
"""

from __future__ import annotations

import math
import random
from collections.abc import Callable, Iterable


def _real(value: float, name: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite real number, not bool")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be a finite real number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _count(value: int, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise TypeError(f"{name} must be a positive integer")
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _step(dt: float) -> float:
    value = _real(dt, "dt")
    if value <= 0:
        raise ValueError("dt must be positive")
    return value


def brownian_increments(n: int, dt: float, *, seed: object = None) -> tuple[float, ...]:
    """Return ``n`` independent pseudo-normal increments, each N(0, dt).

    ``seed`` initializes a *private* ``random.Random`` state, making identical
    calls reproducible without mutating the caller's process-global RNG.
    ``dt`` is a positive, constant time step. Independence and Gaussian laws
    are properties of the ideal probabilistic model, only approximated by PRNG.
    """
    _count(n, "n")
    delta = _step(dt)
    rng = random.Random(seed)
    sigma = math.sqrt(delta)
    return tuple(rng.gauss(0.0, sigma) for _ in range(n))


def brownian_path(n: int, dt: float, *, x0: float = 0.0, seed: object = None) -> tuple[float, ...]:
    """Construct discrete times' Brownian positions, including initial ``x0``.

    Positions correspond to t=0, dt, ..., n*dt; times are not returned.
    """
    start = _real(x0, "x0")
    path = [start]
    for dw in brownian_increments(n, dt, seed=seed):
        position = path[-1] + dw
        if not math.isfinite(position):
            raise ValueError("Brownian path overflowed to a nonfinite value")
        path.append(position)
    return tuple(path)


def quadratic_variation(path: Iterable[float]) -> float:
    """Return sum of adjacent squared increments along a sampled path.

    For standard Brownian motion and refining partitions, this converges in
    probability to elapsed time. One finite realization need not equal time.
    """
    positions = tuple(_real(x, "path entry") for x in path)
    if len(positions) < 2:
        raise ValueError("quadratic variation requires at least two positions")
    variation = math.fsum((right - left) ** 2 for left, right in zip(positions, positions[1:]))
    if not math.isfinite(variation):
        raise ValueError("quadratic variation overflowed")
    return variation


def euler_maruyama(
    drift: Callable[[float, float], float],
    diffusion: Callable[[float, float], float],
    t0: float,
    y0: float,
    dt: float,
    steps: int,
    *,
    increments: Iterable[float] | None = None,
    seed: object = None,
) -> tuple[tuple[float, float], ...]:
    """Solve a scalar Itô SDE by explicit Euler–Maruyama.

    Equation: dY_t = drift(t,Y_t) dt + diffusion(t,Y_t) dW_t.
    At each left endpoint, update ``y += drift*dt + diffusion*dW``.
    Returns ``((t0, y0), ..., (t0 + steps*dt, y_steps))``.

    Give either ``increments`` (user-supplied Brownian ΔW samples) *or* seed.
    Caller-provided increments are useful for repeatable convergence tests and
    coupling simulations. No adaptive refinement, multi-dimensional noise,
    strong convergence estimate, or discontinuity handling is supplied.
    """
    if not callable(drift) or not callable(diffusion):
        raise TypeError("drift and diffusion must be callables")
    count = _count(steps, "steps")
    delta = _step(dt)
    start_t = _real(t0, "t0")
    value = _real(y0, "y0")
    if increments is None:
        noise = brownian_increments(count, delta, seed=seed)
    else:
        if seed is not None:
            raise ValueError("provide increments or seed, not both")
        noise = tuple(_real(dw, "Brownian increment") for dw in increments)
        if len(noise) != count:
            raise ValueError("number of Brownian increments must equal steps")

    trajectory: list[tuple[float, float]] = [(start_t, value)]
    for i, dw in enumerate(noise):
        t = _real(start_t + i * delta, "time")
        a = _real(drift(t, value), "drift value")
        b = _real(diffusion(t, value), "diffusion value")
        value = _real(value + a * delta + b * dw, "updated state")
        next_t = _real(start_t + (i + 1) * delta, "time")
        trajectory.append((next_t, value))
    return tuple(trajectory)


def geometric_brownian_terminal(
    s0: float, mu: float, sigma: float, elapsed: float, w_elapsed: float
) -> float:
    """Analytical GBM solution given one Wiener displacement at elapsed time.

    For dS = μS dt + σS dW, S(T) = S(0) exp[(μ-σ²/2)T + σW(T)].
    This is an *Itô* result; the -σ²/2 correction is essential.
    ``s0`` must be nonnegative. ``w_elapsed`` is supplied by caller so this is
    a deterministic calculation, not a sampler.
    """
    start = _real(s0, "s0")
    growth = _real(mu, "mu")
    volatility = _real(sigma, "sigma")
    horizon = _real(elapsed, "elapsed")
    displacement = _real(w_elapsed, "w_elapsed")
    if start < 0 or horizon < 0:
        raise ValueError("s0 and elapsed must be nonnegative")
    if start == 0:
        return 0.0
    try:
        result = start * math.exp((growth - volatility * volatility / 2) * horizon + volatility * displacement)
    except OverflowError as exc:
        raise ValueError("geometric Brownian solution overflowed") from exc
    return _real(result, "geometric Brownian solution")
