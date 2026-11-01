"""Exact scalar finite-horizon discrete-time linear quadratic control.

For ``x[k+1] = a*x[k] + b*u[k]`` and constant rational coefficients, the
objective is ``sum(q*x[k]**2 + r*u[k]**2, k=0..N-1) + terminal_q*x[N]**2``.
``r > 0``, ``q >= 0`` and ``terminal_q >= 0`` ensure a unique optimal control
at each stage. All arithmetic is exact ``Fraction`` with explicit budgets.

This is not continuous-time optimal control, robust or constrained LQR,
nonlinear control, or a generic Riccati differential-equation solver.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction
from itertools import islice


Rational = int | Fraction
_MAX_HORIZON = 64
_MAX_INPUT_BITS = 4096
_MAX_INTERMEDIATE_BITS = 8192


def _rational(value: Rational, name: str) -> Fraction:
    if type(value) is not int and not isinstance(value, Fraction):
        raise TypeError(f"{name} must be an int or Fraction, not float/bool")
    exact = Fraction(value)
    if max(abs(exact.numerator).bit_length(), exact.denominator.bit_length()) > _MAX_INPUT_BITS:
        raise ValueError(f"{name} exceeds input rational-bit budget")
    return exact


def _checked(value: Fraction, name: str) -> Fraction:
    if max(abs(value.numerator).bit_length(), value.denominator.bit_length()) > _MAX_INTERMEDIATE_BITS:
        raise ValueError(f"{name} exceeds intermediate rational-bit budget")
    return value


@dataclass(frozen=True, slots=True)
class ScalarLQR:
    """The exact backward Riccati solution for a scalar rational problem.

    ``cost_to_go[k]`` is ``P_k`` with optimal value ``V_k(x)=P_k*x*x``.
    ``gains[k]`` is ``K_k`` with optimal feedback ``u_k=-K_k*x_k``.
    There are ``horizon + 1`` cost-to-go coefficients and ``horizon`` gains.
    """

    a: Fraction
    b: Fraction
    q: Fraction
    r: Fraction
    terminal_q: Fraction
    horizon: int
    cost_to_go: tuple[Fraction, ...]
    gains: tuple[Fraction, ...]

    def control(self, state: Rational, stage: int) -> Fraction:
        """The unique exact optimal control at a known stage and state."""
        if type(stage) is not int:
            raise TypeError("stage must be an integer")
        if not 0 <= stage < self.horizon:
            raise ValueError("stage must be within the active control horizon")
        return _checked(-self.gains[stage] * _rational(state, "state"), "control")

    def value(self, state: Rational, stage: int = 0) -> Fraction:
        """Optimal remaining cost from stage ``stage``, including terminal."""
        if type(stage) is not int:
            raise TypeError("stage must be an integer")
        if not 0 <= stage <= self.horizon:
            raise ValueError("stage must be within the inclusive value horizon")
        x = _rational(state, "state")
        return _checked(self.cost_to_go[stage] * x * x, "value")


@dataclass(frozen=True, slots=True)
class ScalarTrajectory:
    """Immutable states, controls, and their exact finite-horizon objective."""

    states: tuple[Fraction, ...]
    controls: tuple[Fraction, ...]
    total_cost: Fraction


def solve_scalar_lqr(
    a: Rational, b: Rational, q: Rational, r: Rational,
    terminal_q: Rational, horizon: int,
) -> ScalarLQR:
    """Solve the unconstrained rational scalar discrete-time LQR by Bellman.

    ``P_N=terminal_q``, ``D_k=r+b²*P_(k+1)``, then
    ``K_k=a*b*P_(k+1)/D_k`` and
    ``P_k=q+a²*P_(k+1)*r/D_k``.

    The positive ``r`` makes all denominators strictly positive and every
    stage control strictly convex. A zero horizon is allowed and gives only
    terminal cost. Horizons exceed 64 or over-budget fractions are rejected.
    """
    if type(horizon) is not int:
        raise TypeError("horizon must be an integer")
    if not 0 <= horizon <= _MAX_HORIZON:
        raise ValueError(f"horizon must be between 0 and {_MAX_HORIZON}")
    aa, bb = _rational(a, "a"), _rational(b, "b")
    qq, rr, terminal = (_rational(q, "q"), _rational(r, "r"),
                        _rational(terminal_q, "terminal_q"))
    if qq < 0 or terminal < 0:
        raise ValueError("state and terminal weights must be nonnegative")
    if rr <= 0:
        raise ValueError("control weight r must be strictly positive")

    costs = [Fraction(0)] * (horizon + 1)
    gains = [Fraction(0)] * horizon
    costs[horizon] = terminal
    for stage in range(horizon - 1, -1, -1):
        following = costs[stage + 1]
        denominator = _checked(rr + bb * bb * following, "Riccati denominator")
        gains[stage] = _checked(aa * bb * following / denominator,
                                "feedback gain")
        costs[stage] = _checked(qq + aa * aa * following * rr / denominator,
                                "cost-to-go coefficient")
    return ScalarLQR(aa, bb, qq, rr, terminal, horizon, tuple(costs), tuple(gains))


def evaluate_trajectory(
    solution: ScalarLQR, initial_state: Rational,
    controls: Iterable[Rational],
) -> ScalarTrajectory:
    """Evaluate **any** supplied finite control sequence against the objective.

    Consumes at most ``horizon+1`` entries, so even an infinite generator
    fails after bounded work. All controls and evolving states must remain
    within the explicit rational-bit budget.
    """
    if not isinstance(solution, ScalarLQR):
        raise TypeError("solution must be a ScalarLQR")
    try:
        iterator = iter(controls)
    except TypeError as exc:
        raise TypeError("controls must be an iterable") from exc
    raw = tuple(islice(iterator, solution.horizon + 1))
    if len(raw) != solution.horizon:
        raise ValueError("control sequence length must exactly equal horizon")
    chosen = tuple(_rational(v, "control") for v in raw)
    states = [_rational(initial_state, "initial_state")]
    cost = Fraction(0)
    for u in chosen:
        x = states[-1]
        cost = _checked(cost + solution.q * x * x + solution.r * u * u,
                        "trajectory cost")
        states.append(_checked(solution.a * x + solution.b * u, "trajectory state"))
    cost = _checked(cost + solution.terminal_q * states[-1] * states[-1],
                    "trajectory cost")
    return ScalarTrajectory(tuple(states), chosen, cost)


def optimal_trajectory(solution: ScalarLQR, initial_state: Rational) -> ScalarTrajectory:
    """Simulate unique exact optimal time-varying feedback gains."""
    if not isinstance(solution, ScalarLQR):
        raise TypeError("solution must be a ScalarLQR")
    state = _rational(initial_state, "initial_state")
    chosen: list[Fraction] = []
    for stage in range(solution.horizon):
        control = solution.control(state, stage)
        chosen.append(control)
        state = _checked(solution.a * state + solution.b * control,
                         "trajectory state")
    answer = evaluate_trajectory(solution, initial_state, chosen)
    if answer.total_cost != solution.value(initial_state):
        raise ArithmeticError("exact Bellman identity failed")
    return answer


__all__ = [
    "ScalarLQR", "ScalarTrajectory", "solve_scalar_lqr",
    "evaluate_trajectory", "optimal_trajectory",
]
