# Finite-horizon scalar optimal control (exact discrete LQR)

The module [`calc.control`](../src/calc/control.py) implements a deliberately
small instance of optimal control: a **scalar, linear, deterministic,
discrete-time** system with an unconstrained quadratic objective. All
coefficients, feedback gains, states and costs are rational numbers, computed
exactly by Python `Fraction` with explicit resource limits.

## Problem and assumptions

For a fixed finite horizon \(N\), choose controls \(u_0,\ldots,u_{N-1}\)
subject to

\[
x_{k+1}=a x_k+b u_k
\]

to minimize

\[
J = \sum_{k=0}^{N-1}(q x_k^2 + r u_k^2)+q_T x_N^2.
\]

The inputs \(a,b\) may be any rationals, while \(q\geq0\),
\(q_T\geq0\), and \(r>0\). Strictly positive \(r\) ensures every
single-stage minimization is strictly convex, hence the optimum control
is unique. The dynamics need not be stable, and \(b\) can be zero.

For the optimal remaining cost \(V_k(x)=P_k x^2\), Bellman's principle
gives the backward scalar Riccati recursion:

\[
\begin{aligned}
P_N &= q_T,\\
D_k &= r+b^2P_{k+1},\\
K_k &= \frac{abP_{k+1}}{D_k},\\
P_k &= q+\frac{a^2P_{k+1}r}{D_k},\\
u_k^* &= -K_k x_k.
\end{aligned}
\]

The denominator is strictly positive. The formula for \(P_k\) is
algebraically equivalent to the usual Riccati subtraction form, but its
nonnegative factors avoid artificial numerical cancellation.

## Executable example

```python
from fractions import Fraction
from calc.control import solve_scalar_lqr, optimal_trajectory, evaluate_trajectory

solution = solve_scalar_lqr(a=1, b=1, q=1, r=1, terminal_q=1, horizon=2)
assert solution.cost_to_go == (Fraction(8, 5), Fraction(3, 2), Fraction(1))
assert solution.gains == (Fraction(3, 5), Fraction(1, 2))

path = optimal_trajectory(solution, initial_state=5)
assert path.states == (Fraction(5), Fraction(2), Fraction(1))
assert path.controls == (Fraction(-3), Fraction(-1))
assert path.total_cost == 40
assert solution.value(5) == 40

# Compare an arbitrary, suboptimal open-loop control sequence.
candidate = evaluate_trajectory(solution, 5, controls=(0, 0))
assert candidate.total_cost > path.total_cost
```

The exact completing-the-square identity independently checks each stage:

\[
qx^2+ru^2+P_{k+1}(ax+bu)^2-P_kx^2
=D_k(u+K_kx)^2\geq0.
\]

The test suite verifies this identity for 80 reproducible rational instances,
checks the hand-derived one- and two-stage results, enumerates 289 two-stage
rational open-loop alternatives, and validates invalid/coefficient-size cases.

## API

- `solve_scalar_lqr(a, b, q, r, terminal_q, horizon)` gives an immutable
  `ScalarLQR` with `cost_to_go` and `gains`.
- `solution.control(state, stage)` evaluates \(-K_k x\).
- `solution.value(state, stage=0)` gives \(P_k x^2\).
- `optimal_trajectory(solution, initial_state)` returns exact states,
  controls and total cost.
- `evaluate_trajectory(solution, initial_state, controls)` evaluates any
  supplied sequence; consuming at most `horizon + 1` values prevents
  infinitely yielding iterables from hanging validation.

Horizon is bounded between 0 and 64, inclusive. Coefficients, states,
controls and intermediate rational growth are explicitly bounded; over-budget
inputs are rejected instead of switching silently to approximate arithmetic.

## What this does **not** implement

This is not Pontryagin's continuous-time maximum principle, a Riccati ODE,
constrained or stochastic optimal control, vector/matrix LQR, nonlinear MPC,
model identification, or a claim of optimality for a real physical plant
outside the exact stated model. A grid enumeration of comparison controls
is only a test witness; optimality follows from the algebraic completed
square and positivity of \(D_k\).

## Independent reference

- Stanford University, EE363, *Linear Quadratic Regulator: Discrete-Time
  Finite Horizon* lecture material, [course lectures](https://ee363.stanford.edu/lectures.html)
  and [discrete-time LQR notes](https://ee363.stanford.edu/archive/lectures/dlqr.pdf).

This implementation and exposition are original; the reference is linked for
mathematical background, not copied as reusable source code.
