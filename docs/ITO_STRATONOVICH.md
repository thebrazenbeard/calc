# Itô and Stratonovich: a finite-grid scalar comparison

**Status:** bounded, instructional implementation in `calc.stochastic_conventions`. This supplements the Itô Euler–Maruyama solver in [`STOCHASTIC.md`](STOCHASTIC.md); neither module implements general stochastic analysis.

## Two conventions for Brownian noise

For one standard Brownian motion (W_t), a scalar Stratonovich SDE is written

\[
dY_t = a_S(t,Y_t)\,dt+b(t,Y_t)\circ dW_t.
\]

With sufficiently regular coefficients, the **equivalent Itô SDE** is

\[
dY_t = \underbrace{\left(a_S(t,Y_t)+\frac12b(t,Y_t)\frac{\partial b}{\partial y}(t,Y_t)\right)}_{a_I(t,Y_t)}\,dt
       + b(t,Y_t)\,dW_t.
\]

The `stratonovich_to_ito_drift` and `ito_to_stratonovich_drift` functions evaluate this relation **at one `(t, y)` point**. The caller provides `diffusion_dx(t, y)`, the derivative of diffusion with respect to the *state* at fixed time. The functions do not estimate derivatives. If (b) does not depend on (y), the correction vanishes.

For example (a_S(t,y)=\mu y), (b(t,y)=\sigma y), and (\partial_y b=\sigma) give (a_I(t,y)=(\mu+\sigma^2/2)y). Confusing these drifts changes the model, even when the diffusion is unchanged.

## Two numerical updates under the **same** increments

For a positive fixed timestep (h), user-supplied increments (\Delta W_k), time (t_k), and state (Y_k):

| Scheme | Convention | Finite-grid update |
|---|---|---|
| Euler–Maruyama, `calc.stochastic.euler_maruyama` | Itô | (Y_{k+1}=Y_k+a_I(t_k,Y_k)h+b(t_k,Y_k)\Delta W_k) |
| Stochastic Heun, `stratonovich_heun` | Stratonovich | Predict (P_k=Y_k+a_S(t_k,Y_k)h+b(t_k,Y_k)\Delta W_k), then set (Y_{k+1}=Y_k+\frac12[a_S(t_k,Y_k)+a_S(t_{k+1},P_k)]h+\frac12[b(t_k,Y_k)+b(t_{k+1},P_k)]\Delta W_k) |

This is **explicit stochastic Heun**: the right-endpoint coefficients use the *predicted* state, not a solved implicit equation. For smooth scalar coefficients with one Brownian driver it is a standard Stratonovich-oriented predictor–corrector. The caller supplies the exact same `increments` sequence to both solvers; one solver does not independently resample randomness. The implemented Heun routine never draws its own random samples.

```python
import math
from calc.stochastic import brownian_increments, euler_maruyama
from calc.stochastic_conventions import stratonovich_heun, stratonovich_to_ito_drift

mu, sigma = 0.1, 0.3
dt, n = 0.001, 1000
noise = brownian_increments(n, dt, seed=29)

a_strat = lambda t, y: mu * y
b = lambda t, y: sigma * y
db_dy = lambda t, y: sigma
a_ito = lambda t, y: stratonovich_to_ito_drift(a_strat, b, db_dy, t, y)

strat_path = stratonovich_heun(a_strat, b, 0, 2, dt, n, increments=noise)
ito_path = euler_maruyama(a_ito, b, 0, 2, dt, n, increments=noise)
exact_terminal = 2 * math.exp(mu * n * dt + sigma * sum(noise))
```

The **exact** common continuous model here has terminal value (Y_T=Y_0\exp(\mu T+\sigma W_T)). The numerically integrated paths are approximations, so they need not coincide step by step or exactly equal the analytic terminal value. `calc.stochastic.geometric_brownian_terminal` provides an independent Itô-form analytic comparator when supplied the converted drift (\mu+\sigma^2/2).

When (b=0), stochastic Heun becomes the **deterministic explicit trapezoidal (Heun) method**—not explicit Euler. With constant drift and diffusion, both schemes reproduce the finite-grid additive solution exactly up to floating-point rounding.

## A discrete witness for Itô's extra term

Consider the function (f(W)=W^2). Given finite increments starting at (W_0=0), the algebraic identity

\[
W_n^2=2\sum_{k=0}^{n-1}W_k\Delta W_k+\sum_{k=0}^{n-1}(\Delta W_k)^2
\]

holds for any values of the increments (up to floating-point roundoff). The first sum is the **left-endpoint** discretization of the Itô integral (\int 2W\,dW). The second is sampled **quadratic variation**. In contrast,

\[
2\sum_{k=0}^{n-1}\frac{W_k+W_{k+1}}2\Delta W_k=W_n^2
\]

is the corresponding midpoint sum. The difference between midpoint and left-point sums is exactly sampled quadratic variation on a finite grid. `brownian_square_integrals(increments)` returns both sums, the quadratic variation, and terminal displacement/square in a frozen `BrownianSquareIntegrals` record.

For true Brownian motion as the partition is refined, quadratic variation tends to elapsed time in the appropriate stochastic sense. Consequently the continuous formulas are (d(W_t^2)=2W_t\,dW_t+dt) for Itô and (d(W_t^2)=2W_t\circ dW_t) for Stratonovich. **One finite simulated path does not have quadratic variation equal to elapsed time exactly.**

## Preconditions, safeguards, and limits

- Scalar real-valued state, one Brownian driver, forward constant positive timestep, positive integer step count up to 1,000,000. Only caller-provided finite increments are accepted; their count must match `steps`.
- Input values, derivative/coefficients, predicted states, outputs, and times must be finite. Updates reject time-grid stagnation (`t + dt == t` at floating-point magnitude) or overflow. Infinite noise generators are capped by `steps + 1` reads before rejection.
- No automatic check proves the given increments are Brownian or independent, or that a supplied `diffusion_dx` is truly the derivative. Smoothness, existence/uniqueness, and convergence conditions remain the caller's mathematical responsibility.
- Neither exact simulation nor error bounds, adaptivity, irregular grids, multidimensional/correlated Brownian noise, Lévy noise, stochastic jumps, general semimartingales, or Malliavin calculus are implemented. This is not a validated financial or scientific simulation engine.

## Primary reference starting points

- David Cohen, [*Lecture notes for mini-course on numerical methods for SDEs*](https://www.math.chalmers.se/~cohend/MT2016/LNminicourseSDEv1.pdf), §4: the two stochastic integral conventions and conversion term.
- [*Numerical Integration of Stochastic Differential Equations: The Heun Algorithm Revisited and the Itô–Stratonovich Calculus*](https://arpi.unipi.it/retrieve/4228425c-3f3f-4b0a-b8c2-c96f0ec191ad/entropy-27-00910.pdf), research article discussing Heun predictor–corrector and convention choices.
- Jonathan Goodman, Courant Institute, [*Stochastic Calculus*, Lesson 3: Itô integral](https://math.nyu.edu/faculty/goodman/teaching/StochCalc2018/notes/Lesson3.pdf), for the role of quadratic variation in the Itô chain rule.

The formulas, examples, implementation, and tests here are independently written and do not reproduce source code or extended text from these references.
