# Stochastic calculus: an executable introductory slice

**Status: limited implementation** — `calc.stochastic` demonstrates discretized Brownian motion and scalar Itô SDEs. It is not a general stochastic-analysis, Malliavin-calculus, or numerical-finance library.

## Model and units

For ideal standard Brownian motion \(W_t\), increments over a time interval \(\Delta t > 0\) are Gaussian with mean zero and variance \(\Delta t\):

\[\Delta W \sim N(0,\Delta t),\qquad W_{t_{k+1}} = W_{t_k}+\Delta W_k.\]

`brownian_increments(n, dt, seed=...)` generates `n` pseudorandom samples using `random.Random.gauss(0, sqrt(dt))`. `brownian_path` accumulates them, starting at user-supplied `x0`. An identical seed reproduces the same samples within a given Python runtime; this is deterministic simulation, **not cryptographic randomness**. The continuous-time Brownian object is idealized; the returned object is only a discrete grid.

For a scalar Itô stochastic differential equation (SDE)

\[dY_t=a(t,Y_t)\,dt+b(t,Y_t)\,dW_t,\]

`euler_maruyama` uses left-endpoint sampling:

\[Y_{k+1}=Y_k+a(t_k,Y_k)\Delta t+b(t_k,Y_k)\Delta W_k.\]

The method has an explicit input: callables `drift(t, y)` and `diffusion(t, y)`, initial time and state, positive `dt`, and positive integer step count. A caller may provide the exact `increments` sequence to compare two numerical schemes under the **same** Brownian realization. Otherwise the solver internally samples seeded increments. Return values are a tuple of `(time, state)` pairs, including the initial condition. Only forward, constant time steps and **scalar** SDEs are supported.

```python
from calc.stochastic import brownian_path, euler_maruyama

w = brownian_path(100, 0.01, seed=2026)
trajectory = euler_maruyama(
    lambda t, y: 0.05*y, lambda t, y: 0.2*y,
    t0=0, y0=1, dt=0.01, steps=100, seed=2026,
)
```

## Two independent mathematical checks

1. If \(b=0\), the Euler–Maruyama update reduces exactly to the deterministic explicit Euler step. For \(Y'=Y\) starting at 1 with ten steps of size 0.1, the terminal numerical value is exactly \(1.1^{10}\) up to floating-point rounding.
2. The geometric Brownian SDE \(dS_t=\mu S_tdt+\sigma S_tdW_t\) has the **Itô** solution

   \[S_T=S_0\exp((\mu-\tfrac12\sigma^2)T+\sigma W_T).\]

   `geometric_brownian_terminal` evaluates this independent analytical comparator from a supplied Wiener displacement. The correction term \(-\sigma^2/2\) distinguishes this Itô interpretation from an incorrectly applied ordinary chain rule.

`quadratic_variation(path)` computes \(\sum_k (W_{t_{k+1}}-W_{t_k})^2\) for a sampled path. Standard Brownian quadratic variation converges toward elapsed time as the mesh refines in the suitable probabilistic sense; an individual finite path is not guaranteed to equal elapsed time.

## Boundaries and numerical warnings

- The solver assumes well-posed drift/diffusion functions and does not establish strong or weak convergence for arbitrary coefficients. Smooth/Lipschitz and growth assumptions matter; violating them can make explicit Euler schemes diverge.
- Pseudorandom finite samples cannot establish true independence or Gaussian distribution. Large-sample tests here check only one reproducible empirical realization.
- The methods do not construct Itô integrals of arbitrary adapted processes, stochastic derivatives, martingale filtrations, stopping times, Lévy jumps, Stratonovich integrals, multi-dimensional correlated noise, or Malliavin derivatives.
- Deterministic exceptions are raised for invalid time steps, nonfinite values, mismatched increments, and nonfinite computed states. Runtime overflow may still be possible with extreme scale choices.
- This package is instructional mathematics, not an audited pricing, forecasting, or risk-management system.

## Reference starting points (source links, not copied implementations)

- MIT OpenCourseWare, *Advanced Stochastic Processes* (2013), [lecture notes](https://ocw.mit.edu/courses/15-070j-advanced-stochastic-processes-fall-2013/resources/lecture-notes/), especially Brownian motion, quadratic variation, and Itô calculus lectures.
- MIT OpenCourseWare, *Topics in Mathematics with Applications in Finance* (2024), [Week 14](https://ocw.mit.edu/courses/18-642-topics-in-mathematics-with-applications-in-finance-fall-2024/pages/week-14/) for Itô stochastic calculus and geometric Brownian motion.
- MIT OpenCourseWare, *Topics in Mathematics with Applications in Finance* (2013), [Lecture 18: Itō Calculus](https://ocw.mit.edu/courses/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/ef2c66c8079ba656210ad1fd4a5e2fa8_MIT18_S096F13_lecnote18.pdf).

These sources explain underlying mathematical objects; the present implementation and tests were written independently.
