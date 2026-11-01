# Scalar root finding and optimization

Status: **implemented as a small standard-library numerical module**, not a general optimization framework.

This module distinguishes **root finding** (`f(x)=0`) from **minimization** (`argmin f(x)`). Both occur in calculus, but their convergence assumptions are different. The algorithms live in `src/calc/optimization.py` and have independently calculated tests in `tests/test_optimization.py`.

## Bracketed bisection

Assume `f` is continuous on `[a,b]` and `f(a)` and `f(b)` have opposite signs. The intermediate value theorem then guarantees at least one root. Each iteration replaces one endpoint with the midpoint, preserving the sign-change bracket. After `n` bisections, the interval width is at most `(b-a)/2^n` in exact arithmetic.

```python
from calc.optimization import bisection

root = bisection(lambda x: x**3 - 2, 1, 2)
assert abs(root - 2**(1/3)) < 1e-10
```

This implementation checks **two independent conditions** before reporting approximate success: the horizontal bracket has become narrow, *and* `abs(f(x)) <= f_tol`. An exactly zero function value can return immediately. A sign-jump `f(x)=-1` below a threshold and `+1` above it has a sign-changing bracket but no root; the solver raises `ConvergenceError` instead of silently presenting the jump location as a root. No finite sample can prove continuity, however, and an arbitrarily tiny discontinuity could still evade this heuristic.

## Safeguarded Newton

Newton's proposed step is `x_new = x - f(x)/f'(x)`. This implementation keeps a valid sign-changing bracket and only accepts Newton steps that land safely inside it; otherwise it bisects. Zero, unusable, or non-finite derivative values cause a bisection fallback. As with bisection, small residual and horizontal accuracy are both required unless the result is an exact evaluated zero.

```python
from calc.optimization import safeguarded_newton

root = safeguarded_newton(
    lambda x: x*x - 2,
    lambda x: 2*x,
    0, 2, x0=1,
)
```

It does **not** discover brackets, prove derivative validity, or isolate all roots. It is slower than specialized production implementations for many workloads.

## Golden-section bounded minimization

The golden ratio `φ=(1+√5)/2` supplies interior points with a reusable spacing ratio `1/φ≈0.618`. For a **unimodal** function, comparing function values at these points narrows the interval containing a minimum. The function returns an estimated *argument* (`x`), not the objective value.

```python
from calc.optimization import golden_section_minimize

argmin = golden_section_minimize(lambda x: (x-3)**2, -1, 10)
assert abs(argmin - 3) < 1e-7
```

It includes endpoints as candidates, allowing a monotone objective to select a boundary. The solver does not establish unimodality; oscillatory or multimodal functions may have a local or wrong minimum. If large constants flatten observable objective values through IEEE-754 rounding, the reported `x` can be less accurate than the requested horizontal tolerance suggests.

## Contract and limitations

- Python 3.11+; standard library only; real and finite input/output values.
- `abs_tol` is a positive horizontal absolute tolerance, `rel_tol` is a nonnegative horizontal relative tolerance for roots, `f_tol` is a nonnegative *absolute* function-residual tolerance. A badly scaled `f` may need a custom `f_tol`.
- `max_iter` is a positive integer; exhausted budget and unverified convergence raise `ConvergenceError` rather than returning an unqualified point.
- Bracketed root solvers allow reversed endpoint order; minimization requires `a < b`.
- Neither gradients, multidimensional optimization, constraints, stochastic objectives, global extrema, nonsmooth subgradient methods, automatic differentiation, nor formal error certification are provided here.

## References and provenance

All code and worked examples here were written for this repository; external material is used as references, not copied.

1. SciPy official [scalar root finding API](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html), with separate documented [bisection](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.bisect.html) and [Newton](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.newton.html) procedures.
2. SciPy official [scalar minimization and root-finding overview](https://docs.scipy.org/doc/scipy/reference/optimize.html) and [optimization tutorial](https://docs.scipy.org/doc/scipy/tutorial/optimize.html). These distinguish bounded minimization, golden search, Brent's method, and bracketed roots.

For production scientific problems, established solver libraries provide additional algorithms, diagnostics, and compatibility checks. This module is an inspectable teaching implementation.
