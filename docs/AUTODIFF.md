# Forward-mode automatic differentiation

**Implementation:** [`src/calc/autodiff.py`](../src/calc/autodiff.py)  
**Status:** tested teaching implementation, not a full differentiation framework.

## Mathematical contract

A dual number is written \(a+b\varepsilon\), with \(\varepsilon^2=0\). The first component is the computed value; the second is its derivative along a seeded input direction.

The product rule follows from dual multiplication:

\[
(a+b\varepsilon)(c+d\varepsilon)=ac+(ad+bc)\varepsilon.
\]

For a differentiable function \(f\colon\mathbb R^n\to\mathbb R^m\), propagating input tangents \(v\) produces a Jacobian-vector product \(J_f(x)v\), subject to floating-point rounding. This is *automatic differentiation* (AD), not finite-difference differentiation: it does not perturb \(x\) by a small step and subtract nearby values.

## Executable examples

With `PYTHONPATH=src` from the repository root:

```python
from calc.autodiff import derivative, directional_derivative, gradient, jacobian, sin, exp

assert abs(derivative(lambda x: x**3, 2) - 12) < 1e-12

# User functions accept a tuple of Dual values for multivariable AD.
g = gradient(lambda v: v[0]**2 + sin(v[1]), (3, 0))
assert g == (6.0, 1.0)

J = jacobian(lambda v: (v[0]*v[1], exp(v[0])), (2, 3))
# Two output coordinates x*y and exp(x): 2 rows, 2 columns.
assert abs(J[0][0] - 3) < 1e-12
assert abs(J[0][1] - 2) < 1e-12

directional_derivative(lambda v: v[0]**2*v[1], (2, 3), (-1, 2))  # -4
```

One scalar-output derivative and a directional derivative each require **one function evaluation**. `gradient` and `jacobian` rerun the user function once per input variable (about \(n\) evaluations). This is efficient when the number of input dimensions is small; reverse mode is usually preferable for a scalar output with very many inputs, but is not implemented here.

## Supported primitives

| Operation | Derivative rule | Domain restriction |
| --- | --- | --- |
| `+`, `-`, `*`, `/` | linearity, product, quotient | division requires nonzero denominator |
| `a**p` | \(p a^{p-1}a'\) | negative `a` only for fixed integral `p`; undefined points reject |
| `a**b` (variable `b`) | \(a^b[b'\log(a)+ba'/a]\) | positive `a` |
| `sin(a)`, `cos(a)` | `cos(a) a'`, `-sin(a) a'` | real arguments in radians |
| `exp(a)` | `exp(a) a'` | real input; finite output |
| `log(a)` | `a'/a` | positive input |
| `sqrt(a)` | `a'/(2 sqrt(a))` | strictly positive input for finite derivative |

## Deliberate limitations

- The module uses **real floating-point numbers**, not arbitrary precision, complex duals, symbolic formulae, or interval-certified derivatives.
- Use its provided `sin`, `cos`, `exp`, `log`, and `sqrt`. Passing `Dual` to ordinary `math.sin`, for example, raises `TypeError` instead of silently discarding derivative information.
- It does not trace arbitrary native functions, C-extension calls, array kernels, comprehensions relying on comparisons of duals, side effects, or branches. Nonsmooth functions and discontinuities need separate domain analysis.
- `gradient` and `jacobian` evaluate user functions repeatedly with different seeds; stateful or random functions can produce inconsistent derivatives.
- A Jacobian result uses **output rows / input columns**. A direction vector is **not automatically normalized**.
- Derivative claims rely on the function being differentiable at the queried point and on the supported operations. An algorithm can return a finite number even when a model's mathematical assumptions fail elsewhere.

## Validation

Run `PYTHONPATH=src python -m unittest discover -s tests -p test_autodiff.py -v`.

Independent analytic oracles include the product, quotient, chain and power rules, a three-input gradient, a three-output Jacobian, and negative/domain cases. These tests verify selected cases, not all possible programs.

## Sources

- MIT OpenCourseWare 18.S096, *Forward Automatic Differentiation via Dual Numbers*, Edelman & Johnson (2023): https://ocw.mit.edu/courses/18-s096-matrix-calculus-for-machine-learning-and-beyond-january-iap-2023/resources/ocw_18s096_lecture05-part2-new_2023jan27_mp4/
- MIT OpenCourseWare 18.S096, *Forward and Reverse-Mode Automatic Differentiation* lecture notes: https://ocw.mit.edu/courses/18-s096-matrix-calculus-for-machine-learning-and-beyond-january-iap-2023/mit18_s096iap23_lec08.pdf
- Baydin et al., *Automatic Differentiation in Machine Learning: A Survey*, JMLR 18(153), 2018: https://jmlr.org/papers/v18/17-468.html

References provide mathematical context. This implementation and examples were written for this repository, not copied from those sources.
