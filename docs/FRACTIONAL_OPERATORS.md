# Fractional operators: Caputo, Riemann–Liouville, Grünwald–Letnikov

**Status:** three narrowly implemented operations, not a general fractional differential equation solver. These operators differ in their treatment of initial data, even when they share the word *derivative*. Source code: [`src/calc/fractional.py`](../src/calc/fractional.py); executable witnesses: [`tests/test_fractional.py`](../tests/test_fractional.py).

## 1. Left-terminal-zero power derivatives

For a real power function $f(t)=t^\beta$, $\beta\ge0$, and $0<\alpha<1$, the Riemann–Liouville derivative at $x>0$ is

\[
{}_0D_x^\alpha x^\beta
=\frac{\Gamma(\beta+1)}{\Gamma(\beta-\alpha+1)}x^{\beta-\alpha}.
\]

For the Caputo derivative of the **nonconstant** power $\beta>0$, the same expression applies under these restrictions. For $f(t)=1$, however,

\[
{}_0^C D_x^\alpha 1 = 0, \qquad
{}_0D_x^\alpha 1 = \frac{x^{-\alpha}}{\Gamma(1-\alpha)}.
\]

That difference is essential when specifying initial-value problems. At $\alpha=1$, both implemented formulas return the ordinary first derivative. The $\alpha=0$ boundary is deliberately excluded because endpoint conventions differ.

At $x=0$, only a **finite right-hand limit** is returned: zero when $\beta>\alpha$, $\Gamma(\beta+1)$ when $\beta=\alpha$, or a domain error for a divergent power. The Caputo constant is always zero, and the ordinary first derivative of a constant is zero.

```python
import math
from calc.fractional import caputo_power_derivative, riemann_liouville_power_derivative

assert caputo_power_derivative(0.5, 0, 4) == 0
assert math.isclose(riemann_liouville_power_derivative(0.5, 0, 4),
                    1 / (2 * math.sqrt(math.pi)))
assert math.isclose(caputo_power_derivative(0.5, 1, 4),
                    4 / math.sqrt(math.pi))
```

The gamma ratio is evaluated with logarithms in IEEE binary floating point. For power exponents at least one million, a first-order asymptotic gamma-ratio expansion avoids catastrophic subtraction of nearly equal `lgamma` values (but introduces an $O(\beta^{-2})$ remainder for fixed noninteger order). It is not an arbitrary-precision symbolic result; extreme powers and exponents can still be rejected if the result is nonfinite.

## 2. Sampled Grünwald–Letnikov approximation

For $0<\alpha<1$, $x>0$, and a uniform grid $h=x/N$, the left-terminal-zero approximation is

\[
G_{h}^{\alpha}f(x)=h^{-\alpha}\sum_{k=0}^{N}(-1)^k
\binom{\alpha}{k}f(x-kh),
\qquad w_{k+1}=w_k\frac{k-\alpha}{k+1},\quad w_0=1.
\]

The recurrence avoids computing generalized binomial coefficients via enormous factorials. A compensated sum and logarithmic final scaling reduce some floating-point overflow and cancellation risks, but **do not certify the error**.

```python
from calc.fractional import grunwald_letnikov_derivative

approximately = grunwald_letnikov_derivative(lambda t: t*t, 1.0, 0.5, steps=2048)
# Analytical RL/Caputo half-derivative of t² at 1: 8/(3 sqrt(pi))
```

The sampled method approximates a **Riemann–Liouville-type** derivative for appropriate regular functions. It does **not** automatically produce the Caputo derivative when $f(0)\ne0$: a constant has nonzero RL derivative but zero Caputo derivative. Convergence depends on smoothness, boundaries, grid resolution, and the operator definition; singular functions and discontinuities need separate analysis.

Implementation constraints: `steps` is an integer from 1 through 200,000; evaluation uses only $[0,x]$. Time and storage are both $O(N)$. The method rejects unrepresentable sample spacing, non-real or nonfinite function outputs, and nonfinite convolution terms. It makes no automatic choice of grid size and supplies no rigorous accuracy estimate.

## Independently derived regression values

| Input | Expected value | Why it matters |
| --- | --- | --- |
| Caputo half derivative of 1 | 0 | Not the RL definition |
| RL half derivative of 1 at 4 | $1/(2\sqrt\pi)$ | Initial-data distinction |
| Half derivative of $t$ at 4 | $4/\sqrt\pi$ | Gamma power rule |
| Half derivative of $t^2$ at 1 | $8/(3\sqrt\pi)$ | Numerical convergence witness |
| Fractional power at $x=0$, $\beta<\alpha$ | Diverges; error | Must not return an invented finite value |

## Sources and provenance

- Francesco Mainardi, *Essentials of Fractional Calculus*, NIST-hosted lecture manuscript: [Mainardi PDF](https://math.nist.gov/~DLozier/SF21/SF21slides/Mainardi.pdf). Definitions, gamma-function power rules and distinctions; no text or code copied.
- NIST Digital Library of Mathematical Functions, [§2.6 Distributional Methods](https://dlmf.nist.gov/2.6). Background on Riemann–Liouville fractional integrals and special-function context.
- *The Grünwald–Letnikov method for fractional differential equations*, *Computers & Mathematics with Applications* (2011), [article record](https://www.sciencedirect.com/science/article/pii/S0898122111002173). Numerical context and finite-difference interpretation.

All code and examples in this repository are independently written. External sources are references, not vendored dependencies; license permissions for reproduction of external material must be evaluated separately.
