# Complex-variable calculus: a bounded numerical implementation

**Status:** Executable educational examples; not a symbolic complex-analysis package or a rigorous holomorphicity/topology prover.

## Mathematical contracts

Let \(f(z)=u(x,y)+iv(x,y)\), where \(z=x+iy\). A holomorphic derivative, when it exists, is the direction-independent limit

\[
f'(z)=\lim_{h\to0,\ h\in\mathbb C} \frac{f(z+h)-f(z)}h.
\]

The `complex_derivative` routine approximates this quotient with **one chosen complex direction** \(d\) and a real step \(h>0\):

\[
D_{h,d}f(z)=\frac{f(z+hd)-f(z-hd)}{2hd},\qquad |d|=1.
\]

This is a central finite difference, not proof the limit exists. In particular, \(f(z)=\overline z\) yields approximately \(+1\) along the real axis and \(-1\) along the imaginary axis.

For a continuously real-differentiable holomorphic function, the Cauchy–Riemann equations are \(u_x=v_y\) and \(u_y=-v_x\). `cauchy_riemann_residual` returns numerical estimates of \((u_x-v_y,\ u_y+v_x)\). Small residuals at **one** sampled point do not certify analyticity in a neighborhood; even exact satisfaction at an isolated point can be insufficient.

For a piecewise continuously differentiable path \(\gamma:[a,b]\to\mathbb C\), the complex contour integral is

\[
\int_\gamma f(z)\,dz = \int_a^b f(\gamma(t))\,\gamma'(t)\,dt.
\]

`contour_integral` applies the composite Simpson rule to the parameter integral. The caller must supply **both** `path` and `path_derivative` rather than relying on a hidden numerical tangent. Panel count must be even. It does not automatically detect singularities between evaluation nodes or provide a certified quadrature error.

The winding-number estimate uses

\[
n(\gamma;c) = \frac{1}{2\pi i}\oint_\gamma\frac{dz}{z-c}.
\]

`winding_number` checks approximate closure, samples to reject encounters with the center, and returns an integer only when the quadrature result lies within the requested tolerance of one. This is a **numerical plausibility check**; paths can hide a winding or pole between samples.

## Reproducible examples

```python
import cmath
import math
from calc.complex_calc import (
    complex_derivative, cauchy_riemann_residual,
    contour_integral, winding_number,
)

point = 1 + 2j
assert abs(complex_derivative(lambda z: z*z, point) - 2*point) < 1e-8
assert max(map(abs, cauchy_riemann_residual(lambda z: z.conjugate(), point))) > 1

circle = lambda t: cmath.exp(1j*t)
tangent = lambda t: 1j*cmath.exp(1j*t)
value = contour_integral(lambda z: 1/z, circle, tangent, 0, 2*math.pi)
assert abs(value - 2j*math.pi) < 1e-8
assert winding_number(circle, tangent, 0) == 1
```

The independent mathematical oracles are \((z^2)'=2z\), the Cauchy–Riemann violation of complex conjugation, and the unit-circle residue \(\oint dz/z=2\pi i\). The test suite also covers orientation reversal, zero winding outside a circle, doubled winding, invalid inputs and sampled singularities.

## Limitations and failure modes

- **Smoothness:** Central finite differences assume the function is regular near the point; a small result cannot distinguish removable singularities, branch cuts, or accidental cancellation.
- **Precision:** The default step is proportional to \(10^{-5}\max(1,|z|)\); rounding and truncation errors have competing scales. It is not tuned to arbitrary user functions.
- **Contour sampling:** Simpson's error estimates presuppose regularity; fast oscillation or an unsampled pole can defeat a finite sampling grid. Check convergence with panel refinements and analyze singularities independently.
- **Topology:** The numerical winding estimate is not a certified homotopy invariant in floating point. It rejects only observed crossings of its center.
- **Callable trust:** Functions and paths supplied as Python callables execute with the caller's privileges. Do not execute untrusted Python expressions.
- **Excluded:** Residue extraction, analytic continuation, Laurent/Taylor series, Riemann surfaces, branch management, interval-certified enclosures, proof-producing complex analysis.

## Primary learning references

1. [MIT OpenCourseWare 18.04 Complex Variables with Applications, lecture notes](https://ocw.mit.edu/courses/18-04-complex-variables-with-applications-spring-2018/resources/lecture-notes/) — complex differentiability, Cauchy–Riemann conditions, line integrals, Cauchy's theorem.
2. [MIT OCW 18.04 Topic 3: Line Integrals and Cauchy's Theorem (PDF)](https://ocw.mit.edu/courses/18-04-complex-variables-with-applications-spring-2018/1968cab6de7c96ead1c83001459a007b_MIT18_04S18_topic3.pdf) — parametrization and contour integrals.
3. [MIT OCW 18.04 Topic 4: Cauchy's Integral Formula (PDF)](https://ocw.mit.edu/courses/18-04-complex-variables-with-applications-spring-2018/b1a90df19f643b974555bbbb93138a48_MIT18_04S18_topic4.pdf) — contour integrals and analytic consequences.

No source implementation or course text was copied. External educational material remains governed by its original license and attribution conditions.
