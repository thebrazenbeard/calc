# Cross-domain mathematical oracle checks (Wave 10)

**Evidence class:** locally executed, model-authored tests with independently
derived mathematical expectations. These checks do **not** constitute formal
verification, an external scientific peer review, or a second implementation
from a different software stack.

**Scope:** selected public APIs in `algebra`, `numerical`, `optimization`,
`collocation`, `spectral`, `fractional`, `intervals`, and `matrix_calc`. This
document records only the cases exercised by
[`tests/test_cross_domain_oracles.py`](../tests/test_cross_domain_oracles.py).

## Fresh reproducible checks

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_cross_domain_oracles.py -v
PYTHONPATH=src python -m unittest discover -s tests -q
```

The focused file contains **20 deterministic test methods** with multiple
independent, bounded cases per method. It uses only the Python standard library.

| Domain | Independent oracle | Meaningful boundary |
| --- | --- | --- |
| Exact polynomial calculus | 60 seeded rational coefficient sets, evaluated using the closed form \(\sum c_k(b^{k+1}-a^{k+1})/(k+1)\) | Signed intervals and antiderivative constants |
| Numerical quadrature | Exact cubic/affine primitives; \(\int_0^\pi\sin x\,dx=2\) | Reversed bounds and panel variation |
| ODE solvers | \(y'=y\Rightarrow y=e^t\), \((y_1,y_2)'=(-y_2,y_1)\Rightarrow(\cos t,\sin t)\) | Euler and RK4 refinement ratios |
| Spectral collocation | A hand-differentiated quintic and \(e^x\) interpolation | Endpoint nodes and degree refinement |
| Symmetric spectral decomposition | Characteristic polynomial of \(\bigl[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\bigr]\), exact Rayleigh extrema | Large diagonal dynamic range \(10^{-300}\) versus \(10^{300}\) |
| Fractional derivatives | \(\Gamma(3/2)=\sqrt\pi/2\), \(\Gamma(5/2)=3\sqrt\pi/4\) | Caputo constant zero versus Riemann–Liouville constant nonzero |
| Interval arithmetic | 70 seeded exact binary64-rational product extrema | Infinity enclosure and division by a zero-containing interval |
| Matrix calculus | Explicit rational Frobenius products and affine matrix derivative coefficients | Central finite difference JVP shapes and coefficients |
| Root finding/minimization | \(\sqrt2\), \(\sqrt[3]2\), and minimizer \(2/7\) | Default tolerances on a smooth non-dyadic root |

## Confirmed defect exposed and repaired

**Medium severity, high confidence:** both the previous bisection and
safeguarded-Newton implementations could raise `ConvergenceError` on the
smooth polynomial \(f(x)=x^2-2\) over `[0,2]` under their **default**
tolerances. Their stopping logic treated a sufficiently narrow bracket as a
terminal error whenever the residual had not yet satisfied `f_tol`, even though
further representable subdivisions would attain the requested residual. The
integration agent fixed the logic in `src/calc/optimization.py` to continue
refinement until the residual criterion is satisfied or a genuine
stagnation/iteration-budget limit occurs. Both default-parameter regression
cases are retained in the oracle test file and pass after the repair.

This is distinct from intentionally rejecting jump discontinuities whose
brackets shrink around a sign change with no actual zero. That rejection
remains an important correctness check in the repository's existing tests.

## Limits of these results

- The expected values arise from basic algebra, differential equations,
  characteristic polynomials, gamma identities and exact `Fraction` arithmetic.
  The tests do not simply call a second code path to produce the answer.
- Fixed-seed tests are repeatable, but no finite set guarantees all inputs or
  all numerical conditioning regimes. No proof of global error bounds is made.
- Numerical convergence ratios exercise benign smooth problems, not stiff
  differential equations, singular integrals, distributions, or arbitrary PDEs.
- Interval checks apply to supported binary64 endpoint operations and do not
  certify transcendental interval functions or complex intervals.
- The tests were authored by a specialist model agent and executed in the
  same Python environment as the implementation. They are **independent
  expectations**, not a fully independent evaluator, provider, or runtime.
- GitHub Actions/cloud CI requires its own remote run evidence; a successful
  local command cannot substitute for that.

## Mathematical source anchors

These are external explanatory sources, **not** copied implementation code:

1. MIT OpenCourseWare, [ES.1803: numerical methods for first-order differential equations](https://ocw.mit.edu/courses/es-1803-differential-equations-spring-2024/mites_1803_s24_topic11.pdf). Describes Euler and RK4, including their characteristic first- and fourth-order convergence behavior.
2. J.-P. Berrut and L. N. Trefethen (2004), [*Barycentric Lagrange Interpolation*, SIAM Review](https://doi.org/10.1137/S0036144502417715). Mathematical basis for the Chebyshev–Lobatto interpolation examples.
3. IEEE Standards Association, [IEEE 1788 interval arithmetic](https://standards.ieee.org/ieee/1788/4431/). Normative background for interval enclosures; the repository does **not** claim full IEEE 1788 conformance.
4. NIST [Digital Library of Mathematical Functions](https://dlmf.nist.gov/5), Gamma Function chapter. Standard definitions and half-integer gamma identities used for fractional power derivative values.

External references remain separately licensed; this repository does not
reproduce their source code.
