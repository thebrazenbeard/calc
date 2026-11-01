# Continuum Calc

An open, source-conscious **atlas of calculus and algebra**: a taxonomy of mathematical ideas, tested executable examples, and a place to study where different calculi connect and diverge.

> **Scope:** No finite catalogue can guarantee it contains *every* type of calculus or algebra. The index is extensible, and an indexed subject is not necessarily implemented. This project separates *implemented*, *documented*, and *planned* work instead of representing a numerical method as a complete theory.

## Where to start

- `src/calc/` — dependency-free Python implementations of selected operations, designed to be auditable and easy to test.
- `tests/` — independently calculated mathematical examples, boundaries, and failure cases.
- `docs/CATALOG.md` — domain map and what belongs in each mathematical family.
- `docs/SOURCES.md` — upstream references, provenance, and license notes.
- `docs/ADVANCED_METHODS.md` — assumptions and limits of advanced examples.
- `docs/REVIEW.md` — adversarial technical review and known gaps.
- `docs/BRANCHES.md` — rules and topic-branch manifest.
- `docs/topics/` — individually researched subject guides and worked examples.
- `src/calc/topic_manifest.json` — machine-readable scope, tests, source IDs, and formal-proof status for each indexed topic.
- `docs/EVIDENCE_MODEL.md` — how catalogue labels map to evidence rather than whole-field claims.
- `docs/REVIEW_WAVE5.md` — hostile numerical review and regression evidence from the subsequent build wave.
- `docs/REVIEW_WAVE6.md` — independent algebraic, stochastic and spectral verification, including repaired edge cases.
- `docs/REVIEW_WAVE7.md` — independent verification of geometric algebra and numerical PDE modules.
- `docs/REVIEW_WAVE8.md` — independent verification of spectral/umbral/categorical/homological/tensor modules.
- `docs/REVIEW_WAVE9.md` — independent verification of finite measure, control and advanced algebra modules.
- `docs/CI_QUALIFICATION.md` — build/test workflow, offline package smoke test, and cloud CI evidence limits.
- `docs/MATH_VERIFICATION.md` — independent mathematical oracles and repaired numerical root regressions.
- `docs/REVIEW_WAVE10.md` — hostile Wave 10 mathematical and evidence-integrity checks.
- `explorer/index.html` — dependency-free interactive visual explorer for selected numerical calculus examples.

The CLI shows the subject map without executing arbitrary mathematical expressions:

```bash
PYTHONPATH=src python -m calc topics
PYTHONPATH=src python -m calc topic fractional
PYTHONPATH=src python -m calc evidence fractional --json
PYTHONPATH=src python -m unittest discover -s tests -v
```

For a browser-based illustration, open `explorer/index.html` directly. It uses
bounded numerical plots and clearly labels approximations; it is not a
mathematical proof system or an arbitrary-expression evaluator.

## Mathematical domains

| Family | Examples | Meaning of coverage |
| --- | --- | --- |
| Foundations | Limits, ordinary differential/integral calculus, complex analysis, measure theory | Core concepts and selected algorithms |
| Multivariable / geometry | Gradients, Jacobians, Hessians, tensor calculus, exterior forms, manifolds | Selected computations; geometry topics catalogued |
| Dynamics | ODE/PDE, calculus of variations, functional analysis, optimal control | Some numerical IVP methods; theory catalogued |
| Generalized calculus | Fractional, stochastic, Itô, Malliavin, time-scale, q-, discrete, nonsmooth | Explicit specialized examples, not general solvers |
| Computational calculus | Quadrature, differentiation, finite elements, automatic differentiation | Selected algorithms and references |
| Algebra | Polynomial, linear, abstract, Boolean, Lie, Clifford, geometric, tropical, homological | Selected exact algebra operations and references |

**Names matter:** Lambda calculus is a formal system of computation; relational calculus is a query formalism. Neither is a variant of real-variable differentiation. The catalogue distinguishes them rather than falsely asserting a common derivative operator.

## Executable building blocks

The Python package is intentionally modular. Each file declares its numerical or algebraic contract instead of claiming a universal solution:

| Module | Included demonstrations |
| --- | --- |
| `algebra` | Exact rational polynomials, matrix operations, linear systems, modular arithmetic |
| `differential` | Derivatives, gradients, Jacobians, Hessians, vector differential operators |
| `numerical` | Trapezoid/Simpson quadrature, scalar/vector Euler and RK4 |
| `advanced` | Fractional power integrals, q-/finite differences, wedges, variational residuals |
| `autodiff` | Forward-mode dual-number differentiation |
| `symbolic` | Restricted expression trees and symbolic differentiation, not arbitrary evaluation |
| `stochastic` | Brownian paths, quadratic variation, Euler–Maruyama, geometric Brownian motion |
| `geometry` | Metrics, coordinate transformations, Christoffel symbols, exterior derivatives |
| `structures` | Finite-group validation, prime fields, Boolean algebra, tropical semirings |
| `complex_calc` | Complex finite-difference diagnostics and numerical contour integrals |
| `pde` | Bounded 1D heat/wave finite-difference solvers with stability checks |
| `optimization` | Bisection, safeguarded Newton and golden-section search |
| `transforms` | Direct DFT/IDFT, convolutions and a finite z-transform |
| `lambda_calc` | Safe capture-avoiding formal lambda-term reduction |
| `matrix_calc` | Exact Frobenius algebra and selected matrix derivatives |
| `grobner` | Small exact rational multivariate polynomial ideals and Buchberger reduction |
| `time_scales` | Finite-grid delta/nabla differences and integrals |
| `fractional` | Explicit Caputo/Riemann–Liouville power derivatives and bounded GL approximation |
| `intervals` | Binary64 interval enclosures with outward rounding |
| `forms` | Exact rational exterior forms, wedge products, Euclidean Hodge stars |
| `nonsmooth` | Exact scalar convex subgradients and proximal examples |
| `discrete` | Rational finite differences, summation by parts and power sums |
| `stochastic_conventions` | Itô–Stratonovich drift conversion and finite-grid Stratonovich Heun |
| `spectral` | Symmetric Jacobi eigendecomposition and residual diagnostics; not a PDE spectral solver |
| `lie` | Exact matrix commutators, Jacobi identity and small Lie subalgebra structure constants |
| `clifford` | Exact rational multivectors for signed diagonal Clifford/geometric algebras |
| `distributions` | Compact rational step functions, Dirac jump derivatives and weak pairings |
| `fem` | Exact-rational 1D piecewise-linear Poisson finite-element assembly and solution |
| `finite_volumes` | Conservative 1D periodic upwind advection with CFL checks |
| `collocation` | Chebyshev–Lobatto interpolation and nodal polynomial differentiation |
| `umbral` | Exact falling-factorial difference operators, Stirling numbers and Newton interpolation |
| `categories` | Finite categories, composition/identity verification, functors, natural transformations |
| `homology` | Exact rational chain complexes, boundary ranks and Betti-number computations |
| `tensors` | Fixed-basis rational tensor products, permutations and contractions |
| `measure` | Exact finite-atom measures, integrals, Radon–Nikodym densities and conditional expectation |
| `control` | Exact scalar discrete-time finite-horizon LQR and feedback policies |
| `free_algebra` | Rational noncommutative word polynomials and formal derivations |
| `operator_algebra` | Rational matrix adjoints, projections and finite star-algebra closure |
| `differential_algebra` | Rational finite-jet polynomials with formal total derivatives |
| `limits` | Exact rational epsilon–delta certificates for bounded polynomial cases; numerical diagnostics are not proofs |
| `integrals` | Exact rational polynomial antiderivatives, definite integrals and elementary FTC identities |

For instance, exact algebra and numerical integration can coexist without conflating their guarantees:

```python
from fractions import Fraction
from calc.algebra import Polynomial
from calc.numerical import simpson

p = Polynomial((1, 0, 1))   # 1 + x^2
assert p.differentiate().coefficients == (Fraction(0), Fraction(2))
assert abs(simpson(lambda x: x*x, 0, 1, n=100) - 1/3) < 1e-12
```

The fractional, stochastic, symbolic, geometry, complex, and PDE modules each cover **specific examples**, not entire theories. A passing test proves its stated contract only. See module documentation for limitations.

## Branch-per-paradigm policy

The repository description promises a branch for each paradigm. `main` is the shared tested foundation and navigable index. Curated `topic/<slug>` branches can contain a focused subject guide or experiment while inheriting the common foundation. Branch content must be independently verified and linked from the index; a branch's existence does not imply a complete implementation of its subject.

## Design rules

1. State mathematical domains, units, assumptions, approximation errors, and limitations explicitly.
2. Prefer standard-library code for basic algorithms; introduce a dependency only when its mathematical role is justified.
3. Verify examples against independent analytic values and test invalid/boundary inputs.
4. Keep external repositories as cited references unless their licenses and provenance support reuse.
5. Keep an honest ledger of unimplemented topics. An encyclopedia of titles is not an encyclopedia of proofs.

## License and external material

This repository is newly created. A project license has **not yet been selected**; do not assume the absence of a license grants third-party reuse rights. External projects remain separately licensed and are attributed in `docs/SOURCES.md`.

