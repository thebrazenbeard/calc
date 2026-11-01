# Independent hostile review: Wave 6

**Date:** 2026-11-01  
**Scope:** `forms`, `nonsmooth`, `discrete`, `spectral`, `stochastic_conventions`, and the corresponding topic-manifest bindings.  
**Owner:** independent Wave 7 review agent; the implementation modules were not edited by this reviewer.

This document records tested examples, defects, and limits; it is **not** a formal proof of the five theories. All tests use deterministic inputs, fixed seeds, and independently constructed expectations. Source code in other repositories has not been copied into this review.

**Review outcome:** 18/18 targeted adversarial tests pass after repairs. The reviewer reproduced three concrete defect classes (drift-conversion overflow, unbounded Hodge dimension, spectral dynamic-range underflow); all now have passing regressions. There is no confirmed unresolved defect inside the tested scope.

## Reproduction

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_hostile_wave6.py -v
PYTHONPATH=src python -m unittest discover -s tests -q
```

## Independent witnesses

| Area | Independent oracle / challenge | Evidence | Confidence |
| --- | --- | --- | --- |
| Exterior algebra | Exact rational pairing identity `alpha ∧ ⋆ beta = <alpha,beta> vol`, orientation reversal, repeated Hodge star, graded commutativity, wedge associativity | Seeded nontrivial forms across dimensions 1–5, all degrees | High within bounded finite-dimensional Euclidean scope |
| Nonsmooth convex calculus | Direct maximum-of-affine construction, supporting-line inequality, endpoint subgradients, one-sided directional derivatives | 90 reproducible randomly seeded affine families with rational coefficients | High for scalar finite convex maxima |
| Proximal calculus | Subgradient stationarity and objective comparison against independently enumerated rational candidates; Moreau envelope recovered through minimizing point | Absolute and hinge examples with exact fractions | High for specified closed-form cases |
| Discrete sums | Newton finite differences of `k^degree`, Faulhaber telescoping up to `n=10^12`, summation-by-parts, discrete fundamental theorem | Exact integer/rational identities, fixed-seed sequences | High for bounded discrete inputs |
| Symmetric eigendecomposition | Closed-form quadratic eigenvalues, determinant and trace identities, Rayleigh extrema, eigenpair residuals, orthonormality | 70 seeded 2×2 and 45 seeded 3×3 symmetric matrices | High for moderate finite matrix scales; *not* an arbitrary-conditioning guarantee |
| Itô/Stratonovich | Analytic state-dependent drift conversion and inverse; pathwise rational midpoint/left-sum/quadratic-variation identities | Dyadic increments against `Fraction`, 45 seeded paths; constant-coefficient stochastic Heun | High for finite-grid algebra; no continuous-process convergence proof |
| Evidence registry | Imported example-test bindings and absence of claimed formal proof | Four new promoted subject records | High for local manifest consistency |

## Confirmed defects and resolution

1. **False overflow in scalar drift conversion — reproduced and fixed.** At `t=y=0`, with Stratonovich drift `-1e308`, diffusion `b(y)=1e308+4y`, and analytic `db/dy=4`, the full Itô drift is finite, `+1e308`, despite the intermediate correction `0.5*b*db/dy` exceeding binary64 range. The original implementation rejected this input. The integration agent added a rare exact-binary64 rational fallback for the full expression; the hostile test now checks both forward and inverse conversion. **Severity:** moderate. **Confidence:** high. No claim is made for arbitrary extended-precision stochastic coefficients.

2. **Unbounded exterior-form dimensions — reproduced and fixed.** A constant rational zero-form in a dimension such as one million was accepted, but Hodge complementation would construct a million-index basis. The integration agent capped ambient dimensions at 16 and sparse coefficient counts. The hostile test rejects the original adversarial construction. **Severity:** moderate availability/resource risk. **Confidence:** high.

3. **Spectral dynamic-range underflow — independently reproduced and fixed.** The original Jacobi solver returned `0` for the smaller eigenvalue of exactly diagonal `diag(1e-300, 1e300)` after uniformly normalizing by the largest entry. The first-coordinate Rayleigh quotient likewise returned `0` instead of `1e-300`. These were **silently incorrect finite results**, even though both true outputs were representable. The integration agent added an exact diagonal fast path and a rational fallback for a prematurely zeroed Rayleigh quotient. The hostile regression now verifies the repaired values. **Severity:** moderate mathematical correctness for extreme dynamic range. **Confidence:** high.

## Unresolved scope and uncertainty

- A finite-dimensional exterior algebra is **not** differential forms on arbitrary smooth manifolds or a full exterior derivative. The independent identities check its claimed constant-coefficient Euclidean operations only.
- Clarke generalized gradients are only checked where they coincide with subgradients of convex finite max-affine functions. No arbitrary nonsmooth/nonconvex calculus is implemented.
- Seeded tests cover small real symmetric matrices, not arbitrary dense spectral PDE methods or robust eigendecomposition under all extreme condition numbers. Numerical residuals alone are diagnostics, not certificates.
- The Brownian tests are exact **finite-grid identities** with caller-supplied increments. They do not prove convergence in probability or pathwise convergence of numerical SDE schemes.
- The registry's `implemented` label means selected tested operations, not an entire theory. It does not certify completeness, upstream license grants, or formal proofs.
- This review adds no external source code. Background mathematical links remain in the adjacent individual topic documents and central source index; those links alone do not confer reuse permission.
