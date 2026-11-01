# Independent hostile review: Wave 7

**Date:** 2026-11-01  
**Target:** `clifford`, `lie`, `fem`, `finite_volumes`, `distributions`  
**Reviewer ownership:** only `tests/test_hostile_wave7.py` and this report; no production modules, registry records, GitHub refs, or external dependencies were modified.

## Result

**No confirmed defect in the tested scope.** The review adds 13 tests that independently derive expected results from defining mathematical identities, exact closed forms, and finite-volume recurrences. This is a bounded **partial** verification, not proof of the complete algebraic or analytic theories. All 13 focused checks pass; the complete shared suite of 654 tests also passed at the review checkpoint. These counts can change as further waves land.

Run:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_hostile_wave7.py -v
PYTHONPATH=src python -m unittest discover -s tests -q
```

## Independent mathematical challenges

| Target | Independent oracle | Results and limits |
| --- | --- | --- |
| Clifford/geometric algebra | Directly reorder generator words under `e_i e_j=-e_j e_i` for `i != j` and `e_i²=g_ii`, then compare every pair of 16 blades under three distinct 4D signatures. Independently sum signed coefficient squares for `scalar(reverse(X)*X)`. | 768 blade-product checks and an exact mixed-signature rational quadratic witness pass. Supports only signed orthogonal metrics and up to eight dimensions. |
| Lie algebra | Derive the six nonzero `sl₂` basis brackets from the standard `H,E,F` relations, independently test the three cyclic `so₃` brackets, reconstruct all constants, and verify `ad_[X,Y]=[ad_X,ad_Y]` on exact rational combinations. | Exact structure constants, coordinate orientation, and adjoint representation agree. Does not generalize to classification of Lie algebras or unbounded matrices. |
| One-dimensional finite elements | Solve `-u''=q` analytically as `u(x)=q(x-a)(b-x)/2` and compare all FEM nodes across 48 randomized nonuniform rational meshes. For 35 more random piecewise-constant loads, independently integrate twice using the exact Green-function boundary condition. Check the weak energy/load equality. | All exact nodal values, zero Galerkin residuals and weak identities pass. This covers homogeneous Dirichlet 1D Poisson with piecewise-constant loads, not arbitrary PDEs, functions or dimensions. |
| Periodic finite volumes | Build four upwind steps from the independent rational convex recurrence with Courant number `1/4`, both signed velocities, 72 seeded initial states. Check mass, nonincreasing periodic total variation, and maximum principle for five velocities. | All finite-grid conservation, dyadic state, and monotonicity checks pass. This does not show convergence for arbitrary discontinuous solutions, nonuniform meshes or nonlinear fluxes. |
| Distributional derivatives | Construct jumps as independently evaluated `f(x+)−f(x−)` on 48 seeded nonoverlapping rational step functions; compare atomic weights and total mass. Independently evaluate a signed four-atom measure against a polynomial, checking weak integration by parts. | All exact jump locations, strengths, pairing and cancellation witnesses pass. Polynomial pairings have a compact smooth cutoff interpretation near support; polynomials alone are not Schwartz test functions. |

## Negative and engineering boundary review

- The implementation and original specialist tests provide explicit finite-size, type, and domain checks. This independent pass stressed positive mathematical invariants more heavily than every invalid-input combination; the latter remain an evidence gap.
- `CliffordAlgebra` requires unit diagonal generator signatures, while `Multivector` rejects floating coefficients and mismatched algebra identities. The sign of a Clifford quadratic form is not necessarily positive.
- The Lie module works over exact rational matrix representations and must verify closure before declaring structure constants. Its finite basis examples are not generic Lie algebra realizations.
- FEM arithmetic is exact over rationals but resource-bounded. An exact rational nodal result does not imply a global continuous-solution error certificate between nodes.
- Upwind advection conserves mass algebraically. Its returned IEEE-754 floats need not preserve an exact real-valued invariant bit-for-bit for all representable inputs. CFL validation and numerical failure cases are tested by original module tests; only the independently derived dyadic cases are asserted exactly here.
- Distributional differentiation is demonstrated only on compactly supported rational step functions. A signed Dirac measure is a finite, atomic distribution; the module does not implement arbitrary distributions, weak Sobolev spaces, or manifolds.

## Source and license checks

The five topic guides cite mathematical references from mathematical lecture notes, MIT OpenCourseWare, Cambridge/LeVeque, FEniCS and Clawpack, and the PyGAE `clifford` project. Those references are linked for background and cross-reading, not treated as an automatic source-code reuse license. No third-party source files or licensed snippets were copied as part of this review. External website availability, all license terms, and the full mathematical claims in every linked work have **not** been independently verified here.

**Verification verdict:** Pass for the enumerated exact examples and independently tested identities; partial for full engineering validation and arbitrary-input correctness. No confirmed unresolved finding in the scope examined.
