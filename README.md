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

The CLI shows the subject map without executing arbitrary mathematical expressions:

```bash
PYTHONPATH=src python -m calc topics
PYTHONPATH=src python -m calc topic fractional
PYTHONPATH=src python -m unittest discover -s tests -v
```

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

