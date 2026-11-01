# Paradigm branch map

The upstream `calc` description promises a distinct Git branch for each paradigm. The shared `main` branch should remain the canonical tested core and the navigable catalogue; `topic/<slug>` branches are **focused research and teaching lenses**, not independent incompatible copies of the library. Their creation, current head, and CI status must be read back from GitHub. This document is a publication contract, not evidence that any remote branch exists.

## Branch rules

1. Start each `topic/<slug>` from one recorded, exact `main` commit. Do not create a branch from another unreviewed subject branch.
2. Put one named paradigm's original guide, mathematical definitions, worked examples, counterexamples, references, and implementation status on that branch. Other topics stay on `main` as linked neighbours.
3. A topic branch may add subject-specific algorithms or tests, but must not silently change shared mathematical contracts. Changes to shared code go through a normal integrated review back into `main`.
4. Record mathematical assumptions and whether results are exact, numerical approximations, or heuristic illustrations. A page without a solver must say so.
5. Run `PYTHONPATH=src python -m unittest discover -s tests -v` on its actual branch head; check the branch's exact GitHub head and workflow status before marking qualified.
6. Link and cite sources without copying upstream notebooks/text. License checks are separate from mathematical review.
7. Retain one stable slug per subject. If a classification changes, keep redirects or an index mapping rather than silently repurposing its identifier.

## First subject guides

Each guide below is an original, source-linked starting point in `docs/topics/`. The suggested branch names are intentionally predictable. Publication status is **not verified here**.

| Intended branch | Guide | Distinguishing contract |
| --- | --- | --- |
| `topic/differential` | [Differential calculus](topics/differential.md) | Rates of change of real functions |
| `topic/integral` | [Integral calculus](topics/integral.md) | Accumulation and integrability |
| `topic/complex-analysis` | [Complex calculus](topics/complex-analysis.md) | Holomorphicity and contour integration |
| `topic/multivariable` | [Multivariable calculus](topics/multivariable.md) | Partial/total derivatives |
| `topic/vector` | [Vector calculus](topics/vector.md) | Divergence/curl and fields |
| `topic/fractional` | [Fractional calculus](topics/fractional.md) | Noninteger-order operators |
| `topic/stochastic` | [Stochastic calculus](topics/stochastic.md) | Itô integrals and quadratic variation |
| `topic/calculus-of-variations` | [Variational calculus](topics/calculus-of-variations.md) | Stationarity of functionals |
| `topic/differential-geometry` | [Differential geometry](topics/differential-geometry.md) | Coordinate-invariant differentiation |
| `topic/exterior` | [Exterior calculus](topics/exterior.md) | Differential forms and wedge |
| `topic/discrete` | [Discrete calculus](topics/discrete.md) | Differences and sums |
| `topic/q-calculus` | [q-calculus](topics/q-calculus.md) | Multiplicative shift derivatives |
| `topic/time-scales` | [Time-scale calculus](topics/time-scales.md) | Continuous/discrete unification |
| `topic/lambda-calculus` | [Lambda calculus](topics/lambda-calculus.md) | Reduction of computational terms |
| `topic/linear-algebra` | [Linear algebra](topics/linear-algebra.md) | Linear maps and matrix invariants |
| `topic/abstract-algebra` | [Abstract algebra](topics/abstract-algebra.md) | Axiomatic algebraic structures |
| `topic/lie-algebra` | [Lie algebra](topics/lie-algebra.md) | Brackets and Jacobi identity |
| `topic/tropical-algebra` | [Tropical algebra](topics/tropical-algebra.md) | Min-plus semiring |
| `topic/boolean-algebra` | [Boolean algebra](topics/boolean-algebra.md) | Logical complement and lattice laws |

The full taxonomy in [CATALOG.md](CATALOG.md) names many more paradigms. Those require their own subject-level sourcing, worked example, counterexample and verification before an empty symbolic branch name could mean anything useful. Branch existence is not a substitute for content.

## Merge and audit criteria

For any branch proposed as complete, the integrator should be able to produce: branch ref and commit SHA, one original worked example, one falsifying edge case, exact references, declared license status, implementation-to-theory mapping, actual test command and exit code, and a human-readable unresolved-gap section. A numerical sample cannot prove a general analytic theorem; a code review is not a formal proof.

