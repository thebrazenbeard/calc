# Independent hostile review: Wave 9

**Date:** 2026-11-01  
**Reviewer:** separate Wave 9 specialist, independently authored tests and findings.  
**Target:** `calc.measure`, `calc.free_algebra`, `calc.operator_algebra`, `calc.control`, and `calc.differential_algebra`.

## Verdict

**No confirmed unresolved defect within the exercised scope.** Fourteen independent adversarial tests passed. The complete shared repository suite also passed **757/757** tests on Python 3.12 at the final review checkpoint. These are examples and mathematical contract checks, not proofs of all inputs or all of the named fields. The five modules deliberately implement finite/bounded subsets of their larger subjects.

## Reproduction

From the project root:

```bash
PYTHONPATH=src python -m unittest tests.test_hostile_wave9 -v
PYTHONPATH=src python -m unittest discover -s tests -q
```

## Requirement and evidence matrix

| Target | Independent mathematical oracle / hostile probe | Result | Confidence and limits |
| --- | --- | --- | --- |
| Finite measures | 60 deterministic seeded Radon–Nikodym density examples, including null reference atoms; exact change-of-measure integration identity | Pass | High for nonnegative rational finite atom measures; no general sigma-algebra or continuum measure theorem |
| Conditional expectation | 75 seeded exact rational partitions; integrals preserved within **each** cell and across the entire space; zero-mass cell convention, duplicate/overlap rejection | Pass | High for finite partition conditional means, not arbitrary conditional expectation |
| Free noncommutative polynomial algebra | Explicit noncommuting word products and a custom generator-image derivation, 55 seeded exact Leibniz and Jacobi identities, invalid generators and powers | Pass | High for small bounded words over Q; no general free-algebra completion or noncommutative Gröbner bases |
| Finite rational matrix operators | 45 seeded exact transpose-of-product and trace-cyclicity identities; 40 rank-one positive-semidefinite constructions each tested on a 49-point rational vector grid; diagonal star-algebra closure and nonunital rejection | Pass | High for the small real rational matrices and predicates tested; not complex/infinite-dimensional C*-algebra certification |
| Scalar finite-horizon LQR | 65 seeded independently completed-square one-step solutions, comparing gain and value against closed formulas and arbitrary nearby controls; multi-stage Bellman value decomposition; zero-horizon boundary | Pass | High for deterministic rational scalar LQR with nonnegative state penalties and strictly positive control penalty; no general control theory |
| Differential polynomial algebra | Explicit formal chain and Leibniz rules; 65 seeded independent monomial-total-derivative evaluation oracles, 55 product Leibniz identities, top-jet truncation and type rejection | Pass | High for the finite bounded jet ring where the formal total derivative is representable; no analytic ODE solution |

## Review method and source handling

- The expected values came from direct exact-rational formulas or mathematical identities, not from implementation helper calls computing expected results.
- The review intentionally probes zero masses, singular absolute continuity, noncommutation, nonunital spans, invalid generator names, finite-horizon boundaries, and missing highest-order jets.
- Inputs and expected values are deterministic. The random seed is fixed for reproducibility.
- Source documents for each module link to mathematical literature or university lecture notes. These are references, not a grant of code/document reuse rights; no third-party source code was copied into the review.
- This reviewer did **not** modify production modules, shared build configuration, the topic manifest, or any GitHub references.

## Unresolved coverage and risk boundaries

1. **Finite versus infinite objects:** finite atoms, bounded words, small real matrices and truncated jet variables cannot establish all results about general measure spaces, noncommutative completions, operator norms or differential fields.
2. **Resource budgets:** several operations bound dimensions, coefficient bits and term counts. Large valid mathematical objects may be intentionally rejected; that is part of the specified finite educational scope.
3. **General optimal control:** unconstrained scalar deterministic discrete LQR is a special case, not arbitrary nonlinear, constrained, stochastic, or continuous-time optimal control.
4. **Formal verification:** independent algebraic property tests are useful regression witnesses but are not machine-checked proofs for every value in an unbounded domain.
5. **Cloud execution:** the reported 757/757 total is a local Python test result. This review did not verify a green GitHub Actions run or any remote deployment.

The review outcome is **pass for the expressly bounded examples**, with the general-theory gaps above left explicitly open.
