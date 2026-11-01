# Independent hostile review: Wave 10

**Date:** 2026-11-01  
**Scope:** exact polynomial limits/continuity and integration, cross-domain mathematical oracles, GitHub Actions contract, static interactive explorer, and evidence-registry integrity.  
**Review files:** `tests/test_hostile_wave10.py` and this document. The reviewer did not edit production or other specialists' files.

## Evidence and outcome

The review uses independently derived rational arithmetic and behavior at public APIs. Model-written tests are **internal executable checks**, not external accreditation or formal mathematical proofs. A Python 3.12 local test run is distinct from a GitHub-hosted CI run.

| Surface | Independent challenge | Evidence and boundary |
| --- | --- | --- |
| Polynomial limits | 60 seeded rational polynomials, independent direct evaluation at five fractions of each reported epsilon-delta bound | Checks a claimed local Lipschitz witness; never infers arbitrary-function limits from finite samples |
| One-sided diagnostics | Jump discontinuity and removable-discontinuity counterexample where the function is deliberately not called at the center | Confirms sampling remains descriptive, not a continuity proof |
| Polynomial integration | 99 exact power-integral cases with reversed endpoints, rational constants, analytically derived substitution and integration-by-parts examples | Exact rational identities for `calc.algebra.Polynomial`, not general elementary integrability |
| Evidence registry | External mutation of nested lookup bindings and entire returned manifest, followed by fresh lookups | Detects accidental poisoning of cached catalogue evidence |
| CI workflow | Static workflow trigger, read-only token and installed-wheel smoke contract; locally parsed workflow and built wheel evidence owned by CI specialist | A source/configuration check is **not** a hosted Actions result |
| Cross-domain oracles | Independent closed forms for numerical, spectral, ODE, fractional, matrix and optimization contracts | Real local executions, limited selected mathematical and numerical inputs |
| Explorer static assets | Parsed HTML labels/unique IDs and SVG accessibility metadata, checked for expression execution, network sinks, remote scripts and CSS resources | Static source checks pass; no actual browser acceptance or assistive-technology test was available |
| Root-solver regression | Smooth \(x^2-2\) bracket, independently asserted residual and bracket containment for both solvers | Confirms both methods can produce a finite root satisfying requested default residual |

## Confirmed defects and repairs

1. **Cached evidence mutation (medium severity, high confidence):** `get_evidence()` originally returned a shallow record copy, allowing a caller to mutate `bindings` and silently alter later lookups. An analogous issue affected the public cached manifest. The coordinator switched to a private cached loader and defensive deep copies; two hostile mutation regressions now pass.
2. **Premature valid-root rejection (medium severity, high confidence):** independent root oracles revealed that bisection and safeguarded Newton could reject a valid smooth root when a bracketing tolerance stopped refinement before an achievable residual. The coordinator corrected root stopping behavior. Cross-domain oracle regressions now pass.

## Residual risks

> A large suite of model-authored tests can establish many bounded examples, but it does not certify general mathematical correctness, asymptotic convergence, or zero security risk. Independent external verification and hosted CI remain desirable.

- Polynomials are exact over rationals only within declared degree/bit budgets; the epsilon-delta construction does not cover arbitrary real functions.
- The numerical one-sided sampler does not prove limits, discontinuities or removable continuity from finitely many observations.
- The explorer is an educational static client, not a general expression evaluator or an interpreter of untrusted user formulas. Browser execution and assistive-technology acceptance require their own evidence.
- GitHub's legacy commit `statuses` and Actions check runs are distinct. Empty legacy statuses alone do **not** prove Actions never ran. The available workflow-run connector filters pull-request-triggered runs; no fresh hosted run was qualified in this review.
- Repository topics marked `implemented` mean selected test-backed operations; every manifest topic still declares `formal_proof: false`.

## Reproduction

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_hostile_wave10.py -v
PYTHONPATH=src python -m unittest discover -s tests -q
```

## Final verification checkpoint

- Hostile review tests: **15/15 pass** with Python 3.12, including security/resource checks and mathematical oracles.
- Full integrated repository tests: **862/862 pass** (`PYTHONPATH=src python -m unittest discover -s tests -q`) at the final reviewer checkpoint.
- The explorer's JavaScript passed source-level safety checks. The implementation specialist also ran a JavaScript syntax check; **a live-browser interaction test was unavailable**, so pointer events, chart rendering and assistive technology remain unqualified.
- GitHub-hosted Actions CI remains **unverified**, despite a working local CI workflow contract and an isolated wheel test reported by the CI specialist.

**Verdict:** No confirmed unresolved defect within the reviewed, bounded Wave 10 scope. This does not elevate the repository to a proof-assistant or unrestricted numerical mathematics library.
