# Independent hostile review — `thebrazenbeard/calc`

Review date: 2026-11-01. Scope: the initial repository build, its scope claims,
mathematical implementations, source licensing, and promised topic-branch layout.
This is a model review, not a formal proof or independent external audit.

## Contract and evidence labels

The user's objective is a deeply populated reference and executable toolkit for
calculus **and** algebra, informed by public repos and external scholarly sources.
The GitHub repository description additionally says each calculus paradigm is
isolated in its **own branch**. These are distinct acceptance conditions.

Evidence classes used below: **RETRIEVED** (live GitHub source/metadata),
**SOURCE_INSPECTION** (local code), **TEST** (executed Python behavior),
**INFERENCE** (reasoned concern), **UNKNOWN** (not independently established).

| Requirement or risk | Best validation | Current finding |
| --- | --- | --- |
| Breadth without pretending literal exhaustivity | Versioned taxonomy with separately tracked covered/planned areas | To inspect after implementation |
| A topic per branch | `git branch -r` or GitHub branch list, checking actual topic content | **OPEN** — initial empty repo has no topic branches |
| Correct algebra and analysis routines | Independent analytic oracles, negative cases, deterministic identities | **PARTIAL** — 7 extra contract tests passed; one original numerical test failed |
| Reproducible installation | Clean `pip install .` and import from outside source tree | To test after implementation |
| Source provenance and rights | Original LICENSE files and attribution/reuse manifest | Partial: Jon Krohn MIT confirmed |
| Clear supported domains | API docstrings and examples naming assumptions and error bounds | To inspect |

## Material falsification targets

1. **High — Scope/representation mismatch (RETRIEVED).** The repository
   description promises independently isolated topic branches, not just directories
   under `main`. A single populated `main` branch cannot verify that promise.
   Distinguish a common core (`main`) from catalogued topics, and publish
   `topic/<slug>` branches pointing to real topic-specific content, or explicitly
   amend the promise with Patrick's consent. The architecture should not quietly
   create dozens of duplicative, diverging Python implementations.
2. **High — Mathematically incorrect results are worse than missing methods
   (INFERENCE).** Numerical code needs analytical cross-checks at normal and
   adversarial scales: zero-length intervals, negative orientation, discontinuities,
   near-singular matrices, singular Jacobians, vanishing step sizes, integer and
   rational coefficients, non-finite inputs. Every public routine should specify
   its domain and error behavior.
3. **High — False 'all calculi' assertion (INFERENCE).** 'Calculus' spans
   mathematically different constructions, e.g. infinitesimal calculus,
   differential geometry, discrete/time-scale, stochastic, functional,
   fractional, variational, nonstandard and the *lambda calculus* of programming
   languages. No finite list is plausibly exhaustive, and a documentation entry
   must not imply a corresponding tested implementation. Mark each topic
   `implemented`, `documented`, `source-indexed`, or `planned`.
4. **Medium — Source/license contamination (RETRIEVED).**
   `jonkrohn/ML-foundations` has a committed MIT `LICENSE` (copyright © 2020 Jon
   Krohn). Its source may be reused only with the applicable license notice in
   copies/substantial portions. Its README references notebooks on introductory
   linear algebra, limits/derivatives, partial derivatives/integrals, statistics
   and optimization; it is a useful teaching source, not authority for every
   advanced calculus. Start with citations and independently written examples.
   Do not conflate MIT licensing of code/notebooks with licensing of linked
   commercial videos, books or third-party datasets.
5. **Medium — Numerical stability and explicit guarantees (INFERENCE).**
   A finite-difference API must define step-size handling; Simpson's method
   requires an even number of panels; adaptive integration needs a work/recursion
   cap and an honest non-convergence response; an ODE solver needs step direction,
   endpoint, vector-dimension and non-finite-output checks.
   **Confirmed counterexample (TEST):** On the initial working tree, evaluating
   `adaptive_simpson(lambda x: cos(8*x), 0, pi)` returned `pi` instead of the
   exact integral `0`. Every sampled dyadic node happens to have value `1`,
   so the Richardson difference is zero while the true error is large. This
   defeated its heuristic error estimate and caused the provided oscillatory
   regression to fail. **A generic finite set of function samples cannot prove
   an error bound for an arbitrary unknown oscillatory integrand.** Possible
   mitigation: add non-dyadic cross-sampling or a frequency-aware method,
   retain the documented non-certification, and distinguish heuristic from
   guaranteed tolerance. Do not merely remove the failing test.
6. **Medium — Misleading algebra generality (INFERENCE).** A matrix solver using
   ordinary floats is not automatically valid for arbitrary rings; exact
   `Fraction` polynomial arithmetic does not imply a complete symbolic algebra
   system. The public API should avoid generic labels like 'solve everything'.
7. **Medium — Installation and test independence (INFERENCE).** A test executed
   from the source working directory might import source files by accident. A
   clean package install and test from a separate directory is stronger evidence.
   Expected test values must come from analytic identities or independently
   derived results, not a second transcription of implementation logic.
8. **Medium — Citation mapping is broader than source support (SOURCE_INSPECTION).**
   In `docs/CATALOG.md`, the real-analysis/geometric-measure rows cite S15,
   an advanced stochastic-processes syllabus, and the process-calculi row
   cites S25, a lambda-calculus/type-theory encyclopedia entry. Those URLs do
   not establish the named results in their respective rows. Likewise the
   stochastic-calculus row includes S10 (fractional derivatives, a different
   subject). Replace broad proxy citations with direct primary references
   when building those subject guides; until then these remain navigation
   leads, not evidence for mathematical definitions or implementations.

## Targeted property tests

- Algebra: polynomial distributivity, coefficient trimming, `d/dx ∫p dx=p`,
  exact rational evaluation, determinant multiplicativity for small nonsingular
  matrices, linear solve residual and singular-matrix rejection.
- Differential: derivative of `sin` = `cos`, mixed partials for a polynomial,
  Hessian symmetry for smooth test functions, gradient directional dot product,
  expected Jacobian shape and orientation.
- Numerical: Simpson integrates cubic polynomials exactly in exact arithmetic
  (within floating tolerance here), integral sign reverses with endpoint order,
  `y'=y` gives `e^t` approximately, RK4 fourth-order convergence, adaptive
  integration signals exhaustion rather than silently returning a false bound.
- Advanced: `Δ_h x²=2hx+h²`, a `q`-derivative of `x²` gives `(1+q)x` when
  admissible, fractional integral of `1` obeys `x^α/Γ(α+1)` on `x≥0`, and
  exterior products satisfy graded commutativity. Check all claimed domains.

## Directly checked external reference

- [Jon Krohn, *ML-foundations*, `README.md` on `master`](https://github.com/jonkrohn/ML-foundations/blob/master/README.md)
  — introductory ML mathematical foundations. Retrieved 2026-11-01.
- [Jon Krohn, *ML-foundations*, `LICENSE` on `master`](https://github.com/jonkrohn/ML-foundations/blob/master/LICENSE)
  — MIT, copyright 2020 Jon Krohn. Retrieved 2026-11-01.

## Review disposition

**PARTIAL, with a confirmed numerical failure / OPEN.** The contractual gap and
reference-license status are documented. I executed:

- `PYTHONPATH=calc/src python -m unittest discover -s calc/tests -p 'test_*.py' -v`
  from the workspace on 2026-11-01: **57 tests, 1 failure**, the oscillatory
  adaptive-Simpson case above; all other 56 passed. This predates adding the
  independent `test_contracts.py` file.
- `PYTHONPATH=calc/src python -m unittest discover -s calc/tests -p 'test_contracts.py' -v`:
  **7 tests, all passed**. These cover polynomial distributivity and Leibniz,
  determinant multiplicativity, brute-force CRT comparisons across moduli 1–8,
  curl-of-gradient, Simpson's cubic exactness, RK4 convergence, and exterior
  graded symmetry/associativity.

Installed distribution, actual branch topology and CI still require readback
after integration. No green status is inferred from another worker's report.
