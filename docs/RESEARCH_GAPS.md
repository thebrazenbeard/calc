# Independent research and coverage audit — 2026-11-01

**Subject:** `thebrazenbeard/calc` local source after Wave 3, the public repository description, its 19 topic branches, `docs/CATALOG.md`, `docs/SOURCES.md`, `docs/BRANCHES.md`, and the user-supplied `jonkrohn/ML-foundations`. This is an independent internal subagent review of the stated objectives, source suitability and branch architecture; it is **not** a formal mathematical proof, legal opinion, external model audit, or comprehensive live GitHub CI qualification.

**Verdict: PARTIAL, with a useful working foundation.** The current work supports a broad *research index* and a growing set of demonstrators. It does not support the repository description's literal claim of an “exhaustive” archive covering “every known branch” with every paradigm in its own branch. Those terms have no finite, generally accepted mathematical boundary.

## Inventory and evidence

| Observation as of this audit | Evidence class | Result |
| --- | --- | --- |
| **78** indexed slugs across seven families | Executed `calc.catalog.TOPICS` inspection | Confirmed locally |
| **31** slugs labeled `implemented`, **47** `documented` | Executed registry inspection | Confirmed locally; “implemented” explicitly means *some executable example*, not complete theory |
| **19** original subject guides and 19 matching intended `topic/<slug>` branches | File inventory and branch-map inspection | Confirmed local guides; three branch `TOPIC.md` files independently read back from GitHub |
| Parent integrator's publication record reports **19/19** branch file readbacks | `docs/PUBLICATION.md` source record | Reported evidence; this reviewer sampled three, not all nineteen |
| **34** bibliography identifiers, many containing multiple URLs | `docs/SOURCES.md` inspection | Confirmed source inventory; not all links/proofs/license terms independently audited |
| `jonkrohn/ML-foundations` at `master` has an MIT license, copyright © 2020 Jon Krohn | Direct GitHub `LICENSE` content readback | Confirmed, source blob `fef341762bd4da13cc23ceba76d15fbbe453fce1` |
| The current repository does not declare its own license | Public GitHub metadata and README source | Confirmed at review; licensing decision remains Patrick's |
| Catalogue unit tests run | Python 3.12 `PYTHONPATH=src python -m unittest discover -s tests -p test_catalog.py -v` | 6/6 PASS, exit 0 |

The other specialists' full-unit-test and wheel-build receipts are separate execution evidence. This review does not elevate them into formal verification. A topic's existence in the index does not imply that an algorithm has been implemented or qualified at every boundary.

## Material falsifiers and gaps

### 1. Public claim exceeds substantiated coverage — Medium severity, high confidence

The GitHub repository description asserts an “exhaustive” archive of *every known branch* and branch isolation for *each paradigm*. Source inspection finds **78** indexed topics and only **19** curated topic branches. The main README correctly warns that no finite enumeration can establish mathematical exhaustivity. These two public descriptions conflict.

> **Hostile challenge:** If a new area of noncommutative geometry, rough integration, synthetic mathematics or a novel typed process calculus appears tomorrow, what theorem or authoritative classification would certify this archive still contains “every” calculus? There is no such supplied completeness criterion.

**Remedy direction:** Retain the ambitious continuously expanding objective while making the public promise falsifiable: “an extensible, independently sourced atlas of mathematical calculi and algebras.” Keep an explicit `INDEXED / GUIDE / DEMONSTRATED / TESTED / FORMAL_PROOF` evidence ladder. Do not reinterpret a module-level example as whole-field qualification.

### 2. Status and guide coverage do not line up — Medium severity, high confidence

The 31 `implemented` slugs do not map one-to-one to worked subject guides. **16 implemented slugs currently have no `docs/topics/<slug>.md` guide**:

`differential-forms`, `ordinary-differential-equations`, `partial-differential-equations`, `ito`, `finite-difference`, `automatic-differentiation`, `numerical-differentiation`, `numerical-integration`, `integral-transforms`, `symbolic-calculus`, `root-finding`, `numerical-optimization`, `polynomial`, `matrix`, `group-theory`, `field-theory`.

Conversely, **four guides have `documented` status**: `time-scales`, `integral`, `lie-algebra`, `vector`. For `vector`, `differential.py` does have divergence and curl examples: this is a likely registry classification omission rather than a missing numerical implementation. For `integral`, numerical quadrature exists, but that does not implement general exact/measure-theoretic integration.

**Remedy direction:** Add a machine-readable topic manifest with `concept_kind` (theory / operation / numerical method / algebraic structure / formal calculus), `guide_path`, `module_api`, `test_ids`, `scope_statement`, `primary_source_ids`, `evidence_level`, and `branch_head`. Generate status tables from this manifest; do not maintain mismatched prose and hand-written status sets.

### 3. Branch isolation is not the same as a qualified branch — Medium severity, high confidence

I independently read `TOPIC.md` from `topic/time-scales`, `topic/lambda-calculus`, and `topic/abstract-algebra`; each points to the exact Wave 2 base `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`. That is good provenance. However Wave 3 added new modules on `main` after those branches split, including `lambda_calc.py` and `complex_calc.py`. **Direct GitHub readback confirmed** `src/calc/lambda_calc.py` exists on `main` while the same path on `topic/lambda-calculus` returns 404. A guide branch therefore does **not** automatically contain the corresponding newest implementation. Parent readback of all 19 subject files is strong evidence of remote file presence; no branch-specific CI readback or convergence proof was supplied.

**Remedy direction:** For each maintained branch, record its exact common-base SHA, current head, the matching API artifact SHA, and whether the API runs on that branch. Use normal reviewed forward integrations from `main` when appropriate. A subject branch should distinguish “research only” from “runnable here.”

### 4. Reference mapping is broad but sometimes indirect — Medium severity, high confidence

The bibliography is valuable as a **discovery index**, but several citations are not direct source specifications:

- `S07` is a differential-forms course **syllabus**, weak support for a specific tensor-covariant-derivative rule.
- `S13` is an MIT **course catalogue**; it proves that subjects are taught, not any particular formula.
- `S15` is a stochastic-processes **syllabus**, not a specialist Malliavin/rough-path or measure-theory proof.
- `S19` includes a pharmacokinetic modeling paper; it demonstrates a time-scale application, but the general delta/nabla calculus should be sourced from a foundational text.
- `S22` is strong philosophical background on constructive/infinitesimal approaches, but synthetic differential geometry needs a dedicated mathematical model reference.
- `S24` is an introductory Itô/SDE lecture source. It explicitly cannot qualify Malliavin calculus, Lévy–Itô jumps, rough paths, or regularity structures.
- `S14` is a broad SageMath API collection, too indirect for a general spectral/holomorphic functional-calculus theorem.

The source file itself acknowledges many of these limitations. Do not silently upgrade them to exact mathematical specifications.

### 5. External licensing is not a single interchangeable permission — Medium severity, high confidence

The cited ML-foundations code is **MIT-licensed**, but copied code would still require the relevant copyright/license notice. MIT OpenCourseWare's [current terms](https://ocw.mit.edu/pages/privacy-and-terms-of-use/) describe **CC BY-NC-SA 4.0** for its covered content; the terms are not interchangeable with MIT, BSD or the unlicensed state of `calc`. [SageMath's licensing documentation](https://doc.sagemath.org/html/en/reference/history_and_license/index.html) describes GPL terms; [SymPy's official licensing explanation](https://docs.sympy.org/dev/citing.html) describes BSD terms. Cite and link these resources; do not bulk-import notebooks, course PDFs or external source files without a per-asset decision.

**Remedy direction:** Patrick should choose a license for original `calc` material. If any upstream content is actually copied or adapted, create a provenance manifest with exact upstream revision, file/asset, SPDX license, notices, transformation, redistribution restrictions, and compatibility review. No copied upstream content was evident in the audited original guides.

## Priority-ordered next research/build waves

Priorities weigh **mathematical centrality**, **cross-domain usefulness**, **falsifiability**, and **implementation feasibility**. They are engineering recommendations, not claims that higher topics are unimportant.

| Priority | Next slice | Observable acceptance / hostile counterexample | Primary starting source |
| --- | --- | --- | --- |
| **P0** | Single machine-readable taxonomy/provenance manifest | 78 registry entries reconcile with guide, source, module, test, branch evidence; unverified never equals implemented | Repository's exact `catalog.py`, `docs/CATALOG.md` |
| **P0** | Branch and CI qualification | Exact head of each curated branch; branch-local import/test; no stale main/branch implementation implication | GitHub branch refs and workflow readback |
| **P0** | Original-work licensing decision and reuse ledger | Explicit project license chosen; imported assets tracked individually | [GitHub licensing guidance](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository) |
| **P1** | Linear algebra and multivariable matrix calculus | Solve/rank/nullspace, eigen/spectral residuals, conditioning, Jacobian of matrix maps with independent analytic witnesses | [MIT 18.06SC](https://ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/), [MIT 18.S096](https://ocw.mit.edu/courses/18-s096-matrix-calculus-for-machine-learning-and-beyond-january-iap-2023/) |
| **P1** | Full discrete/difference operator family | Forward/backward/central, summation by parts, step-size limit, q→1 boundary; explicit singular handling | [SymPy finite differences](https://docs.sympy.org/latest/explanation/special_topics/finite_diff_derivatives.html), [Bohner–Peterson time scales](https://www.math.unl.edu/~apeterson1/sample_book.pdf) |
| **P1** | Symbolic and numerical integrals with assumptions | Safe expression differentiation/integration for a bounded grammar; improper/discontinuous counterexamples; certified-vs-heuristic flags | [SymPy integral documentation](https://docs.sympy.org/latest/modules/integrals/integrals.html), [SciPy quadrature](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html) |
| **P1** | Computational commutative algebra | Multivariate polynomial order; Gröbner basis checked by S-polynomial reduction and ideal-membership examples | [Sage multivariate ideals](https://doc.sagemath.org/html/en/reference/polynomial_rings/sage/rings/polynomial/multi_polynomial_ideal.html) |
| **P2** | Differential forms and tensor geometry | Coordinate covariance and pullback; `d²=0`, torsion-free metric connection; degeneracy/singularity rejection | [MIT differential forms](https://math.mit.edu/~shaoyunb/18.952/18_952_Syllabus.pdf) plus a proof-level textbook |
| **P2** | Fractional operators with explicit conventions | Riemann–Liouville vs Caputo constant derivative difference, endpoint singularity, semigroup assumptions | [Fractional review DOI](https://doi.org/10.1016/j.jcp.2019.03.008) |
| **P2** | Itô vs Stratonovich and stochastic error | Same-noise path coupling; geometric-Brownian error convergence; Itô chain-rule correction under declared assumptions | [Nualart stochastic/Malliavin notes](https://nualart.ku.edu/Course.pdf) |
| **P2** | Nonsmooth/functional calculus | Classical vs Clarke/convex subgradient at `|x|`; Fréchet vs Gâteaux counterexample | [Clarke generalized-gradient paper](https://doi.org/10.1287/moor.1.2.165) |
| **P3** | Malliavin derivatives and jump processes | Wiener functional derivative with domain and Sobolev space declared; Poisson compensator not silently equated to Brownian noise | [Nualart course](https://nualart.ku.edu/Course.pdf) |
| **P3** | Rough integration and singular SPDE | Explicit path regularity, iterated-integral/signature hypotheses; no Euler/standard Itô shortcut | [Lyons foundational paper](https://doi.org/10.4171/rmi/240), [Hairer regularity structures](https://arxiv.org/abs/1303.5113) |
| **P3** | Synthetic foundations and higher algebra | Distinguish hyperreal from nilpotent infinitesimals; specify ambient category/model before calculation | [Anders Kock, *Synthetic Differential Geometry*](https://users-math.au.dk/kock/sdg99.pdf), [Lurie, *Higher Algebra*](https://www.math.ias.edu/~lurie/papers/HA.pdf) |
| **P3** | Formal calculi/proof verification | Typed λ reduction and capture avoidance; relational/process logic separate from numeric differentiation; machine-checked example when feasible | [SEP lambda calculus](https://plato.stanford.edu/entries/lambda-calculus/), [mathlib4](https://github.com/leanprover-community/mathlib4) |

### How to use `jonkrohn/ML-foundations` next

It is a practical **learning-path reference**, not a wholesale implementation plan. Anchor beginner routes to its linear algebra and calculus notebooks (the first four) and optimization notebook 8, then cross-check definitions against MIT/SymPy/SciPy and implement **independently written** examples with exact test oracles. Cite the [repository](https://github.com/jonkrohn/ML-foundations) and its [MIT license](https://github.com/jonkrohn/ML-foundations/blob/master/LICENSE); preserve upstream notices if any actual material is imported. The notebook selection was checked against the repository's documented subject coverage, not audited cell by cell here.

## Evidence and limits

- **Direct tool evidence:** local catalog import/count and 6/6 catalog tests; GitHub remote `ML-foundations/LICENSE`; three sampled topic-branch `TOPIC.md` readbacks; positive `main` and negative `topic/lambda-calculus` reads for `src/calc/lambda_calc.py`; public `calc` metadata; official/source paper pages listed above.
- **Source inspection:** 19 local guides, numbered bibliography, project README, branch policy and publication receipts.
- **Not established by this reviewer:** a complete 19-branch CI matrix; a license selected for `calc`; per-resource terms for every of the 34 source IDs; a proof assistant certification; whole-field correctness for any advanced calculus.
- **No mutations outside this report:** no code, references, branch refs, project metadata or licensing altered.

**Conclusion:** The useful next move is not simply more named folders. Make the 78-topic evidence map machine-readable, qualify which mathematical contract each module actually solves, and promote topics to implementation/proof status only with a falsifiable witness. That lets the atlas expand aggressively without confusing breadth of labels with depth of mathematics.
