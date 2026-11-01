# Research sources, provenance, and reuse boundaries

This is a **linked source inventory**. A citation verifies that a reference is relevant; it does **not** license copying its text, images, notebook cells, proofs, or code into `calc`. Copyright, attribution, and compatibility must be checked per proposed reuse. Unless explicitly noted below, no source code or educational material has been copied into this repository. Dates of this audit: 2026-11-01.

## Core calculus and machine-learning foundations

| ID | Primary resource | Relevant scope | Source/reuse note |
| --- | --- | --- | --- |
| S01 | [MIT OCW: Single Variable Calculus (18.01SC)](https://ocw.mit.edu/courses/18-01sc-single-variable-calculus-fall-2010/) | Limits, derivatives, integrals, series | Link/reference only; OCW content license must be checked per asset. |
| S02 | [MIT OCW: Multivariable Calculus (18.02SC)](https://ocw.mit.edu/courses/18-02sc-multivariable-calculus-fall-2010/) | Multivariable and vector calculus | Link/reference only. |
| S03 | [SymPy: Calculus documentation](https://docs.sympy.org/latest/tutorials/intro-tutorial/calculus.html) | Symbolic differentiation, integration, limits, series, finite differences | Official software manual; do not copy its code into `calc` without dependency/license assessment. |
| S04 | [Jon Krohn: ML-foundations](https://github.com/jonkrohn/ML-foundations) and [license](https://github.com/jonkrohn/ML-foundations/blob/master/LICENSE) | Linear algebra and calculus for ML, plus statistics/optimization. Notebooks 1–4 are directly relevant; [notebook 8](https://github.com/jonkrohn/ML-foundations/blob/master/notebooks/8-optimization.ipynb) adds optimization | **Verified MIT**, copyright © 2020 Jon Krohn. User-suggested source. Attribution/notice required for copying substantial code or text; **link-only here**. |
| S05 | [MIT OCW: Matrix Calculus for Machine Learning and Beyond](https://ocw.mit.edu/courses/18-s096-matrix-calculus-for-machine-learning-and-beyond-january-iap-2023/) | Matrix differentials, optimization, gradients, calculus of variations | Link/reference only. |
| S06 | [SciPy: Integration tutorial](https://docs.scipy.org/doc/scipy/tutorial/integrate.html) and [solve_ivp reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) | Quadrature, ODE solvers, numerical methodology | SciPy project states **BSD** licensing in [its FAQ](https://scipy.org/faq/); dependency is preferable to code copying. |
| S07 | [MIT Mathematics: Theory of Differential Forms (18.952 syllabus)](https://math.mit.edu/~shaoyunb/18.952/18_952_Syllabus.pdf) | Tensors, exterior forms, differential forms, integration | Course syllabus; attribution/reuse terms not assessed. |
| S08 | [SymPy: Finite-difference derivatives](https://docs.sympy.org/latest/explanation/special_topics/finite_diff_derivatives.html) | Finite-difference formulas, order/accuracy | Official conceptual/software reference. |
| S09 | [SciML: DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) | ODE, SDE, DDE, DAE, PDE numerical workflows | Study supported problem categories and tests; no code incorporated. |

## Specialized and advanced analysis

| ID | Primary resource | Relevant scope | Source/reuse note |
| --- | --- | --- | --- |
| S10 | [Review of definitions of fractional derivatives and other operators (2019)](https://doi.org/10.1016/j.jcp.2019.03.008) and [Open-access fractional formalisms review](https://www.mdpi.com/2227-7390/13/22/3643) | Caputo, Riemann–Liouville, Grünwald–Letnikov, Riesz and other definitions | Scholarly references; compare operator definitions before implementation. License of individual paper/figures not assumed. |
| S11 | [NIST Digital Library of Mathematical Functions](https://dlmf.nist.gov/) | Special functions, asymptotic expansions, transforms, reference identities | Government-hosted authoritative reference; check per-section copyright before reuse. |
| S12 | [JAX: Automatic Differentiation](https://docs.jax.dev/en/latest/automatic-differentiation.html) and [advanced AD](https://docs.jax.dev/en/latest/advanced_autodiff.html) | Forward/reverse differentiation of programs, higher derivatives | Official software documentation; algorithmic inspiration only. |
| S13 | [MIT Department of Mathematics course catalogue](https://catalog.mit.edu/subjects/18/) | Tensor, differential geometry, manifolds, Lie groups, geometric analysis | Canonical subject taxonomy, not a code source. |
| S14 | [SageMath reference manual](https://doc.sagemath.org/html/en/reference/) | Operator algebra, functional analysis, symbolic algebra, number theory | Official API reference, license not checked for copying. |
| S15 | [MIT OCW: Advanced Stochastic Processes](https://ocw.mit.edu/courses/15-070j-advanced-stochastic-processes-fall-2013/pages/syllabus/) | Measure-theoretic foundations, stochastic integration | Syllabus-level scope; not a substitute for proof-level implementation references. |
| S16 | [SageMath: Algebras](https://doc.sagemath.org/html/en/reference/algebras/index.html) and [SymPy: calculus](https://docs.sympy.org/latest/tutorials/intro-tutorial/calculus.html) | Algebraic structures and symbolic operators | Official feature/API inventories; no code imported. |
| S17 | [MIT OCW: Calculus of Variations (optimal control)](https://ocw.mit.edu/courses/16-323-principles-of-optimal-control-spring-2008/resources/lec5/) | First variation, Euler–Lagrange, optimal trajectories | Lecture notes link; source license must be verified for redistribution. |
| S18 | [MIT OCW: Convex Optimization](https://ocw.mit.edu/courses/6-253-convex-analysis-and-optimization-spring-2012/) | Convex subgradients and optimization | Link/reference only. |
| S19 | [Time Scale Calculus: multi-dose pharmacokinetic modeling (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11579191/) and [q-calculus paper](https://onlinelibrary.wiley.com/doi/full/10.1155/2018/9089865) | Delta derivatives on time scales, Jackson q-derivatives | Open papers; different operators require separate APIs. |
| S20 | [SageMath: Clifford Algebras](https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/clifford_algebra.html) | Exterior and Clifford operations | Official API/manual. |
| S21 | [Roman & Rota: The umbral calculus (1978)](https://doi.org/10.1016/0001-8708(78)90087-7) | Classical finite operator/umbral calculus | Foundational research paper; no redistribution rights inferred. |
| S22 | [Stanford Encyclopedia of Philosophy: Constructive Mathematics](https://plato.stanford.edu/entries/mathematics-constructive/) and [Continuity and Infinitesimals](https://plato.stanford.edu/entries/continuity/) | Constructive, nonstandard, and smooth-infinitesimal foundations | Scholarly conceptual reference; not an algorithm spec. |
| S23 | [MIT OCW: Complex Variables with Applications](https://ocw.mit.edu/courses/18-04-complex-variables-with-applications-spring-2018/) | Holomorphic derivatives, contour integrals, residues | Link/reference only. |
| S24 | [MIT OCW: Topics in Mathematics with Applications in Finance, Lecture 24](https://ocw.mit.edu/courses/18-642-topics-in-mathematics-with-applications-in-finance-fall-2024/resources/18642-lecture-24-version-2_mp4/) and [Lecture 25](https://ocw.mit.edu/courses/18-642-topics-in-mathematics-with-applications-in-finance-fall-2024/resources/18642-lecture-25-version-2_mp4/) | Itô integrals, stochastic calculus, SDEs | Introductory overview; Malliavin, jump, rough-path details need dedicated primary references before execution. |

## Algebra and formal systems

| ID | Primary resource | Relevant scope | Source/reuse note |
| --- | --- | --- | --- |
| S25 | [Stanford Encyclopedia: Lambda Calculus](https://plato.stanford.edu/entries/lambda-calculus/) and [Type Theory](https://plato.stanford.edu/entries/type-theory/) | λ-calculus, β-reduction, types, computation | Clarifies categorical separation from differential calculus; do not reproduce lengthy prose. |
| S26 | [MIT OCW: Algebra I (18.701)](https://ocw.mit.edu/courses/18-701-algebra-i-fall-2010/) and [MIT catalogue](https://catalog.mit.edu/subjects/18/) | Groups, fields, rings, linear transformations, Lie topics | Link/reference only. |
| S27 | [SageMath: Lie Algebras](https://doc.sagemath.org/html/en/reference/algebras/lie_algebras.html) | Finite/infinite dimensional Lie algebra families | Official API/reference, not an implementation license. |
| S28 | [MIT OCW: Computational Commutative Algebra and Algebraic Geometry](https://ocw.mit.edu/courses/18-704-seminar-in-algebra-and-number-theory-computational-commutative-algebra-and-algebraic-geometry-fall-2008/) | Polynomial ideals, Gröbner bases, varieties | Link/reference only. |
| S29 | [Lean mathlib4](https://github.com/leanprover-community/mathlib4) | Formal statements and machine-checkable proofs | Useful for independently checking future mathematical claims; no proof files copied. |
| S30 | [MIT OCW: Measure and Integration (18.125)](https://ocw.mit.edu/courses/18-125-measure-and-integration-fall-2003/) | Lebesgue measure/integration, real analysis, Lp spaces, weak/measure-based methods | Introductory/graduate curriculum; specialist distributional-derivative sources still needed. |
| S31 | [Stanford Encyclopedia: Natural Deduction](https://plato.stanford.edu/entries/natural-deduction/) | Proof calculi, inference rules and normalization | Conceptual logic reference; source rights not assumed. |
| S32 | [Milner, Parrow & Walker: A Calculus of Mobile Processes, Part I](https://publish.lfcs.inf.ed.ac.uk/reports/89/ECS-LFCS-89-85/) | Original π-calculus of mobile communicating processes | Original institutional paper; no reproduction. |
| S33 | [CMU 15-445/645: Relational Model & Algebra lecture](https://15445.courses.cs.cmu.edu/spring2026/slides/01-relationalmodel.pdf) | Relation algebra versus declarative relational calculus | Official course lecture, link-only. |
| S34 | [Ehrhard & Regnier: The differential λ-calculus](https://doi.org/10.1016/S0304-3975(03)00392-X) | Formal differential λ-calculus, computational derivatives | Original research article; do not conflate with ordinary derivative API. |

## Source-first contribution policy

1. Prefer original explanations, short original examples, independent test oracles, and links to maintained primary sources.
2. **Dependency ≠ copy permission.** If importing third-party code, document version, SPDX identifier, exact license URL, copyright notices, transitive dependencies, and intended distribution mode.
3. A repository's license applies only to content actually covered by it. External notebook media and linked textbooks may have separate rights.
4. For the user-supplied [ML-foundations](https://github.com/jonkrohn/ML-foundations), the first four notebooks and optimization notebook form a valuable learning progression. Use them as a *reference spine* and use independently derived examples in `calc`. It is not a catalog of all modern calculus or algebra.
5. Technical claims need numerical tests, symbolic identities under stated assumptions, or formal proofs; citation to a lesson is evidence of a definition, not evidence that `calc` correctly implements it.
6. New source submissions should record stable URL/DOI, author, precise scope, last review date, license as **verified/unknown**, which concepts they support, and whether content was copied/adapted or merely linked.

### Research gaps that must remain explicitly open

- Dedicated foundational references for Malliavin calculus, rough paths, stochastic PDE, higher categorical and homotopical algebra, synthetic differential geometry, and generalized stochastic integrals.
- Verification of reuse terms for all OCW and SageMath assets before any copying; current entries are merely links.
- Independent formal or numerical witnesses for advanced operations in the `calc` implementation. The taxonomy's presence is not an executable capability claim.
