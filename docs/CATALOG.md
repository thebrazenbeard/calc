# The calculus and algebra atlas

> Status: **research index**, not an assertion that every item below has an executable implementation. The code, API reference, and test suite determine which entries are operational. `Indexed` means a subject is identified and placed in the atlas; `planned` means implementation or a worked lesson is still needed.

There is no mathematically canonical, finite enumeration of **every type** of calculus or algebra. Both names designate families of theories, methods, and formalisms; branches overlap, and research continues to create subfields. This is an extensible atlas, not a misleading claim of completeness. It distinguishes three senses of *calculus*: analysis of change/accumulation, algebraic/combinatorial generalizations of operators, and formal calculi in logic/computer science. The last category is not ordinary differentiation.

## How to navigate

The atlas is organized by **what is being differentiated or manipulated**, not just by the label attached to a course. A theory may occupy several groups: for example, a fractional stochastic differential equation joins fractional operators, probability, and numerical analysis. See [SOURCES.md](SOURCES.md) for the numbered, independently checkable reference links and license notes. The taxonomy itself is original explanatory organization; no third-party educational content is embedded.

Suggested implementation milestones:

1. **Core:** exact elementary algebra, polynomial arithmetic, linear algebra, limits, derivatives, quadrature, finite differences, ODE integration, error estimation.
2. **Intermediate:** multivariable/vector/matrix calculus, higher-order derivatives, root finding, constrained optimization, integral transforms, symbolic differentiation, series.
3. **Advanced:** exterior/tensor/variational/fractional/stochastic/functional/distributional operators, specialized algebras, PDE discretization.
4. **Research:** nonstandard/constructive, nonlinear and noncommutative generalizations, proof-carrying statements, specialized formal calculi. Put domain assumptions and independent numerical or formal witnesses beside each implementation.

## Calculus of real and complex-valued functions

| Family | Specific paradigms and concepts to catalog | Representative question or operator | Entry references |
| --- | --- | --- | --- |
| Precalculus bridge | Functions, composition, inverse functions, trigonometry, exponentials, logarithms, sequences, asymptotic notation | Domain, codomain, continuity prerequisites | [S01], [S02] |
| Single-variable differential | Limits, continuity, derivative from first principles, product/chain/implicit rules, higher derivatives, Taylor theorem | `d f(x)/dx` | [S01], [S03], [S04] |
| Single-variable integral | Riemann/Darboux, improper integrals, fundamental theorem, substitution, integration by parts, numerical quadrature | `integral_a^b f(x) dx` | [S01], [S03], [S04] |
| Sequences and series | Taylor/Maclaurin, Fourier series, power series, asymptotic series, uniform convergence | Convergence and truncation error | [S01], [S03], [S11] |
| Multivariable | Partial/directional derivatives, gradients, Jacobian, Hessian, total differentials, iterated/multiple integrals, constrained extrema | `∇f`, `J_f`, `H_f` | [S02], [S03], [S05] |
| Vector analysis | Divergence, curl, line/surface/volume integrals, Green/Stokes/divergence theorems | `∇·F`, `∇×F` | [S02], [S05] |
| Matrix calculus | Differentials of matrix maps, trace derivatives, tensor-valued Jacobians, reverse-mode gradients | `d tr(XᵀAX)/dX` | [S05], [S12] |
| Real analysis | ε–δ foundations, absolute continuity, bounded variation, Riemann–Stieltjes integration, Lebesgue integration, dominated convergence | Limits versus integrals | [S01], [S30] |
| Complex analysis | Holomorphic/meromorphic derivatives, Cauchy–Riemann equations, contour integrals, residues, analytic continuation | `∮ f(z) dz` | [S23], [S11] |
| Infinite series and transforms | Laplace, Fourier, Mellin, z-transforms, integral equations, Green functions, special functions | Change domains to solve operators | [S11], [S17] |
| Differential equations | First/higher-order ODE, initial/boundary value problems, stiff systems, systems of ODEs, differential-algebraic equations | `y'=f(t,y)` | [S06], [S09] |
| PDEs | Elliptic/parabolic/hyperbolic equations, weak solutions, conservation laws, characteristic curves, Fourier/spectral/FEM/finite-volume methods | Heat, wave, Laplace equations | [S06], [S09] |

## Geometric, tensor, and higher-dimensional change

| Family | Specific paradigms and concepts | Central operator/object | Entry references |
| --- | --- | --- | --- |
| Tensor calculus | Tensor products, index notation, contraction, covariant derivatives, connections, curvature tensors | `∇_μ T` | [S07], [S13] |
| Differential geometry | Manifolds, tangent/cotangent spaces, pushforwards/pullbacks, Lie derivatives, geodesics | Derivative along a manifold | [S07], [S13] |
| Exterior calculus | Differential forms, wedge/exterior product, exterior derivative, interior product, generalized Stokes theorem | `d: Ω^k→Ω^(k+1)` | [S07], [S13], [S20] |
| Hodge and de Rham | Hodge star, codifferential, Laplace–de Rham, de Rham cohomology | `★`, `δ`, `d` | [S07], [S20] |
| Symplectic/contact | Symplectic forms, Poisson brackets, Hamiltonian flows, contact structures | `{f,g}` | [S07], [S13] |
| Lie calculus | Lie group flows, Lie algebra commutators, Lie derivatives, Maurer–Cartan forms | `[X,Y]`, `𝓛_X` | [S13], [S20] |
| Geometric measure | Hausdorff measures, rectifiable currents, varifolds, generalized surface area | Weak geometry | [S30] |
| Discrete differential geometry | Mesh exterior derivatives, discrete Hodge star, discrete curvature, computational geometry | Cochain complexes | [S07], [S20] |
| Geometric/Clifford calculus | Multivectors, geometric product, grade projections, multivector differentiation | Clifford geometric product | [S20] |

## Generalized derivatives, integrals, and operator calculi

| Family | Subtypes and important distinctions | Typical diagnostic / edge | Entry references |
| --- | --- | --- | --- |
| Finite-difference/discrete | Forward, backward, central differences; divided differences; sums; summation by parts | `Δ_h f(x)=f(x+h)-f(x)` (scaled derivatives divide by `h`) | [S03], [S08] |
| Time-scale calculus | Delta/nabla derivatives on closed subsets of ℝ; hybrid continuous/discrete domains | Reduces to ordinary derivatives or differences in special cases | [S19] |
| Quantum `q`-calculus | Jackson q-derivative, q-integral, q-difference equations, `(p,q)` extensions | `D_q f(x)=(f(qx)-f(x))/((q-1)x)` for allowed `q,x` | [S19] |
| Umbral/operator calculus | Shift-invariant delta operators, Sheffer/Appell sequences, operator identities | Abstracts differentiation through linear operators | [S21] |
| Fractional calculus | Riemann–Liouville, Caputo, Grünwald–Letnikov, Weyl, Marchaud, Riesz, Hadamard, Hilfer | Orders α not restricted to integers; definitions are **not interchangeable** | [S10] |
| Nonlocal calculus | Fractional Laplacian, nonlocal gradients/divergences, integral kernels | Boundary conditions affect which fractional operator is meant | [S10] |
| Functional calculus | Operator functions (holomorphic, continuous/Borel spectral calculus) and functional derivatives of maps between function spaces | Avoid conflating operator functional calculus with variational derivative | [S14], [S16] |
| Distributional/weak | Schwartz distributions, weak derivatives, Sobolev spaces, weak PDE formulations | Distributional derivative of a step is a Dirac delta | [S30] |
| Generalized/nonsmooth | Subgradients, Clarke generalized gradients, convex analysis, proximal maps | `|x|` is not classically differentiable at zero | [S18] |
| Automatic differentiation | Forward/reverse mode, Jacobian-vector and vector-Jacobian products, custom rules, higher-order AD | Chain rule through programs; not symbolic simplification or finite differences | [S12] |
| Symbolic/computer algebra | Expression trees, simplification, exact integration, symbolic limits, assumptions, rewriting | Results may be conditional or unevaluated | [S03], [S16] |
| Numerical calculus | Adaptive quadrature, differentiation, interpolation, integration error, root-finding, stability/conditioning | Error estimates must state smoothness and step-size assumptions | [S06], [S08], [S09] |
| Interval/validated calculus | Interval enclosures, automatic bounds, directed rounding, Taylor models | Certified enclosure rather than unqualified floating-point approximation | [S14] |
| Asymptotic calculus | Perturbation, multiple scales, singular perturbations, matched expansions, stationary phase, steepest descent | Approximation regimes and remainder terms | [S11] |

## Probabilistic, variational, and infinite-dimensional change

| Family | Subtypes | Critical mathematical distinction | Entry references |
| --- | --- | --- | --- |
| Stochastic calculus | Itô integral, Itô lemma, Stratonovich integral, semimartingales, stochastic differential equations | An Itô chain rule contains quadratic-variation terms | [S15], [S24] |
| Jump processes | Poisson random measures, Lévy–Itô decompositions, compensated jumps, integro-differential generators | Must specify conventions for jumps/integration | [S24] |
| Malliavin calculus | Derivative of Wiener functionals, Skorokhod divergence, integration by parts on Gaussian space | Infinite-dimensional probabilistic differentiability | [S24] |
| Rough paths/regularity | Controlled rough paths, signatures, rough integration, regularity structures | Paths can be too irregular for ordinary Riemann–Stieltjes integration | [S24] |
| Stochastic PDEs | Random forcing, stochastic heat/wave equations, mild/weak solutions | Noise models and spaces are indispensable | [S09], [S24] |
| Variational calculus | First/second variation, Euler–Lagrange equations, natural boundary conditions, constrained/isoperimetric variations | `δJ[y]=0` for admissible perturbations | [S05], [S17] |
| Optimal control | Pontryagin maximum principle, Hamilton–Jacobi–Bellman, dynamic programming | Optimize trajectories subject to dynamics | [S17] |
| Shape/topological derivatives | Boundary and domain variations, shape optimization, free-boundary problems | Shape itself is the decision variable | [S17] |
| Infinite-dimensional analysis | Banach/Hilbert differentiation, Fréchet vs Gâteaux derivatives, weak/strong convergence | Directional derivative need not be full Fréchet derivative | [S14], [S30] |
| Functional/measure derivatives | Functional gradients, measure-valued gradients, Wasserstein calculus, Lions derivatives | Specify topology/measure space before differentiating | [S14], [S30] |

## Alternative foundations and formal calculi

These are important neighbors, but **not synonyms** for the ordinary integral/derivative. Keep them searchable under separate namespaces or lessons.

| Family | Subtypes | Boundary / caveat | Entry references |
| --- | --- | --- | --- |
| Constructive/intuitionistic analysis | Bishop-style constructive analysis, Brouwerian analysis, computable analysis, constructive integration | Constructive proof obligations differ from classical existence proofs | [S22] |
| Nonstandard/synthetic analysis | Hyperreal infinitesimal nonstandard analysis; smooth infinitesimal analysis using nilpotents | Different logical foundations; do not mix their infinitesimal axioms | [S22] |
| Lambda calculus (CS/logic) | Untyped/typed λ-terms, β/η reduction, Church encodings, evaluation strategies, Curry–Howard links | This is a formal model of computation, **not a derivative operator** | [S25] |
| Logical calculi | Propositional and predicate calculi, sequent/natural deduction, modal, temporal and description logics | Formal inference systems, not calculus of rates of change | [S31] |
| Process calculi | π-calculus, CCS, ambient/process algebra, communicating transitions | Concurrency semantics | [S32] |
| Relational calculus | Tuple/domain relational calculi, predicate queries; distinction from relational algebra | Database query logic | [S33] |
| Combinatory calculus | SKI combinators, term rewriting, categorical combinators | Computation without bound variables | [S25] |
| Differential λ-calculus | Formal derivatives of terms and resource-sensitive computational structure | Related vocabulary; advanced formal system distinct from real analysis | [S34] |

## Algebra atlas

Algebra is a prerequisite for most calculus but also a broad family in its own right. There are multiple classification axes: structures (groups/rings/fields), linearity, commutativity, grading, geometry, computation, and logic. Core teaching sources are [S04], [S16], [S20], and [S26]; specialized Lie and computational algebra references are [S27] and [S28]. Formal proof development can be cross-checked against [S29] without claiming that prose or tests themselves are machine-checked proofs.

| Family | Branches and objects worth independent chapters | Sample invariants / operations |
| --- | --- | --- |
| Elementary | Arithmetic laws, rational expressions, equations and inequalities, systems, factoring, functions, polynomials, radicals, exponents | Domain restrictions, extraneous roots, exact fractions |
| Polynomial/computational | Euclidean division, gcd, square-free factoring, resultants, Gröbner bases, elimination, symbolic summation | `a=bq+r`, `deg r < deg b` |
| Linear | Vectors, matrices, rank/nullity, bases, linear maps, determinants, eigenvalues, SVD, least squares, sparse operators | Rank-nullity theorem, invertibility, numerical conditioning |
| Abstract | Semigroups, monoids, groups, subgroups, homomorphisms, quotient groups, rings, ideals, fields, modules, algebras | Associativity, identity, inverses, structure-preserving maps |
| Commutative | Noetherian rings, localization, ideals, modules, primary decomposition, commutative algebraic geometry | `ab=ba`; localization and ideal structure |
| Noncommutative | Free associative algebras, operator algebras, matrix rings, Weyl algebras, quantum groups | Multiplication order matters |
| Lie | Lie algebras/groups, brackets, representations, enveloping algebras, root systems | Bilinear antisymmetry and Jacobi identity |
| Exterior/graded | Exterior (Grassmann) algebra, wedge product, symmetric algebras, differential graded algebra, superalgebra | `v∧v=0`, graded sign rules |
| Tensor/multilinear | Tensor products, duality, contractions, multilinear forms, tensor decompositions | Basis covariance and tensor rank |
| Geometric/Clifford | Clifford algebra, geometric algebra, quaternions and geometric products | `v²=Q(v)`, geometric versus wedge product |
| Boolean/logical | Boolean algebra, propositional expressions, lattices, Heyting algebras, finite-state logic | Idempotent conjunction/disjunction and complement laws |
| Lattice/order | Partial orders, semilattices, complete lattices, residuated lattices, ordered groups/fields | Meet/join, order-preservation |
| Tropical/idempotent | Min-plus/max-plus semirings, tropical polynomials and geometry | Tropical addition is min/max, not usual `+` |
| Homological | Chain complexes, kernels/cokernels, exact sequences, derived functors, homology/cohomology | `d∘d=0` |
| Homotopical/higher | Operads, simplicial/cubical objects, ∞-categories, higher algebra | Requires explicit model/conventions |
| Representation | Group/ring/Lie representations, character theory, module categories | Structure-preserving maps into linear operators |
| Universal/categorical | Universal algebra, free objects, adjunctions, monads, category theory | Universal properties and morphism equations |
| Arithmetic/number-theoretic | Modular arithmetic, finite fields, algebraic numbers, cyclotomic rings, valuations, p-adic numbers | Exact modular congruences, field extensions |
| Algebraic geometry | Coordinate rings, varieties, schemes, sheaves, divisors, intersection theory | Polynomial equations linked to geometric spaces |
| Operator/functional | Normed/Banach/C*-algebras, spectral theory, functional calculus | Topological/norm compatibility |
| Nonassociative | Jordan, alternative, octonion, Poisson, vertex algebras | Bracketing and identities must be explicit |
| Combinatorial | Incidence/Hopf algebras, symmetric functions, matroids, algebraic combinatorics | Recurrences, coproducts, combinatorial identities |
| Relational | Relational algebra, joins, selections, projections, semiring provenance | Database relation operations |
| Differential algebra | Differential rings/fields, derivations, differential Galois theory | Leibniz rule applied to abstract rings |
| Difference algebra | Rings with endomorphisms, recurrence algebra, Ore operators | Shift operators and twisted products |
| Geometric and probabilistic interfaces | Symplectic/Poisson algebras, stochastic matrices, moment algebras | Geometry and probability constraints |

## Independent witness examples for future executable lessons

An example is a **test target**, not evidence that a feature currently exists.

| Branch | Input | Expected independent result | Pitfall exposed |
| --- | --- | --- | --- |
| Elementary algebra | `(x²-1)/(x-1)` | `x+1` only for `x ≠ 1` | Cancelling a domain restriction |
| Linear algebra | `A=[[1,2],[2,4]]` | `det(A)=0`, rank 1 | Singular systems |
| Single-variable derivative | `f(x)=x³`, `x=2` | `f'(2)=12` | Step-size or symbol handling |
| Integral | `∫_0^1 x² dx` | `1/3` | Quadrature accuracy versus exact result |
| Multivariable | `f=x²y`, at `(2,3)` | `∇f=(12,4)` | Coordinate order |
| Vector calculus | `F=(x,y,z)` | `div F = 3` | Partial derivative handling |
| Finite difference | `Δ_1(x²)` | `2x+1` | Forgetting discrete offset |
| Fractional (Caputo) | `D_C^α 1`, `0<α<1` | `0` | Confusing Caputo and Riemann–Liouville |
| Stochastic (Itô) | `f(W_t)=W_t²` | `d(W_t²)=2W_t dW_t + dt` | Missing quadratic variation |
| Variations | `L=(y')²/2 - V(y)` | `y'' + V'(y) = 0` | Incorrect Euler–Lagrange signs |
| Exterior algebra | `a` one-form | `a∧a=0` | Missing antisymmetry |
| λ-calculus | `(λx.x) y` | β-reduces to `y` | Treating substitution as numeric derivative |

## Contribution quality bar

For every new paradigm, provide: precise definition and assumptions; relation to neighboring paradigms; at least one worked original example; executable API if applicable; domain/error behavior; a counterexample; test or proof witness; computational complexity/precision note if relevant; original source; license status. Separate **mathematical correctness** from **software implementation status**. Avoid implying a numeric approximation equals a theorem or that two similarly named generalized derivatives agree.

**Branch topology proposal:** retain an integrated `main` index/library and add curated, explicitly named `topic/<paradigm>` branches for deep treatments. The repository description's "each paradigm ... in its own branch" is therefore treated as an expansion target, not falsely marked as already accomplished. Avoid duplicating a full independent codebase per branch; pin the shared-library commit and record topic-specific examples/proofs in the branch. This is a proposal for the integrator, not an authorization to create branches.
