# Dedicated branch: differential-geometry

Inherited shared foundation: `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`.

This branch curates one independently sourced paradigm. The document below is an original teaching guide, not proof of a general mathematical solver. The shared tested library remains on `main`.

---

# Differential geometry and tensor calculus

**Definition.** Differential geometry studies derivatives and invariants on manifolds without treating the choice of coordinates as physical structure. A metric \(g_{ij}\) allows inner products of tangent vectors; its Levi–Civita connection has Christoffel symbols

\[\Gamma^k_{ij}=\tfrac12g^{k\ell}(\partial_i g_{j\ell}+\partial_j g_{i\ell}-\partial_\ell g_{ij}).\]

**Worked example.** In flat Euclidean plane polar coordinates \(g=\mathrm{diag}(1,r^2)\), the nonzero Christoffel coefficients are \(\Gamma^r_{\theta\theta}=-r\) and \(\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}=1/r\). Nonzero connection coefficients do not imply nonzero curvature: the underlying plane is flat.

**Failure boundary.** Christoffel symbols are not tensor components, and formulas involving \(g^{-1}\) fail at degenerate metrics or singular coordinate charts such as polar \(r=0\). A numerical differentiation routine cannot establish coordinate invariance or tensor transformation laws by testing one chart.

**Repository relation.** Separate differential-form ideas in [exterior calculus](exterior.md); compare coordinate-dependent [multivariable calculus](multivariable.md) and [linear algebra](linear-algebra.md). A general manifold engine is not claimed.

**Primary starting source:** [MIT Mathematics: Theory of Differential Forms, 18.952 syllabus](https://math.mit.edu/~shaoyunb/18.952/18_952_Syllabus.pdf), a curriculum map rather than a proof-level implementation reference.

