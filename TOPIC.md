# Dedicated branch: exterior

Inherited shared foundation: `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`.

This branch curates one independently sourced paradigm. The document below is an original teaching guide, not proof of a general mathematical solver. The shared tested library remains on `main`.

---

# Exterior calculus and differential forms

**Definition.** The exterior product \(\wedge\) is alternating: for one-forms \(a,b\), \(a\wedge b=-b\wedge a\), hence \(a\wedge a=0\) over characteristic not two. The exterior derivative sends a \(k\)-form to a \((k+1)\)-form and obeys \(d^2=0\) and the graded Leibniz rule.

**Worked example.** In coordinates \((x,y)\), let \(\omega=x\,dy\). Then \(d\omega=dx\wedge dy\). For \(\eta=y\,dx\), \(d\eta=dy\wedge dx=-dx\wedge dy\). Thus \(d(\omega+\eta)=0\). This is an exact sign check independent of numerical quadrature.

**Failure boundary.** The wedge product is not the ordinary scalar multiplication of coefficients: basis indices and antisymmetry determine cancellation. Integrating forms requires orientation and manifold hypotheses; an implementation of wedge multiplication alone is not a Stokes theorem verifier.

**Repository relation.** A narrow homogeneous-form wedge product exists in `calc.advanced`; differential-form exterior differentiation and integration require separate verified APIs. Continue to [vector calculus](vector.md) and [differential geometry](differential-geometry.md).

**Primary starting source:** [MIT Mathematics: Theory of Differential Forms](https://math.mit.edu/~shaoyunb/18.952/18_952_Syllabus.pdf); [SageMath Clifford algebra manual](https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/clifford_algebra.html) for related computation.

