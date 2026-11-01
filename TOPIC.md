# Dedicated branch: multivariable

Inherited shared foundation: `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`.

This branch curates one independently sourced paradigm. The document below is an original teaching guide, not proof of a general mathematical solver. The shared tested library remains on `main`.

---

# Multivariable calculus

**Definition.** For \(f:\mathbb R^n\to\mathbb R^m\), the total derivative \(Df(a)\) is a linear map satisfying \(f(a+h)=f(a)+Df(a)h+o(\lVert h\rVert)\). The Jacobian represents this map in chosen coordinates. Individual partial derivatives can exist without guaranteeing total differentiability.

**Worked example.** With \(f(x,y)=x^2y\),

\[\nabla f=(2xy,x^2),\quad \nabla f(2,3)=(12,4).\]

For a vector map \(F=(xy,x+y)\), the Jacobian is \(\begin{pmatrix}y&x\\1&1\end{pmatrix}\), hence \(J_F(2,3)=\begin{pmatrix}3&2\\1&1\end{pmatrix}\).

**Failure boundary.** The existence of coordinate-wise partial derivatives at a point does not ensure a continuous or differentiable function. Hessians may be asymmetric at points lacking the smoothness needed for equality of mixed partials, and numerical Hessians amplify rounding error.

**Repository relation.** `calc.differential` includes approximations to partials, gradient, Jacobian and Hessian. See [vector calculus](vector.md), [differential geometry](differential-geometry.md), and [linear algebra](linear-algebra.md).

**Primary starting source:** [MIT OCW 18.02SC Multivariable Calculus](https://ocw.mit.edu/courses/18-02sc-multivariable-calculus-fall-2010/).

