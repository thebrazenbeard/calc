# Dedicated branch: calculus-of-variations

Inherited shared foundation: `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`.

This branch curates one independently sourced paradigm. The document below is an original teaching guide, not proof of a general mathematical solver. The shared tested library remains on `main`.

---

# Calculus of variations

**Definition.** Variational calculus optimizes functionals—maps from admissible functions to numbers—rather than merely optimizing a finite-dimensional function. For \(J[y]=\int_a^b L(x,y,y')\,dx\) with fixed endpoints and suitable smoothness, a stationary path obeys the Euler–Lagrange equation

\[\frac{\partial L}{\partial y}-\frac{d}{dx}\frac{\partial L}{\partial y'}=0.\]

**Worked example.** Let \(L=(y')^2/2\). Then \(\partial L/\partial y=0\) and \(\partial L/\partial y'=y'\), hence \(y''=0\). With \(y(0)=0\), \(y(1)=2\), the stationary solution is \(y(x)=2x\). Strict convexity here also establishes it is the minimizer.

**Failure boundary.** Stationarity alone is not necessarily a minimum. Free endpoint problems require natural boundary conditions, constraints require additional multipliers or suitable methods, and low-regularity extremizers can demand weak formulations.

**Repository relation.** `calc.advanced` offers a finite-difference Euler–Lagrange residual illustration, not a functional optimizer. See [multivariable calculus](multivariable.md) and [differential geometry](differential-geometry.md).

**Primary starting source:** [MIT OCW 16.323, Calculus of Variations lecture](https://ocw.mit.edu/courses/16-323-principles-of-optimal-control-spring-2008/resources/lec5/).

