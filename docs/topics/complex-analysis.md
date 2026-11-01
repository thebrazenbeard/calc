# Complex-variable calculus

**Definition.** A function \(f:\mathbb C\to\mathbb C\) is complex differentiable at \(z_0\) when \((f(z_0+h)-f(z_0))/h\) converges as **complex** \(h\to0\), independent of direction. Holomorphicity on a neighbourhood is much stronger than real differentiability in two coordinates. Under suitable regularity, the Cauchy–Riemann equations express the constraint.

**Worked example.** For \(f(z)=z^2\),

\[\frac{(z+h)^2-z^2}{h}=2z+h\to2z.\]

The limit is direction-independent, so \(f\) is entire. In contrast, \(g(z)=\bar z\) gives \(\bar h/h\), which approaches different values along the real and imaginary axes; it is nowhere complex differentiable.

**Failure boundary.** A real Jacobian existing at a point does not imply a complex derivative. Contour-integral theorems require specific domain, analyticity and orientation assumptions. The popular identity \(oint dz/z=2\pi i\) is valid for a positively oriented loop enclosing zero once, not for any arbitrary curve.

**Repository relation.** Currently a reference/guide, not a verified complex-symbolic engine. Compare [multivariable calculus](multivariable.md) for real-coordinate derivatives and [integral calculus](integral.md) for the distinction from complex contour integrals.

**Primary starting source:** [MIT OCW 18.04 Complex Variables with Applications](https://ocw.mit.edu/courses/18-04-complex-variables-with-applications-spring-2018/).

