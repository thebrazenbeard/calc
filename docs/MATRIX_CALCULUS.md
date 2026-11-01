# Matrix calculus: explicit contracts and examples

**Implementation:** [`src/calc/matrix_calc.py`](../src/calc/matrix_calc.py). These are selected real-matrix operations, not an arbitrary tensor-differentiation system or a machine-checked proof library.

## Frobenius pairing and derivatives

For real matrices of the same shape, the Frobenius inner product is

\[
\langle A,B\rangle_F = \operatorname{tr}(A^T B)=\sum_{i,j} A_{ij}B_{ij}.
\]

The gradient is defined relative to this pairing: the first differential of a scalar map \(F\) is \(dF(X)[H]=\langle\nabla_X F,H\rangle_F\).

Two analytic identities are implemented:

1. \(\nabla_X \|X\|_F^2 = 2X\).
2. For \(A\in\mathbb R^{p\times m}\), \(X\in\mathbb R^{m\times n}\), \(B\in\mathbb R^{n\times p}\), \(F(X)=\operatorname{tr}(AXB)\) has gradient \(\nabla_X F=(BA)^T\in\mathbb R^{m\times n}\).

For example, \(A=\begin{pmatrix}1&2\\3&4\end{pmatrix}\) and \(B=\begin{pmatrix}5&6\\7&8\end{pmatrix}\) give \(\nabla_X\operatorname{tr}(AXB)=\begin{pmatrix}23&31\\34&46\end{pmatrix}\). This is independently checked against the derivative of a scalar trace formula in the tests.

`frobenius_inner`, `squared_frobenius_gradient` and `trace_axb_gradient` return exact `fractions.Fraction` entries if *all* inputs are integers or rational fractions. The presence of any float selects finite binary floating-point calculations, whose output is not exact mathematical arithmetic.

## Numerical matrix gradients and Jacobian-vector products

`matrix_gradient(f, X)` approximates each gradient entry by a centered coordinate difference:

\[
(\nabla_X f)_{ij}\approx \frac{f(X+h_{ij}E_{ij})-f(X-h_{ij}E_{ij})}{2h_{ij}},
\qquad h_{ij}=\varepsilon\max(1,|X_{ij}|).
\]

`matrix_jvp(F, X, H)` estimates the directional derivative of a real matrix-to-matrix function:

\[
DF(X)[H]\approx \frac{F(X+hH)-F(X-hH)}{2h}.
\]

The transform receives immutable tuples of finite floats. Its output must stay rectangular with the same shape under both perturbations. A nonzero direction that cannot perturb an input coordinate at available float resolution causes an error; silently returning zero would be misleading.

```python
from calc.matrix_calc import frobenius_inner, matrix_gradient, trace_axb_gradient

assert frobenius_inner(((1, 2),), ((3, 4),)) == 11
assert trace_axb_gradient(((1, 2), (3, 4)), ((5, 6), (7, 8))) == ((23, 31), (34, 46))
approx = matrix_gradient(lambda x: sum(v*v for row in x for v in row), ((1., 2.),))
assert abs(approx[0][0] - 2) < 1e-5
```

### Boundaries and failure behavior

- Real matrices only; complex Wirtinger calculus, generalized spectral derivatives, sparse operators and reverse-mode differentiation are outside scope.
- Dimensions must be nonempty, rectangular and compatible. Boolean and nonfinite inputs are rejected.
- Floating-point underflow, catastrophic cancellation and ill-conditioning may make a finite-difference estimate inaccurate even when the inputs are accepted. A smooth function and adequate step scaling are assumptions, not proven facts.
- Analytic identities are exact over rational entries. Numerical functions convert input to finite floating-point values and do not certify truncation or roundoff error.
- User-provided callables execute normally; do not pass untrusted functions.

### References

- [MIT OpenCourseWare, 18.S096 Matrix Calculus for Machine Learning and Beyond](https://ocw.mit.edu/courses/18-s096-matrix-calculus-for-machine-learning-and-beyond-january-iap-2023/) — notation and gradient/differential context.
- [MIT OpenCourseWare, 18.06SC Linear Algebra](https://ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/) — trace, matrix products and vector-space prerequisites.
- [SciPy numerical differentiation](https://docs.scipy.org/doc/scipy/reference/differentiate.html) — a reference for numerical-differentiation limitations; this module does not incorporate SciPy code.

This text and implementation are original; source links are attribution and study references, not imported code.
