# A bounded, exact finite-element example: 1D Poisson

The `calc.fem` module implements a concrete finite-element problem, not an
arbitrary-dimensional PDE solver:

\[
-u''(x)=f(x),\qquad a<x<b,\qquad u(a)=u(b)=0.
\]

The nodes `a=x₀<x₁<⋯<xₙ=b` are finite *strictly increasing rational*
coordinates, and `f` is constant on each mesh element. Trial and test functions
are continuous piecewise-affine (P1/hat) functions that vanish at both endpoints.

## Weak form and assembly

Multiplying the differential equation by an interior hat function `φᵢ`,
integrating by parts, and using zero boundary values gives

\[
\int_a^b u_h'(x)\,\varphi_i'(x)\,dx
=\int_a^b f(x)\,\varphi_i(x)\,dx.
\]

On element `e=[xᵢ,xᵢ₊₁]` with length `hₑ` and constant load `fₑ`, the exact local
matrices are

\[
K_e=\frac1{h_e}\begin{bmatrix}1&-1\\-1&1\end{bmatrix},
\qquad
b_e=\frac{f_e h_e}{2}\begin{bmatrix}1\\1\end{bmatrix}.
\]

Assembling only *interior* degrees of freedom produces a symmetric positive
definite tridiagonal system. An exact-rational Thomas elimination solves it in
linear numbers of rational operations. This does not mean bit complexity is
linear; large rational denominators can grow, so the implementation has explicit
size/precision limits.

## Executable example

```python
from fractions import Fraction as Q
from calc.fem import assemble_poisson_1d, solve_poisson_1d

nodes = (Q(0), Q(1, 4), Q(1, 2), Q(1))
system = assemble_poisson_1d(nodes, forcing=1)
assert system.diagonal == (Q(8), Q(6))
assert system.upper == (Q(-4),)

sol = system.solve()
assert sol.nodal_values == (Q(0), Q(3, 32), Q(1, 8), Q(0))
assert system.residual(sol.nodal_values[1:-1]) == (Q(0), Q(0))

midpoint = Q(1, 8)
assert sol.evaluate(midpoint) == Q(3, 64)

# Equivalent shorthand:
assert solve_poisson_1d(nodes).nodal_values == sol.nodal_values
```

The exact continuous solution for `f=1` on `[0,1]` is
`u(x)=x(1-x)/2`. For this particular load, P1 FEM has *exact nodal values*;
between nodes, however, the piecewise-affine interpolant differs from the
quadratic by `h²/8` at the midpoint of an element of length `h`. The tests verify
the predicted factor-four error decrease when the mesh width halves. Exact
algebraic residual zero is **not** zero PDE discretization error.

The call `assemble_poisson_1d(nodes, forcing=(f₀, f₁, …))` accepts a separate
constant rational forcing value for each element. The returned object's
`residual(interior_values)` evaluates `K·v−b` exactly. The solution object's
`energy_seminorm_squared()` evaluates `∫|u_h'|²` exactly on the mesh.

## Deliberate boundaries

- Accepts `int` or `Fraction`, but **not binary floats**: exactness refers to
  supplied rational data, not a high-precision reconstruction of measured
  values. Use `Fraction("0.1")` to represent the rational `1/10` explicitly.
- Homogeneous Dirichlet boundary conditions only; no Neumann, Robin, periodic,
  or time-dependent boundary operators.
- One-dimensional mesh; linear P1 basis; piecewise-constant body force;
  `−u''=f` with unit diffusion coefficient. No unstructured 2D/3D meshes,
  nonlinear PDEs, varying coefficients, advection, or adaptive refinement.
- At most 513 mesh nodes and 4096-bit numerators/denominators in checked
  quantities. Excessive input and intermediate sizes raise `ValueError`.
- Exact arithmetic proves identities of the *discrete* linear system; it is
  neither a formal proof of continuous PDE well-posedness nor a certified error
  estimate for arbitrary solutions.

## Sources

- [The FEniCS Book, *Automated Solution of Differential Equations by the Finite
  Element Method*](https://fenicsproject.org/pub/book/book/fenics-book-2011-06-14.pdf)
  — finite-element weak forms, assembly, and Poisson examples.
- [FEniCS/UFL documentation](https://fenics.readthedocs.io/_/downloads/ufl/en/latest/pdf/)
  — representation of variational forms; the implementation here does *not*
  depend on FEniCS.

Formulas and code here were independently written for this repository.
