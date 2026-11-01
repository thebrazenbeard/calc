# Linear algebra

**Definition.** A linear map \(T:V\to W\) preserves sums and scalar multiplication: \(T(u+v)=T(u)+T(v)\), \(T(au)=aT(u)\). Matrices represent finite-dimensional maps after choices of bases. Rank and determinant encode structural properties but have different hypotheses.

**Worked example.** For \(A=\begin{pmatrix}1&2\\2&4\end{pmatrix}\), the second row is twice the first. Therefore \(\operatorname{rank}A=1\), \(\det A=4-4=0\), and the homogeneous system has the nonzero solution \((-2,1)^T\). The matrix is not invertible.

**Failure boundary.** Solving a singular system by ordinary inversion is invalid. For floating-point matrices, algebraic non-singularity and numerical conditioning are different questions; exact rational elimination and robust large-scale numerical linear algebra require different designs.

**Repository relation.** `calc.algebra` implements several exact matrix and rational polynomial operations. Read [multivariable calculus](multivariable.md) for Jacobians and [differential geometry](differential-geometry.md) for metric tensors.

**Primary starting sources:** [MIT OCW Linear Algebra within the mathematics curriculum](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/) and user-suggested [Jon Krohn, ML-foundations (MIT)](https://github.com/jonkrohn/ML-foundations) as attributed optional learning progression.

