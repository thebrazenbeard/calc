# Dedicated branch: lie-algebra

Inherited shared foundation: `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`.

This branch curates one independently sourced paradigm. The document below is an original teaching guide, not proof of a general mathematical solver. The shared tested library remains on `main`.

---

# Lie algebra

**Definition.** A Lie algebra is a vector space with a bilinear bracket \([x,y]\) satisfying alternation \([x,x]=0\) and the Jacobi identity \([x,[y,z]]+[y,[z,x]]+[z,[x,y]]=0\). In characteristic other than two, alternation implies antisymmetry.

**Worked example.** Matrix commutator \([A,B]=AB-BA\) is a Lie bracket on square matrices. For \(A=\begin{pmatrix}1&0\\0&0\end{pmatrix}\) and \(B=\begin{pmatrix}0&1\\0&0\end{pmatrix}\), \(AB=B\) and \(BA=0\), so \([A,B]=B\). Reversing arguments gives \([B,A]=-B\).

**Failure boundary.** Not every arbitrary antisymmetric binary operation obeys Jacobi. A Lie bracket also is not ordinary associative multiplication; the commutator is derived from one associative algebra but many Lie algebras are studied abstractly.

**Repository relation.** Catalogued for future verified bracket/representation APIs. See [linear algebra](linear-algebra.md), [abstract algebra](abstract-algebra.md) and [differential geometry](differential-geometry.md).

**Primary starting source:** [SageMath official Lie Algebras reference](https://doc.sagemath.org/html/en/reference/algebras/lie_algebras.html).

