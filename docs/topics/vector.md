# Vector calculus

**Definition.** For a differentiable vector field \(F=(F_1,F_2,F_3)\), the divergence \(\nabla\cdot F=\sum_i\partial_iF_i\) measures infinitesimal net flux. Curl \(\nabla\times F\) measures infinitesimal circulation in oriented three-dimensional Euclidean coordinates. The geometric integral theorems connect these local operators with boundaries.

**Worked example.** For \(F=(x,y,z)\), \(\nabla\cdot F=1+1+1=3\) and \(\nabla\times F=0\). Over the unit ball, the divergence theorem predicts outward flux \(3\times(4\pi/3)=4\pi\), matching the sphere's area because \(F\cdot n=1\) there.

**Failure boundary.** The textbook three-component cross-product curl is dimension-specific; general manifolds need differential forms and a metric/Hodge operator. Stokes/divergence theorems require appropriate smoothness, orientation and boundary hypotheses.

**Repository relation.** `calc.differential.divergence` and `.curl` provide numerical component-field derivatives. The exact geometric theorems are educational references; see [exterior calculus](exterior.md) and [differential geometry](differential-geometry.md).

**Primary starting source:** [MIT OCW 18.02SC](https://ocw.mit.edu/courses/18-02sc-multivariable-calculus-fall-2010/).

