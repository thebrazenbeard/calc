# Dedicated branch: q-calculus

Inherited shared foundation: `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`.

This branch curates one independently sourced paradigm. The document below is an original teaching guide, not proof of a general mathematical solver. The shared tested library remains on `main`.

---

# q-calculus

**Definition.** The Jackson q-derivative of a function at nonzero \(x\), with \(q\ne1\), is

\[D_q f(x)=\frac{f(qx)-f(x)}{(q-1)x}.\]

It uses a multiplicative rather than additive shift and approaches the ordinary derivative for sufficiently smooth \(f\) as \(q\to1\).

**Worked example.** For \(f(x)=x^n\), \(D_qx^n=[n]_q x^{n-1}\), where \([n]_q=(q^n-1)/(q-1)\). For \(n=2\), \(D_q x^2=(q+1)x\). With \(q=2,x=3\), the exact answer is \(9\), whereas the ordinary derivative gives \(6\).

**Failure boundary.** The displayed formula has a zero denominator at \(x=0\) or \(q=1\). Extensions by limits may exist for specific functions, but code must state whether it accepts those points. Variants of q-integration and \((p,q)\)-operators are not interchangeable.

**Repository relation.** `calc.advanced` has an illustrative Jackson derivative; it is not a generalized q-special-function library. Compare [discrete calculus](discrete.md), [time-scale calculus](time-scales.md) and [differential calculus](differential.md).

**Primary starting source:** [Research article on q-calculus, Hindawi (2018)](https://onlinelibrary.wiley.com/doi/full/10.1155/2018/9089865).

