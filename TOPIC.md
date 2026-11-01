# Dedicated branch: lambda-calculus

Inherited shared foundation: `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`.

This branch curates one independently sourced paradigm. The document below is an original teaching guide, not proof of a general mathematical solver. The shared tested library remains on `main`.

---

# Lambda calculus (formal computation)

**Definition.** Untyped lambda calculus has variables, abstractions \(\lambda x.M\), and applications \(MN\). Beta reduction replaces \((\lambda x.M)N\) by \(M[x:=N]\), using capture-avoiding substitution. This is a formal computational calculus, not differentiation of real functions.

**Worked example.** The identity term \(I=\lambda x.x\) satisfies \(Iy=(\lambda x.x)y\to_\beta y\). For Church numerals, the successor operation maps a representation of \(n\) to a representation of \(n+1\); it does not compute a numerical derivative.

**Failure boundary.** Naive textual replacement is incorrect: substituting \(y\) for \(x\) inside \(\lambda y.x\) must rename the bound \(y\) before insertion to avoid variable capture. Evaluation order also matters for non-normalizing expressions; untyped lambda calculus is not strongly normalizing.

**Repository relation.** Included as a clearly separate formal branch of the word *calculus*. No beta-reduction interpreter or type checker is claimed. See [abstract algebra](abstract-algebra.md) and [symbolic/computational catalog](../CATALOG.md).

**Primary starting source:** [Stanford Encyclopedia of Philosophy: Lambda Calculus](https://plato.stanford.edu/entries/lambda-calculus/).

