# Dedicated branch: abstract-algebra

Inherited shared foundation: `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`.

This branch curates one independently sourced paradigm. The document below is an original teaching guide, not proof of a general mathematical solver. The shared tested library remains on `main`.

---

# Abstract algebra

**Definition.** Abstract algebra studies structures through axioms and maps preserving their operations. A group has an associative binary operation, identity and inverses; a ring has compatible additive and multiplicative structures. Neither a finite group tester nor polynomial arithmetic proves facts about all groups or rings.

**Worked example.** Integers modulo five under addition form a group. The additive inverse of \(2\) is \(3\), because \(2+3\equiv0\pmod5\). Nonzero residues under multiplication also form a group since five is prime: \(2^{-1}=3\) since \(2\cdot3\equiv1\pmod5\).

**Failure boundary.** Modulo four, the nonzero residues do not all have multiplicative inverses: \(2\cdot2\equiv0\pmod4\). A structure checker must test closure, associativity, identity and inverses; a few sampled products cannot certify a large unbounded algebra.

**Repository relation.** `calc.algebra` includes modular arithmetic and exact polynomial operations; more specialized structures belong in separately checked modules. See [linear algebra](linear-algebra.md), [Lie algebra](lie-algebra.md), [Boolean algebra](boolean-algebra.md) and [tropical algebra](tropical-algebra.md).

**Primary starting source:** [MIT OCW 18.701 Algebra I](https://ocw.mit.edu/courses/18-701-algebra-i-fall-2010/).

