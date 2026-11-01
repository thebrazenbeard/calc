# Dedicated branch: boolean-algebra

Inherited shared foundation: `97bc4f0ba576dd1cecf2e9500c6e186ca9dd23d9`.

This branch curates one independently sourced paradigm. The document below is an original teaching guide, not proof of a general mathematical solver. The shared tested library remains on `main`.

---

# Boolean algebra

**Definition.** A Boolean algebra has join \(∨\), meet \(∧\), complements and bottom/top satisfying distributive lattice laws. In the two-element algebra, join is OR, meet is AND and complement is NOT; the structure is not generally the same as arithmetic over the integers.

**Worked example.** The distributive identity \(a\wedge(b\vee c)=(a\wedge b)\vee(a\wedge c)\) holds for every Boolean valuation. For \(a=1,b=0,c=1\), each side is \(1\). The complement law \(a\wedge\neg a=0\) is also independently checked for both values.

**Failure boundary.** Logical XOR corresponds to addition modulo two in a Boolean ring, while OR does not: \(1\vee1=1\), yet \(1+1\equiv0\pmod2\). Conflating the lattice and ring operations breaks formulas.

**Repository relation.** Bridge between [abstract algebra](abstract-algebra.md) and logical calculi. Formal proof engines, satisfiability solvers and BDDs are not implied by these identities; see the [catalogue](../CATALOG.md).

**Primary starting source:** [University of Warwick: Boolean Algebra Axioms](https://www.dcs.warwick.ac.uk/~wmb/fsv/BoolAlgAxioms.html) for direct axioms; [Stanford Encyclopedia: Natural Deduction](https://plato.stanford.edu/entries/natural-deduction/) for the distinct logical-inference perspective.

