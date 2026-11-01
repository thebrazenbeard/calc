# Finite measure theory and conditional expectation

`calc.measure` is a dependency-free, exact-rational laboratory for finite
measure spaces. It implements selected fundamental identities rather than
claiming to solve measure theory on general infinite sets.

## Mathematical model

Let \(\Omega=\{x_1,\ldots,x_n\}\), where the weights
\(\mu_i=\mu(\{x_i\})\in\mathbb{Q}_{\ge 0}\).
Every subset of \(\Omega\) is measurable, and

\[
\mu(A)=\sum_{x_i\in A}\mu_i,
\qquad \int f\,d\mu=\sum_i f(x_i)\mu_i.
\]

All supplied masses and observable values must be exact `int` or `Fraction`;
floats and booleans are deliberately rejected. An atom may have zero mass.
There can be at most 128 labelled atoms and coefficient inputs have a
2048-bit limit. Results are `Fraction` values with a separate size bound.

If \(M=\mu(\Omega)>0\), probability normalization is
\(P(\{x_i\})=\mu_i/M\). The API also provides
\(E[f]=M^{-1}\int f\,d\mu\) and
\(\operatorname{Var}(f)=E[(f-E[f])^2]\).

### Finite Radon–Nikodym density

For two finite measures \(\nu\) and \(\mu\) on **the same atoms**,
\(\nu\ll\mu\) precisely when \(\mu_i=0\Rightarrow\nu_i=0\). Then

\[
\frac{d\nu}{d\mu}(x_i)=
\begin{cases}
\nu_i/\mu_i,&\mu_i>0,\\
0,&\mu_i=0\quad\text{(chosen convention).}
\end{cases}
\]

The second branch is only a convention: densities are unique *almost
everywhere*, not pointwise on null atoms. The implementation rejects
non-absolutely-continuous pairs rather than dividing by zero.

### Conditional expectation onto a finite partition

For a disjoint partition \(\mathcal P=\{C_1,\ldots,C_k\}\) covering
the atom set, conditional expectation is constant on every cell:

\[
E_\mu[f\mid\mathcal P](x)
=\frac{\sum_{y\in C_j}f(y)\mu_y}{\mu(C_j)}
\quad (x\in C_j,\;\mu(C_j)>0).
\]

It is set to zero on a null cell (another almost-sure convention).
The defining identity
\(\int_{C_j}E[f\mid\mathcal P]d\mu=\int_{C_j}f\,d\mu\)
holds exactly. The tests also verify the tower identity for nested finite
partitions and nonnegative variance identities.

## Executable example

```python
from fractions import Fraction as F
from calc.measure import FiniteMeasure

reference = FiniteMeasure({'a': 1, 'b': 3, 'null': 0})
assert reference.total_mass == 4
assert reference.normalize().mass(('b',)) == F(3, 4)

f = {'a': 2, 'b': 6, 'null': 1000}
assert reference.expectation(f) == 5
assert reference.variance(f) == 3

conditional = reference.conditional_expectation(
    f, (('a', 'b'), ('null',)),
)
assert conditional == {'a': F(5), 'b': F(5), 'null': F(0)}
assert reference.integrate(conditional) == reference.integrate(f)

target = FiniteMeasure({'a': F(1, 2), 'b': F(3, 2), 'null': 0})
density = reference.density_of(target)
assert density == {'a': F(1, 2), 'b': F(1, 2), 'null': F(0)}
assert reference.integrate(density) == target.total_mass
```

Run the independently derived rational tests:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_measure.py -v
```

## Boundary of the claim

- Finite labelled atoms use the full power-set sigma-algebra implicitly;
  no arbitrary or infinite measurable space is accepted.
- This implementation does **not** construct Lebesgue measure, perform
  abstract sigma-algebra generation, or establish the general
  Radon–Nikodym theorem.
- All input observables must provide a value on every atom, even a null atom.
  Different values on null atoms yield the same integrals.
- `conditional_expectation` accepts a genuine finite partition, not an
  arbitrary sub-sigma-algebra representation. It chooses zero on null cells.
- `normalize` requires positive total mass. Exact rational outputs still
  have explicit resource ceilings and are not symbolic proof objects.

## Mathematical references and provenance

- [MIT 18.175 Theory of Probability, Lecture 25](https://ocw.mit.edu/courses/18-175-theory-of-probability-spring-2014/f33e55917467aa3cf5afa57ec32f4a6a_MIT18_175S14_Lecture25.pdf) discusses the general Radon–Nikodym theorem; the code implements its elementary finite-atomic special case only.
- [MIT Theory of Probability lecture notes, Section 2](https://dspace.mit.edu/bitstream/handle/1721.1/96865/18-175-fall-2008/contents/lecture-notes/section_2.pdf) discusses absolute continuity and the Radon–Nikodym formula.
- [Mathias Drton, Advanced Probability lecture notes](https://www.mat.univie.ac.at/~mathias/Advanced_Probability.pdf) includes conditional expectation for finite partitions.

These are linked mathematical references. No external source code or
third-party repository text was copied into this implementation.
