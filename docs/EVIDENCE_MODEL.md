# Topic evidence register — source contract v1

The canonical machine-readable register lives at [`src/calc/topic_manifest.json`](../src/calc/topic_manifest.json), alongside the Python package so it can be distributed in the wheel. It covers the **78 slugs present in `calc.catalog` at commit `c4df9e25ce23ef553b988f4b156c0e097b5065f1`**. It is a provenance-oriented map of claims, not a proof that mathematical disciplines have been exhausted or solved.

## What a record establishes

Each `topics[]` record contains:

| Field | Meaning and evidence boundary |
| --- | --- |
| `slug`, `family`, `kind` | Stable lookup identifier, current catalogue family, and one coarse category: `analysis`, `algebra`, `formal`, or `computation`. A cross-disciplinary theory may require more than one classification later. |
| `scope` | A deliberately limited statement of the mathematical subject or the concrete demonstrated operation. The phrase does **not** claim exhaustive capability. |
| `status` | Mirror of `calc.catalog` at the recorded baseline: `implemented` or `documented`. The former means selected executable examples, not whole-field implementation. |
| `evidence_level` | `indexed_only` for topics with no qualified local API; `example_tested` when a concrete API and associated test file are identified. This is not a formal theorem or real-system qualification. |
| `guide` | Path to an original `docs/topics/<slug>.md` guide if one exists; otherwise null. All 19 existing guides are mapped exactly once. |
| `branch` | Conventional `topic/<slug>` Git ref when a guide exists. It is an index entry, **not** a live assertion that the remote ref is current, tested, protected, or even reachable. |
| `source_ids` | Identifiers from `docs/SOURCES.md`. They refer to linked educational/software/academic sources. They do not assert direct proof coverage or permission to reuse source content. |
| `bindings` | Importable Python module plus one or more public classes/functions demonstrating a selected operation. Empty for topics without an example. |
| `test_paths` | Existing local unit-test file(s) associated with the binding. Their presence and module imports are checked; individual mathematical assertions still require direct review and execution. |
| `formal_proof` | False for all current records. A numerical/structural unit test is not a machine-checked proof. |

The manifest also stores the exact catalogue Git base and the intended interpretation of status, sources, and branches. That base SHA refers to **source provenance**, not current GitHub branch status or runtime qualification.

## Independent consistency checks

Run from the repository root:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_manifest_data.py -v
```

The five contract tests verify:

1. Exactly one record for each current catalogue slug, with no duplicate or unknown entries.
2. Catalogue status/family identity, explicit scope, defined evidence vocabulary and false formal-proof claims.
3. Exact coverage of all `docs/topics/*.md` guides, canonical paths and branch names.
4. Every source key is defined in the 34-entry `docs/SOURCES.md` table; no invented references.
5. Every claimed implemented binding has importable public symbols and existing test-file paths.

These checks are structural/source evidence. They cannot establish branch CI, numerical convergence across a function class, source-license clearance, or correctness of every algorithm. The regular mathematical test suite remains independent evidence.

## Promotion protocol for new operations

When a new method is implemented (for example the pending time-scale, interval, or computational-algebra modules), the integrator should update **both** `calc.catalog` and the manifest in one reviewed change:

1. Read the exact Git head and test the new public operation against independently calculated mathematical identities and meaningful counterexamples.
2. Record the narrow `scope`, current `source_ids`, the importable API symbols, and actual test paths. Never mark a complete theory implemented because one operator works.
3. Promote `status` to `implemented` and `evidence_level` to `example_tested` only after focused tests pass; leave `formal_proof` false unless an independently checkable proof artifact exists.
4. Add or revise a subject guide, with domain assumptions, edge cases, uncertainty, attribution, and a mathematical witness.
5. Reconcile a corresponding `topic/<slug>` branch with an exact tested base, preserving its guide and merge ancestry. Verify the branch ref and differential against `main` through GitHub. A branch string in the manifest does not substitute for this readback.
6. Run all manifest checks and the integrated suite, verify built-package inclusion, inspect the diff, publish with ref lease and read back the exact remote effect.

Conversely, if hostile review invalidates an implementation's witness, demote the record or mark the operation qualification unresolved rather than silently preserving a positive status.

## Known limitations of the current bibliography

Several `source_ids` are broad surveys or API inventories, not direct proof-level evidence for specialized claims. Examples include Malliavin calculus (`S15`, `S24`), root finding (`S14`), category theory (`S13`), and tropical algebra (`S14`). Dedicated primary references should replace these proxies as those chapters deepen. The review ledger in [`RESEARCH_GAPS.md`](RESEARCH_GAPS.md) ranks stronger source and theorem-level work.

The user-suggested [ML-foundations](https://github.com/jonkrohn/ML-foundations) is linked as source `S04`; its MIT license was verified for reference purposes. No notebook code, text, figures, or protected educational assets were imported into this manifest. External sources have independent reuse terms. The repository itself still has no selected project license.

## What remains unresolved

- **Exhaustiveness:** There is no finite agreed list of every calculus or algebra. The 78 subjects are an expandable classification, not a completeness certificate.
- **Branch identity:** These data do not themselves resolve current remote heads, CI runs, or ancestry after subsequent changes.
- **Semantic mapping:** Importable symbols and associated tests are a useful audit spine, but independent API-level mathematical claims need explicit invariant-to-test traceability and, where appropriate, formal proofs.
- **Packaging:** The parent integrator owns package-data configuration and wheel smoke tests; this specialist owns the register and structural tests only.

This contract keeps an ambitious research programme reconstructible without confusing reference links, code existence, unit tests, remote effects and mathematical proof.
