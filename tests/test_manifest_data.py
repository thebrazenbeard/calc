"""Source-first consistency checks for the mathematical evidence register."""

from __future__ import annotations

import importlib
import json
from pathlib import Path
import re
import unittest

from calc.catalog import TOPICS


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "src" / "calc" / "topic_manifest.json"


class EvidenceRegistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.registry = json.loads(MANIFEST.read_text(encoding="utf-8"))
        cls.records = cls.registry["topics"]

    def test_schema_and_exact_catalog_coverage(self) -> None:
        self.assertEqual(self.registry["schema_version"], 1)
        self.assertRegex(self.registry["catalog_base_commit"], r"^[a-f0-9]{40}$")
        slugs = [record["slug"] for record in self.records]
        self.assertEqual(len(slugs), len(set(slugs)))
        self.assertEqual(set(slugs), {topic.slug for topic in TOPICS})

    def test_status_and_evidence_boundaries(self) -> None:
        catalogue = {topic.slug: topic for topic in TOPICS}
        for record in self.records:
            with self.subTest(slug=record["slug"]):
                topic = catalogue[record["slug"]]
                self.assertEqual(record["family"], topic.family)
                self.assertEqual(record["status"], topic.status)
                self.assertIn(record["kind"], {"analysis", "algebra", "formal", "computation"})
                self.assertGreater(len(record["scope"].strip()), 12)
                self.assertEqual(record["formal_proof"], False)
                if topic.status == "implemented":
                    self.assertEqual(record["evidence_level"], "example_tested")
                    self.assertTrue(record["bindings"])
                else:
                    self.assertEqual(record["evidence_level"], "indexed_only")
                    self.assertEqual(record["bindings"], [])
                    self.assertEqual(record["test_paths"], [])

    def test_guides_and_branches_are_exactly_mapped(self) -> None:
        guides = {str(p.relative_to(ROOT)) for p in (ROOT / "docs" / "topics").glob("*.md")}
        recorded = {record["guide"] for record in self.records if record["guide"] is not None}
        self.assertEqual(guides, recorded)
        for record in self.records:
            with self.subTest(slug=record["slug"]):
                if record["guide"] is None:
                    self.assertIsNone(record["branch"])
                else:
                    self.assertEqual(record["guide"], f"docs/topics/{record['slug']}.md")
                    self.assertEqual(record["branch"], f"topic/{record['slug']}")

    def test_source_identifiers_exist(self) -> None:
        source_text = (ROOT / "docs" / "SOURCES.md").read_text(encoding="utf-8")
        allowed = set(re.findall(r"^\| (S\d{2}) \|", source_text, re.MULTILINE))
        self.assertEqual(len(allowed), 34)
        for record in self.records:
            with self.subTest(slug=record["slug"]):
                self.assertEqual(len(record["source_ids"]), len(set(record["source_ids"])))
                self.assertTrue(record["source_ids"])
                self.assertLessEqual(len(record["source_ids"]), 4)
                self.assertTrue(set(record["source_ids"]) <= allowed)

    def test_executable_bindings_and_tests_exist(self) -> None:
        for record in self.records:
            if record["status"] != "implemented":
                continue
            with self.subTest(slug=record["slug"]):
                self.assertTrue(record["test_paths"])
                for test_path in record["test_paths"]:
                    self.assertRegex(test_path, r"^tests/test_[a-z0-9_]+\.py$")
                    self.assertTrue((ROOT / test_path).is_file())
                for binding in record["bindings"]:
                    self.assertRegex(binding["module"], r"^calc\.[a-z_]+$")
                    module = importlib.import_module(binding["module"])
                    self.assertTrue(binding["symbols"])
                    for symbol in binding["symbols"]:
                        self.assertFalse(symbol.startswith("_"))
                        self.assertTrue(callable(getattr(module, symbol, None)))


if __name__ == "__main__":
    unittest.main()
