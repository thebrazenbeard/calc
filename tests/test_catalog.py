"""Stable discovery contracts for the subject index."""

import unittest
from contextlib import redirect_stdout
from io import StringIO

from calc.catalog import TOPICS, get_topic, list_topics
from calc.__main__ import main


class CatalogTests(unittest.TestCase):
    def test_stable_unique_slugs(self):
        self.assertEqual(len(TOPICS), len({t.slug for t in TOPICS}))
        self.assertGreaterEqual(len(TOPICS), 60)

    def test_cross_domain_coverage(self):
        for slug in ("ito", "fractional", "lambda-calculus", "clifford-algebra", "linear-algebra", "calculus-of-variations"):
            self.assertTrue(get_topic(slug).family)

    def test_filtering(self):
        algebra = list_topics("algebra-structures")
        self.assertTrue(algebra)
        self.assertTrue(all(t.family == "algebra-structures" for t in algebra))

    def test_unknown_topic(self):
        with self.assertRaises(KeyError):
            get_topic("nonexistent")

    def test_cli(self):
        output = StringIO()
        with redirect_stdout(output):
            self.assertEqual(main(["topic", "ito"]), 0)
        self.assertIn("Itô calculus", output.getvalue())

    def test_implemented_means_a_specific_example(self):
        self.assertEqual(get_topic("automatic-differentiation").status, "implemented")
        self.assertEqual(get_topic("fractional").status, "implemented")
        self.assertEqual(get_topic("malliavin").status, "documented")


if __name__ == "__main__":
    unittest.main()

