"""Public read-only topic evidence CLI contracts."""

from contextlib import redirect_stdout
from io import StringIO
import json
import unittest

from calc.__main__ import main
from calc.evidence import get_evidence, load_evidence


class EvidenceCliTests(unittest.TestCase):
    def test_one_record_matches_packaged_manifest(self):
        self.assertEqual(len(load_evidence()["topics"]), 78)
        record = get_evidence("fractional")
        self.assertEqual(record["slug"], "fractional")
        self.assertIn("S", " ".join(record["source_ids"]))

    def test_json_cli(self):
        out = StringIO()
        with redirect_stdout(out):
            self.assertEqual(main(["evidence", "lambda-calculus", "--json"]), 0)
        data = json.loads(out.getvalue())
        self.assertEqual(data["slug"], "lambda-calculus")
        self.assertIn("formal", data["kind"])

    def test_human_cli(self):
        out = StringIO()
        with redirect_stdout(out):
            self.assertEqual(main(["evidence", "time-scales"]), 0)
        self.assertIn("scope:", out.getvalue())
        self.assertIn("formal proof: False", out.getvalue())

    def test_unknown_slug(self):
        with self.assertRaises(KeyError):
            get_evidence("not-a-calculus")


if __name__ == "__main__":
    unittest.main()
