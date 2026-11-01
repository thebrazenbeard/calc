"""Static CI contract checks; these do not replace a real GitHub Actions run."""

from pathlib import Path
import re
import unittest


WORKFLOW = Path(__file__).resolve().parents[1] / ".github" / "workflows" / "test.yml"


class WorkflowContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = WORKFLOW.read_text(encoding="utf-8")

    def test_push_pull_request_and_manual_trigger(self):
        for event in ("push", "pull_request", "workflow_dispatch"):
            with self.subTest(event=event):
                self.assertRegex(self.text, rf"(?m)^  {re.escape(event)}:\s*$")
        self.assertNotIn("pull_request_target", self.text)

    def test_minimum_token_scope_and_untrusted_checkout(self):
        self.assertRegex(self.text, r"(?m)^permissions:\s*\n  contents: read\s*$")
        self.assertIn("persist-credentials: false", self.text)
        self.assertNotIn("secrets.", self.text)

    def test_supported_python_matrix_and_repository_tests(self):
        for version in ('"3.11"', '"3.12"', '"3.13"'):
            with self.subTest(version=version):
                self.assertIn(version, self.text)
        self.assertIn("PYTHONPATH=src python -m unittest discover -s tests -v", self.text)
        self.assertIn("python -m compileall -q src", self.text)

    def test_wheel_and_isolated_installed_import_smoke(self):
        self.assertIn("python -m pip wheel --no-deps --no-build-isolation", self.text)
        self.assertIn("python -m venv", self.text)
        self.assertIn("--no-index --no-deps", self.text)
        self.assertIn("RUNNER_TEMP", self.text)
        self.assertIn("calc.evidence", self.text)
        self.assertIn("calc-atlas", self.text)


if __name__ == "__main__":
    unittest.main()
