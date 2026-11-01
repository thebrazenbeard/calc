"""Static distribution-contract checks for the zero-dependency explorer.

These do not claim to prove browser rendering or mathematical correctness.
"""

from html.parser import HTMLParser
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
EXPLORER = ROOT / "explorer"


class AssetParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tags = []

    def handle_starttag(self, tag, attributes):
        self.tags.append((tag, dict(attributes)))


class ExplorerAssetContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.html = (EXPLORER / "index.html").read_text(encoding="utf-8")
        cls.script = (EXPLORER / "app.js").read_text(encoding="utf-8")
        cls.styles = (EXPLORER / "styles.css").read_text(encoding="utf-8")
        cls.parser = AssetParser()
        cls.parser.feed(cls.html)

    def test_all_assets_exist_and_are_local(self):
        srcs = [attrs["src"] for tag, attrs in self.parser.tags
                if tag == "script" and "src" in attrs]
        hrefs = [attrs["href"] for tag, attrs in self.parser.tags
                 if tag == "link" and "href" in attrs]
        self.assertEqual(srcs, ["./app.js"])
        self.assertEqual(hrefs, ["./styles.css"])
        self.assertTrue((EXPLORER / "app.js").is_file())
        self.assertTrue((EXPLORER / "styles.css").is_file())

    def test_rejects_dynamic_execution_and_network_dependencies(self):
        self.assertNotRegex(self.script, r"\beval\s*\(")
        self.assertNotRegex(self.script, r"\bnew\s+Function\s*\(")
        self.assertNotRegex(self.script, r"\bfetch\s*\(")
        self.assertNotRegex(self.script, r"\bXMLHttpRequest\b")
        self.assertNotRegex(self.script, r"\binnerHTML\s*=")
        self.assertNotRegex(self.html, r"(?i)<script[^>]+src=[\"']https?://")
        self.assertNotRegex(self.styles, r"(?i)@import\b")

    def test_controls_are_semantically_labelled(self):
        controls = {attrs["id"] for tag, attrs in self.parser.tags
                    if tag in ("input", "select", "button") and "id" in attrs}
        for key in ("function-select", "view-range", "point-slider", "integral-start",
                    "integral-end", "reset-button"):
            with self.subTest(control=key):
                self.assertIn(key, controls)
        labels = {attrs["for"] for tag, attrs in self.parser.tags
                  if tag == "label" and "for" in attrs}
        self.assertTrue({"function-select", "view-range", "point-slider",
                         "integral-start", "integral-end"}.issubset(labels))

    def test_keyboard_and_reduced_motion_contracts(self):
        self.assertIn('event.key === "ArrowLeft"', self.script)
        self.assertIn('event.key === "ArrowRight"', self.script)
        self.assertIn("pointerdown", self.script)
        self.assertIn("focus-visible", self.styles)
        self.assertIn("prefers-reduced-motion: reduce", self.styles)
        self.assertIn('tabindex="0"', self.html)
        self.assertIn("<noscript>", self.html)

    def test_reference_and_estimated_values_are_separate(self):
        for key in ("height-output", "derivative-output", "derivative-estimate",
                    "integral-output", "integral-estimate", "input-error"):
            with self.subTest(output=key):
                self.assertIn(f'id="{key}"', self.html)
        self.assertIn("centralDifference", self.script)
        self.assertIn("simpson", self.script)
        self.assertIn("const slices = 600", self.script)
        self.assertIn("derivative: x =>", self.script)
        self.assertIn("primitive: x =>", self.script)

    def test_bounded_input_and_interval_error_contract(self):
        self.assertIn('max="8"', self.html)
        self.assertIn('min="2"', self.html)
        self.assertIn("a < b", self.script)
        self.assertIn("Number.isFinite(a)", self.script)
        self.assertIn("Invalid interval", self.script)

    def test_no_external_assets_or_generated_content_claims(self):
        self.assertNotRegex(self.html, r"(?i)<iframe\b")
        self.assertNotRegex(self.html, r"(?i)<img[^>]+src=[\"']https?://")
        self.assertIn("No expression execution", self.html)
        self.assertIn("not</em> certified error bounds", self.html)


if __name__ == "__main__":
    unittest.main()
