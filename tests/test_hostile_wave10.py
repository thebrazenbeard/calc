"""Independent Wave 10 adversarial witnesses.

These checks target public contracts and mathematical invariants. Passing
them is not a proof of general numerical convergence or formal correctness.
"""

from __future__ import annotations

import unittest
from fractions import Fraction
import random
import re
from html.parser import HTMLParser
from pathlib import Path

from calc.algebra import Polynomial
from calc.catalog import list_topics
from calc.evidence import get_evidence, load_evidence
from calc.integrals import (definite_integral, integration_by_parts,
                            polynomial_antiderivative, substitution_integral)
from calc.limits import polynomial_limit_certificate, one_sided_samples
from calc.optimization import bisection, safeguarded_newton

ROOT = Path(__file__).resolve().parents[1]


class EvidenceIsolationReview(unittest.TestCase):
    def test_mutating_nested_record_does_not_poison_subsequent_lookup(self):
        original = get_evidence("differential")
        expected = tuple(binding["module"] for binding in original["bindings"])
        original["bindings"].append({"module": "tampered", "symbols": []})
        original["source_ids"].append("FAKE")
        original["test_paths"].clear()
        again = get_evidence("differential")
        self.assertEqual(tuple(binding["module"] for binding in again["bindings"]), expected)
        self.assertNotIn("FAKE", again["source_ids"])
        self.assertTrue(again["test_paths"])

    def test_mutating_public_manifest_does_not_poison_catalogue(self):
        first = load_evidence()
        original_count = len(first["topics"])
        first["topics"].clear()
        first["schema_version"] = -1
        second = load_evidence()
        self.assertEqual(second["schema_version"], 1)
        self.assertEqual(len(second["topics"]), original_count)

    def test_implemented_is_bounded_examples_not_formal_proof(self):
        records = {entry["slug"]: entry for entry in load_evidence()["topics"]}
        topics = list_topics()
        self.assertEqual({topic.slug for topic in topics}, set(records))
        for topic in topics:
            record = records[topic.slug]
            with self.subTest(slug=topic.slug):
                self.assertEqual(record["status"], topic.status)
                self.assertIs(record["formal_proof"], False)
                if topic.status == "implemented":
                    self.assertEqual(record["evidence_level"], "example_tested")
                    self.assertTrue(record["test_paths"])
                    self.assertTrue(record["bindings"])


class PolynomialIntegrationReview(unittest.TestCase):
    def test_power_rule_against_independent_exact_endpoint_formula(self):
        for n in range(0, 33):
            p = Polynomial((0,) * n + (1,))
            for a, b in ((-2, 3), (-1, 1), (Fraction(1, 3), Fraction(5, 2))):
                with self.subTest(power=n, start=a, stop=b):
                    expected = (Fraction(b) ** (n + 1) -
                                Fraction(a) ** (n + 1)) / (n + 1)
                    self.assertEqual(definite_integral(p, a, b), expected)
                    self.assertEqual(definite_integral(p, b, a), -expected)

    def test_primitive_derivative_and_constant_fixed(self):
        p = Polynomial((Fraction(1, 3), -2, 0, Fraction(5, 7)))
        primitive = polynomial_antiderivative(p, constant=Fraction(-11, 13))
        self.assertEqual(primitive.differentiate(), p)
        self.assertEqual(primitive(0), Fraction(-11, 13))

    def test_substitution_integral_with_nonmonotone_inner(self):
        # p(y)=3y²+1, q(x)=x²+x+1; P(y)=y³+y.
        p = Polynomial((1, 0, 3))
        q = Polynomial((1, 1, 1))
        self.assertEqual(substitution_integral(p, q, -1, 2), Fraction(348))
        self.assertEqual(substitution_integral(p, q, 2, -1), Fraction(-348))
        # q(-1)=q(0)=1: oriented integral vanishes despite nonzero integrand.
        self.assertEqual(substitution_integral(p, q, -1, 0), Fraction(0))

    def test_integration_by_parts_nonconstant_with_exact_oracle(self):
        # u=x³+2, v=x²-1 -> integral of 3x²(x²-1) on [0,2] = 56/5.
        u = Polynomial((2, 0, 0, 1))
        v = Polynomial((-1, 0, 1))
        self.assertEqual(integration_by_parts(u, v, 0, 2), Fraction(56, 5))
        self.assertEqual(integration_by_parts(u, v, 2, 0), Fraction(-56, 5))

    def test_exact_input_contract_rejects_inexact_and_bool(self):
        p = Polynomial((1, 2, 3))
        for value in (True, 1.25, "2"):
            with self.subTest(value=value), self.assertRaises(TypeError):
                definite_integral(p, value, 2)


class LimitsContinuityReview(unittest.TestCase):
    def test_seeded_polynomial_epsilon_delta_certificates_against_direct_polynomials(self):
        rng = random.Random(20261101)
        for _ in range(60):
            coefficients = tuple(Fraction(rng.randint(-7, 7), rng.randint(1, 9))
                                 for _ in range(rng.randint(1, 7)))
            c = Fraction(rng.randint(-3, 3), rng.randint(1, 6))
            epsilon = Fraction(rng.randint(1, 7), rng.randint(1, 10))
            radius = Fraction(rng.randint(1, 4), rng.randint(1, 5))
            cert = polynomial_limit_certificate(coefficients, c, epsilon, radius=radius)
            self.assertGreater(cert.delta, 0)
            self.assertLessEqual(cert.delta, radius)
            independent = lambda x: sum((a * x**k for k, a in enumerate(coefficients)),
                                        Fraction(0))
            self.assertEqual(cert.value, independent(c))
            for part in (Fraction(-3, 4), Fraction(-1, 2), Fraction(1, 3),
                         Fraction(2, 3), Fraction(9, 10)):
                x = c + cert.delta * part
                bound = cert.error_bound_at(x)
                error = abs(independent(x) - independent(c))
                self.assertLessEqual(error, bound)
                self.assertLess(error, epsilon)

    def test_one_sided_samples_are_diagnostics_not_continuity_proof(self):
        # A function agreeing on both punctured sides could have any value
        # at the center. Sampling there cannot establish continuity.
        def removable_discontinuity(x):
            if x == 0:
                raise AssertionError("callback must not evaluate at the center")
            return 0

        observations = one_sided_samples(removable_discontinuity, 0, count=12)
        self.assertEqual(observations.left, (0,) * 12)
        self.assertEqual(observations.right, (0,) * 12)
        # And a jump is observable from independent signs without a limit verdict.
        jump = one_sided_samples(lambda x: -2 if x < 0 else 3, 0, count=8)
        self.assertEqual(jump.left, (-2,) * 8)
        self.assertEqual(jump.right, (3,) * 8)

    def test_certificate_rejects_outside_delta_and_invalid_parameters(self):
        c = polynomial_limit_certificate((1, 2, 1), 0, Fraction(1, 100))
        with self.assertRaises(ValueError):
            c.error_bound_at(c.point + c.delta)
        for epsilon in (0, -1):
            with self.subTest(epsilon=epsilon), self.assertRaises(ValueError):
                polynomial_limit_certificate((1, 2), 0, epsilon)
        for invalid in (float("inf"), 0, -1):
            with self.subTest(step=invalid), self.assertRaises(ValueError):
                one_sided_samples(lambda x: x, 0, initial_step=invalid)


class _ExplorerHTML(HTMLParser):
    def __init__(self):
        super().__init__()
        self.attributes: list[tuple[str, dict[str, str]]] = []

    def handle_starttag(self, tag, attrs):
        self.attributes.append((tag, dict(attrs)))


class ExplorerAndWorkflowReview(unittest.TestCase):
    def test_explorer_uses_no_expression_execution_or_network_sink(self):
        script = (ROOT / "explorer" / "app.js").read_text(encoding="utf-8")
        markup = (ROOT / "explorer" / "index.html").read_text(encoding="utf-8")
        css = (ROOT / "explorer" / "styles.css").read_text(encoding="utf-8")
        for forbidden in (r"\beval\s*\(", r"\bnew\s+Function\s*\(",
                          r"\binnerHTML\b", r"\bfetch\s*\(",
                          r"\bXMLHttpRequest\b", r"\bWebSocket\b",
                          r"\bsendBeacon\s*\("):
            with self.subTest(forbidden=forbidden):
                self.assertNotRegex(script, forbidden)
        self.assertNotRegex(markup, r"(?i)<script[^>]+src\s*=\s*['\"]https?://")
        self.assertNotRegex(css, r"(?i)url\s*\(\s*['\"]?https?://")
        self.assertIn("Object.freeze", script)
        self.assertIn("textContent", script)
        self.assertIn("createElementNS", script)

    def test_explorer_controls_have_machine_checkable_accessibility_links(self):
        doc = _ExplorerHTML()
        doc.feed((ROOT / "explorer" / "index.html").read_text(encoding="utf-8"))
        ids = [a["id"] for _, a in doc.attributes if "id" in a]
        self.assertEqual(len(ids), len(set(ids)), "duplicate HTML id breaks label targets")
        label_targets = {a["for"] for tag, a in doc.attributes if tag == "label" and "for" in a}
        required = {"function-select", "view-range", "point-slider", "integral-start", "integral-end"}
        self.assertTrue(required <= set(ids))
        self.assertTrue(required <= label_targets)
        plots = [a for tag, a in doc.attributes if tag == "svg" and a.get("id") == "plot"]
        self.assertEqual(len(plots), 1)
        self.assertEqual(plots[0].get("role"), "img")
        self.assertTrue(plots[0].get("aria-label"))

    def test_workflow_source_contract_and_boundaries(self):
        workflow = (ROOT / ".github" / "workflows" / "test.yml").read_text(encoding="utf-8")
        for trigger in ("push", "pull_request", "workflow_dispatch"):
            self.assertRegex(workflow, rf"(?m)^  {trigger}:\s*$")
        self.assertRegex(workflow, r"(?m)^permissions:\s*\n  contents: read\s*$")
        self.assertIn("persist-credentials: false", workflow)
        self.assertIn("--no-index --no-deps", workflow)
        self.assertIn("calc.evidence", workflow)
        self.assertNotIn("pull_request_target", workflow)
        self.assertNotIn("secrets.", workflow)


class RootSolverRegressionReview(unittest.TestCase):
    def test_smooth_root_reaches_requested_function_residual(self):
        fn = lambda x: x*x - 2.0
        derivative = lambda x: 2.0*x
        for name, result in (("bisection", bisection(fn, 0.0, 2.0)),
                             ("Newton", safeguarded_newton(fn, derivative, 0.0, 2.0))):
            with self.subTest(method=name):
                self.assertGreaterEqual(result, 0.0)
                self.assertLessEqual(result, 2.0)
                self.assertLessEqual(abs(fn(result)), 1e-12)


if __name__ == "__main__":
    unittest.main()
