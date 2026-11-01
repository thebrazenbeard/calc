"""Independent adversarial witnesses for the Wave 5 modules.

The oracles are exact rational arithmetic, hand-derived identities, or the
package/distribution contract. In particular, a test that currently fails is
kept as a regression rather than weakening its mathematical assertion.
"""

from __future__ import annotations

from fractions import Fraction
from importlib import import_module, resources
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tomllib
import unittest

from calc.evidence import get_evidence, load_evidence
from calc.fractional import (
    caputo_power_derivative,
    grunwald_letnikov_derivative,
    riemann_liouville_power_derivative,
)
from calc.grobner import (
    MultiPolynomial,
    ResourceLimitError,
    grobner_basis,
    leading_term,
    normal_form,
    s_polynomial,
    variables,
)
from calc.intervals import Interval
from calc.matrix_calc import (
    frobenius_inner,
    matrix_gradient,
    matrix_jvp,
    squared_frobenius_gradient,
    trace_axb_gradient,
)
from calc.time_scales import FiniteTimeScale


REPO_ROOT = Path(__file__).resolve().parents[1]


class HostileMatrixCalculusTests(unittest.TestCase):
    def test_exact_trace_gradient_orientation_by_hand(self):
        a = ((1, 2, 3), (-1, 4, 0))
        b = ((2, 3), (-1, 5), (7, -2), (0, 1))
        self.assertEqual(
            trace_axb_gradient(a, b),
            ((-1, -6, 9, -1), (16, 18, 6, 4), (6, -3, 21, 0)),
        )

    def test_exact_frobenius_rational_oracle(self):
        a = ((Fraction(1, 3), Fraction(-2, 5)), (Fraction(7, 11), 2))
        b = ((Fraction(5, 7), Fraction(3, 2)), (Fraction(-11, 7), -3))
        expected = Fraction(5, 21) - Fraction(3, 5) - 1 - 6
        self.assertEqual(frobenius_inner(a, b), expected)
        self.assertEqual(squared_frobenius_gradient(a)[0][0], Fraction(2, 3))

    def test_gradient_of_linear_map_not_rejected_on_intermediate_overflow(self):
        # f(+1)=+1e308 and f(-1)=-1e308 are both finite. Their subtraction
        # overflows binary64, but (f(+1)-f(-1))/2 = 1e308 is finite.
        result = matrix_gradient(lambda x: 1e308 * x[0][0], ((0.0,),), step=1.0)
        self.assertEqual(result, ((1e308,),))

    def test_jvp_of_linear_map_not_rejected_on_intermediate_overflow(self):
        result = matrix_jvp(
            lambda x: ((1e308 * x[0][0],),),
            ((0.0,),),
            ((1.0,),),
            step=1.0,
        )
        self.assertEqual(result, ((1e308,),))

    def test_frobenius_cancellation_of_individually_overflowing_products(self):
        # Exact binary64 inputs have products whose sum cancels to zero.
        self.assertEqual(
            frobenius_inner(((1e308, 1e308),), ((2.0, -2.0),)), 0.0
        )

    def test_trace_gradient_cancellation_of_overflowed_dot_products(self):
        # ∇ tr(A X B) = (BA)^T, with BA=[2e308 - 2e308]=[0].
        self.assertEqual(
            trace_axb_gradient(((1e308,), (1e308,)), ((2.0, -2.0),)),
            ((0.0,),),
        )


class HostileGrobnerTests(unittest.TestCase):
    def test_nontrivial_elimination_basis_independently_derived(self):
        # x^2=y and xy=1 imply y^3=1 and x=y^2. Lex x>y thus has
        # reduced basis {x-y^2, y^3-1}, an independent algebraic witness.
        x, y = variables(2)
        basis = grobner_basis((x**2 - y, x * y - 1), order="lex")
        self.assertEqual(basis, (x - y**2, y**3 - 1))
        self.assertEqual(tuple(leading_term(g)[0] for g in basis), ((1, 0), (0, 3)))

    def test_rational_coefficient_gcd_over_q(self):
        (x,) = variables(1)
        b = grobner_basis((2 * x**2 - 2, 3 * x - 3))
        self.assertEqual(b, (x - 1,))
        self.assertEqual((Fraction(2, 3) * x**2 + Fraction(1, 2)).evaluate((3,)),
                         Fraction(13, 2))

    def test_seeded_common_linear_factor_univariate_ideals(self):
        # Over Q[x], <(x-a)(x-b), (x-a)(x-c)> = <x-a> if b != c.
        (x,) = variables(1)
        sample = random.Random(10073)
        for _ in range(35):
            a, b, c = sample.sample(range(-11, 12), 3)
            with self.subTest(root=a, first=b, second=c):
                self.assertEqual(
                    grobner_basis(((x-a)*(x-b), (x-a)*(x-c))),
                    (x-a,),
                )

    def test_grlex_buchberger_pairs_reduce_to_zero(self):
        x, y = variables(2)
        generators = (x**2 - y, x * y - 1)
        basis = grobner_basis(generators, order="grlex")
        for poly in generators:
            self.assertEqual(normal_form(poly, basis, order="grlex"),
                             MultiPolynomial.zero(2))
        for i in range(len(basis)):
            for j in range(i + 1, len(basis)):
                self.assertEqual(
                    normal_form(s_polynomial(basis[i], basis[j], "grlex"),
                                basis, order="grlex"),
                    MultiPolynomial.zero(2),
                )

    def test_declared_polynomial_resource_limits_are_enforced(self):
        with self.assertRaises(ResourceLimitError):
            MultiPolynomial({(513,): 1}, nvars=1)
        (x,) = variables(1)
        with self.assertRaises(ResourceLimitError):
            x**513
        with self.assertRaises(ResourceLimitError):
            grobner_basis((x**2 - 1, x - 1, x**3 - 1), max_pairs=1)


class HostileTimeScalesTests(unittest.TestCase):
    def test_finite_difference_exactly_two_despite_subtraction_overflow(self):
        grid = FiniteTimeScale((0.0, 1e308))
        function = lambda t: -1e308 if t == 0.0 else 1e308
        self.assertEqual(grid.delta_derivative(function, 0.0), 2.0)
        self.assertEqual(grid.nabla_derivative(function, 1e308), 2.0)

    def test_irregular_time_scale_fundamental_theorem_discrete(self):
        grid = FiniteTimeScale((1.0, 2.0, 4.0, 7.0))
        function = lambda t: t * t
        increments = (4 - 1, 16 - 4, 49 - 16)
        self.assertEqual(
            grid.delta_integral(lambda t: grid.delta_derivative(function, t), 1.0, 7.0),
            sum(increments),
        )
        self.assertEqual(
            grid.nabla_integral(lambda t: grid.nabla_derivative(function, t), 1.0, 7.0),
            sum(increments),
        )

    def test_delta_integral_cancellation_of_overflowed_terms(self):
        grid = FiniteTimeScale((0.0, 2.0, 4.0))
        self.assertEqual(
            grid.delta_integral(lambda t: 1e308 if t == 0.0 else -1e308,
                                0.0, 4.0),
            0.0,
        )

    def test_nabla_integral_cancellation_of_overflowed_terms(self):
        grid = FiniteTimeScale((0.0, 2.0, 4.0))
        self.assertEqual(
            grid.nabla_integral(lambda t: 1e308 if t == 2.0 else -1e308,
                                0.0, 4.0),
            0.0,
        )


class HostileFractionalTests(unittest.TestCase):
    def test_half_derivative_gamma_closed_forms(self):
        self.assertAlmostEqual(caputo_power_derivative(0.5, 1, 4),
                               4 / math.sqrt(math.pi), places=12)
        self.assertAlmostEqual(riemann_liouville_power_derivative(0.5, 0, 4),
                               1 / (2 * math.sqrt(math.pi)), places=12)
        self.assertEqual(caputo_power_derivative(0.5, 0, 4), 0.0)
        self.assertAlmostEqual(caputo_power_derivative(1, 3, 2), 12.0, places=12)

    def test_gl_approaches_nonzero_riemann_liouville_constant(self):
        observed = grunwald_letnikov_derivative(lambda t: 1.0,
                                                1.0, 0.5, steps=8192)
        self.assertAlmostEqual(observed, 1 / math.sqrt(math.pi), delta=0.001)
        self.assertNotAlmostEqual(observed, 0.0, delta=0.1)

    def test_zero_endpoint_singularity_is_explicit(self):
        with self.assertRaises(ValueError):
            riemann_liouville_power_derivative(0.5, 0, 0)
        self.assertEqual(caputo_power_derivative(0.5, 2, 0), 0.0)


class HostileIntervalTests(unittest.TestCase):
    @staticmethod
    def _encases(interval: Interval, exact_lower: Fraction, exact_upper: Fraction):
        if math.isfinite(interval.lower):
            assert Fraction.from_float(interval.lower) <= exact_lower
        if math.isfinite(interval.upper):
            assert Fraction.from_float(interval.upper) >= exact_upper

    def test_seeded_interval_arithmetic_exact_rational_extrema(self):
        sample = random.Random(28091)
        candidates = (-1e308, -1e150, -7.5, -0.1, -1e-307, -0.0,
                      0.0, 1e-307, 0.1, 3.25, 1e150, 1e308)
        for _ in range(180):
            l0, l1 = sorted(sample.choices(candidates, k=2))
            r0, r1 = sorted(sample.choices(candidates, k=2))
            left, right = Interval(l0, l1), Interval(r0, r1)
            a, b, c, d = map(Fraction.from_float, (l0, l1, r0, r1))
            with self.subTest(left=left, right=right):
                self._encases(left + right, a + c, b + d)
                self._encases(left - right, a - d, b - c)
                products = (a * c, a * d, b * c, b * d)
                self._encases(left * right, min(products), max(products))
                if not (c <= 0 <= d):
                    quotients = (a / c, a / d, b / c, b / d)
                    self._encases(left / right, min(quotients), max(quotients))

    def test_exact_binary64_point_and_overflow_enclosure(self):
        p = Interval.point(0.1)
        result = p * p
        expected = Fraction.from_float(0.1) ** 2
        self._encases(result, expected, expected)
        self.assertEqual(Interval.point(1e308) * Interval.point(2.0),
                         Interval(sys.float_info.max, math.inf))

    def test_interval_division_excluding_zero_denominator(self):
        with self.assertRaises(ZeroDivisionError):
            Interval(-1.0, 2.0) / Interval(-0.0, 3.0)
        self.assertIsNone(Interval(0.0, 1.0).intersection(Interval(2.0, 3.0)))


class HostileEvidenceManifestTests(unittest.TestCase):
    def test_packaged_resource_is_json_and_explicitly_declared(self):
        raw = resources.files("calc").joinpath("topic_manifest.json").read_text()
        self.assertEqual(json.loads(raw), load_evidence())
        metadata = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
        self.assertIn("topic_manifest.json",
                      metadata["tool"]["setuptools"]["package-data"]["calc"])

    def test_implemented_topic_bindings_and_paths_are_real(self):
        topics = load_evidence()["topics"]
        slugs = [record["slug"] for record in topics]
        self.assertEqual(len(slugs), len(set(slugs)))
        self.assertGreaterEqual(len(slugs), 70)
        source_ids = set()
        for line in (REPO_ROOT / "docs" / "SOURCES.md").read_text().splitlines():
            if line.startswith("| S"):
                source_ids.add(line.split("|")[1].strip())
        for record in topics:
            with self.subTest(slug=record["slug"]):
                self.assertIn(record["evidence_level"], ("indexed_only", "example_tested"))
                self.assertEqual(record["formal_proof"], False)
                self.assertTrue(record["source_ids"])
                self.assertTrue(set(record["source_ids"]) <= source_ids)
                for name in ("guide",):
                    relative_path = record.get(name)
                    if relative_path:
                        path = REPO_ROOT / relative_path
                        self.assertTrue(path.is_file(), path)
                        self.assertTrue(path.resolve().is_relative_to(REPO_ROOT))
                if record["evidence_level"] == "example_tested":
                    self.assertEqual(record["status"], "implemented")
                    self.assertTrue(record["bindings"])
                    self.assertTrue(record["test_paths"])
                    for binding in record["bindings"]:
                        module = import_module(binding["module"])
                        for symbol in binding["symbols"]:
                            self.assertTrue(hasattr(module, symbol), (module, symbol))
                    for relative_path in record["test_paths"]:
                        path = REPO_ROOT / relative_path
                        self.assertTrue(path.is_file(), path)
                        self.assertTrue(path.resolve().is_relative_to(REPO_ROOT))

    def test_evidence_cli_returns_exact_packaged_record(self):
        record = get_evidence("time-scales")
        command = subprocess.run(
            [sys.executable, "-m", "calc", "evidence", "time-scales", "--json"],
            capture_output=True, text=True, check=False, cwd=REPO_ROOT,
        )
        self.assertEqual(command.returncode, 0, command.stderr)
        self.assertEqual(json.loads(command.stdout), record)


if __name__ == "__main__":
    unittest.main()
