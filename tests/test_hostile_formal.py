"""Independent falsifiers for formal lambda and small symbolic calculus.

These tests use mathematical binding invariants and analytic differentiation
identities, not a copy of either implementation's algorithms.
"""

import math
import unittest

from calc.lambda_calc import (
    Abs, App, ReductionLimitExceeded, Var as LambdaVar,
    alpha_equivalent, format_term, free_vars, normalize, parse, substitute,
)
from calc.symbolic import Var as SymbolicVar, log, parse_expression, sin, exp


class HostileLambdaContracts(unittest.TestCase):
    def test_manual_deep_ast_respects_node_budget_before_python_recursion(self):
        """A size limit should reject a giant tree, not leak RecursionError."""
        term = LambdaVar("x")
        for _ in range(1200):
            term = App(LambdaVar("f"), term)
        with self.assertRaises(ReductionLimitExceeded):
            normalize(term, max_nodes=10)
        # An overly generous node budget must not bypass the depth guard.
        with self.assertRaises(ReductionLimitExceeded):
            normalize(term, max_nodes=100_000)

    def test_nested_capture_avoidance_and_free_variable_preservation(self):
        term = parse(r"\y.\z.x y")
        result = substitute(term, "x", LambdaVar("y"))
        self.assertEqual(free_vars(result), frozenset({"y"}))
        self.assertTrue(alpha_equivalent(result, parse(r"\a.\z.y a")))

    def test_alpha_equivalence_distinguishes_binder_depth(self):
        self.assertTrue(alpha_equivalent(parse(r"\x.\y.x y"),
                                         parse(r"\a.\b.a b")))
        self.assertFalse(alpha_equivalent(parse(r"\x.\y.x y"),
                                          parse(r"\a.\b.b a")))

    def test_normal_order_ignores_unneeded_divergent_argument(self):
        omega = r"(\z.z z) (\z.z z)"
        term = parse(r"(\x.\y.y) (" + omega + ")")
        self.assertTrue(alpha_equivalent(normalize(term, max_steps=1), parse(r"\q.q")))

    def test_formatter_round_trips_nested_scopes(self):
        cases = (
            r"(\x.x) (\y.y)",
            r"\x.\y.x (\z.z)",
            r"f (\x.x y) z",
            r"(\x.\y.x y) a",
        )
        for source in cases:
            with self.subTest(source=source):
                original = parse(source)
                rendered = format_term(original)
                self.assertTrue(alpha_equivalent(parse(rendered), original))

    def test_church_numeral_addition_reduces_to_five(self):
        """A nontrivial beta-reduction chain has a known normal form."""
        add_two_three = parse(
            r"(\m.\n.\f.\x.m f (n f x)) "
            r"(\f.\x.f (f x)) "
            r"(\f.\x.f (f (f x)))"
        )
        five = parse(r"\f.\x.f (f (f (f (f x))))")
        self.assertTrue(alpha_equivalent(normalize(add_two_three, max_steps=100), five))


class HostileSymbolicContracts(unittest.TestCase):
    def test_domain_is_not_erased_by_zero_multiplication(self):
        x = SymbolicVar("x")
        for expression, point, error in (
            (0 * log(x), -1.0, ValueError),
            (0 * (1 / x), 0.0, ZeroDivisionError),
            ((1 / x) ** 0, 0.0, ZeroDivisionError),
        ):
            with self.subTest(expression=str(expression)):
                with self.assertRaises(error):
                    expression.simplify().eval({"x": point})

    def test_chain_product_derivative_independent_analytic_witness(self):
        x = SymbolicVar("x")
        expression = exp(x) * sin(x**2)
        derivative = expression.diff("x")
        for point in (-1.25, -0.3, 0.0, 0.65, 1.5):
            expected = math.exp(point) * (math.sin(point**2) + 2*point*math.cos(point**2))
            with self.subTest(point=point):
                self.assertAlmostEqual(derivative.eval({"x": point}), expected, delta=1e-10)

    def test_parser_refuses_dynamic_python_constructs(self):
        for source in (
            "__import__('math')",
            "x.__class__",
            "(lambda x: x)(2)",
            "[x for x in (1,2)]",
            "sin.__call__(x)",
            "x[0]",
        ):
            with self.subTest(source=source), self.assertRaises(ValueError):
                parse_expression(source)

    def test_symbolic_power_still_respects_real_log_domain(self):
        expression = parse_expression("log(x)**2")
        self.assertAlmostEqual(expression.diff("x").eval({"x": math.e}), 2/math.e)
        with self.assertRaises(ValueError):
            expression.eval({"x": -2})


if __name__ == "__main__":
    unittest.main()
