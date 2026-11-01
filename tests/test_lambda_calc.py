"""Independent witnesses for untyped lambda-calculus substitution and beta reduction."""

import unittest

from calc.lambda_calc import (
    Abs, App, LambdaSyntaxError, ReductionLimitExceeded, Var,
    alpha_equivalent, beta_step, format_term, free_vars, normalize, parse,
    substitute,
)


class LambdaSyntaxTests(unittest.TestCase):
    def test_identity_and_application_parsing(self):
        self.assertEqual(parse(r"(\x.x) y"), App(Abs("x", Var("x")), Var("y")))
        self.assertEqual(parse("λx.x"), Abs("x", Var("x")))

    def test_bad_expressions_are_rejected(self):
        for source in ("", "x; y", "x + y", r"\x x", "(x", "x)", r"\1.x",
                       "x.__class__", "__import__('os')", "λx."):
            with self.subTest(source=source):
                with self.assertRaises(LambdaSyntaxError):
                    parse(source)

    def test_resource_limit_for_parser(self):
        with self.assertRaises(LambdaSyntaxError):
            parse("x " * 2200)

    def test_format_round_trip_alpha_equivalence(self):
        term = parse(r"\x.(x (\y.y))")
        self.assertTrue(alpha_equivalent(term, parse(format_term(term))))


class LambdaSemanticsTests(unittest.TestCase):
    def test_free_vars_with_shadowing(self):
        self.assertEqual(free_vars(parse(r"\x.(x (\x.y))")), frozenset({"y"}))

    def test_identity_reduction(self):
        self.assertEqual(normalize(parse(r"(\x.x) y")), Var("y"))

    def test_k_combinator(self):
        self.assertEqual(normalize(parse(r"((\x.\y.x) a) b")), Var("a"))

    def test_alpha_equivalence_respects_binding(self):
        self.assertTrue(alpha_equivalent(parse(r"\x.x"), parse(r"\z.z")))
        self.assertFalse(alpha_equivalent(parse(r"\x.y"), parse(r"\x.x")))
        self.assertFalse(alpha_equivalent(parse(r"\x.\y.x"), parse(r"\x.\y.y")))

    def test_capture_avoiding_substitution(self):
        result = substitute(parse(r"\y.x"), "x", Var("y"))
        self.assertTrue(alpha_equivalent(result, parse(r"\z.y")))
        self.assertEqual(free_vars(result), frozenset({"y"}))

    def test_capture_avoiding_beta_reduction(self):
        result = normalize(parse(r"(\x.\y.x) y"))
        self.assertTrue(alpha_equivalent(result, parse(r"\z.y")))
        self.assertEqual(free_vars(result), frozenset({"y"}))

    def test_nested_capture_and_application(self):
        result = normalize(parse(r"(\x.\y.\z.x y) y"))
        self.assertTrue(alpha_equivalent(result, parse(r"\u.\z.y u")))
        self.assertEqual(free_vars(result), frozenset({"y"}))

    def test_substitution_respects_shadowing(self):
        self.assertEqual(substitute(parse(r"\x.x"), "x", Var("z")), parse(r"\x.x"))

    def test_normal_order_discards_divergent_argument(self):
        omega = r"(\z.z z) (\z.z z)"
        self.assertEqual(normalize(parse(r"(\x.a) (" + omega + ")"), max_steps=1), Var("a"))

    def test_omega_is_bounded(self):
        with self.assertRaises(ReductionLimitExceeded):
            normalize(parse(r"(\x.x x) (\x.x x)"), max_steps=8)

    def test_step_is_outermost_leftmost(self):
        term = parse(r"(\x.x) ((\y.y) z)")
        self.assertEqual(beta_step(term), parse(r"(\y.y) z"))

    def test_reduction_under_abstraction(self):
        self.assertEqual(normalize(parse(r"\x.((\y.y) x)")), parse(r"\x.x"))

    def test_budget_validation_and_normal_form_zero_budget(self):
        self.assertEqual(normalize(Var("x"), max_steps=0), Var("x"))
        with self.assertRaises(ValueError):
            normalize(Var("x"), max_steps=-1)
        with self.assertRaises(ReductionLimitExceeded):
            normalize(parse(r"(\x.x) y"), max_steps=0)

    def test_expanding_redex_refused_before_memory_growth(self):
        source = r"(\x.x x x x x) (\z.z z z z)"
        with self.assertRaisesRegex(ReductionLimitExceeded, "next reduction"):
            normalize(parse(source), max_nodes=22)


if __name__ == "__main__":
    unittest.main()
