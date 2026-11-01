"""Independent mathematical and parser-boundary tests for symbolic expressions."""

import math
import unittest
from fractions import Fraction

from calc.symbolic import Const, Var, cos, exp, log, parse_expression, sin


class SymbolicCalculusTests(unittest.TestCase):
    def test_polynomial_rules(self):
        x = Var("x")
        f = x**3 + 2 * x**2 - 7 * x + 5
        self.assertAlmostEqual(f.eval({"x": 2}), 7)
        self.assertAlmostEqual(f.diff("x").eval({"x": 2}), 13)
        self.assertAlmostEqual(f.diff("x").diff("x").eval({"x": 2}), 16)

    def test_chain_product_and_quotient_rules(self):
        x = Var("x")
        f = sin(x * x) / (1 + x)
        point = 0.4
        expected = ((2 * point * math.cos(point**2) * (1 + point))
                    - math.sin(point**2)) / (1 + point)**2
        self.assertAlmostEqual(f.diff("x").eval({"x": point}), expected, places=12)

    def test_log_exp_and_trigonometry(self):
        x = Var("x")
        self.assertAlmostEqual(log(1 + x*x).diff("x").eval({"x": 2}), 0.8, places=12)
        self.assertAlmostEqual(exp(sin(x)).diff("x").eval({"x": 0}), 1.0, places=12)
        self.assertAlmostEqual(cos(x).diff("x").eval({"x": math.pi / 2}), -1.0, places=12)

    def test_independent_variables(self):
        x, y = Var("x"), Var("y")
        f = x*x + 3*y
        self.assertAlmostEqual(f.diff("y").eval({"x": 4, "y": 10}), 3.0)
        self.assertAlmostEqual(f.diff("z").eval({"x": 4, "y": 10}), 0.0)

    def test_minimal_simplification(self):
        x = Var("x")
        self.assertEqual(((x + 0) * 1).simplify(), x)
        self.assertEqual(((Const(2) + 3) * 4).simplify(), Const(20))
        self.assertEqual((x**1).simplify(), x)
        self.assertEqual((x**0).simplify(), Const(1))
        self.assertEqual(Const(2)**-2, Const(Fraction(1, 4)))
        self.assertIsInstance((Const(2)**-2).value, Fraction)

    def test_parser_math_matches_manually_derived_values(self):
        expression = parse_expression("sin(x**2) + exp(-x) / (1+x)")
        point = 0.6
        expected = 2 * point * math.cos(point**2) - math.exp(-point) * (point + 2) / (1 + point)**2
        self.assertAlmostEqual(expression.diff("x").eval({"x": point}), expected, places=11)
        self.assertAlmostEqual(parse_expression("x**-2").eval({"x": 2}), 0.25)

    def test_parser_refuses_execution_and_dynamic_syntax(self):
        bad = (
            "__import__('os').system('echo unsafe')", "open('secret')", "sin.__class__",
            "lambda x: x", "[x for x in range(3)]", "f'{x}'", "True",
            "x if x else y", "x[0]", "sin(x, y)", "x**10000", "x**x",
            "x.__dict__", "sin(x=2)", "1e309", "[1,2]", "1 < 2",
        )
        for source in bad:
            with self.subTest(source=source), self.assertRaises((ValueError, SyntaxError)):
                parse_expression(source)

    def test_parser_resource_limits(self):
        with self.assertRaises(ValueError):
            parse_expression("x + " * 700 + "x")
        with self.assertRaises(ValueError):
            parse_expression("sin(" * 50 + "x" + ")" * 50)

    def test_domain_input_and_arithmetic_errors(self):
        x = Var("x")
        with self.assertRaises(KeyError):
            x.eval({})
        with self.assertRaises(ValueError):
            x.eval({"x": float("nan")})
        with self.assertRaises(ValueError):
            log(x).eval({"x": -1})
        with self.assertRaises(ZeroDivisionError):
            (1 / x).eval({"x": 0})
        with self.assertRaises(ZeroDivisionError):
            ((1 / x)**0).eval({"x": 0})
        with self.assertRaises(ValueError):
            parse_expression("sin(x", max_nodes=20)

    def test_explicit_power_constraints(self):
        x = Var("x")
        with self.assertRaises(ValueError):
            _ = x**1.5
        with self.assertRaises(ValueError):
            _ = x**200
        with self.assertRaises(TypeError):
            _ = x**True


if __name__ == "__main__":
    unittest.main()
