import unittest
from calc.optimization import ConvergenceError, golden_section_minimize


class IterationBoundaryTests(unittest.TestCase):
    def test_last_permitted_iteration_can_satisfy_tolerance(self):
        for objective, endpoint in ((lambda x: x, 0.0), (lambda x: -x, 1.0)):
            with self.subTest(endpoint=endpoint):
                try:
                    result = golden_section_minimize(objective, 0, 1, abs_tol=0.7, max_iter=1)
                except ConvergenceError as exc:
                    self.fail(f"One reduction is enough to satisfy width 0.7: {exc}")
                self.assertEqual(result, endpoint)

    def test_last_iteration_insufficient_still_refuses(self):
        with self.assertRaises(ConvergenceError):
            golden_section_minimize(lambda x: x, 0, 1, abs_tol=0.6, max_iter=1)
