"""Independent algebraic axioms and counterexamples for finite structures."""

import math
import unittest

from calc.structures import (
    FiniteGroup,
    PrimeField,
    boolean_complement,
    boolean_join,
    boolean_meet,
    boolean_xor,
    tropical_add,
    tropical_multiply,
)


class FiniteGroupTests(unittest.TestCase):
    def test_integers_mod_four(self):
        group = FiniteGroup([[(i + j) % 4 for j in range(4)] for i in range(4)])
        self.assertEqual(group.order, 4)
        self.assertEqual(group.identity, 0)
        self.assertTrue(group.is_abelian)
        self.assertEqual(group.mul(2, 3), 1)
        self.assertEqual(group.inverse(1), 3)
        self.assertEqual(group.element_order(2), 2)
        self.assertEqual(group.cyclic_generators(), (1, 3))

    def test_permutation_group_is_nonabelian(self):
        # S3: permutations of {0,1,2}, composition right-to-left.
        from itertools import permutations

        elements = list(permutations(range(3)))
        identity_index = elements.index((0, 1, 2))
        lookup = {p: i for i, p in enumerate(elements)}
        table = [[lookup[tuple(p[q[k]] for k in range(3))] for q in elements] for p in elements]
        group = FiniteGroup(table)
        self.assertEqual(group.identity, identity_index)
        self.assertFalse(group.is_abelian)
        self.assertEqual(group.order, 6)
        for x in range(6):
            self.assertEqual(group.mul(x, group.inverse(x)), identity_index)

    def test_reject_nonassociative_operation(self):
        # Identity and closure hold, yet associativity fails. This is not a group.
        malformed = [
            [0, 1, 2, 3],
            [1, 0, 3, 2],
            [2, 3, 0, 1],
            [3, 2, 1, 2],
        ]
        with self.assertRaisesRegex(ValueError, "associative"):
            FiniteGroup(malformed)

    def test_reject_no_identity_and_missing_inverse(self):
        with self.assertRaises(ValueError):
            FiniteGroup([[0, 0], [0, 0]])
        with self.assertRaises(ValueError):
            FiniteGroup([[0, 1], [1, 1]])

    def test_malformed_tables_and_invalid_elements(self):
        for table in ([], [[0]], [[0, 1]], [[0, 1], [1, 2]], [[0, 1], [1, True]]):
            if table == [[0]]:
                self.assertEqual(FiniteGroup(table).order, 1)
            else:
                with self.subTest(table=table), self.assertRaises((ValueError, TypeError)):
                    FiniteGroup(table)
        group = FiniteGroup([[0, 1], [1, 0]])
        with self.assertRaises(ValueError):
            group.mul(2, 0)


class PrimeFieldTests(unittest.TestCase):
    def test_field_arithmetic_five(self):
        field = PrimeField(5)
        self.assertEqual(field.add(4, 3), 2)
        self.assertEqual(field.sub(1, 4), 2)
        self.assertEqual(field.mul(4, 4), 1)
        self.assertEqual(field.pow(3, 4), 1)
        self.assertEqual(field.inv(3), 2)
        self.assertEqual(field.div(3, 2), 4)
        self.assertEqual(field.pow(3, -1), 2)

    def test_exhaustive_small_field_axioms(self):
        for prime in (2, 3, 5, 7, 11):
            f = PrimeField(prime)
            for a in range(prime):
                self.assertEqual(f.add(a, 0), a)
                self.assertEqual(f.mul(a, 1), a)
                if a:
                    self.assertEqual(f.mul(a, f.inv(a)), 1)
                for b in range(prime):
                    self.assertEqual(f.add(a, b), f.add(b, a))
                    self.assertEqual(f.mul(a, b), f.mul(b, a))
                    for c in range(prime):
                        self.assertEqual(f.mul(a, f.add(b, c)), f.add(f.mul(a, b), f.mul(a, c)))

    def test_reject_composite_and_unbounded(self):
        for n in (-1, 0, 1, 4, 9, 15, 21, 341, 561, 1105, 2**64):
            with self.subTest(n=n), self.assertRaises(ValueError):
                PrimeField(n)
        with self.assertRaises(TypeError):
            PrimeField(True)
        self.assertEqual(PrimeField(18446744073709551557).p, 18446744073709551557)

    def test_invalid_inputs_and_zero_inverse(self):
        field = PrimeField(7)
        with self.assertRaises(ZeroDivisionError):
            field.inv(0)
        with self.assertRaises(ZeroDivisionError):
            field.div(3, 7)
        with self.assertRaises(TypeError):
            field.add(1.2, 3)
        with self.assertRaises(TypeError):
            field.mul(True, 1)


class BooleanAndTropicalTests(unittest.TestCase):
    def test_boolean_identities_and_de_morgan(self):
        for a in (False, True):
            self.assertEqual(boolean_join(a, False), a)
            self.assertEqual(boolean_meet(a, True), a)
            self.assertFalse(boolean_meet(a, boolean_complement(a)))
            self.assertTrue(boolean_join(a, boolean_complement(a)))
            for b in (False, True):
                self.assertEqual(boolean_complement(boolean_join(a, b)),
                                 boolean_meet(boolean_complement(a), boolean_complement(b)))
                self.assertEqual(boolean_xor(a, b), a != b)

    def test_boolean_rejects_non_bools(self):
        with self.assertRaises(TypeError):
            boolean_join(0, True)
        with self.assertRaises(TypeError):
            boolean_complement(1)
        with self.assertRaises(TypeError):
            boolean_join(True, 1)
        with self.assertRaises(TypeError):
            boolean_meet(False, 1)

    def test_tropical_semiring_identities(self):
        infinity = math.inf
        for x in (-3.5, 0, 4, infinity):
            self.assertEqual(tropical_add(x, infinity), x)
            self.assertEqual(tropical_multiply(x, 0), x)
            self.assertEqual(tropical_multiply(x, infinity), infinity)
        self.assertEqual(tropical_add(3, -2), -2)
        self.assertEqual(tropical_multiply(3, -2), 1)

    def test_tropical_distributive(self):
        for a in (-3, 0, 4, math.inf):
            for b in (-2, 1, math.inf):
                for c in (-8, 2, math.inf):
                    left = tropical_multiply(a, tropical_add(b, c))
                    right = tropical_add(tropical_multiply(a, b), tropical_multiply(a, c))
                    self.assertEqual(left, right)

    def test_tropical_reject_invalid_floats(self):
        for bad in (float("nan"), -math.inf, True, 1j):
            with self.subTest(bad=bad), self.assertRaises((TypeError, ValueError)):
                tropical_add(0, bad)


if __name__ == "__main__":
    unittest.main()
