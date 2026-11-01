"""Independent finite-ring witnesses and boundary tests for Z/nZ."""

from itertools import repeat
from math import gcd
import random
import unittest

from calc.rings import PrincipalIdeal, ResidueRing, chinese_remainder


class ResidueRingTests(unittest.TestCase):
    def test_exact_representative_operations(self):
        r = ResidueRing(12)
        self.assertEqual(r.normalize(-5), 7)
        self.assertEqual(r.add(11, 5), 4)
        self.assertEqual(r.subtract(2, 7), 7)
        self.assertEqual(r.neg(4), 8)
        self.assertEqual(r.multiply(-3, 7), 3)
        self.assertEqual(r.power(5, 0), 1)
        self.assertEqual(r.power(5, 2), 1)

    def test_exhaustive_ring_axioms_and_distributivity(self):
        for n in range(2, 14):
            ring = ResidueRing(n)
            self.assertEqual(ring.add(n - 1, 1), 0)
            for a in range(n):
                self.assertEqual(ring.add(a, 0), a)
                self.assertEqual(ring.multiply(a, 1), a)
                self.assertEqual(ring.add(a, ring.neg(a)), 0)
                for b in range(n):
                    self.assertEqual(ring.add(a, b), ring.add(b, a))
                    self.assertEqual(ring.multiply(a, b), ring.multiply(b, a))
                    for c in range(n):
                        self.assertEqual(ring.add(ring.add(a, b), c), ring.add(a, ring.add(b, c)))
                        self.assertEqual(ring.multiply(ring.multiply(a, b), c), ring.multiply(a, ring.multiply(b, c)))
                        self.assertEqual(ring.multiply(a, ring.add(b, c)),
                                         ring.add(ring.multiply(a, b), ring.multiply(a, c)))

    def test_units_inverses_and_nonzero_zero_divisors(self):
        for n in (4, 6, 8, 9, 12, 15, 17, 25, 30):
            r = ResidueRing(n)
            expected_units = tuple(a for a in range(n) if gcd(a, n) == 1)
            self.assertEqual(r.units(), expected_units)
            for a in range(n):
                self.assertEqual(r.is_unit(a), a in expected_units)
                expected_zero_divisor = a != 0 and any((a * b) % n == 0 for b in range(1, n))
                self.assertEqual(r.is_zero_divisor(a), expected_zero_divisor)
                if a in expected_units:
                    self.assertEqual(r.multiply(a, r.inverse(a)), 1)
                else:
                    with self.assertRaises(ZeroDivisionError):
                        r.inverse(a)

    def test_annihilator_matches_direct_multiplication(self):
        seeded = random.Random(0xAF15)
        for _ in range(80):
            n = seeded.randrange(2, 97)
            a = seeded.randrange(-n * 2, n * 2)
            ring = ResidueRing(n)
            annihilator = ring.annihilator(a)
            expected = tuple(x for x in range(n) if a * x % n == 0)
            self.assertEqual(annihilator.elements(), expected)
        self.assertEqual(ResidueRing(12).annihilator(0).order, 12)
        self.assertEqual(ResidueRing(12).annihilator(5).order, 1)

    def test_validate_all_inputs_and_resource_limits(self):
        for n in (-1, 0, 1, 4097, 2**32):
            with self.subTest(n=n), self.assertRaises(ValueError):
                ResidueRing(n)
        for invalid in (True, 3.0, "3", None):
            with self.subTest(invalid=invalid), self.assertRaises(TypeError):
                ResidueRing(invalid)
        r = ResidueRing(17)
        for invalid in (True, 3.0, "3", None):
            with self.subTest(invalid=invalid), self.assertRaises(TypeError):
                r.add(invalid, 1)
        for invalid in (-1, 1_000_001):
            with self.assertRaises(ValueError):
                r.power(3, invalid)
        with self.assertRaises(TypeError):
            r.power(3, False)
        with self.assertRaises(ValueError):
            r.ideal(*range(129))


class PrincipalIdealTests(unittest.TestCase):
    def test_independent_span_oracle_for_100_small_ideals(self):
        seeded = random.Random(0xEDA1)
        for _ in range(100):
            n = seeded.randrange(2, 120)
            generators = tuple(seeded.randrange(-2 * n, 2 * n) for _ in range(seeded.randrange(0, 5)))
            ring = ResidueRing(n)
            ideal = ring.ideal(*generators)
            # Independent additive-span construction, not the gcd formula.
            span = {0}
            for generator in generators:
                span = {(old + k * generator) % n for old in span for k in range(n)}
            self.assertEqual(set(ideal.elements()), span)
            self.assertEqual(ideal.order, len(span))
            self.assertEqual(ideal.quotient_order, n // len(span))
            for a in range(n):
                self.assertEqual(ideal.contains(a), a in span)

    def test_ideal_sum_intersection_and_product_as_actual_sets(self):
        r = ResidueRing(12)
        i, j = r.ideal(4), r.ideal(6)
        self.assertEqual(i.elements(), (0, 4, 8))
        self.assertEqual(j.elements(), (0, 6))
        a, b = set(i.elements()), set(j.elements())
        self.assertEqual(set(i.sum(j).elements()), {(x + y) % 12 for x in a for y in b})
        self.assertEqual(set(i.intersection(j).elements()), a & b)
        self.assertEqual(set(i.product(j).elements()), {(x * y) % 12 for x in a for y in b})

    def test_quotient_residue_classes_and_homomorphism(self):
        seeded = random.Random(0xCA55)
        for _ in range(70):
            n = seeded.randrange(2, 101)
            i = ResidueRing(n).ideal(seeded.randrange(-n, n))
            a, b = seeded.randrange(-999, 999), seeded.randrange(-999, 999)
            self.assertEqual(i.quotient_class(a + b), i.quotient_add(a, b))
            self.assertEqual(i.quotient_class(a * b), i.quotient_multiply(a, b))
            for element in i.elements():
                self.assertEqual(i.quotient_class(a + element), i.quotient_class(a))
        whole = ResidueRing(12).ideal(5)
        self.assertEqual(whole.quotient_order, 1)
        self.assertEqual(whole.quotient_class(999), 0)
        zero_ideal = ResidueRing(12).ideal()
        self.assertEqual(zero_ideal.elements(), (0,))
        self.assertEqual(zero_ideal.quotient_order, 12)

    def test_invalid_ideal_divisor_and_cross_ring_ops(self):
        r12, r10 = ResidueRing(12), ResidueRing(10)
        for bad in (-1, 0, 5, 13):
            with self.assertRaises(ValueError):
                PrincipalIdeal(r12, bad)
        with self.assertRaises(TypeError):
            PrincipalIdeal(r12, True)
        with self.assertRaises(TypeError):
            PrincipalIdeal("12", 2)
        for operation in ("sum", "intersection", "product"):
            with self.assertRaises(ValueError):
                getattr(r12.ideal(2), operation)(r10.ideal(2))


class ChineseRemainderTests(unittest.TestCase):
    def test_classic_three_moduli(self):
        self.assertEqual(chinese_remainder((2, 3, 2), (3, 5, 7)), (23, 105))
        self.assertEqual(chinese_remainder((-1, 9), (4, 5)), (19, 20))
        self.assertEqual(chinese_remainder((123,), (17,)), (4, 17))

    def test_exhaustive_coprime_pair_congruences(self):
        for m in range(2, 11):
            for n in range(2, 11):
                if gcd(m, n) != 1:
                    continue
                for a in range(m):
                    for b in range(n):
                        result, product = chinese_remainder((a, b), (m, n))
                        self.assertEqual(product, m * n)
                        self.assertTrue(0 <= result < product)
                        self.assertEqual(result % m, a)
                        self.assertEqual(result % n, b)

    def test_resource_and_invalid_rejections(self):
        for bad_residues, bad_moduli in (
            ((), ()), ((1, 2), (3,)), ((1, 2), (4, 6)),
            ((True,), (5,)), ((1,), (1,)), ((1,), (5000,)),
        ):
            with self.subTest(case=(bad_residues, bad_moduli)), self.assertRaises((ValueError, TypeError)):
                chinese_remainder(bad_residues, bad_moduli)
        with self.assertRaises(ValueError):
            chinese_remainder(repeat(1), repeat(2))
        with self.assertRaises(ValueError):
            chinese_remainder((0,) * 10, (101, 103, 107, 109, 113, 127, 131, 137, 139, 149))


if __name__ == "__main__":
    unittest.main()
