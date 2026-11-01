"""Independent finite-module axioms and submodule/quotient oracles."""

from itertools import count, product
import random
import unittest

from calc.modules import FiniteSubmodule, QuotientModule, ResidueFreeModule


class FiniteModuleTests(unittest.TestCase):
    def test_canonical_arithmetic(self):
        m = ResidueFreeModule(6, 2)
        self.assertEqual(m.order, 36)
        self.assertEqual(m.zero, (0, 0))
        self.assertEqual(m.normalize((-1, 13)), (5, 1))
        self.assertEqual(m.add((5, 1), (4, 8)), (3, 3))
        self.assertEqual(m.neg((5, 1)), (1, 5))
        self.assertEqual(m.subtract((2, 2), (5, 1)), (3, 1))
        self.assertEqual(m.scalar_multiply(-2, (2, 4)), (2, 4))
        self.assertEqual(m.scalar_multiply(0, (4, 2)), (0, 0))

    def test_module_axioms_exhaustive_small_spaces(self):
        for modulus, rank in ((2, 2), (3, 2), (4, 1)):
            module = ResidueFreeModule(modulus, rank)
            elements = tuple(product(range(modulus), repeat=rank))
            for x in elements:
                self.assertEqual(module.add(x, module.zero), x)
                self.assertEqual(module.add(x, module.neg(x)), module.zero)
                self.assertEqual(module.scalar_multiply(1, x), x)
                for scalar in range(modulus):
                    for other_scalar in range(modulus):
                        self.assertEqual(module.scalar_multiply(scalar + other_scalar, x),
                                         module.add(module.scalar_multiply(scalar, x),
                                                    module.scalar_multiply(other_scalar, x)))
                        self.assertEqual(module.scalar_multiply(scalar * other_scalar, x),
                                         module.scalar_multiply(scalar, module.scalar_multiply(other_scalar, x)))
                for y in elements:
                    self.assertEqual(module.add(x, y), module.add(y, x))
                    for scalar in range(modulus):
                        self.assertEqual(module.scalar_multiply(scalar, module.add(x, y)),
                                         module.add(module.scalar_multiply(scalar, x),
                                                    module.scalar_multiply(scalar, y)))
                    for z in elements:
                        self.assertEqual(module.add(module.add(x, y), z),
                                         module.add(x, module.add(y, z)))

    def test_exact_independent_generated_span(self):
        m = ResidueFreeModule(6, 2)
        n = m.generated_submodule((2, 0), (0, 3))
        expected = {(x, y) for x in (0, 2, 4) for y in (0, 3)}
        self.assertEqual(set(n.elements), expected)
        self.assertEqual(n.order, 6)
        self.assertEqual(n.index, 6)
        self.assertTrue(n.contains((-2, -3)))
        self.assertFalse(n.contains((1, 0)))

    def test_seeded_generator_span_matches_enumerated_coefficients(self):
        rng = random.Random(0xC0DE5)
        for _ in range(100):
            modulus = rng.randrange(2, 12)
            rank = rng.randrange(1, 4)
            module = ResidueFreeModule(modulus, rank)
            generators = tuple(tuple(rng.randrange(-15, 16) for _ in range(rank))
                               for _ in range(rng.randrange(0, 3)))
            submodule = module.generated_submodule(*generators)
            independent = {
                tuple(sum(a * g[i] for a, g in zip(coefficients, generators)) % modulus
                      for i in range(rank))
                for coefficients in product(range(modulus), repeat=len(generators))
            }
            self.assertEqual(set(submodule.elements), independent)
            self.assertEqual(submodule.order * submodule.index, module.order)

    def test_zero_and_whole_submodules(self):
        m = ResidueFreeModule(5, 2)
        zero = m.generated_submodule()
        self.assertEqual(zero.elements, ((0, 0),))
        self.assertEqual(zero.index, 25)
        whole = m.generated_submodule((1, 0), (0, 1))
        self.assertEqual(whole.order, m.order)
        self.assertEqual(whole.index, 1)
        self.assertEqual(whole.quotient().order, 1)
        self.assertEqual(zero.quotient().order, 25)

    def test_coset_partition_and_operation_well_defined(self):
        m = ResidueFreeModule(4, 2)
        sub = m.generated_submodule((2, 0), (0, 2))
        quotient = sub.quotient()
        self.assertEqual(quotient.order, 4)
        members = [member for coset in quotient.cosets for member in coset]
        self.assertEqual(sorted(members), list(m.elements()))
        self.assertEqual(len(members), len(set(members)))
        self.assertEqual(quotient.zero_class, quotient.class_of(m.zero))
        for x, y in product(m.elements(), repeat=2):
            i, j = quotient.class_of(x), quotient.class_of(y)
            self.assertEqual(quotient.add(i, j), quotient.class_of(m.add(x, y)))
            for scalar in range(m.modulus):
                self.assertEqual(quotient.scalar_multiply(scalar, i),
                                 quotient.class_of(m.scalar_multiply(scalar, x)))
            self.assertEqual(quotient.neg(i), quotient.class_of(m.neg(x)))

    def test_quotient_class_independent_of_representative(self):
        m = ResidueFreeModule(6, 2)
        sub = m.generated_submodule((2, 0), (0, 3))
        q = sub.quotient()
        for x in m.elements():
            for shift in sub.elements:
                self.assertEqual(q.class_of(x), q.class_of(m.add(x, shift)))
        self.assertEqual(q.order * sub.order, m.order)

    def test_zero_module_size_one_coset_group(self):
        m = ResidueFreeModule(2, 8)
        sub = m.generated_submodule((1, 0, 0, 0, 0, 0, 0, 0))
        self.assertEqual(sub.order, 2)
        self.assertEqual(sub.index, 128)
        self.assertEqual(sub.quotient().order, 128)

    def test_invalid_ambient_moduli_and_dimensions(self):
        for args in ((1, 1), (0, 1), (4097, 1), (3, 0), (2, 9), (7, 5)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                ResidueFreeModule(*args)
        for args in ((True, 1), (3.0, 1), (3, False), (3, "2")):
            with self.subTest(args=args), self.assertRaises(TypeError):
                ResidueFreeModule(*args)

    def test_invalid_vectors_scalars_and_bounded_generator_reads(self):
        m = ResidueFreeModule(9, 2)
        for value in ((1,), (1, 2, 3)):
            with self.assertRaises(ValueError):
                m.normalize(value)
        for value in ((True, 2), (1.0, 2), ("a", 2)):
            with self.assertRaises(TypeError):
                m.normalize(value)
        for value in (3.0, True, "3"):
            with self.assertRaises(TypeError):
                m.scalar_multiply(value, (1, 2))
        with self.assertRaises(TypeError):
            m.normalize(None)
        with self.assertRaises(TypeError):
            m.normalize("12")
        with self.assertRaises(ValueError):
            m.normalize(count())  # finite consumption of rank+1 values
        with self.assertRaises(ValueError):
            m.generated_submodule(*([(1, 2)] * 17))

    def test_quotient_rejects_implicit_constructors_and_bad_indices(self):
        module = ResidueFreeModule(4, 1)
        sub = module.generated_submodule((2,))
        with self.assertRaises(TypeError):
            FiniteSubmodule(module, ((0,),))
        with self.assertRaises(TypeError):
            QuotientModule(sub, (((0,),),))
        q = sub.quotient()
        for bad in (-1, q.order):
            with self.assertRaises(ValueError):
                q.add(0, bad)
        with self.assertRaises(TypeError):
            q.neg(True)


if __name__ == "__main__":
    unittest.main()
