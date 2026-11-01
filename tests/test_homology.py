"""Exact, independent chain-complex and rational-homology witnesses."""

from fractions import Fraction
import random
import unittest

from calc.homology import ChainComplex, rational_rank


class RationalMatrixRankTests(unittest.TestCase):
    def test_rank_zero_rectangular_and_dependent(self):
        self.assertEqual(rational_rank(()), 0)
        self.assertEqual(rational_rank(((), (), ())), 0)
        self.assertEqual(rational_rank(((0, 0), (0, 0))), 0)
        self.assertEqual(rational_rank(((1, 2, 3), (2, 4, 6), (-1, -2, -3))), 1)

    def test_rank_full_and_rational(self):
        self.assertEqual(rational_rank(((0, 0, 1), (1, 0, 0), (0, 1, 0))), 3)
        self.assertEqual(rational_rank(((Fraction(1, 3), Fraction(1, 2)),
                                        (Fraction(2, 3), Fraction(2, 3)))), 2)

    def test_matrix_invalid(self):
        with self.assertRaises(ValueError):
            rational_rank(((1,), (1, 2)))
        with self.assertRaises(TypeError):
            rational_rank(((1.0,),))
        with self.assertRaises(TypeError):
            rational_rank(((True,),))
        with self.assertRaises(ValueError):
            rational_rank(((1 << 2050,),))
        with self.assertRaises(TypeError):
            rational_rank((row for row in ((1,),)))


class FiniteHomologyTests(unittest.TestCase):
    def test_oriented_interval(self):
        # One edge from vertex 0 to 1: d(e01)=v1-v0.
        c = ChainComplex((2, 1), {1: ((-1,), (1,))})
        self.assertEqual(c.boundary(1), ((Fraction(-1),), (Fraction(1),)))
        self.assertEqual(c.boundary_rank(1), 1)
        self.assertEqual(c.betti_numbers(), (1, 0))
        self.assertEqual(c.euler_characteristic(), 1)
        self.assertEqual(c.homology_euler_characteristic(), 1)

    def test_circle_from_one_vertex_one_loop(self):
        c = ChainComplex((1, 1), {1: ((0,),)})
        self.assertEqual(c.betti_numbers(), (1, 1))
        self.assertEqual(c.euler_characteristic(), 0)
        self.assertEqual(c.homology_euler_characteristic(), 0)

    def test_triangle_boundary_and_filling(self):
        # Columns e01,e12,e02; d(eij)=vj-vi.
        edge_boundary = ((-1, 0, -1), (1, -1, 0), (0, 1, 1))
        boundary_only = ChainComplex((3, 3), {1: edge_boundary})
        self.assertEqual(boundary_only.betti_numbers(), (1, 1))
        self.assertEqual(boundary_only.euler_characteristic(), 0)
        # Oriented 2-simplex [0,1,2] has boundary [1,2]-[0,2]+[0,1].
        filled = ChainComplex((3, 3, 1),
                              {1: edge_boundary, 2: ((1,), (1,), (-1,))})
        self.assertEqual(filled.betti_numbers(), (1, 0, 0))
        self.assertEqual(filled.boundary_rank(2), 1)
        self.assertEqual(filled.euler_characteristic(), 1)
        self.assertEqual(filled.homology_euler_characteristic(), 1)

    def test_sphere_cellular_chain_complex(self):
        # Cellular model C2=Q, C1=0, C0=Q: no boundary maps.
        c = ChainComplex((1, 0, 1))
        self.assertEqual(c.betti_numbers(), (1, 0, 1))
        self.assertEqual(c.euler_characteristic(), 2)
        self.assertEqual(c.homology_euler_characteristic(), 2)

    def test_two_disconnected_vertices_and_implicit_zero(self):
        c = ChainComplex((2,))
        self.assertEqual(c.betti_numbers(), (2,))
        d = ChainComplex((1, 1))
        self.assertEqual(d.boundary(1), ((Fraction(0),),))
        self.assertEqual(d.betti_numbers(), (1, 1))

    def test_rational_cancellation_chain_condition(self):
        # d1=(1/2, 1/2), d2=(1,-1)^T so d1*d2=0 exactly.
        c = ChainComplex((1, 2, 1), {
            1: ((Fraction(1, 2), Fraction(1, 2)),),
            2: ((1,), (-1,)),
        })
        self.assertEqual(c.betti_numbers(), (0, 0, 0))
        self.assertEqual(c.euler_characteristic(), 0)

    def test_seeded_exact_chain_complexes(self):
        # For d1=(a,b) and d2=(b,-a)^T, composition always vanishes.
        # Both maps have rank 1 for every nonzero rational pair.
        rng = random.Random(100271)
        for _ in range(50):
            a = Fraction(rng.randint(-8, 8), rng.randint(1, 7))
            b = Fraction(rng.randint(-8, 8), rng.randint(1, 7))
            if not a and not b:
                b = Fraction(1)
            c = ChainComplex((1, 2, 1), {
                1: ((a, b),),
                2: ((b,), (-a,)),
            })
            self.assertEqual(c.betti_numbers(), (0, 0, 0))
            self.assertEqual(c.euler_characteristic(), c.homology_euler_characteristic())

    def test_invalid_nonchain(self):
        with self.assertRaisesRegex(ValueError, "chain identity"):
            ChainComplex((1, 1, 1), {1: ((1,),), 2: ((1,),)})

    def test_invalid_dimensions_and_boundary_shape(self):
        for dims in ((), (33,), (1,) * 9, (-1,), (True,), (1.5,)):
            with self.subTest(dims=dims), self.assertRaises((TypeError, ValueError)):
                ChainComplex(dims)
        with self.assertRaises(ValueError):
            ChainComplex((2, 1), {1: ((1,),)})
        with self.assertRaises(ValueError):
            ChainComplex((1, 2), {1: ((1,),)})
        with self.assertRaises(ValueError):
            ChainComplex((1, 2), {0: ((1, 2),)})
        with self.assertRaises(TypeError):
            ChainComplex((1, 1), {True: ((1,),)})
        with self.assertRaises(TypeError):
            ChainComplex((dimension for dimension in (1, 1)))
        with self.assertRaises(TypeError):
            ChainComplex((1, 1), {1: ((value for value in (1,)),)})

    def test_invalid_coefficients(self):
        for value in (1.0, True, "1", 1j, 1 << 2050):
            with self.subTest(value=value), self.assertRaises((TypeError, ValueError)):
                ChainComplex((1, 1), {1: ((value,),)})

    def test_invalid_boundary_degrees(self):
        c = ChainComplex((1, 1))
        self.assertEqual(c.boundary_rank(0), 0)
        self.assertEqual(c.boundary_rank(2), 0)
        for degree in (-1, 3):
            with self.assertRaises(ValueError):
                c.boundary_rank(degree)
        with self.assertRaises(ValueError):
            c.boundary(0)
        with self.assertRaises(TypeError):
            c.boundary_rank(True)

    def test_immutable_snapshots(self):
        original = [[-1], [1]]
        c = ChainComplex((2, 1), {1: original})
        original[0][0] = 55
        self.assertEqual(c.betti_numbers(), (1, 0))
        self.assertEqual(c.boundary(1)[0][0], -1)
        with self.assertRaises((AttributeError, TypeError)):
            c.dimensions = (1,)


if __name__ == "__main__":
    unittest.main()
