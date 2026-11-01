"""Independent identities for finite free associative rational algebras."""

from __future__ import annotations

from fractions import Fraction
from itertools import islice, product
import random
import unittest

from calc.free_algebra import FreeAlgebra, FreePolynomial, ResourceLimitError


class FreeAssociativeAlgebraTests(unittest.TestCase):
    def setUp(self) -> None:
        self.a = FreeAlgebra(("x", "y", "z"))
        self.x = self.a.variable("x")
        self.y = self.a.variable("y")
        self.z = self.a.variable("z")

    def test_empty_word_is_unit_and_zero(self):
        self.assertEqual(self.a.one().terms, (((), Fraction(1)),))
        self.assertFalse(self.a.zero())
        for p in (self.x, self.x*self.y, self.x+self.y*Fraction(2, 3)):
            self.assertEqual(p * self.a.one(), p)
            self.assertEqual(self.a.one() * p, p)
            self.assertEqual(p * self.a.zero(), self.a.zero())

    def test_word_order_is_mathematically_distinct(self):
        xy, yx = self.x * self.y, self.y * self.x
        self.assertNotEqual(xy, yx)
        self.assertEqual(xy.terms[0][0], ("x", "y"))
        self.assertEqual(yx.terms[0][0], ("y", "x"))
        self.assertEqual((self.x + self.y)**2,
                         self.x*self.x + xy + yx + self.y*self.y)

    def test_exact_rationals_and_cancellation(self):
        p = Fraction(2, 7)*self.x*self.y + Fraction(3, 5)*self.y
        self.assertEqual((p - p).terms, ())
        self.assertEqual((p + Fraction(1, 7)*self.x*self.y).as_dict()[("x", "y")],
                         Fraction(3, 7))
        self.assertEqual((self.x + Fraction(1, 2))*(self.x - Fraction(1, 2)),
                         self.x*self.x - Fraction(1, 4))

    def test_seeded_associativity_and_distributivity(self):
        rng = random.Random(0xBEEF)
        words = [(), ("x",), ("y",), ("z",), ("x", "y"), ("y", "x"), ("z", "z")]
        for _ in range(80):
            samples = []
            for _ in range(3):
                coefficients = {word: Fraction(rng.randrange(-3,4),rng.randrange(1,5))
                                for word in rng.sample(words, 4)}
                samples.append(FreePolynomial(self.a, coefficients))
            p, q, r = samples
            self.assertEqual((p*q)*r, p*(q*r))
            self.assertEqual(p*(q+r), p*q+p*r)
            self.assertEqual((p+q)*r, p*r+q*r)

    def test_integer_powers(self):
        self.assertEqual(self.x**0, self.a.one())
        self.assertEqual(self.x**1, self.x)
        self.assertEqual(self.x**5, self.a.monomial(("x",)*5))
        self.assertEqual(self.a.zero()**0, self.a.one())

    def test_free_partial_derivative_each_occurrence(self):
        xyx = self.x*self.y*self.x
        self.assertEqual(xyx.derivative("x"), self.y*self.x+self.x*self.y)
        self.assertEqual(xyx.derivative("y"), self.x*self.x)
        self.assertEqual(self.a.constant(7).derivative("x"), self.a.zero())
        self.assertEqual((self.x**4).derivative("x"), 4*self.x**3)

    def test_free_partial_derivative_leibniz_seeded(self):
        rng = random.Random(19)
        words = [(), ("x",), ("y",), ("x", "x"), ("y", "x"), ("x", "y")]
        for _ in range(70):
            p = FreePolynomial(self.a,{word:rng.randint(-4,4) for word in rng.sample(words,3)})
            q = FreePolynomial(self.a,{word:rng.randint(-4,4) for word in rng.sample(words,3)})
            self.assertEqual((p*q).derivative("x"),p.derivative("x")*q+p*q.derivative("x"))

    def test_arbitrary_generator_images_extend_to_derivation(self):
        images = {"x":self.x*self.y, "y":self.x}
        self.assertEqual((self.x*self.y).derivation(images),
                         self.x*self.y*self.y+self.x*self.x)
        p = self.x+self.y*self.z
        q = self.z*self.x-self.y
        self.assertEqual((p*q).derivation(images),
                         p.derivation(images)*q+p*q.derivation(images))
        self.assertEqual(self.z.derivation(images),self.a.zero())

    def test_commutator_antisymmetric_and_jacobi(self):
        p,q,r=self.x,self.y,self.z
        self.assertEqual(p.commutator(q),-q.commutator(p))
        self.assertEqual(p.commutator(p),self.a.zero())
        jacobi=(p.commutator(q)).commutator(r)+(q.commutator(r)).commutator(p)+(r.commutator(p)).commutator(q)
        self.assertEqual(jacobi,self.a.zero())

    def test_commutator_derivation_identity(self):
        lhs=self.x.commutator(self.y*self.z)
        rhs=self.x.commutator(self.y)*self.z+self.y*self.x.commutator(self.z)
        self.assertEqual(lhs,rhs)

    def test_formal_algebra_identity_not_matrix_accident(self):
        self.assertNotEqual(self.x*self.y*self.x, self.x*self.x*self.y)
        self.assertEqual((self.x+self.y)*(self.x-self.y),
                         self.x*self.x-self.x*self.y+self.y*self.x-self.y*self.y)

    def test_immutable_canonical_sparse_representation(self):
        terms={("y",):1,("x",):0,("x", "z"):Fraction(-2,3)}
        p=FreePolynomial(self.a,terms)
        terms[("y",)]=8
        self.assertEqual(p.terms,((('x','z'),Fraction(-2,3)),(('y',),Fraction(1))))
        d=p.as_dict()
        d[("y",)]=0
        self.assertEqual(p.as_dict()[("y",)],1)
        with self.assertRaises(AttributeError):
            p.terms=()

    def test_bad_generator_metadata(self):
        for names in ((), ("x","x"), ("a b",), ("é",), ("A"*33,), tuple("abcdefghij")):
            with self.subTest(names=names),self.assertRaises(ValueError):
                FreeAlgebra(names)
        with self.assertRaises(ValueError):
            FreeAlgebra(["x","y"])

    def test_bad_terms_or_types(self):
        with self.assertRaises(TypeError):
            self.a.constant(True)
        with self.assertRaises(TypeError):
            self.a.constant(0.5)
        with self.assertRaises(TypeError):
            FreePolynomial(self.a,[(("x",),1)])
        with self.assertRaises(TypeError):
            FreePolynomial(self.a,{"x":1})
        with self.assertRaises(ValueError):
            FreePolynomial(self.a,{("unknown",):1})
        with self.assertRaises(ValueError):
            self.x.derivative("bad")
        with self.assertRaises(TypeError):
            self.x.derivation({"x":1})
        with self.assertRaises(ValueError):
            self.x.derivation({"bad":self.y})

    def test_cross_algebra_is_rejected(self):
        other=FreeAlgebra(("y","x","z"))
        with self.assertRaises(ValueError):
            self.x+other.variable("x")
        with self.assertRaises(ValueError):
            self.x.derivation({"x":other.variable("x")})

    def test_resource_limits_word_coeff_degree_and_term_count(self):
        with self.assertRaises(ResourceLimitError):
            self.a.monomial(("x",)*65)
        with self.assertRaises(ResourceLimitError):
            self.a.constant(1<<4100)
        with self.assertRaises(ValueError):
            self.x**33
        with self.assertRaises(ResourceLimitError):
            self.x**32*self.x**32*self.x
        too_many={word:1 for word in islice(product(("x","y","z"),repeat=9),4097)}
        with self.assertRaises(ResourceLimitError):
            FreePolynomial(self.a,too_many)

    def test_product_pair_limit_and_derivation_limit(self):
        words=list(islice(product(("x","y","z"),repeat=6),320))
        p=FreePolynomial(self.a,{word:1 for word in words})
        with self.assertRaises(ResourceLimitError):
            p*p
        # Each term in a long-word polynomial differentiates into many
        # image terms; reject on operation count before an unbounded loop.
        many=FreePolynomial(self.a,{word:1 for word in words})
        image=FreePolynomial(self.a,{word:1 for word in words})
        with self.assertRaises(ResourceLimitError):
            many.derivation({"x":image,"y":image,"z":image})


if __name__=="__main__":
    unittest.main()
