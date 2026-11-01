"""Independent finite-category axioms and counterexamples."""

import unittest
from itertools import repeat

from calc.categories import (
    Arrow, FiniteCategory, FiniteFunctor, NaturalTransformation,
    monoid_category, poset_category,
)


def chain(*objects):
    return poset_category(objects, [(a, b) for i, a in enumerate(objects)
                                    for b in objects[i:]])


def functor_from_object_map(source, target, obj_map):
    return FiniteFunctor(
        source, target, obj_map,
        {arrow.name: target.hom(obj_map[arrow.source], obj_map[arrow.target])[0]
         for arrow in source.arrows},
    )


def constants_monoid():
    """Actual 2-point transformation monoid under function composition."""
    labels = ("id", "c0", "c1")
    return monoid_category(labels, "id", {(g, f): (f if g == "id" else g)
                                           for g in labels for f in labels})


class CategoryAxioms(unittest.TestCase):
    def test_terminal_singleton(self):
        cat = FiniteCategory(("x",), (Arrow("1", "x", "x"),),
                             {"x": "1"}, {("1", "1"): "1"})
        self.assertEqual(cat.compose("1", "1"), "1")
        self.assertEqual(cat.hom("x", "x"), ("1",))

    def test_poset_thin_category_and_all_composable_triples(self):
        cat = chain("a", "b", "c", "d")
        self.assertEqual(len(cat.arrows), 10)
        for f in cat.arrows:
            self.assertEqual(cat.compose(cat.identity(f.target), f.name), f.name)
            self.assertEqual(cat.compose(f.name, cat.identity(f.source)), f.name)
        for f in cat.arrows:
            for g in cat.arrows:
                if f.target != g.source:
                    continue
                for h in cat.arrows:
                    if g.target == h.source:
                        self.assertEqual(cat.compose(h.name, cat.compose(g.name, f.name)),
                                         cat.compose(cat.compose(h.name, g.name), f.name))
        self.assertEqual(cat.hom("d", "a"), ())
        self.assertRaises(ValueError, cat.compose, cat.hom("a", "b")[0],
                          cat.hom("b", "c")[0])

    def test_monoid_category_associativity_from_function_composition(self):
        cat = constants_monoid()
        self.assertEqual(cat.compose("c0", "c1"), "c0")
        self.assertEqual(cat.compose("c1", "c0"), "c1")
        self.assertEqual(cat.compose("id", "c1"), "c1")

    def test_immutability_of_input_tables(self):
        morphs = {("1", "1"): "1"}
        cat = FiniteCategory(("x",), (Arrow("1", "x", "x"),), {"x": "1"}, morphs)
        morphs[("1", "1")] = "wrong"
        self.assertEqual(cat.compose("1", "1"), "1")
        with self.assertRaises(TypeError):
            cat.composition[("1", "1")] = "wrong"

    def test_invalid_categories(self):
        with self.assertRaisesRegex(ValueError, "unique"):
            FiniteCategory(("x", "x"), (Arrow("1", "x", "x"),),
                           {"x": "1"}, {("1", "1"): "1"})
        with self.assertRaisesRegex(ValueError, "unknown endpoint"):
            FiniteCategory(("x",), (Arrow("f", "y", "x"),),
                           {"x": "f"}, {("f", "f"): "f"})
        with self.assertRaisesRegex(ValueError, "incomplete"):
            FiniteCategory(("x",), (Arrow("1", "x", "x"),), {"x": "1"}, {})
        with self.assertRaisesRegex(ValueError, "wrong result type"):
            FiniteCategory(("a", "b"),
                           (Arrow("ia", "a", "a"), Arrow("ib", "b", "b"),
                            Arrow("f", "a", "b")),
                           {"a": "ia", "b": "ib"},
                           {("ia", "ia"): "ia", ("ib", "ib"): "ib",
                            ("f", "ia"): "f", ("ib", "f"): "ia"})

    def test_invalid_identity_law(self):
        with self.assertRaisesRegex(ValueError, "identity law"):
            FiniteCategory(("x",),
                           (Arrow("e", "x", "x"), Arrow("a", "x", "x")),
                           {"x": "e"},
                           {("e", "e"): "e", ("e", "a"): "e",
                            ("a", "e"): "a", ("a", "a"): "a"})

    def test_nonassociative_unital_magma_is_rejected(self):
        elems = ("e", "a", "b")
        table = {(g, f): (f if g == "e" else g if f == "e" else
                          {("a", "a"): "a", ("a", "b"): "a",
                           ("b", "a"): "b", ("b", "b"): "a"}[(g, f)])
                 for g in elems for f in elems}
        with self.assertRaisesRegex(ValueError, "associativity"):
            monoid_category(elems, "e", table)

    def test_invalid_poset_relations(self):
        with self.assertRaisesRegex(ValueError, "reflexive"):
            poset_category(("a", "b"), {("a", "a")})
        with self.assertRaisesRegex(ValueError, "antisymmetric"):
            poset_category(("a", "b"), {("a", "a"), ("b", "b"),
                                             ("a", "b"), ("b", "a")})
        with self.assertRaisesRegex(ValueError, "transitive"):
            poset_category(("a", "b", "c"), {("a", "a"), ("b", "b"),
                                                   ("c", "c"), ("a", "b"),
                                                   ("b", "c")})

    def test_limits_and_invalid_labels(self):
        with self.assertRaisesRegex(ValueError, "1..32"):
            poset_category(tuple(str(i) for i in range(33)), ())
        with self.assertRaisesRegex(ValueError, "too many morphisms"):
            chain(*tuple(str(i) for i in range(12)))
        with self.assertRaises(ValueError):
            Arrow(" badly ", "x", "y")
        with self.assertRaises(TypeError):
            FiniteCategory("x", (), {}, {})

    def test_infinite_iterables_are_bounded_before_exhaustion(self):
        with self.assertRaisesRegex(ValueError, "1..32"):
            poset_category(repeat("x"), ())
        with self.assertRaisesRegex(ValueError, "too many morphisms"):
            poset_category(("x",), repeat(("x", "x")))
        with self.assertRaisesRegex(ValueError, "exceeds 64"):
            FiniteCategory(("x",), repeat(Arrow("e", "x", "x")),
                           {"x": "e"}, {("e", "e"): "e"})


class FunctorAxioms(unittest.TestCase):
    def test_monotone_poset_maps_are_functors(self):
        src, dst = chain("a", "b", "c"), chain("x", "y")
        functor = functor_from_object_map(src, dst, {"a": "x", "b": "x", "c": "y"})
        self.assertEqual(functor.on_objects["c"], "y")
        for (g, f), gf in src.composition.items():
            self.assertEqual(functor.on_arrows[gf],
                             dst.compose(functor.on_arrows[g], functor.on_arrows[f]))

    def test_bad_endpoint_or_identity_or_composition(self):
        src, dst = chain("a", "b"), chain("x", "y")
        good = functor_from_object_map(src, dst, {"a": "x", "b": "y"})
        broken = dict(good.on_arrows)
        broken[src.hom("a", "b")[0]] = dst.identity("x")
        with self.assertRaisesRegex(ValueError, "endpoints"):
            FiniteFunctor(src, dst, dict(good.on_objects), broken)
        src_monoid = monoid_category(("e", "a"), "e",
                                    {("e", "e"): "e", ("e", "a"): "a",
                                     ("a", "e"): "a", ("a", "a"): "a"})
        group2 = monoid_category(("1", "s"), "1",
                                 {("1", "1"): "1", ("1", "s"): "s",
                                  ("s", "1"): "s", ("s", "s"): "1"})
        with self.assertRaisesRegex(ValueError, "composition"):
            FiniteFunctor(src_monoid, group2, {"•": "•"}, {"e": "1", "a": "s"})
        with self.assertRaisesRegex(ValueError, "identity"):
            FiniteFunctor(src_monoid, group2, {"•": "•"}, {"e": "s", "a": "1"})

    def test_functor_mapping_must_be_total(self):
        src = chain("a", "b")
        dst = chain("x", "y")
        with self.assertRaisesRegex(ValueError, "total"):
            FiniteFunctor(src, dst, {"a": "x"}, {})


class NaturalityAxioms(unittest.TestCase):
    def test_pointwise_order_transformation(self):
        src, dst = chain("a", "b", "c"), chain("x", "y")
        f = functor_from_object_map(src, dst, {"a": "x", "b": "x", "c": "y"})
        g = functor_from_object_map(src, dst, {"a": "x", "b": "y", "c": "y"})
        eta = NaturalTransformation(f, g, {"a": dst.identity("x"),
                                           "b": dst.hom("x", "y")[0],
                                           "c": dst.identity("y")})
        for arrow in src.arrows:
            self.assertEqual(
                dst.compose(g.on_arrows[arrow.name], eta.components[arrow.source]),
                dst.compose(eta.components[arrow.target], f.on_arrows[arrow.name]))

    def test_noncommuting_naturality_square_is_rejected(self):
        src = monoid_category(("e", "a"), "e",
                              {("e", "e"): "e", ("e", "a"): "a",
                               ("a", "e"): "a", ("a", "a"): "a"})
        target = constants_monoid()
        f = FiniteFunctor(src, target, {"•": "•"}, {"e": "id", "a": "c0"})
        g = FiniteFunctor(src, target, {"•": "•"}, {"e": "id", "a": "c1"})
        with self.assertRaisesRegex(ValueError, "naturality"):
            NaturalTransformation(f, g, {"•": "id"})
        eta = NaturalTransformation(f, g, {"•": "c1"})
        self.assertEqual(eta.components["•"], "c1")

    def test_natural_component_endpoint_and_totality(self):
        source = chain("a", "b")
        target = chain("x", "y")
        f = functor_from_object_map(source, target, {"a": "x", "b": "y"})
        with self.assertRaisesRegex(ValueError, "exactly one"):
            NaturalTransformation(f, f, {"a": target.identity("x")})
        with self.assertRaisesRegex(ValueError, "endpoints"):
            NaturalTransformation(f, f, {"a": target.identity("y"),
                                         "b": target.identity("y")})


if __name__ == "__main__":
    unittest.main()
