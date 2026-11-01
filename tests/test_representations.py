"""Independent exact checks of explicitly enumerated group representations."""

import unittest
from fractions import Fraction
from itertools import permutations, repeat

from calc.representations import (
    FiniteGroup, RationalRepresentation, character_pairing, cyclic_group,
    direct_sum, intertwines, regular_representation, trivial_representation,
)


def sign_c2(group):
    return RationalRepresentation(group, {"e": ((1,),), "g": ((-1,),)})


class GroupLawTests(unittest.TestCase):
    def test_trivial_and_cyclic_group_tables(self):
        singleton = cyclic_group(1)
        self.assertEqual(singleton.multiply("e", "e"), "e")
        for n in range(2, 12):
            group = cyclic_group(n)
            self.assertEqual(group.elements[0], "e")
            for i, a in enumerate(group.elements):
                for j, b in enumerate(group.elements):
                    self.assertEqual(group.multiply(a, b), group.elements[(i + j) % n])
                self.assertEqual(group.multiply(a, group.inverse(a)), "e")

    def test_group_rejects_missing_pairs_and_nonclosure(self):
        with self.assertRaises(ValueError):
            FiniteGroup(("e", "a"), {("e", "e"): "e"}, "e")
        table = {(a, b): "e" for a in ("e", "a") for b in ("e", "a")}
        with self.assertRaises(ValueError):
            FiniteGroup(("e", "a"), table, "e")

    def test_group_rejects_nonassociative_full_table_with_inverses(self):
        names = ("e", "a", "b")
        table = {(a, b): (b if a == "e" else a if b == "e" else "e")
                 for a in names for b in names}
        with self.assertRaisesRegex(ValueError, "associative"):
            FiniteGroup(names, table, "e")

    def test_invalid_group_labels_and_bounds(self):
        for bad in (0, True, 17, -1):
            with self.assertRaises(ValueError):
                cyclic_group(bad)
        with self.assertRaises(ValueError):
            FiniteGroup(("e", "e"), {("e", "e"): "e"}, "e")
        with self.assertRaises(TypeError):
            FiniteGroup("e", {}, "e")
        with self.assertRaises(ValueError):
            FiniteGroup(repeat("e"), {}, "e")

    def test_unknown_element_does_not_silently_coerce(self):
        group = cyclic_group(2)
        with self.assertRaises(KeyError):
            group.multiply("g", "missing")
        with self.assertRaises(KeyError):
            group.inverse("missing")


class RepresentationTests(unittest.TestCase):
    def test_regular_c2_is_exact_swap(self):
        group = cyclic_group(2)
        regular = regular_representation(group)
        self.assertEqual(regular.matrix("e"), ((1, 0), (0, 1)))
        self.assertEqual(regular.matrix("g"), ((0, 1), (1, 0)))
        self.assertEqual(regular.dimension, 2)

    def test_c3_regular_is_cyclic_permutation(self):
        group = cyclic_group(3)
        regular = regular_representation(group)
        self.assertEqual(regular.matrix("g"), ((0, 0, 1), (1, 0, 0), (0, 1, 0)))
        self.assertEqual(regular.matrix("g2"), ((0, 1, 0), (0, 0, 1), (1, 0, 0)))

    def test_nonabelian_s3_regular_representation(self):
        # Independent permutation composition generates a noncommutative
        # group rather than relying on a hard-coded Cayley table.
        maps = tuple(permutations(range(3)))
        labels = tuple("".join(map(str, perm)) for perm in maps)
        table = {(labels[i], labels[j]): "".join(str(a[b[k]]) for k in range(3))
                 for i, a in enumerate(maps) for j, b in enumerate(maps)}
        g = FiniteGroup(labels, table, "012")
        self.assertNotEqual(g.multiply("102", "021"), g.multiply("021", "102"))
        regular = regular_representation(g)
        self.assertEqual(regular.dimension, 6)
        self.assertEqual(regular.character()["012"], 6)
        self.assertTrue(all(regular.character()[x] == 0 for x in labels if x != "012"))
        self.assertEqual(character_pairing(regular, trivial_representation(g)), 1)

    def test_regular_character_is_group_order_at_identity_only(self):
        for n in (1, 2, 3, 4, 8):
            g = cyclic_group(n)
            chars = regular_representation(g).character()
            self.assertEqual(chars["e"], n)
            self.assertTrue(all(chars[x] == 0 for x in g.elements[1:]))

    def test_character_pairing_closed_form(self):
        for n in (2, 3, 5):
            g = cyclic_group(n)
            trivial = trivial_representation(g)
            regular = regular_representation(g)
            self.assertEqual(character_pairing(trivial, trivial), 1)
            self.assertEqual(character_pairing(regular, trivial), 1)
            self.assertEqual(character_pairing(regular, regular), n)

    def test_c2_rational_irreducible_character_orthogonality(self):
        g = cyclic_group(2)
        triv, sign = trivial_representation(g), sign_c2(g)
        self.assertEqual(character_pairing(triv, sign), 0)
        self.assertEqual(character_pairing(sign, sign), 1)
        self.assertEqual(character_pairing(regular_representation(g), sign), 1)

    def test_c2_regular_decomposes_via_explicit_intertwiner(self):
        g = cyclic_group(2)
        decomposed = direct_sum(trivial_representation(g), sign_c2(g))
        regular = regular_representation(g)
        self.assertEqual(decomposed.matrix("g"), ((1, 0), (0, -1)))
        # Columns (1,1) and (1,-1) are eigenvectors of the swap action.
        basis_change = ((1, 1), (1, -1))
        self.assertTrue(intertwines(decomposed, regular, basis_change))
        self.assertFalse(intertwines(decomposed, regular, ((1, 0), (0, 1))))
        self.assertEqual(decomposed.character(), regular.character())

    def test_reject_nonhomomorphic_assignment(self):
        g = cyclic_group(2)
        with self.assertRaisesRegex(ValueError, "homomorphism"):
            RationalRepresentation(g, {"e": ((1,),), "g": ((2,),)})
        with self.assertRaisesRegex(ValueError, "identity"):
            RationalRepresentation(g, {"e": ((-1,),), "g": ((1,),)})

    def test_reject_inexact_or_bad_dimensions(self):
        g = cyclic_group(2)
        with self.assertRaises(TypeError):
            RationalRepresentation(g, {"e": ((1,),), "g": ((1.0,),)})
        with self.assertRaises(ValueError):
            RationalRepresentation(g, {"e": ((1, 0), (0, 1)), "g": ((1,),)})
        with self.assertRaises(ValueError):
            RationalRepresentation(g, {"e": ((1,),)})

    def test_rational_coefficient_storage_is_immutable_and_exact(self):
        g = cyclic_group(2)
        x = [[1]]
        y = [[-1]]
        rep = RationalRepresentation(g, {"e": x, "g": y})
        x[0][0] = 5
        y[0][0] = 5
        self.assertEqual(rep.matrix("g"), ((Fraction(-1),),))
        with self.assertRaises(TypeError):
            rep.matrix("g")[0][0] = 4

    def test_invalid_intertwiner_shape_and_group(self):
        a = trivial_representation(cyclic_group(2))
        b = trivial_representation(cyclic_group(3))
        with self.assertRaises(ValueError):
            intertwines(a, b, ((1,),))
        with self.assertRaises(ValueError):
            intertwines(a, a, ((1, 0),))

    def test_direct_sum_honors_dimension_budget(self):
        g = cyclic_group(16)
        regular = regular_representation(g)
        with self.assertRaisesRegex(ValueError, "dimension budget"):
            direct_sum(regular, trivial_representation(g))

    def test_character_and_pairing_do_not_alias_mutable_callers(self):
        g = cyclic_group(2)
        rep = sign_c2(g)
        chars = rep.character()
        chars["g"] = Fraction(1_000_000)
        self.assertEqual(rep.character()["g"], -1)
        self.assertEqual(character_pairing(rep, rep), 1)

    def test_input_overflow_budget_is_rejected(self):
        g = cyclic_group(2)
        with self.assertRaises(ValueError):
            RationalRepresentation(g, {"e": ((1,),), "g": ((1 << 4097,),)})


if __name__ == "__main__":
    unittest.main()
