import unittest

from compare_smiles import compare_output, graph_key, parse, source_assertions


class SmilesComparison(unittest.TestCase):
    def test_graph_identity_detects_isotope_stereo_charge_and_bond_loss(self):
        for source, changed in [("[13CH3]O", "CO"), ("F[C@H](Cl)Br", "F[C@@H](Cl)Br"),
                                ("F/C=C/F", "FC=CF"), ("[NH4+]", "N"), ("C=C", "CC")]:
            expected = graph_key(parse(source))
            self.assertEqual(compare_output(expected, {"Ok": source})["status"], "match")
            self.assertEqual(compare_output(expected, {"Ok": changed})["status"], "difference")

    def test_reference_cleanup_cannot_masquerade_as_an_atom_permutation(self):
        self.assertEqual(source_assertions("F/C=C(/F)F")["double_bond"], 1)
        self.assertEqual(source_assertions("FC=C(F)F")["double_bond"], 0)
        self.assertEqual(source_assertions("F[C@H](Cl)Br")["tetrahedral"], 1)
        self.assertEqual(graph_key(parse("F/C=C(/F)F")), graph_key(parse("FC=C(F)F")))

    def test_only_unmapped_unlabelled_hydrogen_vertices_normalize(self):
        self.assertEqual(graph_key(parse("[H]C")), graph_key(parse("C")))
        self.assertNotEqual(graph_key(parse("[2H]C")), graph_key(parse("C")))
        self.assertNotEqual(graph_key(parse("[H:7]C")), graph_key(parse("C")))
        self.assertEqual(compare_output("C", {"Err": "unavailable"})["status"], "write_error")


if __name__ == "__main__":
    unittest.main()
