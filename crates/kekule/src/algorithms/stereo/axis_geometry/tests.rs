use super::*;
use crate::test_fixtures::read_smiles;

#[test]
fn trigonal_axis_endpoints_use_ligand_and_lone_pair_geometry() {
    for (source, expected) in [
        ("C(=C)(F)Cl", true),
        ("C1(Cl)=CC=CC=C1", true),
        ("N1(C)C=CC=C1", true),
        ("N(C)(C)C(=O)C", true),
        ("[C+](C)(C)C", true),
        ("B(C)(C)C", true),
        ("C1(C)CCCCC1", false),
        ("N1(C)CCCCC1", false),
        ("S(=O)(C)C1=CC=CC=C1", false),
        ("[Se](=O)(C)C1=CC=CC=C1", false),
        ("[S+](C)(C)C1=CC=CC=C1", false),
        ("P(C)(C)C1=CC=CC=C1", false),
        ("N(C)(C)S(=O)C", false),
    ] {
        let mol = read_smiles(source).unwrap();
        assert_eq!(
            atom_is_atropisomeric_sp2_endpoint(&mol, AtomId::new(0), 0),
            expected,
            "{source}"
        );
        assert!(!atom_is_atropisomeric_sp2_endpoint(&mol, AtomId::new(0), 1));
    }
}
