#[path = "support/smiles_contract.rs"]
mod contract;

use kekule::{core::Molecule, smiles};

fn molecule(input: &str) -> Molecule {
    let mut result = smiles::parse_str(input)
        .unwrap()
        .interpret()
        .unwrap()
        .into_molecule()
        .unwrap();
    result.perceive().unwrap();
    result
}

#[test]
fn fuzz_contract_covers_supported_writers_and_projections() {
    for input in [
        "CCO",
        "c1ccncc1",
        "[nH]1cccc1",
        "[NH4+]",
        "[CH3]",
        "[13CH3:7][C@H](O)F",
        "[2H]O[H]",
        "[H:4]OC",
        "F/C=C/Cl",
        "C1CC2CCC1C2",
        "F[C@H](Cl)[C@@H](Br)I |&1:1,3|",
        "F[C@H](Cl)[C@@H](Br)I |o1:1,3|",
        "F[C@H](Cl)Br |r|",
    ] {
        assert!(contract::check_molecule(&molecule(input)) >= 2, "{input}");
    }
}

#[test]
fn fuzz_contract_detects_invalid_output_composition_connectivity_and_stereo_loss() {
    for (source, bad_output) in [
        ("CCO", "C("),
        ("CCO", "[C]"),
        ("CCO", "C.C.O"),
        ("CCO", "COC"),
        ("[13CH3:7]CO", "CCO"),
        ("F/C=C/Cl", "FC=CCl"),
        ("F[C@H](Cl)[C@@H](Br)I |&1:1,3|", "F[C@H](Cl)[C@@H](Br)I"),
    ] {
        let original = molecule(source);
        let canonical = smiles::write_canonical(&original).unwrap();
        assert!(
            std::panic::catch_unwind(|| contract::assert_output(
                &original,
                bad_output,
                Some(&canonical)
            ))
            .is_err(),
            "undetected mutation: {bad_output}"
        );
    }
}
