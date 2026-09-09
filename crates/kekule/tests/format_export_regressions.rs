use kekule::core::{Atom, Element, HydrogenDeclaration, Molecule, MoleculeEditor};
use kekule::descriptors::{molecular_formula, HydrogenCountPolicy};
use kekule::molfile::{self, MolfileWriteVersion};
use kekule::sdf::{self, SdfRecordInterpretation, SdfWriteOptions};
use kekule::smiles;

fn molecule(source: &str) -> Molecule {
    smiles::to_molecules(source).unwrap().pop().unwrap()
}

#[test]
fn isomeric_export_preserves_aromatic_carbon_and_nitrogen_isotopes() {
    for source in ["[13CH]1=CC=CC=C1", "[13cH]1ccccc1", "[15nH]1cccc1"] {
        for perceive_first in [false, true] {
            let mut original = molecule(source);
            if perceive_first {
                original.perceive().unwrap();
            }
            let written = smiles::write_isomeric(&original).unwrap();
            let mut restored = molecule(&written);
            assert_eq!(
                original
                    .atoms()
                    .filter_map(|(_, atom)| atom.isotope)
                    .collect::<Vec<_>>(),
                restored
                    .atoms()
                    .filter_map(|(_, atom)| atom.isotope)
                    .collect::<Vec<_>>(),
                "{source} -> {written}"
            );
            original.perceive().unwrap();
            restored.perceive().unwrap();
            assert_eq!(
                molecular_formula(&original, HydrogenCountPolicy::IncludePerceived).unwrap(),
                molecular_formula(&restored, HydrogenCountPolicy::IncludePerceived).unwrap(),
                "{source} -> {written}"
            );
        }
    }
}

#[test]
fn bracket_metadata_export_preserves_the_perceived_hydrogen_count() {
    for (symbol, map, isotope, charge) in [
        ("C", Some(7), None, 0),
        ("C", None, Some(13), 0),
        ("N", None, None, 1),
        ("O", Some(3), Some(18), 0),
    ] {
        let mut atom = Atom::new(Element::from_symbol(symbol).unwrap());
        atom.atom_map = map;
        atom.isotope = isotope;
        atom.formal_charge = charge;
        let mut editor = MoleculeEditor::new();
        editor.add_atom(atom).unwrap();
        let mut original = editor.finish().unwrap();
        for writer in [smiles::write, smiles::write_isomeric] {
            let error = writer(&original).unwrap_err();
            assert!(error.to_string().contains("hydrogen perception"));
        }
        if map.is_some() || charge != 0 {
            let error = smiles::write_canonical(&original).unwrap_err();
            assert!(error.to_string().contains("hydrogen perception"));
        }
        original.perceive().unwrap();
        for writer in [smiles::write, smiles::write_isomeric] {
            let written = writer(&original).unwrap();
            let mut restored = molecule(&written);
            restored.perceive().unwrap();
            assert_eq!(
                molecular_formula(&original, HydrogenCountPolicy::IncludePerceived).unwrap(),
                molecular_formula(&restored, HydrogenCountPolicy::IncludePerceived).unwrap(),
                "{symbol} -> {written}"
            );
            assert_eq!(restored.atoms().next().unwrap().1.atom_map, map);
        }
    }
}

#[test]
fn canonical_hydrogen_collapse_requires_a_known_parent_count() {
    let mut inferred = molecule("[H]C");
    let error = smiles::write_canonical(&inferred).unwrap_err();
    assert!(error.to_string().contains("hydrogen perception"));
    inferred.perceive().unwrap();
    for original in [inferred, molecule("[H][CH3]")] {
        let written = smiles::write_canonical(&original).unwrap();
        assert_eq!(written, "C");
        let mut restored = molecule(&written);
        restored.perceive().unwrap();
        assert_eq!(
            molecular_formula(&original, HydrogenCountPolicy::IncludePerceived).unwrap(),
            molecular_formula(&restored, HydrogenCountPolicy::IncludePerceived).unwrap()
        );
    }
}

#[test]
fn canonical_bracket_hydrogens_are_materialized_after_explicit_perception() {
    let mut editor = MoleculeEditor::new();
    let mut carbon = Atom::new(Element::from_symbol("C").unwrap());
    carbon.atom_map = Some(7);
    editor.add_atom(carbon).unwrap();
    let mut methane = editor.finish().unwrap();
    assert!(smiles::write_canonical(&methane).is_err());
    methane.perceive().unwrap();
    let written = smiles::write_canonical(&methane).unwrap();
    assert_eq!(written, "[CH4:7]");
    let mut restored = molecule(&written);
    restored.perceive().unwrap();
    assert_eq!(
        molecular_formula(&methane, HydrogenCountPolicy::IncludePerceived).unwrap(),
        molecular_formula(&restored, HydrogenCountPolicy::IncludePerceived).unwrap()
    );
}

#[test]
fn unperceived_organic_and_fixed_bracket_atoms_remain_writeable() {
    for source in [
        "C",
        "N",
        "O",
        "F",
        "[NH4+]",
        "[Na+]",
        "[13CH4]",
        "N[C@H](O)C",
        "F/C=C/F",
    ] {
        let original = molecule(source);
        let written = smiles::write_isomeric(&original).unwrap();
        let restored = molecule(&written);
        assert_eq!(restored.formal_charge(), original.formal_charge());
        assert_eq!(
            restored.stereo_elements().count(),
            original.stereo_elements().count()
        );
    }
}

#[test]
fn minimum_formal_charge_round_trips_without_overflow() {
    let parsed = molecule("[C-128]");
    let mut atom = Atom::new(Element::from_symbol("C").unwrap());
    atom.formal_charge = i8::MIN;
    atom.hydrogens = HydrogenDeclaration::Fixed(0);
    let mut editor = MoleculeEditor::new();
    editor.add_atom(atom).unwrap();
    let constructed = editor.finish().unwrap();
    for original in [parsed, constructed] {
        for writer in [
            smiles::write,
            smiles::write_isomeric,
            smiles::write_canonical,
        ] {
            let written = writer(&original).unwrap();
            assert_eq!(molecule(&written).formal_charge(), i64::from(i8::MIN));
        }
    }
}

#[test]
fn molfile_control_text_in_positional_headers_is_not_a_terminator() {
    let original = molecule("CO");
    for writer in [molfile::write_v2000, molfile::write_v3000] {
        let written = writer(&original).unwrap();
        for header in 0..3 {
            let mut lines = written.lines().collect::<Vec<_>>();
            lines[header] = "M  END";
            let document = molfile::parse_str(&lines.join("\n")).unwrap();
            assert_eq!(
                document.to_molecules().unwrap().pop().unwrap().atom_count(),
                2
            );
        }
    }
}

#[test]
fn sdf_control_text_title_round_trips_or_is_rejected_before_writing() {
    let mol = molfile::write_v2000(&molecule("CO")).unwrap();
    let model = molfile::parse_str(&mol).unwrap().to_model().unwrap();
    for version in [MolfileWriteVersion::V2000, MolfileWriteVersion::V3000] {
        let options = SdfWriteOptions { version };
        let record = SdfRecordInterpretation::new("M  END", model.clone(), Vec::new());
        let written = sdf::write_records(&[record], options).unwrap();
        let document = sdf::parse_str(&written).unwrap();
        assert_eq!(document.records().len(), 1);
        assert_eq!(document.records()[0].title(), "M  END");
        assert_eq!(
            document.records()[0]
                .to_model()
                .unwrap()
                .topology()
                .atom_count(),
            2
        );
        for title in ["$$$$", " $$$$ "] {
            let record = SdfRecordInterpretation::new(title, model.clone(), Vec::new());
            let mut output = Vec::new();
            assert!(sdf::write_records_to(&mut output, &[record], options).is_err());
            assert!(output.is_empty());
        }
    }
}
