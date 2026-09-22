use kekule_openff::ForceField;

const ROSEMARY: &str = include_str!("../data/rosemary.offxml");

#[test]
fn file_loader_uses_literal_paths_and_reports_path_on_failure() {
    let path = std::env::temp_dir().join(format!(
        "kekule custom offxml {}.offxml",
        std::process::id()
    ));
    std::fs::write(&path, ROSEMARY).unwrap();
    let loaded = ForceField::from_file(&path).unwrap();
    assert_eq!(
        loaded.nonbonded_settings(),
        ForceField::rosemary().unwrap().nonbonded_settings()
    );
    std::fs::write(&path, "<SMIRNOFF>").unwrap();
    assert!(ForceField::from_file(&path)
        .unwrap_err()
        .to_string()
        .contains(path.to_str().unwrap()));
    std::fs::write(&path, [0xff, 0xfe]).unwrap();
    assert!(ForceField::from_file(&path).is_err());
    std::fs::remove_file(&path).unwrap();
    assert!(ForceField::from_file(&path)
        .unwrap_err()
        .to_string()
        .contains(path.to_str().unwrap()));
}

#[test]
fn deferred_charging_and_unsupported_semantics_still_fail() {
    for (changed, diagnostic) in [
        (
            ROSEMARY.replace("openff-gnn-am1bcc-1.0.0.pt", "custom.pt"),
            "model_file",
        ),
        (ROSEMARY.replace("7981e7f5", "0981e7f5"), "model_file_hash"),
        (
            ROSEMARY
                .replace("<NAGLCharges", "<ToolkitAM1BCC")
                .replace("</NAGLCharges", "</ToolkitAM1BCC"),
            "ToolkitAM1BCC",
        ),
        (
            ROSEMARY.replace("<vdW version=\"0.4\"", "<vdW version=\"0.5\""),
            "version",
        ),
        (
            ROSEMARY.replace("periodic_method=\"cutoff\"", "periodic_method=\"garbage\""),
            "periodic_method",
        ),
        (
            ROSEMARY.replace(
                "nonperiodic_potential=\"Coulomb\"",
                "nonperiodic_potential=\"garbage\"",
            ),
            "nonperiodic_potential",
        ),
        (
            ROSEMARY.replace(
                "fractional_bondorder_interpolation=\"linear\"",
                "fractional_bondorder_interpolation=\"cubic\"",
            ),
            "interpolation",
        ),
        (
            ROSEMARY.replace("<vdW version", "<vdW method=\"cutoff\" version"),
            "method",
        ),
        (
            ROSEMARY.replace(
                "<Electrostatics version",
                "<Electrostatics method=\"PME\" version",
            ),
            "method",
        ),
        (
            ROSEMARY.replace("<Bonds version", "<Bonds mystery=\"value\" version"),
            "mystery",
        ),
        (
            ROSEMARY.replace("length=\"", "length_bondorder1=\""),
            "fractional-bond-order",
        ),
    ] {
        assert_ne!(changed, ROSEMARY, "test must modify its fixture");
        let message = ForceField::from_offxml(&changed).unwrap_err().to_string();
        assert!(message.contains(diagnostic), "{diagnostic}: {message}");
    }
    let document = roxmltree::Document::parse(ROSEMARY).unwrap();
    let nagl = document
        .descendants()
        .find(|n| n.has_tag_name("NAGLCharges"))
        .unwrap()
        .range();
    let mut missing = ROSEMARY.to_owned();
    missing.replace_range(nagl, "");
    assert!(ForceField::from_offxml(&missing)
        .unwrap_err()
        .to_string()
        .contains("missing handler NAGLCharges"));
}

#[test]
fn malformed_optional_sections_are_not_silently_ignored() {
    for xml in [
        ROSEMARY.replace(
            "<Constraints version=\"0.3\"",
            "<Constraints version=\"999\"",
        ),
        ROSEMARY.replace(
            "<LibraryCharges version=\"0.3\"",
            "<LibraryCharges version=\"999\"",
        ),
        ROSEMARY.replace(
            "<ImproperTorsions version=\"0.3\"",
            "<ImproperTorsions version=\"999\"",
        ),
        ROSEMARY.replace(
            "</Constraints>",
            "</Constraints><Constraints version=\"0.3\"/>",
        ),
        ROSEMARY.replace("<Constraints version", "<Constraints unknown=\"x\" version"),
        ROSEMARY.replace("<SMIRNOFF ", "<SMIRNOFF xmlns=\"urn:other\" "),
    ] {
        assert_ne!(xml, ROSEMARY);
        assert!(ForceField::from_offxml(&xml).is_err());
    }
}

#[test]
fn rule_metadata_is_optional_but_smirks_and_parameter_values_are_required() {
    let xml = ROSEMARY.replace("id=\"b1\"", "name=\"custom bond\"");
    assert_ne!(xml, ROSEMARY);
    let ff = ForceField::from_offxml(&xml).unwrap();
    let molecule = kekule::smiles::to_molecules("[H][H]").unwrap().remove(0);
    // Explicitly append a final anonymous rule to ensure IDs are not assignment keys.
    let xml = xml.replace("</Bonds>", "<Bond smirks=\"[#1:1]-[#1:2]\" length=\".74*angstrom\" k=\"1*kcal/(mol*angstrom**2)\"/></Bonds>");
    let labels = ForceField::from_offxml(&xml)
        .unwrap()
        .label_molecule(&molecule)
        .unwrap();
    assert!(labels.bonds.values().all(|p| p.id.is_empty()));
    assert_eq!(
        ff.nonbonded_settings(),
        ForceField::rosemary().unwrap().nonbonded_settings()
    );
    assert!(ForceField::from_offxml(&ROSEMARY.replace("length=", "missing_length=")).is_err());
}
