use super::*;

fn perceive(molecule: &mut Molecule) {
    molecule.perceive().expect("molecule should be perceived");
}

fn one_smiles(input: &str) -> Result<Molecule, String> {
    let mut molecules = kekule::smiles::to_molecules(input).map_err(|error| error.to_string())?;
    if molecules.len() != 1 {
        return Err(format!(
            "expected one connected molecule, found {}",
            molecules.len()
        ));
    }
    Ok(molecules.pop().expect("component count was checked"))
}

use std::time::{SystemTime, UNIX_EPOCH};

fn disconnected_sdf_record() -> String {
    // Two components with interleaved atoms and bonds in source order.
    "disconnected
  benchmark-regression

  5  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  3  5  1  0  0  0  0
  2  4  1  0  0  0  0
  1  3  1  0  0  0  0
M  CHG  1   2  -1
M  END
>  <ID>
salt-regression

$$$$
"
    .to_owned()
}

#[test]
fn sdf_parse_preserves_disconnected_components_and_source_order() {
    let input = Input {
        path: PathBuf::from("disconnected.sdf"),
        text: disconnected_sdf_record() + &simple_sdf_record("next record"),
    };
    let output = evaluate("io.sdf.v2000.parse", "smoke", &input).unwrap();
    let records = output["records"].as_array().unwrap();
    assert_eq!(records.len(), 2);
    let record = &records[0];
    assert_eq!(record["record_index"], 0);
    assert_eq!(record["status"], "ok");
    assert_eq!(record["title"], "disconnected");
    assert_eq!(record["atom_count"], 5);
    assert_eq!(record["bond_count"], 3);
    assert_eq!(record["properties"], json!({"ID": "salt-regression"}));
    let atoms = record["atoms"].as_array().unwrap();
    assert_eq!(atoms.len(), 5);
    for (index, symbol) in ["C", "O", "N", "C", "F"].iter().enumerate() {
        assert_eq!(atoms[index]["index"], index);
        assert_eq!(atoms[index]["symbol"], *symbol);
        assert_eq!(
            atoms[index]["formal_charge"],
            if index == 1 { -1 } else { 0 }
        );
    }
    assert_eq!(
        record["bonds"],
        json!([
            {"index":0,"begin_atom_index":2,"end_atom_index":4,"bond_type":"SINGLE",
             "is_aromatic":false,"stereo":"STEREONONE","bond_direction":"NONE"},
            {"index":1,"begin_atom_index":1,"end_atom_index":3,"bond_type":"SINGLE",
             "is_aromatic":false,"stereo":"STEREONONE","bond_direction":"NONE"},
            {"index":2,"begin_atom_index":0,"end_atom_index":2,"bond_type":"SINGLE",
             "is_aromatic":false,"stereo":"STEREONONE","bond_direction":"NONE"}
        ])
    );
    assert_eq!(records[1]["record_index"], 1);
    assert_eq!(records[1]["title"], "next record");
    assert_eq!(records[1]["atom_count"], 1);
    assert_eq!(records[1]["atoms"][0]["index"], 0);
}

#[test]
fn single_molecule_features_return_errors_for_disconnected_sdf() {
    let mut input = Input {
        path: PathBuf::from("disconnected.sdf"),
        text: disconnected_sdf_record(),
    };
    for feature in ["descriptor.molecular", "io.sdf.v2000.write"] {
        let error = evaluate(feature, "smoke", &input).unwrap_err();
        assert!(error
            .to_string()
            .contains("expected one connected molecule, found 2"));
    }
    input.text = simple_sdf_record("next record");
    for feature in ["descriptor.molecular", "io.sdf.v2000.write"] {
        let output = evaluate(feature, "smoke", &input).unwrap();
        assert_eq!(output["records"][0]["status"], "ok");
    }
}

#[test]
fn implementation_dispatch_uses_current_molfile_benchmark_ids() {
    let root = temp_workspace_root("mol-benchmark-dispatch");
    let fixture = root.join("fixture.sdf");
    fs::write(&fixture, simple_sdf_record("methane")).expect("fixture should write");

    for benchmark_id in [
        "io.mol.v2000.parse",
        "io.mol.v2000.write",
        "io.mol.v3000.parse",
        "io.mol.v3000.write",
    ] {
        let expected = evaluate_path(benchmark_id, "pubchem-100k", &fixture)
            .expect("benchmark should compare");
        assert_eq!(expected["records"][0]["status"], "ok");
    }

    fs::remove_dir_all(root).ok();
}

#[test]
fn implementation_dispatch_supports_mmcif_document_rows() {
    let root = temp_workspace_root("mmcif-document-dispatch");
    let fixture = root.join("fixture.cif");
    fs::write(
        &fixture,
        r#"data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA CA . ALA ALA A A 1 1 ? 1.00 10.00 1.0 2.0 3.0 1
"#,
    )
    .expect("fixture should write");

    let expected = evaluate_path("io.mmcif.parse", "pdb-1000", &fixture)
        .expect("mmCIF document benchmark should compare");
    let atom_site = &expected["atom_site_rows"];
    assert_eq!(atom_site["status"], "ok");
    assert_eq!(atom_site["row_count"], 1);
    assert_eq!(atom_site["rows"][0]["id"], "1");
    assert_eq!(atom_site["rows"][0]["label_alt_id"], Value::Null);
    assert_eq!(atom_site["rows"][0]["pdbx_PDB_ins_code"], Value::Null);
    assert_eq!(atom_site["rows"][0]["Cartn_z"], "3.0");

    fs::remove_dir_all(root).ok();
}
#[test]
fn implementation_dispatch_supports_hydrogen_transforms() {
    let root = temp_workspace_root("hydrogen-transforms-dispatch");
    let fixture = root.join("fixture.sdf");
    fs::write(&fixture, simple_sdf_record("methane")).expect("fixture should write");

    let expected = evaluate_path("chem.hydrogen-transforms", "pubchem-100k", &fixture)
        .expect("benchmark should compare");
    let record = &expected["records"][0];

    assert_eq!(record["status"], "ok");
    assert_eq!(record["atom_count_after_add"], 5);
    assert_eq!(
        record["added_hydrogens_by_parent"],
        json!([{ "parent_atom_index": 0, "count": 4 }])
    );
    assert_eq!(record["round_trip"]["status"], "ok");

    fs::remove_dir_all(root).ok();
}

#[test]
fn implementation_dispatch_supports_query_benchmark() {
    let root = temp_workspace_root("query-benchmark-dispatch");
    let smarts_fixture = root.join("fixture.smi");
    fs::write(&smarts_fixture, "CCO\nC1=CC=CC=C1\n").expect("fixture should write");

    let parsed = evaluate_path("query.smarts", "pubchem-100k", &smarts_fixture)
        .expect("SMARTS benchmark should compare");
    assert_eq!(parsed["records"][0]["status"], "ok");
    assert_eq!(parsed["records"][0]["atom_count"], 3);
    assert_eq!(parsed["records"][1]["bond_count"], 6);

    let molecule_fixture = root.join("fixture.sdf");
    fs::write(&molecule_fixture, simple_sdf_record("methane")).expect("fixture should write");
    let matched = evaluate_path("algo.substructure.vf2", "pubchem-100k", &molecule_fixture)
        .expect("substructure benchmark should compare");
    assert_eq!(matched["records"][0]["status"], "ok");
    assert_eq!(matched["records"][0]["queries"][0]["smarts"], "[#6]");
    assert_eq!(matched["records"][0]["queries"][0]["matches"], json!([[0]]));

    fs::remove_dir_all(root).ok();
}

#[test]
fn implementation_dispatch_uses_current_isomeric_smiles_benchmark_id() {
    let root = temp_workspace_root("isomeric-smiles-benchmark-dispatch");
    let fixture = root.join("fixture.smi");
    fs::write(
        &fixture,
        [
            "CCO CID:plain",
            "C[C@@H](C(=O)O)N CID:tetrahedral",
            "C(=C\\F)\\F CID:double-bond",
        ]
        .join("\n"),
    )
    .expect("fixture should write");

    let expected = evaluate_path("io.smiles.isomeric", "pubchem-100k", &fixture)
        .expect("benchmark should compare");
    let records = expected["records"]
        .as_array()
        .expect("records should be an array");

    assert_eq!(records.len(), 3);
    assert!(records.iter().all(|record| record["status"] == "ok"));
    assert_eq!(records[0]["input_smiles"], "CCO");
    assert_eq!(records[0]["stereo"]["atom_descriptors"], json!([]));
    assert_eq!(records[0]["stereo"]["bond_descriptors"], json!([]));
    assert!(!records[1]["stereo"]["atom_descriptors"]
        .as_array()
        .expect("atom descriptors should be an array")
        .is_empty());
    assert!(!records[2]["stereo"]["bond_descriptors"]
        .as_array()
        .expect("bond descriptors should be an array")
        .is_empty());

    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_representation_benchmark_preserves_components_and_failed_records() {
    let root = temp_workspace_root("stereo-representation-components");
    let fixture = root.join("fixture.smi");
    fs::write(
        &fixture,
        "[Na+].F[C@H](Cl)Br.F/C=C/F CID:components\nF[C@ CID:invalid\n",
    )
    .expect("fixture should write");
    let output = evaluate_path("stereo.representation", "pubchem-100k", &fixture).unwrap();
    let records = output["records"].as_array().unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0]["record_index"], 0);
    assert_eq!(records[0]["status"], "ok");
    assert_eq!(records[0]["atom_count"], 9);
    assert_eq!(records[0]["bond_count"], 6);
    let elements = records[0]["stereo_elements"].as_array().unwrap();
    assert_eq!(elements.len(), 2);
    assert_eq!(elements[0]["index"], 0);
    assert_eq!(elements[0]["center_atom_index"], 2);
    assert_eq!(
        elements[0]["carriers"],
        json!([
            {"atom_index": 1}, {"atom_index": 3}, {"atom_index": 4}, {"implicit_hydrogen": true}
        ])
    );
    assert_eq!(elements[1]["index"], 1);
    assert_eq!(elements[1]["center_bond_index"], 4);
    assert_eq!(elements[1]["left_atom_index"], 6);
    assert_eq!(elements[1]["right_atom_index"], 7);
    assert_eq!(elements[1]["left_carrier"], json!({"atom_index": 5}));
    assert_eq!(elements[1]["right_carrier"], json!({"atom_index": 8}));
    assert_eq!(records[1]["record_index"], 1);
    assert_eq!(records[1]["status"], "parse_error");
    assert_eq!(records[1]["title"], "CID:invalid");
    fs::remove_dir_all(root).ok();
}

#[test]
fn isomeric_smiles_benchmark_preserves_components_and_failures() {
    let root = temp_workspace_root("isomeric-smiles-components");
    let fixture = root.join("fixture.smi");
    fs::write(
        &fixture,
        [
            "[Na+].F[C@H](Cl)Br.F/C=C/F CID:components",
            "F[C@ CID:invalid",
            "F[C@](Cl)(Br)[CH5] CID:invalid-valence",
            "CCO CID:non-stereo-control",
        ]
        .join("\n"),
    )
    .expect("fixture should write");
    let output = evaluate_path("io.smiles.isomeric", "pubchem-100k", &fixture).unwrap();
    let records = output["records"].as_array().unwrap();
    assert_eq!(records.len(), 4);
    assert_eq!(records[0]["status"], "ok");
    assert_eq!(records[0]["normalized_perceived"]["atom_count"], 9);
    assert_eq!(records[0]["normalized_perceived"]["bond_count"], 6);
    assert_eq!(records[0]["stereo"]["status"], "ok");
    let atoms = records[0]["stereo"]["atom_descriptors"].as_array().unwrap();
    let bonds = records[0]["stereo"]["bond_descriptors"].as_array().unwrap();
    assert_eq!(atoms.len(), 1);
    assert_eq!(atoms[0]["descriptor"], "R");
    assert_eq!(bonds.len(), 1);
    assert_eq!(bonds[0]["descriptor"], "E");
    assert_eq!(records[1]["record_index"], 1);
    assert_eq!(records[1]["status"], "parse_error");
    assert_eq!(records[2]["record_index"], 2);
    assert_eq!(records[2]["status"], "perception_error");
    assert_eq!(records[3]["record_index"], 3);
    assert_eq!(records[3]["status"], "ok");
    assert_eq!(records[3]["input_smiles"], "CCO");
    assert_eq!(records[3]["stereo"]["atom_descriptors"], json!([]));
    assert_eq!(records[3]["stereo"]["bond_descriptors"], json!([]));
    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_benchmark_keeps_source_marks_in_the_document_and_coordinates_in_the_model() {
    let root = temp_workspace_root("stereo-source-and-geometry");
    let smiles_path = root.join("source.smi");
    fs::write(&smiles_path, "F/C=C\\F source\n").unwrap();
    let value = evaluate_path("stereo.representation", "smoke", &smiles_path).unwrap();
    let record = &value["records"][0];
    assert_eq!(record["document"]["source"], "F/C=C\\F");
    assert_eq!(
        record["document"]["stereo_bond_marks"],
        json!([
            {"bond_index": 0, "kind": "directional_up", "source": "smiles"},
            {"bond_index": 2, "kind": "directional_down", "source": "smiles"},
        ])
    );
    assert_eq!(
        record["document"]["stereo_sources"],
        json!([
            {"element_index": 0, "source": "smiles", "specifiedness": "specified"},
        ])
    );
    assert_eq!(record["stereo_elements"][0]["type"], "double_bond");
    assert!(record.get("stereo_bond_marks").is_none());
    assert!(record["stereo_elements"][0].get("source").is_none());

    let molecule = smiles::to_molecules("FC(Cl)(Br)I").unwrap().pop().unwrap();
    let mut orientations = Vec::new();
    for reflection in [1.0, -1.0] {
        let points = [
            [1.0, 1.0, 1.0],
            [0.0, 0.0, 0.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ]
        .map(|p| kekule::geometry::Point3::new(p[0] * reflection, p[1], p[2]))
        .to_vec();
        let positions = kekule::structure::Positions::new(kekule::units::Quantity::new(
            points,
            kekule::units::ANGSTROM,
        ))
        .unwrap();
        let model = kekule::structure::Model::from_molecule(&molecule, &positions).unwrap();
        let text = molfile::write_model(&model, molfile::MolfileWriteOptions::default()).unwrap();
        let path = root.join("geometry.mol");
        fs::write(&path, &text).unwrap();
        let value = evaluate_path("stereo.perception", "smoke", &path).unwrap();
        let record = &value["records"][0];
        assert_eq!(record["status"], "ok");
        assert_eq!(record["report"]["created_element_indices"], json!([0]));
        assert_eq!(record["document"]["source"], text);
        assert_eq!(
            record["document"]["stereo_sources"][0]["source"],
            "coordinates_3d"
        );
        orientations.push(record["stereo_elements"][0]["orientation"].clone());
    }
    assert_ne!(orientations[0], orientations[1]);
    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_cip_benchmark_retains_records_without_descriptors_and_parse_failures() {
    let root = temp_workspace_root("stereo-cip-descriptor-filter");
    let fixture = root.join("fixture.smi");
    fs::write(
        &fixture,
        [
            "CC CID:no-stereo",
            "F[C@ CID:invalid",
            "C[C@H](N)C(=O)O CID:stereo",
        ]
        .join("\n"),
    )
    .expect("fixture should write");

    let expected =
        evaluate_path("stereo.cip", "pubchem-100k", &fixture).expect("benchmark should compare");
    let records = expected["records"]
        .as_array()
        .expect("records should be an array");

    assert_eq!(records.len(), 3);
    assert_eq!(records[0]["title"], "CID:no-stereo");
    assert_eq!(records[0]["status"], "ok");
    assert_eq!(records[0]["atom_count"], 2);
    assert_eq!(records[0]["atom_descriptors"], json!([]));
    assert_eq!(records[0]["bond_descriptors"], json!([]));
    assert_eq!(records[1]["record_index"], 1);
    assert_eq!(records[1]["status"], "parse_error");
    assert_eq!(records[2]["title"], "CID:stereo");
    assert!(!records[2]["atom_descriptors"]
        .as_array()
        .expect("atom descriptors should be an array")
        .is_empty());

    fs::remove_dir_all(root).ok();
}

#[test]
fn smiles_features_attempt_stereo_wildcards_and_invalid_inputs() {
    let input = Input {
        path: PathBuf::from("regression.smi"),
        text: "F[C@H](Cl)Br atom\nF/C=C/F bond\nF\\C=C\\F reverse\n* wildcard\nF[C@ invalid\n"
            .into(),
    };
    for feature in [
        "io.smiles.parse",
        "io.smiles.isomeric",
        "io.smiles.canonical",
    ] {
        let output = evaluate(feature, "smoke", &input).unwrap();
        let records = output["records"].as_array().unwrap();
        assert_eq!(records.len(), 5, "{feature}");
        for (index, record) in records.iter().enumerate() {
            assert_eq!(record["record_index"], index, "{feature}");
            assert_eq!(
                record["status"],
                if index < 3 { "ok" } else { "parse_error" },
                "{feature}: {record}"
            );
        }
    }
    for line in input.text.lines() {
        let record = Input {
            path: input.path.clone(),
            text: line.into(),
        };
        let result = evaluate("io.smiles.write", "smoke", &record);
        if line.contains('@') && !line.ends_with("invalid")
            || line.contains('/')
            || line.contains('\\')
        {
            assert!(result.unwrap_err().to_string().contains("stereochemistry"));
        } else {
            assert_eq!(result.unwrap()["records"][0]["status"], "parse_error");
        }
    }
    for feature in ["stereo.cip", "stereo.representation", "stereo.perception"] {
        let wildcard = Input {
            path: input.path.clone(),
            text: "* wildcard\n".into(),
        };
        let output = evaluate(feature, "smoke", &wildcard).unwrap();
        assert_eq!(output["records"][0]["status"], "parse_error", "{feature}");
    }
}

#[test]
fn stereo_cip_benchmark_preserves_disconnected_record_indices() {
    let root = temp_workspace_root("stereo-cip-components");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "[Na+].C[C@H](N)C(=O)O.[Cl-] CID:salt\n").expect("fixture should write");

    let expected = evaluate_path("stereo.cip", "pubchem-100k", &fixture)
        .expect("disconnected record should compare");
    let records = expected["records"].as_array().expect("record array");
    assert_eq!(records.len(), 1);
    assert_eq!(records[0]["atom_count"], 8);
    assert_eq!(records[0]["bond_count"], 5);
    assert_eq!(
        records[0]["atom_descriptors"],
        json!([
            {"atom_index": 2, "descriptor": "S"}
        ])
    );

    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_cip_benchmark_reports_perception_failure() {
    let root = temp_workspace_root("stereo-cip-perception-failure");
    let fixture = root.join("fixture.smi");
    for smiles in ["F[C@](Cl)(Br)[CH5]", "C(C)(C)(C)(C)C"] {
        fs::write(&fixture, format!("{smiles} CID:invalid-valence\n"))
            .expect("fixture should write");
        let error = evaluate_path("stereo.cip", "pubchem-100k", &fixture)
            .expect_err("failed perception must remain visible, with or without stereo");
        assert!(error.to_string().contains("record 0 perception failed"));
    }

    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_cip_benchmark_retains_isolated_hydrogen_component() {
    let root = temp_workspace_root("stereo-cip-isolated-hydrogen");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "[HH].C[C@H](N)C(=O)O CID:hydrogen-component\n")
        .expect("fixture should write");

    let expected =
        evaluate_path("stereo.cip", "pubchem-100k", &fixture).expect("record should compare");
    assert_eq!(expected["records"][0]["atom_count"], 7);
    assert_eq!(
        expected["records"][0]["atom_descriptors"][0]["atom_index"],
        2
    );
    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_cip_benchmark_uses_rdkit_default_hydrogen_indexing() {
    let root = temp_workspace_root("stereo-cip-rdkit-h-index");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "[H][C@](F)(Cl)Br CID:explicit-h\n").expect("fixture should write");

    let expected =
        evaluate_path("stereo.cip", "pubchem-100k", &fixture).expect("benchmark should compare");
    let records = expected["records"]
        .as_array()
        .expect("records should be an array");

    assert_eq!(records.len(), 1);
    assert_eq!(records[0]["atom_count"], 4);
    assert_eq!(records[0]["bond_count"], 3);
    assert_eq!(records[0]["atom_descriptors"][0]["atom_index"], 0);

    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_cip_benchmark_reads_all_sdf_pack_records() {
    let root = temp_workspace_root("stereo-cip-sdf-pack");
    let fixture = root.join("fixture.sdf");
    fs::write(
        &fixture,
        [
            chiral_wedge_sdf_record("first"),
            chiral_wedge_sdf_record("second"),
        ]
        .join(""),
    )
    .expect("fixture should write");

    let expected =
        evaluate_path("stereo.cip", "pubchem-100k", &fixture).expect("benchmark should compare");
    let records = expected["records"]
        .as_array()
        .expect("records should be an array");

    assert_eq!(records.len(), 2);
    assert_eq!(records[0]["title"], "first");
    assert_eq!(records[1]["title"], "second");
    assert!(records
        .iter()
        .all(|record| record["atom_count"].as_u64() == Some(5)));
    assert!(records.iter().all(|record| !record["atom_descriptors"]
        .as_array()
        .expect("atom descriptors should be an array")
        .is_empty()));

    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_perception_benchmark_records_reference_preparation_errors_per_record() {
    let molecule = one_smiles("C(C)(C)(C)(C)C")
        .expect("pentavalent neutral carbon should remain an interpretation-valid graph");
    let mut record = IndexedSmallRecord {
        record_index: 0,
        title: "invalid neutral-carbon valence".to_owned(),
        molecule,
        sdf_fields: BTreeMap::new(),
    };

    let value = stereo_perception_record_json(&mut record);

    assert_eq!(
        value.get("status").and_then(Value::as_str),
        Some("perception_error")
    );
    assert!(value.get("report").is_none());
}

#[test]
fn smiles_component_benchmarks_preserve_source_record_cardinality() {
    let root = temp_workspace_root("smiles-component-benchmark-cardinality");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "CC.Cl.Cl multi\nC=C connected\n").expect("fixture should write");

    let parsed = evaluate_path("io.smiles.parse", "pubchem-100k", &fixture)
        .expect("parse benchmark should serialize");
    assert_eq!(parsed["records"].as_array().map(Vec::len), Some(2));
    assert_eq!(parsed["records"][0]["status"], "ok");
    assert_eq!(parsed["records"][0]["raw"]["atom_count"], 4);
    assert_eq!(parsed["records"][0]["raw"]["bond_count"], 1);
    assert!(parsed["records"][0].get("normalized_perceived").is_some());
    assert!(parsed["records"][0].get("write_round_trip").is_some());

    let written = evaluate_path("io.smiles.write", "pubchem-100k", &fixture)
        .expect("write benchmark should serialize");
    assert_eq!(written["records"].as_array().map(Vec::len), Some(2));
    assert_eq!(written["records"][0]["status"], "ok");
    assert_eq!(
        written["records"][0]["normalized_perceived"]["atom_count"],
        4
    );
    assert_eq!(
        written["records"][0]["normalized_perceived"]["bond_count"],
        1
    );

    let records =
        read_smiles_records(&Input::read(&fixture).unwrap()).expect("fixture should interpret");
    let reparsed = records[0]
        .components
        .iter()
        .map(|molecule| {
            let text = smiles::write(molecule).expect("component should write");
            let document = smiles::parse_str(&text).expect("written component should parse");
            smiles::interpret(&document)
                .expect("written component should interpret")
                .into_molecule()
                .expect("written component should remain connected")
        })
        .collect::<Vec<_>>();
    assert_eq!(
        written["records"][0]["normalized_perceived"],
        smiles_components_perceived_semantic_json(&reparsed)
    );
    let connected_written = smiles::write(
        records[1]
            .components
            .first()
            .expect("connected record should have one molecule"),
    )
    .expect("connected record should write");
    let connected_document =
        smiles::parse_str(&connected_written).expect("written component should parse");
    let connected_reparsed = smiles::interpret(&connected_document)
        .expect("written component should interpret")
        .into_molecule()
        .expect("written component should remain connected");
    assert_eq!(
        written["records"][1]["normalized_perceived"],
        smiles_perceived_semantic_json(connected_reparsed)
    );

    let skipped = IndexedSmilesRecord {
        record_index: 0,
        status: "ok".to_owned(),
        title: "missing components".to_owned(),
        input_smiles: "CC.Cl.Cl".to_owned(),
        components: Vec::new(),
    };
    let skipped = smiles_write_record_json(&skipped).expect("error record should serialize");
    assert_ne!(skipped["status"], "ok");

    let stereo = evaluate_path("stereo.perception", "pubchem-100k", &fixture)
        .expect("stereo benchmark should serialize");
    assert_eq!(stereo["records"].as_array().map(Vec::len), Some(2));
    assert_eq!(stereo["records"][0]["status"], "ok");
    assert_eq!(stereo["records"][0]["atom_count"], 4);
    assert_eq!(stereo["records"][0]["bond_count"], 1);
    assert!(stereo["records"][0]["report"].get("candidates").is_some());
    assert!(stereo["records"][0]
        .get("source_stereo_element_indices")
        .is_none());

    let stereo_fixture = root.join("stereo.smi");
    fs::write(&stereo_fixture, "F/C=C/F.F/C=C/F directional\n")
        .expect("stereo fixture should write");
    let stereo = evaluate_path("stereo.perception", "pubchem-100k", &stereo_fixture)
        .expect("component stereo benchmark should serialize");
    assert_eq!(stereo["records"].as_array().map(Vec::len), Some(1));
    assert_eq!(
        stereo["records"][0]["report"]["assembled_elements"]
            .as_array()
            .map(Vec::len),
        Some(2)
    );
    assert_eq!(
        stereo["records"][0]["report"]["assembled_elements"][0]["index"],
        0
    );
    assert_eq!(
        stereo["records"][0]["report"]["assembled_elements"][1]["index"],
        1
    );

    fs::remove_dir_all(root).ok();
}

#[test]
fn dssp_comparison_matches_residues_by_source_identity_not_container_order() {
    let mut expected = json!({
        "status": "ok",
        "residues": [
            {"chain_id": "B", "sequence_id": 1, "insertion_code": null, "label_chain_id": "B", "label_sequence_id": 1, "residue_name": "ALA", "sheet": 4, "strand": 8, "ladders": [19, null]},
            {"chain_id": "D", "sequence_id": 1, "insertion_code": null, "label_chain_id": "D", "label_sequence_id": 1, "residue_name": "GLY", "sheet": 7, "strand": 9, "ladders": [21, 19]}
        ]
    });
    let mut actual = json!({
        "status": "ok",
        "residues": [
            {"chain_id": "D", "sequence_id": 1, "insertion_code": null, "label_chain_id": "D", "label_sequence_id": 1, "residue_name": "GLY", "sheet": 12, "strand": 16, "ladders": [31, 30]},
            {"chain_id": "B", "sequence_id": 1, "insertion_code": null, "label_chain_id": "B", "label_sequence_id": 1, "residue_name": "ALA", "sheet": 10, "strand": 15, "ladders": [30, null]}
        ]
    });
    normalize_benchmark_for_comparison_in_place("bio.secondary-structure.dssp", &mut expected);
    normalize_benchmark_for_comparison_in_place("bio.secondary-structure.dssp", &mut actual);
    assert_eq!(expected, actual);
}

#[test]
fn comparison_normalizes_undirected_bonds_and_ring_order() {
    let expected = json!({
        "records": [{
            "bonds": [
                {"index": 0, "begin_atom_index": 5, "end_atom_index": 0, "bond_type": "SINGLE", "stereo": "STEREONONE"}
            ],
            "rings": [[5, 3, 1]]
        }]
    });
    let actual = json!({
        "records": [{
            "bonds": [
                {"index": 7, "begin_atom_index": 0, "end_atom_index": 5, "bond_type": "SINGLE", "stereo": "STEREONONE"}
            ],
            "rings": [[1, 3, 5]]
        }]
    });

    assert_eq!(
        normalize_for_comparison(&expected),
        normalize_for_comparison(&actual)
    );
}

#[test]
fn smiles_semantic_records_assert_topology_and_atom_identity() {
    let single = one_smiles("CC").expect("single bond should parse");
    let double = one_smiles("C=C").expect("double bond should parse");
    assert_ne!(
        smiles_perceived_bonds_json(&single),
        smiles_perceived_bonds_json(&double)
    );

    let aromatic = one_smiles("c1ccccc1").expect("benzene should parse");
    let mut perceived_aromatic = aromatic.clone();
    perceive(&mut perceived_aromatic);
    assert_eq!(
        explicit_valence_json(&perceived_aromatic, AtomId::new(0)),
        3
    );
    let mut aromatic_cyclohexyne = one_smiles("C1=CC#CC=C1").expect("cyclohexyne parses");
    perceive(&mut aromatic_cyclohexyne);
    let alkyne_atoms = aromatic_cyclohexyne
        .bonds()
        .find_map(|(id, bond)| {
            (aromatic_cyclohexyne.bond_is_aromatic(id).ok().flatten() == Some(true)
                && bond.order == BondOrder::Triple)
                .then_some(bond.endpoints())
        })
        .expect("aromaticized triple bond is retained");
    assert_eq!(
        explicit_valence_json(&aromatic_cyclohexyne, alkyne_atoms.0),
        4
    );
    assert_eq!(
        explicit_valence_json(&aromatic_cyclohexyne, alkyne_atoms.1),
        4
    );
    let mut thiophene = one_smiles("c1ccsc1").expect("thiophene parses");
    perceive(&mut thiophene);
    let sulfur_id = thiophene
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "S").then_some(id))
        .expect("sulfur atom");
    assert_eq!(explicit_valence_json(&thiophene, sulfur_id), 2);
    let mut phosphorus_ring = one_smiles("C(F)(F)(F)P1P(P(P(P1C(F)(F)F)C(F)(F)F)C(F)(F)F)C(F)(F)F")
        .expect("phosphorus ring parses");
    perceive(&mut phosphorus_ring);
    for (phosphorus_id, _phosphorus) in phosphorus_ring
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "P")
    {
        assert_eq!(
            phosphorus_ring.atom_is_aromatic(phosphorus_id).unwrap(),
            Some(true)
        );
        assert_eq!(explicit_valence_json(&phosphorus_ring, phosphorus_id), 3);
    }
    let mut phosphinine = one_smiles("C1=CC=PC=C1").expect("phosphinine parses");
    perceive(&mut phosphinine);
    let phosphinine_phosphorus = phosphinine
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "P").then_some(id))
        .expect("phosphinine phosphorus");
    assert_eq!(
        explicit_valence_json(&phosphinine, phosphinine_phosphorus),
        3
    );
    let document = kekule::smiles::parse_str("CN(C)CCO.C1=CC=C2C(=C1)C3=NC4=C5C=CC=CC5=C([N-]4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=N8)N=C2[N-]3.[Cu+2]")
        .expect("anionic macrocycle mixture parses");
    let mut anionic_macrocycle = kekule::smiles::interpret(&document)
        .expect("anionic macrocycle mixture interprets")
        .into_molecules()
        .swap_remove(1);
    perceive(&mut anionic_macrocycle);
    let anionic_nitrogen = anionic_macrocycle
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "N"
                && atom.formal_charge < 0
                && anionic_macrocycle.atom_is_aromatic(id).ok().flatten() == Some(true))
            .then_some(id)
        })
        .expect("anionic aromatic nitrogen");
    assert_eq!(
        explicit_valence_json(&anionic_macrocycle, anionic_nitrogen),
        2
    );
    let mut cyclopentadienyl =
        one_smiles("[CH-]1[C-]=[C-][C-]=[C-]1").expect("cyclopentadienyl anion parses");
    perceive(&mut cyclopentadienyl);
    let anionic_carbon_with_h = cyclopentadienyl
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "C"
                && atom.formal_charge < 0
                && cyclopentadienyl.atom_is_aromatic(id).ok().flatten() == Some(true)
                && atom.hydrogens.explicit_count() > 0)
                .then_some(id)
        })
        .expect("anionic aromatic carbon with explicit hydrogen");
    let anionic_carbon = cyclopentadienyl
        .atom(anionic_carbon_with_h)
        .expect("anionic carbon should exist");
    assert_eq!(
        explicit_valence_json(&cyclopentadienyl, anionic_carbon_with_h)
            + anionic_carbon.hydrogens.explicit_count(),
        3
    );
    let mut substituted_cyclopentadienyl =
        one_smiles("C[C-]1[C-]=[C-][C-]=[C-]1").expect("substituted cyclopentadienyl parses");
    perceive(&mut substituted_cyclopentadienyl);
    let substituted_anionic_carbon = substituted_cyclopentadienyl
        .atoms()
        .find_map(|(id, atom)| {
            let degree = substituted_cyclopentadienyl
                .incident_bonds(id)
                .ok()?
                .count();
            (atom.element.symbol() == "C"
                && atom.formal_charge < 0
                && substituted_cyclopentadienyl
                    .atom_is_aromatic(id)
                    .ok()
                    .flatten()
                    == Some(true)
                && degree == 3)
                .then_some(id)
        })
        .expect("substituted anionic carbon");
    assert_eq!(
        explicit_valence_json(&substituted_cyclopentadienyl, substituted_anionic_carbon,),
        3
    );
    let mut fused_triazine =
        one_smiles("O=[N+]([O-])c2cc(-c1nn5c(=O)c(C=Cc3c(O)ccc4c3cccc4)nnc5s1)ccc2")
            .expect("fused triazine should parse");
    perceive(&mut fused_triazine);
    let tricoordinate_aromatic_nitrogen = fused_triazine
        .atoms()
        .find_map(|(id, atom)| {
            let aromatic_degree = fused_triazine
                .incident_bonds(id)
                .ok()?
                .filter(|(bond, _)| {
                    fused_triazine.bond_is_aromatic(*bond).ok().flatten() == Some(true)
                })
                .count();
            (atom.element.symbol() == "N"
                && fused_triazine.atom_is_aromatic(id).ok().flatten() == Some(true)
                && aromatic_degree >= 3)
                .then_some(id)
        })
        .expect("tri-coordinate aromatic nitrogen");
    assert_eq!(
        explicit_valence_json(&fused_triazine, tricoordinate_aromatic_nitrogen),
        3
    );
    let localized_bonds = smiles_perceived_bonds_json(&aromatic);
    assert_eq!(
        localized_bonds
            .iter()
            .filter(|bond| bond["bond_type"] == "SINGLE" && bond["is_aromatic"] == false)
            .count(),
        3
    );
    assert_eq!(
        localized_bonds
            .iter()
            .filter(|bond| bond["bond_type"] == "DOUBLE" && bond["is_aromatic"] == false)
            .count(),
        3
    );
    assert!(perceived_aromatic
        .bonds()
        .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double)));
    assert!(smiles_perceived_bonds_json(&perceived_aromatic)
        .iter()
        .all(|bond| bond["bond_type"] == "AROMATIC" && bond["is_aromatic"] == true));

    let labeled = one_smiles("[13CH3:7]C").expect("labeled carbon should parse");
    let atoms = smiles_perceived_atoms_json(&labeled);
    assert!(atoms
        .iter()
        .any(|atom| atom["isotope"] == 13 && atom["atom_map"] == 7));
    assert!(atoms.iter().all(|atom| atom["neighbors"].is_array()));
}

#[test]
fn canonical_smiles_records_do_not_prefilter_unsupported_categories() {
    let root = temp_workspace_root("canonical-no-prefilter");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "* CID:example\n").expect("fixture should write");

    let records =
        read_smiles_records(&Input::read(&fixture).unwrap()).expect("records should load");

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].record_index, 0);
    assert_eq!(records[0].status, "parse_error");
    assert_eq!(records[0].input_smiles, "*");
    assert!(records[0].components.is_empty());
}

#[test]
fn canonical_smiles_benchmark_perceives_before_writing() {
    let root = temp_workspace_root("canonical-perceive-before-write");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "C1=CC=CC=C1 CID:benzene\n").expect("fixture should write");

    let records =
        read_smiles_records(&Input::read(&fixture).unwrap()).expect("records should load");
    let item = stereo_smiles_record_json(&records[0], smiles::SmilesWriteMode::Canonical)
        .expect("canonical record should render");

    assert_eq!(item["status"], "ok");
    assert_eq!(item["normalized_perceived"]["atom_count"], 6);
    assert!(item["normalized_perceived"]["atoms"]
        .as_array()
        .unwrap()
        .iter()
        .all(|atom| atom["aromatic"] == true));
}

#[test]
fn canonical_smiles_benchmark_matches_rdkit_parse_status_for_invalid_input() {
    let root = temp_workspace_root("canonical-invalid-input");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "[Cl-](Br)Br CID:invalid\n").expect("fixture should write");

    let records =
        read_smiles_records(&Input::read(&fixture).unwrap()).expect("records should load");
    let item = stereo_smiles_record_json(&records[0], smiles::SmilesWriteMode::Canonical)
        .expect("canonical record should render");

    assert_eq!(item["status"], "parse_error");
}

#[test]
fn smiles_semantics_match_rdkit_aromatic_carbonyl_valence() {
    let molecule =
        one_smiles("CCCCCCCc1cc2c(=O)ccn(O)c2cc1").expect("aromatic carbonyl SMILES should parse");

    let item = smiles_perceived_semantic_json(molecule);
    let atoms = item["atoms"]
        .as_array()
        .expect("perceived atoms should be an array");

    assert!(atoms.iter().any(|atom| {
        atom["symbol"] == "C"
            && atom["aromatic"] == true
            && atom["explicit_valence"] == 4
            && atom["neighbors"].as_array().is_some_and(|neighbors| {
                neighbors.iter().any(|neighbor| {
                    neighbor["bond_type"] == "DOUBLE"
                        && neighbor["atom"]
                            .as_str()
                            .is_some_and(|key| key.starts_with("008|O|0|0|0|0|false|2|"))
                })
            })
    }));
    assert!(!atoms.iter().any(|atom| {
        atom["symbol"] == "C" && atom["aromatic"] == true && atom["explicit_valence"] == 5
    }));
    assert!(atoms.iter().any(|atom| {
        atom["symbol"] == "N"
            && atom["aromatic"] == true
            && atom["explicit_valence"] == 3
            && atom["neighbors"].as_array().is_some_and(|neighbors| {
                neighbors.iter().any(|neighbor| {
                    neighbor["bond_type"] == "SINGLE"
                        && neighbor["atom"]
                            .as_str()
                            .is_some_and(|key| key.starts_with("008|O|0|0|0|1|false|1|"))
                })
            })
    }));
    assert!(!atoms.iter().any(|atom| {
        atom["symbol"] == "N" && atom["aromatic"] == true && atom["explicit_valence"] == 4
    }));
}

#[test]
fn smiles_semantics_match_rdkit_aromatic_nh_no_implicit_flag() {
    let molecule = one_smiles("[nH]1cccc1").expect("aromatic nH SMILES should parse");

    let item = smiles_perceived_semantic_json(molecule);
    let atoms = item["atoms"]
        .as_array()
        .expect("perceived atoms should be an array");

    assert!(atoms.iter().any(|atom| {
        atom["symbol"] == "N"
            && atom["aromatic"] == true
            && atom["explicit_hydrogens"] == 1
            && atom["implicit_hydrogens"] == 0
            && atom["no_implicit_hydrogens"] == false
    }));
}

#[test]
fn smiles_semantics_derive_promoted_aromatic_nh_valence_without_feedback() {
    let molecule = one_smiles("CCOC(=O)C1=C(C(=C(N1)C)C(=O)OC(C)(C)C)C")
        .expect("substituted pyrrole SMILES should parse");

    let item = smiles_perceived_semantic_json(molecule);
    let atoms = item["atoms"]
        .as_array()
        .expect("perceived atoms should be an array");

    assert!(
        atoms.iter().any(|atom| {
            atom["symbol"] == "N"
                && atom["aromatic"] == true
                && atom["explicit_hydrogens"] == 1
                && atom["implicit_hydrogens"] == 0
                && atom["no_implicit_hydrogens"] == false
                && atom["explicit_valence"] == 3
        }),
        "{atoms:#?}"
    );
    assert!(!atoms.iter().any(|atom| {
        atom["symbol"] == "N" && atom["aromatic"] == true && atom["explicit_valence"] == 4
    }));
}

#[test]
fn aromatic_sulfonium_benchmark_uses_trivalent_donor_valence() {
    let mut molecule = one_smiles("C[s+]1cccc1").unwrap();
    molecule.perceive().unwrap();
    let (sulfur, atom) = molecule
        .atoms()
        .find(|(_, atom)| atom.element.symbol() == "S")
        .unwrap();
    assert_eq!(atom.formal_charge, 1);
    assert_eq!(explicit_valence_json(&molecule, sulfur), 3);
    assert_eq!(molecule.implicit_hydrogens(sulfur).unwrap(), Some(0));
}

fn simple_sdf_record(title: &str) -> String {
    format!(
        "{title}
  xtask-test

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"
    )
}

fn chiral_wedge_sdf_record(title: &str) -> String {
    format!(
        "{title}
  xtask-test

  5  4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"
    )
}

fn temp_workspace_root(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time should be available")
        .as_nanos();
    let root = env::temp_dir().join(format!("kekule-xtask-{label}-{}-{nonce}", process::id()));
    fs::create_dir_all(&root).expect("temporary workspace root should create");
    root
}

fn evaluate_path(feature: &str, dataset: &str, path: &Path) -> Result<Value, Box<dyn Error>> {
    evaluate(feature, dataset, &Input::read(path)?)
}

#[test]
fn package_metadata_is_release_consistent() {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let workspace_manifest: toml::Value = toml::from_str(
        &fs::read_to_string(workspace_root.join("Cargo.toml"))
            .expect("workspace manifest should read"),
    )
    .expect("workspace manifest should parse");
    let workspace_package = &workspace_manifest["workspace"]["package"];
    let workspace_version = workspace_package["version"]
        .as_str()
        .expect("workspace package version should be a string");
    assert!(!workspace_version.is_empty());
    assert_eq!(
        workspace_package["repository"].as_str(),
        Some("https://github.com/choutkaj/kekule")
    );

    for (relative_path, package_name, publish) in [
        ("crates/kekule/Cargo.toml", "kekule", true),
        (
            "crates/kekule-potentials/Cargo.toml",
            "kekule-potentials",
            true,
        ),
        ("crates/kekule-traj/Cargo.toml", "kekule-traj", true),
        ("benchmarks/Cargo.toml", "kekule-bench", false),
    ] {
        let manifest: toml::Value = toml::from_str(
            &fs::read_to_string(workspace_root.join(relative_path))
                .expect("package manifest should read"),
        )
        .expect("package manifest should parse");
        assert_eq!(manifest["package"]["name"].as_str(), Some(package_name));
        assert_eq!(
            manifest["package"]["version"]["workspace"].as_bool(),
            Some(true)
        );
        assert_eq!(manifest["package"]["publish"].as_bool(), Some(publish));
    }

    let fuzz_manifest: toml::Value = toml::from_str(
        &fs::read_to_string(workspace_root.join("fuzz/Cargo.toml"))
            .expect("fuzz manifest should read"),
    )
    .expect("fuzz manifest should parse");
    assert_eq!(
        fuzz_manifest["package"]["name"].as_str(),
        Some("kekule-fuzz")
    );
    assert_eq!(fuzz_manifest["package"]["version"].as_str(), Some("0.0.0"));
    assert_eq!(fuzz_manifest["package"]["publish"].as_bool(), Some(false));

    for (relative_path, section, dependency) in [
        ("crates/kekule-traj/Cargo.toml", "dependencies", "kekule"),
        (
            "crates/kekule-potentials/Cargo.toml",
            "dependencies",
            "kekule",
        ),
        (
            "crates/kekule-potentials/Cargo.toml",
            "dev-dependencies",
            "kekule-traj",
        ),
    ] {
        let manifest: toml::Value = toml::from_str(
            &fs::read_to_string(workspace_root.join(relative_path))
                .expect("package manifest should read"),
        )
        .expect("package manifest should parse");
        assert_eq!(
            manifest[section][dependency]["version"].as_str(),
            Some(workspace_version),
            "{relative_path} must require the workspace release version of {dependency}"
        );
    }
}
