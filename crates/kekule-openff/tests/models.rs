//! Public-API tests against independently generated OpenFF observations.
use kekule::{core::Molecule, smiles};
use kekule_openff::{ChargeSource, ForceField, NaglModel, ParameterizedTopology};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{collections::BTreeMap, path::PathBuf};

const ROSEMARY: &str = include_str!("../data/rosemary.offxml");
fn reference() -> Value {
    serde_json::from_reader(flate2::read::GzDecoder::new(
        &include_bytes!("fixtures/models.json.gz")[..],
    ))
    .unwrap()
}
fn xml(model: &Value) -> String {
    let preset = ForceField::rosemary().unwrap();
    ROSEMARY
        .replace(
            preset.charge_model().unwrap().model_file(),
            model["name"].as_str().unwrap(),
        )
        .replace(
            preset.charge_model().unwrap().checkpoint_sha256(),
            model["checkpoint_sha256"].as_str().unwrap(),
        )
}
fn bundles() -> PathBuf {
    std::env::var_os("KEKULE_OPENFF_MODELS")
        .expect("set KEKULE_OPENFF_MODELS to the two-model export directory")
        .into()
}
fn close(actual: f64, expected: f64, atol: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= atol + 1e-12 * expected.abs(),
        "{actual} != {expected}"
    );
}
fn charges(actual: &[f64], expected: &Value, molecule: &Molecule, atol: f64) {
    assert_eq!(actual.len(), expected.as_array().unwrap().len());
    for (i, (_, atom)) in molecule.atoms().enumerate() {
        close(
            actual[i],
            expected[atom.atom_map.unwrap() as usize - 1]
                .as_f64()
                .unwrap(),
            atol,
        );
    }
    close(actual.iter().sum(), molecule.formal_charge() as f64, 1e-10);
}

// Compare every numeric term and its atom-map association, preserving multiplicity.
fn parameters(actual: &ParameterizedTopology, expected: &Value, molecule: &Molecule) {
    type Rows = BTreeMap<Vec<u32>, Vec<BTreeMap<String, f64>>>;
    let maps = |atoms: &[kekule::topology::InstanceAtomId]| {
        atoms
            .iter()
            .map(|a| molecule.atom(a.atom()).unwrap().atom_map.unwrap())
            .collect::<Vec<_>>()
    };
    let mut observed: BTreeMap<&str, Vec<(Vec<u32>, Value)>> = BTreeMap::new();
    observed.insert(
        "Bonds",
        actual
            .bonds()
            .iter()
            .map(|t| {
                (
                    maps(&t.atoms),
                    json!({"k":t.parameter.k.value(), "length":t.parameter.length.value()}),
                )
            })
            .collect(),
    );
    observed.insert(
        "Angles",
        actual
            .angles()
            .iter()
            .map(|t| {
                (
                    maps(&t.atoms),
                    json!({"k":t.parameter.k.value(), "angle":t.parameter.angle.value()}),
                )
            })
            .collect(),
    );
    observed.insert(
        "Constraints",
        actual
            .constraints()
            .iter()
            .map(|t| {
                (
                    maps(&t.atoms),
                    json!({"distance":t.parameter.distance.value()}),
                )
            })
            .collect(),
    );
    observed.insert(
        "vdW",
        molecule
            .atoms()
            .zip(actual.vdw())
            .map(|((_, a), p)| {
                (
                    vec![a.atom_map.unwrap()],
                    json!({"sigma":p.sigma.value(), "epsilon":p.epsilon.value()}),
                )
            })
            .collect(),
    );
    for (name, terms) in [
        ("ProperTorsions", actual.proper_torsions()),
        ("ImproperTorsions", actual.improper_torsions()),
    ] {
        observed.insert(name, terms.iter().flat_map(|t| t.parameter.terms.iter().map(|p| (maps(&t.atoms), json!({"k":p.k.value(),"phase":p.phase.value(),"periodicity":p.periodicity,"idivf":p.idivf})))).collect());
    }
    for (handler, rows) in observed {
        let canonical = |mut atoms: Vec<u32>| {
            if handler == "ImproperTorsions" {
                atoms[1..].sort();
                atoms
            } else {
                let reverse = atoms.iter().rev().copied().collect();
                atoms.min(reverse)
            }
        };
        let collect = |rows: Vec<(Vec<u32>, Value)>| {
            let mut result = Rows::new();
            for (atoms, p) in rows {
                result
                    .entry(canonical(atoms))
                    .or_default()
                    .push(serde_json::from_value(p).unwrap());
            }
            for terms in result.values_mut() {
                terms.sort_by(|a, b| {
                    a.values()
                        .zip(b.values())
                        .map(|(a, b)| a.total_cmp(b))
                        .find(|o| !o.is_eq())
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            }
            result
        };
        let actual = collect(rows);
        let expected = collect(
            expected[handler]
                .as_array()
                .unwrap()
                .iter()
                .map(|r| {
                    (
                        serde_json::from_value(r["maps"].clone()).unwrap(),
                        r["parameters"].clone(),
                    )
                })
                .collect(),
        );
        assert_eq!(
            actual.keys().collect::<Vec<_>>(),
            expected.keys().collect::<Vec<_>>(),
            "{handler}"
        );
        for (key, terms) in actual {
            assert_eq!(terms.len(), expected[&key].len(), "{handler}: {key:?}");
            for (a, b) in terms.iter().zip(&expected[&key]) {
                assert_eq!(a.keys().collect::<Vec<_>>(), b.keys().collect::<Vec<_>>());
                for (name, &value) in a {
                    close(value, b[name], 1e-10);
                }
            }
        }
    }
}

#[test]
fn model_declarations_survive_offxml_loading() {
    let report = reference();
    let models = report["models"].as_array().unwrap();
    assert_eq!(models.len(), 2);
    for model in models {
        let ff = ForceField::from_offxml(&xml(model)).unwrap();
        assert_eq!(
            ff.charge_model().unwrap().model_file(),
            model["name"].as_str().unwrap()
        );
        assert_eq!(
            ff.charge_model().unwrap().checkpoint_sha256(),
            model["checkpoint_sha256"].as_str().unwrap()
        );
    }
    assert_eq!(models[0]["lookup_entries"], 13944);
    assert_eq!(models[1]["lookup_entries"], 0);
}

#[test]
#[ignore = "requires externally exported models; set KEKULE_OPENFF_MODELS (see VALIDATION.md)"]
fn both_models_reproduce_complete_openff_parameterization() {
    let report = reference();
    let mut counts = [0usize; 3];
    for reference in report["models"].as_array().unwrap() {
        let name = reference["name"].as_str().unwrap();
        let directory = bundles().join(name.trim_end_matches(".pt"));
        for (file, expected) in reference["bundle_sha256"].as_object().unwrap() {
            assert_eq!(
                format!(
                    "{:x}",
                    Sha256::digest(std::fs::read(directory.join(file)).unwrap())
                ),
                expected.as_str().unwrap()
            );
        }
        let model = NaglModel::load(&directory).unwrap();
        let ff = ForceField::from_file(directory.join("force-field.offxml")).unwrap();
        assert_eq!(ff.charge_model(), Some(model.identity()));
        assert_eq!(
            model.lookup_entry_count(),
            reference["lookup_entries"].as_u64().unwrap() as usize
        );
        for r in reference["records"].as_array().unwrap() {
            let molecule = smiles::to_molecules(r["smiles"].as_str().unwrap())
                .unwrap()
                .remove(0);
            let before = molecule.clone();
            let p = ff
                .parameterize_molecule(&molecule, &model)
                .unwrap_or_else(|e| panic!("{name} {}: {e}", r["id"]));
            charges(
                p.charges().value(),
                &r["system"]["charges"],
                &molecule,
                5e-5,
            );
            parameters(&p, &r["system"]["parameters"], &molecule);
            match &p.charge_sources()[0] {
                ChargeSource::Library { .. } => counts[0] += 1,
                ChargeSource::Lookup {
                    model: identity, ..
                } => {
                    assert_eq!(identity, model.identity());
                    counts[1] += 1;
                }
                ChargeSource::Inference { model: identity } => {
                    assert_eq!(identity, model.identity());
                    counts[2] += 1;
                }
            }
            if r["inference"]["status"] == "ok" {
                charges(
                    model.infer_charges(&molecule).unwrap().charges.value(),
                    &r["inference"]["charges"],
                    &molecule,
                    1e-6,
                );
                let f = model.atom_features(&molecule).unwrap();
                for (i, (_, atom)) in molecule.atoms().enumerate() {
                    let expected = r["inference"]["features"][atom.atom_map.unwrap() as usize - 1]
                        .as_array()
                        .unwrap();
                    assert_eq!(f[i].len(), expected.len());
                    for (&a, b) in f[i].iter().zip(expected) {
                        close(f64::from(a), b.as_f64().unwrap(), 1e-6);
                    }
                }
            } else {
                assert!(model.infer_charges(&molecule).is_err());
                assert!(model.atom_features(&molecule).is_err());
            }
            if r["assignment"]["status"] == "ok" {
                charges(
                    model.assign_charges(&molecule).unwrap().charges.value(),
                    &r["assignment"]["charges"],
                    &molecule,
                    5e-5,
                );
            } else {
                assert!(model.assign_charges(&molecule).is_err());
            }
            assert_eq!(molecule, before);
        }
        assert_eq!(reference["records"].as_array().unwrap().len(), 66);
        if name != "openff-gnn-am1bcc-1.0.0.pt" {
            let molecule =
                smiles::to_molecules(reference["records"][0]["smiles"].as_str().unwrap())
                    .unwrap()
                    .remove(0);
            assert!(ForceField::rosemary()
                .unwrap()
                .parameterize_molecule(&molecule, &model)
                .unwrap_err()
                .to_string()
                .contains("identity mismatch"));
        }
    }
    assert!(
        counts.iter().all(|&n| n > 0),
        "library, lookup and inference paths must all be exercised"
    );
}

struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let path =
            std::env::temp_dir().join(format!("kekule-nagl-contract-{}", std::process::id()));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
#[ignore = "requires externally exported models; set KEKULE_OPENFF_MODELS (see VALIDATION.md)"]
fn bundle_validation_and_model_binding_fail_before_parameterization() {
    let directory = bundles().join("openff-gnn-am1bcc-0.1.0-rc.2");
    let scratch = Scratch::new();
    let manifest: Value =
        serde_json::from_slice(&std::fs::read(directory.join("model.json")).unwrap()).unwrap();
    std::fs::copy(directory.join("weights.bin"), scratch.0.join("weights.bin")).unwrap();
    let cases = [
        ("/schema", json!(99)),
        ("/preparation", json!("unknown")),
        ("/checkpoint_sha256", json!("bad")),
        ("/weights_sha256", json!("0".repeat(64))),
        ("/config/convolution/architecture", json!("GATConv")),
        (
            "/config/convolution/layers/0/aggregator_type",
            json!("pool"),
        ),
        (
            "/config/convolution/layers/0/hidden_feature_size",
            json!(4096),
        ),
        ("/config/atom_features/0/name", json!("unknown_feature")),
        ("/config/atom_features/0/categories", json!(["C", "C"])),
        (
            "/config/readouts/am1bcc_charges/postprocess",
            json!("compute_partial_charges"),
        ),
        (
            "/tensors/convolution_module.gcn_layers.0.fc_self.weight/shape",
            json!([512, 22]),
        ),
        (
            "/tensors/convolution_module.gcn_layers.0.fc_self.weight/offset",
            json!(usize::MAX),
        ),
        ("/lookup_tables", json!({"wrong_readout":[]})),
    ];
    for (pointer, value) in cases {
        let mut changed = manifest.clone();
        *changed.pointer_mut(pointer).unwrap() = value;
        std::fs::write(
            scratch.0.join("model.json"),
            serde_json::to_vec(&changed).unwrap(),
        )
        .unwrap();
        assert!(NaglModel::load(&scratch.0).is_err(), "{pointer}");
    }
    std::fs::write(
        scratch.0.join("model.json"),
        serde_json::to_vec(&manifest).unwrap(),
    )
    .unwrap();
    let model = NaglModel::load(&scratch.0).unwrap();
    let ff = ForceField::from_file(directory.join("force-field.offxml")).unwrap();
    let molecule = smiles::to_molecules(
        reference()["models"][1]["records"][0]["smiles"]
            .as_str()
            .unwrap(),
    )
    .unwrap()
    .remove(0);
    let expected = ff.parameterize_molecule(&molecule, &model).unwrap();
    for changed in [
        xml(&reference()["models"][1]).replace(model.identity().model_file(), "different.pt"),
        xml(&reference()["models"][1])
            .replace(model.identity().checkpoint_sha256(), &"0".repeat(64)),
    ] {
        assert!(ForceField::from_offxml(&changed)
            .unwrap()
            .parameterize_molecule(&molecule, &model)
            .unwrap_err()
            .to_string()
            .contains("identity mismatch"));
    }
    assert_eq!(
        ff.parameterize_molecule(&molecule, &model)
            .unwrap()
            .charges(),
        expected.charges()
    );
    // Reorder an external model's element vocabulary and matching tensor columns.
    // This represents the same learned function and catches implicit feature ordering.
    let mut permuted = manifest.clone();
    permuted["config"]["atom_features"][0]["categories"]
        .as_array_mut()
        .unwrap()
        .reverse();
    let mut bytes = std::fs::read(directory.join("weights.bin")).unwrap();
    for suffix in ["fc_self.weight", "fc_neigh.weight"] {
        let tensor = &manifest["tensors"][format!("convolution_module.gcn_layers.0.{suffix}")];
        let offset = tensor["offset"].as_u64().unwrap() as usize * 4;
        let rows = tensor["shape"][0].as_u64().unwrap() as usize;
        let columns = tensor["shape"][1].as_u64().unwrap() as usize;
        for row in 0..rows {
            let start = offset + row * columns * 4;
            let original = bytes[start..start + 40].to_vec();
            for column in 0..10 {
                bytes[start + 4 * column..start + 4 * column + 4]
                    .copy_from_slice(&original[4 * (9 - column)..4 * (9 - column) + 4]);
            }
        }
    }
    permuted["weights_sha256"] = json!(format!("{:x}", Sha256::digest(&bytes)));
    std::fs::write(scratch.0.join("weights.bin"), bytes).unwrap();
    std::fs::write(
        scratch.0.join("model.json"),
        serde_json::to_vec(&permuted).unwrap(),
    )
    .unwrap();
    let permuted = NaglModel::load(&scratch.0).unwrap();
    let p = ff.parameterize_molecule(&molecule, &permuted).unwrap();
    for (&a, &b) in p.charges().value().iter().zip(expected.charges().value()) {
        close(a, b, 1e-6);
    }
    // A malformed lookup must fail, rather than silently discard surplus charges.
    std::fs::copy(directory.join("weights.bin"), scratch.0.join("weights.bin")).unwrap();
    let ash: Value = serde_json::from_slice(
        &std::fs::read(bundles().join("openff-gnn-am1bcc-1.0.0/model.json")).unwrap(),
    )
    .unwrap();
    let report = reference();
    let case = report["models"][0]["records"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["id"] == "methane")
        .unwrap();
    let methane = smiles::to_molecules(case["smiles"].as_str().unwrap())
        .unwrap()
        .remove(0);
    let key = model.lookup_identifier(&methane).unwrap();
    let mut entry = ash["lookup_tables"]["am1bcc_charges"]
        .as_array()
        .unwrap()
        .iter()
        .find(|e| e["inchi"] == key)
        .unwrap()
        .clone();
    entry["charges"].as_array_mut().unwrap().push(json!(0.0));
    let mut changed = manifest.clone();
    changed["lookup_tables"] = json!({"am1bcc_charges": [entry]});
    std::fs::write(
        scratch.0.join("model.json"),
        serde_json::to_vec(&changed).unwrap(),
    )
    .unwrap();
    let broken_lookup = NaglModel::load(&scratch.0).unwrap();
    assert!(ff
        .parameterize_molecule(&methane, &broken_lookup)
        .unwrap_err()
        .to_string()
        .contains("lookup charge array"));
    // Match the checksum so truncated float32 storage reaches the decoder.
    let truncated = [0u8; 3];
    changed["weights_sha256"] = json!(format!("{:x}", Sha256::digest(truncated)));
    std::fs::write(
        scratch.0.join("model.json"),
        serde_json::to_vec(&changed).unwrap(),
    )
    .unwrap();
    std::fs::write(scratch.0.join("weights.bin"), truncated).unwrap();
    assert!(NaglModel::load(&scratch.0)
        .unwrap_err()
        .to_string()
        .contains("weights must be little-endian float32"));
    std::fs::remove_file(scratch.0.join("weights.bin")).unwrap();
    assert!(NaglModel::load(&scratch.0).is_err());
    std::fs::File::create(scratch.0.join("model.json"))
        .unwrap()
        .set_len(64 * 1024 * 1024 + 1)
        .unwrap();
    assert!(NaglModel::load(&scratch.0)
        .unwrap_err()
        .to_string()
        .contains("file size"));
}
