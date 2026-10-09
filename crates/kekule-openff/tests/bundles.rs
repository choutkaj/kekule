//! Bundle loading and validation against the bundled Ash model, so the
//! loader's rejection paths run without externally exported files.
#![cfg(feature = "ash")]

use std::io::Read;
use std::path::PathBuf;

use kekule_openff::{ErrorKind, NaglModel};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!("kekule-ash-bundle-{}", std::process::id()));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn write(&self, manifest: &Value, weights: &[u8]) {
        std::fs::write(
            self.0.join("model.json"),
            serde_json::to_vec(manifest).unwrap(),
        )
        .unwrap();
        std::fs::write(self.0.join("weights.bin"), weights).unwrap();
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The bundled weights decoded to the exported `weights.bin` layout.
fn ash_weights() -> Vec<u8> {
    let mut planes = Vec::new();
    flate2::read::ZlibDecoder::new(kekule_openff_ash::WEIGHT_PLANES_ZLIB)
        .read_to_end(&mut planes)
        .unwrap();
    let count = planes.len() / 4;
    let mut weights = vec![0; planes.len()];
    for (plane, bytes) in planes.chunks_exact(count).enumerate() {
        for (value, byte) in bytes.iter().enumerate() {
            weights[value * 4 + plane] = *byte;
        }
    }
    weights
}

#[test]
fn exported_ash_layout_loads_like_the_bundled_model_and_rejects_corruption() {
    let manifest: Value = serde_json::from_slice(kekule_openff_ash::MANIFEST).unwrap();
    let weights = ash_weights();
    assert_eq!(
        format!("{:x}", Sha256::digest(&weights)),
        manifest["weights_sha256"].as_str().unwrap()
    );
    let scratch = Scratch::new();
    scratch.write(&manifest, &weights);
    let loaded = NaglModel::load(&scratch.0).unwrap();
    let bundled = NaglModel::ash().unwrap();
    assert_eq!(loaded.identity(), bundled.identity());
    assert_eq!(
        kekule_openff::diagnostics::lookup_entry_count(&loaded),
        kekule_openff::diagnostics::lookup_entry_count(&bundled)
    );
    let mut ethanol = kekule::smiles::to_molecules("CCO").unwrap().remove(0);
    ethanol.perceive().unwrap();
    ethanol.add_hydrogens().unwrap();
    assert_eq!(
        kekule_openff::diagnostics::infer_charges(&loaded, &ethanol).unwrap(),
        kekule_openff::diagnostics::infer_charges(&bundled, &ethanol).unwrap()
    );

    let cases = [
        ("/schema", json!(1)),
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
            json!([512, 21]),
        ),
        (
            "/tensors/convolution_module.gcn_layers.0.fc_self.weight/offset",
            json!(usize::MAX),
        ),
        ("/lookup_tables", json!({"wrong_readout": []})),
        ("/domain/allowed_elements", json!([0])),
    ];
    for (pointer, value) in cases {
        let mut changed = manifest.clone();
        *changed.pointer_mut(pointer).unwrap() = value;
        scratch.write(&changed, &weights);
        let error = NaglModel::load(&scratch.0).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::Model, "{pointer}: {error}");
        assert_eq!(error.path(), Some(scratch.0.as_path()), "{pointer}");
    }

    // Truncated weights fail on their checksum, before any tensor is read.
    scratch.write(&manifest, &weights[..weights.len() - 4]);
    let error = NaglModel::load(&scratch.0).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::Model);
    assert!(error.to_string().contains("checksum"), "{error}");

    // A lookup entry whose charges do not match its atoms fails on a hit.
    let mut changed = manifest.clone();
    let entries = changed["lookup_tables"]["am1bcc_charges"]
        .as_array_mut()
        .unwrap();
    let methane = kekule_openff::diagnostics::lookup_key(&bundled, &{
        let mut methane = kekule::smiles::to_molecules("C").unwrap().remove(0);
        methane.perceive().unwrap();
        methane.add_hydrogens().unwrap();
        methane
    })
    .unwrap()
    .unwrap();
    let entry = entries
        .iter_mut()
        .find(|entry| entry["inchi"] == methane)
        .unwrap();
    entry["charges"].as_array_mut().unwrap().push(json!(0.0));
    scratch.write(&changed, &weights);
    let broken = NaglModel::load(&scratch.0).unwrap();
    let mut molecule = kekule::smiles::to_molecules("C").unwrap().remove(0);
    molecule.perceive().unwrap();
    molecule.add_hydrogens().unwrap();
    let error = broken.assign_charges(&molecule).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::Model);

    // Lookup entries have no identifier size limit; selection computes none.
    let mut changed = manifest.clone();
    changed["lookup_tables"]["am1bcc_charges"][0]["charges"] = json!(vec![0.0; 5000]);
    scratch.write(&changed, &weights);
    NaglModel::load(&scratch.0).unwrap();
}
