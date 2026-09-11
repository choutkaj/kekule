//! JSON-lines adapter for optional external Molfile stereo interchange checks.
//! Build with `cargo build -p xtask --example molfile_interchange_probe --locked`.

use std::io::{self, BufRead};

use kekule::{core::StereoElementKind, molfile, stereo, structure::Model};
use serde_json::{json, Value};

fn probe(input: Value) -> Result<Value, String> {
    let source = input["molfile"].as_str().ok_or("missing molfile")?;
    let mut model = molfile::parse_str(source)
        .map_err(|e| e.to_string())?
        .to_model()
        .map_err(|e| e.to_string())?;
    if input["clear_double_stereo"].as_bool() == Some(true) {
        let molecules = model.topology().molecules().collect::<Vec<_>>();
        if molecules.len() != 1 {
            return Err("clearing stereo requires one component".into());
        }
        let mut editor = molecules[0].molecule().edit();
        let ids = editor
            .stereo_elements()
            .filter_map(|(id, element)| {
                matches!(element.kind, StereoElementKind::DoubleBond(_)).then_some(id)
            })
            .collect::<Vec<_>>();
        for id in ids {
            editor
                .remove_stereo_element(id)
                .map_err(|e| e.to_string())?;
        }
        let molecule = editor.finish().map_err(|e| e.to_string())?;
        model = Model::from_molecule(&molecule, model.positions()).map_err(|e| e.to_string())?;
    }
    let mut labels = Vec::new();
    let mut atom_offset = 0;
    let mut bond_offset = 0;
    for occurrence in model.topology().molecules() {
        let mut molecule = occurrence.molecule().clone();
        molecule.perceive().map_err(|e| e.to_string())?;
        let assigned = stereo::assign_cip_descriptors(&mut molecule).map_err(|e| e.to_string())?;
        for assignment in assigned.assigned {
            let element = molecule
                .stereo_element(assignment.element)
                .map_err(|e| e.to_string())?;
            let (focus, index) = match &element.kind {
                StereoElementKind::Tetrahedral(stereo) => {
                    ("atom", atom_offset + stereo.center.index())
                }
                StereoElementKind::DoubleBond(stereo) => {
                    ("bond", bond_offset + stereo.bond.index())
                }
                StereoElementKind::Axis(stereo) => ("bond", bond_offset + stereo.axis.index()),
            };
            labels.push((focus, index, format!("{:?}", assignment.descriptor)));
        }
        atom_offset += molecule.atom_count();
        bond_offset += molecule.bond_count();
    }
    labels.sort();
    Ok(json!({
        "status": "ok", "labels": labels,
        "v2000": molfile::write_model_v2000(&model).map_err(|e| e.to_string()),
        "v3000": molfile::write_model_v3000(&model).map_err(|e| e.to_string()),
    }))
}

fn main() {
    for line in io::stdin().lock().lines() {
        let result = line
            .map_err(|e| e.to_string())
            .and_then(|line| serde_json::from_str(&line).map_err(|e| e.to_string()))
            .and_then(probe);
        println!(
            "{}",
            result.unwrap_or_else(|error| json!({"status": "error", "message": error}))
        );
    }
}
