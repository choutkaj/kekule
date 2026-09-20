use super::io::IndexedStereoPerceptionRecord;
use crate::boxed_error;
use kekule::{
    core::{AtomId, Molecule, StereoDescriptor, StereoElementKind, StereoGroupKind},
    stereo,
};
use serde_json::{json, Value};
use std::{collections::BTreeMap, error::Error};
fn offset_object_u64(value: &mut Value, key: &str, offset: u64) {
    if let Some(number) = value.get(key).and_then(Value::as_u64) {
        value[key] = json!(number + offset);
    }
}

pub(crate) fn stereo_cip_record_json(
    record: &mut IndexedStereoPerceptionRecord,
) -> Result<Value, Box<dyn Error>> {
    if record.components.is_empty() {
        return Ok(json!({
            "record_index": record.record_index,
            "status": "parse_error",
            "title": record.title,
        }));
    }
    let mut atom_count = 0;
    let mut bond_count = 0;
    let mut atom_descriptors = Vec::new();
    let mut bond_descriptors = Vec::new();
    for molecule in &mut record.components {
        molecule.perceive().map_err(|error| {
            boxed_error(format!(
                "record {} perception failed: {error:?}",
                record.record_index
            ))
        })?;
        stereo::assign_cip_descriptors(molecule).map_err(|error| {
            boxed_error(format!(
                "record {} CIP assignment failed: {error:?}",
                record.record_index
            ))
        })?;
        let atom_index: BTreeMap<_, _> = molecule
            .atom_ids()
            .enumerate()
            .map(|(index, id)| (id, index as u64))
            .collect();
        let mut atoms = cip_atom_descriptors_json(molecule, &atom_index);
        let mut bonds = cip_bond_descriptors_json(molecule, &atom_index);
        for atom in &mut atoms {
            offset_object_u64(atom, "atom_index", atom_count);
        }
        for bond in &mut bonds {
            offset_object_u64(bond, "begin_atom_index", atom_count);
            offset_object_u64(bond, "end_atom_index", atom_count);
        }
        atom_count += atom_index.len() as u64;
        bond_count += molecule.bond_count();
        atom_descriptors.extend(atoms);
        bond_descriptors.extend(bonds);
    }
    Ok(json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": atom_count,
        "bond_count": bond_count,
        "atom_descriptors": atom_descriptors,
        "bond_descriptors": bond_descriptors,
    }))
}

pub(crate) fn cip_atom_descriptors_json(
    mol: &Molecule,
    atom_index: &BTreeMap<AtomId, u64>,
) -> Vec<Value> {
    let mut descriptors = mol
        .stereo_elements()
        .filter_map(|(id, element)| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => mol
                .cip_descriptor(id)
                .expect("live stereo element")
                .and_then(|descriptor| {
                    let atom_index = *atom_index.get(&stereo.center)?;
                    Some(json!({
                        "atom_index": atom_index,
                        "descriptor": stereo_descriptor_json(descriptor),
                    }))
                }),
            StereoElementKind::Axis(_) | StereoElementKind::DoubleBond(_) => None,
        })
        .collect::<Vec<_>>();
    descriptors.sort_by_key(|value| {
        value
            .get("atom_index")
            .and_then(Value::as_u64)
            .unwrap_or(u64::MAX)
    });
    descriptors
}

pub(crate) fn cip_bond_descriptors_json(
    mol: &Molecule,
    atom_index: &BTreeMap<AtomId, u64>,
) -> Vec<Value> {
    let mut descriptors = mol
        .stereo_elements()
        .filter_map(|(id, element)| match &element.kind {
            StereoElementKind::DoubleBond(stereo) => mol
                .cip_descriptor(id)
                .expect("live stereo element")
                .and_then(|descriptor| {
                    let begin_atom_index = *atom_index.get(&stereo.left)?;
                    let end_atom_index = *atom_index.get(&stereo.right)?;
                    Some(json!({
                        "begin_atom_index": begin_atom_index,
                        "end_atom_index": end_atom_index,
                        "descriptor": stereo_descriptor_json(descriptor),
                    }))
                }),
            StereoElementKind::Axis(stereo) => mol
                .cip_descriptor(id)
                .expect("live stereo element")
                .and_then(|descriptor| {
                    let bond = mol.bond(stereo.axis).expect("live stereo axis");
                    let (begin, end) = bond.endpoints();
                    let begin_atom_index = *atom_index.get(&begin)?;
                    let end_atom_index = *atom_index.get(&end)?;
                    Some(json!({
                        "begin_atom_index": begin_atom_index,
                        "end_atom_index": end_atom_index,
                        "descriptor": stereo_descriptor_json(descriptor),
                    }))
                }),
            StereoElementKind::Tetrahedral(_) => None,
        })
        .collect::<Vec<_>>();
    descriptors.sort_by(|left, right| {
        let left_key = (
            left.get("begin_atom_index")
                .and_then(Value::as_u64)
                .unwrap_or(u64::MAX),
            left.get("end_atom_index")
                .and_then(Value::as_u64)
                .unwrap_or(u64::MAX),
        );
        let right_key = (
            right
                .get("begin_atom_index")
                .and_then(Value::as_u64)
                .unwrap_or(u64::MAX),
            right
                .get("end_atom_index")
                .and_then(Value::as_u64)
                .unwrap_or(u64::MAX),
        );
        left_key.cmp(&right_key).then_with(|| {
            left.get("descriptor")
                .and_then(Value::as_str)
                .unwrap_or("")
                .cmp(
                    right
                        .get("descriptor")
                        .and_then(Value::as_str)
                        .unwrap_or(""),
                )
        })
    });
    descriptors
}

pub(crate) fn stereo_descriptor_json(descriptor: StereoDescriptor) -> &'static str {
    match descriptor {
        StereoDescriptor::R => "R",
        StereoDescriptor::S => "S",
        StereoDescriptor::LowerR => "r",
        StereoDescriptor::LowerS => "s",
        // RDKit's CIPLabeler renders these pseudoasymmetric bond descriptors
        // as lowercase e/z; uppercase E/Z remain distinct descriptors.
        StereoDescriptor::SeqTrans => "e",
        StereoDescriptor::SeqCis => "z",
        StereoDescriptor::E => "E",
        StereoDescriptor::Z => "Z",
        StereoDescriptor::M => "M",
        StereoDescriptor::P => "P",
        StereoDescriptor::LowerM => "m",
        StereoDescriptor::LowerP => "p",
    }
}

pub(crate) fn stereo_group_kind_json(kind: StereoGroupKind) -> &'static str {
    match kind {
        StereoGroupKind::Absolute => "absolute",
        StereoGroupKind::Relative => "relative",
        StereoGroupKind::Racemic => "racemic",
        StereoGroupKind::And => "and",
        StereoGroupKind::Or => "or",
    }
}
