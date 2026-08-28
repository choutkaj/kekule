use crate::*;

use super::chemistry::*;
use super::io::{interpret_smiles, IndexedSmilesRecord};
use super::stereo::*;

pub(crate) fn smiles_write_record_json(
    record: &IndexedSmilesRecord,
) -> Result<Value, Box<dyn Error>> {
    if record.components.is_empty() {
        let mut item = smiles_error_record_json(record);
        if record.status == "ok" {
            item["status"] = json!("write_error");
        }
        return Ok(item);
    }
    let written = record
        .components
        .iter()
        .map(smiles::write)
        .collect::<Result<Vec<_>, _>>()?;
    let reparsed = written
        .iter()
        .map(|text| interpret_smiles(text))
        .collect::<Result<Vec<_>, _>>();
    let Ok(reparsed) = reparsed else {
        return Ok(json!({
            "record_index": record.record_index,
            "status": "write_reparse_error",
            "title": record.title,
            "input_smiles": record.input_smiles,
        }));
    };
    let normalized_perceived = if reparsed.len() == 1 {
        smiles_perceived_semantic_json(reparsed.into_iter().next().expect("one reparsed component"))
    } else {
        smiles_components_perceived_semantic_json(&reparsed)
    };
    Ok(json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "input_smiles": record.input_smiles,
        "normalized_perceived": normalized_perceived,
    }))
}

pub(crate) fn canonical_smiles_record_json(
    record: &IndexedSmilesRecord,
    exact_smiles: bool,
) -> Result<Value, Box<dyn Error>> {
    let Some(molecule) = &record.molecule else {
        return Ok(smiles_error_record_json(record));
    };
    let mut molecule = molecule.clone();
    if molecule.perceive().is_err() {
        return Ok(json!({
            "record_index": record.record_index,
            "status": "parse_error",
            "title": record.title,
            "input_smiles": record.input_smiles,
        }));
    }
    let written = smiles::write_canonical(&molecule)?;
    let reparsed = match interpret_smiles(&written) {
        Ok(reparsed) => reparsed,
        Err(_) => {
            return Ok(json!({
                "record_index": record.record_index,
                "status": "write_reparse_error",
                "title": record.title,
                "input_smiles": record.input_smiles,
                "canonical_smiles": written,
            }));
        }
    };
    let mut item = json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "input_smiles": record.input_smiles,
        "normalized_perceived": smiles_perceived_semantic_json(reparsed),
    });
    if exact_smiles {
        item["canonical_smiles"] = json!(written);
    }
    Ok(item)
}

pub(crate) fn isomeric_smiles_record_json(
    record: &IndexedSmilesRecord,
) -> Result<Value, Box<dyn Error>> {
    let Some(molecule) = &record.molecule else {
        return Ok(smiles_error_record_json(record));
    };
    let mut molecule = molecule.clone();
    if molecule.perceive().is_err() {
        return Ok(json!({
            "record_index": record.record_index,
            "status": "perception_error",
            "title": record.title,
            "input_smiles": record.input_smiles,
        }));
    }
    let written = match smiles::write_isomeric(&molecule) {
        Ok(written) => written,
        Err(error) => {
            return Ok(json!({
                "record_index": record.record_index,
                "status": "write_error",
                "title": record.title,
                "input_smiles": record.input_smiles,
                "message": error.message(),
            }));
        }
    };
    let reparsed = match interpret_smiles(&written) {
        Ok(reparsed) => reparsed,
        Err(_) => {
            return Ok(json!({
                "record_index": record.record_index,
                "status": "write_reparse_error",
                "title": record.title,
                "input_smiles": record.input_smiles,
            }));
        }
    };
    Ok(json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "input_smiles": record.input_smiles,
        "normalized_perceived": smiles_perceived_semantic_json(reparsed.clone()),
        "stereo": smiles_isomeric_stereo_semantic_json(reparsed),
    }))
}

pub(crate) fn isomeric_smiles_record_is_stereo_bearing(record: &IndexedSmilesRecord) -> bool {
    if !record.input_smiles.contains('@')
        && !record.input_smiles.contains('/')
        && !record.input_smiles.contains('\\')
    {
        return false;
    }
    let Some(molecule) = &record.molecule else {
        return false;
    };
    let mut molecule = molecule.clone();
    molecule.perceive().is_ok()
}

pub(crate) fn smiles_parse_record_json(record: &IndexedSmilesRecord) -> Value {
    if record.components.is_empty() {
        return smiles_error_record_json(record);
    }
    let reparsed = record
        .components
        .iter()
        .map(|molecule| {
            smiles::write(molecule)
                .map_err(|_| ())
                .and_then(|text| interpret_smiles(&text).map_err(|_| ()))
        })
        .collect::<Result<Vec<_>, _>>();
    let round_trip = match reparsed {
        Ok(reparsed) => smiles_components_perceived_semantic_json(&reparsed),
        Err(_) => json!({ "status": "write_reparse_error" }),
    };
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "input_smiles": record.input_smiles,
        "raw": smiles_components_raw_semantic_json(&record.components),
        "normalized_perceived": smiles_components_perceived_semantic_json(&record.components),
        "write_round_trip": round_trip,
    })
}

pub(crate) fn smiles_error_record_json(record: &IndexedSmilesRecord) -> Value {
    json!({
        "record_index": record.record_index,
        "status": record.status,
        "title": record.title,
        "input_smiles": record.input_smiles,
    })
}

pub(crate) fn smiles_raw_semantic_json(molecule: &Molecule) -> Value {
    let mol = &molecule;
    json!({
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": basic_atoms_json(mol),
        "bonds": basic_bonds_json(mol),
    })
}

pub(crate) fn smiles_components_raw_semantic_json(components: &[Molecule]) -> Value {
    if let [molecule] = components {
        return smiles_raw_semantic_json(molecule);
    }
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut atom_offset = 0u64;
    let mut bond_offset = 0u64;
    for component in components {
        let mol = component;
        for mut atom in basic_atoms_json(mol) {
            offset_object_u64(&mut atom, "index", atom_offset);
            atoms.push(atom);
        }
        for mut bond in basic_bonds_json(mol) {
            offset_object_u64(&mut bond, "index", bond_offset);
            offset_object_u64(&mut bond, "begin_atom_index", atom_offset);
            offset_object_u64(&mut bond, "end_atom_index", atom_offset);
            bonds.push(bond);
        }
        atom_offset += mol.atom_count() as u64;
        bond_offset += mol.bond_count() as u64;
    }
    json!({
        "atom_count": atom_offset,
        "bond_count": bond_offset,
        "atoms": atoms,
        "bonds": bonds,
    })
}

pub(super) fn offset_object_u64(value: &mut Value, key: &str, offset: u64) {
    if offset == 0 {
        return;
    }
    let Some(number) = value
        .as_object_mut()
        .and_then(|object| object.get_mut(key))
        .and_then(|value| value.as_u64())
    else {
        return;
    };
    value
        .as_object_mut()
        .expect("checked object")
        .insert(key.to_owned(), json!(number + offset));
}

pub(crate) fn smiles_perceived_semantic_json(mut molecule: Molecule) -> Value {
    if molecule.perceive().is_err() {
        return json!({ "status": "perception_error" });
    }
    let mol = &molecule;
    json!({
        "status": "ok",
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": smiles_perceived_atoms_json(mol),
        "bonds": smiles_perceived_bonds_json(mol),
    })
}

pub(crate) fn smiles_components_perceived_semantic_json(components: &[Molecule]) -> Value {
    let mut molecules = components.to_vec();
    if molecules
        .iter_mut()
        .any(|molecule| molecule.perceive().is_err())
    {
        return json!({ "status": "perception_error" });
    }
    let atom_count = molecules
        .iter()
        .map(|molecule| molecule.atom_count())
        .sum::<usize>();
    let bond_count = molecules
        .iter()
        .map(|molecule| molecule.bond_count())
        .sum::<usize>();
    let mut atoms = molecules
        .iter()
        .flat_map(smiles_perceived_atom_entries_json)
        .collect::<Vec<_>>();
    atoms.sort_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then_with(|| left.1.to_string().cmp(&right.1.to_string()))
    });
    let atoms = atoms.into_iter().map(|(_, atom)| atom).collect::<Vec<_>>();
    let mut bonds = molecules
        .iter()
        .flat_map(smiles_perceived_bonds_json)
        .collect::<Vec<_>>();
    bonds.sort_by_key(Value::to_string);
    json!({
        "status": "ok",
        "atom_count": atom_count,
        "bond_count": bond_count,
        "atoms": atoms,
        "bonds": bonds,
    })
}

pub(crate) fn hydrogen_transform_semantic_json(mut molecule: Molecule) -> Value {
    let _ = valence::perceive_valence_with_options(
        &mut molecule,
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    );
    let mol = &molecule;
    let atoms = mol
        .atoms()
        .map(|(id, atom)| {
            let mut neighbors = mol
                .neighbors(id)
                .expect("live atoms have valid adjacency")
                .map(AtomId::index)
                .collect::<Vec<_>>();
            neighbors.sort();
            json!({
                "atom_index": id.index(),
                "atomic_number": atom.element.atomic_number(),
                "symbol": atom.element.symbol(),
                "formal_charge": atom.formal_charge,
                "isotope": atom.isotope,
                "atom_map": atom.atom_map,
                "encoded_hydrogens": usize::from(atom.hydrogens.explicit_count())
                    + usize::from(mol.implicit_hydrogens(id).ok().flatten().unwrap_or(0)),
                "neighbors": neighbors,
            })
        })
        .collect::<Vec<_>>();
    json!({
        "status": "ok",
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": atoms,
    })
}

pub(crate) fn smiles_isomeric_stereo_semantic_json(mut molecule: Molecule) -> Value {
    if molecule.perceive().is_err() {
        return json!({ "status": "perception_error" });
    }
    if stereo::assign_cip_descriptors(&mut molecule).is_err() {
        return json!({ "status": "cip_error" });
    }
    let mol = &molecule;
    json!({
        "status": "ok",
        "atom_descriptors": smiles_cip_atom_descriptor_keys_json(mol),
        "bond_descriptors": smiles_cip_bond_descriptor_keys_json(mol),
    })
}

pub(crate) fn smiles_cip_atom_descriptor_keys_json(mol: &Molecule) -> Vec<Value> {
    let mut descriptors = mol
        .stereo_elements()
        .filter_map(|(id, element)| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => mol
                .cip_descriptor(id)
                .ok()
                .flatten()
                .and_then(|descriptor| {
                    let atom = mol.atom(stereo.center).ok()?;
                    Some(json!({
                        "center_atom": smiles_perceived_atom_key(mol, stereo.center, atom),
                        "descriptor": stereo_descriptor_json(descriptor),
                    }))
                }),
            StereoElementKind::Axis(_) | StereoElementKind::DoubleBond(_) => None,
        })
        .collect::<Vec<_>>();
    descriptors.sort_by_key(|value| value.to_string());
    descriptors
}

pub(crate) fn smiles_cip_bond_descriptor_keys_json(mol: &Molecule) -> Vec<Value> {
    let mut descriptors = mol
        .stereo_elements()
        .filter_map(|(id, element)| match &element.kind {
            StereoElementKind::DoubleBond(stereo) => mol
                .cip_descriptor(id)
                .ok()
                .flatten()
                .and_then(|descriptor| {
                    let left = mol.atom(stereo.left).ok()?;
                    let right = mol.atom(stereo.right).ok()?;
                    let mut endpoint_atoms = [
                        smiles_perceived_atom_key(mol, stereo.left, left),
                        smiles_perceived_atom_key(mol, stereo.right, right),
                    ];
                    endpoint_atoms.sort();
                    Some(json!({
                        "endpoint_atoms": endpoint_atoms,
                        "descriptor": stereo_descriptor_json(descriptor),
                    }))
                }),
            StereoElementKind::Axis(stereo) => {
                mol.cip_descriptor(id)
                    .ok()
                    .flatten()
                    .and_then(|descriptor| {
                        let bond = mol.bond(stereo.axis).ok()?;
                        let (begin, end) = bond.endpoints();
                        let begin_atom = mol.atom(begin).ok()?;
                        let end_atom = mol.atom(end).ok()?;
                        let mut endpoint_atoms = [
                            smiles_perceived_atom_key(mol, begin, begin_atom),
                            smiles_perceived_atom_key(mol, end, end_atom),
                        ];
                        endpoint_atoms.sort();
                        Some(json!({
                            "endpoint_atoms": endpoint_atoms,
                            "descriptor": stereo_descriptor_json(descriptor),
                        }))
                    })
            }
            StereoElementKind::Tetrahedral(_) => None,
        })
        .collect::<Vec<_>>();
    descriptors.sort_by_key(|value| value.to_string());
    descriptors
}

pub(crate) fn smiles_perceived_bonds_json(mol: &Molecule) -> Vec<Value> {
    let mut bonds = mol
        .bonds()
        .map(|(bond_id, bond)| {
            let left = mol.atom(bond.a()).expect("bond endpoint should exist");
            let right = mol.atom(bond.b()).expect("bond endpoint should exist");
            let mut endpoints = [
                smiles_perceived_atom_key(mol, bond.a(), left),
                smiles_perceived_atom_key(mol, bond.b(), right),
            ];
            endpoints.sort();
            json!({
                "endpoint_atoms": endpoints,
                "bond_type": smiles_semantic_bond_type(mol, bond_id, bond),
                "is_aromatic": mol.bond_is_aromatic(bond_id).ok().flatten().unwrap_or(false),
            })
        })
        .collect::<Vec<_>>();
    bonds.sort_by_key(|value| value.to_string());
    bonds
}

pub(crate) fn smiles_perceived_atoms_json(mol: &Molecule) -> Vec<Value> {
    let mut atoms = smiles_perceived_atom_entries_json(mol);
    atoms.sort_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then_with(|| left.1.to_string().cmp(&right.1.to_string()))
    });
    atoms.into_iter().map(|(_, value)| value).collect()
}

pub(crate) fn smiles_perceived_atom_entries_json(mol: &Molecule) -> Vec<(String, Value)> {
    mol
        .atoms()
        .map(|(id, atom)| {
            let (explicit_hydrogens, implicit_hydrogens) =
                smiles_effective_hydrogens(mol, id, atom);
            let no_implicit_hydrogens =
                smiles_effective_no_implicit_hydrogens(mol, id, atom);
            let explicit_valence = explicit_valence_json(mol, id) + explicit_hydrogens;
            let mut neighbors = mol
                .incident_bonds(id)
                .expect("atom should exist")
                .map(|(bond_id, bond)| {
                    let neighbor_id = if bond.a() == id { bond.b() } else { bond.a() };
                    let neighbor = mol.atom(neighbor_id).expect("bond endpoint should exist");
                    json!({
                        "atom": smiles_perceived_atom_key(mol, neighbor_id, neighbor),
                        "bond_type": smiles_semantic_bond_type(mol, bond_id, bond),
                        "is_aromatic": mol.bond_is_aromatic(bond_id).ok().flatten().unwrap_or(false),
                    })
                })
                .collect::<Vec<_>>();
            neighbors.sort_by_key(|value| value.to_string());
            (
                smiles_perceived_atom_key(mol, id, atom),
                json!({
                    "atomic_number": atom.element.atomic_number(),
                    "symbol": atom.element.symbol(),
                    "formal_charge": atom.formal_charge,
                    "isotope": atom.isotope,
                    "explicit_hydrogens": explicit_hydrogens,
                    "implicit_hydrogens": implicit_hydrogens,
                    "no_implicit_hydrogens": no_implicit_hydrogens,
                    "explicit_valence": explicit_valence,
                    "atom_map": atom.atom_map,
                    "aromatic": mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false),
                    "neighbors": neighbors,
                }),
            )
        })
        .collect()
}

pub(crate) fn smiles_perceived_atom_key(mol: &Molecule, id: AtomId, atom: &Atom) -> String {
    let (explicit_hydrogens, implicit_hydrogens) = smiles_effective_hydrogens(mol, id, atom);
    let no_implicit_hydrogens = smiles_effective_no_implicit_hydrogens(mol, id, atom);
    let explicit_valence = explicit_valence_json(mol, id) + explicit_hydrogens;
    format!(
        "{:03}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        atom.element.atomic_number(),
        atom.element.symbol(),
        atom.formal_charge,
        atom.isotope.unwrap_or(0),
        explicit_hydrogens,
        implicit_hydrogens,
        no_implicit_hydrogens,
        explicit_valence,
        atom.atom_map.unwrap_or(0),
        mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false)
    )
}

pub(crate) fn smiles_semantic_bond_type(mol: &Molecule, id: BondId, bond: &Bond) -> &'static str {
    if mol.bond_is_aromatic(id).ok().flatten().unwrap_or(false) {
        "AROMATIC"
    } else {
        bond_order_json(bond.order)
    }
}

pub(crate) fn smiles_effective_hydrogens(mol: &Molecule, id: AtomId, atom: &Atom) -> (u8, u8) {
    let implicit = mol.implicit_hydrogens(id).ok().flatten().unwrap_or(0);
    // Normalize only the reference-facing benchmark record. The molecule
    // retains the represented explicit/perceived implicit split.
    if atom.element.symbol() == "N"
        && mol.atom_is_aromatic(id).ok().flatten() == Some(true)
        && atom.hydrogens.explicit_count() == 0
        && implicit == 1
    {
        (1, 0)
    } else {
        (atom.hydrogens.explicit_count(), implicit)
    }
}

pub(crate) fn smiles_effective_no_implicit_hydrogens(
    mol: &Molecule,
    id: AtomId,
    atom: &Atom,
) -> bool {
    if atom.element.symbol() == "N"
        && mol.atom_is_aromatic(id).ok().flatten() == Some(true)
        && atom.formal_charge == 0
        && (atom.hydrogens.explicit_count() > 0
            || mol.implicit_hydrogens(id).ok().flatten() == Some(1))
    {
        false
    } else {
        !atom.hydrogens.allows_implicit()
    }
}
