//! Direct observations of public molecular APIs. No reference chemistry here.
use super::chemistry::{atom_json, bond_order_json};
use super::io::read_stereo_perception_records_by_suffix;
use crate::*;

pub(super) fn molecular(feature: &str, input: &Input) -> Result<Value, Box<dyn Error>> {
    let mut records = read_stereo_perception_records_by_suffix(input)?;
    let properties = if input.extension().is_some_and(|s| s == "sdf") {
        sdf::parse_str(&input.text)?
            .records()
            .iter()
            .map(|r| {
                r.data_fields()
                    .iter()
                    .map(|f| json!({"name":f.name(),"value":f.value()}))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    } else {
        vec![Vec::new(); records.len()]
    };
    let values = records.iter_mut().enumerate().map(|(index, record)| {
        if record.components.is_empty() { return Err(boxed_error("no parsed molecular components")); }
        let mut components = Vec::new();
        for (molecule, positions) in record.components.iter_mut().zip(&record.positions) {
            molecule.perceive()?;
            if feature == "stereo.perception" {
                if let Some(positions) = positions {
                    let mut editor = molecule.edit();
                    stereo::materialize_coordinate_stereo(&mut editor, positions)?;
                    *molecule = editor.finish()?;
                    molecule.perceive()?;
                }
            }
            let mut value = graph(molecule, positions.as_ref())?;
            if feature == "stereo.perception" {
                let mut candidates = stereo::detect_stereo_candidates(molecule).iter().map(|candidate| match candidate {
                    StereoCandidate::Tetrahedral { center, .. } => json!({"type":"tetrahedral","focus":[center.raw()]}),
                    StereoCandidate::DoubleBond { bond, .. } => {
                        let bond = molecule.bond(*bond).expect("live candidate bond");
                        let mut ends = [bond.a().raw(), bond.b().raw()]; ends.sort();
                        json!({"type":"double_bond","focus":ends})
                    }
                }).collect::<Vec<_>>();
                candidates.sort_by_key(Value::to_string);
                value["candidates"] = json!(candidates);
            }
            components.push(value);
        }
        Ok(json!({"record_index":index,"status":"ok","title":record.title,
            "components":components,"properties":properties.get(index).ok_or("missing record properties")?}))
    }).collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    Ok(json!({"records":values}))
}

pub(super) fn graph(
    mol: &Molecule,
    positions: Option<&kekule::structure::Positions>,
) -> Result<Value, Box<dyn Error>> {
    let index: BTreeMap<_, _> = mol
        .atom_ids()
        .enumerate()
        .map(|(index, id)| (id, index as i64))
        .collect();
    let atoms = mol
        .atoms()
        .map(|(id, atom)| {
            let mut value = atom_json(mol, id, atom);
            value["index"] = json!(index[&id]);
            value["implicit_hydrogens"] = json!(mol.implicit_hydrogens(id)?);
            value["explicit_valence"] = json!(valence::represented_valence(mol, id)?);
            value["no_implicit_hydrogens"] = json!(!atom.hydrogens.allows_implicit());
            if let Some(positions) = positions {
                let point = positions
                    .position_at(id.index())?
                    .value_in(kekule::units::ANGSTROM)?;
                value["coord"] = json!([point.x, point.y, point.z]);
            }
            Ok(value)
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let bonds = mol.bonds().map(|(id,bond)| {
        let mut ends = [index[&bond.a()],index[&bond.b()]];
        if bond.order != BondOrder::Dative { ends.sort(); }
        Ok(json!({"begin_atom_index":ends[0],"end_atom_index":ends[1],
            "bond_type":if mol.bond_is_aromatic(id)? == Some(true) { "AROMATIC" } else { bond_order_json(bond.order) },"is_aromatic":mol.bond_is_aromatic(id)?}))
    }).collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let carrier = |c: &StereoCarrier| match c {
        StereoCarrier::Atom(id) => index[id],
        StereoCarrier::ImplicitHydrogen => -1,
        StereoCarrier::ImplicitLonePair => -2,
    };
    let mut elements = Vec::new();
    let mut focuses = BTreeMap::new();
    for (id, element) in mol.stereo_elements() {
        let value = match &element.kind {
            StereoElementKind::Tetrahedral(s) => {
                let mut carriers = s.carriers.iter().map(carrier).collect::<Vec<_>>();
                let parity = s.orientation.map(|o| {
                    u8::from(o == TetrahedralOrientation::CounterClockwise)
                        ^ permutation_parity(&carriers)
                });
                carriers.sort();
                json!({"type":"tetrahedral","focus":[index[&s.center]],"carriers":carriers,"parity":parity})
            }
            StereoElementKind::DoubleBond(s) => {
                let mut ends = [index[&s.left], index[&s.right]];
                let selected = [carrier(&s.left_carrier), carrier(&s.right_carrier)];
                let canonical = [
                    canonical_carrier(mol, s.left, s.right, &index)?,
                    canonical_carrier(mol, s.right, s.left, &index)?,
                ];
                let parity = s.orientation.map(|o| {
                    u8::from(o == DoubleBondOrientation::Opposite)
                        ^ u8::from(selected[0] != canonical[0])
                        ^ u8::from(selected[1] != canonical[1])
                });
                let mut carriers = canonical;
                if ends[0] > ends[1] {
                    ends.swap(0, 1);
                    carriers.swap(0, 1);
                }
                json!({"type":"double_bond","focus":ends,"carriers":carriers,"parity":parity})
            }
            StereoElementKind::Axis(s) => {
                let bond = mol.bond(s.axis)?;
                let mut ends = [index[&bond.a()], index[&bond.b()]];
                let selected = s.carriers.iter().map(carrier).collect::<Vec<_>>();
                let mut carriers = [
                    canonical_carrier(mol, bond.a(), bond.b(), &index)?,
                    canonical_carrier(mol, bond.b(), bond.a(), &index)?,
                ];
                if selected.len() != 2 {
                    return Err(boxed_error("axis stereo must have two endpoint carriers"));
                }
                let parity = s.orientation.map(|o| {
                    u8::from(o == AxisOrientation::CounterClockwise)
                        ^ u8::from(selected[0] != carriers[0])
                        ^ u8::from(selected[1] != carriers[1])
                });
                if ends[0] > ends[1] {
                    ends.swap(0, 1);
                    carriers.swap(0, 1);
                }
                json!({"type":"axis","focus":ends,"carriers":carriers,"parity":parity})
            }
        };
        focuses.insert(id, json!({"type":value["type"],"focus":value["focus"]}));
        elements.push(value);
    }
    elements.sort_by_key(Value::to_string);
    let mut groups = mol
        .stereo_groups()
        .map(|(_, group)| {
            let mut members = group
                .members
                .iter()
                .map(|id| focuses[id].clone())
                .collect::<Vec<_>>();
            members.sort_by_key(Value::to_string);
            json!({"kind":super::stereo::stereo_group_kind_json(group.kind),"members":members})
        })
        .collect::<Vec<_>>();
    groups.sort_by_key(Value::to_string);
    Ok(
        json!({"atom_count":atoms.len(),"bond_count":bonds.len(),"atoms":atoms,"bonds":bonds,"stereo":elements,"groups":groups}),
    )
}

fn canonical_carrier(
    mol: &Molecule,
    center: AtomId,
    other: AtomId,
    index: &BTreeMap<AtomId, i64>,
) -> Result<i64, Box<dyn Error>> {
    Ok(mol
        .neighbors(center)?
        .filter(|id| *id != other)
        .map(|id| index[&id])
        .min()
        .unwrap_or(-1))
}

fn permutation_parity(values: &[i64]) -> u8 {
    let mut parity = 0;
    for i in 0..values.len() {
        for j in i + 1..values.len() {
            parity ^= u8::from(values[i] > values[j]);
        }
    }
    parity
}

pub(super) fn write(feature: &str, input: &Input) -> Result<Value, Box<dyn Error>> {
    let mut written = Vec::new();
    if feature.starts_with("io.smiles.") {
        for record in super::io::read_smiles_records(input)? {
            if record.status != "ok" || record.components.is_empty() {
                return Err(boxed_error("SMILES parse failed"));
            }
            let mode = match feature {
                "io.smiles.canonical" => smiles::SmilesWriteMode::Canonical,
                "io.smiles.isomeric" => smiles::SmilesWriteMode::Isomeric,
                _ => smiles::SmilesWriteMode::default(),
            };
            let mut components = record.components;
            let mut output = Vec::new();
            for molecule in &mut components {
                molecule.perceive()?;
                output.push(smiles::write_molecule(
                    molecule,
                    smiles::SmilesWriteOptions { mode },
                )?);
            }
            if mode == smiles::SmilesWriteMode::Canonical {
                output.sort();
            }
            let text = output.join(".");
            if feature == "io.smiles.canonical" {
                // Canonical text must be invariant to atom numbering and a fixed
                // point of reading/writing. These checks supplement RDKit identity.
                for molecule in &components {
                    let original =
                        smiles::write_molecule(molecule, smiles::SmilesWriteOptions { mode })?;
                    let permuted = reverse_numbering(molecule)?;
                    if smiles::write_molecule(&permuted, smiles::SmilesWriteOptions { mode })?
                        != original
                    {
                        return Err(boxed_error(
                            "canonical writer depends on atom or bond numbering",
                        ));
                    }
                    let reparsed = smiles::to_molecules(&original)?;
                    if reparsed.len() != 1 {
                        return Err(boxed_error("canonical output changed component count"));
                    }
                    let mut reparsed = reparsed[0].clone();
                    reparsed.perceive()?;
                    if smiles::write_molecule(&reparsed, smiles::SmilesWriteOptions { mode })?
                        != original
                    {
                        return Err(boxed_error("canonical writer has no fixed point"));
                    }
                }
            }
            written.push(json!({"path":"output.smi","text":format!("{} {}",text,record.title)}));
        }
    } else {
        let records = if input.extension().is_some_and(|s| s == "sdf") {
            sdf::interpret(&sdf::parse_str(&input.text)?)?.into_records()
        } else {
            let doc = molfile::parse_str(&input.text)?;
            vec![SdfRecordInterpretation::new(
                doc.header().title(),
                molfile::interpret(&doc)?.into_model(),
                vec![],
            )]
        };
        if feature.starts_with("io.sdf.") {
            written.push(json!({"path":"output.sdf","text":sdf::write_v2000(&records)?}));
        } else {
            let version = if feature.contains("v3000") {
                molfile::MolfileWriteVersion::V3000
            } else {
                molfile::MolfileWriteVersion::V2000
            };
            for record in records {
                written.push(json!({"path":"output.mol","text":molfile::write_model(record.model(), molfile::MolfileWriteOptions { version })?}));
            }
        }
    }
    Ok(json!({"written":written}))
}

// Change only numbering through the public construction API. Retain every atom,
// directed bond, stereo carrier and relation group. This is a writer contract
// check, not an alternative chemistry algorithm or a reference answer.
fn reverse_numbering(source: &Molecule) -> Result<Molecule, Box<dyn Error>> {
    let mut editor = kekule::core::MoleculeEditor::new();
    let mut atoms = BTreeMap::new();
    for (id, atom) in source.atoms().collect::<Vec<_>>().into_iter().rev() {
        atoms.insert(id, editor.add_atom(atom.clone())?);
    }
    let mut bonds = BTreeMap::new();
    for (id, bond) in source.bonds().collect::<Vec<_>>().into_iter().rev() {
        bonds.insert(
            id,
            editor.add_bond(atoms[&bond.a()], atoms[&bond.b()], bond.order)?,
        );
    }
    let remap = |carrier: &mut StereoCarrier| {
        if let StereoCarrier::Atom(atom) = carrier {
            *atom = atoms[atom];
        }
    };
    let mut elements = BTreeMap::new();
    for (id, element) in source.stereo_elements() {
        let mut kind = element.kind.clone();
        match &mut kind {
            StereoElementKind::Tetrahedral(s) => {
                s.center = atoms[&s.center];
                for carrier in &mut s.carriers {
                    remap(carrier);
                }
            }
            StereoElementKind::DoubleBond(s) => {
                s.bond = bonds[&s.bond];
                s.left = atoms[&s.left];
                s.right = atoms[&s.right];
                remap(&mut s.left_carrier);
                remap(&mut s.right_carrier);
            }
            StereoElementKind::Axis(s) => {
                s.axis = bonds[&s.axis];
                for carrier in &mut s.carriers {
                    remap(carrier);
                }
            }
        }
        elements.insert(
            id,
            editor.add_stereo_element(kekule::core::StereoElement::new(kind))?,
        );
    }
    for (_, group) in source.stereo_groups() {
        editor.add_stereo_group(kekule::core::StereoGroup {
            kind: group.kind,
            members: group.members.iter().map(|id| elements[id]).collect(),
        })?;
    }
    let mut molecule = editor.finish()?;
    molecule.perceive()?;
    Ok(molecule)
}
