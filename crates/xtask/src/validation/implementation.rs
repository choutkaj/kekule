use crate::*;

pub(crate) fn implementation_expected(
    feature: &str,
    corpus: &str,
    fixture_path: &Path,
) -> Result<Value, Box<dyn Error>> {
    match feature {
        "io.sdf.v2000.parse" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(sdf_record_json).collect::<Vec<_>>() }))
        }
        "io.sdf.v2000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let molecules = records
                .into_iter()
                .map(|record| record.molecule)
                .collect::<Vec<_>>();
            let written = write_sdf_v2000(&molecules)?;
            let records = read_sdf_v2000_records(&written, SdfParseOptions::default())?
                .into_iter()
                .enumerate()
                .map(|(index, record)| IndexedSmallRecord {
                    record_index: index,
                    title: record.title,
                    molecule: record.molecule,
                })
                .collect::<Vec<_>>();
            Ok(json!({ "records": records.iter().map(sdf_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v2000.parse" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(mol_parse_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v3000.parse" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let records = records
                .into_iter()
                .enumerate()
                .map(|(index, record)| {
                    let written = write_mol_v3000(&record.molecule)?;
                    let molecule = read_mol_v3000_str(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title: molecule_title(&molecule.mol),
                        molecule,
                    })
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
            Ok(json!({ "records": records.iter().map(mol_parse_record_json).collect::<Vec<_>>() }))
        }
        "core.conformers" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(conformer_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v2000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let records = records
                .into_iter()
                .enumerate()
                .map(|(index, record)| {
                    let written = write_mol_v2000(&record.molecule)?;
                    let molecule = read_mol_v2000_str(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title: molecule_title(&molecule.mol),
                        molecule,
                    })
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
            Ok(json!({ "records": records.iter().map(mol_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v3000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let records = records
                .into_iter()
                .enumerate()
                .map(|(index, record)| {
                    let written = write_mol_v3000(&record.molecule)?;
                    let molecule = read_mol_v3000_str(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title: molecule_title(&molecule.mol),
                        molecule,
                    })
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
            Ok(json!({ "records": records.iter().map(mol_record_json).collect::<Vec<_>>() }))
        }
        "io.smiles.parse" => {
            let records = read_smiles_records(fixture_path)?;
            Ok(
                json!({ "records": records.iter().map(smiles_parse_record_json).collect::<Vec<_>>() }),
            )
        }
        "io.smiles.write" => {
            let records = read_smiles_records(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter()
                    .map(smiles_write_record_json)
                    .collect::<Result<Vec<_>, Box<dyn Error>>>()?
            }))
        }
        "io.smiles.canonical" => {
            let records = read_canonical_smiles_records(fixture_path)?;
            let exact_smiles = corpus == "tiny";
            Ok(json!({
                "records": records
                    .iter()
                    .map(|record| canonical_smiles_record_json(record, exact_smiles))
                    .collect::<Result<Vec<_>, Box<dyn Error>>>()?
            }))
        }
        "algo.rings.fast" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(ring_membership_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.rings.sssr" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(ring_set_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.valence.rdkit-like" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(valence_record_json).collect::<Vec<_>>() }),
            )
        }
        "chem.sanitize.rdkit-like" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(sanitized_atom_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.aromaticity.rdkit-like" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(aromaticity_record_json).collect::<Vec<_>>() }),
            )
        }
        "io.mmcif.parse" | "bio.hierarchy.smcra" => {
            let input = fs::read_to_string(fixture_path)?;
            let molecule = read_mmcif_str(&input, MmcifParseOptions::default())?;
            Ok(mmcif_expected_json(&molecule))
        }
        _ => Err(boxed_error(format!(
            "no implementation comparison configured for feature `{feature}`"
        ))),
    }
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedSmallRecord {
    pub(crate) record_index: usize,
    pub(crate) title: String,
    pub(crate) molecule: SmallMolecule,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedSmilesRecord {
    pub(crate) record_index: usize,
    pub(crate) status: String,
    pub(crate) title: String,
    pub(crate) input_smiles: String,
    pub(crate) molecule: Option<SmallMolecule>,
}

pub(crate) fn read_small_records_by_suffix(
    path: &Path,
) -> Result<Vec<IndexedSmallRecord>, Box<dyn Error>> {
    let input = fs::read_to_string(path)?;
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("mol" | "mdl")
    ) {
        let molecule = read_mol_v2000_str(&input)?;
        return Ok(vec![IndexedSmallRecord {
            record_index: 0,
            title: molecule_title(&molecule.mol),
            molecule,
        }]);
    }
    Ok(read_sdf_v2000_records(&input, SdfParseOptions::default())?
        .into_iter()
        .enumerate()
        .map(|(index, record)| small_record(index, record))
        .collect())
}

pub(crate) fn small_record(index: usize, record: SdfRecord) -> IndexedSmallRecord {
    IndexedSmallRecord {
        record_index: index,
        title: record.title,
        molecule: record.molecule,
    }
}

pub(crate) fn read_smiles_records(path: &Path) -> Result<Vec<IndexedSmilesRecord>, Box<dyn Error>> {
    let mut records = Vec::new();
    for (index, raw_line) in fs::read_to_string(path)?.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(2, char::is_whitespace);
        let smiles = parts.next().unwrap_or_default().to_owned();
        let title = parts.next().unwrap_or_default().trim().to_owned();
        if smiles_unsupported_subset_reason(&smiles).is_some() {
            records.push(IndexedSmilesRecord {
                record_index: index,
                status: "unsupported".to_owned(),
                title,
                input_smiles: smiles,
                molecule: None,
            });
            continue;
        }
        let (status, molecule) = match read_smiles_str(&smiles, SmilesParseOptions) {
            Ok(molecule) => ("ok".to_owned(), Some(molecule)),
            Err(_) => ("parse_error".to_owned(), None),
        };
        records.push(IndexedSmilesRecord {
            record_index: index,
            status,
            title,
            input_smiles: smiles,
            molecule,
        });
    }
    Ok(records)
}

pub(crate) fn read_canonical_smiles_records(
    path: &Path,
) -> Result<Vec<IndexedSmilesRecord>, Box<dyn Error>> {
    let mut records = Vec::new();
    for (index, raw_line) in fs::read_to_string(path)?.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(2, char::is_whitespace);
        let smiles = parts.next().unwrap_or_default().to_owned();
        let title = parts.next().unwrap_or_default().trim().to_owned();
        let (status, molecule) = match read_smiles_str(&smiles, SmilesParseOptions) {
            Ok(molecule) => ("ok".to_owned(), Some(molecule)),
            Err(_) => ("parse_error".to_owned(), None),
        };
        records.push(IndexedSmilesRecord {
            record_index: index,
            status,
            title,
            input_smiles: smiles,
            molecule,
        });
    }
    Ok(records)
}

pub(crate) fn smiles_unsupported_subset_reason(smiles: &str) -> Option<&'static str> {
    smiles
        .chars()
        .any(|ch| matches!(ch, '@' | '/' | '\\' | '*'))
        .then_some("unsupported")
}

pub(crate) fn sdf_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": atoms_json(mol),
        "bonds": bonds_json(mol),
        "properties": sdf_properties_json(mol),
    })
}

pub(crate) fn mol_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": atoms_json(mol),
        "bonds": bonds_json(mol),
    })
}

pub(crate) fn conformer_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "conformers": mol.conformers().map(|(_, conformer)| {
            mol.atom_ids()
                .filter_map(|atom_id| {
                    conformer.position(atom_id).map(|point| json!({
                        "atom_index": atom_id.raw(),
                        "x": point.x,
                        "y": point.y,
                        "z": point.z,
                    }))
                })
                .collect::<Vec<_>>()
        }).collect::<Vec<_>>(),
        "atoms": mol.atoms().map(|(id, atom)| conformer_atom_json(id, atom)).collect::<Vec<_>>(),
    })
}

pub(crate) fn mol_parse_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "conformers": conformers_json(mol),
        "atoms": atoms_json(mol),
    })
}

pub(crate) fn conformers_json(mol: &Molecule) -> Vec<Vec<Value>> {
    mol.conformers()
        .map(|(_, conformer)| {
            mol.atom_ids()
                .filter_map(|atom_id| {
                    conformer.position(atom_id).map(|point| {
                        json!({
                            "atom_index": atom_id.raw(),
                            "x": point.x,
                            "y": point.y,
                            "z": point.z,
                        })
                    })
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
}

pub(crate) fn conformer_atom_json(id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "atom_map": atom.atom_map,
        "aromatic": atom.aromatic,
    })
}

pub(crate) fn ring_membership_record_json(record: &mut IndexedSmallRecord) -> Value {
    let membership = perceive_ring_membership(&mut record.molecule.mol);
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_in_ring": mol.atom_ids().map(|id| membership.atom_in_ring(id)).collect::<Vec<_>>(),
        "bond_in_ring": mol.bond_ids().map(|id| membership.bond_in_ring(id)).collect::<Vec<_>>(),
    })
}

pub(crate) fn ring_set_record_json(record: &mut IndexedSmallRecord) -> Value {
    match perceive_ring_set(&mut record.molecule.mol) {
        Ok(ring_set) => json!({
            "record_index": record.record_index,
            "status": "ok",
            "title": record.title,
            "rings": ring_set
                .rings()
                .iter()
                .map(|ring| ring.atoms.iter().map(|atom| atom.raw()).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
        }),
        Err(_) => json!({
            "record_index": record.record_index,
            "status": "resource_error",
            "title": record.title,
        }),
    }
}

pub(crate) fn sanitized_atom_record_json(record: &mut IndexedSmallRecord) -> Value {
    match sanitize_small_molecule(&mut record.molecule, SanitizeOptions::default()) {
        Ok(_) => json!({
            "record_index": record.record_index,
            "status": "ok",
            "title": record.title,
            "atoms": basic_atoms_json(&record.molecule.mol),
        }),
        Err(_) => json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
        }),
    }
}

pub(crate) fn valence_record_json(record: &mut IndexedSmallRecord) -> Value {
    let report = perceive_valence(&mut record.molecule.mol, ValenceModel::RdkitLike);
    if !report.is_ok() {
        return json!({
            "record_index": record.record_index,
            "status": "valence_error",
            "title": record.title,
        });
    }
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atoms": record
            .molecule
            .mol
            .atoms()
            .map(|(id, atom)| valence_atom_json(&record.molecule.mol, id, atom))
            .collect::<Vec<_>>(),
    })
}

pub(crate) fn aromaticity_record_json(record: &mut IndexedSmallRecord) -> Value {
    let status = sanitize_small_molecule(
        &mut record.molecule,
        SanitizeOptions {
            perceive_valence: true,
            perceive_rings: true,
            perceive_aromaticity: false,
        },
    )
    .and_then(|_| {
        perceive_aromaticity(&mut record.molecule.mol, AromaticityModel::RdkitLike)
            .map_err(molecules::prelude::SanitizeError::Aromaticity)
    });
    if status.is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
        });
    }
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_aromatic": mol.atoms().map(|(_, atom)| atom.aromatic).collect::<Vec<_>>(),
        "bond_aromatic": mol.bonds().map(|(_, bond)| bond.aromatic).collect::<Vec<_>>(),
    })
}

pub(crate) fn smiles_write_record_json(
    record: &IndexedSmilesRecord,
) -> Result<Value, Box<dyn Error>> {
    let Some(molecule) = &record.molecule else {
        return Ok(smiles_error_record_json(record));
    };
    let written = write_smiles(molecule, SmilesWriteOptions)?;
    let reparsed = match read_smiles_str(&written, SmilesParseOptions) {
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
        "sanitized": smiles_sanitized_semantic_json(reparsed),
    }))
}

pub(crate) fn canonical_smiles_record_json(
    record: &IndexedSmilesRecord,
    exact_smiles: bool,
) -> Result<Value, Box<dyn Error>> {
    let Some(molecule) = &record.molecule else {
        return Ok(smiles_error_record_json(record));
    };
    let written = write_canonical_smiles(molecule, CanonicalSmilesWriteOptions)?;
    let reparsed = match read_smiles_str(&written, SmilesParseOptions) {
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
        "sanitized": smiles_sanitized_semantic_json(reparsed),
    });
    if exact_smiles {
        item["canonical_smiles"] = json!(written);
    }
    Ok(item)
}

pub(crate) fn smiles_parse_record_json(record: &IndexedSmilesRecord) -> Value {
    let Some(molecule) = &record.molecule else {
        return smiles_error_record_json(record);
    };
    let written = write_smiles(molecule, SmilesWriteOptions);
    let round_trip = match written
        .as_ref()
        .map_err(|_| ())
        .and_then(|text| read_smiles_str(text, SmilesParseOptions).map_err(|_| ()))
    {
        Ok(reparsed) => smiles_sanitized_semantic_json(reparsed),
        Err(_) => json!({ "status": "write_reparse_error" }),
    };
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "input_smiles": record.input_smiles,
        "raw": smiles_raw_semantic_json(molecule),
        "sanitized": smiles_sanitized_semantic_json(molecule.clone()),
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

pub(crate) fn smiles_raw_semantic_json(molecule: &SmallMolecule) -> Value {
    let mol = &molecule.mol;
    json!({
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": basic_atoms_json(mol),
        "bonds": basic_bonds_json(mol),
    })
}

pub(crate) fn smiles_sanitized_semantic_json(mut molecule: SmallMolecule) -> Value {
    match sanitize_small_molecule(&mut molecule, SanitizeOptions::default()) {
        Ok(_) => {
            let mol = &molecule.mol;
            json!({
                "status": "ok",
                "atom_count": mol.atom_count(),
                "bond_count": mol.bond_count(),
                "atoms": smiles_sanitized_atoms_json(mol),
                "bonds": smiles_sanitized_bonds_json(mol),
            })
        }
        Err(_) => json!({ "status": "sanitize_error" }),
    }
}

pub(crate) fn smiles_sanitized_bonds_json(mol: &Molecule) -> Vec<Value> {
    let mut bonds = mol
        .bonds()
        .map(|(_, bond)| {
            let left = mol.atom(bond.a()).expect("bond endpoint should exist");
            let right = mol.atom(bond.b()).expect("bond endpoint should exist");
            let mut endpoints = [
                smiles_sanitized_atom_key(mol, bond.a(), left),
                smiles_sanitized_atom_key(mol, bond.b(), right),
            ];
            endpoints.sort();
            json!({
                "endpoint_atoms": endpoints,
                "bond_type": smiles_semantic_bond_type(bond),
                "is_aromatic": bond.aromatic,
            })
        })
        .collect::<Vec<_>>();
    bonds.sort_by_key(|value| value.to_string());
    bonds
}

pub(crate) fn smiles_sanitized_atoms_json(mol: &Molecule) -> Vec<Value> {
    let mut atoms = mol
        .atoms()
        .map(|(id, atom)| {
            let (explicit_hydrogens, implicit_hydrogens) = smiles_effective_hydrogens(atom);
            let explicit_valence = explicit_valence_json(mol, id) + explicit_hydrogens;
            let mut neighbors = mol
                .incident_bonds(id)
                .expect("atom should exist")
                .map(|(_, bond)| {
                    let neighbor_id = if bond.a() == id { bond.b() } else { bond.a() };
                    let neighbor = mol.atom(neighbor_id).expect("bond endpoint should exist");
                    json!({
                        "atom": smiles_sanitized_atom_key(mol, neighbor_id, neighbor),
                        "bond_type": smiles_semantic_bond_type(bond),
                        "is_aromatic": bond.aromatic,
                    })
                })
                .collect::<Vec<_>>();
            neighbors.sort_by_key(|value| value.to_string());
            (
                smiles_sanitized_atom_key(mol, id, atom),
                json!({
                    "atomic_number": atom.element.atomic_number(),
                    "symbol": atom.element.symbol(),
                    "formal_charge": atom.formal_charge,
                    "isotope": atom.isotope,
                    "explicit_hydrogens": explicit_hydrogens,
                    "implicit_hydrogens": implicit_hydrogens,
                    "no_implicit_hydrogens": atom.no_implicit_hydrogens,
                    "explicit_valence": explicit_valence,
                    "atom_map": atom.atom_map,
                    "aromatic": atom.aromatic,
                    "neighbors": neighbors,
                }),
            )
        })
        .collect::<Vec<_>>();
    atoms.sort_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then_with(|| left.1.to_string().cmp(&right.1.to_string()))
    });
    atoms.into_iter().map(|(_, value)| value).collect()
}

pub(crate) fn smiles_sanitized_atom_key(mol: &Molecule, id: AtomId, atom: &Atom) -> String {
    let (explicit_hydrogens, implicit_hydrogens) = smiles_effective_hydrogens(atom);
    let explicit_valence = explicit_valence_json(mol, id) + explicit_hydrogens;
    format!(
        "{:03}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        atom.element.atomic_number(),
        atom.element.symbol(),
        atom.formal_charge,
        atom.isotope.unwrap_or(0),
        explicit_hydrogens,
        implicit_hydrogens,
        atom.no_implicit_hydrogens,
        explicit_valence,
        atom.atom_map.unwrap_or(0),
        atom.aromatic
    )
}

pub(crate) fn smiles_semantic_bond_type(bond: &Bond) -> &'static str {
    if bond.aromatic {
        "AROMATIC"
    } else {
        bond_order_json(bond.order)
    }
}

pub(crate) fn smiles_effective_hydrogens(atom: &Atom) -> (u8, u8) {
    if atom.element.symbol() == "N"
        && atom.aromatic
        && atom.explicit_hydrogens == 0
        && atom.implicit_hydrogens == Some(1)
    {
        (1, 0)
    } else {
        (
            atom.explicit_hydrogens,
            atom.implicit_hydrogens.unwrap_or(0),
        )
    }
}

pub(crate) fn atoms_json(mol: &Molecule) -> Vec<Value> {
    mol.atoms()
        .map(|(id, atom)| atom_json(id, atom))
        .collect::<Vec<_>>()
}

pub(crate) fn atom_json(id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "atom_map": atom.atom_map,
        "radical": atom.radical.map(radical_json),
        "unpaired_electrons": atom.radical.map(AtomRadical::unpaired_electron_count).unwrap_or(0),
        "aromatic": atom.aromatic,
    })
}

pub(crate) fn basic_atoms_json(mol: &Molecule) -> Vec<Value> {
    mol.atoms()
        .map(|(id, atom)| basic_atom_json(id, atom))
        .collect::<Vec<_>>()
}

pub(crate) fn basic_atom_json(id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "atom_map": atom.atom_map,
        "aromatic": atom.aromatic,
    })
}

pub(crate) fn valence_atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "implicit_hydrogens": atom.implicit_hydrogens.unwrap_or(0),
        "explicit_valence": explicit_valence_json(mol, id) + atom.explicit_hydrogens,
    })
}

pub(crate) fn explicit_valence_json(mol: &Molecule, atom: AtomId) -> u8 {
    mol.incident_bonds(atom)
        .ok()
        .into_iter()
        .flatten()
        .map(|(_, bond)| match bond.order {
            BondOrder::Zero | BondOrder::Dative => 0,
            BondOrder::Single | BondOrder::Aromatic => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
        })
        .sum()
}

pub(crate) fn bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| bond_json(id.raw(), bond))
        .collect::<Vec<_>>()
}

pub(crate) fn bond_json(index: u32, bond: &Bond) -> Value {
    json!({
        "index": index,
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": bond.aromatic,
        "stereo": bond_stereo_json(bond.order, bond.stereo),
        "bond_direction": bond_direction_json(bond.order, bond.stereo),
    })
}

pub(crate) fn basic_bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| basic_bond_json(id.raw(), bond))
        .collect::<Vec<_>>()
}

pub(crate) fn basic_bond_json(index: u32, bond: &Bond) -> Value {
    json!({
        "index": index,
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": bond.aromatic,
        "stereo": legacy_bond_stereo_json(bond.stereo),
    })
}

pub(crate) fn legacy_bond_stereo_json(stereo: Option<BondStereo>) -> &'static str {
    match stereo {
        None | Some(BondStereo::Unspecified) => "STEREONONE",
        Some(BondStereo::E) => "STEREOE",
        Some(BondStereo::Z) => "STEREOZ",
        Some(BondStereo::Up) | Some(BondStereo::Down) | Some(BondStereo::Any) => "STEREOANY",
    }
}

pub(crate) fn radical_json(radical: AtomRadical) -> &'static str {
    match radical {
        AtomRadical::Singlet => "SINGLET",
        AtomRadical::Doublet => "DOUBLET",
        AtomRadical::Triplet => "TRIPLET",
    }
}

pub(crate) fn bond_order_json(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Zero => "ZERO",
        BondOrder::Single => "SINGLE",
        BondOrder::Double => "DOUBLE",
        BondOrder::Triple => "TRIPLE",
        BondOrder::Quadruple => "QUADRUPLE",
        BondOrder::Aromatic => "AROMATIC",
        BondOrder::Dative => "DATIVE",
    }
}

pub(crate) fn bond_stereo_json(order: BondOrder, stereo: Option<BondStereo>) -> &'static str {
    match (order, stereo) {
        (_, None | Some(BondStereo::Unspecified)) => "STEREONONE",
        (_, Some(BondStereo::E)) => "STEREOE",
        (_, Some(BondStereo::Z)) => "STEREOZ",
        (BondOrder::Double, Some(BondStereo::Any)) => "STEREOANY",
        _ => "STEREONONE",
    }
}

pub(crate) fn bond_direction_json(order: BondOrder, stereo: Option<BondStereo>) -> &'static str {
    match (order, stereo) {
        (BondOrder::Single, Some(BondStereo::Up)) => "BEGINWEDGE",
        (BondOrder::Single, Some(BondStereo::Down)) => "BEGINDASH",
        (BondOrder::Single, Some(BondStereo::Any)) => "UNKNOWN",
        _ => "NONE",
    }
}

pub(crate) fn sdf_properties_json(mol: &Molecule) -> Value {
    let mut props = serde_json::Map::new();
    for (key, value) in mol.props() {
        let Some(name) = key.strip_prefix("sdf.field.") else {
            continue;
        };
        if let PropValue::String(text) = value {
            props.insert(name.to_owned(), json!(text));
        }
    }
    Value::Object(props)
}

pub(crate) fn molecule_title(mol: &Molecule) -> String {
    match mol.props().get("sdf.title") {
        Some(PropValue::String(title)) => title.clone(),
        _ => String::new(),
    }
}

pub(crate) fn mmcif_expected_json(molecule: &MacroMolecule) -> Value {
    json!({
        "atom_site_rows": atom_site_rows_json(molecule),
        "structure": structure_json(molecule),
    })
}

pub(crate) fn atom_site_rows_json(molecule: &MacroMolecule) -> Value {
    let rows = molecule
        .hierarchy
        .atom_sites()
        .map(|(site_id, site)| {
            let residue = molecule
                .hierarchy
                .residue(site.residue)
                .expect("residue exists");
            let chain = molecule
                .hierarchy
                .chain(residue.chain)
                .expect("chain exists");
            let model = molecule.hierarchy.model(chain.model).expect("model exists");
            let atom = molecule.mol.atom(site.atom).expect("atom exists");
            let point = first_conformer_point(molecule, site.atom);
            json!({
                "group_PDB": site.metadata.group_pdb.clone(),
                "id": site.metadata.atom_site_id.clone().unwrap_or_else(|| (site_id.raw() + 1).to_string()),
                "type_symbol": site.metadata.type_symbol.clone().unwrap_or_else(|| atom.element.symbol().to_owned()),
                "label_atom_id": site.metadata.label_atom_id,
                "auth_atom_id": site.metadata.auth_atom_id,
                "label_alt_id": site.metadata.label_alt_id,
                "label_comp_id": residue.label_comp_id,
                "auth_comp_id": residue.author_comp_id,
                "label_asym_id": site.metadata.label_asym_id.clone().unwrap_or_else(|| chain.label_id.clone()),
                "auth_asym_id": site.metadata.auth_asym_id.clone().or_else(|| chain.author_id.clone()),
                "label_seq_id": residue.label_seq_id.map(|value| value.to_string()),
                "auth_seq_id": residue.author_seq_id,
                "pdbx_PDB_ins_code": residue.insertion_code,
                "occupancy": site.metadata.occupancy_raw.clone().or_else(|| site.metadata.occupancy.map(|value| format!("{value:.2}"))),
                "B_iso_or_equiv": site.metadata.b_factor_raw.clone().or_else(|| site.metadata.b_factor.map(|value| format!("{value:.2}"))),
                "Cartn_x": site.metadata.cartn_x_raw.clone().or_else(|| point.map(|point| format!("{:.3}", point.x))),
                "Cartn_y": site.metadata.cartn_y_raw.clone().or_else(|| point.map(|point| format!("{:.3}", point.y))),
                "Cartn_z": site.metadata.cartn_z_raw.clone().or_else(|| point.map(|point| format!("{:.3}", point.z))),
                "pdbx_PDB_model_num": model.model_id,
            })
        })
        .collect::<Vec<_>>();
    json!({
        "status": "ok",
        "row_count": rows.len(),
        "rows": rows,
    })
}

pub(crate) fn structure_json(molecule: &MacroMolecule) -> Value {
    json!({
        "status": "ok",
        "models": molecule.hierarchy.models().map(|(model_id, model)| {
            json!({
                "id": model_id.raw(),
                "chains": model.chains.iter().map(|chain_id| {
                    let chain = molecule.hierarchy.chain(*chain_id).expect("chain exists");
                    json!({
                        "id": chain.author_id.clone().unwrap_or_else(|| chain.label_id.clone()),
                        "residues": chain.residues.iter().map(|residue_id| {
                            let residue = molecule.hierarchy.residue(*residue_id).expect("residue exists");
                            json!({
                                "name": residue.name,
                                "hetflag": residue_hetflag_json(molecule, residue),
                                "sequence_id": residue_sequence_json(residue),
                                "insertion_code": residue.insertion_code,
                                "atoms": residue.atom_sites.iter().map(|site_id| {
                                    let site = molecule.hierarchy.atom_site(*site_id).expect("site exists");
                                    let atom = molecule.mol.atom(site.atom).expect("atom exists");
                                    let name = site
                                        .metadata
                                        .label_atom_id
                                        .clone()
                                        .unwrap_or_else(|| atom.element.symbol().to_owned());
                                    let coord = first_conformer_point(molecule, site.atom)
                                        .map(|point| json!([point.x, point.y, point.z]))
                                        .unwrap_or(Value::Null);
                                    json!({
                                        "name": name,
                                        "full_name": name,
                                        "altloc": site.metadata.label_alt_id,
                                        "element": site.metadata.type_symbol.clone().unwrap_or_else(|| atom.element.symbol().to_owned()),
                                        "occupancy": site.metadata.occupancy,
                                        "bfactor": site.metadata.b_factor,
                                        "coord": coord,
                                    })
                                }).collect::<Vec<_>>(),
                            })
                        }).collect::<Vec<_>>(),
                    })
                }).collect::<Vec<_>>(),
            })
        }).collect::<Vec<_>>(),
    })
}

pub(crate) fn first_conformer_point(
    molecule: &MacroMolecule,
    atom: AtomId,
) -> Option<molecules::prelude::Point3> {
    molecule
        .mol
        .first_conformer()
        .and_then(|(_, conformer)| conformer.position(atom))
}

pub(crate) fn residue_sequence_json(residue: &molecules::prelude::Residue) -> Value {
    if let Some(author_seq_id) = &residue.author_seq_id {
        return author_seq_id
            .parse::<i32>()
            .map(Value::from)
            .unwrap_or_else(|_| json!(author_seq_id));
    }
    residue.label_seq_id.map(Value::from).unwrap_or(Value::Null)
}

pub(crate) fn residue_hetflag_json(
    molecule: &MacroMolecule,
    residue: &molecules::prelude::Residue,
) -> Value {
    let is_hetatm = residue.atom_sites.iter().any(|site_id| {
        molecule
            .hierarchy
            .atom_site(*site_id)
            .ok()
            .and_then(|site| site.metadata.group_pdb.as_deref())
            == Some("HETATM")
    });
    if is_hetatm {
        if residue.name == "HOH" {
            json!("W")
        } else {
            json!(format!("H_{}", residue.name))
        }
    } else {
        Value::Null
    }
}
