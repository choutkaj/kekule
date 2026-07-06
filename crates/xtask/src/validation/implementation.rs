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
            let written = sdf::write_v2000(&molecules)?;
            let records = sdf::read_v2000_records(&written, SdfParseOptions::default())?
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
                    let written = molfile::write_v3000(&record.molecule)?;
                    let molecule = molfile::read_v3000_str(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title: molecule_title(molecule.graph()),
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
                    let written = molfile::write_v2000(&record.molecule)?;
                    let molecule = molfile::read_v2000_str(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title: molecule_title(molecule.graph()),
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
                    let written = molfile::write_v3000(&record.molecule)?;
                    let molecule = molfile::read_v3000_str(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title: molecule_title(molecule.graph()),
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
            let exact_smiles = corpus == "smoke";
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
            let molecule = bio::read_mmcif_str(&input, MmcifParseOptions::default())?;
            Ok(mmcif_expected_json(&molecule))
        }
        "stereo.representation" => {
            let records = read_stereo_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(stereo_record_json).collect::<Vec<_>>() }))
        }
        "stereo.perception" => {
            let mut records = read_stereo_records_by_suffix(fixture_path)?;
            for record in &mut records {
                perception::sanitize_with_options(
                    &mut record.molecule,
                    SanitizeOptions {
                        perceive_stereo: false,
                        ..SanitizeOptions::default()
                    },
                )?;
            }
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(stereo_perception_record_json)
                    .collect::<Vec<_>>()
            }))
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
        let molecule = molfile::read_v2000_str(&input)?;
        return Ok(vec![IndexedSmallRecord {
            record_index: 0,
            title: molecule_title(molecule.graph()),
            molecule,
        }]);
    }
    Ok(sdf::read_v2000_records(&input, SdfParseOptions::default())?
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
        let (status, molecule) =
            match smiles::read_str_with_options(&smiles, SmilesParseOptions::default()) {
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

pub(crate) fn read_stereo_records_by_suffix(
    path: &Path,
) -> Result<Vec<IndexedSmallRecord>, Box<dyn Error>> {
    let input = fs::read_to_string(path)?;
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("txt" | "smi" | "smiles")
    ) {
        return Ok(read_smiles_records(path)?
            .into_iter()
            .filter_map(|record| {
                record.molecule.map(|molecule| IndexedSmallRecord {
                    record_index: record.record_index,
                    title: record.title,
                    molecule,
                })
            })
            .collect());
    }
    let molecule = if input.contains("V3000") {
        molfile::read_v3000_str(&input)?
    } else {
        molfile::read_v2000_str(&input)?
    };
    Ok(vec![IndexedSmallRecord {
        record_index: 0,
        title: molecule_title(molecule.graph()),
        molecule,
    }])
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
        let (status, molecule) =
            match smiles::read_str_with_options(&smiles, SmilesParseOptions::default()) {
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
        .any(|ch| matches!(ch, '*'))
        .then_some("unsupported")
}

pub(crate) fn sdf_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = record.molecule.graph();
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
    let mol = record.molecule.graph();
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
    let mol = record.molecule.graph();
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
    let mol = record.molecule.graph();
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "conformers": conformers_json(mol),
        "atoms": atoms_json(mol),
    })
}

pub(crate) fn stereo_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = record.molecule.graph();
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "stereo_elements": stereo_elements_json(mol),
        "stereo_groups": stereo_groups_json(mol),
        "stereo_bond_marks": stereo_bond_marks_json(mol),
    })
}

pub(crate) fn stereo_perception_record_json(record: &mut IndexedSmallRecord) -> Value {
    let report = stereo::perceive_stereo(record.molecule.graph_mut());
    let mol = record.molecule.graph();
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "report": stereo_perception_report_json(&report),
        "stereo_elements": stereo_elements_json(mol),
        "stereo_groups": stereo_groups_json(mol),
        "stereo_bond_marks": stereo_bond_marks_json(mol),
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
    let membership = rings::perceive_ring_membership(record.molecule.graph_mut());
    let mol = record.molecule.graph();
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_in_ring": mol.atom_ids().map(|id| membership.atom_in_ring(id)).collect::<Vec<_>>(),
        "bond_in_ring": mol.bond_ids().map(|id| membership.bond_in_ring(id)).collect::<Vec<_>>(),
    })
}

pub(crate) fn ring_set_record_json(record: &mut IndexedSmallRecord) -> Value {
    match rings::perceive_ring_set(record.molecule.graph_mut()) {
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
    match perception::sanitize_with_options(&mut record.molecule, SanitizeOptions::default()) {
        Ok(_) => json!({
            "record_index": record.record_index,
            "status": "ok",
            "title": record.title,
            "atoms": basic_atoms_json(record.molecule.graph()),
        }),
        Err(_) => json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
        }),
    }
}

pub(crate) fn valence_record_json(record: &mut IndexedSmallRecord) -> Value {
    let report = valence::perceive_valence(record.molecule.graph_mut(), ValenceModel::RdkitLike);
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
            .graph()
            .atoms()
            .map(|(id, atom)| valence_atom_json(record.molecule.graph(), id, atom))
            .collect::<Vec<_>>(),
    })
}

pub(crate) fn aromaticity_record_json(record: &mut IndexedSmallRecord) -> Value {
    let status = perception::sanitize_with_options(
        &mut record.molecule,
        SanitizeOptions {
            perceive_valence: true,
            perceive_rings: true,
            perceive_aromaticity: false,
            perceive_stereo: false,
        },
    )
    .and_then(|_| {
        aromaticity::perceive_aromaticity(record.molecule.graph_mut(), AromaticityModel::RdkitLike)
            .map_err(SanitizeError::Aromaticity)
    });
    if status.is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
        });
    }
    let mol = record.molecule.graph();
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
    let written = smiles::write_with_options(molecule, SmilesWriteOptions::default())?;
    let reparsed = match smiles::read_str_with_options(&written, SmilesParseOptions::default()) {
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
    let mut molecule = molecule.clone();
    if perception::sanitize_with_options(&mut molecule, SanitizeOptions::default()).is_err() {
        return Ok(json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
            "input_smiles": record.input_smiles,
        }));
    }
    let written =
        smiles::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions::default())?;
    let reparsed = match smiles::read_str_with_options(&written, SmilesParseOptions::default()) {
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
    let written = smiles::write_with_options(molecule, SmilesWriteOptions::default());
    let round_trip = match written.as_ref().map_err(|_| ()).and_then(|text| {
        smiles::read_str_with_options(text, SmilesParseOptions::default()).map_err(|_| ())
    }) {
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
    let mol = molecule.graph();
    json!({
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": basic_atoms_json(mol),
        "bonds": basic_bonds_json(mol),
    })
}

pub(crate) fn smiles_sanitized_semantic_json(mut molecule: SmallMolecule) -> Value {
    match perception::sanitize_with_options(&mut molecule, SanitizeOptions::default()) {
        Ok(_) => {
            let mol = molecule.graph();
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
            let no_implicit_hydrogens = smiles_effective_no_implicit_hydrogens(atom);
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
                    "no_implicit_hydrogens": no_implicit_hydrogens,
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
    let no_implicit_hydrogens = smiles_effective_no_implicit_hydrogens(atom);
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

pub(crate) fn smiles_effective_no_implicit_hydrogens(atom: &Atom) -> bool {
    if atom.element.symbol() == "N"
        && atom.aromatic
        && atom.formal_charge == 0
        && (atom.explicit_hydrogens > 0 || atom.implicit_hydrogens == Some(1))
    {
        false
    } else {
        atom.no_implicit_hydrogens
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
    let atom_record = mol.atom(atom).ok();
    let bonds = mol
        .incident_bonds(atom)
        .ok()
        .into_iter()
        .flatten()
        .map(|(_, bond)| bond)
        .collect::<Vec<_>>();
    let has_non_aromatic_bond = bonds.iter().any(|bond| !bond.aromatic);
    let has_non_aromatic_multiple_bond = bonds.iter().any(|bond| {
        !bond.aromatic
            && matches!(
                bond.order,
                BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
            )
    });
    let aromatic_bond_count = bonds.iter().filter(|bond| bond.aromatic).count();
    let doubled: u8 = bonds
        .into_iter()
        .map(|bond| {
            if bond.aromatic {
                return aromatic_bond_valence_twice(
                    atom_record,
                    has_non_aromatic_bond,
                    has_non_aromatic_multiple_bond,
                    aromatic_bond_count,
                );
            }
            match bond.order {
                BondOrder::Zero | BondOrder::Dative => 0,
                BondOrder::Single | BondOrder::Aromatic => 2,
                BondOrder::Double => 4,
                BondOrder::Triple => 6,
                BondOrder::Quadruple => 8,
            }
        })
        .sum();
    doubled / 2
}

fn aromatic_bond_valence_twice(
    atom: Option<&Atom>,
    has_non_aromatic_bond: bool,
    has_non_aromatic_multiple_bond: bool,
    aromatic_bond_count: usize,
) -> u8 {
    let Some(atom) = atom else {
        return 2;
    };
    if atom.aromatic && has_non_aromatic_multiple_bond {
        return 2;
    }
    match atom.element.symbol() {
        "C" if atom.formal_charge < 0 && atom.explicit_hydrogens > 0 => 2,
        "O" | "S" | "Se" | "Te" if atom.formal_charge == 0 && atom.explicit_hydrogens == 0 => 2,
        "N" if atom.formal_charge < 0 => 2,
        "N" if atom.formal_charge == 0 && atom.explicit_hydrogens > 0 => 2,
        "N" if atom.formal_charge == 0 && has_non_aromatic_bond => 2,
        "N" if atom.formal_charge == 0 && aromatic_bond_count >= 3 => 2,
        _ => 3,
    }
}

pub(crate) fn bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| bond_json(id.raw(), bond, mol.stereo_bond_mark(id)))
        .collect::<Vec<_>>()
}

pub(crate) fn bond_json(index: u32, bond: &Bond, stereo: Option<&StereoBondMark>) -> Value {
    json!({
        "index": index,
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": bond.aromatic,
        "stereo": bond_stereo_json(bond.order, stereo),
        "bond_direction": bond_direction_json(bond.order, stereo),
    })
}

pub(crate) fn basic_bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| basic_bond_json(id.raw(), bond, mol.stereo_bond_mark(id)))
        .collect::<Vec<_>>()
}

pub(crate) fn basic_bond_json(index: u32, bond: &Bond, stereo: Option<&StereoBondMark>) -> Value {
    json!({
        "index": index,
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": bond.aromatic,
        "stereo": legacy_bond_stereo_json(stereo),
    })
}

pub(crate) fn legacy_bond_stereo_json(stereo: Option<&StereoBondMark>) -> &'static str {
    match stereo.map(|mark| mark.kind) {
        None => "STEREONONE",
        Some(
            StereoBondMarkKind::DirectionalUp
            | StereoBondMarkKind::DirectionalDown
            | StereoBondMarkKind::WedgeUp
            | StereoBondMarkKind::WedgeDown
            | StereoBondMarkKind::WedgeEither
            | StereoBondMarkKind::DoubleBondEither,
        ) => "STEREOANY",
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

pub(crate) fn bond_stereo_json(order: BondOrder, stereo: Option<&StereoBondMark>) -> &'static str {
    match (order, stereo.map(|mark| mark.kind)) {
        (_, None) => "STEREONONE",
        (BondOrder::Double, Some(StereoBondMarkKind::DoubleBondEither)) => "STEREOANY",
        _ => "STEREONONE",
    }
}

pub(crate) fn bond_direction_json(
    order: BondOrder,
    stereo: Option<&StereoBondMark>,
) -> &'static str {
    match (order, stereo.map(|mark| mark.kind)) {
        (
            BondOrder::Single,
            Some(StereoBondMarkKind::DirectionalUp | StereoBondMarkKind::WedgeUp),
        ) => "BEGINWEDGE",
        (
            BondOrder::Single,
            Some(StereoBondMarkKind::DirectionalDown | StereoBondMarkKind::WedgeDown),
        ) => "BEGINDASH",
        (BondOrder::Single, Some(StereoBondMarkKind::WedgeEither)) => "UNKNOWN",
        _ => "NONE",
    }
}

pub(crate) fn stereo_perception_report_json(report: &StereoPerceptionReport) -> Value {
    json!({
        "is_ok": report.is_ok(),
        "candidates": report
            .candidates
            .iter()
            .map(stereo_candidate_json)
            .collect::<Vec<_>>(),
        "issues": report
            .issues
            .iter()
            .map(stereo_perception_issue_json)
            .collect::<Vec<_>>(),
        "assembled_elements": report
            .assembled_elements
            .iter()
            .enumerate()
            .map(|(index, element)| stereo_element_json(index as u32, element))
            .collect::<Vec<_>>(),
        "created_element_indices": report
            .created_elements
            .iter()
            .map(|id| id.raw())
            .collect::<Vec<_>>(),
    })
}

pub(crate) fn stereo_candidate_json(candidate: &StereoCandidate) -> Value {
    match candidate {
        StereoCandidate::Tetrahedral { center, carriers } => json!({
            "type": "tetrahedral",
            "center_atom_index": center.raw(),
            "carriers": carriers.iter().map(stereo_carrier_json).collect::<Vec<_>>(),
        }),
        StereoCandidate::DoubleBond {
            bond,
            left,
            right,
            left_carriers,
            right_carriers,
        } => json!({
            "type": "double_bond",
            "center_bond_index": bond.raw(),
            "left_atom_index": left.raw(),
            "right_atom_index": right.raw(),
            "left_carriers": left_carriers.iter().map(stereo_carrier_json).collect::<Vec<_>>(),
            "right_carriers": right_carriers.iter().map(stereo_carrier_json).collect::<Vec<_>>(),
        }),
    }
}

pub(crate) fn stereo_perception_issue_json(issue: &StereoPerceptionIssue) -> Value {
    match issue {
        StereoPerceptionIssue::MissingStereoAtom { element, atom } => json!({
            "type": "missing_stereo_atom",
            "element_index": element.raw(),
            "atom_index": atom.raw(),
        }),
        StereoPerceptionIssue::MissingStereoBond { element, bond } => json!({
            "type": "missing_stereo_bond",
            "element_index": element.raw(),
            "bond_index": bond.raw(),
        }),
        StereoPerceptionIssue::InvalidTetrahedralCarrierCount {
            element,
            center,
            carrier_count,
        } => json!({
            "type": "invalid_tetrahedral_carrier_count",
            "element_index": element.raw(),
            "center_atom_index": center.raw(),
            "carrier_count": carrier_count,
        }),
        StereoPerceptionIssue::DuplicateTetrahedralCarrier {
            element,
            center,
            carrier,
        } => json!({
            "type": "duplicate_tetrahedral_carrier",
            "element_index": element.raw(),
            "center_atom_index": center.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoPerceptionIssue::TetrahedralCarrierNotAdjacent {
            element,
            center,
            carrier,
        } => json!({
            "type": "tetrahedral_carrier_not_adjacent",
            "element_index": element.raw(),
            "center_atom_index": center.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoPerceptionIssue::TetrahedralHydrogenCarrierUnavailable { element, center } => json!({
            "type": "tetrahedral_hydrogen_carrier_unavailable",
            "element_index": element.raw(),
            "center_atom_index": center.raw(),
        }),
        StereoPerceptionIssue::InvalidDoubleBondOrder {
            element,
            bond,
            order,
        } => json!({
            "type": "invalid_double_bond_order",
            "element_index": element.raw(),
            "bond_index": bond.raw(),
            "bond_order": bond_order_json(*order),
        }),
        StereoPerceptionIssue::DoubleBondFocusMismatch {
            element,
            bond,
            left,
            right,
        } => json!({
            "type": "double_bond_focus_mismatch",
            "element_index": element.raw(),
            "bond_index": bond.raw(),
            "left_atom_index": left.raw(),
            "right_atom_index": right.raw(),
        }),
        StereoPerceptionIssue::DoubleBondCarrierIsFocusAtom {
            element,
            endpoint,
            carrier,
        } => json!({
            "type": "double_bond_carrier_is_focus_atom",
            "element_index": element.raw(),
            "endpoint_atom_index": endpoint.raw(),
            "carrier_atom_index": carrier.raw(),
        }),
        StereoPerceptionIssue::DoubleBondCarrierNotAdjacent {
            element,
            endpoint,
            carrier,
        } => json!({
            "type": "double_bond_carrier_not_adjacent",
            "element_index": element.raw(),
            "endpoint_atom_index": endpoint.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoPerceptionIssue::DoubleBondHydrogenCarrierUnavailable { element, endpoint } => {
            json!({
                "type": "double_bond_hydrogen_carrier_unavailable",
                "element_index": element.raw(),
                "endpoint_atom_index": endpoint.raw(),
            })
        }
        StereoPerceptionIssue::UnsupportedAxisElement { element } => json!({
            "type": "unsupported_axis_element",
            "element_index": element.raw(),
        }),
        StereoPerceptionIssue::AmbiguousTetrahedralWedgeMarks { center, mark_count } => json!({
            "type": "ambiguous_tetrahedral_wedge_marks",
            "center_atom_index": center.raw(),
            "mark_count": mark_count,
        }),
        StereoPerceptionIssue::UnassembledTetrahedralBondMark { bond, kind } => json!({
            "type": "unassembled_tetrahedral_bond_mark",
            "bond_index": bond.raw(),
            "kind": stereo_bond_mark_kind_json(*kind),
        }),
        StereoPerceptionIssue::AmbiguousDirectionalBondMarks {
            double_bond,
            endpoint,
            mark_count,
        } => json!({
            "type": "ambiguous_directional_bond_marks",
            "double_bond_index": double_bond.raw(),
            "endpoint_atom_index": endpoint.raw(),
            "mark_count": mark_count,
        }),
        StereoPerceptionIssue::UnpairedDirectionalBondMark { bond } => json!({
            "type": "unpaired_directional_bond_mark",
            "bond_index": bond.raw(),
        }),
        StereoPerceptionIssue::UnsupportedSourceBondMark { bond, kind } => json!({
            "type": "unsupported_source_bond_mark",
            "bond_index": bond.raw(),
            "kind": stereo_bond_mark_kind_json(*kind),
        }),
        StereoPerceptionIssue::CouldNotCreateElement { message } => json!({
            "type": "could_not_create_element",
            "message": message,
        }),
    }
}

pub(crate) fn stereo_elements_json(mol: &Molecule) -> Vec<Value> {
    mol.stereo_elements()
        .map(|(id, element)| stereo_element_json(id.raw(), element))
        .collect()
}

pub(crate) fn stereo_element_json(index: u32, element: &StereoElement) -> Value {
    let mut object = serde_json::Map::new();
    object.insert("index".to_owned(), json!(index));
    object.insert(
        "specifiedness".to_owned(),
        json!(stereo_specifiedness_json(element.specifiedness)),
    );
    object.insert(
        "source".to_owned(),
        json!(stereo_source_json(element.source)),
    );
    if let Some(group) = element.group {
        object.insert("group_index".to_owned(), json!(group.raw()));
    }
    if let Some(descriptor) = element.descriptor {
        object.insert(
            "descriptor".to_owned(),
            json!(stereo_descriptor_json(descriptor)),
        );
    }
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            object.insert("type".to_owned(), json!("tetrahedral"));
            object.insert("center_atom_index".to_owned(), json!(stereo.center.raw()));
            object.insert(
                "carriers".to_owned(),
                Value::Array(
                    stereo
                        .carriers
                        .iter()
                        .map(stereo_carrier_json)
                        .collect::<Vec<_>>(),
                ),
            );
            object.insert(
                "orientation".to_owned(),
                json!(tetrahedral_orientation_json(stereo.orientation)),
            );
        }
        StereoElementKind::DoubleBond(stereo) => {
            object.insert("type".to_owned(), json!("double_bond"));
            object.insert("center_bond_index".to_owned(), json!(stereo.bond.raw()));
            object.insert("left_atom_index".to_owned(), json!(stereo.left.raw()));
            object.insert("right_atom_index".to_owned(), json!(stereo.right.raw()));
            object.insert(
                "left_carrier".to_owned(),
                stereo_carrier_json(&stereo.left_carrier),
            );
            object.insert(
                "right_carrier".to_owned(),
                stereo_carrier_json(&stereo.right_carrier),
            );
            object.insert(
                "orientation".to_owned(),
                json!(double_bond_orientation_json(stereo.orientation)),
            );
        }
        StereoElementKind::Axis(stereo) => {
            object.insert("type".to_owned(), json!("axis"));
            object.insert("axis_bond_index".to_owned(), json!(stereo.axis.raw()));
            object.insert(
                "carriers".to_owned(),
                Value::Array(
                    stereo
                        .carriers
                        .iter()
                        .map(stereo_carrier_json)
                        .collect::<Vec<_>>(),
                ),
            );
            object.insert(
                "orientation".to_owned(),
                json!(axis_orientation_json(stereo.orientation)),
            );
        }
    }
    Value::Object(object)
}

pub(crate) fn stereo_groups_json(mol: &Molecule) -> Vec<Value> {
    mol.stereo_groups()
        .map(|(id, group)| stereo_group_json(id.raw(), group))
        .collect()
}

pub(crate) fn stereo_group_json(index: u32, group: &StereoGroup) -> Value {
    json!({
        "index": index,
        "kind": stereo_group_kind_json(group.kind),
        "members": group.members.iter().map(|member| member.raw()).collect::<Vec<_>>(),
    })
}

pub(crate) fn stereo_bond_marks_json(mol: &Molecule) -> Vec<Value> {
    let mut marks = mol
        .stereo_bond_marks()
        .map(|mark| {
            json!({
                "bond_index": mark.bond.raw(),
                "kind": stereo_bond_mark_kind_json(mark.kind),
                "source": stereo_source_json(mark.source),
            })
        })
        .collect::<Vec<_>>();
    marks.sort_by_key(|value| {
        value
            .get("bond_index")
            .and_then(Value::as_u64)
            .unwrap_or(u64::MAX)
    });
    marks
}

pub(crate) fn stereo_carrier_json(carrier: &StereoCarrier) -> Value {
    match carrier {
        StereoCarrier::Atom(atom) => json!({ "atom_index": atom.raw() }),
        StereoCarrier::ImplicitHydrogen => json!({ "implicit_hydrogen": true }),
    }
}

pub(crate) fn stereo_specifiedness_json(specifiedness: StereoSpecifiedness) -> &'static str {
    match specifiedness {
        StereoSpecifiedness::Specified => "specified",
        StereoSpecifiedness::Unknown => "unknown",
        StereoSpecifiedness::Unspecified => "unspecified",
        StereoSpecifiedness::InvalidCleared => "invalid_cleared",
    }
}

pub(crate) fn stereo_source_json(source: StereoSource) -> &'static str {
    match source {
        StereoSource::Smiles => "smiles",
        StereoSource::MolfileV2000 => "molfile_v2000",
        StereoSource::MolfileV3000 => "molfile_v3000",
        StereoSource::Coordinates2D => "coordinates_2d",
        StereoSource::Coordinates3D => "coordinates_3d",
        StereoSource::Reaction => "reaction",
        StereoSource::User => "user",
    }
}

pub(crate) fn stereo_descriptor_json(descriptor: StereoDescriptor) -> &'static str {
    match descriptor {
        StereoDescriptor::R => "R",
        StereoDescriptor::S => "S",
        StereoDescriptor::LowerR => "r",
        StereoDescriptor::LowerS => "s",
        StereoDescriptor::E => "E",
        StereoDescriptor::Z => "Z",
        StereoDescriptor::M => "M",
        StereoDescriptor::P => "P",
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

pub(crate) fn tetrahedral_orientation_json(orientation: TetrahedralOrientation) -> &'static str {
    match orientation {
        TetrahedralOrientation::Clockwise => "clockwise",
        TetrahedralOrientation::CounterClockwise => "counter_clockwise",
    }
}

pub(crate) fn double_bond_orientation_json(orientation: DoubleBondOrientation) -> &'static str {
    match orientation {
        DoubleBondOrientation::Together => "together",
        DoubleBondOrientation::Opposite => "opposite",
    }
}

pub(crate) fn axis_orientation_json(orientation: AxisOrientation) -> &'static str {
    match orientation {
        AxisOrientation::Clockwise => "clockwise",
        AxisOrientation::CounterClockwise => "counter_clockwise",
    }
}

pub(crate) fn stereo_bond_mark_kind_json(kind: StereoBondMarkKind) -> &'static str {
    match kind {
        StereoBondMarkKind::DirectionalUp => "directional_up",
        StereoBondMarkKind::DirectionalDown => "directional_down",
        StereoBondMarkKind::WedgeUp => "wedge_up",
        StereoBondMarkKind::WedgeDown => "wedge_down",
        StereoBondMarkKind::WedgeEither => "wedge_either",
        StereoBondMarkKind::DoubleBondEither => "double_bond_either",
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
        .hierarchy()
        .atom_sites()
        .map(|(site_id, site)| {
            let residue = molecule
                .hierarchy()
                .residue(site.residue)
                .expect("residue exists");
            let chain = molecule
                .hierarchy()
                .chain(residue.chain)
                .expect("chain exists");
            let model = molecule.hierarchy().model(chain.model).expect("model exists");
            let atom = molecule.graph().atom(site.atom).expect("atom exists");
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
        "models": molecule.hierarchy().models().map(|(model_id, model)| {
            json!({
                "id": model_id.raw(),
                "chains": model.chains.iter().map(|chain_id| {
                    let chain = molecule.hierarchy().chain(*chain_id).expect("chain exists");
                    json!({
                        "id": chain.author_id.clone().unwrap_or_else(|| chain.label_id.clone()),
                        "residues": chain.residues.iter().map(|residue_id| {
                            let residue = molecule.hierarchy().residue(*residue_id).expect("residue exists");
                            json!({
                                "name": residue.name,
                                "hetflag": residue_hetflag_json(molecule, residue),
                                "sequence_id": residue_sequence_json(residue),
                                "insertion_code": residue.insertion_code,
                                "atoms": residue.atom_sites.iter().map(|site_id| {
                                    let site = molecule.hierarchy().atom_site(*site_id).expect("site exists");
                                    let atom = molecule.graph().atom(site.atom).expect("atom exists");
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

pub(crate) fn first_conformer_point(molecule: &MacroMolecule, atom: AtomId) -> Option<Point3> {
    molecule
        .graph()
        .first_conformer()
        .and_then(|(_, conformer)| conformer.position(atom))
}

pub(crate) fn residue_sequence_json(residue: &Residue) -> Value {
    if let Some(author_seq_id) = &residue.author_seq_id {
        return author_seq_id
            .parse::<i32>()
            .map(Value::from)
            .unwrap_or_else(|_| json!(author_seq_id));
    }
    residue.label_seq_id.map(Value::from).unwrap_or(Value::Null)
}

pub(crate) fn residue_hetflag_json(molecule: &MacroMolecule, residue: &Residue) -> Value {
    let is_hetatm = residue.atom_sites.iter().any(|site_id| {
        molecule
            .hierarchy()
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
