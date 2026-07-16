use crate::*;

pub(crate) fn implementation_expected(
    feature: &str,
    corpus: &str,
    fixture_path: &Path,
) -> Result<Value, Box<dyn Error>> {
    match feature {
        "bio.secondary-structure.dssp" => dssp_record_json(fixture_path),
        "io.sdf.v2000.parse" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(sdf_record_json).collect::<Vec<_>>() }))
        }
        "io.sdf.v2000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let records = records
                .into_iter()
                .map(|record| {
                    let fields = record
                        .sdf_fields
                        .into_iter()
                        .map(|(name, value)| SdfDataField::new(name, value))
                        .collect();
                    SdfRecord::new(record.title, record.molecule, fields)
                })
                .collect::<Vec<_>>();
            let written = sdf::write_v2000(&records)?;
            let records = interpret_sdf(&written)?
                .into_iter()
                .enumerate()
                .map(|(index, record)| small_record(index, record))
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
                    let title = record.title;
                    let written = molfile::write_v3000(&record.molecule)?;
                    let molecule = interpret_molfile(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title,
                        molecule,
                        sdf_fields: BTreeMap::new(),
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
                    let title = record.title;
                    let written = molfile::write_v2000(&record.molecule)?;
                    let molecule = interpret_molfile(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title,
                        molecule,
                        sdf_fields: BTreeMap::new(),
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
                    let title = record.title;
                    let written = molfile::write_v3000(&record.molecule)?;
                    let molecule = interpret_molfile(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title,
                        molecule,
                        sdf_fields: BTreeMap::new(),
                    })
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
            Ok(json!({ "records": records.iter().map(mol_record_json).collect::<Vec<_>>() }))
        }
        "io.smiles.parse" => {
            let records = read_nonisomeric_smiles_records(fixture_path)?;
            Ok(
                json!({ "records": records.iter().map(smiles_parse_record_json).collect::<Vec<_>>() }),
            )
        }
        "io.smiles.write" => {
            let records = read_nonisomeric_smiles_records(fixture_path)?;
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
        "io.smiles.isomeric" => {
            let records = read_canonical_smiles_records(fixture_path)?;
            let stereo_only = corpus != "smoke";
            Ok(json!({
                "records": records
                    .iter()
                    .filter(|record| {
                        !stereo_only || isomeric_smiles_record_is_stereo_bearing(record)
                    })
                    .map(isomeric_smiles_record_json)
                    .collect::<Result<Vec<_>, Box<dyn Error>>>()?
            }))
        }
        "query.smarts" => Ok(json!({
            "records": smarts_query_records_json(fixture_path)?
        })),
        "algo.substructure.vf2" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(substructure_record_json)
                    .collect::<Vec<_>>()
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
        "chem.hydrogen-normalization" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(hydrogen_normalization_record_json)
                    .collect::<Vec<_>>()
            }))
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
        "algo.canonical-ranking" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(canonical_ranking_record_json)
                    .collect::<Vec<_>>()
            }))
        }
        "stereo.representation" => {
            let records = read_stereo_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(stereo_record_json).collect::<Vec<_>>() }))
        }
        "stereo.perception" => {
            let mut records = read_stereo_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(stereo_perception_record_json)
                    .collect::<Vec<_>>()
            }))
        }
        "stereo.cip" => {
            let mut records = read_stereo_records_by_suffix(fixture_path)?;
            let remove_plain_hydrogens = matches!(
                fixture_path.extension().and_then(|ext| ext.to_str()),
                Some("txt" | "smi" | "smiles")
            );
            Ok(json!({
                "records": records
                    .iter_mut()
                    .filter_map(|record| stereo_cip_record_json(record, remove_plain_hydrogens))
                    .collect::<Vec<_>>()
            }))
        }
        _ => Err(boxed_error(format!(
            "no implementation comparison configured for feature `{feature}`"
        ))),
    }
}

fn dssp_record_json(fixture_path: &Path) -> Result<Value, Box<dyn Error>> {
    let input = fs::read_to_string(fixture_path)?;
    let document = mmcif::parse_str(&input, MmcifParseOptions::default())?;
    let interpretation = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            model_selection: MmcifModelSelection::First,
            ..MmcifInterpretOptions::default()
        },
    )?;
    let result = match dssp::assign(interpretation.model(), dssp::DsspOptions::default()) {
        Ok(result) => result,
        Err(dssp::DsspError::NoAnalyzableProteinResidues) => {
            return Ok(json!({
                "status": "no_analyzable_residues",
                "residues": [],
            }));
        }
        Err(error) => return Err(Box::new(error)),
    };
    let residues = result.residues().collect::<Vec<_>>();
    let residues_by_key = residues
        .iter()
        .map(|residue| (residue.key(), *residue))
        .collect::<BTreeMap<_, _>>();
    let records = residues
        .iter()
        .map(|residue| {
            let source = residue.source();
            let sequence_id = source
                .author_sequence_id
                .as_deref()
                .and_then(|value| value.parse::<i32>().ok())
                .or(source.label_sequence_id);
            json!({
                "chain_id": source.chain_author_id.as_ref().unwrap_or(&source.chain_label_id),
                "sequence_id": sequence_id,
                "insertion_code": source.insertion_code,
                "label_chain_id": source.chain_label_id,
                "label_sequence_id": source.label_sequence_id,
                "residue_name": source.residue_name,
                "residue_one_letter": dssp_residue_letter(&source.residue_name),
                "secondary_structure": residue.secondary_structure().code().to_string(),
                "chain_break": dssp_chain_break_json(residue.chain_break()),
                "phi_degrees": residue.phi_degrees(),
                "psi_degrees": residue.psi_degrees(),
                "tco": residue.tco(),
                "kappa_degrees": residue.kappa_degrees(),
                "alpha_degrees": residue.alpha_degrees(),
                "helix_positions": residue.helix_positions().map(dssp_helix_position_json),
                "sheet": residue.sheet(),
                "strand": residue.strand(),
                "ladders": residue.beta_partners().map(|partner| partner.map(|partner| partner.ladder)),
                "beta_parallel": residue.beta_partners().map(|partner| partner.map(|partner| partner.parallel)),
                "acceptors": residue.acceptors().map(|bond| dssp_bond_json(bond, &residues_by_key)),
                "donors": residue.donors().map(|bond| dssp_bond_json(bond, &residues_by_key)),
            })
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "status": "ok",
        "residues": records,
    }))
}

fn dssp_chain_break_json(chain_break: dssp::DsspChainBreak) -> &'static str {
    match chain_break {
        dssp::DsspChainBreak::None => "none",
        dssp::DsspChainBreak::NewChain => "new_chain",
        dssp::DsspChainBreak::Gap => "gap",
        _ => "unknown",
    }
}

fn dssp_helix_position_json(position: dssp::DsspHelixPosition) -> &'static str {
    match position {
        dssp::DsspHelixPosition::None => "none",
        dssp::DsspHelixPosition::Start => "start",
        dssp::DsspHelixPosition::End => "end",
        dssp::DsspHelixPosition::StartAndEnd => "start_and_end",
        dssp::DsspHelixPosition::Middle => "middle",
        _ => "unknown",
    }
}

fn dssp_bond_json(
    bond: Option<dssp::DsspHydrogenBond>,
    residues: &BTreeMap<dssp::DsspResidueKey, &dssp::DsspResidue>,
) -> Value {
    let Some(bond) = bond else {
        return Value::Null;
    };
    let source = residues[&bond.partner].source();
    json!({
        "partner_chain_id": source.chain_author_id.as_ref().unwrap_or(&source.chain_label_id),
        "partner_sequence_id": dssp_sequence_id(source),
        "partner_insertion_code": source.insertion_code,
        "energy_kcal_per_mol": bond.energy_kcal_per_mol,
    })
}

fn dssp_sequence_id(source: &dssp::DsspResidueSource) -> Option<i32> {
    source
        .author_sequence_id
        .as_deref()
        .and_then(|value| value.parse::<i32>().ok())
        .or(source.label_sequence_id)
}

fn dssp_residue_letter(name: &str) -> char {
    match name.to_ascii_uppercase().as_str() {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        _ => 'X',
    }
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedSmallRecord {
    pub(crate) record_index: usize,
    pub(crate) title: String,
    pub(crate) molecule: SmallMolecule,
    pub(crate) sdf_fields: BTreeMap<String, String>,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedSmilesRecord {
    pub(crate) record_index: usize,
    pub(crate) status: String,
    pub(crate) title: String,
    pub(crate) input_smiles: String,
    pub(crate) molecule: Option<SmallMolecule>,
}

const BOUNDED_SUBSTRUCTURE_QUERIES: &[&str] = &[
    "[#6]",
    "[!#6]",
    "A",
    "a",
    "[C,N]",
    "[C,H]",
    "[H,D]",
    "[!H]",
    "[#6]-[#8]",
    "C=O",
    "[O;H1]",
    "[#8;+0]",
    "[#6,#7;H1]",
    "[#6;R]",
    "[R0]",
    "C@C",
    "C!@C",
    "c1ccccc1",
];

fn smarts_query_records_json(path: &Path) -> Result<Vec<Value>, Box<dyn Error>> {
    let mut records = Vec::new();
    for (record_index, raw_line) in fs::read_to_string(path)?.lines().enumerate() {
        let smarts = raw_line.trim();
        if smarts.is_empty() || smarts.starts_with('#') {
            continue;
        }
        match query::parse_smarts(smarts) {
            Ok(graph) => records.push(json!({
                "record_index": record_index,
                "status": "ok",
                "smarts": smarts,
                "atom_count": graph.atom_count(),
                "bond_count": graph.bond_count(),
            })),
            Err(_) => records.push(json!({
                "record_index": record_index,
                "status": "parse_error",
                "smarts": smarts,
                "atom_count": Value::Null,
                "bond_count": Value::Null,
            })),
        }
    }
    Ok(records)
}

fn substructure_record_json(record: &mut IndexedSmallRecord) -> Value {
    if perception::sanitize(&mut record.molecule).is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
            "queries": [],
        });
    }
    let mut queries = Vec::new();
    for smarts in BOUNDED_SUBSTRUCTURE_QUERIES {
        let graph =
            query::parse_smarts(smarts).expect("checked-in bounded validation SMARTS must parse");
        let matches = substructure::find_substructure_matches(record.molecule.graph(), &graph)
            .expect("sanitized validation molecule must satisfy query prerequisites");
        let mut atom_sets = matches
            .into_iter()
            .map(|query_match| {
                let mut atoms = query_match
                    .atoms()
                    .iter()
                    .map(|atom| atom.raw())
                    .collect::<Vec<_>>();
                atoms.sort_unstable();
                atoms
            })
            .collect::<Vec<_>>();
        atom_sets.sort_unstable();
        atom_sets.dedup();
        queries.push(json!({
            "smarts": smarts,
            "matches": atom_sets,
        }));
    }
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "queries": queries,
    })
}

pub(crate) fn read_small_records_by_suffix(
    path: &Path,
) -> Result<Vec<IndexedSmallRecord>, Box<dyn Error>> {
    let input = fs::read_to_string(path)?;
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("mol" | "mdl")
    ) {
        let document = molfile::parse_str(&input)?;
        let title = document.header().title().to_owned();
        let molecule = molfile::interpret(&document)?.into_molecule();
        return Ok(vec![IndexedSmallRecord {
            record_index: 0,
            title,
            molecule,
            sdf_fields: BTreeMap::new(),
        }]);
    }
    Ok(interpret_sdf(&input)?
        .into_iter()
        .enumerate()
        .map(|(index, record)| small_record(index, record))
        .collect())
}

pub(crate) fn small_record(index: usize, record: SdfRecord) -> IndexedSmallRecord {
    let title = record.title().to_owned();
    let sdf_fields = record
        .data_fields()
        .iter()
        .map(|field| (field.name().to_owned(), field.value().to_owned()))
        .collect();
    IndexedSmallRecord {
        record_index: index,
        title,
        molecule: record.into_molecule(),
        sdf_fields,
    }
}

fn interpret_molfile(input: &str) -> Result<SmallMolecule, Box<dyn Error>> {
    let document = molfile::parse_str(input)?;
    Ok(molfile::interpret(&document)?.into_molecule())
}

fn interpret_sdf(input: &str) -> Result<Vec<SdfRecord>, Box<dyn Error>> {
    let document = sdf::parse_str(input, SdfParseOptions::default())?;
    Ok(sdf::interpret(&document)?.into_records())
}

fn interpret_smiles(input: &str) -> Result<SmallMolecule, Box<dyn Error>> {
    let document = smiles::parse_str(input)?;
    Ok(smiles::interpret(&document)?.into_molecule())
}

pub(crate) fn read_smiles_records(path: &Path) -> Result<Vec<IndexedSmilesRecord>, Box<dyn Error>> {
    read_smiles_records_with_filter(path, |smiles| smiles.contains('*'))
}

pub(crate) fn read_nonisomeric_smiles_records(
    path: &Path,
) -> Result<Vec<IndexedSmilesRecord>, Box<dyn Error>> {
    read_smiles_records_with_filter(path, |smiles| {
        smiles_unsupported_subset_reason(smiles).is_some()
    })
}

fn read_smiles_records_with_filter(
    path: &Path,
    unsupported: impl Fn(&str) -> bool,
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
        if unsupported(&smiles) {
            records.push(IndexedSmilesRecord {
                record_index: index,
                status: "unsupported".to_owned(),
                title,
                input_smiles: smiles,
                molecule: None,
            });
            continue;
        }
        let (status, molecule) = match interpret_smiles(&smiles) {
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
                    sdf_fields: BTreeMap::new(),
                })
            })
            .collect());
    }
    if !matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("mol" | "mdl")
    ) {
        return read_small_records_by_suffix(path);
    }
    let document = molfile::parse_str(&input)?;
    let title = document.header().title().to_owned();
    let molecule = molfile::interpret(&document)?.into_molecule();
    Ok(vec![IndexedSmallRecord {
        record_index: 0,
        title,
        molecule,
        sdf_fields: BTreeMap::new(),
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
        let (status, molecule) = match interpret_smiles(&smiles) {
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
    let mol = record.molecule.graph();
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": atoms_json(mol),
        "bonds": bonds_json(mol),
        "properties": record.sdf_fields,
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
                    conformer.position(atom_id).map(|point| {
                        let point = point.value();
                        json!({
                            "atom_index": atom_id.raw(),
                            "x": point.x,
                            "y": point.y,
                            "z": point.z,
                        })
                    })
                })
                .collect::<Vec<_>>()
        }).collect::<Vec<_>>(),
        "atoms": mol.atoms().map(|(id, atom)| conformer_atom_json(mol, id, atom)).collect::<Vec<_>>(),
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
    let sanitize = perception::sanitize_with_options(
        &mut record.molecule,
        SanitizeOptions {
            perceive_stereo: false,
            ..SanitizeOptions::default()
        },
    );
    if sanitize.is_err() {
        let mol = record.molecule.graph();
        return json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
            "atom_count": mol.atom_count(),
            "bond_count": mol.bond_count(),
        });
    }
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

pub(crate) fn stereo_cip_record_json(
    record: &mut IndexedSmallRecord,
    remove_plain_hydrogens: bool,
) -> Option<Value> {
    let sanitize = perception::sanitize_with_options(
        &mut record.molecule,
        SanitizeOptions {
            perceive_stereo: false,
            ..SanitizeOptions::default()
        },
    );
    if sanitize.is_err() {
        return None;
    }
    let perception_report = stereo::perceive_stereo_with_options(
        record.molecule.graph_mut(),
        stereo::StereoPerceptionOptions {
            assign_coordinates: false,
            ..stereo::StereoPerceptionOptions::default()
        },
    );
    if !perception_report.is_ok() {
        return None;
    }
    stereo::assign_cip_descriptors(record.molecule.graph_mut());
    let mol = record.molecule.graph();
    let atom_index = rdkit_default_atom_index(mol, remove_plain_hydrogens);
    let atom_descriptors = cip_atom_descriptors_json(mol, &atom_index);
    let bond_descriptors = cip_bond_descriptors_json(mol, &atom_index);
    if atom_descriptors.is_empty() && bond_descriptors.is_empty() {
        return None;
    }
    Some(json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": atom_index.len(),
        "bond_count": rdkit_default_bond_count(mol, &atom_index),
        "atom_descriptors": atom_descriptors,
        "bond_descriptors": bond_descriptors,
    }))
}

pub(crate) fn cip_atom_descriptors_json(
    mol: &Molecule,
    atom_index: &BTreeMap<AtomId, u32>,
) -> Vec<Value> {
    let mut descriptors = mol
        .stereo_elements()
        .filter_map(|(id, element)| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => mol
                .cip_descriptor(id)
                .ok()
                .flatten()
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
    atom_index: &BTreeMap<AtomId, u32>,
) -> Vec<Value> {
    let mut descriptors = mol
        .stereo_elements()
        .filter_map(|(id, element)| match &element.kind {
            StereoElementKind::DoubleBond(stereo) => mol
                .cip_descriptor(id)
                .ok()
                .flatten()
                .and_then(|descriptor| {
                    let begin_atom_index = *atom_index.get(&stereo.left)?;
                    let end_atom_index = *atom_index.get(&stereo.right)?;
                    Some(json!({
                        "begin_atom_index": begin_atom_index,
                        "end_atom_index": end_atom_index,
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
                        let begin_atom_index = *atom_index.get(&begin)?;
                        let end_atom_index = *atom_index.get(&end)?;
                        Some(json!({
                            "begin_atom_index": begin_atom_index,
                            "end_atom_index": end_atom_index,
                            "descriptor": stereo_descriptor_json(descriptor),
                        }))
                    })
            }
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

fn rdkit_default_atom_index(mol: &Molecule, remove_plain_hydrogens: bool) -> BTreeMap<AtomId, u32> {
    let mut index = BTreeMap::new();
    for (atom_id, atom) in mol.atoms() {
        if remove_plain_hydrogens && rdkit_default_removes_hydrogen(atom) {
            continue;
        }
        index.insert(atom_id, index.len() as u32);
    }
    index
}

fn rdkit_default_bond_count(mol: &Molecule, atom_index: &BTreeMap<AtomId, u32>) -> usize {
    mol.bonds()
        .filter(|(_, bond)| {
            atom_index.contains_key(&bond.a()) && atom_index.contains_key(&bond.b())
        })
        .count()
}

fn rdkit_default_removes_hydrogen(atom: &Atom) -> bool {
    atom.element.symbol() == "H"
        && atom.isotope.is_none()
        && atom.formal_charge == 0
        && atom.radical.is_none()
        && atom.atom_map.is_none()
        && atom.props.is_empty()
}

pub(crate) fn conformers_json(mol: &Molecule) -> Vec<Vec<Value>> {
    mol.conformers()
        .map(|(_, conformer)| {
            mol.atom_ids()
                .filter_map(|atom_id| {
                    conformer.position(atom_id).map(|point| {
                        let point = point.value();
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

pub(crate) fn conformer_atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "atom_map": atom.atom_map,
        "aromatic": mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false),
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
    let report = valence::perceive_valence_with_options(
        record.molecule.graph_mut(),
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    );
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

pub(crate) fn hydrogen_normalization_record_json(record: &mut IndexedSmallRecord) -> Value {
    if perception::sanitize_with_options(&mut record.molecule, SanitizeOptions::default()).is_err()
    {
        return json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
        });
    }
    let added = match hydrogens::add_hydrogens(&mut record.molecule) {
        Ok(report) => report,
        Err(_) => {
            return json!({
                "record_index": record.record_index,
                "status": "add_error",
                "title": record.title,
            });
        }
    };
    let atom_count_after_add = record.molecule.atom_count();
    let mut added_by_parent = BTreeMap::<usize, usize>::new();
    for entry in added.added {
        *added_by_parent.entry(entry.parent.index()).or_default() += 1;
    }

    if !valence::perceive_valence(record.molecule.graph_mut(), ValenceModel::RdkitLike).is_ok() {
        return json!({
            "record_index": record.record_index,
            "status": "add_error",
            "title": record.title,
        });
    }
    if hydrogens::remove_hydrogens(&mut record.molecule).is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "remove_error",
            "title": record.title,
        });
    }

    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count_after_add": atom_count_after_add,
        "added_hydrogens_by_parent": added_by_parent
            .into_iter()
            .map(|(parent_atom_index, count)| json!({
                "parent_atom_index": parent_atom_index,
                "count": count,
            }))
            .collect::<Vec<_>>(),
        "round_trip": hydrogen_normalized_semantic_json(record.molecule.clone()),
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
        "atom_aromatic": mol.atoms().map(|(id, _)| mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false)).collect::<Vec<_>>(),
        "bond_aromatic": mol.bonds().map(|(id, _)| mol.bond_is_aromatic(id).ok().flatten().unwrap_or(false)).collect::<Vec<_>>(),
    })
}

pub(crate) fn canonical_ranking_record_json(record: &mut IndexedSmallRecord) -> Value {
    let options = SanitizeOptions {
        perceive_valence: true,
        perceive_rings: true,
        perceive_aromaticity: true,
        perceive_stereo: false,
    };
    if perception::sanitize_with_options(&mut record.molecule, options).is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
        });
    }
    let ranking = canon::atom_ranking(record.molecule.graph());
    let mut classes = BTreeMap::<u32, Vec<usize>>::new();
    for (atom, rank) in ranking.iter() {
        classes.entry(rank).or_default().push(atom.index());
    }
    let mut classes = classes.into_values().collect::<Vec<_>>();
    classes.sort();
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "classes": classes,
    })
}

pub(crate) fn smiles_write_record_json(
    record: &IndexedSmilesRecord,
) -> Result<Value, Box<dyn Error>> {
    let Some(molecule) = &record.molecule else {
        return Ok(smiles_error_record_json(record));
    };
    let written = smiles::write_with_options(molecule, SmilesWriteOptions::default())?;
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
            "status": "parse_error",
            "title": record.title,
            "input_smiles": record.input_smiles,
        }));
    }
    let written =
        smiles::write_canonical_with_options(&molecule, CanonicalSmilesWriteOptions::default())?;
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
        "sanitized": smiles_sanitized_semantic_json(reparsed),
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
    if perception::sanitize_with_options(&mut molecule, SanitizeOptions::default()).is_err() {
        return Ok(json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
            "input_smiles": record.input_smiles,
        }));
    }
    let written = match smiles::write_isomeric_with_options(&molecule, Default::default()) {
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
        "sanitized": smiles_sanitized_semantic_json(reparsed.clone()),
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
    perception::sanitize_with_options(&mut molecule, SanitizeOptions::default()).is_ok()
}

pub(crate) fn smiles_parse_record_json(record: &IndexedSmilesRecord) -> Value {
    let Some(molecule) = &record.molecule else {
        return smiles_error_record_json(record);
    };
    let written = smiles::write_with_options(molecule, SmilesWriteOptions::default());
    let round_trip = match written
        .as_ref()
        .map_err(|_| ())
        .and_then(|text| interpret_smiles(text).map_err(|_| ()))
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

pub(crate) fn hydrogen_normalized_semantic_json(mut molecule: SmallMolecule) -> Value {
    let _ = valence::perceive_valence_with_options(
        molecule.graph_mut(),
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    );
    let mol = molecule.graph();
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
                "encoded_hydrogens": usize::from(atom.explicit_hydrogens)
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

pub(crate) fn smiles_isomeric_stereo_semantic_json(mut molecule: SmallMolecule) -> Value {
    if perception::sanitize_with_options(&mut molecule, SanitizeOptions::default()).is_err() {
        return json!({ "status": "sanitize_error" });
    }
    stereo::assign_cip_descriptors(molecule.graph_mut());
    let mol = molecule.graph();
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
                        "center_atom": smiles_sanitized_atom_key(mol, stereo.center, atom),
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
                        smiles_sanitized_atom_key(mol, stereo.left, left),
                        smiles_sanitized_atom_key(mol, stereo.right, right),
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
                            smiles_sanitized_atom_key(mol, begin, begin_atom),
                            smiles_sanitized_atom_key(mol, end, end_atom),
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

pub(crate) fn smiles_sanitized_bonds_json(mol: &Molecule) -> Vec<Value> {
    let mut bonds = mol
        .bonds()
        .map(|(bond_id, bond)| {
            let left = mol.atom(bond.a()).expect("bond endpoint should exist");
            let right = mol.atom(bond.b()).expect("bond endpoint should exist");
            let mut endpoints = [
                smiles_sanitized_atom_key(mol, bond.a(), left),
                smiles_sanitized_atom_key(mol, bond.b(), right),
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

pub(crate) fn smiles_sanitized_atoms_json(mol: &Molecule) -> Vec<Value> {
    let mut atoms = mol
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
                        "atom": smiles_sanitized_atom_key(mol, neighbor_id, neighbor),
                        "bond_type": smiles_semantic_bond_type(mol, bond_id, bond),
                        "is_aromatic": mol.bond_is_aromatic(bond_id).ok().flatten().unwrap_or(false),
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
                    "aromatic": mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false),
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
    if atom.element.symbol() == "N"
        && mol.atom_is_aromatic(id).ok().flatten() == Some(true)
        && atom.explicit_hydrogens == 0
        && implicit == 1
    {
        (1, 0)
    } else {
        (atom.explicit_hydrogens, implicit)
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
        && (atom.explicit_hydrogens > 0 || mol.implicit_hydrogens(id).ok().flatten() == Some(1))
    {
        false
    } else {
        atom.no_implicit_hydrogens
    }
}

pub(crate) fn atoms_json(mol: &Molecule) -> Vec<Value> {
    mol.atoms()
        .map(|(id, atom)| atom_json(mol, id, atom))
        .collect::<Vec<_>>()
}

pub(crate) fn atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
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
        "aromatic": mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false),
    })
}

pub(crate) fn basic_atoms_json(mol: &Molecule) -> Vec<Value> {
    mol.atoms()
        .map(|(id, atom)| basic_atom_json(mol, id, atom))
        .collect::<Vec<_>>()
}

pub(crate) fn basic_atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "atom_map": atom.atom_map,
        "aromatic": mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false),
    })
}

pub(crate) fn valence_atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "implicit_hydrogens": mol.implicit_hydrogens(id).ok().flatten().unwrap_or(0),
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
        .collect::<Vec<_>>();
    let has_non_aromatic_bond = bonds
        .iter()
        .any(|(id, _)| mol.bond_is_aromatic(*id).ok().flatten() != Some(true));
    let has_non_aromatic_multiple_bond = bonds.iter().any(|(id, bond)| {
        mol.bond_is_aromatic(*id).ok().flatten() != Some(true)
            && matches!(
                bond.order,
                BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
            )
    });
    let has_marked_aromatic_high_order_bond = bonds.iter().any(|(id, bond)| {
        mol.bond_is_aromatic(*id).ok().flatten() == Some(true)
            && matches!(bond.order, BondOrder::Triple | BondOrder::Quadruple)
    });
    let aromatic_bond_count = bonds
        .iter()
        .filter(|(id, _)| mol.bond_is_aromatic(*id).ok().flatten() == Some(true))
        .count();
    let doubled: u8 = bonds
        .into_iter()
        .map(|(id, bond)| {
            if mol.bond_is_aromatic(id).ok().flatten() == Some(true) {
                if has_marked_aromatic_high_order_bond {
                    return match bond.order {
                        BondOrder::Triple => 6,
                        BondOrder::Quadruple => 8,
                        _ => 2,
                    };
                }
                return aromatic_bond_valence_twice(
                    atom_record,
                    mol.atom_is_aromatic(atom).ok().flatten() == Some(true),
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
    atom_aromatic: bool,
    has_non_aromatic_bond: bool,
    has_non_aromatic_multiple_bond: bool,
    aromatic_bond_count: usize,
) -> u8 {
    let Some(atom) = atom else {
        return 2;
    };
    if atom_aromatic && has_non_aromatic_multiple_bond {
        return 2;
    }
    match atom.element.symbol() {
        "C" if atom.formal_charge < 0
            && (atom.explicit_hydrogens > 0
                || has_non_aromatic_bond
                || aromatic_bond_count >= 3) =>
        {
            2
        }
        "P" | "As" | "Sb"
            if atom.formal_charge == 0
                && atom.explicit_hydrogens == 0
                && (has_non_aromatic_bond || aromatic_bond_count >= 3) =>
        {
            2
        }
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
        .map(|(id, bond)| bond_json(mol, id, bond, mol.stereo_bond_mark(id)))
        .collect::<Vec<_>>()
}

pub(crate) fn bond_json(
    mol: &Molecule,
    id: BondId,
    bond: &Bond,
    stereo: Option<&StereoBondMark>,
) -> Value {
    json!({
        "index": id.raw(),
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": mol.bond_is_aromatic(id).ok().flatten().unwrap_or(false),
        "stereo": bond_stereo_json(bond.order, stereo),
        "bond_direction": bond_direction_json(bond.order, stereo),
    })
}

pub(crate) fn basic_bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| basic_bond_json(mol, id, bond, mol.stereo_bond_mark(id)))
        .collect::<Vec<_>>()
}

pub(crate) fn basic_bond_json(
    mol: &Molecule,
    id: BondId,
    bond: &Bond,
    stereo: Option<&StereoBondMark>,
) -> Value {
    json!({
        "index": id.raw(),
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": mol.bond_is_aromatic(id).ok().flatten().unwrap_or(false),
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
        AtomRadical::Quartet => "QUARTET",
        AtomRadical::Quintet => "QUINTET",
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
            .map(|(index, element)| stereo_element_json(index as u32, element, None))
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
        StereoPerceptionIssue::InvalidAxisCarrierCount {
            element,
            axis,
            carrier_count,
        } => json!({
            "type": "invalid_axis_carrier_count",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier_count": carrier_count,
        }),
        StereoPerceptionIssue::AxisCarrierIsFocusAtom {
            element,
            axis,
            carrier,
        } => json!({
            "type": "axis_carrier_is_focus_atom",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier_atom_index": carrier.raw(),
        }),
        StereoPerceptionIssue::AxisCarrierNotAdjacent {
            element,
            axis,
            carrier,
        } => json!({
            "type": "axis_carrier_not_adjacent",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier": stereo_carrier_json(carrier),
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
        .map(|(id, element)| {
            stereo_element_json(id.raw(), element, mol.cip_descriptor(id).ok().flatten())
        })
        .collect()
}

pub(crate) fn stereo_element_json(
    index: u32,
    element: &StereoElement,
    descriptor: Option<StereoDescriptor>,
) -> Value {
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
    if let Some(descriptor) = descriptor {
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
        StereoCarrier::ImplicitLonePair => json!({ "implicit_lone_pair": true }),
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
        StereoDescriptor::SeqTrans => "seqTrans",
        StereoDescriptor::SeqCis => "seqCis",
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
