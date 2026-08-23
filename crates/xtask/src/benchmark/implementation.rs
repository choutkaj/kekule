use crate::*;

pub(crate) fn implementation_expected(
    benchmark: &str,
    _corpus: &str,
    fixture_path: &Path,
) -> Result<Value, Box<dyn Error>> {
    match benchmark {
        "bio.secondary-structure.dssp" => dssp_record_json(fixture_path),
        "io.mmcif.parse" => mmcif_document_json(fixture_path),
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
                    SdfRecordInterpretation::new(record.title, vec![record.molecule], fields)
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
        "descriptor.molecular" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(molecular_descriptor_record_json)
                    .collect::<Vec<_>>()
            }))
        }
        "descriptor.rotatable-bonds.rdkit-strict" => {
            if matches!(
                fixture_path
                    .extension()
                    .and_then(|extension| extension.to_str()),
                Some("smi" | "smiles" | "txt")
            ) {
                let records = read_smiles_records(fixture_path)?;
                return Ok(json!({
                    "records": records
                        .iter()
                        .map(rotatable_bond_smiles_record_json)
                        .collect::<Vec<_>>()
                }));
            }
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter()
                    .map(rotatable_bond_record_json)
                    .collect::<Vec<_>>()
            }))
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
            let exact_smiles = false;
            Ok(json!({
                "records": records
                    .iter()
                    .map(|record| canonical_smiles_record_json(record, exact_smiles))
                    .collect::<Result<Vec<_>, Box<dyn Error>>>()?
            }))
        }
        "io.smiles.isomeric" => {
            let records = read_canonical_smiles_records(fixture_path)?;
            let stereo_only = true;
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
        "chem.hydrogen-transforms" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(hydrogen_transform_record_json)
                    .collect::<Vec<_>>()
            }))
        }
        "chem.perception.default" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(default_perception_atom_record_json).collect::<Vec<_>>() }),
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
            let mut records = read_stereo_perception_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(stereo_perception_group_record_json)
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
            "no implementation comparison configured for benchmark `{benchmark}`"
        ))),
    }
}

const MMCIF_ATOM_SITE_FIELDS: &[&str] = &[
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "auth_atom_id",
    "label_alt_id",
    "label_comp_id",
    "auth_comp_id",
    "label_asym_id",
    "auth_asym_id",
    "label_seq_id",
    "auth_seq_id",
    "pdbx_PDB_ins_code",
    "occupancy",
    "B_iso_or_equiv",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "pdbx_PDB_model_num",
];

fn mmcif_document_json(fixture_path: &Path) -> Result<Value, Box<dyn Error>> {
    let input = fs::read_to_string(fixture_path)?;
    let document = mmcif::parse_str(&input, MmcifParseOptions::default())?;
    let table = document
        .blocks()
        .iter()
        .find_map(|block| block.loop_with_tag("_atom_site.id"));
    let row_count = table.map_or(0, |table| table.row_count());
    let rows = (0..row_count)
        .map(|row_index| {
            let mut row = serde_json::Map::new();
            for field in MMCIF_ATOM_SITE_FIELDS {
                let tag = format!("_atom_site.{field}");
                let value = table
                    .and_then(|table| table.value(row_index, &tag))
                    .and_then(|value| value.optional_text());
                row.insert((*field).to_owned(), json!(value));
            }
            Value::Object(row)
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "atom_site_rows": {
            "status": "ok",
            "row_count": row_count,
            "rows": rows,
        }
    }))
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
    let result = match dssp::assign(interpretation.model().view(), dssp::DsspOptions::default()) {
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
    residues: &BTreeMap<kekule::topology::InstanceResidueId, &dssp::DsspResidue>,
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
    pub(crate) molecule: Molecule,
    pub(crate) sdf_fields: BTreeMap<String, String>,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedSmilesRecord {
    pub(crate) record_index: usize,
    pub(crate) status: String,
    pub(crate) title: String,
    pub(crate) input_smiles: String,
    pub(crate) molecule: Option<Molecule>,
    pub(crate) components: Vec<Molecule>,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedStereoPerceptionRecord {
    pub(crate) record_index: usize,
    pub(crate) title: String,
    pub(crate) components: Vec<Molecule>,
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
    if record.molecule.perceive().is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "perception_error",
            "title": record.title,
            "queries": [],
        });
    }
    let mut queries = Vec::new();
    for smarts in BOUNDED_SUBSTRUCTURE_QUERIES {
        let graph =
            query::parse_smarts(smarts).expect("checked-in bounded benchmark SMARTS must parse");
        let matches = substructure::find_substructure_matches(&record.molecule, &graph)
            .expect("perceived benchmark molecule must satisfy query prerequisites");
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
        let molecule = exactly_one_molecule(molfile::interpret(&document)?.to_molecules())?;
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

pub(crate) fn small_record(index: usize, record: SdfRecordInterpretation) -> IndexedSmallRecord {
    let title = record.title().to_owned();
    let sdf_fields = record
        .data_fields()
        .iter()
        .map(|field| (field.name().to_owned(), field.value().to_owned()))
        .collect();
    IndexedSmallRecord {
        record_index: index,
        title,
        molecule: exactly_one_molecule(record.to_molecules())
            .expect("small-record benchmark requires one connected component"),
        sdf_fields,
    }
}

fn interpret_molfile(input: &str) -> Result<Molecule, Box<dyn Error>> {
    let document = molfile::parse_str(input)?;
    exactly_one_molecule(molfile::interpret(&document)?.to_molecules())
}

fn interpret_sdf(input: &str) -> Result<Vec<SdfRecordInterpretation>, Box<dyn Error>> {
    let document = sdf::parse_str(input, SdfParseOptions::default())?;
    Ok(sdf::interpret(&document)?.to_records())
}

fn interpret_smiles(input: &str) -> Result<Molecule, Box<dyn Error>> {
    let document = smiles::parse_str(input)?;
    exactly_one_molecule(smiles::interpret(&document)?.to_molecules())
}

fn exactly_one_molecule(mut molecules: Vec<Molecule>) -> Result<Molecule, Box<dyn Error>> {
    if molecules.len() != 1 {
        return Err(std::io::Error::other(format!(
            "expected one connected molecule, found {}",
            molecules.len()
        ))
        .into());
    }
    Ok(molecules.pop().expect("component count was checked"))
}

fn interpret_smiles_components(input: &str) -> Result<Vec<Molecule>, Box<dyn Error>> {
    let document = smiles::parse_str(input)?;
    Ok(smiles::interpret(&document)?.to_molecules())
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
                components: Vec::new(),
            });
            continue;
        }
        let (status, molecule, components) = match interpret_smiles_components(&smiles) {
            Ok(components) => {
                let molecule = (components.len() == 1).then(|| components[0].clone());
                ("ok".to_owned(), molecule, components)
            }
            Err(_) => ("parse_error".to_owned(), None, Vec::new()),
        };
        records.push(IndexedSmilesRecord {
            record_index: index,
            status,
            title,
            input_smiles: smiles,
            molecule,
            components,
        });
    }
    Ok(records)
}

pub(crate) fn read_stereo_perception_records_by_suffix(
    path: &Path,
) -> Result<Vec<IndexedStereoPerceptionRecord>, Box<dyn Error>> {
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("txt" | "smi" | "smiles")
    ) {
        return Ok(read_smiles_records(path)?
            .into_iter()
            .map(|record| IndexedStereoPerceptionRecord {
                record_index: record.record_index,
                title: record.title,
                components: record.components,
            })
            .collect());
    }
    Ok(read_stereo_records_by_suffix(path)?
        .into_iter()
        .map(|record| IndexedStereoPerceptionRecord {
            record_index: record.record_index,
            title: record.title,
            components: vec![record.molecule],
        })
        .collect())
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
    let molecule = exactly_one_molecule(molfile::interpret(&document)?.to_molecules())?;
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
            components: molecule.iter().cloned().collect(),
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
    let mol = &record.molecule;
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
    let mol = &record.molecule;
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
    let mol = &record.molecule;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "conformers": Vec::<Value>::new(),
        "atoms": mol.atoms().map(|(id, atom)| conformer_atom_json(mol, id, atom)).collect::<Vec<_>>(),
    })
}

pub(crate) fn molecular_descriptor_record_json(record: &mut IndexedSmallRecord) -> Value {
    if record.molecule.perceive().is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "perception_error",
            "title": record.title,
        });
    }
    let policy = kekule::descriptors::HydrogenCountPolicy::IncludePerceived;
    let result = (|| {
        let formula = kekule::descriptors::molecular_formula(&record.molecule, policy)?;
        let average = kekule::descriptors::average_mass(&record.molecule, policy)?;
        let monoisotopic = kekule::descriptors::monoisotopic_mass(&record.molecule, policy)?;
        Ok::<_, kekule::descriptors::MolecularDescriptorError>((
            formula,
            *average.value(),
            *monoisotopic.value(),
        ))
    })();
    let Ok((formula, average_mass_da, monoisotopic_mass_da)) = result else {
        return json!({
            "record_index": record.record_index,
            "status": "descriptor_error",
            "title": record.title,
        });
    };
    let terms = formula
        .terms()
        .map(|(element, isotope, count)| {
            json!({
                "element": element.symbol(),
                "isotope": isotope,
                "count": count,
            })
        })
        .collect::<Vec<_>>();
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "formula": {
            "terms": terms,
            "formal_charge": formula.formal_charge(),
        },
        "average_mass_da": average_mass_da,
        "monoisotopic_mass_da": monoisotopic_mass_da,
    })
}

pub(crate) fn rotatable_bond_record_json(record: &IndexedSmallRecord) -> Value {
    let molecule = &record.molecule;
    let detected = kekule::rotatable_bonds::detect(
        molecule,
        kekule::rotatable_bonds::RotatableBondOptions::STRICT,
    );
    let bonds = detected
        .bond_ids()
        .iter()
        .copied()
        .map(|bond_id| {
            let bond = molecule
                .bond(bond_id)
                .expect("rotatable-bond detector returns live bond IDs");
            json!({
                "begin_atom_index": bond.a().raw(),
                "end_atom_index": bond.b().raw(),
            })
        })
        .collect::<Vec<_>>();
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "count": detected.len(),
        "bonds": bonds,
    })
}

pub(crate) fn rotatable_bond_smiles_record_json(record: &IndexedSmilesRecord) -> Value {
    if record.status != "ok" {
        return json!({
            "record_index": record.record_index,
            "status": record.status,
            "title": record.title,
        });
    }

    let mut atom_offset = 0usize;
    let mut bonds = Vec::new();
    for component in &record.components {
        let molecule = component;
        let detected = kekule::rotatable_bonds::detect(
            molecule,
            kekule::rotatable_bonds::RotatableBondOptions::STRICT,
        );
        bonds.extend(detected.bond_ids().iter().copied().map(|bond_id| {
            let bond = molecule
                .bond(bond_id)
                .expect("rotatable-bond detector returns live bond IDs");
            json!({
                "begin_atom_index": atom_offset + bond.a().raw() as usize,
                "end_atom_index": atom_offset + bond.b().raw() as usize,
            })
        }));
        atom_offset += molecule.atom_count();
    }
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "count": bonds.len(),
        "bonds": bonds,
    })
}

pub(crate) fn mol_parse_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule;
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
    let mol = &record.molecule;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "stereo_elements": stereo_elements_json(mol),
        "stereo_groups": stereo_groups_json(mol),
    })
}

pub(crate) fn stereo_perception_record_json(record: &mut IndexedSmallRecord) -> Value {
    let source_stereo_elements = record
        .molecule
        .stereo_elements()
        .map(|(id, _)| id)
        .collect::<Vec<_>>();
    if record.molecule.perceive().is_err() {
        let mol = &record.molecule;
        return json!({
            "record_index": record.record_index,
            "status": "perception_error",
            "title": record.title,
            "atom_count": mol.atom_count(),
            "bond_count": mol.bond_count(),
        });
    }
    let candidates = stereo::detect_stereo_candidates(&record.molecule);
    let conformer =
        kekule::core::Conformer::new(kekule::units::ANGSTROM).expect("angstrom is a length unit");
    let mut editor = record.molecule.edit();
    let result = stereo::materialize_coordinate_stereo(&mut editor, &conformer);
    if result.is_ok() {
        record.molecule = editor
            .finish()
            .expect("coordinate stereo materialization preserves publication invariants");
    }
    let mol = &record.molecule;
    match result {
        Ok(report) => json!({
            "record_index": record.record_index,
            "status": "ok",
            "title": record.title,
            "atom_count": mol.atom_count(),
            "bond_count": mol.bond_count(),
            "report": stereo_perception_benchmark_report_json(
                mol,
                &source_stereo_elements,
                &candidates,
                &report,
            ),
            "stereo_elements": stereo_elements_json(mol),
            "stereo_groups": stereo_groups_json(mol),
        }),
        Err(error) => json!({
            "record_index": record.record_index,
            "status": "perception_error",
            "title": record.title,
            "atom_count": mol.atom_count(),
            "bond_count": mol.bond_count(),
            "candidates": candidates.iter().map(stereo_candidate_json).collect::<Vec<_>>(),
            "source_stereo_element_indices": source_stereo_elements
                .iter()
                .map(|id| id.raw())
                .collect::<Vec<_>>(),
            "error": coordinate_stereo_error_json(&error),
            "stereo_elements": stereo_elements_json(mol),
            "stereo_groups": stereo_groups_json(mol),
        }),
    }
}

pub(crate) fn stereo_perception_group_record_json(
    record: &mut IndexedStereoPerceptionRecord,
) -> Value {
    if record.components.is_empty() {
        return json!({
            "record_index": record.record_index,
            "status": "parse_error",
            "title": record.title,
            "atom_count": 0,
            "bond_count": 0,
        });
    }

    let mut atom_count = 0u64;
    let mut bond_count = 0u64;
    let mut element_count = 0u64;
    let mut group_count = 0u64;
    let mut assembled_count = 0u64;
    let mut assembled_elements = Vec::new();
    let mut candidates = Vec::new();
    let mut created_element_indices = Vec::new();
    let mut issues = Vec::new();
    let mut stereo_elements = Vec::new();
    let mut stereo_groups = Vec::new();

    for component in &record.components {
        let mut component_record = IndexedSmallRecord {
            record_index: record.record_index,
            title: record.title.clone(),
            molecule: component.clone(),
            sdf_fields: BTreeMap::new(),
        };
        let mut value = stereo_perception_record_json(&mut component_record);
        if value.get("status").and_then(Value::as_str) != Some("ok") {
            return json!({
                "record_index": record.record_index,
                "status": value.get("status").cloned().unwrap_or_else(|| json!("perception_error")),
                "title": record.title,
                "atom_count": record.components.iter().map(|molecule| molecule.atom_count()).sum::<usize>(),
                "bond_count": record.components.iter().map(|molecule| molecule.bond_count()).sum::<usize>(),
            });
        }
        let object = value
            .as_object_mut()
            .expect("stereo record must be an object");
        let component_atom_count = object
            .get("atom_count")
            .and_then(Value::as_u64)
            .unwrap_or_default();
        let component_bond_count = object
            .get("bond_count")
            .and_then(Value::as_u64)
            .unwrap_or_default();
        let mut report = match object.remove("report") {
            Some(Value::Object(report)) => report,
            _ => panic!("successful stereo record must contain a report"),
        };

        let mut component_assembled = take_array(&mut report, "assembled_elements");
        for element in &mut component_assembled {
            offset_object_u64(element, "index", assembled_count);
            offset_stereo_references(element, atom_count, bond_count, element_count, group_count);
        }
        assembled_count += component_assembled.len() as u64;
        assembled_elements.extend(component_assembled);

        let mut component_candidates = take_array(&mut report, "candidates");
        for candidate in &mut component_candidates {
            offset_stereo_references(
                candidate,
                atom_count,
                bond_count,
                element_count,
                group_count,
            );
        }
        candidates.extend(component_candidates);

        for created in take_array(&mut report, "created_element_indices") {
            if let Some(index) = created.as_u64() {
                created_element_indices.push(json!(index + element_count));
            }
        }
        let mut component_issues = take_array(&mut report, "issues");
        for issue in &mut component_issues {
            offset_stereo_references(issue, atom_count, bond_count, element_count, group_count);
        }
        issues.extend(component_issues);

        let mut component_elements = take_array(object, "stereo_elements");
        for element in &mut component_elements {
            offset_object_u64(element, "index", element_count);
            offset_stereo_references(element, atom_count, bond_count, element_count, group_count);
        }
        let component_element_count = component_elements.len() as u64;
        stereo_elements.extend(component_elements);

        let mut component_groups = take_array(object, "stereo_groups");
        for group in &mut component_groups {
            offset_object_u64(group, "index", group_count);
            if let Some(members) = group
                .as_object_mut()
                .and_then(|object| object.get_mut("members"))
                .and_then(Value::as_array_mut)
            {
                for member in members {
                    if let Some(index) = member.as_u64() {
                        *member = json!(index + element_count);
                    }
                }
            }
        }
        let component_group_count = component_groups.len() as u64;
        stereo_groups.extend(component_groups);

        atom_count += component_atom_count;
        bond_count += component_bond_count;
        element_count += component_element_count;
        group_count += component_group_count;
    }

    candidates.sort_by_key(
        |candidate| match candidate.get("type").and_then(Value::as_str) {
            Some("tetrahedral") => (
                0u8,
                candidate
                    .get("center_atom_index")
                    .and_then(Value::as_u64)
                    .unwrap_or(u64::MAX),
            ),
            Some("double_bond") => (
                1u8,
                candidate
                    .get("center_bond_index")
                    .and_then(Value::as_u64)
                    .unwrap_or(u64::MAX),
            ),
            _ => (u8::MAX, u64::MAX),
        },
    );

    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": atom_count,
        "bond_count": bond_count,
        "report": {
            "is_ok": issues.is_empty(),
            "candidates": candidates,
            "issues": issues,
            "assembled_elements": assembled_elements,
            "created_element_indices": created_element_indices,
        },
        "stereo_elements": stereo_elements,
        "stereo_groups": stereo_groups,
    })
}

pub(crate) fn stereo_cip_record_json(
    record: &mut IndexedSmallRecord,
    remove_plain_hydrogens: bool,
) -> Option<Value> {
    if record.molecule.perceive().is_err() {
        return None;
    }
    if stereo::validate_stereo(&record.molecule).is_err() {
        return None;
    }
    stereo::assign_cip_descriptors(&mut record.molecule).ok()?;
    let mol = &record.molecule;
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
    atom_index: &BTreeMap<AtomId, u64>,
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
    atom_index: &BTreeMap<AtomId, u64>,
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

fn rdkit_default_atom_index(mol: &Molecule, remove_plain_hydrogens: bool) -> BTreeMap<AtomId, u64> {
    let mut index = BTreeMap::new();
    let retained = mol
        .atoms()
        .filter(|(_, atom)| !remove_plain_hydrogens || !rdkit_default_removes_hydrogen(atom));
    for (dense_index, (atom_id, _)) in (0u64..).zip(retained) {
        index.insert(atom_id, dense_index);
    }
    index
}

fn rdkit_default_bond_count(mol: &Molecule, atom_index: &BTreeMap<AtomId, u64>) -> usize {
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
    let _ = mol;
    Vec::new()
}

pub(crate) fn conformer_atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.hydrogens.explicit_count(),
        "atom_map": atom.atom_map,
        "aromatic": mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false),
    })
}

pub(crate) fn ring_membership_record_json(record: &mut IndexedSmallRecord) -> Value {
    let membership = rings::perceive_ring_membership(&mut record.molecule);
    let mol = &record.molecule;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_in_ring": mol.atom_ids().map(|id| membership.atom_in_ring(id)).collect::<Vec<_>>(),
        "bond_in_ring": mol.bond_ids().map(|id| membership.bond_in_ring(id)).collect::<Vec<_>>(),
    })
}

pub(crate) fn ring_set_record_json(record: &mut IndexedSmallRecord) -> Value {
    match rings::perceive_ring_set(&mut record.molecule) {
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

pub(crate) fn default_perception_atom_record_json(record: &mut IndexedSmallRecord) -> Value {
    if record.molecule.perceive().is_ok() {
        json!({
            "record_index": record.record_index,
            "status": "ok",
            "title": record.title,
            "atoms": basic_atoms_json(&record.molecule),
        })
    } else {
        json!({
            "record_index": record.record_index,
            "status": "perception_error",
            "title": record.title,
        })
    }
}

pub(crate) fn valence_record_json(record: &mut IndexedSmallRecord) -> Value {
    let result = valence::perceive_valence_with_options(
        &mut record.molecule,
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    );
    if result.is_err() {
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

            .atoms()
            .map(|(id, atom)| valence_atom_json(&record.molecule, id, atom))
            .collect::<Vec<_>>(),
    })
}

pub(crate) fn hydrogen_transform_record_json(record: &mut IndexedSmallRecord) -> Value {
    if record.molecule.perceive().is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "perception_error",
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

    if valence::perceive_valence(&mut record.molecule, ValenceModel::RdkitLike).is_err() {
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
        "round_trip": hydrogen_transform_semantic_json(record.molecule.clone()),
    })
}

pub(crate) fn aromaticity_record_json(record: &mut IndexedSmallRecord) -> Value {
    if record.molecule.perceive().is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "perception_error",
            "title": record.title,
        });
    }
    let mol = &record.molecule;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_aromatic": mol.atoms().map(|(id, _)| mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false)).collect::<Vec<_>>(),
        "bond_aromatic": mol.bonds().map(|(id, _)| mol.bond_is_aromatic(id).ok().flatten().unwrap_or(false)).collect::<Vec<_>>(),
    })
}

pub(crate) fn canonical_ranking_record_json(record: &mut IndexedSmallRecord) -> Value {
    if record.molecule.perceive().is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "perception_error",
            "title": record.title,
        });
    }
    let ranking = canon::atom_ranking(&record.molecule);
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

fn offset_object_u64(value: &mut Value, key: &str, offset: u64) {
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
        "explicit_hydrogens": atom.hydrogens.explicit_count(),
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
        "explicit_hydrogens": atom.hydrogens.explicit_count(),
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
        "explicit_hydrogens": atom.hydrogens.explicit_count(),
        "implicit_hydrogens": mol.implicit_hydrogens(id).ok().flatten().unwrap_or(0),
        "explicit_valence": explicit_valence_json(mol, id) + atom.hydrogens.explicit_count(),
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
    // The RDKit semantic record treats a pyrrolic donor H as explicit after
    // RDKit's prepared state. Kekule keeps an inferred H in Perception, so derive
    // the comparable bond-valence contribution without rewriting the atom.
    let has_aromatic_nitrogen_hydrogen = atom_record.is_some_and(|atom_record| {
        atom_record.element.symbol() == "N"
            && atom_record.formal_charge == 0
            && mol.atom_is_aromatic(atom).ok().flatten() == Some(true)
            && (atom_record.hydrogens.explicit_count() > 0
                || mol.implicit_hydrogens(atom).ok().flatten() == Some(1))
    });
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
                    has_aromatic_nitrogen_hydrogen,
                );
            }
            match bond.order {
                BondOrder::Zero | BondOrder::Dative => 0,
                BondOrder::Single => 2,
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
    has_aromatic_nitrogen_hydrogen: bool,
) -> u8 {
    let Some(atom) = atom else {
        return 2;
    };
    if atom_aromatic && has_non_aromatic_multiple_bond {
        return 2;
    }
    match atom.element.symbol() {
        "C" if atom.formal_charge < 0
            && (atom.hydrogens.explicit_count() > 0
                || has_non_aromatic_bond
                || aromatic_bond_count >= 3) =>
        {
            2
        }
        "P" | "As" | "Sb"
            if atom.formal_charge == 0
                && atom.hydrogens.explicit_count() == 0
                && (has_non_aromatic_bond || aromatic_bond_count >= 3) =>
        {
            2
        }
        "O" | "S" | "Se" | "Te"
            if atom.formal_charge == 0 && atom.hydrogens.explicit_count() == 0 =>
        {
            2
        }
        "N" if atom.formal_charge < 0 => 2,
        "N" if atom.formal_charge == 0 && has_aromatic_nitrogen_hydrogen => 2,
        "N" if atom.formal_charge == 0 && has_non_aromatic_bond => 2,
        "N" if atom.formal_charge == 0 && aromatic_bond_count >= 3 => 2,
        _ => 3,
    }
}

pub(crate) fn bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| bond_json(mol, id, bond))
        .collect::<Vec<_>>()
}

pub(crate) fn bond_json(mol: &Molecule, id: BondId, bond: &Bond) -> Value {
    json!({
        "index": id.raw(),
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": mol.bond_is_aromatic(id).ok().flatten().unwrap_or(false),
        "stereo": "STEREONONE",
        "bond_direction": "NONE",
    })
}

pub(crate) fn basic_bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| basic_bond_json(mol, id, bond))
        .collect::<Vec<_>>()
}

pub(crate) fn basic_bond_json(mol: &Molecule, id: BondId, bond: &Bond) -> Value {
    json!({
        "index": id.raw(),
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": mol.bond_is_aromatic(id).ok().flatten().unwrap_or(false),
        "stereo": "STEREONONE",
    })
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
        BondOrder::Dative => "DATIVE",
    }
}

pub(crate) fn stereo_perception_benchmark_report_json(
    mol: &Molecule,
    source_stereo_elements: &[StereoElementId],
    candidates: &[StereoCandidate],
    report: &CoordinateStereoMaterializationReport,
) -> Value {
    let assembled_elements = source_stereo_elements
        .iter()
        .filter_map(|id| {
            mol.stereo_element(*id)
                .ok()
                .map(|element| stereo_element_json(u64::from(id.raw()), element, None))
        })
        .collect::<Vec<_>>();
    json!({
        "is_ok": true,
        "candidates": candidates.iter().map(stereo_candidate_json).collect::<Vec<_>>(),
        "issues": [],
        "assembled_elements": assembled_elements,
        "created_element_indices": report
            .created_elements
            .iter()
            .map(|id| id.raw())
            .collect::<Vec<_>>(),
    })
}

fn take_array(object: &mut serde_json::Map<String, Value>, key: &str) -> Vec<Value> {
    match object.remove(key) {
        Some(Value::Array(values)) => values,
        _ => Vec::new(),
    }
}

fn offset_stereo_references(
    value: &mut Value,
    atom_offset: u64,
    bond_offset: u64,
    element_offset: u64,
    group_offset: u64,
) {
    match value {
        Value::Array(values) => {
            for value in values {
                offset_stereo_references(
                    value,
                    atom_offset,
                    bond_offset,
                    element_offset,
                    group_offset,
                );
            }
        }
        Value::Object(object) => {
            for (key, value) in object {
                let offset = if key == "atom_index" || key.ends_with("_atom_index") {
                    atom_offset
                } else if key == "bond_index" || key.ends_with("_bond_index") {
                    bond_offset
                } else if key == "element_index" || key.ends_with("_element_index") {
                    element_offset
                } else if key == "group_index" || key.ends_with("_group_index") {
                    group_offset
                } else {
                    offset_stereo_references(
                        value,
                        atom_offset,
                        bond_offset,
                        element_offset,
                        group_offset,
                    );
                    continue;
                };
                if let Some(index) = value.as_u64() {
                    *value = json!(index + offset);
                }
            }
        }
        _ => {}
    }
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

pub(crate) fn coordinate_stereo_error_json(error: &CoordinateStereoError) -> Value {
    let issues = match error {
        CoordinateStereoError::InvalidStereo(error) => error
            .issues
            .iter()
            .map(|issue| {
                json!({
                    "type": "invalid_stereo",
                    "issue": stereo_validation_issue_json(issue),
                })
            })
            .collect::<Vec<_>>(),
        CoordinateStereoError::CouldNotCreateElement(error) => vec![json!({
            "type": "could_not_create_element",
            "error": format!("{error:?}"),
        })],
        _ => vec![json!({
            "type": "coordinate_stereo_error",
            "error": format!("{error:?}"),
        })],
    };
    json!({ "issues": issues })
}

pub(crate) fn stereo_validation_issue_json(issue: &StereoValidationIssue) -> Value {
    match issue {
        StereoValidationIssue::MissingStereoAtom { element, atom } => json!({
            "type": "missing_stereo_atom",
            "element_index": element.raw(),
            "atom_index": atom.raw(),
        }),
        StereoValidationIssue::MissingStereoBond { element, bond } => json!({
            "type": "missing_stereo_bond",
            "element_index": element.raw(),
            "bond_index": bond.raw(),
        }),
        StereoValidationIssue::InvalidTetrahedralCarrierCount {
            element,
            center,
            carrier_count,
        } => json!({
            "type": "invalid_tetrahedral_carrier_count",
            "element_index": element.raw(),
            "center_atom_index": center.raw(),
            "carrier_count": carrier_count,
        }),
        StereoValidationIssue::DuplicateTetrahedralCarrier {
            element,
            center,
            carrier,
        } => json!({
            "type": "duplicate_tetrahedral_carrier",
            "element_index": element.raw(),
            "center_atom_index": center.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoValidationIssue::TetrahedralCarrierNotAdjacent {
            element,
            center,
            carrier,
        } => json!({
            "type": "tetrahedral_carrier_not_adjacent",
            "element_index": element.raw(),
            "center_atom_index": center.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoValidationIssue::InvalidDoubleBondOrder {
            element,
            bond,
            order,
        } => json!({
            "type": "invalid_double_bond_order",
            "element_index": element.raw(),
            "bond_index": bond.raw(),
            "bond_order": bond_order_json(*order),
        }),
        StereoValidationIssue::DoubleBondFocusMismatch {
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
        StereoValidationIssue::DoubleBondCarrierIsFocusAtom {
            element,
            endpoint,
            carrier,
        } => json!({
            "type": "double_bond_carrier_is_focus_atom",
            "element_index": element.raw(),
            "endpoint_atom_index": endpoint.raw(),
            "carrier_atom_index": carrier.raw(),
        }),
        StereoValidationIssue::DoubleBondCarrierNotAdjacent {
            element,
            endpoint,
            carrier,
        } => json!({
            "type": "double_bond_carrier_not_adjacent",
            "element_index": element.raw(),
            "endpoint_atom_index": endpoint.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoValidationIssue::UnsupportedDoubleBondCarrier {
            element,
            endpoint,
            carrier,
        } => json!({
            "type": "unsupported_double_bond_carrier",
            "element_index": element.raw(),
            "endpoint_atom_index": endpoint.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoValidationIssue::InvalidAxisCarrierCount {
            element,
            axis,
            carrier_count,
        } => json!({
            "type": "invalid_axis_carrier_count",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier_count": carrier_count,
        }),
        StereoValidationIssue::AxisCarrierIsFocusAtom {
            element,
            axis,
            carrier,
        } => json!({
            "type": "axis_carrier_is_focus_atom",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier_atom_index": carrier.raw(),
        }),
        StereoValidationIssue::AxisCarrierNotAdjacent {
            element,
            axis,
            carrier,
        } => json!({
            "type": "axis_carrier_not_adjacent",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoValidationIssue::UnsupportedAxisCarrier {
            element,
            axis,
            carrier,
        } => json!({
            "type": "unsupported_axis_carrier",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
    }
}

pub(crate) fn stereo_elements_json(mol: &Molecule) -> Vec<Value> {
    mol.stereo_elements()
        .map(|(id, element)| {
            stereo_element_json(
                u64::from(id.raw()),
                element,
                mol.cip_descriptor(id).ok().flatten(),
            )
        })
        .collect()
}

pub(crate) fn stereo_element_json(
    index: u64,
    element: &StereoElement,
    descriptor: Option<StereoDescriptor>,
) -> Value {
    let mut object = serde_json::Map::new();
    object.insert("index".to_owned(), json!(index));
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
                json!(stereo.orientation.map(tetrahedral_orientation_json)),
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
                json!(stereo.orientation.map(double_bond_orientation_json)),
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
                json!(stereo.orientation.map(axis_orientation_json)),
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

pub(crate) fn stereo_carrier_json(carrier: &StereoCarrier) -> Value {
    match carrier {
        StereoCarrier::Atom(atom) => json!({ "atom_index": atom.raw() }),
        StereoCarrier::ImplicitHydrogen => json!({ "implicit_hydrogen": true }),
        StereoCarrier::ImplicitLonePair => json!({ "implicit_lone_pair": true }),
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
