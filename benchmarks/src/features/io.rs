use crate::*;

use super::chemistry::{atom_json, atoms_json, bond_json, bonds_json};

#[derive(Debug, Clone)]
pub(crate) struct IndexedSmallRecord {
    pub(crate) record_index: usize,
    pub(crate) title: String,
    pub(crate) molecule: Molecule,
    pub(crate) sdf_fields: BTreeMap<String, String>,
}

pub(super) fn round_trip_molfiles(
    records: Vec<IndexedSmallRecord>,
    write: fn(&Molecule) -> Result<String, molfile::MolWriteError>,
) -> Result<Vec<IndexedSmallRecord>, Box<dyn Error>> {
    records
        .into_iter()
        .enumerate()
        .map(|(index, record)| {
            let written = write(&record.molecule)?;
            let molecule = interpret_molfile(&written)?;
            Ok(IndexedSmallRecord {
                record_index: index,
                title: record.title,
                molecule,
                sdf_fields: BTreeMap::new(),
            })
        })
        .collect()
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedSmilesRecord {
    pub(crate) record_index: usize,
    pub(crate) status: String,
    pub(crate) title: String,
    pub(crate) input_smiles: String,
    pub(crate) components: Vec<Molecule>,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedStereoPerceptionRecord {
    pub(crate) record_index: usize,
    pub(crate) title: String,
    pub(crate) components: Vec<Molecule>,
    pub(crate) positions: Vec<Option<kekule::structure::Positions>>,
    pub(crate) document: Value,
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

pub(super) fn smarts_query_records_json(path: &Input) -> Result<Vec<Value>, Box<dyn Error>> {
    let mut records = Vec::new();
    for (record_index, raw_line) in path.text.clone().lines().enumerate() {
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

pub(super) fn substructure_record_json(record: &mut IndexedSmallRecord) -> Value {
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
    path: &Input,
) -> Result<Vec<IndexedSmallRecord>, Box<dyn Error>> {
    let input = path.text.clone();
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("mol" | "mdl")
    ) {
        let document = molfile::parse_str(&input)?;
        let title = document.header().title().to_owned();
        let molecule = exactly_one_molecule(molfile::interpret(&document)?.into_molecules())?;
        return Ok(vec![IndexedSmallRecord {
            record_index: 0,
            title,
            molecule,
            sdf_fields: BTreeMap::new(),
        }]);
    }
    interpret_sdf(&input)?
        .into_iter()
        .enumerate()
        .map(|(index, record)| small_record(index, record))
        .collect()
}

pub(crate) fn small_record(
    index: usize,
    record: SdfRecordInterpretation,
) -> Result<IndexedSmallRecord, Box<dyn Error>> {
    let title = record.title().to_owned();
    let sdf_fields = record
        .data_fields()
        .iter()
        .map(|field| (field.name().to_owned(), field.value().to_owned()))
        .collect();
    Ok(IndexedSmallRecord {
        record_index: index,
        title,
        molecule: exactly_one_molecule(record.into_molecules())?,
        sdf_fields,
    })
}

pub(super) fn sdf_parse_records_json(path: &Input) -> Result<Vec<Value>, Box<dyn Error>> {
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("mol" | "mdl")
    ) {
        return Ok(read_small_records_by_suffix(path)?
            .iter()
            .map(sdf_record_json)
            .collect());
    }
    interpret_sdf(&path.text)?
        .iter()
        .enumerate()
        .map(|(record_index, record)| {
            let molecules = record.molecules().collect::<Vec<_>>();
            let reports = record.report().molfile_components();
            if molecules.len() != reports.len() {
                return Err(boxed_error("missing SDF component correspondence"));
            }

            // Components can be interleaved in the source. Compare in source order,
            // keeping local molecule IDs separate from record-wide reference IDs.
            let mut atoms = Vec::new();
            for (component, (molecule, report)) in molecules.iter().zip(reports).enumerate() {
                for mapping in report.atom_mappings() {
                    let id = mapping.atom();
                    atoms.push((
                        mapping.source_line(),
                        component,
                        id,
                        atom_json(molecule, id, molecule.atom(id)?),
                    ));
                }
            }
            atoms.sort_by_key(|(line, ..)| *line);
            let mut indices = BTreeMap::new();
            for (index, (_, component, id, atom)) in atoms.iter_mut().enumerate() {
                indices.insert((*component, *id), index);
                atom["index"] = json!(index);
            }

            let mut bonds = Vec::new();
            for (component, (molecule, report)) in molecules.iter().zip(reports).enumerate() {
                for mapping in report.bond_mappings() {
                    let id = mapping.bond();
                    let bond = molecule.bond(id)?;
                    let mut value = bond_json(molecule, id, bond);
                    for (field, atom) in
                        [("begin_atom_index", bond.a()), ("end_atom_index", bond.b())]
                    {
                        value[field] = json!(indices
                            .get(&(component, atom))
                            .ok_or_else(|| boxed_error("missing SDF atom correspondence"))?);
                    }
                    bonds.push((mapping.source_line(), value));
                }
            }
            bonds.sort_by_key(|(line, _)| *line);
            for (index, (_, bond)) in bonds.iter_mut().enumerate() {
                bond["index"] = json!(index);
            }
            let properties = record
                .data_fields()
                .iter()
                .map(|field| (field.name(), field.value()))
                .collect::<BTreeMap<_, _>>();
            Ok(json!({
                "record_index": record_index, "status": "ok", "title": record.title(),
                "atom_count": atoms.len(), "bond_count": bonds.len(),
                "atoms": atoms.into_iter().map(|(_, _, _, atom)| atom).collect::<Vec<_>>(),
                "bonds": bonds.into_iter().map(|(_, bond)| bond).collect::<Vec<_>>(),
                "properties": properties,
            }))
        })
        .collect()
}

pub(super) fn interpret_molfile(input: &str) -> Result<Molecule, Box<dyn Error>> {
    let document = molfile::parse_str(input)?;
    exactly_one_molecule(molfile::interpret(&document)?.into_molecules())
}

pub(super) fn interpret_sdf(input: &str) -> Result<Vec<SdfRecordInterpretation>, Box<dyn Error>> {
    let document = sdf::parse_str(input)?;
    Ok(sdf::interpret(&document)?.into_records())
}

pub(super) fn zero_coordinate_model(
    molecule: &Molecule,
) -> Result<kekule::structure::Model, Box<dyn Error>> {
    let positions = kekule::structure::Positions::new(kekule::units::Quantity::new(
        vec![kekule::geometry::Point3::default(); molecule.atom_count()],
        kekule::units::ANGSTROM,
    ))?;
    Ok(kekule::structure::Model::from_molecule(
        molecule, &positions,
    )?)
}

pub(super) fn interpret_smiles(input: &str) -> Result<Molecule, Box<dyn Error>> {
    let document = smiles::parse_str(input)?;
    exactly_one_molecule(smiles::interpret(&document)?.into_molecules())
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
    Ok(smiles::interpret(&document)?.into_molecules())
}

pub(crate) fn read_smiles_records(
    path: &Input,
) -> Result<Vec<IndexedSmilesRecord>, Box<dyn Error>> {
    let mut records = Vec::new();
    for (index, raw_line) in path.text.clone().lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(2, char::is_whitespace);
        let smiles = parts.next().unwrap_or_default().to_owned();
        let title = parts.next().unwrap_or_default().trim().to_owned();
        let (status, components) = match interpret_smiles_components(&smiles) {
            Ok(components) => ("ok".to_owned(), components),
            Err(_) => ("parse_error".to_owned(), Vec::new()),
        };
        records.push(IndexedSmilesRecord {
            record_index: index,
            status,
            title,
            input_smiles: smiles,
            components,
        });
    }
    Ok(records)
}

pub(crate) fn read_stereo_perception_records_by_suffix(
    path: &Input,
) -> Result<Vec<IndexedStereoPerceptionRecord>, Box<dyn Error>> {
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("txt" | "smi" | "smiles")
    ) {
        return read_smiles_records(path)?
            .into_iter()
            .map(|record| {
                let mut marks = Vec::new();
                if let Ok(document) = smiles::parse_str(&record.input_smiles) {
                    if let Ok(interpretation) = smiles::interpret(&document) {
                        let mut bond_offset = 0usize;
                        for component in interpretation.components() {
                            for mapping in component.report().bond_mappings() {
                                let kind = match document.source().as_bytes().get(mapping.source_offset()) {
                                    Some(b'/') => Some("directional_up"),
                                    Some(b'\\') => Some("directional_down"),
                                    _ => None,
                                };
                                if let Some(kind) = kind {
                                    marks.push(json!({"bond_index": bond_offset + mapping.bond().index(), "kind": kind, "source": "smiles"}));
                                }
                            }
                            bond_offset += component.molecule().bond_count();
                        }
                    }
                }
                Ok(IndexedStereoPerceptionRecord {
                    record_index: record.record_index,
                    title: record.title,
                    positions: vec![None; record.components.len()],
                    components: record.components,
                    document: json!({"format": "smiles", "source": record.input_smiles, "stereo_bond_marks": marks}),
                })
            })
            .collect();
    }
    let input = path.text.clone();
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("mol" | "mdl")
    ) {
        return Ok(vec![stereo_molfile_record(
            0,
            &molfile::parse_str(&input)?,
        )?]);
    }
    sdf::parse_str(&input)?
        .records()
        .iter()
        .enumerate()
        .map(|(index, record)| stereo_molfile_record(index, record.molfile()))
        .collect()
}

fn stereo_molfile_record(
    record_index: usize,
    document: &molfile::MolfileDocument,
) -> Result<IndexedStereoPerceptionRecord, Box<dyn Error>> {
    let format = match document.version() {
        molfile::MolfileVersion::V2000 => "molfile_v2000",
        molfile::MolfileVersion::V3000 => "molfile_v3000",
    };
    let mut marks = Vec::new();
    for (index, line) in document.bond_records().iter().enumerate() {
        let code = if format == "molfile_v2000" {
            line.text().get(9..12).unwrap_or("").trim()
        } else {
            line.text()
                .split_whitespace()
                .find_map(|part| part.strip_prefix("CFG="))
                .unwrap_or("0")
        };
        let kind = match (format, code) {
            (_, "0" | "") => None,
            ("molfile_v2000", "1") | ("molfile_v3000", "1") => Some("wedge_up"),
            ("molfile_v2000", "6") | ("molfile_v3000", "3") => Some("wedge_down"),
            ("molfile_v2000", "4") | ("molfile_v3000", "2") => Some("unknown"),
            ("molfile_v2000", "3") => Some("either_double"),
            _ => {
                return Err(boxed_error(format!(
                    "unrecognized stereo source code {code}"
                )))
            }
        };
        if let Some(kind) = kind {
            marks.push(json!({"bond_index": index, "kind": kind, "source": format}));
        }
    }
    let interpretation = molfile::interpret(document)?;
    let model = interpretation.model();
    let mut positions = Vec::new();
    for instance in model.topology().molecules() {
        let points = instance
            .molecule()
            .atom_ids()
            .map(|atom| {
                model
                    .position(kekule::topology::InstanceAtomId::new(instance.id(), atom))
                    .map(|position| position.into_value())
                    .map_err(Box::new)
            })
            .collect::<Result<Vec<_>, _>>()?;
        positions.push(Some(kekule::structure::Positions::new(
            kekule::units::Quantity::new(points, kekule::units::ANGSTROM),
        )?));
    }
    Ok(IndexedStereoPerceptionRecord {
        record_index,
        title: document.header().title().to_owned(),
        positions,
        components: interpretation.into_molecules(),
        document: json!({"format": format, "source": document.source(), "stereo_bond_marks": marks}),
    })
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
