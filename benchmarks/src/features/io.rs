use crate::{boxed_error, dataset::Input};
use kekule::{core::Molecule, molfile, query, sdf, smiles, substructure};
use serde_json::{json, Value};
use std::error::Error;

#[derive(Debug, Clone)]
pub(crate) struct IndexedSmallRecord {
    pub(crate) record_index: usize,
    pub(crate) title: String,
    pub(crate) molecule: Molecule,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedStereoPerceptionRecord {
    pub(crate) record_index: usize,
    pub(crate) title: String,
    pub(crate) components: Vec<Molecule>,
    pub(crate) positions: Vec<Option<kekule::structure::Positions>>,
}

const BOUNDED_SUBSTRUCTURE_QUERIES: &str = include_str!("../../queries.smarts");

pub(super) fn smarts_query_records_json(path: &Input) -> Result<Vec<Value>, Box<dyn Error>> {
    let mut records = Vec::new();
    for (record_index, raw_line) in path.text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (smarts, title) = line.split_once(char::is_whitespace).unwrap_or((line, ""));
        let graph = query::parse_smarts(smarts)?;
        records.push(
            json!({"record_index":record_index,"status":"ok","smarts":smarts,
            "title":title.trim(),"atom_count":graph.atom_count(),"bond_count":graph.bond_count(),
            "behavior":super::query::observe(&graph)?}),
        );
    }
    Ok(records)
}

pub(super) fn substructure_record_json(
    record: &mut IndexedSmallRecord,
) -> Result<Value, Box<dyn Error>> {
    record.molecule.perceive()?;
    let mut queries = Vec::new();
    for smarts in BOUNDED_SUBSTRUCTURE_QUERIES.lines() {
        let graph =
            query::parse_smarts(smarts).expect("checked-in bounded benchmark SMARTS must parse");
        let matches = substructure::find_substructure_matches_with_options(
            &record.molecule,
            &graph,
            substructure::SubstructureMatchOptions {
                max_matches: usize::MAX,
                uniquify: false,
                ..Default::default()
            },
        )?;
        let mut atom_sets = matches
            .into_iter()
            .map(|query_match| {
                query_match
                    .atoms()
                    .iter()
                    .map(|atom| atom.raw())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        atom_sets.sort_unstable();
        queries.push(json!({
            "smarts": smarts,
            "matches": atom_sets,
        }));
    }
    Ok(json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "queries": queries,
    }))
}

pub(crate) fn read_small_records_by_suffix(
    path: &Input,
) -> Result<Vec<IndexedSmallRecord>, Box<dyn Error>> {
    let mut values = Vec::new();
    for record in read_stereo_perception_records_by_suffix(path)? {
        if record.components.is_empty() {
            return Err(boxed_error("molecular input did not parse"));
        }
        for molecule in record.components {
            values.push(IndexedSmallRecord {
                record_index: values.len(),
                title: record.title.clone(),
                molecule,
            });
        }
    }
    Ok(values)
}

fn interpret_smiles_components(input: &str) -> Result<Vec<Molecule>, Box<dyn Error>> {
    let document = smiles::parse_str(input)?;
    Ok(smiles::interpret(&document)?.into_molecules())
}

pub(crate) fn read_smiles_records(
    path: &Input,
) -> Result<Vec<IndexedStereoPerceptionRecord>, Box<dyn Error>> {
    let mut records = Vec::new();
    for (index, raw_line) in path.text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(2, char::is_whitespace);
        let smiles = parts.next().unwrap_or_default().to_owned();
        let title = parts.next().unwrap_or_default().trim().to_owned();
        if title.starts_with('|') {
            return Err(boxed_error(
                "CXSMILES extensions are not supported by Kekule; refusing to discard them",
            ));
        }
        let components = interpret_smiles_components(&smiles)?;
        records.push(IndexedStereoPerceptionRecord {
            record_index: index,
            title,
            positions: vec![None; components.len()],
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
        return read_smiles_records(path);
    }
    let input = &path.text;
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("mol" | "mdl")
    ) {
        return Ok(vec![stereo_molfile_record(0, &molfile::parse_str(input)?)?]);
    }
    sdf::parse_str(input)?
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
                    .map(|position| {
                        position
                            .value_in(kekule::units::ANGSTROM)
                            .expect("length units")
                    })
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
    })
}
