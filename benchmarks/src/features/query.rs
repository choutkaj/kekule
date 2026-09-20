//! General SMARTS behavioral observations against provenance-pinned molecules.
use crate::boxed_error;
use kekule::{
    query, smiles,
    substructure::{self, PreparedTopologyTarget},
    topology::{InstanceAtomId, Topology},
};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    error::Error,
    sync::{Arc, OnceLock},
};

pub(crate) const CONTRACT: &str = include_str!("../../query-smarts.json");
struct Target {
    id: String,
    topology: Arc<Topology>,
    indices: BTreeMap<InstanceAtomId, usize>,
}
static TARGETS: OnceLock<Result<Vec<Target>, String>> = OnceLock::new();
static PREPARED: OnceLock<Vec<PreparedTopologyTarget<'static>>> = OnceLock::new();

fn targets() -> Result<&'static Vec<Target>, Box<dyn Error>> {
    TARGETS
        .get_or_init(|| {
            let build = || -> Result<Vec<Target>, Box<dyn Error>> {
                let contract: Value = serde_json::from_str(CONTRACT)?;
                contract["targets"]
                    .as_array()
                    .ok_or("missing target panel")?
                    .iter()
                    .map(|row| {
                        let mut molecules =
                            smiles::to_molecules(row["smiles"].as_str().ok_or("target SMILES")?)?;
                        for molecule in &mut molecules {
                            molecule.perceive()?;
                        }
                        let topology = Arc::new(Topology::from_molecules(&molecules)?);
                        let indices = topology
                            .molecules()
                            .flat_map(|m| m.atoms().map(|(id, _)| id))
                            .enumerate()
                            .map(|(i, id)| (id, i))
                            .collect();
                        Ok(Target {
                            id: row["id"].as_str().ok_or("target ID")?.into(),
                            topology,
                            indices,
                        })
                    })
                    .collect()
            };
            build().map_err(|e| e.to_string())
        })
        .as_ref()
        .map_err(|e| boxed_error(format!("SMARTS target preparation failed: {e}")))
}

pub(super) fn observe(graph: &query::QueryGraph) -> Result<Value, Box<dyn Error>> {
    let targets = targets()?;
    let prepared = PREPARED.get_or_init(|| {
        targets
            .iter()
            .map(|t| PreparedTopologyTarget::new(&t.topology))
            .collect()
    });
    let contract: Value = serde_json::from_str(CONTRACT)?;
    let options = substructure::SubstructureMatchOptions {
        max_matches: contract["max_matches"].as_u64().ok_or("match limit")? as usize,
        max_search_states: contract["max_search_states"]
            .as_u64()
            .ok_or("search limit")? as usize,
        max_candidate_pairs: contract["max_candidate_pairs"]
            .as_u64()
            .ok_or("candidate limit")? as usize,
        uniquify: false,
        ..Default::default()
    };
    let mut tags = graph.tagged_atoms().collect::<Vec<_>>();
    tags.sort();
    let mut observations = Vec::new();
    for (target, prepared) in targets.iter().zip(prepared) {
        let matches = prepared
            .find_matches_complete(graph, options)
            .map_err(|e| boxed_error(format!("SMARTS target {}: {e}", target.id)))?;
        let mut mappings = matches
            .iter()
            .map(|m| {
                m.atoms()
                    .iter()
                    .map(|a| target.indices[a])
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        mappings.sort();
        let tagged = mappings
            .iter()
            .map(|m| tags.iter().map(|(_, a)| m[a.index()]).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        observations.push(json!({"target":target.id,"matches":mappings,"tagged_matches":tagged}));
    }
    Ok(
        json!({"version":2,"tags":tags.iter().map(|(tag,atom)|json!({"tag":tag,"query_atom":atom.index()})).collect::<Vec<_>>(),"targets":observations}),
    )
}
