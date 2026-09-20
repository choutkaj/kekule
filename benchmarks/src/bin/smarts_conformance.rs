//! JSONL protocol for optional external SMARTS conformance comparisons.
use kekule::{
    core::AromaticityModel, hydrogens, perception::aromaticity, query, smiles, substructure::*,
    topology::Topology,
};
use serde_json::{json, Value};
use std::io::{self, BufRead};
use std::sync::Arc;

fn evaluate(row: &Value) -> Value {
    let Some(pattern) = row["smarts"].as_str() else {
        return json!({"status":"invalid_request"});
    };
    let query = match query::parse_smarts(pattern) {
        Ok(q) => q,
        Err(e) => {
            return json!({"status":"parse_error","kind":format!("{:?}",e.kind()),"span":[e.span().start,e.span().end],"message":e.message()})
        }
    };
    if row["smiles"].is_null() {
        return json!({"status":"ok","atom_count":query.atom_count(),"bond_count":query.bond_count()});
    }
    let source = row["smiles"].as_str().unwrap();
    let build = || -> Result<Arc<Topology>, Box<dyn std::error::Error>> {
        let mut molecules = smiles::to_molecules(source)?;
        for m in &mut molecules {
            m.perceive()?;
            if row["explicit_hydrogens"] == true {
                hydrogens::add_hydrogens(m)?;
                m.perceive()?;
            }
            if row["mdl"] == true {
                aromaticity::perceive_aromaticity(m, AromaticityModel::Mdl)?;
            }
        }
        Ok(Arc::new(Topology::from_molecules(&molecules)?))
    };
    let topology = match build() {
        Ok(t) => t,
        Err(e) => return json!({"status":"target_error","message":e.to_string()}),
    };
    let indices = topology
        .molecules()
        .flat_map(|m| m.atoms().map(|(id, _)| id))
        .enumerate()
        .map(|(i, id)| (id, i))
        .collect::<std::collections::BTreeMap<_, _>>();
    let options = SubstructureMatchOptions {
        max_matches: 1_000_000,
        max_search_states: 10_000_000,
        max_candidate_pairs: 10_000_000,
        uniquify: false,
        ..Default::default()
    };
    let matches = match find_topology_substructure_matches_complete(&topology, &query, options) {
        Ok(m) => m,
        Err(e) => {
            return json!({"status":"match_error","kind":format!("{e:?}"),"message":e.to_string()})
        }
    };
    let mut mappings = matches
        .iter()
        .map(|m| m.atoms().iter().map(|a| indices[a]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    mappings.sort();
    let mut tags = query.tagged_atoms().collect::<Vec<_>>();
    tags.sort();
    let tagged = mappings
        .iter()
        .map(|m| tags.iter().map(|(_, a)| m[a.index()]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let aromatic_atoms = topology
        .molecules()
        .flat_map(|m| {
            m.molecule()
                .atom_ids()
                .map(move |a| m.molecule().perception().atom_is_aromatic(a).unwrap())
        })
        .collect::<Vec<_>>();
    json!({"status":"ok","atom_count":query.atom_count(),"bond_count":query.bond_count(),"matches":mappings,"tagged_matches":tagged,"tags":tags.iter().map(|(t,_)|*t).collect::<Vec<_>>(),"target_aromatic_atoms":aromatic_atoms})
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let row: Value = serde_json::from_str(&line?)?;
        println!("{}", evaluate(&row));
    }
    Ok(())
}
