//! Indexed JSONL probe for focused differential resonance regressions.
//! Corpus benchmarks use the separately registered, provenance-checked features.
use kekule::{core::BondOrder, perception::resonance::*, smiles};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    error::Error,
    io::{self, BufRead},
};

fn evaluate(v: &Value) -> Result<Value, Box<dyn Error>> {
    let mut molecules = if let Some(path) = v["molfile"].as_str() {
        kekule::molfile::interpret(&kekule::molfile::parse_str(&std::fs::read_to_string(
            path,
        )?)?)?
        .into_molecules()
    } else {
        smiles::to_molecules(v["smiles"].as_str().ok_or("missing smiles")?)?
    };
    if molecules.len() != 1 {
        return Err("probe requires one connected molecule".into());
    }
    let mol = &mut molecules[0];
    mol.perceive()?;
    let flags = ResonanceFlags::from_bits(v["flags"].as_u64().unwrap_or(0).try_into()?)
        .ok_or("invalid flags")?;
    let options = ResonanceOptions {
        flags,
        max_structures: v["max_structures"].as_u64().unwrap_or(1000).try_into()?,
        ..Default::default()
    };
    let result = enumerate_resonance(mol, options)?;
    let mut structures = Vec::new();
    for c in result.contributors() {
        let charges: Vec<_> = c.formal_charges().map(|(_, q)| q).collect();
        let orders: BTreeMap<_, _> = c.bond_orders().collect();
        let mut bonds: Vec<_> = mol
            .bonds()
            .map(|(id, b)| {
                let mut ends = [b.a().index(), b.b().index()];
                ends.sort();
                (
                    ends,
                    match orders[&id] {
                        BondOrder::Zero => 0,
                        BondOrder::Single => 1,
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                        BondOrder::Quadruple => 4,
                        BondOrder::Dative => 5,
                    },
                )
            })
            .collect();
        bonds.sort();
        structures.push((charges, bonds));
    }
    structures.sort();
    Ok(json!({"structures":structures}))
}
fn main() -> Result<(), Box<dyn Error>> {
    for line in io::stdin().lock().lines() {
        let v: Value = serde_json::from_str(&line?)?;
        let out = match evaluate(&v) {
            Ok(v) => v,
            Err(e) => json!({"error":e.to_string()}),
        };
        println!("{out}");
    }
    Ok(())
}
