use super::io::IndexedSmallRecord;
use kekule::{core::BondOrder, perception::resonance::*};
use serde_json::{json, Value};
use std::{collections::BTreeMap, error::Error};

pub(super) fn record(
    record: &mut IndexedSmallRecord,
    feature: &str,
) -> Result<Value, Box<dyn Error>> {
    record.molecule.perceive()?;
    let mol = &record.molecule;
    let mut value = json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
    });
    let mut bonds: Vec<_> = mol
        .bonds()
        .map(|(id, b)| {
            let mut ends = [b.a().index(), b.b().index()];
            ends.sort();
            (ends, id)
        })
        .collect();
    bonds.sort_by_key(|&(ends, _)| ends);
    if feature == "algo.conjugation.rdkit-like" {
        let membership: Vec<_> = bonds
            .iter()
            .map(|&(ends, id)| {
                json!({
                    "atoms": ends,
                    "conjugated": mol.perception().bond_is_conjugated(id).unwrap(),
                })
            })
            .collect();
        value["bonds"] = json!(membership);
        return Ok(value);
    }
    perceive_resonance(&mut record.molecule)?;
    let mol = &record.molecule;
    let groups = mol.perception().resonance_state().unwrap().groups();
    let mut memberships: Vec<_> = groups
        .iter()
        .map(|g| {
            let atoms: Vec<_> = g.atoms.iter().map(|a| a.index()).collect();
            let group_bonds: Vec<_> = bonds
                .iter()
                .filter(|(_, id)| g.bonds.contains(id))
                .map(|&(ends, _)| ends)
                .collect();
            (atoms, group_bonds)
        })
        .collect();
    memberships.sort();
    value["groups"] = json!(memberships
        .into_iter()
        .map(|(atoms, bonds)| json!({"atoms":atoms,"bonds":bonds}))
        .collect::<Vec<_>>());
    if feature == "algo.resonance.groups" {
        return Ok(value);
    }
    value["bonds"] = json!(bonds.iter().map(|(ends, _)| ends).collect::<Vec<_>>());
    let mut profiles = Vec::new();
    for flags in 0..32 {
        let options = ResonanceOptions {
            flags: ResonanceFlags::from_bits(flags).unwrap(),
            ..Default::default()
        };
        let result = enumerate_resonance(mol, options)?;
        let mut structures = Vec::new();
        for c in result.contributors() {
            let charges: Vec<_> = c.formal_charges().map(|(_, q)| q).collect();
            let orders: BTreeMap<_, _> = c.bond_orders().collect();
            let orders: Vec<_> = bonds
                .iter()
                .map(|&(_, id)| match orders[&id] {
                    BondOrder::Zero => 0,
                    BondOrder::Single => 1,
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    BondOrder::Quadruple => 4,
                    BondOrder::Dative => 5,
                })
                .collect();
            structures.push((charges, orders));
        }
        structures.sort();
        let structures: Vec<_> = structures
            .into_iter()
            .map(|(charges, orders)| json!({"charges": charges, "orders": orders}))
            .collect();
        profiles.push(json!({
            "flags": flags,
            "max_structures": 1000,
            "count": structures.len(),
            "structures": structures,
        }));
    }
    value["profiles"] = json!(profiles);
    Ok(value)
}
