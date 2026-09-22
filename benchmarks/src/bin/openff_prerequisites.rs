//! Read-only graph/feature observer for the optional OpenFF prerequisite audit.
//! Inputs are fully explicit, uniquely mapped molecules; outputs use map labels,
//! never an assumed agreement between SMILES traversal and reference atom order.
use kekule::{
    core::{AromaticityModel, RingBasisModel},
    perception::aromaticity,
    smiles,
};
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::io::{self, BufRead};

fn observe(source: &str) -> Result<Value, Box<dyn std::error::Error>> {
    let mut molecules = smiles::to_molecules(source)?;
    if molecules.len() != 1 {
        return Err("expected one connected molecule".into());
    }
    let mut molecule = molecules.remove(0);
    molecule.perceive()?;
    let mut labels = BTreeSet::new();
    for (id, atom) in molecule.atoms() {
        let label = atom.atom_map.ok_or("missing atom map")?;
        if label == 0 || !labels.insert(label) {
            return Err("atom maps must be positive and unique".into());
        }
        if molecule.implicit_hydrogens(id)? != Some(0) {
            return Err("all hydrogens must be explicit".into());
        }
    }
    if molecule.perception().ring_basis_model() != Some(RingBasisModel::FiguerasSssrLike) {
        return Err("unsupported ring basis".into());
    }
    aromaticity::perceive_aromaticity(&mut molecule, AromaticityModel::Mdl)?;
    let mut atoms = Vec::new();
    for (id, atom) in molecule.atoms() {
        let ring_sizes: Vec<_> = (3..=6)
            .map(|size| {
                molecule
                    .ring_set()
                    .unwrap()
                    .rings()
                    .iter()
                    .any(|ring| ring.atoms.len() == size && ring.atoms.contains(&id))
            })
            .collect();
        atoms.push(json!({
            "map": atom.atom_map.unwrap(), "element": atom.element.atomic_number(),
            "formal_charge": atom.formal_charge, "degree": molecule.neighbors(id)?.count(),
            "rings_3_6": ring_sizes,
            "mdl_aromatic": molecule.perception().atom_is_aromatic(id).unwrap()
        }));
    }
    atoms.sort_by_key(|a| a["map"].as_u64().unwrap());
    Ok(json!({"status":"ok", "atoms":atoms}))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let output = match serde_json::from_str::<Value>(&line?) {
            Ok(row) => match row["smiles"].as_str() {
                Some(s) => observe(s)
                    .unwrap_or_else(|e| json!({"status":"error", "message":e.to_string()})),
                None => json!({"status":"error", "message":"missing smiles"}),
            },
            Err(e) => json!({"status":"error", "message":e.to_string()}),
        };
        println!("{output}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maps_bind_observations_across_traversal_order() {
        assert_eq!(
            observe("[H:2][O:1][H:3]").unwrap(),
            observe("[H:3][O:1][H:2]").unwrap()
        );
        assert_eq!(observe("[H:2][O:1][H:3]").unwrap()["atoms"][0]["degree"], 2);
    }

    #[test]
    fn incomplete_or_ambiguous_identity_is_rejected() {
        for source in ["O", "[OH2:1]", "[H:1][O:1][H:3]", "[Na+:1].[Cl-:2]"] {
            assert!(observe(source).is_err(), "{source}");
        }
    }
}
