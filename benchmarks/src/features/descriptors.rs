use crate::*;

use super::chemistry::*;
use super::io::IndexedSmallRecord;

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

pub(crate) fn ring_membership_record_json(record: &mut IndexedSmallRecord) -> Value {
    let membership = rings::perceive_ring_membership(&mut record.molecule);
    let mol = &record.molecule;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_in_ring": mol.atom_ids().map(|id| membership.atom_in_ring(id)).collect::<Vec<_>>(),
        "bond_in_ring": bond_values(mol, |id| json!(membership.bond_in_ring(id))),
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
            "graph": super::strict::graph(&record.molecule, None).unwrap_or_else(|e| json!({"status":"error","message":e.to_string()})),
            "valence": record.molecule.atoms().map(|(id, atom)| valence_atom_json(&record.molecule, id, atom)).collect::<Vec<_>>(),
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
    if record.molecule.perceive().is_err() {
        return json!({"status":"perception_error"});
    }
    let added_graph = super::strict::graph(&record.molecule, None)
        .unwrap_or_else(|e| json!({"status":"error","message":e.to_string()}));
    if hydrogens::remove_hydrogens(&mut record.molecule).is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "remove_error",
            "title": record.title,
        });
    }

    if record.molecule.perceive().is_err() {
        return json!({"status":"perception_error"});
    }
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count_after_add": atom_count_after_add,
        "added_graph": added_graph,
        "added_hydrogens_by_parent": added_by_parent
            .into_iter()
            .map(|(parent_atom_index, count)| json!({
                "parent_atom_index": parent_atom_index,
                "count": count,
            }))
            .collect::<Vec<_>>(),
        "round_trip": super::strict::graph(&record.molecule, None).unwrap_or_else(|e| json!({"status":"error","message":e.to_string()})),
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
        "bond_aromatic": bond_values(mol, |id| json!(mol.bond_is_aromatic(id).expect("live bond"))),
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

fn bond_values(mol: &Molecule, value: impl Fn(kekule::core::BondId) -> Value) -> Vec<Value> {
    let mut bonds = mol
        .bonds()
        .map(|(id, bond)| {
            let mut ends = [bond.a().raw(), bond.b().raw()];
            if bond.order != BondOrder::Dative {
                ends.sort();
            }
            (
                ends,
                json!({"begin_atom_index":ends[0],"end_atom_index":ends[1],"value":value(id)}),
            )
        })
        .collect::<Vec<_>>();
    bonds.sort_by_key(|(ends, _)| *ends);
    bonds.into_iter().map(|(_, value)| value).collect()
}
