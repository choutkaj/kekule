//! Small-molecule tasks. Each returns facts keyed by source atom indices.
//!
//! Every task starts from the same preparation: default perception, removal of
//! non-stereogenic assertions, and stereo materialized from coordinates when
//! the source has them, so downstream tasks see the chemistry `parse` reports.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use kekule::core::{
    AromaticityModel, AxisOrientation, BondOrder, DoubleBondOrientation, StereoCarrier,
    StereoDescriptor, StereoElementId, StereoElementKind, StereoGroupKind, TetrahedralOrientation,
};
use kekule::descriptors::{molecular_formula, HydrogenCountPolicy};
use kekule::perception::{aromaticity, resonance, rings};
use kekule::rotatable_bonds::{self, RotatableBondOptions};
use kekule::stereo::{self, StereoCandidate};
use kekule::substructure::{PreparedTopologyTarget, SubstructureMatchOptions};
use kekule::topology::Topology;
use kekule::{canon, molfile, query, sdf, smiles};
use serde_json::{json, Value};

use crate::input::{self, Component};
use crate::{Facts, Failure};

/// Matches per query and target, the same bound the reference observer uses.
const MAX_MATCHES: usize = 100_000;

pub fn observe(task: &str, format: &str, text: &str, options: &Value) -> Result<Facts, Failure> {
    if task == "write" {
        return write(format, text, options);
    }
    let mut components = input::read(format, text)?;
    for component in &mut components {
        prepare(component)?;
    }
    let mut facts = Facts::new();
    match task {
        "parse" => {
            parse(&components, &mut facts)?;
            if options["coordinates"].as_bool() == Some(true) {
                coordinates(&components, &mut facts)?;
            }
        }
        "rings" => ring_facts(&mut components, &mut facts)?,
        "aromaticity" => aromaticity_facts(&mut components, options, &mut facts)?,
        "conjugation" => conjugation(&components, &mut facts)?,
        "resonance" => resonance_facts(&mut components, &mut facts)?,
        "stereo.candidates" => candidates(&components, &mut facts)?,
        "stereo.cip" => cip(&mut components, &mut facts)?,
        "symmetry" => symmetry(&components, &mut facts),
        "formula" => formula(&components, &mut facts)?,
        "rotatable" => rotatable(&components, &mut facts)?,
        "smarts" => smarts(components, options, &mut facts)?,
        other => {
            return Err(Failure::new(
                "request",
                format!("unknown task small.{other}"),
            ))
        }
    }
    Ok(facts)
}

fn prepare(component: &mut Component) -> Result<(), Failure> {
    let molecule = &mut component.molecule;
    molecule
        .perceive()
        .map_err(|e| Failure::new("perception", e))?;
    let mut editor = molecule.edit();
    stereo::cleanup_stereo(&mut editor, Default::default())
        .map_err(|e| Failure::new("stereo", e))?;
    if let Some(positions) = &component.positions {
        stereo::materialize_coordinate_stereo(&mut editor, positions)
            .map_err(|e| Failure::new("stereo", e))?;
        stereo::cleanup_stereo(&mut editor, Default::default())
            .map_err(|e| Failure::new("stereo", e))?;
    }
    *molecule = editor.finish().map_err(|e| Failure::new("stereo", e))?;
    molecule
        .perceive()
        .map_err(|e| Failure::new("perception", e))
}

/// Bond endpoints as source indices; dative bonds keep their direction.
fn ends(component: &Component, bond: &kekule::core::Bond) -> [usize; 2] {
    let mut ends = [component.atom(bond.a()), component.atom(bond.b())];
    if bond.order != BondOrder::Dative {
        ends.sort_unstable();
    }
    ends
}

fn order(order: BondOrder) -> Value {
    match order {
        BondOrder::Zero => json!(0),
        BondOrder::Single => json!(1),
        BondOrder::Double => json!(2),
        BondOrder::Triple => json!(3),
        BondOrder::Quadruple => json!(4),
        BondOrder::Dative => json!("dative"),
    }
}

fn min_atom(component: &Component) -> usize {
    component.source.iter().copied().min().unwrap_or(0)
}

fn parse(components: &[Component], facts: &mut Facts) -> Result<(), Failure> {
    for component in components {
        let molecule = &component.molecule;
        for (id, atom) in molecule.atoms() {
            let i = component.atom(id);
            let hydrogens = molecule
                .total_hydrogens(id)
                .map_err(|e| Failure::new("perception", e))?;
            let radicals = atom.radical.map(|r| r.electron_count()).unwrap_or(0);
            facts.push(json!(["atom", i, "element", atom.element.symbol()]));
            facts.push(json!(["atom", i, "isotope", atom.isotope]));
            facts.push(json!(["atom", i, "charge", atom.formal_charge]));
            facts.push(json!(["atom", i, "radicals", radicals]));
            facts.push(json!(["atom", i, "hydrogens", hydrogens]));
            facts.push(json!(["atom", i, "map", atom.atom_map]));
        }
        for (_, bond) in molecule.bonds() {
            let [a, b] = ends(component, bond);
            facts.push(json!(["bond", a, b, "order", order(bond.order)]));
        }
        let mut atoms = component.source.clone();
        atoms.sort_unstable();
        facts.push(json!(["component", min_atom(component), "atoms", atoms]));
        stereo_facts(component, facts)?;
    }
    Ok(())
}

/// Source coordinates in Å, for comparing written files with their source.
fn coordinates(components: &[Component], facts: &mut Facts) -> Result<(), Failure> {
    for component in components {
        let Some(positions) = &component.positions else {
            continue;
        };
        for id in component.molecule.atom_ids() {
            let point = positions
                .position_at(id.index())
                .map_err(|e| Failure::new("parse", e))?
                .value_in(kekule::units::ANGSTROM)
                .map_err(|e| Failure::new("parse", e))?;
            facts.push(json!([
                "atom",
                component.atom(id),
                "xyz",
                [point.x, point.y, point.z]
            ]));
        }
    }
    Ok(())
}

fn permutation_parity(values: &[i64]) -> u8 {
    let mut parity = 0;
    for i in 0..values.len() {
        for j in i + 1..values.len() {
            parity ^= u8::from(values[i] > values[j]);
        }
    }
    parity
}

/// Lowest-numbered neighbour of `center` other than `other`, as in RDKit's
/// double-bond and atropisomer reference carriers.
fn reference_carrier(
    component: &Component,
    center: kekule::core::AtomId,
    other: kekule::core::AtomId,
) -> Result<i64, Failure> {
    Ok(component
        .molecule
        .neighbors(center)
        .map_err(|e| Failure::new("stereo", e))?
        .filter(|id| *id != other)
        .map(|id| component.atom(id) as i64)
        .min()
        .unwrap_or(-1))
}

/// Stereo as canonical parity over source indices. Carriers -1 and -2 are an
/// implicit hydrogen and an implicit lone pair.
fn stereo_facts(component: &Component, facts: &mut Facts) -> Result<(), Failure> {
    let molecule = &component.molecule;
    let carrier = |c: &StereoCarrier| match c {
        StereoCarrier::Atom(id) => component.atom(*id) as i64,
        StereoCarrier::ImplicitHydrogen => -1,
        StereoCarrier::ImplicitLonePair => -2,
    };
    let mut members: BTreeMap<StereoElementId, Value> = BTreeMap::new();
    for (id, element) in molecule.stereo_elements() {
        let (kind, focus, carriers, parity) = match &element.kind {
            StereoElementKind::Tetrahedral(s) => {
                let mut carriers = s.carriers.iter().map(carrier).collect::<Vec<_>>();
                let parity = s.orientation.map(|o| {
                    u8::from(o == TetrahedralOrientation::CounterClockwise)
                        ^ permutation_parity(&carriers)
                });
                carriers.sort_unstable();
                (
                    "tetrahedral",
                    vec![component.atom(s.center)],
                    carriers,
                    parity,
                )
            }
            StereoElementKind::DoubleBond(s) => {
                let selected = [carrier(&s.left_carrier), carrier(&s.right_carrier)];
                let mut reference = [
                    reference_carrier(component, s.left, s.right)?,
                    reference_carrier(component, s.right, s.left)?,
                ];
                let parity = s.orientation.map(|o| {
                    u8::from(o == DoubleBondOrientation::Opposite)
                        ^ u8::from(selected[0] != reference[0])
                        ^ u8::from(selected[1] != reference[1])
                });
                let mut focus = [component.atom(s.left), component.atom(s.right)];
                if focus[0] > focus[1] {
                    focus.swap(0, 1);
                    reference.swap(0, 1);
                }
                ("double_bond", focus.to_vec(), reference.to_vec(), parity)
            }
            StereoElementKind::Axis(s) => {
                let bond = molecule
                    .bond(s.axis)
                    .map_err(|e| Failure::new("stereo", e))?;
                let selected = s.carriers.iter().map(carrier).collect::<Vec<_>>();
                let mut reference = [
                    reference_carrier(component, bond.a(), bond.b())?,
                    reference_carrier(component, bond.b(), bond.a())?,
                ];
                if selected.len() != 2 {
                    return Err(Failure::new("stereo", "axis without two endpoint carriers"));
                }
                let parity = s.orientation.map(|o| {
                    u8::from(o == AxisOrientation::CounterClockwise)
                        ^ u8::from(selected[0] != reference[0])
                        ^ u8::from(selected[1] != reference[1])
                });
                let mut focus = [component.atom(bond.a()), component.atom(bond.b())];
                if focus[0] > focus[1] {
                    focus.swap(0, 1);
                    reference.swap(0, 1);
                }
                ("axis", focus.to_vec(), reference.to_vec(), parity)
            }
        };
        facts.push(json!(["stereo", kind, focus, "carriers", carriers]));
        facts.push(json!(["stereo", kind, focus, "parity", parity]));
        members.insert(id, json!([kind, focus]));
    }
    for (_, group) in molecule.stereo_groups() {
        let mut group_members = group
            .members
            .iter()
            .map(|id| members[id].clone())
            .collect::<Vec<_>>();
        group_members.sort_by_key(Value::to_string);
        facts.push(json!([
            "group",
            group_members,
            "kind",
            group_kind(group.kind)
        ]));
    }
    Ok(())
}

fn group_kind(kind: StereoGroupKind) -> &'static str {
    match kind {
        StereoGroupKind::Absolute => "absolute",
        StereoGroupKind::Relative => "relative",
        StereoGroupKind::Racemic => "racemic",
        StereoGroupKind::And => "and",
        StereoGroupKind::Or => "or",
    }
}

fn ring_facts(components: &mut [Component], facts: &mut Facts) -> Result<(), Failure> {
    let mut sizes = BTreeMap::<usize, usize>::new();
    for component in components.iter_mut() {
        let membership = rings::perceive_ring_membership(&mut component.molecule);
        let ring_set = rings::perceive_ring_set(&mut component.molecule)
            .map_err(|e| Failure::new("rings", e))?;
        for ring in ring_set.rings() {
            *sizes.entry(ring.atoms.len()).or_default() += 1;
        }
        let molecule = &component.molecule;
        for id in molecule.atom_ids() {
            facts.push(json!([
                "atom",
                component.atom(id),
                "in_ring",
                membership.atom_in_ring(id)
            ]));
        }
        for (id, bond) in molecule.bonds() {
            let [a, b] = ends(component, bond);
            facts.push(json!([
                "bond",
                a,
                b,
                "in_ring",
                membership.bond_in_ring(id)
            ]));
        }
    }
    for (size, count) in sizes {
        facts.push(json!(["rings", size, "count", count]));
    }
    Ok(())
}

fn aromaticity_facts(
    components: &mut [Component],
    options: &Value,
    facts: &mut Facts,
) -> Result<(), Failure> {
    let model = match options["model"].as_str() {
        Some("rdkit") | None => AromaticityModel::RdkitLike,
        Some("mdl") => AromaticityModel::Mdl,
        Some(other) => return Err(Failure::new("request", format!("unknown model {other:?}"))),
    };
    for component in components.iter_mut() {
        if model != AromaticityModel::RdkitLike {
            aromaticity::perceive_aromaticity(&mut component.molecule, model)
                .map_err(|e| Failure::new("aromaticity", e))?;
        }
        let molecule = &component.molecule;
        for id in molecule.atom_ids() {
            let aromatic = molecule.atom_is_aromatic(id).ok().flatten();
            facts.push(json!(["atom", component.atom(id), "aromatic", aromatic]));
        }
        for (id, bond) in molecule.bonds() {
            let [a, b] = ends(component, bond);
            let aromatic = molecule.bond_is_aromatic(id).ok().flatten();
            facts.push(json!(["bond", a, b, "aromatic", aromatic]));
        }
    }
    Ok(())
}

fn conjugation(components: &[Component], facts: &mut Facts) -> Result<(), Failure> {
    for component in components {
        let molecule = &component.molecule;
        for (id, bond) in molecule.bonds() {
            let [a, b] = ends(component, bond);
            let conjugated = molecule.bond_is_conjugated(id).ok().flatten();
            facts.push(json!(["bond", a, b, "conjugated", conjugated]));
        }
    }
    Ok(())
}

/// Conjugated groups and default-option contributors, per component. A
/// contributor is identified by its nonzero charges and non-single bonds.
fn resonance_facts(components: &mut [Component], facts: &mut Facts) -> Result<(), Failure> {
    for component in components.iter_mut() {
        resonance::perceive_resonance(&mut component.molecule)
            .map_err(|e| Failure::new("resonance", e))?;
        let component = &*component;
        let molecule = &component.molecule;
        let state = molecule
            .perception()
            .resonance_state()
            .ok_or_else(|| Failure::new("resonance", "no resonance state after perception"))?;
        for group in state.groups() {
            let mut atoms = group
                .atoms
                .iter()
                .map(|&a| component.atom(a))
                .collect::<Vec<_>>();
            atoms.sort_unstable();
            let mut bonds = group
                .bonds
                .iter()
                .map(|&id| molecule.bond(id).map(|bond| ends(component, bond)))
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| Failure::new("resonance", e))?;
            bonds.sort_unstable();
            facts.push(json!(["resonance", "group", atoms[0], "atoms", atoms]));
            facts.push(json!(["resonance", "group", atoms[0], "bonds", bonds]));
        }
        let structures = resonance::enumerate_resonance(molecule, Default::default())
            .map_err(|e| Failure::new("resonance", e))?;
        let mut contributors = BTreeSet::new();
        for contributor in structures.contributors() {
            let mut charges = contributor
                .formal_charges()
                .filter(|&(_, q)| q != 0)
                .map(|(id, q)| (component.atom(id), q))
                .collect::<Vec<_>>();
            charges.sort_unstable();
            let mut bonds = Vec::new();
            for (id, bond_order) in contributor.bond_orders() {
                if bond_order != BondOrder::Single {
                    let bond = molecule
                        .bond(id)
                        .map_err(|e| Failure::new("resonance", e))?;
                    bonds.push((ends(component, bond), order(bond_order).to_string()));
                }
            }
            bonds.sort_unstable();
            let charges = charges
                .iter()
                .map(|(i, q)| format!("{i}:{q}"))
                .collect::<Vec<_>>();
            let bonds = bonds
                .iter()
                .map(|([a, b], o)| format!("{a}-{b}:{o}"))
                .collect::<Vec<_>>();
            contributors.insert(format!("{}|{}", charges.join(","), bonds.join(",")));
        }
        let first = min_atom(component);
        facts.push(json!([
            "resonance",
            "contributors",
            first,
            "count",
            structures.contributors().len()
        ]));
        facts.push(json!([
            "resonance",
            "contributors",
            first,
            "set",
            contributors
        ]));
    }
    Ok(())
}

fn candidates(components: &[Component], facts: &mut Facts) -> Result<(), Failure> {
    for component in components {
        let molecule = &component.molecule;
        for candidate in
            stereo::detect_stereo_candidates(molecule).map_err(|e| Failure::new("stereo", e))?
        {
            match candidate {
                StereoCandidate::Tetrahedral { center, .. } => {
                    facts.push(json!([
                        "candidate",
                        "tetrahedral",
                        [component.atom(center)],
                        "present",
                        true
                    ]));
                }
                StereoCandidate::DoubleBond { bond, .. } => {
                    let bond = molecule.bond(bond).map_err(|e| Failure::new("stereo", e))?;
                    facts.push(json!([
                        "candidate",
                        "double_bond",
                        ends(component, bond),
                        "present",
                        true
                    ]));
                }
            }
        }
    }
    Ok(())
}

fn descriptor(descriptor: StereoDescriptor) -> &'static str {
    match descriptor {
        StereoDescriptor::R => "R",
        StereoDescriptor::S => "S",
        StereoDescriptor::LowerR => "r",
        StereoDescriptor::LowerS => "s",
        // RDKit's CIPLabeler writes the pseudoasymmetric double-bond
        // descriptors seqTrans and seqCis as lowercase e and z.
        StereoDescriptor::SeqTrans => "e",
        StereoDescriptor::SeqCis => "z",
        StereoDescriptor::E => "E",
        StereoDescriptor::Z => "Z",
        StereoDescriptor::M => "M",
        StereoDescriptor::P => "P",
        StereoDescriptor::LowerM => "m",
        StereoDescriptor::LowerP => "p",
    }
}

/// CIP labels by focus: atoms for centres, sorted endpoints for double bonds
/// and atropisomeric axes alike, as RDKit labels both on the bond.
fn cip(components: &mut [Component], facts: &mut Facts) -> Result<(), Failure> {
    for component in components.iter_mut() {
        stereo::assign_cip_descriptors(&mut component.molecule)
            .map_err(|e| Failure::new("cip", e))?;
        let component = &*component;
        let molecule = &component.molecule;
        for (id, element) in molecule.stereo_elements() {
            let Some(label) = molecule
                .cip_descriptor(id)
                .map_err(|e| Failure::new("cip", e))?
            else {
                continue;
            };
            let (kind, focus) = match &element.kind {
                StereoElementKind::Tetrahedral(s) => ("atom", vec![component.atom(s.center)]),
                StereoElementKind::DoubleBond(s) => {
                    let mut focus = vec![component.atom(s.left), component.atom(s.right)];
                    focus.sort_unstable();
                    ("bond", focus)
                }
                StereoElementKind::Axis(s) => {
                    let bond = molecule.bond(s.axis).map_err(|e| Failure::new("cip", e))?;
                    ("bond", ends(component, bond).to_vec())
                }
            };
            facts.push(json!(["cip", kind, focus, "label", descriptor(label)]));
        }
    }
    Ok(())
}

/// Atom equivalence classes within each component.
fn symmetry(components: &[Component], facts: &mut Facts) {
    for component in components {
        let ranking = canon::atom_ranking(&component.molecule);
        let mut classes = BTreeMap::<u32, Vec<usize>>::new();
        for (atom, rank) in ranking.iter() {
            classes.entry(rank).or_default().push(component.atom(atom));
        }
        for mut atoms in classes.into_values() {
            atoms.sort_unstable();
            facts.push(json!(["class", atoms[0], "atoms", atoms]));
        }
    }
}

fn formula(components: &[Component], facts: &mut Facts) -> Result<(), Failure> {
    for component in components {
        let formula = molecular_formula(&component.molecule, HydrogenCountPolicy::IncludePerceived)
            .map_err(|e| Failure::new("formula", e))?;
        let first = min_atom(component);
        for (element, isotope, count) in formula.terms() {
            facts.push(json!([
                "formula",
                first,
                "count",
                element.symbol(),
                isotope,
                count
            ]));
        }
        facts.push(json!(["formula", first, "charge", formula.formal_charge()]));
    }
    Ok(())
}

fn rotatable(components: &[Component], facts: &mut Facts) -> Result<(), Failure> {
    for component in components {
        let molecule = &component.molecule;
        let detected = rotatable_bonds::detect(molecule, RotatableBondOptions::STRICT)
            .map_err(|e| Failure::new("rotatable", e))?;
        for &id in detected.bond_ids() {
            let bond = molecule
                .bond(id)
                .map_err(|e| Failure::new("rotatable", e))?;
            let [a, b] = ends(component, bond);
            facts.push(json!(["bond", a, b, "rotatable", true]));
        }
    }
    Ok(())
}

/// Every query in `options.queries` against this record as one target. Each
/// match lists source atom indices in query atom order.
fn smarts(components: Vec<Component>, options: &Value, facts: &mut Facts) -> Result<(), Failure> {
    let queries = options["queries"]
        .as_array()
        .ok_or_else(|| Failure::new("request", "smarts needs options.queries"))?;
    let sources = components
        .iter()
        .map(|component| component.source.clone())
        .collect::<Vec<_>>();
    let topology = Arc::new(
        Topology::from_molecules(components.into_iter().map(|c| c.molecule))
            .map_err(|e| Failure::new("smarts", e))?,
    );
    let target = PreparedTopologyTarget::new(&topology);
    let match_options = SubstructureMatchOptions {
        max_matches: MAX_MATCHES,
        uniquify: false,
        ..Default::default()
    };
    for (index, smarts) in queries.iter().enumerate() {
        let smarts = smarts.as_str().unwrap_or_default();
        let graph = match query::parse_smarts(smarts) {
            Ok(graph) => graph,
            Err(_) => {
                facts.push(json!(["query", index, "error", "parse"]));
                continue;
            }
        };
        let matches = match target.find_matches_with_options(&graph, match_options) {
            Ok(matches) => matches,
            Err(_) => {
                facts.push(json!(["query", index, "error", "match"]));
                continue;
            }
        };
        let mut mappings = matches
            .iter()
            .map(|found| {
                found
                    .atoms()
                    .iter()
                    .map(|atom| sources[atom.molecule().index()][atom.atom().index()])
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        mappings.sort_unstable();
        facts.push(json!(["query", index, "count", mappings.len()]));
        facts.push(json!(["query", index, "matches", mappings]));
    }
    Ok(())
}

/// Kekule's serialization of the input; the reference reads it back.
fn write(format: &str, text: &str, options: &Value) -> Result<Facts, Failure> {
    let target = options["format"].as_str().unwrap_or_default();
    let written = match target {
        "smiles-isomeric" | "smiles-canonical" => {
            let mut components = input::read(format, text)?;
            for component in &mut components {
                prepare(component)?;
            }
            let topology = Topology::from_molecules(components.into_iter().map(|c| c.molecule))
                .map_err(|e| Failure::new("write", e))?;
            let mode = if target == "smiles-canonical" {
                smiles::SmilesWriteMode::Canonical
            } else {
                smiles::SmilesWriteMode::Isomeric
            };
            smiles::write(&topology, smiles::SmilesWriteOptions { mode })
                .map_err(|e| Failure::new("write", e))?
        }
        "molfile-v2000" | "molfile-v3000" | "sdf" => {
            let model = match format {
                "mol" => molfile::parse_str(text)
                    .map_err(|e| Failure::new("parse", e))?
                    .interpret()
                    .map_err(|e| Failure::new("parse", e))?
                    .into_model(),
                "sdf" => {
                    let document = sdf::parse_str(text).map_err(|e| Failure::new("parse", e))?;
                    let [record] = document.records() else {
                        return Err(Failure::new("request", "expected one SDF record"));
                    };
                    record
                        .interpret()
                        .map_err(|e| Failure::new("parse", e))?
                        .into_model()
                }
                other => {
                    return Err(Failure::new(
                        "request",
                        format!("{target} cannot read {other}"),
                    ))
                }
            };
            if target == "sdf" {
                let record = sdf::SdfRecordInterpretation::new("", model, vec![]);
                sdf::write(
                    &[record],
                    sdf::SdfWriteOptions {
                        version: sdf::MolfileWriteVersion::V2000,
                    },
                )
                .map_err(|e| Failure::new("write", e))?
            } else {
                let version = if target == "molfile-v3000" {
                    molfile::MolfileWriteVersion::V3000
                } else {
                    molfile::MolfileWriteVersion::V2000
                };
                molfile::write(&model, molfile::MolfileWriteOptions { version })
                    .map_err(|e| Failure::new("write", e))?
            }
        }
        other => {
            return Err(Failure::new(
                "request",
                format!("unknown write format {other:?}"),
            ))
        }
    };
    Ok(vec![json!(["written", target, written])])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn facts(task: &str, text: &str, options: Value) -> Vec<String> {
        observe(task, "smiles", text, &options)
            .unwrap()
            .iter()
            .map(Value::to_string)
            .collect()
    }

    #[test]
    fn parse_keys_atoms_by_source_position_across_components() {
        let observed = facts("parse", "[Na+].C[C@H](N)C(=O)[O-]", Value::Null);
        let expected = [
            r#"["atom",0,"element","Na"]"#,
            r#"["atom",0,"isotope",null]"#,
            r#"["atom",0,"charge",1]"#,
            r#"["atom",0,"radicals",0]"#,
            r#"["atom",0,"hydrogens",0]"#,
            r#"["atom",0,"map",null]"#,
            r#"["component",0,"atoms",[0]]"#,
        ];
        assert_eq!(&observed[..7], expected);
        assert!(observed.contains(&r#"["component",1,"atoms",[1,2,3,4,5,6]]"#.to_owned()));
        assert!(observed.contains(&r#"["atom",2,"hydrogens",1]"#.to_owned()));
        assert!(observed.contains(&r#"["bond",4,5,"order",2]"#.to_owned()));
        assert!(
            observed.contains(&r#"["stereo","tetrahedral",[2],"carriers",[-1,1,3,4]]"#.to_owned())
        );
        assert!(observed.contains(&r#"["stereo","tetrahedral",[2],"parity",1]"#.to_owned()));
    }

    #[test]
    fn rings_report_membership_and_ring_sizes() {
        let observed = facts("rings", "C1CC1C", Value::Null);
        assert_eq!(
            observed,
            [
                r#"["atom",0,"in_ring",true]"#,
                r#"["atom",1,"in_ring",true]"#,
                r#"["atom",2,"in_ring",true]"#,
                r#"["atom",3,"in_ring",false]"#,
                r#"["bond",0,1,"in_ring",true]"#,
                r#"["bond",1,2,"in_ring",true]"#,
                r#"["bond",0,2,"in_ring",true]"#,
                r#"["bond",2,3,"in_ring",false]"#,
                r#"["rings",3,"count",1]"#,
            ]
        );
    }

    #[test]
    fn cip_labels_centres_by_source_atom() {
        assert_eq!(
            facts("stereo.cip", "C[C@H](N)C(=O)O", Value::Null),
            [r#"["cip","atom",[1],"label","S"]"#]
        );
    }

    #[test]
    fn smarts_reports_matches_and_unparsable_queries() {
        let observed = facts(
            "smarts",
            "CCO",
            json!({"queries": ["[OX2H]", "C~C", "[bad"]}),
        );
        assert_eq!(
            observed,
            [
                r#"["query",0,"count",1]"#,
                r#"["query",0,"matches",[[2]]]"#,
                r#"["query",1,"count",2]"#,
                r#"["query",1,"matches",[[0,1],[1,0]]]"#,
                r#"["query",2,"error","parse"]"#,
            ]
        );
    }

    #[test]
    fn write_returns_the_serialized_text() {
        assert_eq!(
            facts("write", "F/C=C/Cl", json!({"format": "smiles-canonical"})),
            [r#"["written","smiles-canonical","F/C=C/Cl"]"#]
        );
    }
}
