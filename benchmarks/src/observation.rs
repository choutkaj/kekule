//! Schema validation at both observation boundaries. These types deliberately
//! carry no chemistry algorithms: they require measured fields and their types.
#![allow(dead_code)]
use crate::{boxed_error, features};
use serde::Deserialize;
use serde_json::Value;
use std::{collections::BTreeMap, error::Error};

#[derive(Deserialize)]
#[serde(rename_all = "snake_case")]
enum Success {
    Ok,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Records<T> {
    records: Vec<T>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GraphRecord {
    record_index: usize,
    status: Success,
    title: String,
    components: Vec<Graph>,
    properties: Vec<Property>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Property {
    name: String,
    value: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Atom {
    index: usize,
    atomic_number: u16,
    symbol: String,
    formal_charge: i32,
    #[serde(deserialize_with = "required_option")]
    isotope: Option<u16>,
    explicit_hydrogens: u32,
    #[serde(deserialize_with = "required_option")]
    atom_map: Option<u32>,
    #[serde(deserialize_with = "required_option")]
    radical: Option<String>,
    unpaired_electrons: u32,
    aromatic: bool,
    implicit_hydrogens: u32,
    explicit_valence: u32,
    no_implicit_hydrogens: bool,
    coord: Option<[f64; 3]>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Bond {
    begin_atom_index: usize,
    end_atom_index: usize,
    bond_type: String,
    is_aromatic: bool,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Focus {
    #[serde(rename = "type")]
    kind: String,
    focus: Vec<usize>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Stereo {
    #[serde(rename = "type")]
    kind: String,
    focus: Vec<usize>,
    carriers: Vec<i64>,
    #[serde(deserialize_with = "required_option")]
    parity: Option<u8>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Group {
    kind: String,
    members: Vec<Focus>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Graph {
    atom_count: usize,
    bond_count: usize,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    stereo: Vec<Stereo>,
    groups: Vec<Group>,
    candidates: Option<Vec<Focus>>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Identity {
    record_index: usize,
    status: Success,
    title: String,
    identity: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Query {
    record_index: usize,
    status: Success,
    title: String,
    smarts: String,
    atom_count: usize,
    bond_count: usize,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Matches {
    smarts: String,
    matches: Vec<Vec<usize>>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Search {
    record_index: usize,
    status: Success,
    title: String,
    queries: Vec<Matches>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BondValue<T> {
    begin_atom_index: usize,
    end_atom_index: usize,
    value: T,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RingMembership {
    record_index: usize,
    status: Success,
    title: String,
    atom_in_ring: Vec<bool>,
    bond_in_ring: Vec<BondValue<bool>>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RingSet {
    record_index: usize,
    status: Success,
    title: String,
    rings: Vec<Vec<usize>>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Aromaticity {
    record_index: usize,
    status: Success,
    title: String,
    atom_aromatic: Vec<bool>,
    bond_aromatic: Vec<BondValue<bool>>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ValenceAtom {
    index: usize,
    atomic_number: u16,
    symbol: String,
    formal_charge: i32,
    explicit_hydrogens: u32,
    implicit_hydrogens: u32,
    explicit_valence: u32,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Valence {
    record_index: usize,
    status: Success,
    title: String,
    atoms: Vec<ValenceAtom>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BasicAtom {
    index: usize,
    atomic_number: u16,
    symbol: String,
    formal_charge: i32,
    #[serde(deserialize_with = "required_option")]
    isotope: Option<u16>,
    explicit_hydrogens: u32,
    #[serde(deserialize_with = "required_option")]
    atom_map: Option<u32>,
    aromatic: bool,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Perception {
    record_index: usize,
    status: Success,
    title: String,
    atoms: Vec<BasicAtom>,
    graph: Graph,
    valence: Vec<ValenceAtom>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Classes {
    record_index: usize,
    status: Success,
    title: String,
    classes: Vec<Vec<usize>>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Term {
    element: String,
    #[serde(deserialize_with = "required_option")]
    isotope: Option<u16>,
    count: u64,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Formula {
    terms: Vec<Term>,
    formal_charge: i32,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MolecularDescriptor {
    record_index: usize,
    status: Success,
    title: String,
    formula: Formula,
    average_mass_da: f64,
    monoisotopic_mass_da: f64,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Endpoints {
    begin_atom_index: usize,
    end_atom_index: usize,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Rotatable {
    record_index: usize,
    status: Success,
    title: String,
    count: usize,
    bonds: Vec<Endpoints>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct AddedHydrogen {
    parent_atom_index: usize,
    count: usize,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Hydrogens {
    record_index: usize,
    status: Success,
    title: String,
    atom_count_after_add: usize,
    added_graph: Graph,
    added_hydrogens_by_parent: Vec<AddedHydrogen>,
    round_trip: Graph,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct AtomDescriptor {
    atom_index: usize,
    descriptor: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BondDescriptor {
    begin_atom_index: usize,
    end_atom_index: usize,
    descriptor: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Cip {
    record_index: usize,
    status: Success,
    title: String,
    atom_count: usize,
    bond_count: usize,
    atom_descriptors: Vec<AtomDescriptor>,
    bond_descriptors: Vec<BondDescriptor>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Block {
    name: String,
    values: BTreeMap<String, Vec<String>>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Cif {
    blocks: Vec<Block>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Partner {
    partner_chain_id: String,
    partner_sequence_id: i32,
    #[serde(deserialize_with = "required_option")]
    partner_insertion_code: Option<String>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct HydrogenBond {
    partner_chain_id: String,
    partner_sequence_id: i32,
    #[serde(deserialize_with = "required_option")]
    partner_insertion_code: Option<String>,
    energy_kcal_per_mol: f64,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Residue {
    chain_id: String,
    sequence_id: i32,
    #[serde(deserialize_with = "required_option")]
    insertion_code: Option<String>,
    label_chain_id: String,
    #[serde(deserialize_with = "required_option")]
    label_sequence_id: Option<i32>,
    residue_name: String,
    residue_one_letter: String,
    secondary_structure: String,
    chain_break: String,
    #[serde(deserialize_with = "required_option")]
    phi_degrees: Option<f64>,
    #[serde(deserialize_with = "required_option")]
    psi_degrees: Option<f64>,
    #[serde(deserialize_with = "required_option")]
    omega_degrees: Option<f64>,
    #[serde(deserialize_with = "required_option")]
    tco: Option<f64>,
    #[serde(deserialize_with = "required_option")]
    kappa_degrees: Option<f64>,
    #[serde(deserialize_with = "required_option")]
    alpha_degrees: Option<f64>,
    helix_positions: [String; 4],
    #[serde(deserialize_with = "required_option")]
    sheet: Option<usize>,
    #[serde(deserialize_with = "required_option")]
    strand: Option<usize>,
    ladders: [Option<usize>; 2],
    beta_parallel: [Option<bool>; 2],
    beta_partners: [Option<Partner>; 2],
    acceptors: [Option<HydrogenBond>; 2],
    donors: [Option<HydrogenBond>; 2],
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Dssp {
    status: Success,
    residues: Vec<Residue>,
}

fn read<'a, T: Deserialize<'a>>(value: &'a Value) -> Result<(), Box<dyn Error>> {
    T::deserialize(value)?;
    Ok(())
}
pub(crate) fn validate(feature: &str, value: &Value) -> Result<(), Box<dyn Error>> {
    graph_dimensions(value, feature == "stereo.perception")?;
    let collection = match feature {
        "io.mmcif.parse" => "blocks",
        "bio.secondary-structure.dssp" => "residues",
        _ => "records",
    };
    if value
        .get(collection)
        .and_then(Value::as_array)
        .is_none_or(Vec::is_empty)
    {
        return Err(boxed_error(format!("missing or empty {collection}")));
    }
    match feature {
        "io.mmcif.parse" => read::<Cif>(value),
        "bio.secondary-structure.dssp" => read::<Dssp>(value),
        "query.smarts" => read::<Records<Query>>(value),
        "algo.substructure.vf2" => read::<Records<Search>>(value),
        "algo.rings.fast" => read::<Records<RingMembership>>(value),
        "algo.rings.sssr" => read::<Records<RingSet>>(value),
        "algo.aromaticity.rdkit-like" => read::<Records<Aromaticity>>(value),
        "algo.valence.rdkit-like" => read::<Records<Valence>>(value),
        "algo.canonical-ranking" => read::<Records<Classes>>(value),
        "chem.perception.default" => read::<Records<Perception>>(value),
        "chem.hydrogen-transforms" => read::<Records<Hydrogens>>(value),
        "descriptor.molecular" => read::<Records<MolecularDescriptor>>(value),
        "descriptor.rotatable-bonds.rdkit-strict" => read::<Records<Rotatable>>(value),
        "stereo.cip" => read::<Records<Cip>>(value),
        _ if features::is_writer(feature) && feature.starts_with("io.smiles.") => {
            read::<Records<Identity>>(value)
        }
        _ => read::<Records<GraphRecord>>(value),
    }
}

fn graph_dimensions(value: &Value, require_candidates: bool) -> Result<(), Box<dyn Error>> {
    match value {
        Value::Object(fields) => {
            if let Some(atoms) = fields.get("atoms").and_then(Value::as_array) {
                if let Some(count) = fields.get("atom_count") {
                    if atoms.is_empty() || count.as_u64() != Some(atoms.len() as u64) {
                        return Err(boxed_error("invalid graph atom count"));
                    }
                    let bonds = fields
                        .get("bonds")
                        .and_then(Value::as_array)
                        .ok_or("missing graph bonds")?;
                    if fields.get("bond_count").and_then(Value::as_u64) != Some(bonds.len() as u64)
                    {
                        return Err(boxed_error("invalid graph bond count"));
                    }
                    if require_candidates && !fields.get("candidates").is_some_and(Value::is_array)
                    {
                        return Err(boxed_error("missing stereo candidates"));
                    }
                }
            }
            if fields
                .get("components")
                .and_then(Value::as_array)
                .is_some_and(Vec::is_empty)
            {
                return Err(boxed_error("empty molecular components"));
            }
            for child in fields.values() {
                graph_dimensions(child, require_candidates)?;
            }
        }
        Value::Array(values) => {
            for child in values {
                graph_dimensions(child, require_candidates)?;
            }
        }
        _ => (),
    }
    Ok(())
}
pub(crate) fn finite(value: f64) -> Result<(), Box<dyn Error>> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(boxed_error("nonfinite scientific observation"))
    }
}

fn required_option<'de, D: serde::Deserializer<'de>, T: Deserialize<'de>>(
    deserializer: D,
) -> Result<Option<T>, D::Error> {
    Option::deserialize(deserializer)
}
