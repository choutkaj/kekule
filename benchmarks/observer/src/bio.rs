//! Biomolecule tasks over mmCIF. Atom sites are keyed by `_atom_site.id`.

use std::collections::BTreeMap;
use std::sync::Arc;

use kekule::core::BondOrder;
use kekule::dssp::{self, DsspHydrogenBond, DsspOptions};
use kekule::mmcif::{
    self, MmcifDocument, MmcifEnsembleInterpretOptions, MmcifEntry, MmcifInterpretOptions,
    MmcifModelSelection, MmcifValue,
};
use kekule::topology::{AtomSelection, InstanceAtomId, ResidueClass, ResidueId};
use kekule::units::{ANGSTROM, SQUARE_ANGSTROM};
use serde_json::{json, Value};

use crate::{Facts, Failure};

pub fn observe(task: &str, format: &str, text: &str, _options: &Value) -> Result<Facts, Failure> {
    if format != "mmcif" {
        return Err(Failure::new(
            "request",
            format!("bio tasks read mmCIF, not {format:?}"),
        ));
    }
    let document = mmcif::parse_str(text).map_err(|e| Failure::new("parse", e))?;
    match task {
        "cif.syntax" => syntax(&document),
        "hierarchy" => hierarchy(&document),
        "connectivity" => connectivity(&document),
        "classification" => classification(&document),
        "dssp" => secondary_structure(&document),
        "superposition" => superposition(&document),
        other => Err(Failure::new("request", format!("unknown task bio.{other}"))),
    }
}

/// FNV-1a over values separated by 0x1f; the Python references use the same.
fn fingerprint(values: impl Iterator<Item = String>) -> String {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for (index, value) in values.enumerate() {
        let separator: &[u8] = if index == 0 { &[] } else { &[0x1f] };
        for &byte in separator.iter().chain(value.as_bytes()) {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0100_0000_01b3);
        }
    }
    format!("{hash:016x}")
}

/// Unquoted `?` and `.` are null values, distinct from the quoted text.
fn token(value: &MmcifValue) -> String {
    if value.is_missing() {
        format!("\u{0}{}", value.text())
    } else {
        value.text().to_owned()
    }
}

/// Every item value, and each loop column's length and fingerprint.
fn syntax(document: &MmcifDocument) -> Result<Facts, Failure> {
    let mut facts = Facts::new();
    for (index, block) in document.blocks().iter().enumerate() {
        facts.push(json!(["block", index, "name", block.name()]));
        for entry in block.entries() {
            match entry {
                MmcifEntry::Item(item) => {
                    let tag = item.tag().to_ascii_lowercase();
                    facts.push(json!(["item", index, tag, "value", token(item.value())]));
                }
                MmcifEntry::Loop(table) => {
                    for tag in table.tags() {
                        let column = (0..table.row_count())
                            .map(|row| table.value(row, tag).map(token).unwrap_or_default());
                        let tag = tag.to_ascii_lowercase();
                        facts.push(json!(["column", index, tag, "rows", table.row_count()]));
                        facts.push(json!([
                            "column",
                            index,
                            tag,
                            "fingerprint",
                            fingerprint(column)
                        ]));
                    }
                }
            }
        }
    }
    Ok(facts)
}

fn first_model(document: &MmcifDocument) -> Result<mmcif::MmcifInterpretation, Failure> {
    mmcif::interpret(
        document,
        MmcifInterpretOptions {
            model_selection: MmcifModelSelection::First,
            ..Default::default()
        },
    )
    .map_err(|e| Failure::new("interpret", e))
}

/// Site identity for every atom an interpretation kept.
fn site_ids(interpretation: &mmcif::MmcifInterpretation) -> BTreeMap<InstanceAtomId, String> {
    interpretation
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .map(|atom| (atom.atom(), atom.atom_site_id().unwrap_or("?").to_owned()))
        .collect()
}

/// The atom sites Kekule kept from the first model: identity, element,
/// occupancy, B in Å² and coordinates in Å.
fn hierarchy(document: &MmcifDocument) -> Result<Facts, Failure> {
    let interpretation = first_model(document)?;
    let model = interpretation.model();
    let mut facts = Facts::new();
    for instance in interpretation.report().instances() {
        for atom in instance.atoms() {
            let view = model
                .atom(atom.atom())
                .ok_or_else(|| Failure::new("interpret", "reported atom missing from model"))?;
            let point = view
                .position()
                .value_in(ANGSTROM)
                .map_err(|e| Failure::new("units", e))?;
            let b_factor = view
                .b_factor()
                .map(|b| b.value_in(SQUARE_ANGSTROM))
                .transpose()
                .map_err(|e| Failure::new("units", e))?;
            facts.push(json!([
                "site",
                atom.atom_site_id().unwrap_or("?"),
                [
                    atom.label_asym_id(),
                    atom.auth_asym_id(),
                    atom.label_sequence_id(),
                    atom.author_sequence_id(),
                    atom.insertion_code(),
                    atom.label_component_id(),
                    atom.label_atom_name(),
                    view.atom().element.symbol().to_ascii_uppercase(),
                    atom.selected_alternate_location(),
                    view.occupancy(),
                    b_factor,
                    [point.x, point.y, point.z],
                ]
            ]));
        }
    }
    Ok(facts)
}

/// Covalent bonds between kept sites, with orders, and the molecule partition.
fn connectivity(document: &MmcifDocument) -> Result<Facts, Failure> {
    let interpretation = first_model(document)?;
    let ids = site_ids(&interpretation);
    let topology = interpretation.topology();
    let mut facts = Facts::new();
    for instance in topology.molecules() {
        let molecule = instance.molecule();
        let id = |atom| ids[&InstanceAtomId::new(instance.id(), atom)].clone();
        for (_, bond) in molecule.bonds() {
            let mut ends = [id(bond.a()), id(bond.b())];
            if bond.order != BondOrder::Dative {
                ends.sort_by_key(|id| site_order(id));
            }
            let order = match bond.order {
                BondOrder::Zero => json!(0),
                BondOrder::Single => json!(1),
                BondOrder::Double => json!(2),
                BondOrder::Triple => json!(3),
                BondOrder::Quadruple => json!(4),
                BondOrder::Dative => json!("dative"),
            };
            facts.push(json!(["bond", ends[0], ends[1], "order", order]));
        }
        let mut sites = molecule.atom_ids().map(id).collect::<Vec<_>>();
        sites.sort_by_key(|id| site_order(id));
        facts.push(json!(["molecule", sites[0], "sites", sites]));
    }
    Ok(facts)
}

fn residue_class(class: ResidueClass) -> &'static str {
    match class {
        ResidueClass::AminoAcid => "amino-acid",
        ResidueClass::DnaNucleotide => "dna",
        ResidueClass::RnaNucleotide => "rna",
        ResidueClass::Carbohydrate => "carbohydrate",
        ResidueClass::Water => "water",
        ResidueClass::Ion => "ion",
        _ => "other",
    }
}

/// The class of each residue, keyed `chain/number+insertion/component` by
/// author identifiers, which do not depend on which altloc rows were kept.
fn classification(document: &MmcifDocument) -> Result<Facts, Failure> {
    let interpretation = first_model(document)?;
    let model = interpretation.model();
    let mut residues: BTreeMap<ResidueId, (String, &'static str)> = BTreeMap::new();
    for instance in interpretation.report().instances() {
        for atom in instance.atoms() {
            let Some(residue) = model.atom(atom.atom()).and_then(|view| view.residue()) else {
                continue;
            };
            residues.entry(residue.id()).or_insert_with(|| {
                let key = format!(
                    "{}/{}{}/{}",
                    atom.auth_asym_id().unwrap_or(atom.asym_id()),
                    atom.author_sequence_id().unwrap_or("?"),
                    atom.insertion_code().unwrap_or(""),
                    atom.auth_component_id().unwrap_or(atom.component_id()),
                );
                (key, residue_class(residue.class()))
            });
        }
    }
    let mut facts = residues
        .into_values()
        .map(|(key, class)| json!(["residue", key, "class", class]))
        .collect::<Facts>();
    facts.sort_by_key(Value::to_string);
    Ok(facts)
}

/// Numeric site IDs order numerically; others after them, lexically.
fn site_order(id: &str) -> (u8, u64, String) {
    match id.parse::<u64>() {
        Ok(number) => (0, number, String::new()),
        Err(_) => (1, 0, id.to_owned()),
    }
}

fn partner(
    bond: Option<DsspHydrogenBond>,
    residues: &BTreeMap<ResidueId, &dssp::DsspResidue>,
) -> (Value, Value) {
    match bond {
        Some(bond) => {
            let source = residues[&bond.partner].source();
            (
                json!(format!(
                    "{}:{}",
                    source.chain_label_id,
                    source
                        .label_sequence_id
                        .map_or("?".to_owned(), |s| s.to_string())
                )),
                json!(bond.energy_kcal_per_mol),
            )
        }
        None => (Value::Null, Value::Null),
    }
}

/// DSSP per residue of a single-model, single-conformer input, keyed by
/// label chain and sequence number as mkdssp reports them.
fn secondary_structure(document: &MmcifDocument) -> Result<Facts, Failure> {
    let interpretation = first_model(document)?;
    let result = match dssp::assign(
        interpretation.model().as_model_view(),
        DsspOptions::default(),
    ) {
        Ok(result) => result,
        Err(dssp::DsspError::NoAnalyzableProteinResidues) => return Ok(Facts::new()),
        Err(error) => return Err(Failure::new("dssp", error)),
    };
    let residues = result
        .residues()
        .map(|r| (r.key(), r))
        .collect::<BTreeMap<_, _>>();
    let mut facts = Facts::new();
    for residue in result.residues() {
        let source = residue.source();
        let code = residue.secondary_structure().code();
        let [acceptor_1, acceptor_2] = residue.acceptors().map(|b| partner(b, &residues));
        let [donor_1, donor_2] = residue.donors().map(|b| partner(b, &residues));
        facts.push(json!([
            "dssp",
            source.chain_label_id,
            source.label_sequence_id,
            [
                if code == ' ' {
                    ".".to_owned()
                } else {
                    code.to_string()
                },
                residue.phi_degrees(),
                residue.psi_degrees(),
                residue.kappa_degrees(),
                residue.alpha_degrees(),
                residue.tco(),
                acceptor_1.0,
                acceptor_1.1,
                acceptor_2.0,
                acceptor_2.1,
                donor_1.0,
                donor_1.1,
                donor_2.0,
                donor_2.1,
            ]
        ]));
    }
    Ok(facts)
}

/// Cα RMSD of every coordinate model after fitting it onto the first.
fn superposition(document: &MmcifDocument) -> Result<Facts, Failure> {
    let interpretation =
        mmcif::interpret_ensemble(document, MmcifEnsembleInterpretOptions::default())
            .map_err(|e| Failure::new("interpret", e))?;
    let ensemble = interpretation.ensemble();
    let topology: Arc<_> = ensemble.shared_topology();
    let first = &interpretation.reports()[0];
    let alpha = first
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .filter(|atom| atom.label_atom_name() == Some("CA"))
        .filter(|atom| {
            topology
                .atom(atom.atom())
                .is_some_and(|view| view.atom().element.symbol() == "C")
        })
        .map(|atom| atom.atom())
        .collect::<Vec<_>>();
    let count = alpha.len();
    let selection =
        AtomSelection::from_atoms(&topology, alpha).map_err(|e| Failure::new("selection", e))?;
    let rmsd = ensemble
        .aligned_rmsd(0usize, &selection, &selection)
        .map_err(|e| Failure::new("superposition", e))?
        .value_in(ANGSTROM)
        .map_err(|e| Failure::new("units", e))?;
    let mut facts = vec![json!(["alpha", "count", count])];
    for (report, value) in interpretation.reports().iter().zip(rmsd) {
        facts.push(json!(["model", report.selected_model(), "rmsd", value]));
    }
    Ok(facts)
}
