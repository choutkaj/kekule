use std::collections::{BTreeMap, BTreeSet};

use crate::chemistry::localize_source_aromatic_bonds;
use crate::core::{AtomId, BondOrder, Molecule};
use crate::structure::{AtomData, Ensemble, EnsembleMember, Model, ModelBuilder, Positions};
use crate::topology::{InstanceAtomId, Topology};
use crate::units::{Quantity, ANGSTROM};
use std::sync::Arc;

use super::mmcif_interpret as raw;
use super::{MmcifDataBlock, MmcifDocument, MmcifLoopTable, MmcifValue};

/// Interpreted mmCIF model after authoritative component and polymer
/// connectivity has been materialized into the canonical molecule graphs.
#[derive(Debug, Clone)]
pub struct MmcifInterpretation {
    model: Model,
    report: raw::MmcifInterpretationReport,
}

impl MmcifInterpretation {
    #[cfg(test)]
    pub fn report(&self) -> &raw::MmcifInterpretationReport {
        &self.report
    }

    #[cfg(test)]
    pub fn topology(&self) -> &Topology {
        self.model.topology()
    }

    pub fn into_parts(self) -> (Model, raw::MmcifInterpretationReport) {
        (self.model, self.report)
    }
}

/// Interprets mmCIF syntax first, then completes authoritative covalent
/// connectivity without using coordinate-distance guesses.
pub fn interpret_mmcif(
    document: &MmcifDocument,
    options: raw::MmcifInterpretOptions,
) -> Result<MmcifInterpretation, raw::MmcifInterpretError> {
    let interpretation = raw::interpret_mmcif(document, options)?;
    let (model, mut report) = interpretation.into_parts();
    let (model, pending) = rebuild_model_with_connectivity(document, model, &report)?;
    report.template_bonds_pending = pending;
    Ok(MmcifInterpretation { model, report })
}

#[derive(Debug, Clone)]
pub struct MmcifEnsembleInterpretation {
    ensemble: Ensemble,
    reports: Vec<raw::MmcifInterpretationReport>,
}

impl MmcifEnsembleInterpretation {
    pub fn into_parts(self) -> (Ensemble, Vec<raw::MmcifInterpretationReport>) {
        (self.ensemble, self.reports)
    }
}

/// Completes the shared topology of a parsed mmCIF ensemble once, then rebinds
/// every member to the repaired topology without changing dense atom order.
pub fn interpret_mmcif_ensemble(
    document: &MmcifDocument,
    options: raw::MmcifEnsembleInterpretOptions,
) -> Result<MmcifEnsembleInterpretation, raw::MmcifEnsembleInterpretError> {
    let interpretation = raw::interpret_mmcif_ensemble(document, options)?;
    let (source, mut reports) = interpretation.into_parts();
    let first = source
        .member(0)
        .ok_or(raw::MmcifEnsembleInterpretError::EmptyModelSelection)?;
    let model_id = reports
        .first()
        .and_then(raw::MmcifInterpretationReport::selected_model)
        .unwrap_or("<unknown>")
        .to_owned();
    let prototype = Model::with_atom_data(
        source.shared_topology(),
        first.positions().clone(),
        first.cell().copied(),
        first.atom_data().clone(),
    )
    .map_err(|error| raw::MmcifEnsembleInterpretError::Model {
        model_id: model_id.clone(),
        error: interpret_error(error),
    })?;
    let first_report = reports
        .first()
        .ok_or(raw::MmcifEnsembleInterpretError::EmptyModelSelection)?;
    let (prototype, pending) = rebuild_model_with_connectivity(document, prototype, first_report)
        .map_err(|error| raw::MmcifEnsembleInterpretError::Model {
        model_id: model_id.clone(),
        error,
    })?;
    let topology = prototype.shared_topology();
    let mut ensemble = Ensemble::new(Arc::clone(&topology));

    for member in source.members() {
        let positions = Positions::new(&topology, member.positions().values())
            .map_err(raw::MmcifEnsembleInterpretError::Position)?;
        let mut rebuilt = EnsembleMember::new(positions);
        rebuilt.set_cell(member.cell().copied());
        rebuilt
            .set_weight(member.weight())
            .map_err(|error| raw::MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
        rebuilt
            .set_atom_data(
                rebound_atom_data(member.atom_data(), &topology).map_err(|error| {
                    raw::MmcifEnsembleInterpretError::Model {
                        model_id: model_id.clone(),
                        error,
                    }
                })?,
            )
            .map_err(|error| raw::MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
        rebuilt.props_mut().clone_from(member.props());
        ensemble
            .push(rebuilt)
            .map_err(|error| raw::MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
    }

    for report in &mut reports {
        report.template_bonds_pending = pending;
    }
    Ok(MmcifEnsembleInterpretation { ensemble, reports })
}

#[derive(Debug, Clone)]
struct ComponentBond {
    atom_1: String,
    atom_2: String,
    order: ComponentBondOrder,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ComponentBondOrder {
    Localized(BondOrder),
    Aromatic,
}

impl ComponentBondOrder {
    const fn staged_order(self) -> BondOrder {
        match self {
            Self::Localized(order) => order,
            Self::Aromatic => BondOrder::Single,
        }
    }
}

#[derive(Debug, Clone)]
struct BranchSchemeSite {
    entity_id: String,
    asym_id: String,
    number: i32,
    component_id: String,
    author_sequence_id: Option<String>,
}

#[derive(Debug, Clone)]
struct BranchLink {
    entity_id: String,
    number_1: i32,
    atom_1: String,
    number_2: i32,
    atom_2: String,
    order: ComponentBondOrder,
}

#[derive(Debug, Default)]
struct ConnectivityCatalog {
    component_bonds: BTreeMap<String, Vec<ComponentBond>>,
    polymer_types: BTreeMap<String, String>,
    nonstandard_polymer_linkage: BTreeSet<String>,
    branch_sites: Vec<BranchSchemeSite>,
    branch_links: Vec<BranchLink>,
}

impl ConnectivityCatalog {
    fn from_block(block: &MmcifDataBlock) -> Result<Self, raw::MmcifInterpretError> {
        let mut catalog = Self::default();
        if let Some(table) = block.loop_with_tag("_chem_comp_bond.comp_id") {
            for row in 0..table.row_count() {
                let comp_id = required(table, row, "_chem_comp_bond.comp_id")?;
                let atom_1 = required(table, row, "_chem_comp_bond.atom_id_1")?;
                let atom_2 = required(table, row, "_chem_comp_bond.atom_id_2")?;
                let order = component_bond_order(
                    optional(table, row, "_chem_comp_bond.value_order").unwrap_or("sing"),
                    table,
                    row,
                )?;
                catalog
                    .component_bonds
                    .entry(normalized(comp_id))
                    .or_default()
                    .push(ComponentBond {
                        atom_1: atom_1.to_owned(),
                        atom_2: atom_2.to_owned(),
                        order,
                    });
            }
        }
        if let Some(table) = block.loop_with_tag("_entity_poly.entity_id") {
            for row in 0..table.row_count() {
                let entity = required(table, row, "_entity_poly.entity_id")?;
                let kind = required(table, row, "_entity_poly.type")?;
                catalog
                    .polymer_types
                    .insert(entity.to_owned(), kind.to_owned());
                if optional(table, row, "_entity_poly.nstd_linkage").is_some_and(is_yes) {
                    catalog
                        .nonstandard_polymer_linkage
                        .insert(entity.to_owned());
                }
            }
        }
        if let Some(table) = block.loop_with_tag("_pdbx_branch_scheme.entity_id") {
            for row in 0..table.row_count() {
                let entity_id = required(table, row, "_pdbx_branch_scheme.entity_id")?;
                let asym_id = required(table, row, "_pdbx_branch_scheme.asym_id")?;
                let number = required_i32(table, row, "_pdbx_branch_scheme.num")?;
                let component_id = required(table, row, "_pdbx_branch_scheme.mon_id")?;
                let author_sequence_id = optional(table, row, "_pdbx_branch_scheme.pdb_seq_num")
                    .or_else(|| optional(table, row, "_pdbx_branch_scheme.auth_seq_num"))
                    .map(str::to_owned);
                catalog.branch_sites.push(BranchSchemeSite {
                    entity_id: entity_id.to_owned(),
                    asym_id: asym_id.to_owned(),
                    number,
                    component_id: component_id.to_owned(),
                    author_sequence_id,
                });
            }
        }
        if let Some(table) = block.loop_with_tag("_pdbx_entity_branch_link.entity_id") {
            for row in 0..table.row_count() {
                let entity_id = required(table, row, "_pdbx_entity_branch_link.entity_id")?;
                let number_1 = required_i32(
                    table,
                    row,
                    "_pdbx_entity_branch_link.entity_branch_list_num_1",
                )?;
                let atom_1 = required(table, row, "_pdbx_entity_branch_link.atom_id_1")?;
                let number_2 = required_i32(
                    table,
                    row,
                    "_pdbx_entity_branch_link.entity_branch_list_num_2",
                )?;
                let atom_2 = required(table, row, "_pdbx_entity_branch_link.atom_id_2")?;
                let order = component_bond_order(
                    optional(table, row, "_pdbx_entity_branch_link.value_order").unwrap_or("sing"),
                    table,
                    row,
                )?;
                catalog.branch_links.push(BranchLink {
                    entity_id: entity_id.to_owned(),
                    number_1,
                    atom_1: atom_1.to_owned(),
                    number_2,
                    atom_2: atom_2.to_owned(),
                    order,
                });
            }
        }
        Ok(catalog)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ResidueKey {
    asym_id: String,
    entity_id: Option<String>,
    label_sequence_id: Option<i32>,
    author_sequence_id: Option<String>,
    insertion_code: Option<String>,
    occurrence: Option<usize>,
    component_id: String,
}

#[derive(Debug, Clone)]
struct ResidueAtoms {
    key: ResidueKey,
    atoms: BTreeMap<String, AtomId>,
}

fn rebuild_model_with_connectivity(
    document: &MmcifDocument,
    source: Model,
    report: &raw::MmcifInterpretationReport,
) -> Result<(Model, usize), raw::MmcifInterpretError> {
    let block = document
        .block(report.data_block())
        .ok_or_else(|| interpret_error("interpreted mmCIF data block is unavailable"))?;
    let catalog = ConnectivityCatalog::from_block(block)?;
    let topology = source.topology();
    let mut builder = ModelBuilder::new();
    let mut pending = 0usize;

    for (instance_id, instance) in topology.instances() {
        let definition = topology
            .definition_for_instance(instance_id)
            .map_err(interpret_error)?;
        let provenance = report
            .instances()
            .iter()
            .find(|provenance| provenance.molecule() == instance_id)
            .ok_or_else(|| interpret_error("mmCIF instance provenance is incomplete"))?;
        let positions = definition
            .graph()
            .atom_ids()
            .map(|atom| {
                source
                    .position(InstanceAtomId::new(instance_id, atom))
                    .and_then(|position| position.into_unit(ANGSTROM).map_err(Into::into))
                    .map(|position| position.into_value())
                    .map_err(interpret_error)
            })
            .collect::<Result<Vec<_>, _>>()?;

        let new_instance = if let Some(molecule) = definition.macro_molecule() {
            let mut molecule = molecule.clone();
            let mut source_aromatic_bonds = BTreeSet::new();
            apply_instance_connectivity(
                molecule.graph_mut_unchecked_connectedness(),
                provenance,
                &catalog,
                &mut source_aromatic_bonds,
            )?;
            localize_source_aromatic_bonds(
                molecule.graph_mut_unchecked_connectedness(),
                &source_aromatic_bonds,
            )
            .map_err(interpret_error)?;
            let connected = molecule.graph().is_connected();
            if molecule.graph().atom_count() > 1 && !connected {
                pending += 1;
            }
            let definition = builder
                .add_macro_molecule_definition_unchecked_connectedness(&molecule)
                .map_err(interpret_error)?;
            builder
                .add_instance(
                    definition,
                    Quantity::new(positions, ANGSTROM),
                    instance.metadata().clone(),
                )
                .map_err(interpret_error)?
        } else if let Some(molecule) = definition.small_molecule() {
            let mut molecule = molecule.clone();
            let mut source_aromatic_bonds = BTreeSet::new();
            apply_instance_connectivity(
                molecule.graph_mut(),
                provenance,
                &catalog,
                &mut source_aromatic_bonds,
            )?;
            localize_source_aromatic_bonds(molecule.graph_mut(), &source_aromatic_bonds)
                .map_err(interpret_error)?;
            if molecule.graph().atom_count() > 1 && !molecule.graph().is_connected() {
                pending += 1;
            }
            let definition = builder
                .add_small_molecule_definition_unchecked_connectedness(&molecule)
                .map_err(interpret_error)?;
            builder
                .add_instance(
                    definition,
                    Quantity::new(positions, ANGSTROM),
                    instance.metadata().clone(),
                )
                .map_err(interpret_error)?
        } else {
            return Err(interpret_error(
                "mmCIF topology definition has no molecule payload",
            ));
        };

        if new_instance != instance_id {
            return Err(interpret_error(
                "mmCIF connectivity rebuild changed molecule instance order",
            ));
        }
    }

    let mut model = builder.build().map_err(interpret_error)?;
    model.set_cell(source.cell().copied());
    model
        .set_atom_data(rebound_atom_data(
            source.atom_data(),
            &model.shared_topology(),
        )?)
        .map_err(interpret_error)?;
    Ok((model, pending))
}

fn apply_instance_connectivity(
    graph: &mut Molecule,
    provenance: &raw::MmcifInstanceProvenance,
    catalog: &ConnectivityCatalog,
    source_aromatic_bonds: &mut BTreeSet<crate::core::BondId>,
) -> Result<(), raw::MmcifInterpretError> {
    let residues = residue_atoms(provenance);
    for residue in residues.values() {
        if let Some(template) = catalog
            .component_bonds
            .get(&normalized(&residue.key.component_id))
        {
            for bond in template {
                let (Some(&left), Some(&right)) = (
                    residue.atoms.get(&bond.atom_1),
                    residue.atoms.get(&bond.atom_2),
                ) else {
                    continue;
                };
                add_bond_if_missing(graph, left, right, bond.order, source_aromatic_bonds)?;
            }
        }
    }
    apply_branch_links(graph, &residues, catalog, source_aromatic_bonds)?;
    apply_polymer_links(graph, &residues, catalog, source_aromatic_bonds)
}

fn residue_atoms(provenance: &raw::MmcifInstanceProvenance) -> BTreeMap<ResidueKey, ResidueAtoms> {
    let mut residues = BTreeMap::new();
    for atom in provenance.atoms() {
        let key = ResidueKey {
            asym_id: atom.asym_id().to_owned(),
            entity_id: atom.entity_id().map(str::to_owned),
            label_sequence_id: atom.label_sequence_id(),
            author_sequence_id: atom.author_sequence_id().map(str::to_owned),
            insertion_code: atom.insertion_code().map(str::to_owned),
            occurrence: atom.occurrence(),
            component_id: atom.component_id().to_owned(),
        };
        residues
            .entry(key.clone())
            .or_insert_with(|| ResidueAtoms {
                key,
                atoms: BTreeMap::new(),
            })
            .atoms
            .insert(atom.atom_name().to_owned(), atom.atom().atom());
    }
    residues
}

fn apply_branch_links(
    graph: &mut Molecule,
    residues: &BTreeMap<ResidueKey, ResidueAtoms>,
    catalog: &ConnectivityCatalog,
    source_aromatic_bonds: &mut BTreeSet<crate::core::BondId>,
) -> Result<(), raw::MmcifInterpretError> {
    for link in &catalog.branch_links {
        let left = resolve_branch_residue(residues, catalog, &link.entity_id, link.number_1);
        let right = resolve_branch_residue(residues, catalog, &link.entity_id, link.number_2);
        let (Some(left), Some(right)) = (left, right) else {
            continue;
        };
        let (Some(&left_atom), Some(&right_atom)) =
            (left.atoms.get(&link.atom_1), right.atoms.get(&link.atom_2))
        else {
            continue;
        };
        add_bond_if_missing(
            graph,
            left_atom,
            right_atom,
            link.order,
            source_aromatic_bonds,
        )?;
    }
    Ok(())
}

fn resolve_branch_residue<'a>(
    residues: &'a BTreeMap<ResidueKey, ResidueAtoms>,
    catalog: &ConnectivityCatalog,
    entity_id: &str,
    number: i32,
) -> Option<&'a ResidueAtoms> {
    catalog
        .branch_sites
        .iter()
        .filter(|site| site.entity_id == entity_id && site.number == number)
        .find_map(|site| {
            residues.values().find(|residue| {
                residue.key.asym_id == site.asym_id
                    && residue.key.entity_id.as_deref() == Some(site.entity_id.as_str())
                    && normalized(&residue.key.component_id) == normalized(&site.component_id)
                    && site.author_sequence_id.as_deref().is_none_or(|expected| {
                        residue.key.author_sequence_id.as_deref() == Some(expected)
                    })
            })
        })
}

fn apply_polymer_links(
    graph: &mut Molecule,
    residues: &BTreeMap<ResidueKey, ResidueAtoms>,
    catalog: &ConnectivityCatalog,
    source_aromatic_bonds: &mut BTreeSet<crate::core::BondId>,
) -> Result<(), raw::MmcifInterpretError> {
    let mut chains = BTreeMap::<(String, String), Vec<&ResidueAtoms>>::new();
    for residue in residues.values() {
        let (Some(entity), Some(_)) = (
            residue.key.entity_id.as_ref(),
            residue.key.label_sequence_id,
        ) else {
            continue;
        };
        chains
            .entry((residue.key.asym_id.clone(), entity.clone()))
            .or_default()
            .push(residue);
    }
    for ((_, entity), chain) in &mut chains {
        chain.sort_by_key(|residue| residue.key.label_sequence_id);
        let Some(polymer_type) = catalog.polymer_types.get(entity) else {
            continue;
        };
        // The flag states that at least one monomer-to-monomer linkage differs
        // from the linkage implied by the polymer type, but does not identify
        // that pair. Remain conservative rather than inventing a standard bond.
        if catalog.nonstandard_polymer_linkage.contains(entity) {
            continue;
        }
        for pair in chain.windows(2) {
            let left = pair[0];
            let right = pair[1];
            let (Some(left_seq), Some(right_seq)) =
                (left.key.label_sequence_id, right.key.label_sequence_id)
            else {
                continue;
            };
            // label_seq_id includes unmodelled residues, so only adjacent values
            // prove that the two observed residues are consecutive in sequence.
            if right_seq != left_seq.saturating_add(1) {
                continue;
            }
            if is_polypeptide(polymer_type) {
                if let (Some(&carbonyl), Some(&nitrogen)) =
                    (left.atoms.get("C"), right.atoms.get("N"))
                {
                    add_bond_if_missing(
                        graph,
                        carbonyl,
                        nitrogen,
                        ComponentBondOrder::Localized(BondOrder::Single),
                        source_aromatic_bonds,
                    )?;
                }
            } else if is_nucleic_acid(polymer_type) {
                let oxygen = left.atoms.get("O3'").or_else(|| left.atoms.get("O3*"));
                if let (Some(&oxygen), Some(&phosphorus)) = (oxygen, right.atoms.get("P")) {
                    add_bond_if_missing(
                        graph,
                        oxygen,
                        phosphorus,
                        ComponentBondOrder::Localized(BondOrder::Single),
                        source_aromatic_bonds,
                    )?;
                }
            }
        }
    }
    Ok(())
}

fn add_bond_if_missing(
    graph: &mut Molecule,
    left: AtomId,
    right: AtomId,
    source_order: ComponentBondOrder,
    source_aromatic_bonds: &mut BTreeSet<crate::core::BondId>,
) -> Result<(), raw::MmcifInterpretError> {
    let existing = graph.bond_between(left, right).map_err(interpret_error)?;
    let Some(existing) = existing else {
        let bond = graph
            .add_bond(left, right, source_order.staged_order())
            .map_err(interpret_error)?;
        if source_order == ComponentBondOrder::Aromatic {
            source_aromatic_bonds.insert(bond);
        }
        return Ok(());
    };
    let existing_order = graph.bond(existing).map_err(interpret_error)?.order;
    let consistent = match source_order {
        ComponentBondOrder::Aromatic => source_aromatic_bonds.contains(&existing),
        ComponentBondOrder::Localized(order) => {
            !source_aromatic_bonds.contains(&existing) && existing_order == order
        }
    };
    if !consistent {
        let source_order = match source_order {
            ComponentBondOrder::Localized(order) => format!("{order:?}"),
            ComponentBondOrder::Aromatic => "Aromatic".to_owned(),
        };
        return Err(interpret_error(format!(
            "conflicting authoritative mmCIF bond evidence for one atom pair: {existing_order:?} versus {source_order}"
        )));
    }
    Ok(())
}

fn is_yes(value: &str) -> bool {
    matches!(value.to_ascii_lowercase().as_str(), "y" | "yes")
}

fn is_polypeptide(value: &str) -> bool {
    value.to_ascii_lowercase().contains("polypeptide")
}

fn is_nucleic_acid(value: &str) -> bool {
    let value = value.to_ascii_lowercase();
    value.contains("polyribonucleotide") || value.contains("polydeoxyribonucleotide")
}

fn component_bond_order(
    value: &str,
    table: &MmcifLoopTable,
    row: usize,
) -> Result<ComponentBondOrder, raw::MmcifInterpretError> {
    match value.to_ascii_lowercase().as_str() {
        "sing" | "poly" => Ok(ComponentBondOrder::Localized(BondOrder::Single)),
        "doub" | "pi" => Ok(ComponentBondOrder::Localized(BondOrder::Double)),
        "trip" => Ok(ComponentBondOrder::Localized(BondOrder::Triple)),
        "quad" => Ok(ComponentBondOrder::Localized(BondOrder::Quadruple)),
        "arom" | "delo" => Ok(ComponentBondOrder::Aromatic),
        other => Err(row_error(
            table,
            row,
            format!("unsupported component bond order `{other}`"),
        )),
    }
}

fn rebound_atom_data(
    source: &AtomData,
    target: &Arc<Topology>,
) -> Result<AtomData, raw::MmcifInterpretError> {
    let mut atom_data = AtomData::new(target);
    if let Some(occupancies) = source.occupancies() {
        atom_data
            .set_occupancies(occupancies)
            .map_err(interpret_error)?;
    }
    if let Some(b_factors) = source.b_factors() {
        atom_data
            .set_b_factors(b_factors)
            .map_err(interpret_error)?;
    }
    Ok(atom_data)
}

fn normalized(value: &str) -> String {
    value.to_ascii_uppercase()
}

fn required<'a>(
    table: &'a MmcifLoopTable,
    row: usize,
    tag: &str,
) -> Result<&'a str, raw::MmcifInterpretError> {
    let value = table
        .value(row, tag)
        .ok_or_else(|| row_error(table, row, format!("missing required {tag}")))?;
    value
        .optional_text()
        .ok_or_else(|| row_error(table, row, format!("missing required {tag}")))
}

fn required_i32(
    table: &MmcifLoopTable,
    row: usize,
    tag: &str,
) -> Result<i32, raw::MmcifInterpretError> {
    required(table, row, tag)?
        .parse::<i32>()
        .map_err(|_| row_error(table, row, format!("invalid integer {tag}")))
}

fn optional<'a>(table: &'a MmcifLoopTable, row: usize, tag: &str) -> Option<&'a str> {
    table.value(row, tag).and_then(MmcifValue::optional_text)
}

fn row_error(
    table: &MmcifLoopTable,
    row: usize,
    message: impl Into<String>,
) -> raw::MmcifInterpretError {
    raw::MmcifInterpretError {
        line: table
            .row(row)
            .and_then(|values| values.first())
            .map(MmcifValue::line),
        message: message.into(),
    }
}

fn interpret_error(error: impl std::fmt::Display) -> raw::MmcifInterpretError {
    raw::MmcifInterpretError {
        line: None,
        message: error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::MoleculeInstanceId;

    const PEPTIDE: &str = r#"
data_peptide
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.nstd_linkage
1 'polypeptide(L)' no
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
A 1 1 GLY
A 1 2 GLY
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
GLY N CA sing
GLY CA C sing
GLY C O doub
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY A 1 1 0.0 0.0 0.0
ATOM 2 C CA GLY A 1 1 1.4 0.0 0.0
ATOM 3 C C GLY A 1 1 2.8 0.0 0.0
ATOM 4 O O GLY A 1 1 3.8 0.0 0.0
ATOM 5 N N GLY A 1 2 4.1 0.0 0.0
ATOM 6 C CA GLY A 1 2 5.5 0.0 0.0
ATOM 7 C C GLY A 1 2 6.9 0.0 0.0
ATOM 8 O O GLY A 1 2 7.9 0.0 0.0
"#;

    const BRANCH: &str = r#"
data_branch
loop_
_entity.id
_entity.type
2 branched
loop_
_struct_asym.id
_struct_asym.entity_id
B 2
loop_
_pdbx_entity_branch_list.entity_id
_pdbx_entity_branch_list.comp_id
_pdbx_entity_branch_list.num
_pdbx_entity_branch_list.hetero
2 NAG 1 n
2 GAL 2 n
loop_
_pdbx_branch_scheme.asym_id
_pdbx_branch_scheme.entity_id
_pdbx_branch_scheme.mon_id
_pdbx_branch_scheme.num
_pdbx_branch_scheme.pdb_seq_num
B 2 NAG 1 10
B 2 GAL 2 11
loop_
_pdbx_entity_branch_link.link_id
_pdbx_entity_branch_link.entity_id
_pdbx_entity_branch_link.entity_branch_list_num_1
_pdbx_entity_branch_link.comp_id_1
_pdbx_entity_branch_link.atom_id_1
_pdbx_entity_branch_link.entity_branch_list_num_2
_pdbx_entity_branch_link.comp_id_2
_pdbx_entity_branch_link.atom_id_2
_pdbx_entity_branch_link.value_order
1 2 1 NAG O4 2 GAL C1 sing
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
NAG C4 O4 sing
GAL C1 O5 sing
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C4 NAG B 2 . 10 0.0 0.0 0.0
HETATM 2 O O4 NAG B 2 . 10 1.4 0.0 0.0
HETATM 3 C C1 GAL B 2 . 11 2.8 0.0 0.0
HETATM 4 O O5 GAL B 2 . 11 4.2 0.0 0.0
"#;

    #[test]
    fn component_templates_and_polymer_links_complete_peptide_graph() {
        let document = super::super::mmcif_document::parse_mmcif_str(
            PEPTIDE,
            super::super::mmcif_document::MmcifParseOptions::default(),
        )
        .expect("parse peptide");
        let interpretation = interpret_mmcif(&document, raw::MmcifInterpretOptions::default())
            .expect("interpret peptide");
        let graph = interpretation
            .topology()
            .graph_for_instance(MoleculeInstanceId::new(0))
            .expect("peptide instance");
        assert_eq!(graph.bond_count(), 7);
        assert_eq!(
            graph
                .bonds()
                .filter(|(_, bond)| bond.order == BondOrder::Single)
                .count(),
            5
        );
        assert_eq!(
            graph
                .bonds()
                .filter(|(_, bond)| bond.order == BondOrder::Double)
                .count(),
            2
        );
        assert!(graph.is_connected());
        assert_eq!(interpretation.report().template_bonds_pending(), 0);
    }

    #[test]
    fn sequence_gap_is_not_bridged_by_a_fake_polymer_bond() {
        let input = PEPTIDE
            .replace("A 1 2 GLY", "A 1 3 GLY")
            .replace("GLY A 1 2 4.1", "GLY A 1 3 4.1")
            .replace("GLY A 1 2 5.5", "GLY A 1 3 5.5")
            .replace("GLY A 1 2 6.9", "GLY A 1 3 6.9")
            .replace("GLY A 1 2 7.9", "GLY A 1 3 7.9");
        let document = super::super::mmcif_document::parse_mmcif_str(
            &input,
            super::super::mmcif_document::MmcifParseOptions::default(),
        )
        .expect("parse gapped peptide");
        let interpretation = interpret_mmcif(&document, raw::MmcifInterpretOptions::default())
            .expect("interpret gapped peptide");
        let graph = interpretation
            .topology()
            .graph_for_instance(MoleculeInstanceId::new(0))
            .expect("peptide instance");
        assert_eq!(graph.connected_components().len(), 2);
        assert_eq!(interpretation.report().template_bonds_pending(), 1);
    }

    #[test]
    fn branch_link_connects_branched_macromolecule() {
        let document = super::super::mmcif_document::parse_mmcif_str(
            BRANCH,
            super::super::mmcif_document::MmcifParseOptions::default(),
        )
        .expect("parse branch");
        let interpretation = interpret_mmcif(&document, raw::MmcifInterpretOptions::default())
            .expect("interpret branch");
        let graph = interpretation
            .topology()
            .graph_for_instance(MoleculeInstanceId::new(0))
            .expect("branch instance");
        assert_eq!(graph.bond_count(), 3);
        assert!(graph
            .bonds()
            .all(|(_, bond)| bond.order == BondOrder::Single));
        assert!(graph.is_connected());
        assert_eq!(interpretation.report().template_bonds_pending(), 0);
    }
}
