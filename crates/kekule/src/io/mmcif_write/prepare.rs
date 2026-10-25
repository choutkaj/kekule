//! Validate represented chemistry and prepare ordered atom/connection rows.
use super::super::element_symbol::structure_element_symbol;
use super::entities::{AssignedEntityKind, AsymRow, AtomEntityAssignment, EntityPlan, EntityRow};
use super::implied::mark_implied;
use super::one_based_serial;
use super::{fail, invalid_model, MmcifAtomField, MmcifSampleField, MmcifWriteIssue};
use crate::diagnostics::FormatError;
use crate::geometry::Point3;
use crate::io::UnrepresentablePolicy;
use crate::molecule::{Atom, AtomId, BondOrder, ImplicitHydrogens};
use crate::structure::ModelView;
use crate::topology::{
    Hierarchy, InstanceAtomId, InstanceBondId, MoleculeDefinition, MoleculeInstance, ResidueId,
};
use crate::units::ANGSTROM;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub(super) struct AtomRow {
    pub(super) atom: InstanceAtomId,
    pub(super) residue: Option<ResidueId>,
    pub(super) entity_id: String,
    pub(super) asym_id: String,
    pub(super) group_pdb: String,
    pub(super) type_symbol: String,
    pub(super) label_atom_id: String,
    pub(super) label_alt_id: Option<String>,
    pub(super) label_comp_id: String,
    pub(super) label_seq_id: Option<i32>,
    pub(super) insertion_code: Option<String>,
    pub(super) position: Point3,
    pub(super) occupancy: Option<f64>,
    pub(super) b_factor: Option<f64>,
    pub(super) formal_charge: i8,
    pub(super) auth_seq_id: Option<String>,
    pub(super) auth_comp_id: String,
    pub(super) auth_asym_id: String,
    pub(super) auth_atom_id: String,
}

/// One bond of the model, written as a `_struct_conn` row unless the reader
/// derives it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct ConnectionRow {
    pub(super) left: InstanceAtomId,
    pub(super) right: InstanceAtomId,
    pub(super) order: BondOrder,
    /// How the reader derives the bond, or `None` when the block states it.
    pub(super) implied: Option<Implied>,
}

/// Evidence from which an mmCIF reader derives a bond the block leaves out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Implied {
    /// The bundled component template of the bond's residue, which rebuilds
    /// all of that residue's bonds.
    Template,
    /// Polymer linkage between residues adjacent in sequence, within bonding
    /// distance in the written coordinates.
    PolymerLink,
}

#[derive(Debug)]
pub(super) struct PreparedModel {
    pub(super) entities: Vec<EntityRow>,
    pub(super) asyms: Vec<AsymRow>,
    pub(super) atoms: Vec<AtomRow>,
    pub(super) connections: Vec<ConnectionRow>,
    /// Atoms whose fixed implicit-hydrogen counts the rows leave out.
    pub(super) omitted_hydrogens: usize,
}

pub(super) fn prepare_model(
    model: ModelView<'_>,
    plan: EntityPlan,
    implicit_hydrogens: UnrepresentablePolicy,
    coordinate_precision: usize,
) -> Result<PreparedModel, FormatError> {
    let EntityPlan {
        entities,
        asyms,
        atoms: assignments,
    } = plan;

    let mut atoms = Vec::new();
    let mut omitted_hydrogens = 0;
    let hierarchy = model.topology().hierarchy();
    for (id, molecule) in model.topology().instance_records() {
        let definition = model
            .topology()
            .instance(id)
            .expect("listed instances are topology molecules")
            .definition();
        omitted_hydrogens += validate_graph_chemistry(molecule, definition, implicit_hydrogens)?;
        if definition.molecule().atom_ids().any(|atom| {
            hierarchy
                .atom_site_for_atom(molecule.qualify_atom(atom))
                .is_some()
        }) {
            collect_macro_rows(
                model,
                molecule,
                definition,
                hierarchy,
                &assignments,
                &mut atoms,
            )?;
        } else {
            collect_small_rows(model, molecule, definition, &assignments, &mut atoms)?;
        }
    }

    // Rows follow dense atom order, so the interpreted block keeps this
    // topology's order and coordinate files written in it stay aligned.
    atoms.sort_by_key(|row| {
        model
            .topology()
            .atom_index(row.atom)
            .expect("prepared atom-site row references a topology atom")
    });
    number_unsequenced_polymers(&mut atoms, &entities, &asyms, hierarchy);

    let atom_indexes = atoms
        .iter()
        .enumerate()
        .map(|(index, row)| (row.atom, index))
        .collect::<BTreeMap<_, _>>();
    validate_atom_identities(&atoms)?;
    let mut connections = Vec::new();
    for bond in model.topology().bonds() {
        let order = supported_bond_order(bond.id(), bond.order)?;
        let [left, right] = bond.atoms().map(|atom| atom.id());
        connections.push(ConnectionRow {
            left,
            right,
            order,
            implied: None,
        });
    }
    let elements = atoms
        .iter()
        .map(|row| {
            model
                .topology()
                .atom(row.atom)
                .map(|atom| atom.atom().element)
                .ok_or_else(|| invalid_model(format!("invalid topology atom: {}", row.atom)))
        })
        .collect::<Result<Vec<_>, _>>()?;
    mark_implied(
        &atoms,
        &elements,
        &mut connections,
        &atom_indexes,
        coordinate_precision,
    );
    // A polymer link another ensemble sample stretches beyond bonding
    // distance is stated after all, so its partners are checked too.
    let selectors = ConnectionSelectors::new(&atoms);
    for connection in connections
        .iter()
        .filter(|connection| connection.implied != Some(Implied::Template))
    {
        selectors.validate(connection.left, &atoms, &atom_indexes)?;
        selectors.validate(connection.right, &atoms, &atom_indexes)?;
    }
    Ok(PreparedModel {
        entities,
        asyms,
        atoms,
        connections,
        omitted_hydrogens,
    })
}

/// Fails on chemistry mmCIF cannot carry; returns how many atoms have a
/// fixed implicit-hydrogen count that `implicit_hydrogens` omits.
fn validate_graph_chemistry(
    molecule: &MoleculeInstance,
    definition: &MoleculeDefinition,
    implicit_hydrogens: UnrepresentablePolicy,
) -> Result<usize, FormatError> {
    let mut omitted_hydrogens = 0;
    let molecule_definition = definition.molecule();
    if molecule_definition.stereo_elements().next().is_some()
        || molecule_definition.stereo_groups().next().is_some()
    {
        return Err(fail(MmcifWriteIssue::UnsupportedStereo(molecule.id())));
    }
    for (atom_id, atom) in molecule_definition.atoms() {
        let atom_id = molecule.qualify_atom(atom_id);
        if structure_element_symbol(atom.element, atom.isotope.map(std::num::NonZero::get))
            .is_none()
        {
            return Err(fail(MmcifWriteIssue::UnsupportedAtomField {
                atom: atom_id,
                field: MmcifAtomField::Isotope,
            }));
        }
        if atom.radical.is_some() {
            return Err(fail(MmcifWriteIssue::UnsupportedAtomField {
                atom: atom_id,
                field: MmcifAtomField::Radical,
            }));
        }
        if atom.hydrogens != ImplicitHydrogens::default() {
            if implicit_hydrogens == UnrepresentablePolicy::Fail {
                return Err(fail(MmcifWriteIssue::UnsupportedAtomField {
                    atom: atom_id,
                    field: MmcifAtomField::Hydrogens,
                }));
            }
            omitted_hydrogens += 1;
        }
        if atom.atom_map.is_some() {
            return Err(fail(MmcifWriteIssue::UnsupportedAtomField {
                atom: atom_id,
                field: MmcifAtomField::AtomMap,
            }));
        }
        if !(-8..=8).contains(&atom.formal_charge) {
            return Err(fail(MmcifWriteIssue::FormalChargeOutOfRange {
                atom: atom_id,
                charge: atom.formal_charge,
            }));
        }
    }
    Ok(omitted_hydrogens)
}

fn collect_macro_rows(
    model: ModelView<'_>,
    molecule: &MoleculeInstance,
    definition: &MoleculeDefinition,
    hierarchy: &Hierarchy,
    assignments: &BTreeMap<InstanceAtomId, AtomEntityAssignment>,
    rows: &mut Vec<AtomRow>,
) -> Result<(), FormatError> {
    for (atom_id, atom) in definition.molecule().atoms() {
        let qualified = molecule.qualify_atom(atom_id);
        let site = hierarchy
            .atom_site_for_atom(qualified)
            .ok_or_else(|| fail(MmcifWriteIssue::MissingAtomSite(qualified)))?;
        let residue = hierarchy
            .require_residue(site.residue)
            .map_err(|error| invalid_hierarchy("invalid topology hierarchy").with_source(error))?;
        let chain = hierarchy
            .require_chain(residue.chain)
            .map_err(|error| invalid_hierarchy("invalid topology hierarchy").with_source(error))?;
        let assignment = assignments
            .get(&qualified)
            .ok_or_else(|| fail(MmcifWriteIssue::MissingAtomSite(qualified)))?;
        if site
            .metadata
            .label_asym_id
            .as_deref()
            .is_some_and(|value| value != chain.label_id)
        {
            return Err(fail(MmcifWriteIssue::InconsistentAtomSite {
                atom: qualified,
                field: MmcifAtomField::LabelAsymId,
            }));
        }
        if site
            .metadata
            .type_symbol
            .as_deref()
            .is_some_and(|value| !value.eq_ignore_ascii_case(type_symbol(atom)))
        {
            return Err(fail(MmcifWriteIssue::InconsistentAtomSite {
                atom: qualified,
                field: MmcifAtomField::TypeSymbol,
            }));
        }
        let group_pdb = assignment.kind.default_group_pdb().to_owned();
        let label_atom_id = site
            .metadata
            .label_atom_id
            .as_ref()
            .or(site.metadata.auth_atom_id.as_ref())
            .cloned()
            .unwrap_or_else(|| generated_atom_name(atom.element.symbol(), atom_id));
        let label_comp_id = residue
            .label_comp_id
            .clone()
            .unwrap_or_else(|| residue.name.clone());
        rows.push(AtomRow {
            atom: qualified,
            residue: Some(site.residue()),
            entity_id: assignment.entity_id.clone(),
            asym_id: assignment.asym_id.clone(),
            group_pdb,
            type_symbol: type_symbol(atom).to_owned(),
            label_atom_id: label_atom_id.clone(),
            label_alt_id: None,
            label_comp_id: label_comp_id.clone(),
            label_seq_id: residue.label_seq_id,
            insertion_code: residue.insertion_code.clone(),
            position: model
                .position(qualified)
                .ok_or_else(|| invalid_model("inconsistent source model"))?
                .value_in(ANGSTROM)
                .map_err(|error| invalid_model("inconsistent source model").with_source(error))?,
            occupancy: model.occupancy(dense_index(model, qualified)?),
            b_factor: model
                .b_factor(dense_index(model, qualified)?)
                .map(|value| value.value_in(crate::units::SQUARE_ANGSTROM))
                .transpose()
                .map_err(|error| invalid_model("inconsistent source model").with_source(error))?,
            formal_charge: atom.formal_charge,
            auth_seq_id: residue
                .author_seq_id
                .clone()
                .or_else(|| residue.label_seq_id.map(|value| value.to_string())),
            auth_comp_id: residue.author_comp_id.clone().unwrap_or(label_comp_id),
            auth_asym_id: site
                .metadata
                .auth_asym_id
                .clone()
                .or_else(|| chain.author_id.clone())
                .unwrap_or_else(|| chain.label_id.clone()),
            auth_atom_id: site.metadata.auth_atom_id.clone().unwrap_or(label_atom_id),
        });
    }
    Ok(())
}

fn collect_small_rows(
    model: ModelView<'_>,
    molecule: &MoleculeInstance,
    definition: &MoleculeDefinition,
    assignments: &BTreeMap<InstanceAtomId, AtomEntityAssignment>,
    rows: &mut Vec<AtomRow>,
) -> Result<(), FormatError> {
    for (atom_id, atom) in definition.molecule().atoms() {
        let qualified = molecule.qualify_atom(atom_id);
        let assignment = assignments
            .get(&qualified)
            .ok_or_else(|| fail(MmcifWriteIssue::MissingAtomSite(qualified)))?;
        let component_id = if assignment.kind == AssignedEntityKind::Water {
            "HOH"
        } else {
            "MOL"
        };
        let atom_name = if assignment.kind == AssignedEntityKind::Water {
            water_atom_name(definition.molecule(), atom_id)
        } else {
            generated_atom_name(atom.element.symbol(), atom_id)
        };
        rows.push(AtomRow {
            atom: qualified,
            residue: None,
            entity_id: assignment.entity_id.clone(),
            asym_id: assignment.asym_id.clone(),
            group_pdb: assignment.kind.default_group_pdb().to_owned(),
            type_symbol: type_symbol(atom).to_owned(),
            label_atom_id: atom_name.clone(),
            label_alt_id: None,
            label_comp_id: component_id.to_owned(),
            label_seq_id: None,
            insertion_code: None,
            position: model
                .position(qualified)
                .ok_or_else(|| invalid_model("inconsistent source model"))?
                .value_in(ANGSTROM)
                .map_err(|error| invalid_model("inconsistent source model").with_source(error))?,
            occupancy: model.occupancy(dense_index(model, qualified)?),
            b_factor: model
                .b_factor(dense_index(model, qualified)?)
                .map(|value| value.value_in(crate::units::SQUARE_ANGSTROM))
                .transpose()
                .map_err(|error| invalid_model("inconsistent source model").with_source(error))?,
            formal_charge: atom.formal_charge,
            auth_seq_id: None,
            auth_comp_id: component_id.to_owned(),
            auth_asym_id: assignment.asym_id.clone(),
            auth_atom_id: atom_name,
        });
    }
    Ok(())
}

/// The Chemical Component Dictionary names of water atoms: `O`, `H1`, `H2`.
fn water_atom_name(molecule: &crate::molecule::Molecule, atom: AtomId) -> String {
    let hydrogen = molecule
        .atoms()
        .filter(|(_, candidate)| candidate.element.atomic_number() == 1)
        .position(|(id, _)| id == atom);
    match hydrogen {
        Some(index) => format!("H{}", index + 1),
        None => "O".to_owned(),
    }
}

/// The atom's `type_symbol`; validation has already rejected isotopes the
/// format cannot name.
fn type_symbol(atom: &Atom) -> &'static str {
    structure_element_symbol(atom.element, atom.isotope.map(std::num::NonZero::get))
        .expect("validated atoms have a structure element symbol")
}

fn generated_atom_name(symbol: &str, atom: AtomId) -> String {
    format!("{symbol}{}", one_based_serial(atom.raw()))
}

fn supported_bond_order(bond: InstanceBondId, order: BondOrder) -> Result<BondOrder, FormatError> {
    match order {
        BondOrder::Single
        | BondOrder::Double
        | BondOrder::Triple
        | BondOrder::Quadruple
        | BondOrder::Dative => Ok(order),
        BondOrder::Zero => Err(fail(MmcifWriteIssue::UnsupportedBondOrder { bond, order })),
    }
}

/// Rows by label asym and atom name, among which a `_struct_conn` partner
/// selector must match exactly one: the one with its component, and its
/// label sequence number, author sequence number, and insertion code where
/// it states them, as the reader matches partners.
struct ConnectionSelectors<'a> {
    by_name: BTreeMap<(&'a str, &'a str), Vec<&'a AtomRow>>,
}

impl<'a> ConnectionSelectors<'a> {
    fn new(rows: &'a [AtomRow]) -> Self {
        let mut by_name = BTreeMap::<_, Vec<_>>::new();
        for row in rows {
            by_name
                .entry((row.asym_id.as_str(), row.label_atom_id.as_str()))
                .or_default()
                .push(row);
        }
        Self { by_name }
    }

    fn validate(
        &self,
        atom: InstanceAtomId,
        rows: &[AtomRow],
        indexes: &BTreeMap<InstanceAtomId, usize>,
    ) -> Result<(), FormatError> {
        let row = &rows[*indexes
            .get(&atom)
            .ok_or_else(|| invalid_model(format!("missing atom row for {atom}")))?];
        let stated = |field: &Option<String>, candidate: &Option<String>| {
            field.is_none() || field == candidate
        };
        let matches = self
            .by_name
            .get(&(row.asym_id.as_str(), row.label_atom_id.as_str()))
            .map_or(0, |candidates| {
                candidates
                    .iter()
                    .filter(|candidate| {
                        candidate.label_comp_id == row.label_comp_id
                            && row
                                .label_seq_id
                                .is_none_or(|sequence| candidate.label_seq_id == Some(sequence))
                            && stated(&row.auth_seq_id, &candidate.auth_seq_id)
                            && stated(&row.insertion_code, &candidate.insertion_code)
                    })
                    .count()
            });
        if matches != 1 {
            return Err(fail(MmcifWriteIssue::AmbiguousConnectionSelector(atom)));
        }
        Ok(())
    }
}

/// mmCIF sequences every polymer residue by `label_seq_id`. Sources that
/// number residues by author alone, as PDB does, leave a polymer asym without
/// one; its residues are numbered from 1 in chain order. Unobserved residues
/// are unknown, so the numbers count observed residues only.
fn number_unsequenced_polymers(
    rows: &mut [AtomRow],
    entities: &[EntityRow],
    asyms: &[AsymRow],
    hierarchy: &Hierarchy,
) {
    let polymers = entities
        .iter()
        .filter(|entity| entity.kind == AssignedEntityKind::Polymer)
        .map(|entity| entity.id.as_str())
        .collect::<BTreeSet<_>>();
    let chain_position = |residue: ResidueId| {
        let chain = hierarchy
            .residue(residue)
            .expect("atom-site rows reference live residues")
            .chain();
        hierarchy
            .chain(chain)
            .expect("published residues reference live chains")
            .residues()
            .iter()
            .position(|candidate| *candidate == residue)
            .expect("a chain lists its residues")
    };
    for asym in asyms
        .iter()
        .filter(|asym| polymers.contains(asym.entity_id.as_str()))
    {
        let in_asym = |row: &&mut AtomRow| row.asym_id == asym.id;
        if rows
            .iter_mut()
            .filter(in_asym)
            .any(|row| row.label_seq_id.is_some())
        {
            continue;
        }
        let mut residues = rows
            .iter_mut()
            .filter(in_asym)
            .filter_map(|row| row.residue)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        residues.sort_by_key(|residue| chain_position(*residue));
        let numbers = residues.into_iter().zip(1..).collect::<BTreeMap<_, i32>>();
        for row in rows.iter_mut().filter(in_asym) {
            row.label_seq_id = row.residue.map(|residue| numbers[&residue]);
        }
    }
}

fn validate_atom_identities(rows: &[AtomRow]) -> Result<(), FormatError> {
    let mut identities = BTreeSet::new();
    for row in rows {
        let residue = if let Some(sequence) = row.label_seq_id {
            format!(
                "label:{sequence}:{}",
                row.insertion_code.as_deref().unwrap_or("")
            )
        } else if let Some(sequence) = &row.auth_seq_id {
            format!(
                "auth:{sequence}:{}",
                row.insertion_code.as_deref().unwrap_or("")
            )
        } else {
            format!("unsequenced:{:?}", row.residue)
        };
        if !identities.insert((row.asym_id.clone(), residue, row.label_atom_id.clone())) {
            return Err(fail(MmcifWriteIssue::DuplicateAtomIdentity(row.atom)));
        }
    }
    Ok(())
}

fn invalid_hierarchy(message: impl Into<String>) -> FormatError {
    invalid_model(message)
}

pub(super) fn validate_ensemble_sample(
    model: &PreparedModel,
    candidate: &PreparedModel,
    sample: usize,
) -> Result<(), FormatError> {
    if candidate.entities != model.entities {
        return Err(fail(MmcifWriteIssue::IncompatibleEnsembleSample {
            sample,
            field: MmcifSampleField::Entities,
        }));
    }
    if candidate.asyms != model.asyms {
        return Err(fail(MmcifWriteIssue::IncompatibleEnsembleSample {
            sample,
            field: MmcifSampleField::StructuralAsymmetry,
        }));
    }
    let bonds = |model: &PreparedModel| {
        model
            .connections
            .iter()
            .map(|connection| (connection.left, connection.right, connection.order))
            .collect::<Vec<_>>()
    };
    if bonds(candidate) != bonds(model) {
        return Err(fail(MmcifWriteIssue::IncompatibleEnsembleSample {
            sample,
            field: MmcifSampleField::Connectivity,
        }));
    }
    if candidate.atoms.len() != model.atoms.len()
        || candidate
            .atoms
            .iter()
            .zip(&model.atoms)
            .any(|(left, right)| !same_static_atom_row(left, right))
    {
        return Err(fail(MmcifWriteIssue::IncompatibleEnsembleSample {
            sample,
            field: MmcifSampleField::AtomSites,
        }));
    }
    Ok(())
}

fn same_static_atom_row(left: &AtomRow, right: &AtomRow) -> bool {
    left.atom == right.atom
        && left.residue == right.residue
        && left.entity_id == right.entity_id
        && left.asym_id == right.asym_id
        && left.group_pdb == right.group_pdb
        && left.type_symbol == right.type_symbol
        && left.label_atom_id == right.label_atom_id
        && left.label_alt_id == right.label_alt_id
        && left.label_comp_id == right.label_comp_id
        && left.label_seq_id == right.label_seq_id
        && left.insertion_code == right.insertion_code
        && left.formal_charge == right.formal_charge
        && left.auth_seq_id == right.auth_seq_id
        && left.auth_comp_id == right.auth_comp_id
        && left.auth_asym_id == right.auth_asym_id
        && left.auth_atom_id == right.auth_atom_id
}

fn dense_index(
    model: ModelView<'_>,
    atom: InstanceAtomId,
) -> Result<crate::topology::TopologyAtomIndex, FormatError> {
    model
        .topology()
        .atom_index(atom)
        .ok_or_else(|| invalid_model(format!("invalid topology atom: {atom}")))
}

#[cfg(test)]
mod capacity_tests {
    use super::one_based_serial;

    #[test]
    fn one_based_serials_widen_before_incrementing() {
        assert_eq!(one_based_serial(0), 1);
        assert_eq!(one_based_serial(u32::MAX), u64::from(u32::MAX) + 1);
    }
}
