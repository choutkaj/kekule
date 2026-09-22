//! Validate represented chemistry and prepare ordered atom/connection rows.
use super::entities::{AsymRow, AtomEntityAssignment, EntityKind, EntityPlan, EntityRow};
use super::one_based_serial;
use super::MmcifWriteError;
use crate::core::{AtomId, BondOrder, HydrogenDeclaration};
use crate::geometry::Point3;
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct ConnectionRow {
    pub(super) left: InstanceAtomId,
    pub(super) right: InstanceAtomId,
    pub(super) order: BondOrder,
}

#[derive(Debug)]
pub(super) struct PreparedModel {
    pub(super) entities: Vec<EntityRow>,
    pub(super) asyms: Vec<AsymRow>,
    pub(super) atoms: Vec<AtomRow>,
    pub(super) connections: Vec<ConnectionRow>,
}

pub(super) fn prepare_model(
    model: ModelView<'_>,
    plan: EntityPlan,
) -> Result<PreparedModel, MmcifWriteError> {
    let EntityPlan {
        entities,
        asyms,
        atoms: assignments,
    } = plan;

    let mut atoms = Vec::new();
    let hierarchy = model.topology().hierarchy();
    for (id, molecule) in model.topology().instances() {
        let definition = model
            .topology()
            .definition_for_instance(id)
            .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?;
        validate_graph_chemistry(molecule, definition)?;
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

    let hierarchy_order = model
        .topology()
        .hierarchy()
        .chains()
        .flat_map(|(_, chain)| chain.residues().iter().copied())
        .map(|residue| {
            model
                .topology()
                .hierarchy()
                .residue(residue)
                .expect("published hierarchy chain references a live residue")
        })
        .flat_map(|residue| residue.atom_sites().iter().copied())
        .map(|site| {
            model
                .topology()
                .hierarchy()
                .atom_site(site)
                .expect("published hierarchy residue references a live atom site")
        })
        .enumerate()
        .map(|(index, site)| (site.atom(), index))
        .collect::<BTreeMap<_, _>>();
    atoms.sort_by_key(|row| {
        hierarchy_order
            .get(&row.atom)
            .copied()
            .unwrap_or(usize::MAX)
    });

    let atom_indexes = atoms
        .iter()
        .enumerate()
        .map(|(index, row)| (row.atom, index))
        .collect::<BTreeMap<_, _>>();
    validate_atom_identities(&atoms)?;
    let mut connections = Vec::new();
    for (bond_id, bond) in model.topology().bonds() {
        let order = supported_bond_order(bond_id, bond.order)?;
        let molecule = model
            .topology()
            .instance(bond_id.molecule())
            .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?;
        let left = molecule.qualify_atom(bond.a());
        let right = molecule.qualify_atom(bond.b());
        validate_connection_selector(left, &atoms, &atom_indexes)?;
        validate_connection_selector(right, &atoms, &atom_indexes)?;
        connections.push(ConnectionRow { left, right, order });
    }
    Ok(PreparedModel {
        entities,
        asyms,
        atoms,
        connections,
    })
}

fn validate_graph_chemistry(
    molecule: &MoleculeInstance,
    definition: &MoleculeDefinition,
) -> Result<(), MmcifWriteError> {
    let molecule_definition = definition.molecule();
    if molecule_definition.stereo_elements().next().is_some()
        || molecule_definition.stereo_groups().next().is_some()
    {
        return Err(MmcifWriteError::UnsupportedStereo(molecule.id()));
    }
    for (atom_id, atom) in molecule_definition.atoms() {
        let atom_id = molecule.qualify_atom(atom_id);
        if atom.isotope.is_some() {
            return Err(MmcifWriteError::UnsupportedAtomField {
                atom: atom_id,
                field: "isotope",
            });
        }
        if atom.radical.is_some() {
            return Err(MmcifWriteError::UnsupportedAtomField {
                atom: atom_id,
                field: "radical",
            });
        }
        if atom.hydrogens != HydrogenDeclaration::default() {
            return Err(MmcifWriteError::UnsupportedAtomField {
                atom: atom_id,
                field: "hydrogens",
            });
        }
        if atom.atom_map.is_some() {
            return Err(MmcifWriteError::UnsupportedAtomField {
                atom: atom_id,
                field: "atom_map",
            });
        }
        if !(-8..=8).contains(&atom.formal_charge) {
            return Err(MmcifWriteError::FormalChargeOutOfRange {
                atom: atom_id,
                charge: atom.formal_charge,
            });
        }
    }
    Ok(())
}

fn collect_macro_rows(
    model: ModelView<'_>,
    molecule: &MoleculeInstance,
    definition: &MoleculeDefinition,
    hierarchy: &Hierarchy,
    assignments: &BTreeMap<InstanceAtomId, AtomEntityAssignment>,
    rows: &mut Vec<AtomRow>,
) -> Result<(), MmcifWriteError> {
    for (atom_id, atom) in definition.molecule().atoms() {
        let qualified = molecule.qualify_atom(atom_id);
        let site = hierarchy
            .atom_site_for_atom(qualified)
            .ok_or(MmcifWriteError::MissingAtomSite(qualified))?;
        let residue = hierarchy
            .residue(site.residue)
            .map_err(|error| invalid_hierarchy(error.to_string()))?;
        let chain = hierarchy
            .chain(residue.chain)
            .map_err(|error| invalid_hierarchy(error.to_string()))?;
        let assignment = assignments
            .get(&qualified)
            .ok_or(MmcifWriteError::MissingAtomSite(qualified))?;
        if assignment.asym_id != chain.label_id {
            return Err(MmcifWriteError::InconsistentAtomSite {
                atom: qualified,
                field: "entity/asymmetry assignment",
            });
        }
        if site
            .metadata
            .label_asym_id
            .as_deref()
            .is_some_and(|value| value != chain.label_id)
        {
            return Err(MmcifWriteError::InconsistentAtomSite {
                atom: qualified,
                field: "label_asym_id",
            });
        }
        if site
            .metadata
            .type_symbol
            .as_deref()
            .is_some_and(|value| !value.eq_ignore_ascii_case(atom.element.symbol()))
        {
            return Err(MmcifWriteError::InconsistentAtomSite {
                atom: qualified,
                field: "type_symbol",
            });
        }
        let group_pdb = normalized_group_pdb(qualified, None, assignment.kind.default_group_pdb())?;
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
            asym_id: chain.label_id.clone(),
            group_pdb,
            type_symbol: atom.element.symbol().to_owned(),
            label_atom_id: label_atom_id.clone(),
            label_alt_id: None,
            label_comp_id: label_comp_id.clone(),
            label_seq_id: residue.label_seq_id,
            insertion_code: residue.insertion_code.clone(),
            position: model
                .position(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?
                .value_in(ANGSTROM)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            occupancy: model
                .occupancy(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            b_factor: model
                .b_factor(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?
                .map(|value| value.value_in(crate::units::SQUARE_ANGSTROM))
                .transpose()
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
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
) -> Result<(), MmcifWriteError> {
    for (atom_id, atom) in definition.molecule().atoms() {
        let qualified = molecule.qualify_atom(atom_id);
        let assignment = assignments
            .get(&qualified)
            .ok_or(MmcifWriteError::MissingAtomSite(qualified))?;
        let component_id = if assignment.kind == EntityKind::Water {
            "HOH"
        } else {
            "MOL"
        };
        let atom_name = generated_atom_name(atom.element.symbol(), atom_id);
        rows.push(AtomRow {
            atom: qualified,
            residue: None,
            entity_id: assignment.entity_id.clone(),
            asym_id: assignment.asym_id.clone(),
            group_pdb: normalized_group_pdb(qualified, None, assignment.kind.default_group_pdb())?,
            type_symbol: atom.element.symbol().to_owned(),
            label_atom_id: atom_name.clone(),
            label_alt_id: None,
            label_comp_id: component_id.to_owned(),
            label_seq_id: None,
            insertion_code: None,
            position: model
                .position(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?
                .value_in(ANGSTROM)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            occupancy: model
                .occupancy(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            b_factor: model
                .b_factor(qualified)
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?
                .map(|value| value.value_in(crate::units::SQUARE_ANGSTROM))
                .transpose()
                .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?,
            formal_charge: atom.formal_charge,
            auth_seq_id: None,
            auth_comp_id: component_id.to_owned(),
            auth_asym_id: assignment.asym_id.clone(),
            auth_atom_id: atom_name,
        });
    }
    Ok(())
}

fn normalized_group_pdb(
    atom: InstanceAtomId,
    value: Option<&str>,
    default: &str,
) -> Result<String, MmcifWriteError> {
    let value = value.unwrap_or(default);
    if value.eq_ignore_ascii_case("ATOM") {
        Ok("ATOM".to_owned())
    } else if value.eq_ignore_ascii_case("HETATM") {
        Ok("HETATM".to_owned())
    } else {
        Err(MmcifWriteError::InvalidGroupPdb {
            atom,
            value: value.to_owned(),
        })
    }
}

fn generated_atom_name(symbol: &str, atom: AtomId) -> String {
    format!("{symbol}{}", one_based_serial(atom.raw()))
}

fn supported_bond_order(
    bond: InstanceBondId,
    order: BondOrder,
) -> Result<BondOrder, MmcifWriteError> {
    match order {
        BondOrder::Single | BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple => {
            Ok(order)
        }
        BondOrder::Zero | BondOrder::Dative => {
            Err(MmcifWriteError::UnsupportedBondOrder { bond, order })
        }
    }
}

fn validate_connection_selector(
    atom: InstanceAtomId,
    rows: &[AtomRow],
    indexes: &BTreeMap<InstanceAtomId, usize>,
) -> Result<(), MmcifWriteError> {
    let row =
        rows.get(*indexes.get(&atom).ok_or_else(|| {
            MmcifWriteError::InvalidModel(format!("missing atom row for {atom}"))
        })?)
        .expect("atom row index is valid");
    let matches = rows
        .iter()
        .filter(|candidate| {
            candidate.asym_id == row.asym_id
                && candidate.label_atom_id == row.label_atom_id
                && row
                    .label_seq_id
                    .is_none_or(|sequence| candidate.label_seq_id == Some(sequence))
        })
        .count();
    if matches != 1 {
        return Err(MmcifWriteError::AmbiguousConnectionSelector(atom));
    }
    Ok(())
}

fn validate_atom_identities(rows: &[AtomRow]) -> Result<(), MmcifWriteError> {
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
            return Err(MmcifWriteError::DuplicateAtomIdentity(row.atom));
        }
    }
    Ok(())
}

fn invalid_hierarchy(message: impl Into<String>) -> MmcifWriteError {
    MmcifWriteError::InvalidHierarchy {
        message: message.into(),
    }
}

pub(super) fn validate_ensemble_member(
    model: &PreparedModel,
    candidate: &PreparedModel,
    member: usize,
) -> Result<(), MmcifWriteError> {
    if candidate.entities != model.entities {
        return Err(MmcifWriteError::IncompatibleEnsembleMember {
            member,
            field: "entity representation",
        });
    }
    if candidate.asyms != model.asyms {
        return Err(MmcifWriteError::IncompatibleEnsembleMember {
            member,
            field: "structural asymmetry representation",
        });
    }
    if candidate.connections != model.connections {
        return Err(MmcifWriteError::IncompatibleEnsembleMember {
            member,
            field: "connectivity representation",
        });
    }
    if candidate.atoms.len() != model.atoms.len()
        || candidate
            .atoms
            .iter()
            .zip(&model.atoms)
            .any(|(left, right)| !same_static_atom_row(left, right))
    {
        return Err(MmcifWriteError::IncompatibleEnsembleMember {
            member,
            field: "atom-site representation",
        });
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

#[cfg(test)]
mod capacity_tests {
    use super::one_based_serial;

    #[test]
    fn one_based_serials_widen_before_incrementing() {
        assert_eq!(one_based_serial(0), 1);
        assert_eq!(one_based_serial(u32::MAX), u64::from(u32::MAX) + 1);
    }
}
