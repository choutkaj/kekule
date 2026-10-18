//! Exact reconstruction of stored topology state.
//!
//! Native deserialization restores a published topology from its parts. The
//! parts are untrusted, so publication goes through [`TopologyBuilder`] and
//! the hierarchy storage validator; only class values are restored as stored
//! instead of being inferred again, so inference rules can change without
//! changing a stored topology.

use crate::molecule::Molecule;
use crate::properties::TopologyProperties;

use super::{
    Hierarchy, HierarchyBuilder, InstanceAtomId, MoleculeClass, MoleculeDefinitionId, Topology,
    TopologyBuildError, TopologyBuilder,
};

/// One stored definition: its molecule, published class, and whether the
/// class was assigned explicitly.
#[derive(Debug, Clone)]
pub(crate) struct RawDefinition {
    pub(crate) molecule: Molecule,
    pub(crate) class: MoleculeClass,
    pub(crate) explicit: bool,
}

/// The stored parts of one published topology.
#[derive(Debug, Clone)]
pub(crate) struct TopologyParts {
    pub(crate) definitions: Vec<RawDefinition>,
    pub(crate) instances: Vec<MoleculeDefinitionId>,
    /// `None` keeps the builder's instance order.
    pub(crate) atom_order: Option<Vec<InstanceAtomId>>,
    /// Residue classes, explicit or inferred, are stored on its residues.
    pub(crate) hierarchy: Hierarchy,
    pub(crate) properties: TopologyProperties,
}

impl Topology {
    /// Publishes stored parts as a new topology layout with the same semantic
    /// IDs, dense order, hierarchy, classes, and properties.
    pub(crate) fn from_parts(parts: TopologyParts) -> Result<Self, TopologyBuildError> {
        let TopologyParts {
            definitions,
            instances,
            atom_order,
            hierarchy,
            properties,
        } = parts;
        let mut builder = TopologyBuilder::new();
        builder.reserve_definitions(definitions.len())?;
        builder.reserve_instances(instances.len())?;
        for definition in definitions {
            let id = builder.add_molecule_definition(definition.molecule)?;
            builder.set_restored_molecule_class(id, definition.class, definition.explicit)?;
        }
        for definition in instances {
            builder.add_instance(definition)?;
        }
        if let Some(order) = atom_order {
            builder.set_atom_order(order)?;
        }
        // `Hierarchy::from_raw_parts` already guarantees consistent storage;
        // publication checks that every site names a staged atom.
        *builder.hierarchy_mut() = HierarchyBuilder::from_storage(hierarchy);
        builder.install_properties(properties);
        builder.build_restored()
    }
}
