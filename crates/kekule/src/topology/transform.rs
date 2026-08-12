//! Immutable whole-instance topology transformations.
//!
//! These operations preserve complete molecule definitions and instances. They
//! do not split molecules, infer correspondence, or transfer coordinate state
//! implicitly.

use std::fmt;
use std::sync::Arc;

use super::{
    InstanceAtomId, InstanceBondId, MoleculeDefinitionPayload, MoleculeInstanceId, SelectionError,
    Topology, TopologyBuildError, TopologyBuilder, TopologyEditResult, TopologyMapping,
    TopologyMappingError,
};

/// Retains complete molecule instances in source topology order.
///
/// Duplicate identifiers are treated as one membership request. A request
/// containing every source instance preserves the source `Arc<Topology>`. Empty
/// results and invalid identifiers are rejected before target construction.
///
/// # Examples
///
/// ```
/// use kekule::core::{Atom, Element, Molecule};
/// use kekule::small::SmallMolecule;
/// use kekule::topology::{
///     transform, MoleculeInstanceMetadata, MoleculeRole, TopologyBuilder,
/// };
/// use std::sync::Arc;
///
/// let mut water_builder = Molecule::builder();
/// water_builder.add_atom(Atom::new(Element::from_symbol("O").unwrap()))?;
/// let water = SmallMolecule::from_graph(water_builder.build()?);
///
/// let mut builder = TopologyBuilder::new();
/// let definition = builder.add_small_molecule_definition(&water)?;
/// let mut metadata = MoleculeInstanceMetadata::default();
/// metadata.insert_role(MoleculeRole::Solvent);
/// let first = builder.add_instance(definition, metadata.clone())?;
/// builder.add_instance(definition, metadata)?;
/// let source = Arc::new(builder.build()?);
///
/// let edit = transform::retain_instances(&source, [first])?;
/// assert_eq!(edit.topology().instance_count(), 1);
/// assert_eq!(edit.topology().definition_count(), 1);
/// assert!(!Arc::ptr_eq(&source, &edit.mapping().target_arc()));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn retain_instances(
    topology: &Arc<Topology>,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<TopologyEditResult, TopologyTransformError> {
    let retained = validate_instances(topology, instances)?;
    retain_normalized(topology, &retained)
}

/// Removes complete molecule instances while preserving filtered source order.
///
/// Duplicate identifiers are harmless. Removing no instances preserves the
/// source `Arc<Topology>`; removing every instance is rejected.
pub fn remove_instances(
    topology: &Arc<Topology>,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<TopologyEditResult, TopologyTransformError> {
    let removed = validate_instances(topology, instances)?;
    let retained = InstanceMembership {
        members: removed
            .members
            .into_iter()
            .map(|removed| !removed)
            .collect(),
        len: topology.instance_count() - removed.len,
    };
    retain_normalized(topology, &retained)
}

struct InstanceMembership {
    members: Vec<bool>,
    len: usize,
}

impl InstanceMembership {
    fn contains(&self, instance: MoleculeInstanceId) -> bool {
        self.members[instance.index()]
    }
}

fn validate_instances(
    topology: &Topology,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<InstanceMembership, TopologyTransformError> {
    let mut normalized = InstanceMembership {
        members: vec![false; topology.instance_count()],
        len: 0,
    };
    for instance in instances {
        if topology.instance(instance).is_err() {
            return Err(TopologyTransformError::InvalidSourceInstance(instance));
        }
        if !normalized.members[instance.index()] {
            normalized.members[instance.index()] = true;
            normalized.len += 1;
        }
    }
    Ok(normalized)
}

fn retain_normalized(
    topology: &Arc<Topology>,
    retained: &InstanceMembership,
) -> Result<TopologyEditResult, TopologyTransformError> {
    if retained.len == 0 {
        return Err(TopologyTransformError::EmptyTargetTopology);
    }
    if retained.len == topology.instance_count() {
        let mapping = TopologyMapping::between_identical_layouts(topology, topology)?;
        return Ok(TopologyEditResult::new(Arc::clone(topology), mapping)?);
    }

    let mut referenced_definitions = vec![false; topology.definition_count()];
    for (instance_id, instance) in topology.instances() {
        if retained.contains(instance_id) {
            referenced_definitions[instance.definition().index()] = true;
        }
    }
    let retained_definition_count = referenced_definitions
        .iter()
        .filter(|referenced| **referenced)
        .count();

    let mut builder = TopologyBuilder::new();
    builder.reserve_definitions(retained_definition_count)?;
    builder.reserve_instances(retained.len)?;

    let mut definition_targets = vec![None; topology.definition_count()];
    let mut definitions = Vec::with_capacity(retained_definition_count);
    for (source_id, definition) in topology
        .definitions()
        .filter(|(id, _)| referenced_definitions[id.index()])
    {
        let target_id = match definition.payload() {
            MoleculeDefinitionPayload::Small(molecule) => {
                builder.add_small_molecule_definition(molecule)?
            }
            MoleculeDefinitionPayload::Macro(molecule) => {
                builder.add_macro_molecule_definition(molecule)?
            }
        };
        definition_targets[source_id.index()] = Some(target_id);
        definitions.push((source_id, target_id));
    }

    let mut instances = Vec::with_capacity(retained.len);
    for (source_id, instance) in topology
        .instances()
        .filter(|(id, _)| retained.contains(*id))
    {
        let target_definition = definition_targets[instance.definition().index()]
            .expect("retained instance has a retained definition");
        let target_id = builder.add_instance(target_definition, instance.metadata().clone())?;
        instances.push((source_id, target_id));
    }

    let target = Arc::new(builder.build()?);
    let atoms = instances
        .iter()
        .flat_map(|(source_instance, target_instance)| {
            topology
                .graph_for_instance(*source_instance)
                .expect("retained source instance was validated")
                .atom_ids()
                .map(|atom| {
                    (
                        InstanceAtomId::new(*source_instance, atom),
                        InstanceAtomId::new(*target_instance, atom),
                    )
                })
        });
    let bonds = instances
        .iter()
        .flat_map(|(source_instance, target_instance)| {
            topology
                .graph_for_instance(*source_instance)
                .expect("retained source instance was validated")
                .bond_ids()
                .map(|bond| {
                    (
                        InstanceBondId::new(*source_instance, bond),
                        InstanceBondId::new(*target_instance, bond),
                    )
                })
        });
    let mapping = TopologyMapping::from_pairs(
        topology,
        &target,
        definitions,
        instances.iter().copied(),
        atoms,
        bonds,
    )?;
    Ok(TopologyEditResult::new(target, mapping)?)
}

/// Failure to construct an immutable whole-instance topology subset.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TopologyTransformError {
    /// A requested instance does not exist in the source topology.
    InvalidSourceInstance(MoleculeInstanceId),
    /// The requested membership would produce a topology with no instances.
    EmptyTargetTopology,
    /// The filtered target topology could not be constructed.
    TopologyBuild(TopologyBuildError),
    /// Complete checked edit lineage could not be constructed.
    Mapping(TopologyMappingError),
}

impl fmt::Display for TopologyTransformError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidSourceInstance(instance) => {
                write!(formatter, "invalid source molecule instance: {instance}")
            }
            Self::EmptyTargetTopology => {
                formatter.write_str("topology transformation would remove every instance")
            }
            Self::TopologyBuild(error) => {
                write!(formatter, "cannot build target topology: {error}")
            }
            Self::Mapping(error) => write!(formatter, "cannot build topology lineage: {error}"),
        }
    }
}

impl std::error::Error for TopologyTransformError {}

impl From<TopologyBuildError> for TopologyTransformError {
    fn from(error: TopologyBuildError) -> Self {
        Self::TopologyBuild(error)
    }
}

impl From<TopologyMappingError> for TopologyTransformError {
    fn from(error: TopologyMappingError) -> Self {
        Self::Mapping(error)
    }
}

/// Policy for source atoms removed while remapping an atom selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RemovedSelectionPolicy {
    /// Reject the remap and identify the first removed selected atom.
    Error,
    /// Explicitly discard removed selected atoms.
    Drop,
}

/// Failure to remap a topology-bound atom selection.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SelectionRemapError {
    /// The selection is not bound to the supplied source topology.
    SourceTopologyMismatch,
    /// The mapping is not sourced from the supplied source topology.
    MappingSourceMismatch,
    /// The mapping does not target the supplied target topology.
    MappingTargetMismatch,
    /// Strict policy encountered a selected atom removed by the edit.
    RemovedSelectedAtom(InstanceAtomId),
    /// The mapped target selection could not be constructed.
    Selection(SelectionError),
}

impl fmt::Display for SelectionRemapError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SourceTopologyMismatch => formatter
                .write_str("atom selection does not belong to the supplied source topology"),
            Self::MappingSourceMismatch => {
                formatter.write_str("topology mapping does not match the supplied source topology")
            }
            Self::MappingTargetMismatch => {
                formatter.write_str("topology mapping does not match the supplied target topology")
            }
            Self::RemovedSelectedAtom(atom) => {
                write!(formatter, "selected source atom {atom} was removed")
            }
            Self::Selection(error) => {
                write!(formatter, "cannot construct target selection: {error}")
            }
        }
    }
}

impl std::error::Error for SelectionRemapError {}
