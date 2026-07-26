//! Immutable whole-instance topology transformations.
//!
//! These operations preserve complete molecule definitions and instances. They
//! do not split molecules, infer correspondence, or transfer coordinate state
//! implicitly.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use super::{
    InstanceAtomId, InstanceBondId, MoleculeDefinitionId, MoleculeDefinitionPayload,
    MoleculeInstanceId, SelectionError, Topology, TopologyBuildError, TopologyBuilder,
    TopologyEditResult, TopologyMapping, TopologyMappingError,
};

/// Retains complete molecule instances in source topology order.
///
/// Duplicate identifiers are treated as one membership request. A request
/// containing every source instance preserves exact topology identity. Empty
/// results and invalid identifiers are rejected before target construction.
///
/// # Examples
///
/// ```
/// use molecular::core::{Atom, Element, Molecule};
/// use molecular::small::SmallMolecule;
/// use molecular::topology::{
///     transform, MoleculeInstanceMetadata, MoleculeRole, TopologyBuilder,
/// };
///
/// let mut water = Molecule::new();
/// water.add_atom(Atom::new(Element::from_symbol("O").unwrap()))?;
/// let water = SmallMolecule::from_graph(water);
///
/// let mut builder = TopologyBuilder::new();
/// let definition = builder.add_small_molecule_definition(&water)?;
/// let mut metadata = MoleculeInstanceMetadata::default();
/// metadata.insert_role(MoleculeRole::Solvent);
/// let first = builder.add_instance(definition, metadata.clone())?;
/// builder.add_instance(definition, metadata)?;
/// let source = builder.build()?;
///
/// let edit = transform::retain_instances(&source, [first])?;
/// assert_eq!(edit.topology().instance_count(), 1);
/// assert_eq!(edit.topology().definition_count(), 1);
/// assert!(!edit.topology().same_identity(&source));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn retain_instances(
    topology: &Topology,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<TopologyEditResult, TopologyTransformError> {
    let retained = validate_instances(topology, instances)?;
    retain_normalized(topology, &retained)
}

/// Removes complete molecule instances while preserving filtered source order.
///
/// Duplicate identifiers are harmless. Removing no instances preserves exact
/// topology identity; removing every instance is rejected.
pub fn remove_instances(
    topology: &Topology,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<TopologyEditResult, TopologyTransformError> {
    let removed = validate_instances(topology, instances)?;
    let retained = topology
        .instances()
        .map(|(id, _)| id)
        .filter(|id| !removed.contains(id))
        .collect::<BTreeSet<_>>();
    retain_normalized(topology, &retained)
}

fn validate_instances(
    topology: &Topology,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<BTreeSet<MoleculeInstanceId>, TopologyTransformError> {
    let instances = instances.into_iter().collect::<BTreeSet<_>>();
    for instance in &instances {
        if topology.instance(*instance).is_err() {
            return Err(TopologyTransformError::InvalidSourceInstance(*instance));
        }
    }
    Ok(instances)
}

fn retain_normalized(
    topology: &Topology,
    retained: &BTreeSet<MoleculeInstanceId>,
) -> Result<TopologyEditResult, TopologyTransformError> {
    if retained.is_empty() {
        return Err(TopologyTransformError::EmptyTargetTopology);
    }
    if retained.len() == topology.instance_count() {
        let mapping = TopologyMapping::between_identical_layouts(topology, topology)?;
        return Ok(TopologyEditResult::new(topology.clone(), mapping)?);
    }

    let referenced_definitions = topology
        .instances()
        .filter(|(id, _)| retained.contains(id))
        .map(|(_, instance)| instance.definition())
        .collect::<BTreeSet<_>>();

    let mut builder = TopologyBuilder::new();
    builder.reserve_definitions(referenced_definitions.len())?;
    builder.reserve_instances(retained.len())?;

    let mut definitions = BTreeMap::<MoleculeDefinitionId, MoleculeDefinitionId>::new();
    for (source_id, definition) in topology
        .definitions()
        .filter(|(id, _)| referenced_definitions.contains(id))
    {
        let target_id = match definition.payload() {
            MoleculeDefinitionPayload::Small(molecule) => {
                builder.add_small_molecule_definition(molecule)?
            }
            MoleculeDefinitionPayload::Macro(molecule) => {
                builder.add_macro_molecule_definition(molecule)?
            }
        };
        definitions.insert(source_id, target_id);
    }

    let mut instances = Vec::with_capacity(retained.len());
    for (source_id, instance) in topology.instances().filter(|(id, _)| retained.contains(id)) {
        let target_definition = definitions[&instance.definition()];
        let target_id = builder.add_instance(target_definition, instance.metadata().clone())?;
        instances.push((source_id, target_id));
    }

    let target = builder.build()?;
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
