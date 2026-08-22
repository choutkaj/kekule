//! Immutable whole-instance topology transformations.
//!
//! These operations preserve complete molecule definitions and instances. They
//! do not split molecules, infer correspondence, or transfer coordinate state
//! implicitly.

use std::fmt;
use std::sync::Arc;

use super::{MoleculeInstanceId, Topology, TopologyBuildError, TopologyBuilder};

/// Retains complete molecule instances in source topology order.
///
/// Duplicate identifiers are treated as one membership request. A request
/// containing every source instance preserves the source `Arc<Topology>`. Empty
/// results and invalid identifiers are rejected before target construction.
///
/// # Examples
///
/// ```
/// use kekule::core::{Atom, Element, MoleculeEditor};
/// use kekule::topology::{transform, TopologyBuilder};
/// use std::sync::Arc;
///
/// let mut water_builder = MoleculeEditor::new();
/// water_builder.add_atom(Atom::new(Element::from_symbol("O").unwrap()))?;
/// let water = water_builder.finish()?;
///
/// let mut builder = TopologyBuilder::new();
/// let definition = builder.add_molecule_definition(&water)?;
/// let first = builder.add_instance(definition)?;
/// builder.add_instance(definition)?;
/// let source = Arc::new(builder.build()?);
///
/// let target = transform::retain_instances(&source, [first])?;
/// assert_eq!(target.instance_count(), 1);
/// assert_eq!(target.definition_count(), 1);
/// assert!(!Arc::ptr_eq(&source, &target));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn retain_instances(
    topology: &Arc<Topology>,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<Arc<Topology>, TopologyTransformError> {
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
) -> Result<Arc<Topology>, TopologyTransformError> {
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
) -> Result<Arc<Topology>, TopologyTransformError> {
    if retained.len == 0 {
        return Err(TopologyTransformError::EmptyTargetTopology);
    }
    if retained.len == topology.instance_count() {
        return Ok(Arc::clone(topology));
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
    for (source_id, definition) in topology
        .definitions()
        .filter(|(id, _)| referenced_definitions[id.index()])
    {
        let target_id = builder.add_molecule_definition(definition.molecule())?;
        definition_targets[source_id.index()] = Some(target_id);
    }

    for (_, instance) in topology
        .instances()
        .filter(|(id, _)| retained.contains(*id))
    {
        let target_definition = definition_targets[instance.definition().index()]
            .expect("retained instance has a retained definition");
        builder.add_instance(target_definition)?;
    }

    Ok(Arc::new(builder.build()?))
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
        }
    }
}

impl std::error::Error for TopologyTransformError {}

impl From<TopologyBuildError> for TopologyTransformError {
    fn from(error: TopologyBuildError) -> Self {
        Self::TopologyBuild(error)
    }
}
