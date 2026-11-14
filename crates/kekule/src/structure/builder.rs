use crate::DraftError;
use std::fmt;
use std::sync::Arc;

use crate::math::Point3;
use crate::molecule::Molecule;
use crate::topology::{
    HierarchyBuilder, InstanceAtomId, InstanceBondId, MoleculeClass, MoleculeDefinitionId,
    MoleculeInstanceId, Topology, TopologyAtomIndex, TopologyBuildError, TopologyBuilder,
};
use crate::units::Quantity;

use super::draft_geometry::DraftGeometry;
use super::state::{geometry_mutators, geometry_readers, site_mutators, site_readers};
use super::{DenseArrayError, GeometryError, Positions, Structure, StructureError};

/// Convenience builder that assembles topology and one complete structure.
///
/// Each call to [`Self::add_molecule`] adds one molecule occurrence and its
/// dense positions. [`Self::build`] publishes the coordinate-free topology and
/// the matching structure transactionally. Staged geometry state follows staged
/// dense atom order; use [`Self::atom_index`] to address it by atom ID.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StructureBuilder {
    topology: TopologyBuilder,
    state: DraftGeometry,
    extending_structure: bool,
}

impl StructureBuilder {
    /// An empty builder with no molecules or positions.
    pub fn new() -> Self {
        Self::default()
    }

    /// A builder that extends `topology` and its geometry `state`, as
    /// [`Structure::into_builder`] returns.
    pub(super) fn extending(topology: Arc<Topology>, state: DraftGeometry) -> Self {
        Self {
            topology: TopologyBuilder::from_shared(topology),
            state,
            extending_structure: true,
        }
    }

    /// The staged topology, for read-only inspection.
    pub fn topology_builder(&self) -> &TopologyBuilder {
        &self.topology
    }

    /// A restricted static-annotation view. It cannot replace topology or add
    /// instances independently of their coordinates.
    ///
    /// ```compile_fail
    /// use kekule::{structure::StructureBuilder, topology::TopologyBuilder};
    /// let mut builder = StructureBuilder::new();
    /// *builder.topology_builder_mut() = TopologyBuilder::new();
    /// ```
    /// ```compile_fail
    /// use kekule::{structure::StructureBuilder, topology::MoleculeDefinitionId};
    /// let mut builder = StructureBuilder::new();
    /// builder.topology_builder_mut().add_instance(MoleculeDefinitionId::new(0));
    /// ```
    pub fn topology_builder_mut(&mut self) -> StructureTopologyMut<'_> {
        StructureTopologyMut {
            topology: &mut self.topology,
        }
    }

    /// The staged hierarchy.
    pub fn hierarchy(&self) -> &HierarchyBuilder {
        self.topology.hierarchy()
    }

    /// The staged hierarchy, for adding chains, residues, and atom sites and
    /// assigning residue classes; see [`TopologyBuilder::hierarchy_mut`].
    pub fn hierarchy_mut(&mut self) -> &mut HierarchyBuilder {
        self.topology.hierarchy_mut()
    }

    /// Number of staged atoms.
    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

    /// Number of staged bonds.
    pub fn bond_count(&self) -> usize {
        self.topology.bond_count()
    }

    /// Staged qualified atom IDs in current dense order.
    pub fn atom_ids(&self) -> &[InstanceAtomId] {
        self.topology.atom_ids()
    }

    /// Staged qualified bond IDs in dense order.
    pub fn bond_ids(&self) -> &[InstanceBondId] {
        self.topology.bond_ids()
    }

    /// Staged dense index of one atom, valid until the atom order changes.
    pub fn atom_index(&self, atom: InstanceAtomId) -> Option<TopologyAtomIndex> {
        self.topology
            .atom_index(atom)
            .map(|index| TopologyAtomIndex::new(index as u32))
    }

    /// The staged position of one atom, or `None` when `atom` is not in the
    /// staged topology.
    pub fn position(&self, atom: InstanceAtomId) -> Option<Quantity<Point3>> {
        self.state.positions.get(self.atom_index(atom)?.index())
    }

    /// Replaces the staged position of one atom, accepting any length unit.
    ///
    /// The builder is unchanged when an error is returned.
    ///
    /// # Errors
    ///
    /// Returns [`StructureBuildError::Structure`] wrapping [`StructureError::InvalidAtomId`]
    /// when `atom` is not staged, or [`StructureError::Geometry`] when the unit is
    /// not a length or the converted position is not finite.
    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), StructureBuildError> {
        let index = self.checked_atom_index(atom)?;
        Ok(self
            .set_position_at(index, position)
            .map_err(StructureError::from)?)
    }

    /// Applies a sparse coordinate batch atomically; repeated atoms use the
    /// last value.
    pub fn set_atom_positions(
        &mut self,
        values: impl IntoIterator<Item = (InstanceAtomId, Quantity<Point3>)>,
    ) -> Result<(), StructureBuildError> {
        let staged = values
            .into_iter()
            .map(|(atom, position)| {
                let index = self.checked_atom_index(atom)?.index();
                let point = position
                    .into_unit(crate::units::CANONICAL_LENGTH_UNIT)
                    .map_err(DenseArrayError::from)?
                    .into_value();
                if !point.is_finite() {
                    return Err(DenseArrayError::NonFinite { index }.into());
                }
                Ok((index, point))
            })
            .collect::<Result<Vec<_>, StructureBuildError>>()?;
        Ok(self.state.positions.set_canonical_batch(&staged)?)
    }

    /// Checks that [`Self::build`] would succeed, without consuming the builder.
    ///
    /// # Errors
    ///
    /// Returns the error [`Self::build`] would return.
    pub fn validate(&self) -> Result<(), StructureBuildError> {
        self.clone().build().map(|_| ())
    }

    /// Builds the structure, or returns this builder for repair with the error.
    pub fn try_build(self) -> Result<Structure, DraftError<Self, StructureBuildError>> {
        self.clone()
            .build()
            .map_err(|error| DraftError::new(error, self))
    }

    /// Replaces the dense atom order with a permutation of every staged atom.
    ///
    /// Staged structure atom state and topology atom properties move with
    /// their atoms. See [`TopologyBuilder::set_atom_order`]. A rejected order
    /// leaves the builder unchanged.
    pub fn set_atom_order(
        &mut self,
        order: impl IntoIterator<Item = InstanceAtomId>,
    ) -> Result<(), StructureBuildError> {
        let previous = self.topology.reorder_atoms(order)?;
        self.state.reorder_atoms(&previous);
        Ok(())
    }

    /// Stages a reusable definition; the molecule moves into the topology.
    pub fn add_molecule_definition(
        &mut self,
        molecule: Molecule,
    ) -> Result<MoleculeDefinitionId, StructureBuildError> {
        Ok(self.topology.add_molecule_definition(molecule)?)
    }

    /// Overrides automatic classification for one staged molecule definition.
    pub fn set_molecule_class(
        &mut self,
        definition: MoleculeDefinitionId,
        class: MoleculeClass,
    ) -> Result<(), StructureBuildError> {
        Ok(self.topology.set_molecule_class(definition, class)?)
    }

    /// Clears one staged definition's explicit class, so publication infers
    /// it afresh.
    pub fn clear_molecule_class(
        &mut self,
        definition: MoleculeDefinitionId,
    ) -> Result<(), StructureBuildError> {
        Ok(self.topology.clear_molecule_class(definition)?)
    }

    /// Adds an instance with dense positions in its definition's atom order.
    pub fn add_instance(
        &mut self,
        definition: MoleculeDefinitionId,
        positions: &Positions,
    ) -> Result<MoleculeInstanceId, StructureBuildError> {
        let molecule = self.topology.require_definition(definition)?.molecule();
        let added_bonds = molecule.bond_count();
        validate_position_count(molecule.atom_count(), positions.len())?;
        self.reserve(positions.len())?;
        let instance = self.topology.add_instance(definition)?;
        self.extend(positions, added_bonds);
        Ok(instance)
    }

    /// Adds a molecule and dense positions in that molecule's atom order.
    pub fn add_molecule(
        &mut self,
        molecule: Molecule,
        positions: &Positions,
    ) -> Result<MoleculeInstanceId, StructureBuildError> {
        validate_position_count(molecule.atom_count(), positions.len())?;
        self.reserve(positions.len())?;
        let bonds = molecule.bond_count();
        let instance = self.topology.add_molecule(molecule)?;
        self.extend(positions, bonds);
        Ok(instance)
    }

    /// Validates and publishes the staged topology and structure.
    pub fn build(self) -> Result<Structure, StructureBuildError> {
        let topology = Arc::new(self.topology.build()?);
        Ok(Structure::from_state(topology, self.state)?)
    }

    fn checked_atom_index(
        &self,
        atom: InstanceAtomId,
    ) -> Result<TopologyAtomIndex, StructureBuildError> {
        self.atom_index(atom)
            .ok_or_else(|| StructureError::InvalidAtomId(atom).into())
    }

    fn reserve(&mut self, additional: usize) -> Result<(), StructureBuildError> {
        self.state
            .try_reserve_atoms(additional)
            .map_err(|_| StructureBuildError::CapacityOverflow)
    }

    fn extend(&mut self, positions: &Positions, added_bonds: usize) {
        self.state.extend(positions, added_bonds);
        if self.extending_structure {
            self.state.properties.owner_mut().clear();
        }
    }

    fn storage_ref(&self) -> &DraftGeometry {
        &self.state
    }

    fn storage(&mut self) -> &mut DraftGeometry {
        &mut self.state
    }

    geometry_readers!();
    site_readers!();
    geometry_mutators!();
    site_mutators!();
}

fn validate_position_count(expected: usize, actual: usize) -> Result<(), StructureBuildError> {
    if actual != expected {
        return Err(StructureBuildError::InstancePositionCountMismatch { expected, actual });
    }
    Ok(())
}

/// A structure builder edit or the final build was rejected.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum StructureBuildError {
    /// An added instance's positions do not match its molecule's atom count.
    InstancePositionCountMismatch {
        /// The molecule's atom count.
        expected: usize,
        /// The number of positions supplied.
        actual: usize,
    },
    /// Reserving per-atom state exceeded capacity.
    CapacityOverflow,
    /// A topology edit or the topology build failed; see the source.
    Topology(TopologyBuildError),
    /// A chain, residue or atom-site edit failed; see the source.
    Hierarchy(crate::topology::HierarchyError),
    /// Final structure validation failed after topology publication.
    Structure(Box<StructureError>),
}

impl fmt::Display for StructureBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InstancePositionCountMismatch { expected, actual } => write!(
                formatter,
                "The molecule has {}, but {} given.",
                crate::diagnostics::counted(*expected, "atom", "atoms"),
                crate::diagnostics::counted(*actual, "position was", "positions were")
            ),
            Self::CapacityOverflow => {
                formatter.write_str("The structure exceeds the addressable number of atoms.")
            }
            Self::Topology(_) => formatter.write_str("Cannot build the topology."),
            Self::Hierarchy(_) => formatter.write_str("Cannot build the hierarchy."),
            Self::Structure(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for StructureBuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Topology(e) => Some(e),
            Self::Hierarchy(e) => Some(e),
            Self::Structure(e) => e.source(),
            _ => None,
        }
    }
}

impl From<DenseArrayError> for StructureBuildError {
    fn from(error: DenseArrayError) -> Self {
        StructureError::from(error).into()
    }
}

impl From<GeometryError> for StructureBuildError {
    fn from(error: GeometryError) -> Self {
        StructureError::from(error).into()
    }
}

/// Restricted mutable access to static structure-builder annotations. Each method
/// consumes the view, so references remain bound to the original builder borrow.
#[derive(Debug)]
pub struct StructureTopologyMut<'a> {
    topology: &'a mut TopologyBuilder,
}

impl<'a> StructureTopologyMut<'a> {
    /// Length-preserving mutable access to staged static topology annotations.
    pub fn properties_mut(self) -> crate::properties::TopologyPropertiesMut<'a> {
        self.topology.properties_mut()
    }
}

impl From<StructureError> for StructureBuildError {
    fn from(error: StructureError) -> Self {
        Self::Structure(Box::new(error))
    }
}

impl From<TopologyBuildError> for StructureBuildError {
    fn from(error: TopologyBuildError) -> Self {
        Self::Topology(error)
    }
}

impl From<crate::topology::HierarchyError> for StructureBuildError {
    fn from(error: crate::topology::HierarchyError) -> Self {
        Self::Hierarchy(error)
    }
}
