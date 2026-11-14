use crate::diagnostics::{CapacityExceeded, Countable, LengthMismatch};
use crate::draft::{Draft, DraftError};
use crate::fixed_u32_id::FixedId;
use std::fmt;
use std::sync::Arc;

use crate::math::Point3;
use crate::molecule::Molecule;
use crate::topology::{
    HierarchyBuilder, InstanceAtomId, MoleculeClass, MoleculeDefinitionId, MoleculeInstanceId,
    Topology, TopologyAtomIndex, TopologyBuildError, TopologyBuilder,
};
use crate::units::Quantity;

use super::{
    DenseArrayError, Geometry, GeometryData, GeometryError, GeometryMut, Positions, Structure,
    StructureError,
};

/// Convenience builder that assembles topology and one complete structure.
///
/// Each call to [`Self::add_molecule`] adds one molecule occurrence and its
/// dense positions. [`Self::build`] publishes the coordinate-free topology and
/// the matching structure transactionally. The staged [`Geometry`] follows
/// staged dense atom order: [`Self::geometry`] reads it and
/// [`Self::geometry_mut`] edits it; use [`Self::atom_index`] to address it by
/// atom ID. Adding atoms clears staged velocities and forces, which cover the
/// earlier atoms only.
#[derive(Debug, Clone, PartialEq)]
pub struct StructureBuilder {
    topology: TopologyBuilder,
    geometry: Geometry,
    extending_structure: bool,
}

impl Default for StructureBuilder {
    fn default() -> Self {
        Self {
            topology: TopologyBuilder::default(),
            geometry: Geometry::from_data(GeometryData::default()),
            extending_structure: false,
        }
    }
}

impl Draft for StructureBuilder {
    type Published = Structure;
    type Error = StructureBuildError;
    fn into_published(self) -> Result<Structure, StructureBuildError> {
        self.build()
    }
}

impl StructureBuilder {
    /// An empty builder with no molecules or positions.
    pub fn new() -> Self {
        Self::default()
    }

    /// A builder that extends `topology` and its bound `geometry`, as
    /// [`Structure::into_builder`] returns.
    pub(super) fn extending(topology: Arc<Topology>, geometry: Geometry) -> Self {
        Self {
            topology: TopologyBuilder::from_shared(topology),
            geometry,
            extending_structure: true,
        }
    }

    /// The staged geometry, in staged dense atom order.
    pub fn geometry(&self) -> &Geometry {
        &self.geometry
    }

    /// Dimension-preserving mutable access to the staged geometry.
    pub fn geometry_mut(&mut self) -> GeometryMut<'_> {
        GeometryMut::new(&mut self.geometry)
    }

    /// The staged topology.
    pub fn topology(&self) -> &TopologyBuilder {
        &self.topology
    }

    /// The staged topology edits that add no atoms. It cannot replace the
    /// topology or add instances without their coordinates.
    ///
    /// ```compile_fail
    /// use kekule::{structure::StructureBuilder, topology::TopologyBuilder};
    /// let mut builder = StructureBuilder::new();
    /// *builder.topology_mut() = TopologyBuilder::new();
    /// ```
    /// ```compile_fail
    /// use kekule::{structure::StructureBuilder, topology::MoleculeDefinitionId};
    /// let mut builder = StructureBuilder::new();
    /// builder.topology_mut().add_instance(MoleculeDefinitionId::new(0));
    /// ```
    pub fn topology_mut(&mut self) -> StructureBuilderTopologyMut<'_> {
        StructureBuilderTopologyMut {
            topology: &mut self.topology,
        }
    }

    /// Staged dense index of one atom, valid until the atom order changes.
    pub fn atom_index(&self, atom: InstanceAtomId) -> Option<TopologyAtomIndex> {
        self.topology
            .atom_index(atom)
            .map(TopologyAtomIndex::from_index)
    }

    /// The staged position of one atom, or `None` when `atom` is not in the
    /// staged topology.
    pub fn position(&self, atom: InstanceAtomId) -> Option<Quantity<Point3>> {
        self.geometry
            .positions()
            .get(self.atom_index(atom)?.index())
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
            .geometry
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
        if staged.is_empty() {
            return Ok(());
        }
        Ok(self
            .geometry
            .data_in_place()
            .positions
            .set_canonical_batch(&staged)?)
    }

    /// The error [`Self::build`] would return, leaving the builder unchanged.
    pub fn validate(&self) -> Result<(), StructureBuildError> {
        self.check()
    }

    /// Builds the structure, or returns this builder for repair with the error.
    pub fn try_build(self) -> Result<Structure, DraftError<Self, StructureBuildError>> {
        self.attempt()
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
        self.geometry.data_in_place().reorder_atoms(&previous);
        Ok(())
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
        Ok(Structure::from_parts(topology, self.geometry)?)
    }

    fn checked_atom_index(
        &self,
        atom: InstanceAtomId,
    ) -> Result<TopologyAtomIndex, StructureBuildError> {
        self.atom_index(atom)
            .ok_or_else(|| StructureError::InvalidAtomId(atom).into())
    }

    fn reserve(&mut self, additional: usize) -> Result<(), StructureBuildError> {
        self.geometry
            .data_in_place()
            .try_reserve_atoms(additional)
            .map_err(|_| {
                StructureBuildError::CapacityExceeded(CapacityExceeded::new(Countable::Atom))
            })
    }

    fn extend(&mut self, positions: &Positions, added_bonds: usize) {
        let data = self.geometry.data_in_place();
        data.extend(positions, added_bonds);
        if self.extending_structure {
            data.properties.owner_mut().clear();
        }
    }
}

fn validate_position_count(expected: usize, actual: usize) -> Result<(), StructureBuildError> {
    if actual != expected {
        return Err(StructureBuildError::InstancePositionCountMismatch(
            LengthMismatch::new(Countable::Position, expected, actual),
        ));
    }
    Ok(())
}

/// A structure builder edit or the final build was rejected.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum StructureBuildError {
    /// An added instance's positions do not match its molecule's atom count.
    InstancePositionCountMismatch(LengthMismatch),
    /// The structure cannot hold more atoms or residues.
    CapacityExceeded(CapacityExceeded),
    /// A topology edit or the topology build failed; see the source.
    Topology(TopologyBuildError),
    /// A chain, residue or atom-site edit failed; see the source.
    Hierarchy(crate::topology::HierarchyError),
    /// Final structure validation failed after topology publication.
    Structure(StructureError),
}

impl fmt::Display for StructureBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InstancePositionCountMismatch(error) => error.fmt(formatter),
            Self::CapacityExceeded(error) => error.fmt(formatter),
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

/// The edits of a [`StructureBuilder`]'s staged topology that add no atoms:
/// definitions, classes, hierarchy and static annotations. Accessors consume
/// the view, so what they return stays bound to the builder borrow.
#[derive(Debug)]
pub struct StructureBuilderTopologyMut<'a> {
    topology: &'a mut TopologyBuilder,
}

impl<'a> StructureBuilderTopologyMut<'a> {
    /// Length-preserving mutable access to staged static topology annotations.
    pub fn properties_mut(self) -> crate::properties::TopologyPropertiesMut<'a> {
        self.topology.properties_mut()
    }
    /// See [`TopologyBuilder::hierarchy_mut`].
    pub fn hierarchy_mut(self) -> &'a mut HierarchyBuilder {
        self.topology.hierarchy_mut()
    }
    /// Stages a reusable definition without an instance; see
    /// [`StructureBuilder::add_instance`].
    pub fn add_molecule_definition(
        &mut self,
        molecule: Molecule,
    ) -> Result<MoleculeDefinitionId, TopologyBuildError> {
        self.topology.add_molecule_definition(molecule)
    }
    /// See [`TopologyBuilder::set_molecule_class`].
    pub fn set_molecule_class(
        &mut self,
        definition: MoleculeDefinitionId,
        class: MoleculeClass,
    ) -> Result<(), TopologyBuildError> {
        self.topology.set_molecule_class(definition, class)
    }
    /// See [`TopologyBuilder::clear_molecule_class`].
    pub fn clear_molecule_class(
        &mut self,
        definition: MoleculeDefinitionId,
    ) -> Result<(), TopologyBuildError> {
        self.topology.clear_molecule_class(definition)
    }
}

impl From<StructureError> for StructureBuildError {
    fn from(error: StructureError) -> Self {
        Self::Structure(error)
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
