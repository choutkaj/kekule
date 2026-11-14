use crate::diagnostics::Cancellation;
use std::fmt;
use std::sync::Arc;

use crate::math::{PeriodicCell, Point3, RigidTransform};
use crate::molecule::Molecule;
use crate::properties::{GeometryProperties, TopologyPropertiesMut};
use crate::topology::TopologySubsetError;
use crate::topology::{
    AtomSelection, InstanceAtomId, Topology, TopologyAtomIndex, TopologyPerceptionError,
};
use crate::units::Quantity;

use super::builder::{StructureBuildError, StructureBuilder};
use super::draft_geometry::DraftGeometry;
use super::scalars::AtomScalars;
use super::state::{
    geometry_mutators, geometry_readers, site_mutators, site_readers, subset_rows, transform_cell,
    transform_positions, validate_binding,
};
use super::view::{StructureAtomView, StructureKind, StructureView};
use super::{DenseArrayError, GeometryError, Positions};

/// One concrete structure of one topology.
///
/// A structure owns a shared [`Topology`] and one geometry of it, in the
/// topology's dense atom order: a position for every atom, an optional
/// [`PeriodicCell`], optional crystallographic occupancies and B-factors, and
/// geometry properties. The topology stays coordinate-free and may be
/// shared with other structures, geometries, trajectory frames, selections,
/// and prepared potentials.
///
/// Read and change geometry state through the structure's own
/// dimension-preserving methods (`structure.positions()`, `structure.set_cell(..)`,
/// `structure.properties_mut()`, ...). Static annotations live on
/// [`Self::topology`].
///
/// Use [`Self::from_molecule`] for a single connected molecule,
/// [`Self::builder`] to assemble several molecules and their positions, or
/// [`Self::new`] when a topology and correctly ordered positions exist.
#[derive(Debug, Clone)]
pub struct Structure {
    pub(super) topology: Arc<Topology>,
    pub(super) positions: Positions,
    pub(super) cell: Option<PeriodicCell>,
    pub(super) occupancies: AtomScalars,
    pub(super) b_factors: AtomScalars,
    pub(super) properties: GeometryProperties,
}

impl PartialEq for Structure {
    fn eq(&self, other: &Self) -> bool {
        self.topology.shares_layout(&other.topology)
            && self.positions == other.positions
            && self.cell == other.cell
            && self.occupancies == other.occupancies
            && self.b_factors == other.b_factors
            && self.properties == other.properties
    }
}

impl Structure {
    /// Binds dense positions to a topology, without a cell, occupancies,
    /// B-factors, or properties.
    ///
    /// Positions must follow the topology's dense atom order.
    pub fn new(
        topology: impl Into<Arc<Topology>>,
        positions: Positions,
    ) -> Result<Self, StructureError> {
        let properties = GeometryProperties::new(positions.len());
        Self::from_state(
            topology.into(),
            DraftGeometry {
                positions,
                cell: None,
                occupancies: AtomScalars::default(),
                b_factors: AtomScalars::default(),
                properties,
            },
        )
    }

    /// Validates staged state against `topology` and binds its bond rows.
    pub(super) fn from_state(
        topology: Arc<Topology>,
        state: DraftGeometry,
    ) -> Result<Self, StructureError> {
        validate_binding(state.positions.len(), &state.properties, &topology)?;
        let DraftGeometry {
            positions,
            cell,
            occupancies,
            b_factors,
            mut properties,
        } = state;
        properties.bind_bond_rows(topology.bond_count());
        Ok(Self {
            topology,
            positions,
            cell,
            occupancies,
            b_factors,
            properties,
        })
    }

    /// Moves the geometry state out as draft state, keeping the topology.
    pub(super) fn into_state(self) -> (Arc<Topology>, DraftGeometry) {
        (
            self.topology,
            DraftGeometry {
                positions: self.positions,
                cell: self.cell,
                occupancies: self.occupancies,
                b_factors: self.b_factors,
                properties: self.properties,
            },
        )
    }

    /// Returns this structure with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.cell = Some(cell);
        self
    }

    /// An empty [`StructureBuilder`]; equivalent to [`StructureBuilder::new`].
    pub fn builder() -> StructureBuilder {
        StructureBuilder::new()
    }

    /// Moves a structure into coordinated append-oriented construction state.
    pub fn into_builder(self) -> StructureBuilder {
        let (topology, state) = self.into_state();
        StructureBuilder::extending(topology, state)
    }

    /// Copies this structure into append-oriented construction state; see
    /// [`Self::into_builder`].
    pub fn to_builder(&self) -> StructureBuilder {
        self.clone().into_builder()
    }

    /// Builds a single-molecule structure from dense positions in molecule atom order.
    pub fn from_molecule(
        molecule: Molecule,
        positions: &Positions,
    ) -> Result<Self, StructureBuildError> {
        let mut builder = StructureBuilder::new();
        builder.add_molecule(molecule, positions)?;
        builder.build()
    }

    /// The coordinate-free topology of this structure.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// A shared handle to the topology; cloning the [`Arc`] copies no topology
    /// data and keeps the same layout identity.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Number of atoms, equal to the topology's atom count.
    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// Applies a rigid transform to the whole structure transactionally:
    /// positions move and the cell rotates, as for geometries and
    /// trajectory frames. Use [`Self::transform_atoms`] to move a selection
    /// while keeping the cell.
    ///
    /// ```
    /// use kekule::{smiles, math::{Matrix3, PeriodicCell, RigidTransform, Vector3},
    ///     structure::{Structure, Positions}, units::{Quantity, NANOMETER}};
    /// let mut structure = Structure::new(smiles::to_topology("O")?, Positions::zeros(1))?;
    /// let cell = PeriodicCell::orthorhombic(
    ///     Quantity::new(Vector3::new(1.0, 2.0, 3.0), NANOMETER), [true; 3])?;
    /// structure.set_cell(Some(cell));
    /// // A quarter turn about z.
    /// let turn = Matrix3::from_columns(
    ///     Vector3::new(0.0, 1.0, 0.0), Vector3::new(-1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0));
    /// let still = Quantity::new(Vector3::zero(), NANOMETER);
    /// structure.apply_rigid_transform(&RigidTransform::new(turn, still)?)?;
    /// let [a, b, _] = structure.cell().unwrap().vectors().into_value();
    /// assert_eq!((a, b), (Vector3::new(0.0, 1.0, 0.0), Vector3::new(-2.0, 0.0, 0.0)));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), GeometryError> {
        let cell = transform_cell(self.cell, *transform)?;
        let positions = transform_positions(&self.positions, *transform)?;
        self.cell = cell;
        self.positions = positions;
        Ok(())
    }

    /// Borrows this structure as the common kernel input without copying.
    pub fn as_structure_view(&self) -> StructureView<'_> {
        StructureView::bound(
            &self.topology,
            &self.positions,
            self.cell.as_ref(),
            &self.occupancies,
            &self.b_factors,
            &self.properties,
            StructureKind::Standalone,
        )
    }

    /// Installs default perception through a new topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`], perceiving each definition once.
    /// Structure state is kept without copying. Failure leaves the structure
    /// unchanged. The new snapshot shares the original layout, so selections,
    /// buffers, and prepared potentials bound before perception remain usable.
    ///
    /// ```
    /// use kekule::{smiles, structure::{Structure, Positions}};
    ///
    /// let topology = smiles::to_topology("c1ccccc1")?;
    /// let mut structure = Structure::new(topology, Positions::zeros(6))?;
    /// structure.perceive()?;
    /// assert!(structure.topology().instances().all(|m| m.molecule().perception().rings().is_some()));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    /// [`Self::perceive`], stopping partway once `cancellation` is requested.
    ///
    /// Delegates to [`Topology::perceived_cancellable`]. A cancelled call
    /// leaves the structure unchanged; a completed one installs exactly what
    /// [`Self::perceive`] would.
    ///
    /// ```
    /// use kekule::{diagnostics::Cancellation, perception::PerceptionError};
    /// use kekule::{smiles, structure::{Structure, Positions}};
    ///
    /// let topology = smiles::to_topology("c1ccccc1")?;
    /// let mut structure = Structure::new(topology, Positions::zeros(6))?;
    /// let error = structure.perceive_cancellable(&Cancellation::new(|| true)).unwrap_err();
    /// assert!(matches!(error.source, PerceptionError::Cancelled(_)));
    /// assert!(structure.topology().instances().all(|m| m.molecule().perception().rings().is_none()));
    /// structure.perceive_cancellable(&Cancellation::new(|| false))?;
    /// assert!(structure.topology().instances().all(|m| m.molecule().perception().rings().is_some()));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// As [`Topology::perceived_cancellable`].
    pub fn perceive_cancellable(
        &mut self,
        cancellation: &Cancellation,
    ) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived_cancellable(cancellation)?);
        Ok(())
    }

    /// Length-preserving mutable access to the topology's static annotations.
    ///
    /// A topology shared with other owners is first replaced by a snapshot of
    /// it ([`Arc::make_mut`]), so those owners are unaffected. Either way the
    /// layout is kept, so selections, buffers, and prepared potentials remain
    /// usable. See [`Topology::properties_mut`].
    pub fn topology_properties_mut(&mut self) -> TopologyPropertiesMut<'_> {
        Arc::make_mut(&mut self.topology).properties_mut()
    }

    /// Publishes an induced topology subset and projects geometry state
    /// onto it. Structure owner properties are dropped.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, StructureError> {
        let subset = crate::topology::Topology::subset(self.topology.clone(), selection)?;
        let (atoms, bonds) = subset_rows(&subset);
        let state = DraftGeometry {
            positions: self.positions.select_indices(&atoms)?,
            cell: self.cell,
            occupancies: self.occupancies.select(&atoms),
            b_factors: self.b_factors.select(&atoms),
            properties: self
                .properties
                .project(&atoms, &bonds)
                .map_err(GeometryError::from)?,
        };
        Self::from_state(subset.shared_topology(), state)
    }

    /// The position of one atom, or `None` when `atom` is not in the topology.
    pub fn position(&self, atom: InstanceAtomId) -> Option<Quantity<Point3>> {
        self.as_structure_view().position(atom)
    }

    /// One atom with its topology identity and geometry state; see
    /// [`StructureAtomView`].
    pub fn atom(&self, atom: InstanceAtomId) -> Option<StructureAtomView<'_>> {
        self.as_structure_view().atom(atom)
    }

    /// One atom by dense index, or `None` when `index` is out of range.
    pub fn atom_at(&self, index: TopologyAtomIndex) -> Option<StructureAtomView<'_>> {
        self.as_structure_view().atom_at(index)
    }

    /// Atoms with their geometry state, in dense order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = StructureAtomView<'_>> + '_ {
        self.as_structure_view().atoms()
    }

    /// Replaces the position of one atom, accepting any length unit.
    ///
    /// The structure is unchanged when an error is returned.
    ///
    /// # Errors
    ///
    /// Returns [`StructureError::InvalidAtomId`] when `atom` is not in the topology
    /// and [`StructureError::Geometry`] when the unit is not a length or the
    /// converted position is not finite.
    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), StructureError> {
        let index = atom_index(&self.topology, atom)?;
        Ok(self.set_position_at(index, position)?)
    }

    /// Applies a sparse coordinate batch atomically; repeated atoms use the
    /// last value. Only the batch is staged, never the complete array.
    pub fn set_atom_positions(
        &mut self,
        values: impl IntoIterator<Item = (InstanceAtomId, Quantity<Point3>)>,
    ) -> Result<(), StructureError> {
        let staged = values
            .into_iter()
            .map(|(atom, position)| {
                let index = atom_index(&self.topology, atom)?.index();
                let point = position
                    .into_unit(crate::units::CANONICAL_LENGTH_UNIT)
                    .map_err(DenseArrayError::from)?
                    .into_value();
                if !point.is_finite() {
                    return Err(DenseArrayError::NonFinite { index }.into());
                }
                Ok((index, point))
            })
            .collect::<Result<Vec<_>, StructureError>>()?;
        Ok(self.positions.set_canonical_batch(&staged)?)
    }

    fn storage_ref(&self) -> &Self {
        self
    }

    fn storage(&mut self) -> &mut Self {
        self
    }

    geometry_readers!();
    site_readers!();
    geometry_mutators!();
    site_mutators!();
}

pub(super) fn atom_index(
    topology: &Topology,
    atom: InstanceAtomId,
) -> Result<TopologyAtomIndex, StructureError> {
    topology
        .atom_index(atom)
        .ok_or(StructureError::InvalidAtomId(atom))
}

/// Why a structure operation was rejected.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum StructureError {
    /// The atom is not part of the structure's topology.
    InvalidAtomId(InstanceAtomId),
    /// Geometry state is invalid; see the source.
    Geometry(GeometryError),
    /// Subsetting the topology failed; see the source.
    Subset(TopologySubsetError),
}

impl fmt::Display for StructureError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(atom) => write!(
                formatter,
                "{} is not in the structure's topology.",
                crate::diagnostics::Capitalized(atom)
            ),
            Self::Geometry(error) => error.fmt(formatter),
            Self::Subset(_) => formatter.write_str("Cannot subset the structure's topology."),
        }
    }
}

impl std::error::Error for StructureError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Geometry(error) => error.source(),
            Self::Subset(error) => Some(error),
            _ => None,
        }
    }
}

impl From<GeometryError> for StructureError {
    fn from(error: GeometryError) -> Self {
        Self::Geometry(error)
    }
}

impl From<DenseArrayError> for StructureError {
    fn from(error: DenseArrayError) -> Self {
        Self::Geometry(error.into())
    }
}

impl From<TopologySubsetError> for StructureError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Subset(error)
    }
}
