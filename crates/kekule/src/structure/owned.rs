use crate::diagnostics::Cancellation;
use std::fmt;
use std::sync::Arc;

use crate::math::{PeriodicCell, Point3};
use crate::molecule::Molecule;
use crate::properties::TopologyPropertiesMut;
use crate::topology::TopologySubsetError;
use crate::topology::{
    AtomSelection, InstanceAtomId, Topology, TopologyAtomIndex, TopologyPerceptionError,
};
use crate::units::Quantity;

use super::builder::{StructureBuildError, StructureBuilder};
use super::geometry_validation::is_bound;
use super::view::{StructureAtomView, StructureKind, StructureView};
use super::{subset_rows, DenseArrayError, Geometry, GeometryError, GeometryMut, Positions};

/// One concrete structure of one topology.
///
/// A structure is a shared [`Topology`] and one [`Geometry`] of it, in the
/// topology's dense atom order: a position for every atom, an optional
/// [`PeriodicCell`], optional crystallographic occupancies and B-factors,
/// velocities, forces, and geometry properties. The topology stays
/// coordinate-free and may be shared with other structures, structure sets,
/// trajectory frames, selections, and prepared potentials.
///
/// Geometry state reads through `Deref` (`structure.positions()`,
/// `structure.cell()`, ...) and changes through the guard
/// [`Self::geometry_mut`] (`structure.geometry_mut().set_cell(..)`), the same
/// [`GeometryMut`] a structure set lends, which keeps the geometry bound to
/// the topology. [`Self::into_parts`] and [`Self::from_parts`] move the
/// topology and the geometry out and back in. Static annotations live on
/// [`Self::topology`].
///
/// Use [`Self::from_molecule`] for a single connected molecule,
/// [`Self::builder`] to assemble several molecules and their positions, or
/// [`Self::new`] when a topology and correctly ordered positions exist.
///
/// ```
/// use kekule::{smiles, math::{Matrix3, PeriodicCell, RigidTransform, Vector3},
///     structure::{Structure, Positions}, units::{Quantity, NANOMETER}};
/// let mut structure = Structure::new(smiles::to_topology("O")?, Positions::zeros(1))?;
/// let cell = PeriodicCell::orthorhombic(
///     Quantity::new(Vector3::new(1.0, 2.0, 3.0), NANOMETER), [true; 3])?;
/// structure.geometry_mut().set_cell(Some(cell));
/// // A quarter turn about z moves the positions and rotates the cell.
/// let turn = Matrix3::from_columns(
///     Vector3::new(0.0, 1.0, 0.0), Vector3::new(-1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0));
/// let still = Quantity::new(Vector3::zero(), NANOMETER);
/// structure.geometry_mut().apply_rigid_transform(&RigidTransform::new(turn, still)?)?;
/// let [a, b, _] = structure.cell().unwrap().vectors().into_value();
/// assert_eq!((a, b), (Vector3::new(0.0, 1.0, 0.0), Vector3::new(-2.0, 0.0, 0.0)));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct Structure {
    pub(super) topology: Arc<Topology>,
    pub(super) geometry: Geometry,
}

/// Structures are equal when their topologies share a layout and their
/// geometries are equal.
impl PartialEq for Structure {
    fn eq(&self, other: &Self) -> bool {
        self.topology.shares_layout(&other.topology) && self.geometry == other.geometry
    }
}

impl std::ops::Deref for Structure {
    type Target = Geometry;

    fn deref(&self) -> &Self::Target {
        &self.geometry
    }
}

impl Structure {
    /// Binds dense positions to a topology, without a cell, occupancies,
    /// B-factors, velocities, forces, or properties.
    ///
    /// Positions must follow the topology's dense atom order.
    pub fn new(
        topology: impl Into<Arc<Topology>>,
        positions: Positions,
    ) -> Result<Self, StructureError> {
        Self::from_parts(topology, Geometry::new(positions))
    }

    /// A structure of `topology` and `geometry`, binding the geometry as a
    /// structure set binds the geometries it takes ([`Geometry::bind`]).
    ///
    /// # Errors
    ///
    /// [`StructureError::Geometry`] when the geometry does not have a
    /// position per atom of `topology`, or has bond rows of another topology.
    pub fn from_parts(
        topology: impl Into<Arc<Topology>>,
        mut geometry: Geometry,
    ) -> Result<Self, StructureError> {
        let topology = topology.into();
        geometry.bind(&topology)?;
        Ok(Self::bound(topology, geometry))
    }

    /// A structure of a geometry already bound to `topology`.
    pub(crate) fn bound(topology: Arc<Topology>, geometry: Geometry) -> Self {
        debug_assert!(is_bound(
            geometry.atom_count(),
            geometry.properties(),
            &topology
        ));
        Self { topology, geometry }
    }

    /// Moves the topology and the geometry out without copying. The
    /// geometry stays bound to the topology, so [`Self::from_parts`] takes
    /// them back.
    pub fn into_parts(self) -> (Arc<Topology>, Geometry) {
        (self.topology, self.geometry)
    }

    /// Dimension-preserving mutable access to the geometry. A geometry whose
    /// storage is shared with a clone is copied when it is first written.
    pub fn geometry_mut(&mut self) -> GeometryMut<'_> {
        GeometryMut::new(&mut self.geometry)
    }

    /// Returns this structure with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.geometry.set_cell(Some(cell));
        self
    }

    /// An empty [`StructureBuilder`]; equivalent to [`StructureBuilder::new`].
    pub fn builder() -> StructureBuilder {
        StructureBuilder::new()
    }

    /// Moves a structure into coordinated append-oriented construction state.
    pub fn into_builder(self) -> StructureBuilder {
        StructureBuilder::extending(self.topology, self.geometry)
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

    /// Borrows this structure as the common kernel input without copying.
    pub fn as_structure_view(&self) -> StructureView<'_> {
        self.geometry
            .structure_view(&self.topology, StructureKind::Standalone)
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
        let geometry = self.geometry.projected(&atoms, &bonds)?;
        Self::from_parts(subset.shared_topology(), geometry)
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
        Ok(self.geometry.set_position_at(index, position)?)
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
        if staged.is_empty() {
            return Ok(());
        }
        Ok(self
            .geometry
            .data_in_place()
            .positions
            .set_canonical_batch(&staged)?)
    }
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
