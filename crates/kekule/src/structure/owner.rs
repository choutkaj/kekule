use crate::diagnostics::Cancellation;
use crate::DraftError;
use std::fmt;
use std::sync::Arc;

use crate::math::{PeriodicCell, Point3, RigidTransform};
use crate::molecule::Molecule;
use crate::properties::{GeometryProperties, TopologyPropertiesMut};
use crate::topology::TopologySubsetError;
use crate::topology::{
    AtomSelection, HierarchyBuilder, InstanceAtomId, InstanceBondId, MoleculeClass,
    MoleculeDefinitionId, MoleculeInstanceId, Topology, TopologyAtomIndex, TopologyBuildError,
    TopologyBuilder, TopologyPerceptionError,
};
use crate::units::Quantity;

use super::draft_geometry::DraftGeometry;
use super::scalars::{b_factor_quantities, b_factor_quantity, AtomScalars};
use super::state::{
    geometry_mutators, geometry_readers, is_bound, site_mutators, site_readers, subset_rows,
    transform_cell, transform_positions, validate_binding,
};
use super::{DenseArrayError, GeometryError, Positions, SetContext};

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
        StructureBuilder {
            topology: TopologyBuilder::from_shared(topology),
            state,
            extending_structure: true,
        }
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

pub(crate) mod sealed {
    /// Restricts [`super::AsStructureView`] to Kekule's structure owners and views.
    pub trait Sealed {}
}

/// Anything that can be borrowed as a [`StructureView`] without copying: structures,
/// views, geometries, and trajectory frames.
///
/// Each implementor also has an inherent `as_structure_view` method, so the trait
/// only needs importing for generic code. This trait is sealed; other
/// structure-like types, such as a [`crate::traj::FrameBuffer`], are passed
/// through their own `as_structure_view()`.
pub trait AsStructureView: sealed::Sealed {
    /// Borrows the structure as a [`StructureView`] without copying.
    fn as_structure_view(&self) -> StructureView<'_>;
}

impl sealed::Sealed for Structure {}
impl sealed::Sealed for StructureView<'_> {}

impl AsStructureView for Structure {
    fn as_structure_view(&self) -> StructureView<'_> {
        Structure::as_structure_view(self)
    }
}

impl AsStructureView for StructureView<'_> {
    fn as_structure_view(&self) -> StructureView<'_> {
        *self
    }
}

impl<'a> From<&'a Structure> for StructureView<'a> {
    fn from(structure: &'a Structure) -> Self {
        structure.as_structure_view()
    }
}

/// What a [`StructureView`] borrows: a structure, or one geometry of a structure
/// set or trajectory frame source with what its source records about it.
///
/// Kernels that only need coordinates ignore it. Writers and reductions read
/// it so that a geometry's place, velocities, and forces are
/// written, rejected, or deliberately omitted, never silently dropped.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum StructureKind<'a> {
    /// A [`Structure`], or a view materialized as one: a structure held on its
    /// own, with no set around it.
    Standalone,
    /// A geometry of a [`StructureSet`](super::StructureSet), a trajectory frame,
    /// or a frame buffer, with its [`SetContext`].
    InSet(SetContext<'a>),
}

/// A borrowed topology plus the geometry state bound to it: the common
/// input contract of coordinate-dependent kernels.
///
/// Structures, geometries, trajectory frames, and frame buffers all borrow as a
/// view, and [`Self::kind`] says whether it is a structure or a geometry.
/// Accessors return borrows with the view's full lifetime. A view of a
/// trajectory frame has no occupancies or B-factors. Views preserve exact
/// topology identity.
#[derive(Debug, Clone, Copy)]
pub struct StructureView<'a> {
    topology: &'a Arc<Topology>,
    positions: &'a Positions,
    cell: Option<&'a PeriodicCell>,
    occupancies: &'a AtomScalars,
    b_factors: &'a AtomScalars,
    properties: &'a GeometryProperties,
    kind: StructureKind<'a>,
}

impl<'a> StructureView<'a> {
    /// Only owners that already enforce binding may use this.
    pub(crate) fn bound(
        topology: &'a Arc<Topology>,
        positions: &'a Positions,
        cell: Option<&'a PeriodicCell>,
        occupancies: &'a AtomScalars,
        b_factors: &'a AtomScalars,
        properties: &'a GeometryProperties,
        kind: StructureKind<'a>,
    ) -> Self {
        debug_assert!(is_bound(positions.len(), properties, topology));
        Self {
            topology,
            positions,
            cell,
            occupancies,
            b_factors,
            properties,
            kind,
        }
    }

    /// Whether this view borrows a structure or a geometry, with the geometry's
    /// place, velocities, and forces.
    ///
    /// ```
    /// use kekule::{smiles, structure::{Structure, Positions, StructureKind}};
    /// let structure = Structure::new(smiles::to_topology("O")?, Positions::zeros(1))?;
    /// assert_eq!(structure.as_structure_view().kind(), StructureKind::Standalone);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub const fn kind(self) -> StructureKind<'a> {
        self.kind
    }

    /// The topology the viewed structure belongs to.
    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub(crate) const fn topology_arc(self) -> &'a Arc<Topology> {
        self.topology
    }

    /// A shared handle to the bound topology.
    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    /// Number of atoms, equal to the topology's atom count.
    pub fn atom_count(self) -> usize {
        self.positions.len()
    }

    /// Positions in dense topology atom order.
    pub fn positions(self) -> &'a Positions {
        self.positions
    }

    /// The periodic cell, or `None` for a non-periodic structure.
    pub fn cell(self) -> Option<&'a PeriodicCell> {
        self.cell
    }

    /// Structure annotations: owner values plus one row per dense atom and bond.
    pub fn properties(self) -> &'a GeometryProperties {
        self.properties
    }

    /// Occupancy of one atom; `None` when it was never recorded, which is
    /// always the case for a trajectory frame, or when `atom` is out of range.
    pub fn occupancy(self, atom: TopologyAtomIndex) -> Option<f64> {
        self.occupancies.get(atom)
    }

    /// Every occupancy in dense atom order, or `None` when none was recorded.
    pub fn occupancies(self) -> Option<&'a [Option<f64>]> {
        self.occupancies.values()
    }

    /// Isotropic B-factor of one atom; `None` when it was never recorded,
    /// which is always the case for a trajectory frame, or when `atom` is out
    /// of range.
    pub fn b_factor(self, atom: TopologyAtomIndex) -> Option<Quantity<f64>> {
        b_factor_quantity(self.b_factors.get(atom))
    }

    /// Every B-factor in dense atom order and square nanometres, or `None`
    /// when none was recorded.
    pub fn b_factors(self) -> Option<Quantity<&'a [Option<f64>]>> {
        b_factor_quantities(self.b_factors)
    }

    pub(crate) fn occupancy_column(self) -> &'a AtomScalars {
        self.occupancies
    }

    pub(crate) fn b_factor_column(self) -> &'a AtomScalars {
        self.b_factors
    }

    /// The position of one atom, or `None` when `atom` is not in the topology.
    pub fn position(self, atom: InstanceAtomId) -> Option<Quantity<Point3>> {
        self.positions.get(self.topology.atom_index(atom)?.index())
    }

    /// One atom with its topology identity and geometry state.
    pub fn atom(self, atom: InstanceAtomId) -> Option<StructureAtomView<'a>> {
        let atom = self.topology.atom(atom)?;
        Some(StructureAtomView {
            atom,
            structure: self,
        })
    }

    /// One atom by dense index.
    pub fn atom_at(self, index: TopologyAtomIndex) -> Option<StructureAtomView<'a>> {
        let atom = self.topology.atom_at(index)?;
        Some(StructureAtomView {
            atom,
            structure: self,
        })
    }

    /// Atoms with their geometry state, in dense order.
    pub fn atoms(self) -> impl ExactSizeIterator<Item = StructureAtomView<'a>> + 'a {
        self.topology.atoms().map(move |atom| StructureAtomView {
            atom,
            structure: self,
        })
    }

    /// Materializes this borrowed structure as an owned structure, dropping the
    /// set context reported by [`Self::kind`]. A view of a
    /// trajectory frame yields a structure without occupancies or B-factors.
    pub fn to_structure(self) -> Structure {
        Structure {
            topology: Arc::clone(self.topology),
            positions: self.positions.clone(),
            cell: self.cell.copied(),
            occupancies: self.occupancies.clone(),
            b_factors: self.b_factors.clone(),
            properties: self.properties.clone(),
        }
    }
}

/// One atom of a structure: its topology view (through `Deref`) plus the
/// coordinate-dependent state stored for it.
///
/// ```
/// use kekule::{smiles, structure::{Structure, Positions}};
/// let structure = Structure::new(smiles::to_topology("CO")?, Positions::zeros(2))?;
/// let oxygen = structure.atoms().find(|atom| atom.element().symbol() == "O").unwrap();
/// assert_eq!(oxygen.occupancy(), None);
/// assert_eq!(oxygen.position().into_value(), kekule::math::Point3::origin());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct StructureAtomView<'a> {
    atom: crate::topology::AtomView<'a>,
    structure: StructureView<'a>,
}

impl<'a> std::ops::Deref for StructureAtomView<'a> {
    type Target = crate::topology::AtomView<'a>;

    fn deref(&self) -> &Self::Target {
        &self.atom
    }
}

impl<'a> StructureAtomView<'a> {
    /// The topology view with the structure's full lifetime.
    pub const fn topology_atom(self) -> crate::topology::AtomView<'a> {
        self.atom
    }

    /// The atom's stored position in the canonical length unit (nm).
    pub fn position(self) -> Quantity<Point3> {
        self.structure
            .positions()
            .get(self.atom.index().index())
            .expect("bound structures cover every atom")
    }

    /// The atom's occupancy, or `None` when it was never recorded.
    pub fn occupancy(self) -> Option<f64> {
        self.structure.occupancy(self.atom.index())
    }

    /// The atom's isotropic B-factor in nm², or `None` when it was never
    /// recorded.
    pub fn b_factor(self) -> Option<Quantity<f64>> {
        self.structure.b_factor(self.atom.index())
    }

    /// One geometry annotation of this atom. Static annotations are
    /// available through [`crate::topology::AtomView::property`].
    pub fn geometry_property(
        self,
        key: &crate::properties::PropertyKey,
    ) -> Option<crate::properties::PropertyValue> {
        self.structure
            .properties()
            .atoms()
            .value(key, self.atom.index())
            .expect("bound structures cover every atom")
    }
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
