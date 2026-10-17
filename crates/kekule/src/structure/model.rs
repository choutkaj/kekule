use crate::DraftError;
use std::fmt;
use std::sync::Arc;

use crate::geometry::{PeriodicCell, Point3};
use crate::molecule::Molecule;
use crate::properties::{StructureProperties, TopologyPropertiesMut};
use crate::topology::TopologySubsetError;
use crate::topology::{
    AtomSelection, Hierarchy, InstanceAtomId, InstanceBondId, MoleculeClass, MoleculeDefinitionId,
    MoleculeInstanceId, ResidueClass, ResidueId, Topology, TopologyAtomIndex, TopologyBuildError,
    TopologyBuilder, TopologyPerceptionError,
};
use crate::units::Quantity;

use super::model_state::ModelState;
use super::scalars::{b_factor_quantities, b_factor_quantity, AtomScalars};
use super::state::{
    is_bound, site_mutators, site_readers, structure_mutators, structure_readers, subset_rows,
    validate_binding,
};
use super::{DenseArrayError, Positions, StructureError};

/// One concrete structure of one topology.
///
/// A model owns a shared [`Topology`] and one structure of it, in the
/// topology's dense atom order: a position for every atom, an optional
/// [`PeriodicCell`], optional crystallographic occupancies and B-factors, and
/// structure properties. The topology stays coordinate-free and may be
/// shared with other models, ensemble members, trajectory frames, selections,
/// and prepared potentials.
///
/// Read and change structure state through the model's own
/// dimension-preserving methods (`model.positions()`, `model.set_cell(..)`,
/// `model.properties_mut()`, ...). Static annotations live on
/// [`Self::topology`].
///
/// Use [`Self::from_molecule`] for a single connected molecule,
/// [`Self::builder`] to assemble several molecules and their positions, or
/// [`Self::new`] when a topology and correctly ordered positions exist.
#[derive(Debug, Clone)]
pub struct Model {
    pub(super) topology: Arc<Topology>,
    pub(super) positions: Positions,
    pub(super) cell: Option<PeriodicCell>,
    pub(super) occupancies: AtomScalars,
    pub(super) b_factors: AtomScalars,
    pub(super) properties: StructureProperties,
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology)
            && self.positions == other.positions
            && self.cell == other.cell
            && self.occupancies == other.occupancies
            && self.b_factors == other.b_factors
            && self.properties == other.properties
    }
}

impl Model {
    /// Binds dense positions to a topology, without a cell, occupancies,
    /// B-factors, or properties.
    ///
    /// Positions must follow the topology's dense atom order.
    pub fn new(
        topology: impl Into<Arc<Topology>>,
        positions: Positions,
    ) -> Result<Self, ModelError> {
        let properties = StructureProperties::new(positions.len());
        Self::from_state(
            topology.into(),
            ModelState {
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
        state: ModelState,
    ) -> Result<Self, ModelError> {
        validate_binding(state.positions.len(), &state.properties, &topology)?;
        let ModelState {
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

    /// Moves the structure state out as draft state, keeping the topology.
    pub(super) fn into_state(self) -> (Arc<Topology>, ModelState) {
        (
            self.topology,
            ModelState {
                positions: self.positions,
                cell: self.cell,
                occupancies: self.occupancies,
                b_factors: self.b_factors,
                properties: self.properties,
            },
        )
    }

    /// Returns this model with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.cell = Some(cell);
        self
    }

    pub fn builder() -> ModelBuilder {
        ModelBuilder::new()
    }

    /// Moves a model into coordinated append-oriented construction state.
    pub fn into_builder(self) -> ModelBuilder {
        let (topology, state) = self.into_state();
        ModelBuilder {
            topology: TopologyBuilder::from_shared(topology),
            state,
            extending_model: true,
        }
    }

    pub fn to_builder(&self) -> ModelBuilder {
        self.clone().into_builder()
    }

    /// Builds a single-molecule model from dense positions in molecule atom order.
    pub fn from_molecule(
        molecule: Molecule,
        positions: &Positions,
    ) -> Result<Self, ModelBuildError> {
        let mut builder = ModelBuilder::new();
        builder.add_molecule(molecule, positions)?;
        builder.build()
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// Borrows this model as the common kernel input without copying.
    pub fn as_model_view(&self) -> ModelView<'_> {
        ModelView::bound(
            &self.topology,
            &self.positions,
            self.cell.as_ref(),
            &self.occupancies,
            &self.b_factors,
            &self.properties,
        )
    }

    /// Installs default perception through a new topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`], perceiving each definition once.
    /// Structure state is kept without copying. Failure leaves the model
    /// unchanged. The new snapshot shares the original layout, so selections,
    /// buffers, and prepared potentials bound before perception remain usable.
    ///
    /// ```
    /// use kekule::{smiles, structure::{Model, Positions}};
    ///
    /// let topology = smiles::to_topology("c1ccccc1")?;
    /// let mut model = Model::new(topology, Positions::zeros(6))?;
    /// model.perceive()?;
    /// assert!(model.topology().molecules().all(|m| m.molecule().perception().has_rings()));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
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

    /// Publishes an induced topology subset and projects structure state
    /// onto it. Structure owner properties are dropped.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, ModelError> {
        let subset = self.topology.subset(selection)?;
        let (atoms, bonds) = subset_rows(&subset);
        let state = ModelState {
            positions: self.positions.select_indices(&atoms)?,
            cell: self.cell,
            occupancies: self.occupancies.select(&atoms),
            b_factors: self.b_factors.select(&atoms),
            properties: self
                .properties
                .project(&atoms, &bonds)
                .map_err(StructureError::from)?,
        };
        Self::from_state(subset.shared_topology(), state)
    }

    pub fn position(&self, atom: InstanceAtomId) -> Result<Quantity<Point3>, ModelError> {
        self.as_model_view().position(atom)
    }

    /// One atom with its topology identity and structure state; see
    /// [`ModelAtomView`].
    pub fn atom(&self, atom: InstanceAtomId) -> Option<ModelAtomView<'_>> {
        self.as_model_view().atom(atom)
    }

    pub fn atom_at(&self, index: TopologyAtomIndex) -> Option<ModelAtomView<'_>> {
        self.as_model_view().atom_at(index)
    }

    /// Atoms with their structure state, in dense order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = ModelAtomView<'_>> + '_ {
        self.as_model_view().atoms()
    }

    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), ModelError> {
        let index = atom_index(&self.topology, atom)?;
        Ok(self.set_position_at(index, position)?)
    }

    /// Applies a sparse coordinate batch atomically; repeated atoms use the
    /// last value. Only the batch is staged, never the complete array.
    pub fn set_atom_positions(
        &mut self,
        values: impl IntoIterator<Item = (InstanceAtomId, Quantity<Point3>)>,
    ) -> Result<(), ModelError> {
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
            .collect::<Result<Vec<_>, ModelError>>()?;
        Ok(self.positions.set_canonical_batch(&staged)?)
    }

    fn storage_ref(&self) -> &Self {
        self
    }

    fn storage(&mut self) -> &mut Self {
        self
    }

    structure_readers!();
    site_readers!();
    structure_mutators!();
    site_mutators!();
}

pub(super) fn atom_index(
    topology: &Topology,
    atom: InstanceAtomId,
) -> Result<TopologyAtomIndex, ModelError> {
    topology
        .atom_index(atom)
        .ok_or(ModelError::InvalidAtomId(atom))
}

/// Anything that can be borrowed as a [`ModelView`] without copying: models,
/// views, ensemble members, trajectory frames, and frame buffers.
///
/// Each implementor also has an inherent `as_model_view` method, so the trait
/// only needs importing for generic code.
pub trait AsModelView {
    fn as_model_view(&self) -> ModelView<'_>;
}

impl AsModelView for Model {
    fn as_model_view(&self) -> ModelView<'_> {
        Model::as_model_view(self)
    }
}

impl AsModelView for ModelView<'_> {
    fn as_model_view(&self) -> ModelView<'_> {
        *self
    }
}

impl<'a> From<&'a Model> for ModelView<'a> {
    fn from(model: &'a Model) -> Self {
        model.as_model_view()
    }
}

/// A borrowed topology plus the structure state bound to it: the common
/// input contract of coordinate-dependent kernels.
///
/// Models, ensemble members, trajectory frames, and frame buffers all borrow
/// as a view. Accessors return borrows with the view's full lifetime. A view of
/// a trajectory frame has no occupancies or B-factors. Views preserve exact
/// topology identity.
#[derive(Debug, Clone, Copy)]
pub struct ModelView<'a> {
    topology: &'a Arc<Topology>,
    positions: &'a Positions,
    cell: Option<&'a PeriodicCell>,
    occupancies: &'a AtomScalars,
    b_factors: &'a AtomScalars,
    properties: &'a StructureProperties,
}

impl<'a> ModelView<'a> {
    /// Only owners that already enforce binding may use this.
    pub(crate) fn bound(
        topology: &'a Arc<Topology>,
        positions: &'a Positions,
        cell: Option<&'a PeriodicCell>,
        occupancies: &'a AtomScalars,
        b_factors: &'a AtomScalars,
        properties: &'a StructureProperties,
    ) -> Self {
        debug_assert!(is_bound(positions.len(), properties, topology));
        Self {
            topology,
            positions,
            cell,
            occupancies,
            b_factors,
            properties,
        }
    }

    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub(crate) const fn topology_arc(self) -> &'a Arc<Topology> {
        self.topology
    }

    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    pub fn atom_count(self) -> usize {
        self.positions.len()
    }

    pub fn positions(self) -> &'a Positions {
        self.positions
    }

    pub fn cell(self) -> Option<&'a PeriodicCell> {
        self.cell
    }

    pub fn properties(self) -> &'a StructureProperties {
        self.properties
    }

    /// Occupancy of one atom; `None` when it was never recorded, which is
    /// always the case for a trajectory frame.
    pub fn occupancy(self, atom: TopologyAtomIndex) -> Result<Option<f64>, StructureError> {
        self.occupancies.get(atom, self.positions.len())
    }

    /// Every occupancy in dense atom order, or `None` when none was recorded.
    pub fn occupancies(self) -> Option<&'a [Option<f64>]> {
        self.occupancies.values()
    }

    /// Isotropic B-factor of one atom; `None` when it was never recorded,
    /// which is always the case for a trajectory frame.
    pub fn b_factor(
        self,
        atom: TopologyAtomIndex,
    ) -> Result<Option<Quantity<f64>>, StructureError> {
        Ok(b_factor_quantity(
            self.b_factors.get(atom, self.positions.len())?,
        ))
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

    pub fn position(self, atom: InstanceAtomId) -> Result<Quantity<Point3>, ModelError> {
        let index = atom_index(self.topology, atom)?;
        Ok(self.positions.position_at(index.index())?)
    }

    /// One atom with its topology identity and structure state.
    pub fn atom(self, atom: InstanceAtomId) -> Option<ModelAtomView<'a>> {
        let atom = self.topology.atom(atom)?;
        Some(ModelAtomView { atom, model: self })
    }

    /// One atom by dense index.
    pub fn atom_at(self, index: TopologyAtomIndex) -> Option<ModelAtomView<'a>> {
        let atom = self.topology.atom_at(index)?;
        Some(ModelAtomView { atom, model: self })
    }

    /// Atoms with their structure state, in dense order.
    pub fn atoms(self) -> impl ExactSizeIterator<Item = ModelAtomView<'a>> + 'a {
        self.topology
            .atoms()
            .map(move |atom| ModelAtomView { atom, model: self })
    }

    /// Materializes this borrowed structure as an owned model. A view of a
    /// trajectory frame yields a model without occupancies or B-factors.
    pub fn to_model(self) -> Model {
        Model {
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
/// use kekule::{smiles, structure::{Model, Positions}};
/// let model = Model::new(smiles::to_topology("CO")?, Positions::zeros(2))?;
/// let oxygen = model.atoms().find(|atom| atom.element.symbol() == "O").unwrap();
/// assert_eq!(oxygen.occupancy(), None);
/// assert_eq!(oxygen.position().into_value(), kekule::geometry::Point3::origin());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct ModelAtomView<'a> {
    atom: crate::topology::AtomView<'a>,
    model: ModelView<'a>,
}

impl<'a> std::ops::Deref for ModelAtomView<'a> {
    type Target = crate::topology::AtomView<'a>;

    fn deref(&self) -> &Self::Target {
        &self.atom
    }
}

impl<'a> ModelAtomView<'a> {
    /// The topology view with the model's full lifetime.
    pub const fn topology_atom(self) -> crate::topology::AtomView<'a> {
        self.atom
    }

    pub fn position(self) -> Quantity<Point3> {
        self.model
            .positions()
            .position_at(self.atom.index().index())
            .expect("bound structures cover every atom")
    }

    pub fn occupancy(self) -> Option<f64> {
        self.model
            .occupancy(self.atom.index())
            .expect("bound structures cover every atom")
    }

    pub fn b_factor(self) -> Option<Quantity<f64>> {
        self.model
            .b_factor(self.atom.index())
            .expect("bound structures cover every atom")
    }

    /// One structure annotation of this atom. Static annotations are
    /// available through [`crate::topology::AtomView::property`].
    pub fn structure_property(
        self,
        key: &crate::properties::PropertyKey,
    ) -> Option<crate::properties::PropertyValue> {
        self.model
            .properties()
            .atoms()
            .value(key, self.atom.index())
            .expect("bound structures cover every atom")
    }
}

/// Why a model operation was rejected.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelError {
    /// The atom is not part of the model's topology.
    InvalidAtomId(InstanceAtomId),
    /// Structure state is invalid; see the source.
    Structure(StructureError),
    /// Subsetting the topology failed; see the source.
    Subset(TopologySubsetError),
}

impl fmt::Display for ModelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(atom) => write!(formatter, "invalid topology atom: {atom}"),
            Self::Structure(_) => formatter.write_str("invalid model structure"),
            Self::Subset(_) => formatter.write_str("cannot subset model topology"),
        }
    }
}

impl std::error::Error for ModelError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structure(error) => Some(error),
            Self::Subset(error) => Some(error),
            _ => None,
        }
    }
}

impl From<StructureError> for ModelError {
    fn from(error: StructureError) -> Self {
        Self::Structure(error)
    }
}

impl From<DenseArrayError> for ModelError {
    fn from(error: DenseArrayError) -> Self {
        Self::Structure(error.into())
    }
}

impl From<TopologySubsetError> for ModelError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Subset(error)
    }
}

/// Convenience builder that assembles topology and one complete model.
///
/// Each call to [`Self::add_molecule`] adds one molecule occurrence and its
/// dense positions. [`Self::build`] publishes the coordinate-free topology and
/// the matching model transactionally. Staged structure state follows staged
/// dense atom order; use [`Self::atom_index`] to address it by atom ID.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ModelBuilder {
    topology: TopologyBuilder,
    state: ModelState,
    extending_model: bool,
}

impl ModelBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn topology_builder(&self) -> &TopologyBuilder {
        &self.topology
    }

    /// A restricted static-annotation view. It cannot replace topology or add
    /// instances independently of their coordinates.
    ///
    /// ```compile_fail
    /// use kekule::{structure::ModelBuilder, topology::TopologyBuilder};
    /// let mut builder = ModelBuilder::new();
    /// *builder.topology_builder_mut() = TopologyBuilder::new();
    /// ```
    /// ```compile_fail
    /// use kekule::{structure::ModelBuilder, topology::MoleculeDefinitionId};
    /// let mut builder = ModelBuilder::new();
    /// builder.topology_builder_mut().add_instance(MoleculeDefinitionId::new(0));
    /// ```
    pub fn topology_builder_mut(&mut self) -> ModelTopologyMut<'_> {
        ModelTopologyMut {
            topology: &mut self.topology,
        }
    }

    pub fn hierarchy(&self) -> &Hierarchy {
        self.topology.hierarchy()
    }

    pub fn hierarchy_mut(&mut self) -> &mut Hierarchy {
        self.topology.hierarchy_mut()
    }

    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

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

    pub fn position(&self, atom: InstanceAtomId) -> Result<Quantity<Point3>, ModelBuildError> {
        let index = self.checked_atom_index(atom)?;
        Ok(self.state.positions.position_at(index.index())?)
    }

    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), ModelBuildError> {
        let index = self.checked_atom_index(atom)?;
        Ok(self
            .set_position_at(index, position)
            .map_err(ModelError::from)?)
    }

    /// Applies a sparse coordinate batch atomically; repeated atoms use the
    /// last value.
    pub fn set_atom_positions(
        &mut self,
        values: impl IntoIterator<Item = (InstanceAtomId, Quantity<Point3>)>,
    ) -> Result<(), ModelBuildError> {
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
            .collect::<Result<Vec<_>, ModelBuildError>>()?;
        Ok(self.state.positions.set_canonical_batch(&staged)?)
    }

    pub fn validate(&self) -> Result<(), ModelBuildError> {
        self.clone().build().map(|_| ())
    }

    /// Builds the model, or returns this builder for repair with the error.
    pub fn try_build(self) -> Result<Model, DraftError<Self, ModelBuildError>> {
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
    ) -> Result<(), ModelBuildError> {
        let previous = self.topology.reorder_atoms(order)?;
        self.state.reorder_atoms(&previous);
        Ok(())
    }

    /// Stages a reusable definition; the molecule moves into the topology.
    pub fn add_molecule_definition(
        &mut self,
        molecule: Molecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self.topology.add_molecule_definition(molecule)?)
    }

    /// Overrides automatic classification for one staged molecule definition.
    pub fn set_molecule_class(
        &mut self,
        definition: MoleculeDefinitionId,
        class: MoleculeClass,
    ) -> Result<(), ModelBuildError> {
        Ok(self.topology.set_molecule_class(definition, class)?)
    }

    /// Overrides automatic classification for one staged hierarchy residue.
    pub fn set_residue_class(
        &mut self,
        residue: ResidueId,
        class: ResidueClass,
    ) -> Result<(), ModelBuildError> {
        Ok(self.topology.set_residue_class(residue, class)?)
    }

    /// Adds an instance with dense positions in its definition's atom order.
    pub fn add_instance(
        &mut self,
        definition: MoleculeDefinitionId,
        positions: &Positions,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        let molecule = self.topology.definition(definition)?.molecule();
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
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        validate_position_count(molecule.atom_count(), positions.len())?;
        self.reserve(positions.len())?;
        let bonds = molecule.bond_count();
        let instance = self.topology.add_molecule(molecule)?;
        self.extend(positions, bonds);
        Ok(instance)
    }

    /// Validates and publishes the staged topology and model.
    pub fn build(self) -> Result<Model, ModelBuildError> {
        let topology = Arc::new(self.topology.build()?);
        Ok(Model::from_state(topology, self.state)?)
    }

    fn checked_atom_index(
        &self,
        atom: InstanceAtomId,
    ) -> Result<TopologyAtomIndex, ModelBuildError> {
        self.atom_index(atom)
            .ok_or_else(|| ModelError::InvalidAtomId(atom).into())
    }

    fn reserve(&mut self, additional: usize) -> Result<(), ModelBuildError> {
        self.state
            .try_reserve_atoms(additional)
            .map_err(|_| ModelBuildError::CapacityOverflow)
    }

    fn extend(&mut self, positions: &Positions, added_bonds: usize) {
        self.state.extend(positions, added_bonds);
        if self.extending_model {
            self.state.properties.owner_mut().clear();
        }
    }

    fn storage_ref(&self) -> &ModelState {
        &self.state
    }

    fn storage(&mut self) -> &mut ModelState {
        &mut self.state
    }

    structure_readers!();
    site_readers!();
    structure_mutators!();
    site_mutators!();
}

fn validate_position_count(expected: usize, actual: usize) -> Result<(), ModelBuildError> {
    if actual != expected {
        return Err(ModelBuildError::InstancePositionCountMismatch { expected, actual });
    }
    Ok(())
}

/// A model builder edit or the final build was rejected.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelBuildError {
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
    /// Final model validation failed after topology publication.
    Model(Box<ModelError>),
}

impl fmt::Display for ModelBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InstancePositionCountMismatch { expected, actual } => write!(
                formatter,
                "definition instance requires {expected} positions, but received {actual}"
            ),
            Self::CapacityOverflow => {
                formatter.write_str("model construction exceeds coordinate capacity")
            }
            Self::Topology(_) => formatter.write_str("cannot build topology"),
            Self::Hierarchy(_) => formatter.write_str("cannot build hierarchy"),
            Self::Model(_) => formatter.write_str("cannot build model"),
        }
    }
}

impl std::error::Error for ModelBuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Topology(e) => Some(e),
            Self::Hierarchy(e) => Some(e),
            Self::Model(e) => Some(e),
            _ => None,
        }
    }
}

impl From<DenseArrayError> for ModelBuildError {
    fn from(error: DenseArrayError) -> Self {
        ModelError::from(error).into()
    }
}

impl From<StructureError> for ModelBuildError {
    fn from(error: StructureError) -> Self {
        ModelError::from(error).into()
    }
}

/// Restricted mutable access to static model-builder annotations. Each method
/// consumes the view, so references remain bound to the original builder borrow.
#[derive(Debug)]
pub struct ModelTopologyMut<'a> {
    topology: &'a mut TopologyBuilder,
}

impl<'a> ModelTopologyMut<'a> {
    pub fn hierarchy_mut(self) -> &'a mut Hierarchy {
        self.topology.hierarchy_mut()
    }

    /// Length-preserving mutable access to staged static topology annotations.
    pub fn properties_mut(self) -> crate::properties::TopologyPropertiesMut<'a> {
        self.topology.properties_mut()
    }
}

impl From<ModelError> for ModelBuildError {
    fn from(error: ModelError) -> Self {
        Self::Model(Box::new(error))
    }
}

impl From<TopologyBuildError> for ModelBuildError {
    fn from(error: TopologyBuildError) -> Self {
        Self::Topology(error)
    }
}

impl From<crate::topology::HierarchyError> for ModelBuildError {
    fn from(error: crate::topology::HierarchyError) -> Self {
        Self::Hierarchy(error)
    }
}
