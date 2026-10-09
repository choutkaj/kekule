use std::fmt;
use std::sync::Arc;

use crate::core::Molecule;
use crate::geometry::{PeriodicCell, Point3};
use crate::topology::transform::TopologySubsetError;
use crate::topology::{
    AtomSelection, Hierarchy, InstanceAtomId, InstanceBondId, MoleculeClass, MoleculeDefinitionId,
    MoleculeInstanceId, ResidueClass, ResidueId, Topology, TopologyAtomIndex, TopologyBuildError,
    TopologyBuilder, TopologyPerceptionError,
};
use crate::units::Quantity;

use super::{Conformation, ConformationError, ConformationMut, PositionError, Positions};

/// One concrete realization of one immutable topology.
///
/// A model owns a shared [`Topology`] and one bound [`Conformation`]: a
/// position for every topology atom, an optional [`PeriodicCell`], optional
/// occupancies and B-factors, and realization properties. The topology stays
/// coordinate-free and may be shared with other models, ensemble members,
/// trajectory frames, selections, and prepared potentials.
///
/// Conformation state reads through `Deref` (`model.positions()`,
/// `model.cell()`, `model.properties()`); mutate it through
/// [`Self::conformation_mut`]. Static annotations live on
/// [`Self::topology`].
///
/// Use [`Self::from_molecule`] for a single connected molecule,
/// [`Self::builder`] to assemble several molecules and their positions, or
/// [`Self::new`] when a topology and correctly ordered positions exist.
#[derive(Debug, Clone)]
pub struct Model {
    pub(super) topology: Arc<Topology>,
    pub(super) conformation: Conformation,
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.conformation == other.conformation
    }
}

impl std::ops::Deref for Model {
    type Target = Conformation;

    fn deref(&self) -> &Self::Target {
        &self.conformation
    }
}

impl Model {
    /// Binds a conformation (or plain [`Positions`]) to a topology.
    ///
    /// Positions must follow the topology's dense atom order. Bond property
    /// rows of a detached conformation are allocated here.
    pub fn new(
        topology: impl Into<Arc<Topology>>,
        conformation: impl Into<Conformation>,
    ) -> Result<Self, ModelError> {
        let topology = topology.into();
        let mut conformation = conformation.into();
        conformation.validate_for(&topology)?;
        conformation.bind(&topology);
        Ok(Self {
            topology,
            conformation,
        })
    }

    pub fn builder() -> ModelBuilder {
        ModelBuilder::new()
    }

    /// Moves a model into coordinated append-oriented construction state.
    pub fn into_builder(self) -> ModelBuilder {
        ModelBuilder {
            topology: TopologyBuilder::from_shared(self.topology),
            conformation: self.conformation,
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

    pub fn conformation(&self) -> &Conformation {
        &self.conformation
    }

    /// Borrows this model as the common kernel input without copying.
    pub fn as_model_view(&self) -> ModelView<'_> {
        ModelView::bound(&self.topology, &self.conformation)
    }

    /// Dimension-preserving mutable access to the conformation.
    pub fn conformation_mut(&mut self) -> ConformationMut<'_> {
        ConformationMut::new(&mut self.conformation)
    }

    /// Replaces the conformation after validating it against the topology,
    /// returning the previous one.
    pub fn replace_conformation(
        &mut self,
        conformation: impl Into<Conformation>,
    ) -> Result<Conformation, ModelError> {
        let mut conformation = conformation.into();
        conformation.validate_for(&self.topology)?;
        conformation.bind(&self.topology);
        Ok(std::mem::replace(&mut self.conformation, conformation))
    }

    /// Consumes the model without copying either part.
    pub fn into_parts(self) -> (Arc<Topology>, Conformation) {
        (self.topology, self.conformation)
    }

    /// Installs default perception through a new topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`], perceiving each definition once.
    /// The conformation is kept without copying. Failure leaves the model
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

    /// Publishes an induced topology subset and projects the conformation onto
    /// it. Realization owner properties are dropped.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, ModelError> {
        let subset = self.topology.subset(selection)?;
        let correspondence = subset.correspondence();
        let atoms = correspondence
            .source_atom_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let bonds = correspondence
            .source_bond_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let conformation = self.conformation.project(&atoms, &bonds)?;
        Self::new(subset.shared_topology(), conformation)
    }

    pub fn position(&self, atom: InstanceAtomId) -> Result<Quantity<Point3>, ModelError> {
        self.as_model_view().position(atom)
    }

    /// One atom with its topology identity and realization state; see
    /// [`ModelAtomView`].
    pub fn atom(&self, atom: InstanceAtomId) -> Option<ModelAtomView<'_>> {
        self.as_model_view().atom(atom)
    }

    pub fn atom_at(&self, index: TopologyAtomIndex) -> Option<ModelAtomView<'_>> {
        self.as_model_view().atom_at(index)
    }

    /// Atoms with their realization state, in dense order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = ModelAtomView<'_>> + '_ {
        self.as_model_view().atoms()
    }

    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), ModelError> {
        let index = atom_index(&self.topology, atom)?;
        Ok(self.conformation.set_position(index, position)?)
    }

    /// Applies a sparse coordinate batch atomically; repeated atoms use the
    /// last value. Only the batch is staged, never the complete array.
    pub fn set_atom_positions(
        &mut self,
        values: impl IntoIterator<Item = (InstanceAtomId, Quantity<Point3>)>,
    ) -> Result<(), ModelError> {
        let staged = values
            .into_iter()
            .enumerate()
            .map(|(order, (atom, position))| {
                let index = atom_index(&self.topology, atom)?.index();
                let point = position
                    .into_unit(crate::units::CANONICAL_LENGTH_UNIT)
                    .map_err(PositionError::from)?
                    .into_value();
                if !point.is_finite() {
                    return Err(PositionError::NonFinitePosition { index: order }.into());
                }
                Ok((index, point))
            })
            .collect::<Result<Vec<_>, ModelError>>()?;
        Ok(self
            .conformation
            .positions_mut()
            .set_canonical_batch(&staged)?)
    }
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

/// A borrowed topology plus bound conformation: the common input contract of
/// coordinate-dependent kernels.
///
/// Conformation state reads through `Deref`. The inherent accessors below
/// return borrows with the view's full lifetime. Views preserve exact
/// topology identity.
#[derive(Debug, Clone, Copy)]
pub struct ModelView<'a> {
    topology: &'a Arc<Topology>,
    conformation: &'a Conformation,
}

impl std::ops::Deref for ModelView<'_> {
    type Target = Conformation;

    fn deref(&self) -> &Self::Target {
        self.conformation
    }
}

impl<'a> ModelView<'a> {
    /// Borrows a conformation that already has every row `topology` requires.
    pub fn new(
        topology: &'a Arc<Topology>,
        conformation: &'a Conformation,
    ) -> Result<Self, ConformationError> {
        conformation.validate_for(topology)?;
        if !conformation.is_bound_to(topology) {
            return Err(ConformationError::BondCountMismatch {
                expected: topology.bond_count(),
                actual: conformation.properties().bonds().len(),
            });
        }
        Ok(Self::bound(topology, conformation))
    }

    /// Only owners that already enforce binding may use this.
    pub(crate) fn bound(topology: &'a Arc<Topology>, conformation: &'a Conformation) -> Self {
        debug_assert!(conformation.is_bound_to(topology));
        Self {
            topology,
            conformation,
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

    pub fn conformation(self) -> &'a Conformation {
        self.conformation
    }

    pub fn positions(self) -> &'a Positions {
        self.conformation.positions()
    }

    pub fn cell(self) -> Option<&'a PeriodicCell> {
        self.conformation.cell()
    }

    pub fn properties(self) -> &'a crate::properties::RealizationProperties {
        self.conformation.properties()
    }

    pub fn position(self, atom: InstanceAtomId) -> Result<Quantity<Point3>, ModelError> {
        let index = atom_index(self.topology, atom)?;
        Ok(self.conformation.positions().position_at(index.index())?)
    }

    /// One atom with its topology identity and realization state.
    pub fn atom(self, atom: InstanceAtomId) -> Option<ModelAtomView<'a>> {
        let atom = self.topology.atom(atom)?;
        Some(ModelAtomView {
            atom,
            conformation: self.conformation,
        })
    }

    /// One atom by dense index.
    pub fn atom_at(self, index: TopologyAtomIndex) -> Option<ModelAtomView<'a>> {
        let atom = self.topology.atom_at(index)?;
        Some(ModelAtomView {
            atom,
            conformation: self.conformation,
        })
    }

    /// Atoms with their realization state, in dense order.
    pub fn atoms(self) -> impl ExactSizeIterator<Item = ModelAtomView<'a>> + 'a {
        let conformation = self.conformation;
        self.topology
            .atoms()
            .map(move |atom| ModelAtomView { atom, conformation })
    }

    /// Materializes this borrowed realization as an owned model.
    pub fn to_model(self) -> Model {
        Model {
            topology: Arc::clone(self.topology),
            conformation: self.conformation.clone(),
        }
    }
}

/// One atom of a realization: its topology view (through `Deref`) plus the
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
    conformation: &'a Conformation,
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
        self.conformation
            .positions()
            .position_at(self.atom.index().index())
            .expect("bound conformations cover every atom")
    }

    pub fn occupancy(self) -> Option<f64> {
        self.conformation
            .occupancy(self.atom.index())
            .expect("bound conformations cover every atom")
    }

    pub fn b_factor(self) -> Option<Quantity<f64>> {
        self.conformation
            .b_factor(self.atom.index())
            .expect("bound conformations cover every atom")
    }

    /// One realization annotation of this atom. Static annotations are
    /// available through [`crate::topology::AtomView::property`].
    pub fn realization_property(
        self,
        key: &crate::properties::PropertyKey,
    ) -> Option<crate::properties::PropertyValue> {
        self.conformation
            .properties()
            .atoms()
            .value(key, self.atom.index())
            .expect("bound conformations cover every atom")
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelError {
    InvalidAtomId(InstanceAtomId),
    InvalidBondId(InstanceBondId),
    Conformation(ConformationError),
    Subset(TopologySubsetError),
}

impl fmt::Display for ModelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(atom) => write!(formatter, "invalid topology atom: {atom}"),
            Self::InvalidBondId(bond) => write!(formatter, "invalid topology bond: {bond}"),
            Self::Conformation(error) => write!(formatter, "invalid model conformation: {error}"),
            Self::Subset(error) => write!(formatter, "cannot subset model topology: {error}"),
        }
    }
}

impl std::error::Error for ModelError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Conformation(error) => Some(error),
            Self::Subset(error) => Some(error),
            _ => None,
        }
    }
}

impl From<ConformationError> for ModelError {
    fn from(error: ConformationError) -> Self {
        Self::Conformation(error)
    }
}

impl From<PositionError> for ModelError {
    fn from(error: PositionError) -> Self {
        Self::Conformation(error.into())
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
/// the matching model transactionally. The staged conformation follows staged
/// dense atom order; use [`Self::atom_index`] to address it by atom ID.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelBuilder {
    topology: TopologyBuilder,
    conformation: Conformation,
    extending_model: bool,
}

impl Default for ModelBuilder {
    fn default() -> Self {
        Self {
            topology: TopologyBuilder::default(),
            conformation: Conformation::new(Positions::default()),
            extending_model: false,
        }
    }
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

    /// Staged conformation in staged dense atom order.
    pub fn conformation(&self) -> &Conformation {
        &self.conformation
    }

    pub fn conformation_mut(&mut self) -> ConformationMut<'_> {
        ConformationMut::new(&mut self.conformation)
    }

    pub fn position(&self, atom: InstanceAtomId) -> Result<Quantity<Point3>, ModelBuildError> {
        let index = self.checked_atom_index(atom)?;
        Ok(self.conformation.positions().position_at(index.index())?)
    }

    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), ModelBuildError> {
        let index = self.checked_atom_index(atom)?;
        Ok(self
            .conformation
            .set_position(index, position)
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
            .enumerate()
            .map(|(order, (atom, position))| {
                let index = self.checked_atom_index(atom)?.index();
                let point = position
                    .into_unit(crate::units::CANONICAL_LENGTH_UNIT)
                    .map_err(PositionError::from)?
                    .into_value();
                if !point.is_finite() {
                    return Err(PositionError::NonFinitePosition { index: order }.into());
                }
                Ok((index, point))
            })
            .collect::<Result<Vec<_>, ModelBuildError>>()?;
        Ok(self
            .conformation
            .positions_mut()
            .set_canonical_batch(&staged)?)
    }

    pub fn validate(&self) -> Result<(), ModelBuildError> {
        self.clone().build().map(|_| ())
    }

    pub fn try_build(self) -> Result<Model, ModelBuilderError> {
        self.clone().build().map_err(|error| ModelBuilderError {
            error,
            builder: Box::new(self),
        })
    }

    /// Replaces the dense atom order with a permutation of every staged atom.
    ///
    /// Conformation atom state and staged topology atom properties move with
    /// their atoms. See [`TopologyBuilder::set_atom_order`]. A rejected order
    /// leaves the builder unchanged.
    pub fn set_atom_order(
        &mut self,
        order: impl IntoIterator<Item = InstanceAtomId>,
    ) -> Result<(), ModelBuildError> {
        let previous = self.topology.reorder_atoms(order)?;
        self.conformation.reorder_atoms(&previous);
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
        Ok(Model::new(topology, self.conformation)?)
    }

    fn checked_atom_index(
        &self,
        atom: InstanceAtomId,
    ) -> Result<TopologyAtomIndex, ModelBuildError> {
        self.atom_index(atom)
            .ok_or_else(|| ModelError::InvalidAtomId(atom).into())
    }

    fn reserve(&mut self, additional: usize) -> Result<(), ModelBuildError> {
        self.conformation
            .try_reserve_atoms(additional)
            .map_err(|_| ModelBuildError::CapacityOverflow)
    }

    fn extend(&mut self, positions: &Positions, added_bonds: usize) {
        self.conformation.extend(positions, added_bonds);
        if self.extending_model {
            self.conformation
                .properties_storage_mut()
                .owner_mut()
                .clear();
        }
    }
}

fn validate_position_count(expected: usize, actual: usize) -> Result<(), ModelBuildError> {
    if actual != expected {
        return Err(ModelBuildError::InstancePositionCountMismatch { expected, actual });
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelBuildError {
    InstancePositionCountMismatch {
        expected: usize,
        actual: usize,
    },
    CapacityOverflow,
    Topology(TopologyBuildError),
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
            Self::Topology(error) => write!(formatter, "cannot build topology: {error}"),
            Self::Hierarchy(error) => write!(formatter, "cannot build hierarchy: {error}"),
            Self::Model(error) => write!(formatter, "cannot build model: {error}"),
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

impl From<PositionError> for ModelBuildError {
    fn from(error: PositionError) -> Self {
        ModelError::from(error).into()
    }
}

impl From<ConformationError> for ModelBuildError {
    fn from(error: ConformationError) -> Self {
        ModelError::from(error).into()
    }
}

/// Failed model construction retaining its complete original builder for repair.
#[derive(Debug)]
pub struct ModelBuilderError {
    error: ModelBuildError,
    builder: Box<ModelBuilder>,
}

impl ModelBuilderError {
    pub fn error(&self) -> &ModelBuildError {
        &self.error
    }
    pub fn builder(&self) -> &ModelBuilder {
        &self.builder
    }
    pub fn into_builder(self) -> ModelBuilder {
        *self.builder
    }
}

impl fmt::Display for ModelBuilderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.error.fmt(f)
    }
}

impl std::error::Error for ModelBuilderError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.error)
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
