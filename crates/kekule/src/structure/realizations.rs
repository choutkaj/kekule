use std::fmt;
use std::sync::Arc;

use crate::core::Molecule;
use crate::geometry::RigidTransform;
use crate::properties::OwnerProperties;
use crate::topology::transform::TopologySubsetError;
use crate::topology::{
    AtomSelection, Topology, TopologyBuildError, TopologyBuilder, TopologyPerceptionError,
};

use super::{
    AsModelView, Conformation, ConformationError, ConformationMut, Model, ModelView, Positions,
};

pub(crate) mod sealed {
    use super::*;

    /// Crate-private payload behavior behind [`super::Realization`].
    pub trait Payload: Clone {
        fn conformation_storage(&mut self) -> &mut Conformation;
        /// Projects every dense row through a subset correspondence.
        fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, ConformationError>;
        /// Applies a rigid transform to every coordinate-frame-dependent value.
        fn transformed(&self, transform: RigidTransform) -> Result<Self, ConformationError>;
    }
}

/// One realization payload stored in [`Realizations`]: a [`Conformation`]
/// plus payload-specific state.
///
/// Implemented by [`super::EnsembleMember`] and [`super::TrajectoryFrame`].
/// Payload state reads like a conformation through `Deref`.
pub trait Realization:
    sealed::Payload
    + fmt::Debug
    + PartialEq
    + From<Conformation>
    + std::ops::Deref<Target = Conformation>
{
}

/// A finite, stable-order collection of realizations of one shared topology.
///
/// Use the [`super::Ensemble`] alias for non-temporal sets (conformers,
/// alternate experimental models) and [`super::Trajectory`] when order is
/// temporal. Every item is validated against the shared topology on insertion,
/// so items can be read as [`ModelView`]s without further checks.
#[derive(Debug, Clone)]
pub struct Realizations<P> {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    items: Vec<P>,
}

impl<P: Realization> Realizations<P> {
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self {
            topology: topology.into(),
            properties: OwnerProperties::new(),
            items: Vec::new(),
        }
    }

    /// Validates and collects items in order. Fails on the first invalid item.
    pub fn from_items(
        topology: impl Into<Arc<Topology>>,
        items: impl IntoIterator<Item = P>,
    ) -> Result<Self, RealizationError> {
        let mut collection = Self::new(topology);
        for item in items {
            collection.push(item)?;
        }
        Ok(collection)
    }

    /// Collects models that share one topology layout, keeping the first
    /// model's snapshot. Model conformations move without copying.
    pub fn from_models(models: impl IntoIterator<Item = Model>) -> Result<Self, RealizationError> {
        let mut models = models.into_iter();
        let first = models.next().ok_or(RealizationError::EmptySource)?;
        let (topology, conformation) = first.into_parts();
        let mut collection = Self::new(topology);
        collection.items.push(P::from(conformation));
        for model in models {
            if !model.topology().shares_layout(&collection.topology) {
                return Err(RealizationError::TopologyMismatch);
            }
            let (_, conformation) = model.into_parts();
            collection.items.push(P::from(conformation));
        }
        Ok(collection)
    }

    /// Builds a single-molecule collection from dense positions in molecule
    /// atom order.
    pub fn from_molecule_positions(
        molecule: Molecule,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<Self, RealizationError> {
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(molecule)?;
        builder.add_instance(definition)?;
        Self::from_items(
            builder.build()?,
            positions
                .into_iter()
                .map(|positions| P::from(Conformation::new(positions))),
        )
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Collection-level annotations.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    pub fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// One topology-bound item by stable collection index.
    pub fn get(&self, index: usize) -> Option<RealizationView<'_, P>> {
        self.items.get(index).map(|item| RealizationView {
            topology: &self.topology,
            item,
        })
    }

    /// A dimension-preserving editor for one item.
    ///
    /// Bind it with `let mut item = collection.get_mut(index).unwrap()` for
    /// several edits. Whole-item replacement goes through [`Self::replace`].
    ///
    /// ```compile_fail,E0594
    /// use kekule::structure::{Ensemble, EnsembleMember, Positions};
    /// fn overwrite(ensemble: &mut Ensemble) {
    ///     *ensemble.get_mut(0).unwrap() = EnsembleMember::new(Positions::zeros(1));
    /// }
    /// ```
    pub fn get_mut(&mut self, index: usize) -> Option<RealizationMut<'_, P>> {
        self.items
            .get_mut(index)
            .map(|item| RealizationMut { item })
    }

    /// Topology-bound items in stable collection order.
    pub fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = RealizationView<'_, P>> + DoubleEndedIterator {
        let topology = &self.topology;
        self.items
            .iter()
            .map(move |item| RealizationView { topology, item })
    }

    /// Validates and appends one item, binding its bond property rows.
    pub fn push(&mut self, mut item: P) -> Result<(), RealizationError> {
        self.prepare(&mut item)?;
        self.items.push(item);
        Ok(())
    }

    /// Validates and replaces one item, returning the previous payload. The
    /// index is checked first; failure leaves the collection unchanged.
    pub fn replace(&mut self, index: usize, mut item: P) -> Result<P, RealizationError> {
        self.check_index(index)?;
        self.prepare(&mut item)?;
        Ok(std::mem::replace(&mut self.items[index], item))
    }

    /// Removes and returns one item, shifting later items down.
    pub fn remove(&mut self, index: usize) -> Result<P, RealizationError> {
        self.check_index(index)?;
        Ok(self.items.remove(index))
    }

    /// Replaces the positions of every item transactionally, in collection
    /// order. Other item state is kept.
    pub fn replace_positions(
        &mut self,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<(), RealizationError> {
        let positions = positions.into_iter().collect::<Vec<_>>();
        if positions.len() != self.items.len() {
            return Err(RealizationError::ItemCountMismatch {
                expected: self.items.len(),
                actual: positions.len(),
            });
        }
        for values in &positions {
            if values.len() != self.topology.atom_count() {
                return Err(ConformationError::AtomCountMismatch {
                    expected: self.topology.atom_count(),
                    actual: values.len(),
                }
                .into());
            }
        }
        for (item, positions) in self.items.iter_mut().zip(positions) {
            *item.conformation_storage().positions_mut() = positions;
        }
        Ok(())
    }

    /// Copies items in the requested order, sharing the topology. Ranges and
    /// strides can be written as `(start..end).step_by(stride)`; duplicates and
    /// reordering are allowed. Collection properties are kept.
    pub fn select(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, RealizationError> {
        let items = indices
            .into_iter()
            .map(|index| {
                self.check_index(index)?;
                Ok(self.items[index].clone())
            })
            .collect::<Result<Vec<_>, RealizationError>>()?;
        Ok(Self {
            topology: self.shared_topology(),
            properties: self.properties.clone(),
            items,
        })
    }

    /// Publishes one induced topology subset and projects every item onto it.
    /// Collection and item owner properties are dropped.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, RealizationError> {
        let subset = self.topology.subset(selection)?;
        let atoms = subset
            .correspondence()
            .source_atom_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let bonds = subset
            .correspondence()
            .source_bond_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let items = self
            .items
            .iter()
            .map(|item| item.project(&atoms, &bonds))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology: subset.shared_topology(),
            properties: OwnerProperties::new(),
            items,
        })
    }

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once, independent of item count,
    /// and keeps every item and the collection properties without copying.
    /// Failure leaves the collection unchanged. The new snapshot shares the
    /// original layout, so selections, buffers, and prepared potentials stay
    /// usable.
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    /// Consumes the collection, transferring topology, collection properties,
    /// and items without copying.
    pub fn into_parts(self) -> (Arc<Topology>, OwnerProperties, Vec<P>) {
        (self.topology, self.properties, self.items)
    }

    pub fn into_items(self) -> Vec<P> {
        self.items
    }

    /// Applies one transform per item transactionally.
    pub(crate) fn transform_items(
        &mut self,
        transforms: &[RigidTransform],
    ) -> Result<(), ConformationError> {
        debug_assert_eq!(transforms.len(), self.items.len());
        let items = self
            .items
            .iter()
            .zip(transforms)
            .map(|(item, transform)| item.transformed(*transform))
            .collect::<Result<Vec<_>, _>>()?;
        self.items = items;
        Ok(())
    }

    pub(crate) fn map_items<Q: Realization>(self, map: impl FnMut(P) -> Q) -> Realizations<Q> {
        Realizations {
            topology: self.topology,
            properties: self.properties,
            items: self.items.into_iter().map(map).collect(),
        }
    }

    fn check_index(&self, index: usize) -> Result<(), RealizationError> {
        if index >= self.items.len() {
            return Err(RealizationError::IndexOutOfRange {
                index,
                len: self.items.len(),
            });
        }
        Ok(())
    }

    fn prepare(&self, item: &mut P) -> Result<(), RealizationError> {
        item.validate_for(&self.topology)?;
        item.conformation_storage().bind(&self.topology);
        Ok(())
    }
}

impl<'a, P: Realization> IntoIterator for &'a Realizations<P> {
    type Item = RealizationView<'a, P>;
    type IntoIter = RealizationIter<'a, P>;

    fn into_iter(self) -> Self::IntoIter {
        RealizationIter {
            topology: &self.topology,
            items: self.items.iter(),
        }
    }
}

/// Iterator over topology-bound items of a [`Realizations`] collection.
#[derive(Debug, Clone)]
pub struct RealizationIter<'a, P> {
    topology: &'a Arc<Topology>,
    items: std::slice::Iter<'a, P>,
}

impl<'a, P> Iterator for RealizationIter<'a, P> {
    type Item = RealizationView<'a, P>;

    fn next(&mut self) -> Option<Self::Item> {
        let topology = self.topology;
        self.items
            .next()
            .map(|item| RealizationView { topology, item })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.items.size_hint()
    }
}

impl<P> ExactSizeIterator for RealizationIter<'_, P> {}

/// One topology-bound realization borrowed from a collection or buffer.
///
/// Reads payload state through `Deref` and converts to a [`ModelView`]
/// through [`AsModelView`] without copying.
#[derive(Debug)]
pub struct RealizationView<'a, P> {
    topology: &'a Arc<Topology>,
    item: &'a P,
}

impl<P> Clone for RealizationView<'_, P> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<P> Copy for RealizationView<'_, P> {}

impl<P> std::ops::Deref for RealizationView<'_, P> {
    type Target = P;

    fn deref(&self) -> &Self::Target {
        self.item
    }
}

impl<'a, P: Realization> RealizationView<'a, P> {
    /// Binds a payload that already has every row `topology` requires.
    pub fn new(topology: &'a Arc<Topology>, item: &'a P) -> Result<Self, ConformationError> {
        item.validate_for(topology)?;
        if !item.is_bound_to(topology) {
            return Err(ConformationError::BondCountMismatch {
                expected: topology.bond_count(),
                actual: item.properties().bonds().len(),
            });
        }
        Ok(Self { topology, item })
    }

    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    /// The borrowed payload with the view's full lifetime.
    pub fn payload(self) -> &'a P {
        self.item
    }

    /// Borrows the conformation as the common kernel input, with the view's
    /// full lifetime.
    pub fn as_model_view(self) -> ModelView<'a> {
        ModelView::bound(self.topology, self.item)
    }

    /// Materializes the conformation as an owned model, dropping
    /// payload-specific state.
    pub fn to_model(self) -> Model {
        self.as_model_view().to_model()
    }
}

impl<P: Realization> AsModelView for RealizationView<'_, P> {
    fn as_model_view(&self) -> ModelView<'_> {
        RealizationView::as_model_view(*self)
    }
}

/// Dimension-preserving mutable access to one item of a [`Realizations`]
/// collection.
///
/// Reads use `Deref`. There is deliberately no `DerefMut`; replace a whole
/// item through [`Realizations::replace`].
#[derive(Debug)]
pub struct RealizationMut<'a, P> {
    pub(super) item: &'a mut P,
}

impl<P> std::ops::Deref for RealizationMut<'_, P> {
    type Target = P;

    fn deref(&self) -> &Self::Target {
        self.item
    }
}

impl<P: Realization> RealizationMut<'_, P> {
    pub fn conformation_mut(&mut self) -> ConformationMut<'_> {
        ConformationMut::new(self.item.conformation_storage())
    }

    /// Applies a rigid transform transactionally: positions move; cells,
    /// velocities, and forces rotate.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), ConformationError> {
        *self.item = self.item.transformed(*transform)?;
        Ok(())
    }
}

/// Failure of a realization collection operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum RealizationError {
    /// A constructor received no items.
    EmptySource,
    IndexOutOfRange {
        index: usize,
        len: usize,
    },
    /// An input belongs to a different topology layout.
    TopologyMismatch,
    /// A collection-wide replacement has the wrong number of items.
    ItemCountMismatch {
        expected: usize,
        actual: usize,
    },
    MissingWeight {
        member: usize,
    },
    ZeroTotalWeight,
    MissingTime {
        frame: usize,
    },
    NonMonotonicTime {
        frame: usize,
    },
    Conformation(ConformationError),
    Subset(TopologySubsetError),
    TopologyBuild(TopologyBuildError),
}

impl fmt::Display for RealizationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySource => formatter.write_str("realization source is empty"),
            Self::IndexOutOfRange { index, len } => {
                write!(
                    formatter,
                    "index {index} is out of range for {len} realizations"
                )
            }
            Self::TopologyMismatch => {
                formatter.write_str("realization belongs to a different topology layout")
            }
            Self::ItemCountMismatch { expected, actual } => write!(
                formatter,
                "collection requires {expected} replacement items, but received {actual}"
            ),
            Self::MissingWeight { member } => write!(formatter, "member {member} has no weight"),
            Self::ZeroTotalWeight => {
                formatter.write_str("weights must have a positive finite total")
            }
            Self::MissingTime { frame } => write!(formatter, "frame {frame} has no time"),
            Self::NonMonotonicTime { frame } => {
                write!(formatter, "time decreases at frame {frame}")
            }
            Self::Conformation(error) => write!(formatter, "invalid realization: {error}"),
            Self::Subset(error) => write!(formatter, "cannot subset realizations: {error}"),
            Self::TopologyBuild(error) => write!(formatter, "cannot build topology: {error}"),
        }
    }
}

impl std::error::Error for RealizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Conformation(error) => Some(error),
            Self::Subset(error) => Some(error),
            Self::TopologyBuild(error) => Some(error),
            _ => None,
        }
    }
}

impl From<ConformationError> for RealizationError {
    fn from(error: ConformationError) -> Self {
        Self::Conformation(error)
    }
}

impl From<TopologySubsetError> for RealizationError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Subset(error)
    }
}

impl From<TopologyBuildError> for RealizationError {
    fn from(error: TopologyBuildError) -> Self {
        Self::TopologyBuild(error)
    }
}
