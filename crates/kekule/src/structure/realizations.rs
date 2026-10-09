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

/// One realization payload stored in an [`super::Ensemble`] or
/// [`super::Trajectory`]: a [`Conformation`] plus payload-specific state.
///
/// Implemented by [`super::EnsembleMember`] and [`super::TrajectoryFrame`].
/// Payload state reads like a conformation through `Deref`.
pub trait Realization:
    sealed::Payload + fmt::Debug + PartialEq + std::ops::Deref<Target = Conformation>
{
}

/// Crate-private storage shared by the distinct public collection types.
///
/// It owns one topology, collection owner properties, and validated payloads.
/// [`super::Ensemble`] and [`super::Trajectory`] expose it through their own
/// documented APIs; scientific behavior that differs between them, such as
/// weights or time, lives on those types.
#[derive(Debug, Clone)]
pub(crate) struct RealizationStore<P> {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    items: Vec<P>,
}

impl<P: Realization> RealizationStore<P> {
    pub(crate) fn new(topology: Arc<Topology>) -> Self {
        Self {
            topology,
            properties: OwnerProperties::new(),
            items: Vec::new(),
        }
    }

    pub(crate) fn from_items(
        topology: Arc<Topology>,
        items: impl IntoIterator<Item = P>,
    ) -> Result<Self, RealizationError> {
        let mut store = Self::new(topology);
        for item in items {
            store.push(item)?;
        }
        Ok(store)
    }

    /// Collects model conformations sharing one layout, keeping the first
    /// model's snapshot.
    pub(crate) fn from_models(
        models: impl IntoIterator<Item = Model>,
        mut payload: impl FnMut(Conformation) -> P,
    ) -> Result<Self, RealizationError> {
        let mut models = models.into_iter();
        let first = models.next().ok_or(RealizationError::EmptySource)?;
        let (topology, conformation) = first.into_parts();
        let mut store = Self::new(topology);
        store.items.push(payload(conformation));
        for model in models {
            if !model.topology().shares_layout(&store.topology) {
                return Err(RealizationError::TopologyMismatch);
            }
            let (_, conformation) = model.into_parts();
            store.items.push(payload(conformation));
        }
        Ok(store)
    }

    pub(crate) fn topology(&self) -> &Topology {
        &self.topology
    }

    pub(crate) fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub(crate) fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    pub(crate) fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }

    pub(crate) fn len(&self) -> usize {
        self.items.len()
    }

    pub(crate) fn get(&self, index: usize) -> Option<RealizationView<'_, P>> {
        self.items.get(index).map(|item| RealizationView {
            topology: &self.topology,
            item,
        })
    }

    pub(crate) fn get_mut(&mut self, index: usize) -> Option<RealizationMut<'_, P>> {
        self.items
            .get_mut(index)
            .map(|item| RealizationMut { item })
    }

    pub(crate) fn iter(&self) -> RealizationIter<'_, P> {
        RealizationIter {
            topology: &self.topology,
            items: self.items.iter(),
        }
    }

    pub(crate) fn push(&mut self, mut item: P) -> Result<(), RealizationError> {
        self.prepare(&mut item)?;
        self.items.push(item);
        Ok(())
    }

    pub(crate) fn replace(&mut self, index: usize, mut item: P) -> Result<P, RealizationError> {
        self.check_index(index)?;
        self.prepare(&mut item)?;
        Ok(std::mem::replace(&mut self.items[index], item))
    }

    pub(crate) fn remove(&mut self, index: usize) -> Result<P, RealizationError> {
        self.check_index(index)?;
        Ok(self.items.remove(index))
    }

    pub(crate) fn replace_positions(
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

    pub(crate) fn select(
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

    pub(crate) fn subset(&self, selection: &AtomSelection) -> Result<Self, RealizationError> {
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

    pub(crate) fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    pub(crate) fn into_parts(self) -> (Arc<Topology>, OwnerProperties, Vec<P>) {
        (self.topology, self.properties, self.items)
    }

    pub(crate) fn items_mut(&mut self) -> &mut [P] {
        &mut self.items
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

    pub(crate) fn map_items<Q: Realization>(self, map: impl FnMut(P) -> Q) -> RealizationStore<Q> {
        RealizationStore {
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

/// Builds a single-molecule topology for collection constructors.
pub(crate) fn single_molecule_topology(
    molecule: Molecule,
) -> Result<Arc<Topology>, RealizationError> {
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule)?;
    builder.add_instance(definition)?;
    Ok(Arc::new(builder.build()?))
}

/// The collection-independent public API of [`super::Ensemble`] and
/// [`super::Trajectory`], each documented in its own terms. `$item` names one
/// payload in prose ("member", "frame").
macro_rules! realization_collection {
    ($collection:ident, $payload:ty, $view:ty, $editor:ty, $item:literal) => {
        impl $collection {
            /// An empty collection bound to one shared topology.
            pub fn new(topology: impl Into<std::sync::Arc<$crate::topology::Topology>>) -> Self {
                Self {
                    store: $crate::structure::realizations::RealizationStore::new(topology.into()),
                }
            }

            #[doc = concat!(
                "Validates and collects each ", $item, " in order, binding its bond ",
                "property rows. Fails on the first invalid ", $item, "."
            )]
            pub fn from_items(
                topology: impl Into<std::sync::Arc<$crate::topology::Topology>>,
                items: impl IntoIterator<Item = $payload>,
            ) -> Result<Self, $crate::structure::RealizationError> {
                Ok(Self {
                    store: $crate::structure::realizations::RealizationStore::from_items(
                        topology.into(),
                        items,
                    )?,
                })
            }

            pub fn topology(&self) -> &$crate::topology::Topology {
                self.store.topology()
            }

            pub fn shared_topology(&self) -> std::sync::Arc<$crate::topology::Topology> {
                self.store.shared_topology()
            }

            /// Collection-level annotations.
            pub fn properties(&self) -> &$crate::properties::OwnerProperties {
                self.store.properties()
            }

            pub fn properties_mut(&mut self) -> &mut $crate::properties::OwnerProperties {
                self.store.properties_mut()
            }

            pub fn len(&self) -> usize {
                self.store.len()
            }

            pub fn is_empty(&self) -> bool {
                self.store.len() == 0
            }

            #[doc = concat!("One topology-bound ", $item, " by stable index.")]
            pub fn get(&self, index: usize) -> Option<$view> {
                self.store.get(index)
            }

            #[doc = concat!(
                "A dimension-preserving editor for one ", $item, ". Whole-", $item,
                " replacement goes through [`Self::replace`]."
            )]
            pub fn get_mut(&mut self, index: usize) -> Option<$editor> {
                self.store.get_mut(index)
            }

            #[doc = concat!("Topology-bound ", $item, "s in stable order.")]
            pub fn iter(&self) -> $crate::structure::RealizationIter<'_, $payload> {
                self.store.iter()
            }

            #[doc = concat!("Validates and appends one ", $item, ", binding its bond property rows.")]
            pub fn push(&mut self, item: $payload) -> Result<(), $crate::structure::RealizationError> {
                self.store.push(item)
            }

            #[doc = concat!(
                "Validates and replaces one ", $item, ", returning the previous one. ",
                "The index is checked first; failure leaves the collection unchanged."
            )]
            pub fn replace(
                &mut self,
                index: usize,
                item: $payload,
            ) -> Result<$payload, $crate::structure::RealizationError> {
                self.store.replace(index, item)
            }

            #[doc = concat!("Removes and returns one ", $item, ", shifting later ones down.")]
            pub fn remove(
                &mut self,
                index: usize,
            ) -> Result<$payload, $crate::structure::RealizationError> {
                self.store.remove(index)
            }

            #[doc = concat!(
                "Replaces the positions of every ", $item, " transactionally, in order. ",
                "Other state is kept."
            )]
            pub fn replace_positions(
                &mut self,
                positions: impl IntoIterator<Item = $crate::structure::Positions>,
            ) -> Result<(), $crate::structure::RealizationError> {
                self.store.replace_positions(positions)
            }

            #[doc = concat!(
                "Copies ", $item, "s in the requested order, sharing the topology. ",
                "Ranges and strides can be written as `(start..end).step_by(stride)`; ",
                "duplicates and reordering are allowed. Collection properties are kept."
            )]
            pub fn select(
                &self,
                indices: impl IntoIterator<Item = usize>,
            ) -> Result<Self, $crate::structure::RealizationError> {
                Ok(Self {
                    store: self.store.select(indices)?,
                })
            }

            #[doc = concat!(
                "Publishes one induced topology subset and projects every ", $item,
                " onto it. Collection and per-", $item, " owner properties are dropped."
            )]
            pub fn subset(
                &self,
                selection: &$crate::topology::AtomSelection,
            ) -> Result<Self, $crate::structure::RealizationError> {
                Ok(Self {
                    store: self.store.subset(selection)?,
                })
            }

            #[doc = concat!(
                "Installs default perception through one new shared topology snapshot.\n\n",
                "Delegates to [`Topology::perceived`](crate::topology::Topology::perceived) ",
                "once and keeps every ", $item, " and the collection properties without ",
                "copying. Failure leaves the collection unchanged. The new snapshot shares ",
                "the original layout, so selections, buffers, and prepared potentials stay usable."
            )]
            pub fn perceive(&mut self) -> Result<(), $crate::topology::TopologyPerceptionError> {
                self.store.perceive()
            }

            #[doc = concat!(
                "Consumes the collection, transferring topology, collection properties, ",
                "and every ", $item, " without copying."
            )]
            pub fn into_parts(
                self,
            ) -> (
                std::sync::Arc<$crate::topology::Topology>,
                $crate::properties::OwnerProperties,
                Vec<$payload>,
            ) {
                self.store.into_parts()
            }

            pub fn into_items(self) -> Vec<$payload> {
                self.store.into_parts().2
            }
        }

        impl<'a> IntoIterator for &'a $collection {
            type Item = $crate::structure::RealizationView<'a, $payload>;
            type IntoIter = $crate::structure::RealizationIter<'a, $payload>;

            fn into_iter(self) -> Self::IntoIter {
                self.store.iter()
            }
        }
    };
}
pub(crate) use realization_collection;

/// Iterator over topology-bound items of an [`super::Ensemble`] or
/// [`super::Trajectory`].
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

/// Dimension-preserving mutable access to one item of an [`super::Ensemble`]
/// or [`super::Trajectory`].
///
/// Reads use `Deref`. There is deliberately no `DerefMut`; replace a whole
/// item through the collection's `replace`.
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

/// Failure of an [`super::Ensemble`] or [`super::Trajectory`] operation.
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
