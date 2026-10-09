//! Superposition and RMSD over every member of an [`Ensemble`] or frame of a
//! [`Trajectory`].

use crate::structure::{
    AsModelView, Ensemble, ModelView, Realization, RealizationStore, Trajectory,
};
use crate::topology::TopologyAtomIndex;
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};

use super::{
    check_fit_count, check_measurement_count, kabsch_pairs, measure_pairs, AlignmentError,
    AlignmentOptions, FitAtoms, NormalizedWeights, RigidAlignment, Weighting,
};

/// The reference of a collection-wide fit or RMSD.
///
/// Convert from a `usize` item index of the collection itself, or from any
/// borrowed model view (`&model`, `&frame_view`, ...) for an independent
/// reference. Pair atoms of an independent topology with an
/// [`super::AtomCorrespondence`].
#[derive(Debug, Clone, Copy)]
pub enum Reference<'a> {
    /// An item of the collection being measured.
    Index(usize),
    /// An external view.
    View(ModelView<'a>),
}

impl From<usize> for Reference<'_> {
    fn from(index: usize) -> Self {
        Self::Index(index)
    }
}

impl<'a, T: AsModelView + ?Sized> From<&'a T> for Reference<'a> {
    fn from(view: &'a T) -> Self {
        Self::View(view.as_model_view())
    }
}

/// Options for fitting on one set of pairs and measuring another.
#[derive(Debug, Clone, Copy, Default)]
pub struct AlignedRmsdOptions<'a> {
    /// Fit weighting and the periodic policy applied to both steps.
    pub fit: AlignmentOptions<'a>,
    /// Weights of the measured pairs, in pair order.
    pub measurement_weighting: Weighting<'a>,
}

/// Complete record of one successful collection superposition.
#[derive(Debug, Clone, PartialEq)]
pub struct SuperpositionReport {
    reference: Option<usize>,
    alignments: Vec<RigidAlignment>,
}

impl SuperpositionReport {
    /// The collection index of the reference, or `None` for an external one.
    pub const fn reference_index(&self) -> Option<usize> {
        self.reference
    }

    /// One applied alignment per item, in collection order.
    pub fn alignments(&self) -> &[RigidAlignment] {
        &self.alignments
    }

    pub fn alignment(&self, index: usize) -> Option<&RigidAlignment> {
        self.alignments.get(index)
    }

    pub fn len(&self) -> usize {
        self.alignments.len()
    }

    pub fn is_empty(&self) -> bool {
        self.alignments.is_empty()
    }
}

type Pairs = Vec<(TopologyAtomIndex, TopologyAtomIndex)>;

impl<P: Realization> RealizationStore<P> {
    /// [`Self::superpose`] with explicit fit weights and periodic policy.
    pub(crate) fn superpose_with_options<'a>(
        &mut self,
        reference: impl Into<Reference<'a>>,
        atoms: impl Into<FitAtoms<'a>>,
        options: AlignmentOptions<'_>,
    ) -> Result<SuperpositionReport, AlignmentError> {
        let reference = reference.into();
        let reference_index = match reference {
            Reference::Index(index) => Some(index),
            Reference::View(_) => None,
        };
        let alignments = {
            let reference = self.reference_view(reference)?;
            let atoms = atoms.into();
            let pairs = self.pairs(reference, atoms)?;
            check_fit_count(pairs.len())?;
            let weights = NormalizedWeights::new(options.weighting, pairs.len())?;
            self.iter()
                .enumerate()
                .map(|(index, item)| {
                    kabsch_pairs(
                        item.as_model_view(),
                        reference,
                        &pairs,
                        weights,
                        options.periodic_policy,
                    )
                    .map_err(|source| item_error(index, source))
                })
                .collect::<Result<Vec<_>, _>>()?
        };
        let transforms = alignments
            .iter()
            .map(RigidAlignment::transform)
            .collect::<Vec<_>>();
        self.transform_items(&transforms)
            .map_err(AlignmentError::Transform)?;
        Ok(SuperpositionReport {
            reference: reference_index,
            alignments,
        })
    }

    /// [`Self::rmsd`] with explicit weights and periodic policy.
    pub(crate) fn rmsd_with_options<'a>(
        &self,
        reference: impl Into<Reference<'a>>,
        atoms: impl Into<FitAtoms<'a>>,
        options: AlignmentOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
        let reference = self.reference_view(reference.into())?;
        let pairs = self.pairs(reference, atoms.into())?;
        check_measurement_count(pairs.len())?;
        let weights = NormalizedWeights::new(options.weighting, pairs.len())?;
        let values = self
            .iter()
            .enumerate()
            .map(|(index, item)| {
                measure_pairs(
                    item.as_model_view(),
                    reference,
                    &pairs,
                    weights,
                    options.periodic_policy,
                    None,
                )
                .map_err(|source| item_error(index, source))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Quantity::new(values, CANONICAL_LENGTH_UNIT))
    }

    /// [`Self::aligned_rmsd`] with separate fit and measurement weights.
    pub(crate) fn aligned_rmsd_with_options<'a>(
        &self,
        reference: impl Into<Reference<'a>>,
        fit: impl Into<FitAtoms<'a>>,
        measurement: impl Into<FitAtoms<'a>>,
        options: AlignedRmsdOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
        let reference = self.reference_view(reference.into())?;
        let fit = self.pairs(reference, fit.into())?;
        let measurement = self.pairs(reference, measurement.into())?;
        check_fit_count(fit.len())?;
        check_measurement_count(measurement.len())?;
        let fit_weights = NormalizedWeights::new(options.fit.weighting, fit.len())?;
        let measurement_weights =
            NormalizedWeights::new(options.measurement_weighting, measurement.len())?;
        let policy = options.fit.periodic_policy;
        let measure = |moving: ModelView<'_>| {
            let alignment = kabsch_pairs(moving, reference, &fit, fit_weights, policy)?;
            measure_pairs(
                moving,
                reference,
                &measurement,
                measurement_weights,
                policy,
                Some(alignment.transform()),
            )
        };
        let values = self
            .iter()
            .enumerate()
            .map(|(index, item)| {
                measure(item.as_model_view()).map_err(|source| item_error(index, source))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Quantity::new(values, CANONICAL_LENGTH_UNIT))
    }

    fn reference_view<'s>(
        &'s self,
        reference: Reference<'s>,
    ) -> Result<ModelView<'s>, AlignmentError> {
        match reference {
            Reference::Index(index) => self.get(index).map(|item| item.as_model_view()).ok_or(
                AlignmentError::ReferenceOutOfRange {
                    index,
                    len: self.len(),
                },
            ),
            Reference::View(view) => Ok(view),
        }
    }

    /// Validates the pairing once against the collection topology.
    fn pairs(
        &self,
        reference: ModelView<'_>,
        atoms: FitAtoms<'_>,
    ) -> Result<Pairs, AlignmentError> {
        match self.iter().next() {
            Some(item) => atoms.pairs(item.as_model_view(), reference),
            // Without items, validate against the collection topology itself.
            None => match atoms {
                FitAtoms::Selection(selection) => {
                    if !self.topology().shares_layout(reference.topology()) {
                        return Err(AlignmentError::TopologyMismatch);
                    }
                    if !selection.topology().shares_layout(self.topology()) {
                        return Err(AlignmentError::SelectionTopologyMismatch);
                    }
                    Ok(selection
                        .indices()
                        .iter()
                        .map(|index| (*index, *index))
                        .collect())
                }
                FitAtoms::Correspondence(correspondence) => {
                    correspondence
                        .ensure_compatible(self.topology(), reference.topology())
                        .map_err(AlignmentError::Correspondence)?;
                    Ok(correspondence.index_pairs().to_vec())
                }
            },
        }
    }
}

/// Public superposition and RMSD for one collection type, documented with
/// that type's item name.
macro_rules! collection_alignment {
    ($collection:ident, $item:literal) => {
        impl $collection {
            #[doc = concat!(
                "Fits every ", $item, " onto `reference` and applies each transform in place.\n\n",
                "`reference` is an index into this collection or any borrowed view ",
                "(`&model`, `&other.get(0).unwrap()`, ...); `atoms` is an ",
                "[`AtomSelection`](crate::topology::AtomSelection) of this layout or an ",
                "[`AtomCorrespondence`](super::AtomCorrespondence) to an independent reference. ",
                "Positions move and cells rotate, as do trajectory velocities and forces; ",
                "all other state is kept. The fit is transactional: on failure nothing ",
                "changes. Coordinates are fitted as stored; no imaging or unwrapping is ",
                "performed. Clone first to keep the original."
            )]
            pub fn superpose<'a>(
                &mut self,
                reference: impl Into<Reference<'a>>,
                atoms: impl Into<FitAtoms<'a>>,
            ) -> Result<SuperpositionReport, AlignmentError> {
                self.store_mut()
                    .superpose_with_options(reference, atoms, AlignmentOptions::default())
            }

            /// [`Self::superpose`] with explicit fit weights and periodic policy.
            pub fn superpose_with_options<'a>(
                &mut self,
                reference: impl Into<Reference<'a>>,
                atoms: impl Into<FitAtoms<'a>>,
                options: AlignmentOptions<'_>,
            ) -> Result<SuperpositionReport, AlignmentError> {
                self.store_mut()
                    .superpose_with_options(reference, atoms, options)
            }

            #[doc = concat!(
                "Direct RMSD from every ", $item, " to `reference`, measured as stored.\n\n",
                "Never fits or changes the collection. Values are in ",
                "[`CANONICAL_LENGTH_UNIT`], one per ", $item, " in order."
            )]
            pub fn rmsd<'a>(
                &self,
                reference: impl Into<Reference<'a>>,
                atoms: impl Into<FitAtoms<'a>>,
            ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
                self.store()
                    .rmsd_with_options(reference, atoms, AlignmentOptions::default())
            }

            /// [`Self::rmsd`] with explicit weights and periodic policy.
            pub fn rmsd_with_options<'a>(
                &self,
                reference: impl Into<Reference<'a>>,
                atoms: impl Into<FitAtoms<'a>>,
                options: AlignmentOptions<'_>,
            ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
                self.store().rmsd_with_options(reference, atoms, options)
            }

            #[doc = concat!(
                "Fits every ", $item, " on `fit` pairs and measures RMSD over `measurement` ",
                "pairs, without changing the collection or copying coordinates.\n\n",
                "For example, fit a protein backbone and measure ligand motion."
            )]
            pub fn aligned_rmsd<'a>(
                &self,
                reference: impl Into<Reference<'a>>,
                fit: impl Into<FitAtoms<'a>>,
                measurement: impl Into<FitAtoms<'a>>,
            ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
                self.store().aligned_rmsd_with_options(
                    reference,
                    fit,
                    measurement,
                    AlignedRmsdOptions::default(),
                )
            }

            /// [`Self::aligned_rmsd`] with separate fit and measurement weights.
            pub fn aligned_rmsd_with_options<'a>(
                &self,
                reference: impl Into<Reference<'a>>,
                fit: impl Into<FitAtoms<'a>>,
                measurement: impl Into<FitAtoms<'a>>,
                options: AlignedRmsdOptions<'_>,
            ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
                self.store()
                    .aligned_rmsd_with_options(reference, fit, measurement, options)
            }
        }
    };
}

collection_alignment!(Ensemble, "member");
collection_alignment!(Trajectory, "frame");

fn item_error(index: usize, source: AlignmentError) -> AlignmentError {
    AlignmentError::Item {
        index,
        source: Box::new(source),
    }
}
