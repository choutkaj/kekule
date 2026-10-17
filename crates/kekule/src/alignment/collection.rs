//! Superposition and RMSD over every member of an [`Ensemble`] or frame of a
//! [`Trajectory`].

use crate::geometry::RigidTransform;
use crate::structure::{AsModelView, Ensemble, ModelView, Trajectory};
use crate::topology::{AtomSelection, Topology};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};

use super::{
    best_kabsch, best_measure, check_fit_count, check_measurement_count, AlignmentError,
    AlignmentOptions, FitAtoms, NormalizedWeights, Pairs, RigidAlignment, Weighting,
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
    fit: AlignmentOptions<'a>,
    measurement_weighting: Weighting<'a>,
}

impl<'a> AlignedRmsdOptions<'a> {
    /// Fit weighting and the periodic policy applied to both steps.
    pub const fn fit(&self) -> AlignmentOptions<'a> {
        self.fit
    }

    /// Sets [`Self::fit`].
    #[must_use]
    pub const fn with_fit(mut self, fit: AlignmentOptions<'a>) -> Self {
        self.fit = fit;
        self
    }

    /// Weights of the measured pairs, in pair order.
    pub const fn measurement_weighting(&self) -> Weighting<'a> {
        self.measurement_weighting
    }

    /// Sets [`Self::measurement_weighting`].
    #[must_use]
    pub const fn with_measurement_weighting(
        mut self,
        measurement_weighting: Weighting<'a>,
    ) -> Self {
        self.measurement_weighting = measurement_weighting;
        self
    }
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

/// Fits every item onto the reference without changing anything.
///
/// `items` are the borrowed views of one collection over `topology`, in
/// collection order. The caller applies the returned alignments.
fn fit_items<'a>(
    topology: &Topology,
    items: &[ModelView<'a>],
    reference: Reference<'a>,
    atoms: FitAtoms<'a>,
    options: AlignmentOptions<'_>,
) -> Result<SuperpositionReport, AlignmentError> {
    let reference_index = match reference {
        Reference::Index(index) => Some(index),
        Reference::View(_) => None,
    };
    let reference = reference_view(items, reference)?;
    let sets = item_pair_sets(topology, items, reference, atoms)?;
    check_fit_count(atoms.len())?;
    let weights = NormalizedWeights::new(options.weighting, atoms.len())?;
    let alignments = items
        .iter()
        .enumerate()
        .map(|(index, item)| {
            best_kabsch(*item, reference, &sets, weights, options.periodic_policy)
                .map_err(|source| item_error(index, source))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(SuperpositionReport {
        reference: reference_index,
        alignments,
    })
}

/// Direct RMSD from every item to the reference.
fn item_rmsd<'a>(
    topology: &Topology,
    items: &[ModelView<'a>],
    reference: Reference<'a>,
    atoms: FitAtoms<'a>,
    options: AlignmentOptions<'_>,
) -> Result<Quantity<Vec<f64>>, AlignmentError> {
    let reference = reference_view(items, reference)?;
    let sets = item_pair_sets(topology, items, reference, atoms)?;
    check_measurement_count(atoms.len())?;
    let weights = NormalizedWeights::new(options.weighting, atoms.len())?;
    let values = items
        .iter()
        .enumerate()
        .map(|(index, item)| {
            best_measure(
                *item,
                reference,
                &sets,
                weights,
                options.periodic_policy,
                None,
            )
            .map_err(|source| item_error(index, source))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Quantity::new(values, CANONICAL_LENGTH_UNIT))
}

/// RMSD over `measurement` pairs of every item after fitting it on `fit`.
fn item_aligned_rmsd<'a>(
    topology: &Topology,
    items: &[ModelView<'a>],
    reference: Reference<'a>,
    fit: FitAtoms<'a>,
    measurement: FitAtoms<'a>,
    options: AlignedRmsdOptions<'_>,
) -> Result<Quantity<Vec<f64>>, AlignmentError> {
    let reference = reference_view(items, reference)?;
    let fit_sets = item_pair_sets(topology, items, reference, fit)?;
    let measurement_sets = item_pair_sets(topology, items, reference, measurement)?;
    check_fit_count(fit.len())?;
    check_measurement_count(measurement.len())?;
    let fit_weights = NormalizedWeights::new(options.fit.weighting, fit.len())?;
    let measurement_weights =
        NormalizedWeights::new(options.measurement_weighting, measurement.len())?;
    let policy = options.fit.periodic_policy;
    // The best fit fixes the transform; the measurement then takes its own
    // best permutation under that transform.
    let measure = |moving: ModelView<'_>| {
        let alignment = best_kabsch(moving, reference, &fit_sets, fit_weights, policy)?;
        best_measure(
            moving,
            reference,
            &measurement_sets,
            measurement_weights,
            policy,
            Some(alignment.transform()),
        )
    };
    let values = items
        .iter()
        .enumerate()
        .map(|(index, item)| measure(*item).map_err(|source| item_error(index, source)))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Quantity::new(values, CANONICAL_LENGTH_UNIT))
}

fn reference_view<'a>(
    items: &[ModelView<'a>],
    reference: Reference<'a>,
) -> Result<ModelView<'a>, AlignmentError> {
    match reference {
        Reference::Index(index) => {
            items
                .get(index)
                .copied()
                .ok_or(AlignmentError::ReferenceOutOfRange {
                    index,
                    len: items.len(),
                })
        }
        Reference::View(view) => Ok(view),
    }
}

/// Validates the pairing once against the collection topology.
fn item_pair_sets(
    topology: &Topology,
    items: &[ModelView<'_>],
    reference: ModelView<'_>,
    atoms: FitAtoms<'_>,
) -> Result<Vec<Pairs>, AlignmentError> {
    match items.first() {
        Some(item) => atoms.pair_sets(*item, reference),
        // Without items, validate against the collection topology itself.
        None => match atoms {
            FitAtoms::Selection(selection) => {
                check_empty_selection(topology, selection, reference)?;
                Ok(vec![selection
                    .indices()
                    .iter()
                    .map(|index| (*index, *index))
                    .collect()])
            }
            FitAtoms::Correspondence(correspondence) => {
                correspondence
                    .ensure_compatible(topology, reference.topology())
                    .map_err(AlignmentError::Correspondence)?;
                Ok(vec![correspondence.index_pairs().to_vec()])
            }
            FitAtoms::Symmetric(symmetric) => {
                check_empty_selection(topology, symmetric.selection(), reference)?;
                Ok(symmetric.pair_sets())
            }
        },
    }
}

fn check_empty_selection(
    topology: &Topology,
    selection: &AtomSelection,
    reference: ModelView<'_>,
) -> Result<(), AlignmentError> {
    if !topology.shares_layout(reference.topology()) {
        return Err(AlignmentError::TopologyMismatch);
    }
    if !selection.topology().shares_layout(topology) {
        return Err(AlignmentError::SelectionTopologyMismatch);
    }
    Ok(())
}

fn transforms(report: &SuperpositionReport) -> Vec<RigidTransform> {
    report
        .alignments
        .iter()
        .map(RigidAlignment::transform)
        .collect()
}

fn item_error(index: usize, source: AlignmentError) -> AlignmentError {
    AlignmentError::Item {
        index,
        source: Box::new(source),
    }
}

impl Ensemble {
    /// Fits every member onto `reference` and applies each transform in place.
    ///
    /// `reference` is an index into this ensemble or any borrowed view
    /// (`&model`, `&other.get(0).unwrap()`, ...); `atoms` is an
    /// [`AtomSelection`](crate::topology::AtomSelection) of this layout or an
    /// [`AtomCorrespondence`](super::AtomCorrespondence) to an independent
    /// reference. Positions move and cells rotate; all other state is kept.
    /// The fit is transactional: on failure nothing changes. Coordinates are
    /// fitted as stored; no imaging or unwrapping is performed. Clone first to
    /// keep the original. `options` selects the fit weights, in pair order, and
    /// the periodic policy.
    pub fn superpose<'a>(
        &mut self,
        reference: impl Into<Reference<'a>>,
        atoms: impl Into<FitAtoms<'a>>,
        options: AlignmentOptions<'_>,
    ) -> Result<SuperpositionReport, AlignmentError> {
        let report = {
            let items = self.iter().map(|m| m.as_model_view()).collect::<Vec<_>>();
            fit_items(
                self.topology(),
                &items,
                reference.into(),
                atoms.into(),
                options,
            )?
        };
        self.transform_members(&transforms(&report))
            .map_err(AlignmentError::Transform)?;
        Ok(report)
    }

    /// Direct RMSD from every member to `reference`, measured as stored.
    ///
    /// Never fits or changes the ensemble. Values are in
    /// [`CANONICAL_LENGTH_UNIT`], one per member in order; weight them by the
    /// member weights for ensemble statistics. `options` selects the pair
    /// weights and the periodic policy.
    pub fn rmsd<'a>(
        &self,
        reference: impl Into<Reference<'a>>,
        atoms: impl Into<FitAtoms<'a>>,
        options: AlignmentOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
        let items = self.iter().map(|m| m.as_model_view()).collect::<Vec<_>>();
        item_rmsd(
            self.topology(),
            &items,
            reference.into(),
            atoms.into(),
            options,
        )
    }

    /// Fits every member on `fit` pairs and measures RMSD over `measurement`
    /// pairs, without changing the ensemble or copying coordinates.
    ///
    /// For example, fit a protein backbone and measure ligand motion. `options`
    /// holds the fit options and the measurement weights.
    pub fn aligned_rmsd<'a>(
        &self,
        reference: impl Into<Reference<'a>>,
        fit: impl Into<FitAtoms<'a>>,
        measurement: impl Into<FitAtoms<'a>>,
        options: AlignedRmsdOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
        let items = self.iter().map(|m| m.as_model_view()).collect::<Vec<_>>();
        item_aligned_rmsd(
            self.topology(),
            &items,
            reference.into(),
            fit.into(),
            measurement.into(),
            options,
        )
    }
}

impl Trajectory {
    /// Fits every frame onto `reference` and applies each transform in place.
    ///
    /// `reference` is an index into this trajectory or any borrowed view
    /// (`&model`, `&other.get(0).unwrap()`, ...); `atoms` is an
    /// [`AtomSelection`](crate::topology::AtomSelection) of this layout or an
    /// [`AtomCorrespondence`](super::AtomCorrespondence) to an independent
    /// reference. Positions move; cells, velocities, and forces rotate; all
    /// other state is kept. The fit is transactional: on failure nothing
    /// changes. Coordinates are fitted as stored; no imaging or unwrapping is
    /// performed. Clone first to keep the original.
    pub fn superpose<'a>(
        &mut self,
        reference: impl Into<Reference<'a>>,
        atoms: impl Into<FitAtoms<'a>>,
        options: AlignmentOptions<'_>,
    ) -> Result<SuperpositionReport, AlignmentError> {
        let report = {
            let items = self.iter().map(|f| f.as_model_view()).collect::<Vec<_>>();
            fit_items(
                self.topology(),
                &items,
                reference.into(),
                atoms.into(),
                options,
            )?
        };
        self.transform_frames(&transforms(&report))
            .map_err(AlignmentError::Transform)?;
        Ok(report)
    }

    /// Direct RMSD from every frame to `reference`, measured as stored.
    ///
    /// Never fits or changes the trajectory. Values are in
    /// [`CANONICAL_LENGTH_UNIT`], one per frame in temporal order. `options`
    /// selects the pair weights and the periodic policy.
    pub fn rmsd<'a>(
        &self,
        reference: impl Into<Reference<'a>>,
        atoms: impl Into<FitAtoms<'a>>,
        options: AlignmentOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
        let items = self.iter().map(|f| f.as_model_view()).collect::<Vec<_>>();
        item_rmsd(
            self.topology(),
            &items,
            reference.into(),
            atoms.into(),
            options,
        )
    }

    /// Fits every frame on `fit` pairs and measures RMSD over `measurement`
    /// pairs, without changing the trajectory or copying coordinates.
    ///
    /// For example, fit a protein backbone and measure ligand motion. `options`
    /// holds the fit options and the measurement weights.
    pub fn aligned_rmsd<'a>(
        &self,
        reference: impl Into<Reference<'a>>,
        fit: impl Into<FitAtoms<'a>>,
        measurement: impl Into<FitAtoms<'a>>,
        options: AlignedRmsdOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
        let items = self.iter().map(|f| f.as_model_view()).collect::<Vec<_>>();
        item_aligned_rmsd(
            self.topology(),
            &items,
            reference.into(),
            fit.into(),
            measurement.into(),
            options,
        )
    }
}
