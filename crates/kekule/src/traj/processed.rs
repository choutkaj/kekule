//! Samples processed as they are read: [`Processed`] and its
//! [`ProcessingSteps`].

use std::fmt;
use std::sync::{Arc, OnceLock};

use crate::alignment::{self, AlignmentOptions};
use crate::geometry::{Matrix3, PeriodicCell, Point3, RigidTransform, Vector3};
use crate::structure::{
    subset_rows, Distribution, Model, ModelError, Order, Positions, Sample, SampleRead,
    SampleReadError, SampleView, StructureError,
};
use crate::topology::{AtomSelection, Topology, TopologySubsetError};
use crate::units::Quantity;

use super::periodic::MoleculeImager;

/// The steps a [`Processed`] reader runs on every sample as it is read.
/// Each is optional; they run in this order:
///
/// 1. Wrap ([`Self::with_wrap`]): every molecule made whole and imaged
///    around anchor atoms, as [`super::periodic::image_molecules`] does; or
///    make whole ([`Self::with_whole`]), as
///    [`super::periodic::make_molecules_whole`] does. Samples without a cell
///    skip this step.
/// 2. Center ([`Self::with_center`]): the sample translated so that the
///    centroid of the centered atoms lies at the center of its cell, or,
///    without a cell, where it lies in the reader's sample 0 as stored.
/// 3. Fit ([`Self::with_fit`]): the sample superposed onto a reference
///    through the fitted atoms, as
///    [`FrameSuperposer`](super::analysis::FrameSuperposer) does with
///    default options; positions move, and the cell, velocities, and forces
///    rotate.
/// 4. Smooth ([`Self::with_smoothing`]): positions averaged over a window of
///    samples centered on the sample, fewer at the ends, each processed by
///    steps 1–3 and, with a cell, imaged onto the sample first, as
///    [`MoleculeImager::image_onto`] does, so that molecules crossing the
///    cell between samples are averaged whole. Everything but the positions
///    is the sample's own.
/// 5. Subset ([`Self::with_subset`]): the sample projected onto the atoms of
///    a topology subset, which becomes the reader's topology.
///
/// ```
/// use kekule::topology::AtomSelection;
/// use kekule::traj::{FitReference, ProcessingSteps};
///
/// let topology = kekule::smiles::to_topology("CCO")?;
/// let all = AtomSelection::all(topology);
/// let steps = ProcessingSteps::default()
///     .with_whole(true)
///     .with_fit(all, FitReference::Sample(0))
///     .with_smoothing(5)?;
/// assert_eq!(steps.smoothing(), Some(5));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Default)]
pub struct ProcessingSteps {
    whole: bool,
    wrap: Option<AtomSelection>,
    center: Option<AtomSelection>,
    fit: Option<(AtomSelection, FitReference)>,
    smoothing: Option<usize>,
    subset: Option<AtomSelection>,
}

/// What [`ProcessingSteps::with_fit`] superposes the samples onto.
#[derive(Debug, Clone, PartialEq)]
pub enum FitReference {
    /// Sample `n` of the same reader, after the steps before fitting.
    Sample(usize),
    /// These positions, one per atom of the reader's topology.
    Positions(Positions),
}

/// One of the [`ProcessingSteps`], which [`ProcessingError`] names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ProcessingStep {
    /// Wrapping around anchor atoms.
    Wrap,
    /// Centering on atoms.
    Center,
    /// Fitting onto a reference.
    Fit,
    /// Projecting onto an atom subset.
    Subset,
}

impl fmt::Display for ProcessingStep {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Wrap => "wrapping",
            Self::Center => "centering",
            Self::Fit => "fitting",
            Self::Subset => "the subset",
        })
    }
}

impl ProcessingSteps {
    /// Makes every molecule whole across the cell's boundaries; wrapping
    /// does so too.
    #[must_use]
    pub fn with_whole(mut self, whole: bool) -> Self {
        self.whole = whole;
        self
    }

    /// Images every molecule whole around the molecules of `anchors`, whose
    /// combined center moves to the center of the cell.
    #[must_use]
    pub fn with_wrap(mut self, anchors: AtomSelection) -> Self {
        self.wrap = Some(anchors);
        self
    }

    /// Centers every sample on the centroid of `atoms`.
    #[must_use]
    pub fn with_center(mut self, atoms: AtomSelection) -> Self {
        self.center = Some(atoms);
        self
    }

    /// Superposes every sample onto `reference` through `atoms`.
    #[must_use]
    pub fn with_fit(mut self, atoms: AtomSelection, reference: FitReference) -> Self {
        self.fit = Some((atoms, reference));
        self
    }

    /// Averages positions over `window` samples centered on each one.
    ///
    /// # Errors
    ///
    /// [`ProcessingError::InvalidWindow`] when `window` is not odd.
    pub fn with_smoothing(mut self, window: usize) -> Result<Self, ProcessingError> {
        if window.is_multiple_of(2) {
            return Err(ProcessingError::InvalidWindow { window });
        }
        self.smoothing = Some(window);
        Ok(self)
    }

    /// Projects every sample onto the atoms of `atoms`, after the other
    /// steps.
    #[must_use]
    pub fn with_subset(mut self, atoms: AtomSelection) -> Self {
        self.subset = Some(atoms);
        self
    }

    /// The same steps without smoothing.
    #[must_use]
    pub fn without_smoothing(mut self) -> Self {
        self.smoothing = None;
        self
    }

    /// Whether molecules are made whole without wrapping.
    pub const fn is_whole(&self) -> bool {
        self.whole
    }

    /// The anchor atoms of wrapping.
    pub const fn wrap(&self) -> Option<&AtomSelection> {
        self.wrap.as_ref()
    }

    /// The centered atoms.
    pub const fn center(&self) -> Option<&AtomSelection> {
        self.center.as_ref()
    }

    /// The fitted atoms and their reference.
    pub fn fit(&self) -> Option<(&AtomSelection, &FitReference)> {
        self.fit
            .as_ref()
            .map(|(atoms, reference)| (atoms, reference))
    }

    /// The smoothing window.
    pub const fn smoothing(&self) -> Option<usize> {
        self.smoothing
    }

    /// The atoms of the subset.
    pub const fn subset(&self) -> Option<&AtomSelection> {
        self.subset.as_ref()
    }

    /// Whether there is no step, so samples are read as stored.
    pub const fn is_empty(&self) -> bool {
        !self.whole
            && self.wrap.is_none()
            && self.center.is_none()
            && self.fit.is_none()
            && self.smoothing.is_none()
            && self.subset.is_none()
    }
}

/// Failure to set up a [`Processed`] reader.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ProcessingError {
    /// A smoothing window is an odd number of samples.
    InvalidWindow {
        /// The rejected window.
        window: usize,
    },
    /// The atoms of a step belong to another topology than the reader's.
    SelectionTopologyMismatch {
        /// The step.
        step: ProcessingStep,
    },
    /// The atoms of a step are none.
    EmptySelection {
        /// The step.
        step: ProcessingStep,
    },
    /// The reference sample does not exist.
    ReferenceOutOfRange {
        /// The reference sample index.
        index: usize,
        /// The number of samples.
        len: usize,
    },
    /// The reference positions do not fit the reader's topology; see the
    /// source.
    Reference(Box<ModelError>),
    /// The subset cannot be made; see the source.
    Subset(TopologySubsetError),
    /// A step cannot process the first sample or the reference; see the
    /// source.
    Sample(SampleReadError),
}

impl fmt::Display for ProcessingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidWindow { window } => write!(
                formatter,
                "A smoothing window is an odd number of samples, not {window}."
            ),
            Self::SelectionTopologyMismatch { step } => {
                write!(formatter, "The atoms of {step} belong to another topology.")
            }
            Self::EmptySelection { step } => {
                write!(formatter, "The atoms of {step} are none.")
            }
            Self::ReferenceOutOfRange { index, len } => write!(
                formatter,
                "The reference sample {index} does not exist. There are {}.",
                crate::diagnostics::counted(*len, "sample", "samples")
            ),
            Self::Reference(_) => {
                formatter.write_str("The reference positions do not fit the topology.")
            }
            Self::Subset(error) => fmt::Display::fmt(error, formatter),
            Self::Sample(error) => fmt::Display::fmt(error, formatter),
        }
    }
}

impl std::error::Error for ProcessingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Reference(error) => Some(error.as_ref()),
            Self::Subset(error) => std::error::Error::source(error),
            Self::Sample(error) => std::error::Error::source(error),
            _ => None,
        }
    }
}

impl From<SampleReadError> for ProcessingError {
    fn from(error: SampleReadError) -> Self {
        Self::Sample(error)
    }
}

/// The samples of another reader with [`ProcessingSteps`] run on each as it
/// is read.
///
/// Reading sample `i` reads sample `i` of the reader below, and with
/// smoothing the samples of its window, and processes them; nothing is
/// kept, so a host caches samples in front of this reader when it reads
/// them again. The order and the distribution are those of the reader
/// below. [`Self::unsmoothed`] reads the samples without smoothing, which is
/// what analyses read: smoothing is for display. To cache the samples that
/// smoothing averages, a host puts its cache in front of a reader with the
/// other steps and smooths through another `Processed` reader in front of
/// the cache.
///
/// The reader below can be owned, borrowed (`&R`), or shared (`Arc<R>`),
/// such as an `Arc<dyn SampleRead>`.
///
/// ```
/// use kekule::geometry::Point3;
/// use kekule::structure::{Positions, Sample, SampleRead, SampleSet};
/// use kekule::topology::AtomSelection;
/// use kekule::traj::{FitReference, Processed, ProcessingSteps};
/// use kekule::units::{Quantity, NANOMETER};
///
/// let topology = kekule::smiles::to_topology("CCO")?;
/// let mut set = SampleSet::new(topology.clone());
/// for x in [0.0, 1.0] {
///     let points = [Point3::new(x, 0.0, 0.0), Point3::new(x + 0.1, 0.0, 0.0), Point3::new(x, 0.1, 0.0)];
///     let positions = Positions::new(Quantity::new(points, NANOMETER))?;
///     set.push(Sample::new(positions), None)?;
/// }
/// let all = AtomSelection::all(topology);
/// let steps = ProcessingSteps::default().with_fit(all, FitReference::Sample(0));
/// let fitted = Processed::new(&set, steps)?;
/// let sample = fitted.read(1)?;
/// let x = sample.positions().values().value()[0].x;
/// assert!(x.abs() < 1.0e-12);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct Processed<R> {
    inner: R,
    steps: ProcessingSteps,
    /// The topology of the processed samples: the subset's, or the inner
    /// reader's.
    topology: Arc<Topology>,
    /// Plans molecules for making whole, wrapping, and smoothing.
    imager: Option<MoleculeImager>,
    /// The molecules wrapping images around.
    anchors: Option<Vec<bool>>,
    /// Where centering puts the centered atoms in samples without a cell:
    /// their centroid in sample 0 of the reader below.
    rest: OnceLock<Vector3>,
    /// The structure fitting superposes onto.
    reference: Option<Model>,
    /// The dense atom and bond rows the subset keeps.
    subset: Option<(Vec<usize>, Vec<usize>)>,
}

impl<R: SampleRead> Processed<R> {
    /// The samples of `inner` processed by `steps`. The steps are checked
    /// against the reader's topology, the reference made, and sample 0
    /// processed, so steps that cannot process the samples fail here.
    ///
    /// # Errors
    ///
    /// [`ProcessingError`] for atoms of another topology or no atoms, a
    /// reference sample past the last or reference positions of another
    /// atom count, a subset that cannot be made, and a failure to read or
    /// process the reference or sample 0.
    pub fn new(inner: R, steps: ProcessingSteps) -> Result<Self, ProcessingError> {
        let source = Arc::clone(inner.topology());
        for (step, atoms) in [
            (ProcessingStep::Wrap, steps.wrap.as_ref()),
            (ProcessingStep::Center, steps.center.as_ref()),
            (
                ProcessingStep::Fit,
                steps.fit.as_ref().map(|(atoms, _)| atoms),
            ),
            (ProcessingStep::Subset, steps.subset.as_ref()),
        ] {
            if let Some(atoms) = atoms {
                check_atoms(&source, atoms, step)?;
            }
        }
        let imager = (steps.whole || steps.wrap.is_some() || steps.smoothing.is_some())
            .then(|| MoleculeImager::new(Arc::clone(&source)));
        let anchors = match (&imager, &steps.wrap) {
            (Some(imager), Some(anchors)) => Some(
                imager
                    .anchor_groups(anchors)
                    .expect("the anchors were checked"),
            ),
            _ => None,
        };
        let (topology, subset) = match &steps.subset {
            Some(atoms) => {
                let subset = Topology::subset(Arc::clone(&source), atoms)
                    .map_err(ProcessingError::Subset)?;
                let rows = subset_rows(&subset);
                (subset.shared_topology(), Some(rows))
            }
            None => (Arc::clone(&source), None),
        };
        let mut processed = Self {
            inner,
            steps,
            topology,
            imager,
            anchors,
            rest: OnceLock::new(),
            reference: None,
            subset,
        };

        if let Some((_, reference)) = &processed.steps.fit {
            processed.reference = Some(processed.reference_model(reference)?);
        }
        if !processed.inner.is_empty() {
            processed.read(0)?;
        }
        Ok(processed)
    }

    /// The reader below.
    pub const fn inner(&self) -> &R {
        &self.inner
    }

    /// The steps.
    pub const fn steps(&self) -> &ProcessingSteps {
        &self.steps
    }

    /// The reader below and the steps.
    pub fn into_parts(self) -> (R, ProcessingSteps) {
        (self.inner, self.steps)
    }

    /// The samples without smoothing: every other step runs, the subset
    /// included. Without smoothing, they are this reader's samples.
    pub const fn unsmoothed(&self) -> Unsmoothed<'_, R> {
        Unsmoothed { processed: self }
    }

    /// The model fitting superposes onto.
    fn reference_model(&self, reference: &FitReference) -> Result<Model, ProcessingError> {
        match reference {
            FitReference::Sample(index) => {
                let len = self.inner.len();
                if *index >= len {
                    return Err(ProcessingError::ReferenceOutOfRange { index: *index, len });
                }
                let mut sample = self.read_inner(*index)?;
                self.prepare(*index, &mut sample)?;
                Ok(self.view(*index, &sample).as_model_view().to_model())
            }
            FitReference::Positions(positions) => {
                Model::new(Arc::clone(self.inner.topology()), positions.clone())
                    .map_err(|error| ProcessingError::Reference(Box::new(error)))
            }
        }
    }

    /// Sample `index` of the reader below, checked to be bound to its
    /// topology, as every step expects.
    fn read_inner(&self, index: usize) -> Result<Sample, SampleReadError> {
        let sample = self.inner.read(index)?;
        SampleView::new(&self.inner, index, &sample)?;
        Ok(sample)
    }

    /// A view of `sample`, sample `index` of the reader below.
    fn view<'a>(&'a self, index: usize, sample: &'a Sample) -> SampleView<'a> {
        SampleView::of_reader(&self.inner, index, sample)
    }

    /// Wraps or makes whole, then centers, sample `index` of the reader
    /// below.
    fn prepare(&self, index: usize, sample: &mut Sample) -> Result<(), SampleReadError> {
        if let Some(imager) = &self.imager {
            if sample.cell().is_some() && (self.steps.whole || self.anchors.is_some()) {
                let positions = imager
                    .frame_positions(
                        index,
                        self.view(index, sample).as_model_view(),
                        self.anchors.as_deref(),
                    )
                    .map_err(|source| SampleReadError::Periodic {
                        sample: index,
                        source: Box::new(source),
                    })?;
                sample.replace_checked_positions(positions);
            }
        }
        if let Some(center) = &self.steps.center {
            let target = match sample.cell() {
                Some(cell) => cell_center(cell),
                None => self.rest(center)?,
            };
            let shift = target - centroid(sample.positions().values().value(), center);
            let structure = |source| SampleReadError::Structure {
                sample: index,
                source: Box::new(source),
            };
            let transform = RigidTransform::from_canonical(Matrix3::identity(), shift)
                .map_err(|_| structure(StructureError::NonFiniteCoordinate))?;
            sample
                .apply_rigid_transform(&transform)
                .map_err(structure)?;
        }
        Ok(())
    }

    /// The centroid of the centered atoms in sample 0 of the reader below.
    fn rest(&self, center: &AtomSelection) -> Result<Vector3, SampleReadError> {
        if let Some(rest) = self.rest.get() {
            return Ok(*rest);
        }
        let first = self.read_inner(0)?;
        let rest = centroid(first.positions().values().value(), center);
        Ok(*self.rest.get_or_init(|| rest))
    }

    /// Sample `index` of the reader below after steps 1–3, before the
    /// subset.
    fn stepped(&self, index: usize) -> Result<Sample, SampleReadError> {
        let mut sample = self.read_inner(index)?;
        self.prepare(index, &mut sample)?;
        if let (Some((atoms, _)), Some(reference)) = (&self.steps.fit, &self.reference) {
            let fitted = alignment::fit(
                self.view(index, &sample).as_model_view(),
                reference.as_model_view(),
                atoms,
                AlignmentOptions::default(),
            )
            .map_err(|source| SampleReadError::Alignment {
                sample: index,
                source: Box::new(source),
            })?;
            sample
                .apply_rigid_transform(&fitted.transform())
                .map_err(|source| SampleReadError::Structure {
                    sample: index,
                    source: Box::new(source),
                })?;
        }
        Ok(sample)
    }

    /// Sample `index` after steps 1–3, averaged over the window of `half`
    /// samples on each side.
    fn smoothed(&self, index: usize, half: usize) -> Result<Sample, SampleReadError> {
        let mut center = self.stepped(index)?;
        let last = self.inner.len().saturating_sub(1);
        let window = index.saturating_sub(half)..=(index + half).min(last);
        if window.start() == window.end() {
            return Ok(center);
        }
        let mut sum: Vec<Vector3> = center
            .positions()
            .values()
            .value()
            .iter()
            .map(|point| *point - Point3::origin())
            .collect();
        let mut count = 1.0;
        let imager = self.imager.as_ref().expect("smoothing plans molecules");
        for other in window {
            if other == index {
                continue;
            }
            let sample = self.stepped(other)?;
            if sample.cell().is_some() {
                let imaged = imager
                    .reference_positions(
                        other,
                        self.view(other, &sample).as_model_view(),
                        self.view(index, &center).as_model_view(),
                    )
                    .map_err(|source| SampleReadError::Periodic {
                        sample: other,
                        source: Box::new(source),
                    })?;
                add(&mut sum, imaged.values().value());
            } else {
                add(&mut sum, sample.positions().values().value());
            }
            count += 1.0;
        }
        let average: Vec<Point3> = sum.iter().map(|v| Point3::origin() + *v / count).collect();
        let positions = Positions::from_vec(Quantity::new(average, crate::units::NANOMETER))
            .map_err(|source| SampleReadError::Structure {
                sample: index,
                source: Box::new(StructureError::Array(source)),
            })?;
        center.replace_checked_positions(positions);
        Ok(center)
    }

    /// Projects sample `index` onto the subset.
    fn finish(&self, index: usize, sample: Sample) -> Result<Sample, SampleReadError> {
        match &self.subset {
            Some((atoms, bonds)) => {
                sample
                    .projected(atoms, bonds)
                    .map_err(|source| SampleReadError::Structure {
                        sample: index,
                        source: Box::new(source),
                    })
            }
            None => Ok(sample),
        }
    }
}

impl<R: SampleRead> SampleRead for Processed<R> {
    fn topology(&self) -> &Arc<Topology> {
        &self.topology
    }

    fn len(&self) -> usize {
        self.inner.len()
    }

    fn order(&self) -> Option<&Order> {
        self.inner.order()
    }

    fn distribution(&self) -> Option<&Distribution> {
        self.inner.distribution()
    }

    fn peek(&self, _index: usize) -> Option<SampleView<'_>> {
        None
    }

    fn read(&self, index: usize) -> Result<Sample, SampleReadError> {
        let sample = match self.steps.smoothing {
            Some(window) => self.smoothed(index, window / 2)?,
            None => self.stepped(index)?,
        };
        self.finish(index, sample)
    }
}

/// The samples of a [`Processed`] reader without smoothing, from
/// [`Processed::unsmoothed`].
#[derive(Debug)]
pub struct Unsmoothed<'a, R> {
    processed: &'a Processed<R>,
}

impl<R> Clone for Unsmoothed<'_, R> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<R> Copy for Unsmoothed<'_, R> {}

impl<R: SampleRead> SampleRead for Unsmoothed<'_, R> {
    fn topology(&self) -> &Arc<Topology> {
        self.processed.topology()
    }

    fn len(&self) -> usize {
        self.processed.len()
    }

    fn order(&self) -> Option<&Order> {
        self.processed.order()
    }

    fn distribution(&self) -> Option<&Distribution> {
        self.processed.distribution()
    }

    fn peek(&self, _index: usize) -> Option<SampleView<'_>> {
        None
    }

    fn read(&self, index: usize) -> Result<Sample, SampleReadError> {
        let sample = self.processed.stepped(index)?;
        self.processed.finish(index, sample)
    }
}

/// Checks that the atoms of `step` are some atoms of `topology`.
fn check_atoms(
    topology: &Topology,
    atoms: &AtomSelection,
    step: ProcessingStep,
) -> Result<(), ProcessingError> {
    if !topology.shares_layout(atoms.topology()) {
        return Err(ProcessingError::SelectionTopologyMismatch { step });
    }
    if atoms.indices().is_empty() {
        return Err(ProcessingError::EmptySelection { step });
    }
    Ok(())
}

/// The centroid of the `atoms` of `points`, from the origin.
fn centroid(points: &[Point3], atoms: &AtomSelection) -> Vector3 {
    let indices = atoms.indices();
    let mut sum = Vector3::zero();
    for index in indices {
        sum += points[index.index()] - Point3::origin();
    }
    sum / indices.len() as f64
}

/// The center of a cell, from the origin.
fn cell_center(cell: &PeriodicCell) -> Vector3 {
    let [a, b, c] = cell.vectors().into_value();
    (a + b + c) * 0.5
}

/// Adds `points` to `sum`, atom by atom.
fn add(sum: &mut [Vector3], points: &[Point3]) {
    for (total, point) in sum.iter_mut().zip(points) {
        *total += *point - Point3::origin();
    }
}
