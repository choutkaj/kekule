//! Superposition and RMSD over every geometry of a [`StructureSet`], and RMSD
//! over listed geometries of any [`StructureSetRead`], through one implementation.

use crate::math::RigidTransform;
use crate::parallel::{Progress, Run};
use crate::structure::{
    visit_geometry, AsStructureView, Geometry, GeometryView, StructureSet, StructureSetRead,
    StructureView,
};
use crate::topology::{AtomSelection, Topology};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};

use super::{
    best_kabsch, best_measure, check_fit_count, check_measurement_count, AlignmentError,
    AlignmentOptions, FitAtoms, NormalizedWeights, Pairs, PeriodicPolicy, RigidAlignment,
    Weighting,
};

/// The reference of a collection-wide fit or RMSD.
///
/// Convert from a `usize` item index of the collection itself, or from any
/// borrowed structure view (`&structure`, `&frame_view`, ...) for an independent
/// reference. Pair atoms of an independent topology with an
/// [`super::AtomCorrespondence`].
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum Reference<'a> {
    /// An item of the collection being measured.
    Index(usize),
    /// An external view.
    View(StructureView<'a>),
}

impl From<usize> for Reference<'_> {
    fn from(index: usize) -> Self {
        Self::Index(index)
    }
}

impl<'a, T: AsStructureView + ?Sized> From<&'a T> for Reference<'a> {
    fn from(view: &'a T) -> Self {
        Self::View(view.as_structure_view())
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

    /// The alignment applied to item `index`, or `None` when `index` is out of
    /// range.
    pub fn alignment(&self, index: usize) -> Option<&RigidAlignment> {
        self.alignments.get(index)
    }

    /// Number of aligned items.
    pub fn len(&self) -> usize {
        self.alignments.len()
    }

    /// Whether the collection had no items.
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
    items: &[StructureView<'a>],
    reference: Reference<'a>,
    atoms: FitAtoms<'a>,
    options: AlignmentOptions<'_>,
) -> Result<SuperpositionReport, AlignmentError> {
    let reference_index = match reference {
        Reference::Index(index) => Some(index),
        Reference::View(_) => None,
    };
    let reference = reference_view(items, reference)?;
    let sets = item_pair_sets(topology, items.first().copied(), reference, atoms)?;
    check_fit_count(atoms.pair_count())?;
    let weights = NormalizedWeights::new(options.weighting, atoms.pair_count())?;
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

/// What an RMSD measures: direct RMSD over one pairing, or RMSD over one
/// pairing after a fit over another.
#[derive(Clone, Copy)]
enum Measured<'a, 'w> {
    Direct(FitAtoms<'a>, AlignmentOptions<'w>),
    Aligned(FitAtoms<'a>, FitAtoms<'a>, AlignedRmsdOptions<'w>),
}

/// An RMSD to one reference, validated once and measured for many items.
enum Measurement<'r, 'w> {
    Direct {
        reference: StructureView<'r>,
        sets: Vec<Pairs>,
        weights: NormalizedWeights<'w>,
        policy: PeriodicPolicy,
    },
    Aligned {
        reference: StructureView<'r>,
        fit_sets: Vec<Pairs>,
        fit_weights: NormalizedWeights<'w>,
        measurement_sets: Vec<Pairs>,
        measurement_weights: NormalizedWeights<'w>,
        policy: PeriodicPolicy,
    },
}

impl<'r, 'w> Measurement<'r, 'w> {
    /// Validates the pairing against the first item, or the collection's
    /// topology without items, and the weights.
    fn new(
        topology: &Topology,
        first: Option<StructureView<'_>>,
        reference: StructureView<'r>,
        measured: Measured<'_, 'w>,
    ) -> Result<Self, AlignmentError> {
        match measured {
            Measured::Direct(atoms, options) => {
                let sets = item_pair_sets(topology, first, reference, atoms)?;
                check_measurement_count(atoms.pair_count())?;
                Ok(Self::Direct {
                    reference,
                    sets,
                    weights: NormalizedWeights::new(options.weighting, atoms.pair_count())?,
                    policy: options.periodic_policy,
                })
            }
            Measured::Aligned(fit, measurement, options) => {
                let fit_sets = item_pair_sets(topology, first, reference, fit)?;
                let measurement_sets = item_pair_sets(topology, first, reference, measurement)?;
                check_fit_count(fit.pair_count())?;
                check_measurement_count(measurement.pair_count())?;
                Ok(Self::Aligned {
                    reference,
                    fit_sets,
                    fit_weights: NormalizedWeights::new(options.fit.weighting, fit.pair_count())?,
                    measurement_sets,
                    measurement_weights: NormalizedWeights::new(
                        options.measurement_weighting,
                        measurement.pair_count(),
                    )?,
                    policy: options.fit.periodic_policy,
                })
            }
        }
    }

    /// The RMSD of `moving` to the reference.
    fn measure(&self, moving: StructureView<'_>) -> Result<f64, AlignmentError> {
        match self {
            Self::Direct {
                reference,
                sets,
                weights,
                policy,
            } => best_measure(moving, *reference, sets, *weights, *policy, None),
            // The best fit fixes the transform; the measurement then takes
            // its own best permutation under that transform.
            Self::Aligned {
                reference,
                fit_sets,
                fit_weights,
                measurement_sets,
                measurement_weights,
                policy,
            } => {
                let alignment = best_kabsch(moving, *reference, fit_sets, *fit_weights, *policy)?;
                best_measure(
                    moving,
                    *reference,
                    measurement_sets,
                    *measurement_weights,
                    *policy,
                    Some(alignment.transform()),
                )
            }
        }
    }
}

/// The RMSD of every item to the reference.
fn item_values<'a>(
    topology: &Topology,
    items: &[StructureView<'a>],
    reference: Reference<'a>,
    measured: Measured<'_, '_>,
) -> Result<Quantity<Vec<f64>>, AlignmentError> {
    let reference = reference_view(items, reference)?;
    let measurement = Measurement::new(topology, items.first().copied(), reference, measured)?;
    let values = items
        .iter()
        .enumerate()
        .map(|(index, item)| {
            measurement
                .measure(*item)
                .map_err(|source| item_error(index, source))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Quantity::new(values, CANONICAL_LENGTH_UNIT))
}

/// Direct RMSD from every item to the reference.
fn item_rmsd<'a>(
    topology: &Topology,
    items: &[StructureView<'a>],
    reference: Reference<'a>,
    atoms: FitAtoms<'a>,
    options: AlignmentOptions<'_>,
) -> Result<Quantity<Vec<f64>>, AlignmentError> {
    item_values(topology, items, reference, Measured::Direct(atoms, options))
}

/// RMSD over `measurement` pairs of every item after fitting it on `fit`.
fn item_aligned_rmsd<'a>(
    topology: &Topology,
    items: &[StructureView<'a>],
    reference: Reference<'a>,
    fit: FitAtoms<'a>,
    measurement: FitAtoms<'a>,
    options: AlignedRmsdOptions<'_>,
) -> Result<Quantity<Vec<f64>>, AlignmentError> {
    item_values(
        topology,
        items,
        reference,
        Measured::Aligned(fit, measurement, options),
    )
}

/// The RMSD of each of `geometries` of `reader` to the reference, each lent
/// or read in turn, polling `run`'s cancellation before each and reporting
/// the geometries measured.
fn geometry_values(
    reader: &dyn StructureSetRead,
    geometries: impl IntoIterator<Item = usize>,
    reference: Reference<'_>,
    measured: Measured<'_, '_>,
    run: Run<'_>,
) -> Result<Quantity<Vec<f64>>, AlignmentError> {
    let read_error = |error| AlignmentError::Read(Box::new(error));
    let held: Geometry;
    let reference = match reference {
        Reference::View(view) => view,
        Reference::Index(index) => {
            let len = reader.len();
            if index >= len {
                return Err(AlignmentError::ReferenceOutOfRange { index, len });
            }
            match reader.peek(index) {
                Some(view) => view.as_structure_view(),
                None => {
                    held = reader.read(index).map_err(read_error)?;
                    GeometryView::new(reader, index, &held)
                        .map_err(read_error)?
                        .as_structure_view()
                }
            }
        }
    };
    let topology = reader.topology();
    let geometries = geometries.into_iter().collect::<Vec<_>>();
    if geometries.is_empty() {
        run.check().map_err(AlignmentError::Cancelled)?;
        Measurement::new(topology, None, reference, measured)?;
        run.report(Progress::of(0, 0));
        return Ok(Quantity::new(Vec::new(), CANONICAL_LENGTH_UNIT));
    }
    let mut measurement = None;
    let mut values = Vec::with_capacity(geometries.len());
    for (done, &index) in geometries.iter().enumerate() {
        run.report(Progress::of(done, geometries.len()));
        run.check().map_err(AlignmentError::Cancelled)?;
        let value = visit_geometry(reader, index, read_error, |moving| {
            let measurement = match &measurement {
                Some(measurement) => measurement,
                None => measurement.insert(Measurement::new(
                    topology,
                    Some(moving),
                    reference,
                    measured,
                )?),
            };
            measurement
                .measure(moving)
                .map_err(|source| item_error(index, source))
        })?;
        values.push(value);
    }
    run.report(Progress::of(values.len(), values.len()));
    Ok(Quantity::new(values, CANONICAL_LENGTH_UNIT))
}

/// Direct RMSD from each of `geometries` of `reader`, in the listed order, to
/// `reference`, measured as stored, as [`StructureSet::rmsd`] measures every
/// geometry of a set in memory.
///
/// `reference` is a geometry index of `reader` or any borrowed view. Each
/// geometry is lent when the reader lends it and read otherwise, one at a
/// time, so memory stays that of one geometry. Ranges and strides are written
/// as `(start..end).step_by(stride)`. The geometries are measured on the
/// calling thread; to spread them over threads, split the list and measure
/// each part.
///
/// ```
/// use kekule::alignment::geometry_rmsd;
/// use kekule::math::Point3;
/// use kekule::structure::{Positions, Geometry, StructureSet};
/// use kekule::topology::AtomSelection;
/// use kekule::units::{Quantity, NANOMETER};
///
/// let topology = kekule::smiles::to_topology("CCO")?;
/// let mut set = StructureSet::new(topology.clone());
/// for x in [0.0, 0.1, 0.2, 0.3] {
///     let points = [Point3::new(x, 0.0, 0.0), Point3::new(0.1, 0.0, 0.0), Point3::new(0.0, 0.1, 0.0)];
///     set.push(Geometry::new(Positions::new(Quantity::new(points, NANOMETER))?), None)?;
/// }
/// let atoms = AtomSelection::all(topology);
/// let every_other = geometry_rmsd(&set, (0..set.len()).step_by(2), 0, &atoms, Default::default())?;
/// let all = set.rmsd(0, &atoms, Default::default())?;
/// assert_eq!(every_other.value(), &[all.value()[0], all.value()[2]]);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Errors
///
/// [`AlignmentError::ReferenceOutOfRange`] for a reference index past the
/// last geometry, [`AlignmentError::Read`] when a geometry cannot be read,
/// [`AlignmentError::Item`] wrapping the failure of the first failing
/// geometry, named by its geometry index, and otherwise the errors of
/// [`rmsd`](super::rmsd) for the pairing and weights.
pub fn geometry_rmsd<'a>(
    reader: &dyn StructureSetRead,
    geometries: impl IntoIterator<Item = usize>,
    reference: impl Into<Reference<'a>>,
    atoms: impl Into<FitAtoms<'a>>,
    options: AlignmentOptions<'_>,
) -> Result<Quantity<Vec<f64>>, AlignmentError> {
    geometry_rmsd_with(
        reader,
        geometries,
        reference,
        atoms,
        options,
        Run::default(),
    )
}

/// [`geometry_rmsd`] under a [`Run`], measuring on the calling thread.
/// Progress counts the geometries measured.
///
/// # Errors
///
/// As [`geometry_rmsd`], and [`AlignmentError::Cancelled`] once `run`'s
/// cancellation is requested.
pub fn geometry_rmsd_with<'a>(
    reader: &dyn StructureSetRead,
    geometries: impl IntoIterator<Item = usize>,
    reference: impl Into<Reference<'a>>,
    atoms: impl Into<FitAtoms<'a>>,
    options: AlignmentOptions<'_>,
    run: Run<'_>,
) -> Result<Quantity<Vec<f64>>, AlignmentError> {
    let measured = Measured::Direct(atoms.into(), options);
    geometry_values(reader, geometries, reference.into(), measured, run)
}

/// RMSD over `measurement` pairs of each of `geometries` of `reader`, in the
/// listed order, after fitting it onto `reference` over `fit` pairs, as
/// [`StructureSet::aligned_rmsd`] measures every geometry of a set in memory;
/// geometries are reached as [`geometry_rmsd`] reaches them.
///
/// # Errors
///
/// As [`geometry_rmsd`], with the errors of
/// [`aligned_rmsd`](super::aligned_rmsd) for the pairings and weights.
pub fn geometry_aligned_rmsd<'a>(
    reader: &dyn StructureSetRead,
    geometries: impl IntoIterator<Item = usize>,
    reference: impl Into<Reference<'a>>,
    fit: impl Into<FitAtoms<'a>>,
    measurement: impl Into<FitAtoms<'a>>,
    options: AlignedRmsdOptions<'_>,
) -> Result<Quantity<Vec<f64>>, AlignmentError> {
    geometry_aligned_rmsd_with(
        reader,
        geometries,
        reference,
        fit,
        measurement,
        options,
        Run::default(),
    )
}

/// [`geometry_aligned_rmsd`] under a [`Run`], measuring on the calling
/// thread. Progress counts the geometries measured.
///
/// # Errors
///
/// As [`geometry_aligned_rmsd`], and [`AlignmentError::Cancelled`] once
/// `run`'s cancellation is requested.
pub fn geometry_aligned_rmsd_with<'a>(
    reader: &dyn StructureSetRead,
    geometries: impl IntoIterator<Item = usize>,
    reference: impl Into<Reference<'a>>,
    fit: impl Into<FitAtoms<'a>>,
    measurement: impl Into<FitAtoms<'a>>,
    options: AlignedRmsdOptions<'_>,
    run: Run<'_>,
) -> Result<Quantity<Vec<f64>>, AlignmentError> {
    let measured = Measured::Aligned(fit.into(), measurement.into(), options);
    geometry_values(reader, geometries, reference.into(), measured, run)
}

fn reference_view<'a>(
    items: &[StructureView<'a>],
    reference: Reference<'a>,
) -> Result<StructureView<'a>, AlignmentError> {
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

/// Validates the pairing once against the first item, or without items
/// against the collection topology.
fn item_pair_sets(
    topology: &Topology,
    first: Option<StructureView<'_>>,
    reference: StructureView<'_>,
    atoms: FitAtoms<'_>,
) -> Result<Vec<Pairs>, AlignmentError> {
    match first {
        Some(item) => atoms.pair_sets(item, reference),
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
    reference: StructureView<'_>,
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

impl StructureSet {
    /// Fits every geometry onto `reference` and applies each transform in place.
    ///
    /// `reference` is a geometry index of this set or any borrowed view
    /// (`&structure`, `&other.get(0).unwrap()`, ...); `atoms` is an
    /// [`AtomSelection`](crate::topology::AtomSelection) or
    /// [`SymmetricSelection`](super::SymmetricSelection) of this layout, or an
    /// [`AtomCorrespondence`](super::AtomCorrespondence) to an independent
    /// reference. Positions move; cells, velocities, and forces rotate; all
    /// other state, the order among it, is kept. The fit
    /// is transactional: on failure nothing changes. Coordinates are fitted as
    /// stored; no imaging or unwrapping is performed. Clone first to keep the
    /// original; a clone shares coordinates until this fit replaces them.
    /// `options` selects the fit weights, in pair order, and the periodic
    /// policy.
    ///
    /// # Errors
    ///
    /// [`AlignmentError::ReferenceOutOfRange`] for a reference index beyond
    /// the set, [`AlignmentError::Item`] wrapping the failure of the first
    /// failing geometry, and otherwise the errors of [`fit`](super::fit) for
    /// the pairing and weights. [`AlignmentError::Transform`] when a fitted
    /// transform cannot be applied to a geometry.
    pub fn superpose<'a>(
        &mut self,
        reference: impl Into<Reference<'a>>,
        atoms: impl Into<FitAtoms<'a>>,
        options: AlignmentOptions<'_>,
    ) -> Result<SuperpositionReport, AlignmentError> {
        let report = {
            let items = self
                .iter()
                .map(|s| s.as_structure_view())
                .collect::<Vec<_>>();
            fit_items(
                self.topology(),
                &items,
                reference.into(),
                atoms.into(),
                options,
            )?
        };
        self.transform_geometries(&transforms(&report))
            .map_err(AlignmentError::Transform)?;
        Ok(report)
    }

    /// Direct RMSD from every geometry to `reference`, measured as stored.
    ///
    /// Never fits or changes the set. Values are in [`CANONICAL_LENGTH_UNIT`],
    /// one per geometry in set order: plot them against the order's values,
    /// or weight an ensemble's samples by
    /// [`Distribution::sample_probabilities`](crate::structure::Distribution::sample_probabilities). `options` selects the pair
    /// weights and the periodic policy.
    ///
    /// # Errors
    ///
    /// [`AlignmentError::ReferenceOutOfRange`] for a reference index beyond
    /// the set, [`AlignmentError::Item`] wrapping the failure of the first
    /// failing geometry, and otherwise the errors of [`rmsd`](super::rmsd) for
    /// the pairing and weights.
    pub fn rmsd<'a>(
        &self,
        reference: impl Into<Reference<'a>>,
        atoms: impl Into<FitAtoms<'a>>,
        options: AlignmentOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
        let items = self
            .iter()
            .map(|s| s.as_structure_view())
            .collect::<Vec<_>>();
        item_rmsd(
            self.topology(),
            &items,
            reference.into(),
            atoms.into(),
            options,
        )
    }

    /// Fits every geometry on `fit` pairs and measures RMSD over `measurement`
    /// pairs, without changing the set or copying coordinates.
    ///
    /// For example, fit a protein backbone and measure ligand motion, or fit
    /// the atoms a scan holds fixed and measure the ones it moves. `options`
    /// holds the fit options and the measurement weights. Values are in
    /// [`CANONICAL_LENGTH_UNIT`], one per geometry in set order.
    ///
    /// # Errors
    ///
    /// [`AlignmentError::ReferenceOutOfRange`] for a reference index beyond
    /// the set, [`AlignmentError::Item`] wrapping the failure of the first
    /// failing geometry, and otherwise the errors of
    /// [`aligned_rmsd`](super::aligned_rmsd) for the pairing and weights.
    pub fn aligned_rmsd<'a>(
        &self,
        reference: impl Into<Reference<'a>>,
        fit: impl Into<FitAtoms<'a>>,
        measurement: impl Into<FitAtoms<'a>>,
        options: AlignedRmsdOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, AlignmentError> {
        let items = self
            .iter()
            .map(|s| s.as_structure_view())
            .collect::<Vec<_>>();
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
