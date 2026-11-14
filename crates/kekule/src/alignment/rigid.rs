//! Rigid superposition of paired atoms: fitting, RMSD, the weighted Kabsch
//! kernel, and their options and errors.

use crate::diagnostics::{
    Cancelled, Countable, IndexOutOfRange, LengthMismatch, LimitExceeded, TopologyMismatch,
};
use std::fmt;

use crate::math::{Matrix3, Point3, RigidTransform, RigidTransformError, Vector3};
use crate::structure::{GeometryError, StructureSetReadError, StructureView};
use crate::topology::{AtomSelection, InstanceAtomId, TopologyAtomIndex, TopologyPerceptionError};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};

use super::{
    AlignedRmsdOptions, AmbiguousLabel, AtomCorrespondence, CorrespondenceSide,
    LabelCorrespondence, SymmetricSelection,
};

/// Moving/reference topology atom index pairs.
pub(crate) type Pairs = Vec<(TopologyAtomIndex, TopologyAtomIndex)>;

const MIN_FIT_ATOMS: usize = 3;

/// Variance-relative rank threshold used for centered selected geometry.
///
/// Rank two requires `lambda_2 > lambda_1 * RANK_RELATIVE_TOLERANCE`, where
/// the lambdas are the two largest eigenvalues of the weighted scatter matrix.
const RANK_RELATIVE_TOLERANCE: f64 = 1.0e-12;
const JACOBI_RELATIVE_TOLERANCE: f64 = 64.0 * f64::EPSILON;
const JACOBI_MAX_SWEEPS: usize = 24;

/// The atoms paired between a moving and a reference view.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum FitAtoms<'a> {
    /// Each selected atom pairs with itself; both views and the selection share
    /// one topology layout. Pairs follow sorted selection order.
    Selection(&'a AtomSelection),
    /// Explicit ordered pairs between two topologies, in pair order.
    Correspondence(&'a AtomCorrespondence),
    /// A selection paired under every symmetry permutation, keeping the one
    /// with the smallest RMSD; both views and the selection share one
    /// topology layout. Pairs follow sorted selection order on the reference
    /// side.
    Symmetric(&'a SymmetricSelection),
}

impl<'a> From<&'a AtomSelection> for FitAtoms<'a> {
    fn from(selection: &'a AtomSelection) -> Self {
        Self::Selection(selection)
    }
}

impl<'a> From<&'a AtomCorrespondence> for FitAtoms<'a> {
    fn from(correspondence: &'a AtomCorrespondence) -> Self {
        Self::Correspondence(correspondence)
    }
}

impl<'a> From<&'a LabelCorrespondence> for FitAtoms<'a> {
    fn from(correspondence: &'a LabelCorrespondence) -> Self {
        Self::Correspondence(correspondence.correspondence())
    }
}

impl<'a> From<&'a SymmetricSelection> for FitAtoms<'a> {
    fn from(selection: &'a SymmetricSelection) -> Self {
        Self::Symmetric(selection)
    }
}

impl FitAtoms<'_> {
    /// Number of atom pairs.
    pub fn pair_count(self) -> usize {
        match self {
            Self::Selection(selection) => selection.indices().len(),
            Self::Correspondence(correspondence) => correspondence.len(),
            Self::Symmetric(selection) => selection.len(),
        }
    }

    /// Whether there are no atom pairs.
    pub fn is_empty(self) -> bool {
        self.pair_count() == 0
    }

    /// Validates the pairing against both views and returns its candidate
    /// pair sets: one, or one per symmetry permutation with the identity
    /// first. Every set has the same length and reference order.
    pub(crate) fn pair_sets(
        self,
        moving: StructureView<'_>,
        reference: StructureView<'_>,
    ) -> Result<Vec<Pairs>, AlignmentError> {
        match self {
            Self::Selection(selection) => {
                check_selection(selection, moving, reference)?;
                Ok(vec![selection
                    .indices()
                    .iter()
                    .map(|index| (*index, *index))
                    .collect()])
            }
            Self::Correspondence(correspondence) => {
                correspondence.ensure_compatible(moving.topology(), reference.topology())?;
                Ok(vec![correspondence.index_pairs().to_vec()])
            }
            Self::Symmetric(symmetric) => {
                check_selection(symmetric.selection(), moving, reference)?;
                Ok(symmetric.pair_sets())
            }
        }
    }
}

fn check_selection(
    selection: &AtomSelection,
    moving: StructureView<'_>,
    reference: StructureView<'_>,
) -> Result<(), AlignmentError> {
    if !moving.topology().shares_layout(reference.topology()) {
        return Err(AlignmentError::TopologyMismatch(
            TopologyMismatch::Reference,
        ));
    }
    selection
        .ensure_compatible(moving.topology_arc())
        .map_err(|_| AlignmentError::TopologyMismatch(TopologyMismatch::Selection))
}

/// Fits `moving` onto `reference` with the Kabsch algorithm.
///
/// The fit minimizes `sum(w_i * |R x_i + t - y_i|^2)` subject to a proper
/// right-handed rotation, over at least three pairs spanning two dimensions.
/// `options` selects the weights, which follow pair order, and the periodic
/// policy; the default weights are uniform. The result maps moving coordinates
/// into the reference coordinate system, and its RMSD is the weighted post-fit
/// value in [`CANONICAL_LENGTH_UNIT`].
///
/// # Errors
///
/// - [`AlignmentError::TopologyMismatch`] when a selection's views, the
///   selection itself or a correspondence do not share one topology layout.
/// - [`AlignmentError::WeightCountMismatch`],
///   [`AlignmentError::NonFiniteWeight`], or
///   [`AlignmentError::NonPositiveWeight`] for invalid explicit weights.
/// - [`AlignmentError::PeriodicCoordinates`] when either view has a cell under
///   [`PeriodicPolicy::RejectPeriodic`].
/// - [`AlignmentError::InsufficientSelectedAtoms`] for fewer than three fit
///   pairs, and [`AlignmentError::DegenerateGeometry`] when the paired
///   coordinates do not span two dimensions.
/// - [`AlignmentError::NumericalFailure`] or
///   [`AlignmentError::InvalidRigidTransform`] when the solution is not finite.
pub fn fit<'a>(
    moving: StructureView<'_>,
    reference: StructureView<'_>,
    atoms: impl Into<FitAtoms<'a>>,
    options: AlignmentOptions<'_>,
) -> Result<RigidAlignment, AlignmentError> {
    let atoms = atoms.into();
    let sets = atoms.pair_sets(moving, reference)?;
    let weights = NormalizedWeights::new(options.weighting, atoms.pair_count())?;
    best_kabsch(moving, reference, &sets, weights, options.periodic_policy)
}

/// Direct RMSD over paired atoms, measured exactly as stored without fitting.
///
/// `options` selects the weights, which follow pair order, and the periodic
/// policy. The value is in [`CANONICAL_LENGTH_UNIT`].
///
/// # Errors
///
/// - [`AlignmentError::TopologyMismatch`] when a selection's views, the
///   selection itself or a correspondence do not share one topology layout.
/// - [`AlignmentError::WeightCountMismatch`],
///   [`AlignmentError::NonFiniteWeight`], or
///   [`AlignmentError::NonPositiveWeight`] for invalid explicit weights.
/// - [`AlignmentError::PeriodicCoordinates`] when either view has a cell under
///   [`PeriodicPolicy::RejectPeriodic`].
/// - [`AlignmentError::EmptySelection`] when there are no pairs.
/// - [`AlignmentError::NumericalFailure`] when the value is not finite.
pub fn rmsd<'a>(
    moving: StructureView<'_>,
    reference: StructureView<'_>,
    atoms: impl Into<FitAtoms<'a>>,
    options: AlignmentOptions<'_>,
) -> Result<Quantity<f64>, AlignmentError> {
    let atoms = atoms.into();
    let sets = atoms.pair_sets(moving, reference)?;
    let weights = NormalizedWeights::new(options.weighting, atoms.pair_count())?;
    let value = best_measure(
        moving,
        reference,
        &sets,
        weights,
        options.periodic_policy,
        None,
    )?;
    Ok(Quantity::new(value, CANONICAL_LENGTH_UNIT))
}

/// Fits `moving` onto `reference` on `fit` pairs and measures the RMSD over
/// `measurement` pairs under that fit, without changing either view.
///
/// For example, fit a protein backbone and measure ligand motion. The best fit
/// fixes the transform; a symmetric measurement then takes its own best
/// permutation under it. `options` holds the fit options, whose periodic
/// policy applies to both steps, and the measurement weights. The value is in
/// [`CANONICAL_LENGTH_UNIT`].
///
/// # Errors
///
/// Any error of [`fit`] for the fit pairs, and of [`rmsd`] for the
/// measurement pairs.
pub fn aligned_rmsd<'a>(
    moving: StructureView<'_>,
    reference: StructureView<'_>,
    fit: impl Into<FitAtoms<'a>>,
    measurement: impl Into<FitAtoms<'a>>,
    options: AlignedRmsdOptions<'_>,
) -> Result<Quantity<f64>, AlignmentError> {
    let fit = fit.into();
    let measurement = measurement.into();
    let fit_sets = fit.pair_sets(moving, reference)?;
    let measurement_sets = measurement.pair_sets(moving, reference)?;
    check_fit_count(fit.pair_count())?;
    check_measurement_count(measurement.pair_count())?;
    let fit_weights = NormalizedWeights::new(options.fit().weighting, fit.pair_count())?;
    let measurement_weights =
        NormalizedWeights::new(options.measurement_weighting(), measurement.pair_count())?;
    let policy = options.fit().periodic_policy;
    let alignment = best_kabsch(moving, reference, &fit_sets, fit_weights, policy)?;
    let value = best_measure(
        moving,
        reference,
        &measurement_sets,
        measurement_weights,
        policy,
        Some(alignment.transform()),
    )?;
    Ok(Quantity::new(value, CANONICAL_LENGTH_UNIT))
}

/// The positions of one structure as the fitting and measuring kernels
/// read them: a structure's own, or the fitted atoms of one loaded alone,
/// which pairs then index.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Coordinates<'a> {
    /// Positions in the canonical length unit.
    pub(crate) points: &'a [Point3],
    /// Whether the structure has a periodic cell.
    pub(crate) periodic: bool,
}

impl<'a> From<StructureView<'a>> for Coordinates<'a> {
    fn from(view: StructureView<'a>) -> Self {
        Self {
            points: *view.payload().positions().values().value(),
            periodic: view.cell().is_some(),
        }
    }
}

/// The fit with the smallest RMSD over candidate pair sets; ties keep the
/// earliest set. Any failing set fails the fit.
pub(crate) fn best_kabsch<'m, 'r>(
    moving: impl Into<Coordinates<'m>>,
    reference: impl Into<Coordinates<'r>>,
    sets: &[Pairs],
    weights: NormalizedWeights<'_>,
    policy: PeriodicPolicy,
) -> Result<RigidAlignment, AlignmentError> {
    best_coordinates_kabsch(moving.into(), reference.into(), sets, weights, policy)
}

/// [`best_kabsch`] on coordinates, kept out of line: inlined into its
/// callers, such as the pairwise matrix's loop, the kernel compiles to
/// slower code and a fit takes about 40 % longer.
#[inline(never)]
fn best_coordinates_kabsch(
    moving: Coordinates<'_>,
    reference: Coordinates<'_>,
    sets: &[Pairs],
    weights: NormalizedWeights<'_>,
    policy: PeriodicPolicy,
) -> Result<RigidAlignment, AlignmentError> {
    let mut best: Option<RigidAlignment> = None;
    for pairs in sets {
        let alignment = kabsch_pairs(moving, reference, pairs, weights, policy)?;
        if best
            .as_ref()
            .is_none_or(|best| alignment.rmsd.value() < best.rmsd.value())
        {
            best = Some(alignment);
        }
    }
    best.ok_or(AlignmentError::EmptySelection)
}

/// The smallest measured RMSD over candidate pair sets.
pub(crate) fn best_measure<'m, 'r>(
    moving: impl Into<Coordinates<'m>>,
    reference: impl Into<Coordinates<'r>>,
    sets: &[Pairs],
    weights: NormalizedWeights<'_>,
    policy: PeriodicPolicy,
    transform: Option<RigidTransform>,
) -> Result<f64, AlignmentError> {
    let (moving, reference) = (moving.into(), reference.into());
    let mut best = f64::INFINITY;
    for pairs in sets {
        best = best.min(measure_pairs(
            moving, reference, pairs, weights, policy, transform,
        )?);
    }
    if sets.is_empty() {
        return Err(AlignmentError::EmptySelection);
    }
    Ok(best)
}

fn check_periodic(
    moving: Coordinates<'_>,
    reference: Coordinates<'_>,
    policy: PeriodicPolicy,
) -> Result<(), AlignmentError> {
    let moving_periodic = moving.periodic;
    let reference_periodic = reference.periodic;
    if policy == PeriodicPolicy::RejectPeriodic && (moving_periodic || reference_periodic) {
        return Err(AlignmentError::PeriodicCoordinates {
            moving: moving_periodic,
            reference: reference_periodic,
        });
    }
    Ok(())
}

/// RMSD requires at least one pair.
pub(crate) fn check_measurement_count(pairs: usize) -> Result<(), AlignmentError> {
    if pairs == 0 {
        return Err(AlignmentError::EmptySelection);
    }
    Ok(())
}

/// A determined rigid fit requires at least [`MIN_FIT_ATOMS`] pairs.
pub(crate) fn check_fit_count(pairs: usize) -> Result<(), AlignmentError> {
    if pairs < MIN_FIT_ATOMS {
        return Err(AlignmentError::InsufficientSelectedAtoms {
            selected: pairs,
            minimum: MIN_FIT_ATOMS,
        });
    }
    Ok(())
}

/// Weighted RMSD over validated pairs, optionally after moving coordinates
/// are transformed.
fn measure_pairs(
    moving: Coordinates<'_>,
    reference: Coordinates<'_>,
    pairs: &[(TopologyAtomIndex, TopologyAtomIndex)],
    weights: NormalizedWeights<'_>,
    policy: PeriodicPolicy,
    transform: Option<RigidTransform>,
) -> Result<f64, AlignmentError> {
    check_measurement_count(pairs.len())?;
    check_periodic(moving, reference, policy)?;
    let moving_positions = moving.points;
    let reference_positions = reference.points;
    let mut squared_residual = CompensatedSum::default();
    let mut weight_sum = CompensatedSum::default();
    for (pair_index, (moving_index, reference_index)) in pairs.iter().enumerate() {
        let mut point = moving_positions[moving_index.index()];
        if let Some(transform) = transform {
            point = transform.transform_point(point);
        }
        let residual = point - reference_positions[reference_index.index()];
        let weight = weights.at(pair_index);
        squared_residual.add(weight * residual.norm_squared());
        weight_sum.add(weight);
    }
    let mean_squared_residual = squared_residual.value() / weight_sum.value();
    if !mean_squared_residual.is_finite() || mean_squared_residual < 0.0 {
        return Err(AlignmentError::NumericalFailure);
    }
    Ok(mean_squared_residual.sqrt())
}

fn kabsch_pairs(
    moving: Coordinates<'_>,
    reference: Coordinates<'_>,
    pairs: &[(TopologyAtomIndex, TopologyAtomIndex)],
    weights: NormalizedWeights<'_>,
    policy: PeriodicPolicy,
) -> Result<RigidAlignment, AlignmentError> {
    let selected_atom_count = pairs.len();
    check_fit_count(selected_atom_count)?;
    check_periodic(moving, reference, policy)?;

    let moving_positions = moving.points;
    let reference_positions = reference.points;

    let mut moments = AlignmentMoments::default();
    for (pair_index, (moving_index, reference_index)) in pairs.iter().enumerate() {
        moments.add(
            moving_positions[moving_index.index()] - Point3::origin(),
            reference_positions[reference_index.index()] - Point3::origin(),
            weights.at(pair_index),
        );
    }
    let moments = moments.finish()?;

    ensure_rank_two(moments.moving_scatter, AlignmentGeometry::Moving)?;
    ensure_rank_two(moments.reference_scatter, AlignmentGeometry::Reference)?;

    let rotation = proper_rotation(moments.cross_covariance)?;
    let translation =
        moments.reference_centroid - rotation.transform_vector(moments.moving_centroid);
    let transform = RigidTransform::from_canonical(rotation, translation)
        .map_err(AlignmentError::InvalidRigidTransform)?;

    let mut weighted_squared_residual = CompensatedSum::default();
    for (pair_index, (moving_index, reference_index)) in pairs.iter().enumerate() {
        let moving_centered =
            (moving_positions[moving_index.index()] - Point3::origin()) - moments.moving_centroid;
        let reference_centered = (reference_positions[reference_index.index()] - Point3::origin())
            - moments.reference_centroid;
        let residual = rotation.transform_vector(moving_centered) - reference_centered;
        weighted_squared_residual.add(weights.at(pair_index) * residual.norm_squared());
    }
    let mean_squared_residual = weighted_squared_residual.value() / moments.weight_sum;
    if !mean_squared_residual.is_finite() || mean_squared_residual < 0.0 {
        return Err(AlignmentError::NumericalFailure);
    }

    Ok(RigidAlignment {
        transform,
        rmsd: Quantity::new(mean_squared_residual.sqrt(), CANONICAL_LENGTH_UNIT),
        selected_atom_count,
    })
}

/// Options shared by fitting and RMSD measurement.
#[derive(Debug, Clone, Copy, Default)]
pub struct AlignmentOptions<'a> {
    pub(super) weighting: Weighting<'a>,
    pub(super) periodic_policy: PeriodicPolicy,
}

impl<'a> AlignmentOptions<'a> {
    /// Per-pair weighting.
    pub const fn weighting(&self) -> Weighting<'a> {
        self.weighting
    }

    /// Sets [`Self::weighting`].
    #[must_use]
    pub const fn with_weighting(mut self, weighting: Weighting<'a>) -> Self {
        self.weighting = weighting;
        self
    }

    /// Handling of views carrying periodic cells.
    pub const fn periodic_policy(&self) -> PeriodicPolicy {
        self.periodic_policy
    }

    /// Sets [`Self::periodic_policy`].
    #[must_use]
    pub const fn with_periodic_policy(mut self, periodic_policy: PeriodicPolicy) -> Self {
        self.periodic_policy = periodic_policy;
        self
    }
}

/// Per-pair weighting for fitting and RMSD.
#[derive(Debug, Clone, Copy, Default)]
#[non_exhaustive]
pub enum Weighting<'a> {
    /// Give every pair equal weight.
    #[default]
    Uniform,
    /// Positive finite weights in pair order.
    Explicit(&'a [f64]),
}

/// Policy for views carrying periodic cells.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum PeriodicPolicy {
    /// Reject either input when it carries a periodic cell.
    RejectPeriodic,
    /// Ignore cells and use the stored Cartesian coordinates directly.
    ///
    /// This performs no imaging, wrapping, unwrapping, minimum-image
    /// correction, or molecule reconstruction.
    #[default]
    UseStoredCoordinates,
}

/// One successful proper rigid fit from moving to reference coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct RigidAlignment {
    transform: RigidTransform,
    rmsd: Quantity<f64>,
    selected_atom_count: usize,
}

impl RigidAlignment {
    /// Returns the transform mapping moving coordinates into reference
    /// coordinates.
    pub const fn transform(&self) -> RigidTransform {
        self.transform
    }

    /// Returns post-fit weighted RMSD in [`CANONICAL_LENGTH_UNIT`].
    pub const fn rmsd(&self) -> Quantity<f64> {
        self.rmsd
    }

    /// Returns the number of selected atom correspondences used by the fit.
    pub const fn selected_atom_count(&self) -> usize {
        self.selected_atom_count
    }
}

/// The geometry whose numerical rank was insufficient for a determined fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AlignmentGeometry {
    /// The selected moving coordinates are rank deficient.
    Moving,
    /// The selected reference coordinates are rank deficient.
    Reference,
    /// The moving/reference cross-covariance cannot determine two axes.
    CrossCovariance,
}

/// A structured rigid-alignment failure.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum AlignmentError {
    /// The reference does not share the moving structure's topology layout
    /// ([`TopologyMismatch::Reference`]), the atom selection does not share
    /// it ([`TopologyMismatch::Selection`]), or a structure does not share
    /// the layout of its side of an explicit correspondence
    /// ([`TopologyMismatch::Moving`] or [`TopologyMismatch::Reference`]).
    TopologyMismatch(TopologyMismatch),
    /// Dense-order pairing ([`AtomCorrespondence::from_same_layout`]) was
    /// requested for topologies with different layouts.
    LayoutMismatch,
    /// A pair names an atom absent from that side's topology.
    InvalidAtom {
        /// The side whose atom is absent.
        side: CorrespondenceSide,
        /// The position of the pair in the input.
        pair_index: usize,
        /// The absent atom.
        atom: InstanceAtomId,
    },
    /// A pair repeats an atom already paired on the same side.
    DuplicateAtom {
        /// The side with the repeated atom.
        side: CorrespondenceSide,
        /// The position of the repeating pair in the input.
        pair_index: usize,
        /// The repeated atom.
        atom: InstanceAtomId,
    },
    /// Hierarchy labels name more than one selected atom on a side
    /// ([`AtomCorrespondence::from_labels`]).
    AmbiguousLabels {
        /// The side whose selection repeats the labels.
        side: CorrespondenceSide,
        /// Every repeated label of that side, in the dense order of its first
        /// atom.
        labels: Vec<AmbiguousLabel>,
    },
    /// The selection admits more permutations ([`Resource::Permutations`]) or
    /// the symmetry search examined more candidate mappings
    /// ([`Resource::SearchStates`]) than [`SymmetryLimits`](super::SymmetryLimits) allows.
    ///
    /// [`Resource::Permutations`]: crate::diagnostics::Resource::Permutations
    /// [`Resource::SearchStates`]: crate::diagnostics::Resource::SearchStates
    SymmetryLimitExceeded(LimitExceeded),
    /// Perceiving the inferred hydrogen counts a symmetry search needs failed.
    Perception(TopologyPerceptionError),
    /// An atom's hydrogen count is unknown even after perception, so the
    /// symmetry search cannot compare it.
    UnknownHydrogenCount {
        /// The atom with the unknown count.
        atom: TopologyAtomIndex,
    },
    /// Fewer than three atoms were selected.
    InsufficientSelectedAtoms {
        /// Actual selected atom count.
        selected: usize,
        /// Required minimum selected atom count.
        minimum: usize,
    },
    /// Selected geometry does not span two dimensions at the documented
    /// scale-relative tolerance.
    DegenerateGeometry {
        /// The rank-deficient geometry.
        geometry: AlignmentGeometry,
    },
    /// The explicit weight count does not match selection length.
    WeightCountMismatch(LengthMismatch),
    /// One explicit weight is NaN or infinite.
    NonFiniteWeight {
        /// Zero-based position in selection order.
        selection_index: usize,
    },
    /// One explicit weight is zero or negative.
    NonPositiveWeight {
        /// Zero-based position in selection order.
        selection_index: usize,
    },
    /// Periodic coordinates were rejected by policy.
    PeriodicCoordinates {
        /// Whether the moving structure has a cell.
        moving: bool,
        /// Whether the reference structure has a cell.
        reference: bool,
    },
    /// Fixed-size numerical accumulation or decomposition did not produce a
    /// finite, converged solution.
    NumericalFailure,
    /// The existing rigid-transform invariant rejected the fitted result.
    InvalidRigidTransform(RigidTransformError),
    /// RMSD requires at least one atom pair.
    EmptySelection,
    /// A collection reference index does not exist.
    ReferenceOutOfRange(IndexOutOfRange),
    /// One collection item failed.
    Item {
        /// The index of the failed item.
        index: usize,
        /// The item's failure.
        source: Box<AlignmentError>,
    },
    /// A fitted transform could not be applied to one item's state.
    Transform(GeometryError),
    /// A geometry cannot be read; see the source.
    Read(StructureSetReadError),
    /// The measurement was cancelled.
    Cancelled(Cancelled),
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch(mismatch) => mismatch.fmt(formatter),
            Self::LayoutMismatch => formatter.write_str(
                "The moving and reference topologies have different layouts. Give explicit atom pairs.",
            ),
            Self::InvalidAtom {
                side,
                pair_index,
                atom,
            } => write!(
                formatter,
                "Atom pair {pair_index} names {atom}. The {side} topology has no such atom.",
            ),
            Self::DuplicateAtom {
                side,
                pair_index,
                atom,
            } => write!(formatter, "Atom pair {pair_index} repeats the {side} {atom}."),
            Self::AmbiguousLabels { side, labels } => {
                match labels.len() {
                    1 => write!(formatter, "A {side} label names more than one selected atom:")?,
                    count => write!(
                        formatter,
                        "{count} {side} labels name more than one selected atom:"
                    )?,
                }
                for (position, label) in labels.iter().take(3).enumerate() {
                    let separator = if position == 0 { " " } else { "; " };
                    write!(formatter, "{separator}{label}")?;
                }
                if labels.len() > 3 {
                    write!(formatter, "; and {} more", labels.len() - 3)?;
                }
                formatter.write_str(".")
            }
            Self::SymmetryLimitExceeded(limit) => write!(
                formatter,
                "The symmetry search exceeded its {}.",
                limit.describe()
            ),
            Self::Perception(_) => formatter.write_str("Cannot perceive hydrogen counts."),
            Self::UnknownHydrogenCount { atom } => {
                write!(formatter, "The hydrogen count of {atom} is unknown.")
            }
            Self::InsufficientSelectedAtoms { selected, minimum } => write!(
                formatter,
                "The selection has {}. Rigid alignment needs at least {minimum}.",
                crate::diagnostics::counted(*selected, "atom", "atoms")
            ),
            Self::DegenerateGeometry { geometry } => write!(
                formatter,
                "The {geometry} geometry does not span two dimensions."
            ),
            Self::WeightCountMismatch(error) => error.fmt(formatter),
            Self::NonFiniteWeight { selection_index } => write!(
                formatter,
                "The alignment weight at selection index {selection_index} is not finite."
            ),
            Self::NonPositiveWeight { selection_index } => write!(
                formatter,
                "The alignment weight at selection index {selection_index} is not positive."
            ),
            Self::PeriodicCoordinates { moving, reference } => {
                let structures = match (moving, reference) {
                    (true, true) => "Both structures have",
                    (true, false) => "The moving structure has",
                    (false, true) => "The reference structure has",
                    (false, false) => {
                        return formatter
                            .write_str("Alignment does not support periodic coordinates.")
                    }
                };
                write!(
                    formatter,
                    "{structures} a periodic cell. Alignment does not support periodic coordinates."
                )
            }
            Self::NumericalFailure => {
                formatter.write_str("The rigid alignment did not reach a finite solution.")
            }
            Self::InvalidRigidTransform(_) => {
                formatter.write_str("The fitted rigid transform is invalid.")
            }
            Self::EmptySelection => {
                formatter.write_str("The selection is empty. RMSD needs at least one atom pair.")
            }
            Self::ReferenceOutOfRange(error) => error.fmt(formatter),
            Self::Item { index, .. } => write!(formatter, "Cannot align structure {index}."),
            Self::Transform(_) => formatter.write_str("Cannot apply the fitted transform."),
            Self::Read(error) => fmt::Display::fmt(error, formatter),
            Self::Cancelled(error) => fmt::Display::fmt(error, formatter),
        }
    }
}

impl fmt::Display for AlignmentGeometry {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Moving => formatter.write_str("selected moving"),
            Self::Reference => formatter.write_str("selected reference"),
            Self::CrossCovariance => formatter.write_str("selected cross-covariance"),
        }
    }
}

impl std::error::Error for AlignmentError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Perception(error) => Some(error),
            Self::InvalidRigidTransform(error) => Some(error),
            Self::Item { source, .. } => Some(source.as_ref()),
            Self::Transform(error) => Some(error),
            Self::Read(error) => std::error::Error::source(error),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum NormalizedWeights<'a> {
    Uniform,
    Explicit { values: &'a [f64], maximum: f64 },
}

impl<'a> NormalizedWeights<'a> {
    pub(crate) fn new(
        weighting: Weighting<'a>,
        selected_atom_count: usize,
    ) -> Result<Self, AlignmentError> {
        match weighting {
            Weighting::Uniform => Ok(Self::Uniform),
            Weighting::Explicit(values) => {
                if values.len() != selected_atom_count {
                    return Err(AlignmentError::WeightCountMismatch(LengthMismatch::new(
                        Countable::Weight,
                        selected_atom_count,
                        values.len(),
                    )));
                }
                let mut maximum = 0.0_f64;
                for (selection_index, weight) in values.iter().copied().enumerate() {
                    if !weight.is_finite() {
                        return Err(AlignmentError::NonFiniteWeight { selection_index });
                    }
                    if weight <= 0.0 {
                        return Err(AlignmentError::NonPositiveWeight { selection_index });
                    }
                    maximum = maximum.max(weight);
                }
                Ok(Self::Explicit { values, maximum })
            }
        }
    }

    pub(crate) fn at(self, selection_index: usize) -> f64 {
        match self {
            Self::Uniform => 1.0,
            Self::Explicit { values, maximum } => values[selection_index] / maximum,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        if self.sum.abs() >= value.abs() {
            self.correction += (self.sum - next) + value;
        } else {
            self.correction += (value - next) + self.sum;
        }
        self.sum = next;
    }

    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct CompensatedVector {
    values: [CompensatedSum; 3],
}

impl CompensatedVector {
    fn add(&mut self, value: Vector3) {
        for (sum, value) in self.values.iter_mut().zip(components(value)) {
            sum.add(value);
        }
    }

    fn value(self) -> Vector3 {
        let [x, y, z] = self.values.map(CompensatedSum::value);
        Vector3::new(x, y, z)
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct CompensatedMatrix {
    values: [[CompensatedSum; 3]; 3],
}

impl CompensatedMatrix {
    fn add_outer(&mut self, left: Vector3, right: Vector3, scale: f64) {
        for (row, left_value) in components(left).into_iter().enumerate() {
            for (column, right_value) in components(right).into_iter().enumerate() {
                self.values[row][column].add(scale * left_value * right_value);
            }
        }
    }

    fn value(self) -> Matrix {
        self.values.map(|row| row.map(CompensatedSum::value))
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct AlignmentMoments {
    weight_sum: CompensatedSum,
    moving_centroid: CompensatedVector,
    reference_centroid: CompensatedVector,
    moving_scatter: CompensatedMatrix,
    reference_scatter: CompensatedMatrix,
    cross_covariance: CompensatedMatrix,
}

impl AlignmentMoments {
    fn add(&mut self, moving: Vector3, reference: Vector3, weight: f64) {
        let old_weight_sum = self.weight_sum.value();
        self.weight_sum.add(weight);
        let weight_sum = self.weight_sum.value();
        let moving_centroid = self.moving_centroid.value();
        let reference_centroid = self.reference_centroid.value();
        let moving_delta = moving - moving_centroid;
        let reference_delta = reference - reference_centroid;
        let ratio = weight / weight_sum;
        self.moving_centroid.add(moving_delta * ratio);
        self.reference_centroid.add(reference_delta * ratio);

        let scatter_scale = weight * old_weight_sum / weight_sum;
        self.moving_scatter
            .add_outer(moving_delta, moving_delta, scatter_scale);
        self.reference_scatter
            .add_outer(reference_delta, reference_delta, scatter_scale);
        self.cross_covariance
            .add_outer(moving_delta, reference_delta, scatter_scale);
    }

    fn finish(self) -> Result<FinishedMoments, AlignmentError> {
        let result = FinishedMoments {
            weight_sum: self.weight_sum.value(),
            moving_centroid: self.moving_centroid.value(),
            reference_centroid: self.reference_centroid.value(),
            moving_scatter: self.moving_scatter.value(),
            reference_scatter: self.reference_scatter.value(),
            cross_covariance: self.cross_covariance.value(),
        };
        if !result.weight_sum.is_finite()
            || result.weight_sum <= 0.0
            || !result.moving_centroid.is_finite()
            || !result.reference_centroid.is_finite()
            || !matrix_is_finite(result.moving_scatter)
            || !matrix_is_finite(result.reference_scatter)
            || !matrix_is_finite(result.cross_covariance)
        {
            return Err(AlignmentError::NumericalFailure);
        }
        Ok(result)
    }
}

#[derive(Debug, Clone, Copy)]
struct FinishedMoments {
    weight_sum: f64,
    moving_centroid: Vector3,
    reference_centroid: Vector3,
    moving_scatter: Matrix,
    reference_scatter: Matrix,
    cross_covariance: Matrix,
}

type Matrix = [[f64; 3]; 3];
type DominantSingularBases = ([Vector3; 2], [Vector3; 2], [f64; 3]);

fn ensure_rank_two(scatter: Matrix, geometry: AlignmentGeometry) -> Result<(), AlignmentError> {
    let eigenvalues = symmetric_eigenvalues(scatter)?;
    let largest = eigenvalues[0];
    let second = eigenvalues[1];
    if largest <= 0.0 || second <= largest * RANK_RELATIVE_TOLERANCE {
        return Err(AlignmentError::DegenerateGeometry { geometry });
    }
    Ok(())
}

fn symmetric_eigenvalues(mut matrix: Matrix) -> Result<[f64; 3], AlignmentError> {
    if !matrix_is_finite(matrix) {
        return Err(AlignmentError::NumericalFailure);
    }
    for _ in 0..JACOBI_MAX_SWEEPS {
        let (p, q, off_diagonal) = largest_off_diagonal(matrix);
        let scale = matrix_scale(matrix);
        if off_diagonal <= JACOBI_RELATIVE_TOLERANCE * scale {
            let mut eigenvalues = [matrix[0][0], matrix[1][1], matrix[2][2]];
            eigenvalues.sort_by(|left, right| right.total_cmp(left));
            if eigenvalues.into_iter().all(f64::is_finite) {
                return Ok(eigenvalues);
            }
            return Err(AlignmentError::NumericalFailure);
        }
        apply_symmetric_jacobi_rotation(&mut matrix, p, q);
        if !matrix_is_finite(matrix) {
            return Err(AlignmentError::NumericalFailure);
        }
    }
    Err(AlignmentError::NumericalFailure)
}

fn apply_symmetric_jacobi_rotation(matrix: &mut Matrix, p: usize, q: usize) {
    let app = matrix[p][p];
    let aqq = matrix[q][q];
    let apq = matrix[p][q];
    let scale = app.abs().max(aqq.abs()).max(apq.abs());
    let zeta = ((aqq / scale) - (app / scale)) / (2.0 * (apq / scale));
    let tangent = if zeta >= 0.0 {
        1.0 / (zeta + (1.0 + zeta * zeta).sqrt())
    } else {
        -1.0 / (-zeta + (1.0 + zeta * zeta).sqrt())
    };
    let cosine = 1.0 / (1.0 + tangent * tangent).sqrt();
    let sine = tangent * cosine;

    for k in [0, 1, 2] {
        if k == p || k == q {
            continue;
        }
        let akp = matrix[k][p];
        let akq = matrix[k][q];
        let new_kp = cosine * akp - sine * akq;
        let new_kq = sine * akp + cosine * akq;
        matrix[k][p] = new_kp;
        matrix[p][k] = new_kp;
        matrix[k][q] = new_kq;
        matrix[q][k] = new_kq;
    }
    matrix[p][p] = cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq;
    matrix[q][q] = sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq;
    matrix[p][q] = 0.0;
    matrix[q][p] = 0.0;
}

fn proper_rotation(cross_covariance: Matrix) -> Result<Matrix3, AlignmentError> {
    let (left, right, singular_values) = dominant_singular_bases(cross_covariance)?;
    if singular_values[0] <= 0.0
        || singular_values[1] <= singular_values[0] * RANK_RELATIVE_TOLERANCE
    {
        return Err(AlignmentError::DegenerateGeometry {
            geometry: AlignmentGeometry::CrossCovariance,
        });
    }

    let (left_first, left_second) = orthonormal_pair(left[0], left[1])?;
    let (right_first, right_second) = orthonormal_pair(right[0], right[1])?;
    let left_basis = Matrix3::from_columns(left_first, left_second, left_first.cross(left_second));
    let right_basis =
        Matrix3::from_columns(right_first, right_second, right_first.cross(right_second));
    Ok(right_basis.multiply(left_basis.transpose()))
}

/// Returns the two dominant left/right singular vectors and all singular
/// values. A fixed-size one-sided Jacobi SVD avoids squaring the covariance
/// condition number through normal equations.
fn dominant_singular_bases(
    cross_covariance: Matrix,
) -> Result<DominantSingularBases, AlignmentError> {
    if !matrix_is_finite(cross_covariance) {
        return Err(AlignmentError::NumericalFailure);
    }
    let mut work = cross_covariance;
    let mut right_vectors = identity_matrix();
    for _ in 0..JACOBI_MAX_SWEEPS {
        let mut changed = false;
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            let first = matrix_column(work, p);
            let second = matrix_column(work, q);
            let alpha = fused_dot(first, first);
            let beta = fused_dot(second, second);
            let coupling = fused_dot(first, second);
            if !alpha.is_finite() || !beta.is_finite() || !coupling.is_finite() {
                return Err(AlignmentError::NumericalFailure);
            }
            let comparison_scale = alpha.sqrt() * beta.sqrt();
            if coupling.abs() <= JACOBI_RELATIVE_TOLERANCE * comparison_scale {
                continue;
            }
            let scale = alpha.max(beta).max(coupling.abs());
            let zeta = ((beta / scale) - (alpha / scale)) / (2.0 * (coupling / scale));
            let tangent = if zeta >= 0.0 {
                1.0 / (zeta + (1.0 + zeta * zeta).sqrt())
            } else {
                -1.0 / (-zeta + (1.0 + zeta * zeta).sqrt())
            };
            let cosine = 1.0 / (1.0 + tangent * tangent).sqrt();
            let sine = tangent * cosine;
            rotate_matrix_columns(&mut work, p, q, cosine, sine);
            rotate_matrix_columns(&mut right_vectors, p, q, cosine, sine);
            changed = true;
        }
        if !changed {
            break;
        }
    }
    if !one_sided_svd_converged(work) {
        return Err(AlignmentError::NumericalFailure);
    }

    let mut order = [0_usize, 1, 2];
    order.sort_by(|left, right| {
        norm_squared(matrix_column(work, *right))
            .total_cmp(&norm_squared(matrix_column(work, *left)))
    });
    let singular_values = order.map(|column| norm_squared(matrix_column(work, column)).sqrt());
    if !singular_values.into_iter().all(f64::is_finite) {
        return Err(AlignmentError::NumericalFailure);
    }
    if singular_values[0] == 0.0 || singular_values[1] == 0.0 {
        return Err(AlignmentError::DegenerateGeometry {
            geometry: AlignmentGeometry::CrossCovariance,
        });
    }
    let left = [
        matrix_column(work, order[0]) * (1.0 / singular_values[0]),
        matrix_column(work, order[1]) * (1.0 / singular_values[1]),
    ];
    let right = [
        matrix_column(right_vectors, order[0]),
        matrix_column(right_vectors, order[1]),
    ];
    Ok((left, right, singular_values))
}

fn one_sided_svd_converged(matrix: Matrix) -> bool {
    let mut order = [0_usize, 1, 2];
    order.sort_by(|left, right| {
        norm_squared(matrix_column(matrix, *right))
            .total_cmp(&norm_squared(matrix_column(matrix, *left)))
    });
    let first = matrix_column(matrix, order[0]);
    let second = matrix_column(matrix, order[1]);
    let scale = norm_squared(first).sqrt() * norm_squared(second).sqrt();
    fused_dot(first, second).abs() <= 8.0 * JACOBI_RELATIVE_TOLERANCE * scale
}

fn rotate_matrix_columns(matrix: &mut Matrix, p: usize, q: usize, cosine: f64, sine: f64) {
    for row in matrix {
        let first = row[p];
        let second = row[q];
        row[p] = cosine * first - sine * second;
        row[q] = sine * first + cosine * second;
    }
}

fn orthonormal_pair(first: Vector3, second: Vector3) -> Result<(Vector3, Vector3), AlignmentError> {
    let first = normalize(first).ok_or(AlignmentError::NumericalFailure)?;
    let second = second - first * fused_dot(first, second);
    let second = normalize(second).ok_or(AlignmentError::NumericalFailure)?;
    Ok((first, second))
}

fn normalize(value: Vector3) -> Option<Vector3> {
    let norm = norm_squared(value).sqrt();
    if !norm.is_finite() || norm == 0.0 {
        return None;
    }
    Some(value * (1.0 / norm))
}

fn largest_off_diagonal(matrix: Matrix) -> (usize, usize, f64) {
    [(0, 1), (0, 2), (1, 2)]
        .into_iter()
        .map(|(p, q)| (p, q, matrix[p][q].abs()))
        .max_by(|left, right| left.2.total_cmp(&right.2))
        .expect("fixed non-empty pair list")
}

fn matrix_scale(matrix: Matrix) -> f64 {
    matrix
        .into_iter()
        .flatten()
        .map(f64::abs)
        .fold(0.0, f64::max)
}

fn matrix_column(matrix: Matrix, column: usize) -> Vector3 {
    Vector3::new(matrix[0][column], matrix[1][column], matrix[2][column])
}

const fn identity_matrix() -> Matrix {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

const fn components(value: Vector3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

/// The dot product with fused multiply-adds, which the singular value
/// decomposition keeps; [`Vector3::dot`] rounds each product.
fn fused_dot(left: Vector3, right: Vector3) -> f64 {
    left.x
        .mul_add(right.x, left.y.mul_add(right.y, left.z * right.z))
}

fn norm_squared(value: Vector3) -> f64 {
    fused_dot(value, value)
}

fn matrix_is_finite(matrix: Matrix) -> bool {
    matrix.into_iter().flatten().all(f64::is_finite)
}

#[cfg(test)]
mod tests;
