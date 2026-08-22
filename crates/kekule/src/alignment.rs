//! Same-topology, selection-based rigid molecular alignment.
//!
//! The returned transform always maps moving coordinates into reference
//! coordinates:
//!
//! ```text
//! aligned = result.transform().transform_point(moving)
//! ```
//!
//! Alignment is read-only. It does not image periodic coordinates or
//! materialize transformed canonical coordinates.
//!
//! # Examples
//!
//! ```
//! use kekule::alignment::kabsch;
//! use kekule::core::{Atom, BondOrder, Element, MoleculeEditor};
//! use kekule::geometry::Point3;
//! use kekule::structure::{Model, Positions};
//! use kekule::topology::{AtomSelection, TopologyBuilder};
//! use kekule::units::{Quantity, ANGSTROM};
//! use std::sync::Arc;
//!
//! let mut graph = MoleculeEditor::new();
//! let mut previous = None;
//! for _ in 0..3 {
//!     let atom = graph.add_atom(Atom::new(Element::from_symbol("C").unwrap()))?;
//!     if let Some(parent) = previous {
//!         graph.add_bond(parent, atom, BondOrder::Single)?;
//!     }
//!     previous = Some(atom);
//! }
//! let molecule = graph.finish()?;
//! let mut builder = TopologyBuilder::new();
//! let definition = builder.add_molecule_definition(&molecule)?;
//! builder.add_instance(definition)?;
//! let topology = Arc::new(builder.build()?);
//! let moving_points = [
//!     Point3::new(0.0, 0.0, 0.0),
//!     Point3::new(1.0, 0.0, 0.0),
//!     Point3::new(0.0, 1.0, 0.0),
//! ];
//! let reference_points = moving_points.map(|point| Point3::new(
//!     point.x + 2.0,
//!     point.y - 1.0,
//!     point.z + 0.5,
//! ));
//! let moving = Model::new(
//!     Arc::clone(&topology),
//!     Positions::new(
//!         &topology,
//!         Quantity::new(moving_points, ANGSTROM),
//!     )?,
//! )?;
//! let reference = Model::new(
//!     Arc::clone(&topology),
//!     Positions::new(
//!         &topology,
//!         Quantity::new(reference_points, ANGSTROM),
//!     )?,
//! )?;
//! let selection = AtomSelection::from_atoms(
//!     &topology,
//!     topology.atom_ids().iter().copied(),
//! )?;
//!
//! let result = kabsch(moving.view(), reference.view(), &selection)?;
//! let aligned = result.transform().transform_point(moving_points[1]);
//! assert!((aligned.x - reference_points[1].x).abs() < 1.0e-12);
//! assert!(result.rmsd().to_value() < 1.0e-12);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::fmt;
use std::sync::Arc;

use crate::geometry::{Matrix3, Point3, RigidTransform, RigidTransformError, Vector3};
use crate::structure::ModelView;
use crate::topology::AtomSelection;
use crate::units::{Quantity, MODEL_LENGTH_UNIT};

const MIN_SELECTED_ATOMS: usize = 3;

/// Variance-relative rank threshold used for centered selected geometry.
///
/// Rank two requires `lambda_2 > lambda_1 * RANK_RELATIVE_TOLERANCE`, where
/// the lambdas are the two largest eigenvalues of the weighted scatter matrix.
const RANK_RELATIVE_TOLERANCE: f64 = 1.0e-12;
const JACOBI_RELATIVE_TOLERANCE: f64 = 64.0 * f64::EPSILON;
const JACOBI_MAX_SWEEPS: usize = 24;

/// Fits `moving` onto `reference` with uniform weights.
///
/// This is equivalent to [`kabsch_with_options`] with
/// [`KabschOptions::default`]. The result maps moving coordinates into the
/// reference coordinate system.
pub fn kabsch(
    moving: ModelView<'_>,
    reference: ModelView<'_>,
    selection: &AtomSelection,
) -> Result<RigidAlignment, AlignmentError> {
    kabsch_with_options(moving, reference, selection, KabschOptions::default())
}

/// Fits `moving` onto `reference` with explicit weighting and periodic policy.
///
/// Correspondence follows the selection's sorted dense-index order. The two
/// views and selection must share one `Arc<Topology>` allocation. Explicit weights
/// use selection order, not complete-topology order.
///
/// The fit minimizes `sum(w_i * |R x_i + t - y_i|^2)` subject to a proper
/// right-handed rotation. The returned RMSD is the weighted post-fit value in
/// [`MODEL_LENGTH_UNIT`].
pub fn kabsch_with_options(
    moving: ModelView<'_>,
    reference: ModelView<'_>,
    selection: &AtomSelection,
    options: KabschOptions<'_>,
) -> Result<RigidAlignment, AlignmentError> {
    if !Arc::ptr_eq(moving.topology_arc(), reference.topology_arc()) {
        return Err(AlignmentError::TopologyMismatch);
    }
    selection
        .ensure_compatible(moving.topology_arc())
        .map_err(|_| AlignmentError::SelectionTopologyMismatch)?;

    let selected_atom_count = selection.indices().len();
    if selected_atom_count < MIN_SELECTED_ATOMS {
        return Err(AlignmentError::InsufficientSelectedAtoms {
            selected: selected_atom_count,
            minimum: MIN_SELECTED_ATOMS,
        });
    }

    let weights = NormalizedWeights::new(options.weighting, selected_atom_count)?;
    let moving_periodic = moving.cell().is_some();
    let reference_periodic = reference.cell().is_some();
    if options.periodic_policy == PeriodicAlignmentPolicy::RejectPeriodic
        && (moving_periodic || reference_periodic)
    {
        return Err(AlignmentError::PeriodicCoordinates {
            moving: moving_periodic,
            reference: reference_periodic,
        });
    }

    let moving_positions = moving.positions().values();
    let reference_positions = reference.positions().values();
    let moving_positions = moving_positions.value();
    let reference_positions = reference_positions.value();

    let mut moments = AlignmentMoments::default();
    for (selection_index, dense_index) in selection.indices().iter().copied().enumerate() {
        moments.add(
            point_components(moving_positions[dense_index.index()]),
            point_components(reference_positions[dense_index.index()]),
            weights.at(selection_index),
        );
    }
    let moments = moments.finish()?;

    ensure_rank_two(moments.moving_scatter, AlignmentGeometry::Moving)?;
    ensure_rank_two(moments.reference_scatter, AlignmentGeometry::Reference)?;

    let rotation = proper_rotation(moments.cross_covariance)?;
    let moving_centroid = components_vector(moments.moving_centroid);
    let reference_centroid = components_vector(moments.reference_centroid);
    let translation = reference_centroid - rotation.transform_vector(moving_centroid);
    let transform = RigidTransform::new(rotation, translation)
        .map_err(AlignmentError::InvalidRigidTransform)?;

    let mut weighted_squared_residual = CompensatedSum::default();
    for (selection_index, dense_index) in selection.indices().iter().copied().enumerate() {
        let moving_centered = point_components(moving_positions[dense_index.index()]);
        let reference_centered = point_components(reference_positions[dense_index.index()]);
        let moving_centered = subtract(moving_centered, moments.moving_centroid);
        let reference_centered = subtract(reference_centered, moments.reference_centroid);
        let residual = rotation.transform_vector(components_vector(moving_centered))
            - components_vector(reference_centered);
        weighted_squared_residual.add(weights.at(selection_index) * residual.norm_squared());
    }
    let mean_squared_residual = weighted_squared_residual.value() / moments.weight_sum;
    if !mean_squared_residual.is_finite() || mean_squared_residual < 0.0 {
        return Err(AlignmentError::NumericalFailure);
    }

    Ok(RigidAlignment {
        transform,
        rmsd: Quantity::new(mean_squared_residual.sqrt(), MODEL_LENGTH_UNIT),
        selected_atom_count,
    })
}

/// Options for [`kabsch_with_options`].
#[derive(Debug, Clone, Copy, Default)]
pub struct KabschOptions<'a> {
    /// Per-selected-atom weighting policy.
    pub weighting: AlignmentWeighting<'a>,
    /// Handling of models carrying periodic cells.
    pub periodic_policy: PeriodicAlignmentPolicy,
}

/// Per-selected-atom weighting for rigid alignment.
#[derive(Debug, Clone, Copy, Default)]
pub enum AlignmentWeighting<'a> {
    /// Give every selected atom equal weight.
    #[default]
    Uniform,
    /// Use positive finite weights in selection order.
    Explicit(&'a [f64]),
}

/// Policy for models carrying periodic cells.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum PeriodicAlignmentPolicy {
    /// Reject either input when it carries a periodic cell.
    #[default]
    RejectPeriodic,
    /// Ignore cells and fit the stored Cartesian coordinates directly.
    ///
    /// This performs no imaging, wrapping, unwrapping, minimum-image
    /// correction, or molecule reconstruction.
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

    /// Returns post-fit weighted RMSD in [`MODEL_LENGTH_UNIT`].
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
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum AlignmentError {
    /// Moving and reference views do not share one topology allocation.
    TopologyMismatch,
    /// The atom selection belongs to another topology allocation.
    SelectionTopologyMismatch,
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
    WeightCountMismatch {
        /// Selected atom count.
        expected: usize,
        /// Supplied weight count.
        actual: usize,
    },
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
        /// Whether the moving model has a cell.
        moving: bool,
        /// Whether the reference model has a cell.
        reference: bool,
    },
    /// Fixed-size numerical accumulation or decomposition did not produce a
    /// finite, converged solution.
    NumericalFailure,
    /// The existing rigid-transform invariant rejected the fitted result.
    InvalidRigidTransform(RigidTransformError),
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => formatter.write_str(
                "moving and reference models belong to different topology allocations",
            ),
            Self::SelectionTopologyMismatch => {
                formatter.write_str("atom selection belongs to a different topology allocation")
            }
            Self::InsufficientSelectedAtoms { selected, minimum } => write!(
                formatter,
                "rigid alignment requires at least {minimum} selected atoms, but received {selected}"
            ),
            Self::DegenerateGeometry { geometry } => {
                write!(formatter, "{geometry} geometry is rank deficient")
            }
            Self::WeightCountMismatch { expected, actual } => write!(
                formatter,
                "rigid alignment requires {expected} weights, but received {actual}"
            ),
            Self::NonFiniteWeight { selection_index } => write!(
                formatter,
                "alignment weight at selection index {selection_index} is not finite"
            ),
            Self::NonPositiveWeight { selection_index } => write!(
                formatter,
                "alignment weight at selection index {selection_index} is not strictly positive"
            ),
            Self::PeriodicCoordinates {
                moving,
                reference,
            } => write!(
                formatter,
                "periodic alignment is rejected by policy (moving cell: {moving}, reference cell: {reference})"
            ),
            Self::NumericalFailure => {
                formatter.write_str("rigid alignment numerical solution failed")
            }
            Self::InvalidRigidTransform(error) => {
                write!(formatter, "fitted rigid transform is invalid: {error}")
            }
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
            Self::InvalidRigidTransform(error) => Some(error),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum NormalizedWeights<'a> {
    Uniform,
    Explicit { values: &'a [f64], maximum: f64 },
}

impl<'a> NormalizedWeights<'a> {
    fn new(
        weighting: AlignmentWeighting<'a>,
        selected_atom_count: usize,
    ) -> Result<Self, AlignmentError> {
        match weighting {
            AlignmentWeighting::Uniform => Ok(Self::Uniform),
            AlignmentWeighting::Explicit(values) => {
                if values.len() != selected_atom_count {
                    return Err(AlignmentError::WeightCountMismatch {
                        expected: selected_atom_count,
                        actual: values.len(),
                    });
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

    fn at(self, selection_index: usize) -> f64 {
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
    fn add(&mut self, value: [f64; 3]) {
        for (sum, value) in self.values.iter_mut().zip(value) {
            sum.add(value);
        }
    }

    fn value(self) -> [f64; 3] {
        self.values.map(CompensatedSum::value)
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct CompensatedMatrix {
    values: [[CompensatedSum; 3]; 3],
}

impl CompensatedMatrix {
    fn add_outer(&mut self, left: [f64; 3], right: [f64; 3], scale: f64) {
        for (row, left_value) in left.into_iter().enumerate() {
            for (column, right_value) in right.into_iter().enumerate() {
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
    fn add(&mut self, moving: [f64; 3], reference: [f64; 3], weight: f64) {
        let old_weight_sum = self.weight_sum.value();
        self.weight_sum.add(weight);
        let weight_sum = self.weight_sum.value();
        let moving_centroid = self.moving_centroid.value();
        let reference_centroid = self.reference_centroid.value();
        let moving_delta = subtract(moving, moving_centroid);
        let reference_delta = subtract(reference, reference_centroid);
        let ratio = weight / weight_sum;
        self.moving_centroid.add(scale(moving_delta, ratio));
        self.reference_centroid.add(scale(reference_delta, ratio));

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
            || !vector_is_finite(result.moving_centroid)
            || !vector_is_finite(result.reference_centroid)
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
    moving_centroid: [f64; 3],
    reference_centroid: [f64; 3],
    moving_scatter: Matrix,
    reference_scatter: Matrix,
    cross_covariance: Matrix,
}

type Matrix = [[f64; 3]; 3];
type DominantSingularBases = ([[f64; 3]; 2], [[f64; 3]; 2], [f64; 3]);

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
    let left_third = cross(left_first, left_second);
    let right_third = cross(right_first, right_second);
    let left_basis = Matrix3::from_columns(
        components_vector(left_first),
        components_vector(left_second),
        components_vector(left_third),
    );
    let right_basis = Matrix3::from_columns(
        components_vector(right_first),
        components_vector(right_second),
        components_vector(right_third),
    );
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
            let alpha = dot(first, first);
            let beta = dot(second, second);
            let coupling = dot(first, second);
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
        scale(matrix_column(work, order[0]), 1.0 / singular_values[0]),
        scale(matrix_column(work, order[1]), 1.0 / singular_values[1]),
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
    dot(first, second).abs() <= 8.0 * JACOBI_RELATIVE_TOLERANCE * scale
}

fn rotate_matrix_columns(matrix: &mut Matrix, p: usize, q: usize, cosine: f64, sine: f64) {
    for row in matrix {
        let first = row[p];
        let second = row[q];
        row[p] = cosine * first - sine * second;
        row[q] = sine * first + cosine * second;
    }
}

fn orthonormal_pair(
    first: [f64; 3],
    second: [f64; 3],
) -> Result<([f64; 3], [f64; 3]), AlignmentError> {
    let first = normalize(first).ok_or(AlignmentError::NumericalFailure)?;
    let second = subtract(second, scale(first, dot(first, second)));
    let second = normalize(second).ok_or(AlignmentError::NumericalFailure)?;
    Ok((first, second))
}

fn normalize(value: [f64; 3]) -> Option<[f64; 3]> {
    let norm = norm_squared(value).sqrt();
    if !norm.is_finite() || norm == 0.0 {
        return None;
    }
    Some(scale(value, 1.0 / norm))
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

fn matrix_column(matrix: Matrix, column: usize) -> [f64; 3] {
    [matrix[0][column], matrix[1][column], matrix[2][column]]
}

const fn identity_matrix() -> Matrix {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

fn point_components(point: Point3) -> [f64; 3] {
    [point.x, point.y, point.z]
}

fn components_vector(components: [f64; 3]) -> Vector3 {
    Vector3::new(components[0], components[1], components[2])
}

fn subtract(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn scale(value: [f64; 3], factor: f64) -> [f64; 3] {
    [value[0] * factor, value[1] * factor, value[2] * factor]
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0].mul_add(right[0], left[1].mul_add(right[1], left[2] * right[2]))
}

fn norm_squared(value: [f64; 3]) -> f64 {
    dot(value, value)
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn vector_is_finite(value: [f64; 3]) -> bool {
    value.into_iter().all(f64::is_finite)
}

fn matrix_is_finite(matrix: Matrix) -> bool {
    matrix.into_iter().flatten().all(f64::is_finite)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Atom, BondOrder, Element};
    use crate::geometry::PeriodicCell;
    use crate::structure::{Model, Positions};
    use crate::topology::{AtomSelection, Topology, TopologyBuilder};
    use crate::units::{Quantity, ANGSTROM, NANOMETER};

    fn topology(atom_count: usize) -> Arc<Topology> {
        let mut graph = crate::core::MoleculeEditor::new();
        let mut previous = None;
        for _ in 0..atom_count {
            let atom = graph
                .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
                .unwrap();
            if let Some(parent) = previous {
                graph.add_bond(parent, atom, BondOrder::Single).unwrap();
            }
            previous = Some(atom);
        }
        let molecule = graph.finish().unwrap();
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(&molecule).unwrap();
        builder.add_instance(definition).unwrap();
        Arc::new(builder.build().unwrap())
    }

    fn model(topology: &Arc<Topology>, points: &[Point3]) -> Model {
        let positions = Positions::new(topology, Quantity::new(points, ANGSTROM)).unwrap();
        Model::new(Arc::clone(topology), positions).unwrap()
    }

    fn models(moving: &[Point3], reference: &[Point3]) -> (Arc<Topology>, Model, Model) {
        assert_eq!(moving.len(), reference.len());
        let topology = topology(moving.len());
        let moving = model(&topology, moving);
        let reference = model(&topology, reference);
        (topology, moving, reference)
    }

    fn all(topology: &Arc<Topology>) -> AtomSelection {
        AtomSelection::from_atoms(topology, topology.atom_ids().iter().copied()).unwrap()
    }

    fn fixture_points() -> [Point3; 5] {
        [
            Point3::new(-1.2, 0.4, 2.1),
            Point3::new(0.7, -1.5, 0.3),
            Point3::new(2.4, 0.8, -0.6),
            Point3::new(-0.3, 2.2, 1.1),
            Point3::new(1.3, 1.7, 2.8),
        ]
    }

    fn non_axis_transform() -> RigidTransform {
        let axis = Vector3::new(1.0, -2.0, 3.0);
        let axis = axis / axis.norm();
        let angle = 0.731_f64;
        let cosine = angle.cos();
        let sine = angle.sin();
        let one_minus_cosine = 1.0 - cosine;
        let rotation = Matrix3::from_columns(
            Vector3::new(
                cosine + axis.x * axis.x * one_minus_cosine,
                axis.y * axis.x * one_minus_cosine + axis.z * sine,
                axis.z * axis.x * one_minus_cosine - axis.y * sine,
            ),
            Vector3::new(
                axis.x * axis.y * one_minus_cosine - axis.z * sine,
                cosine + axis.y * axis.y * one_minus_cosine,
                axis.z * axis.y * one_minus_cosine + axis.x * sine,
            ),
            Vector3::new(
                axis.x * axis.z * one_minus_cosine + axis.y * sine,
                axis.y * axis.z * one_minus_cosine - axis.x * sine,
                cosine + axis.z * axis.z * one_minus_cosine,
            ),
        );
        RigidTransform::new(rotation, Vector3::new(2.3, -4.1, 0.8)).unwrap()
    }

    fn transform_points(points: &[Point3], transform: RigidTransform) -> Vec<Point3> {
        points
            .iter()
            .copied()
            .map(|point| transform.transform_point(point))
            .collect()
    }

    fn assert_close(left: f64, right: f64, tolerance: f64) {
        assert!(
            (left - right).abs() <= tolerance,
            "{left} is not within {tolerance} of {right}"
        );
    }

    fn assert_point_close(left: Point3, right: Point3, tolerance: f64) {
        assert_close(left.x, right.x, tolerance);
        assert_close(left.y, right.y, tolerance);
        assert_close(left.z, right.z, tolerance);
    }

    fn assert_transform_close(left: RigidTransform, right: RigidTransform, tolerance: f64) {
        for (left, right) in left
            .rotation()
            .columns()
            .into_iter()
            .zip(right.rotation().columns())
        {
            assert_close(left.x, right.x, tolerance);
            assert_close(left.y, right.y, tolerance);
            assert_close(left.z, right.z, tolerance);
        }
        let left = left.translation();
        let right = right.translation();
        assert_close(left.x, right.x, tolerance);
        assert_close(left.y, right.y, tolerance);
        assert_close(left.z, right.z, tolerance);
    }

    #[test]
    fn identity_fit_reports_model_units_and_preserves_inputs() {
        let points = fixture_points();
        let (topology, moving, reference) = models(&points, &points);
        let selection = all(&topology);
        let moving_before = moving.positions().values().value().to_vec();
        let reference_before = reference.positions().values().value().to_vec();
        let selection_before = selection.clone();

        let result = kabsch(moving.view(), reference.view(), &selection).unwrap();

        assert_transform_close(result.transform(), RigidTransform::identity(), 1.0e-12);
        assert_close(result.rmsd().to_value(), 0.0, 1.0e-14);
        assert_eq!(result.rmsd().unit(), MODEL_LENGTH_UNIT);
        assert_eq!(result.selected_atom_count(), points.len());
        assert_eq!(
            moving.positions().values().to_value(),
            moving_before.as_slice()
        );
        assert_eq!(
            reference.positions().values().to_value(),
            reference_before.as_slice()
        );
        assert_eq!(selection, selection_before);
    }

    #[test]
    fn translation_fit_maps_moving_to_reference_not_the_inverse() {
        let moving = fixture_points();
        let expected =
            RigidTransform::new(Matrix3::identity(), Vector3::new(3.2, -1.7, 4.4)).unwrap();
        let reference = transform_points(&moving, expected);
        let (topology, moving_model, reference_model) = models(&moving, &reference);

        let result = kabsch(moving_model.view(), reference_model.view(), &all(&topology)).unwrap();

        assert_transform_close(result.transform(), expected, 1.0e-12);
        for (moving, reference) in moving.into_iter().zip(reference) {
            assert_point_close(
                result.transform().transform_point(moving),
                reference,
                1.0e-11,
            );
            assert!((result.transform().transform_point(reference).x - moving.x).abs() > 1.0);
        }
    }

    #[test]
    fn asymmetric_non_axis_rotation_and_translation_are_recovered() {
        let moving = fixture_points();
        let expected = non_axis_transform();
        let reference = transform_points(&moving, expected);
        let (topology, moving, reference_model) = models(&moving, &reference);

        let result = kabsch(moving.view(), reference_model.view(), &all(&topology)).unwrap();

        assert_transform_close(result.transform(), expected, 2.0e-12);
        assert!(result.rmsd().to_value() < 2.0e-12);
    }

    #[test]
    fn subset_fit_ignores_unselected_atoms() {
        let moving = fixture_points();
        let expected = non_axis_transform();
        let mut reference = transform_points(&moving, expected);
        reference[3] = Point3::new(100.0, -80.0, 40.0);
        reference[4] = Point3::new(-70.0, 60.0, -30.0);
        let (topology, moving_model, reference_model) = models(&moving, &reference);
        let selection =
            AtomSelection::from_atoms(&topology, topology.atom_ids()[..3].iter().copied()).unwrap();

        let result = kabsch(moving_model.view(), reference_model.view(), &selection).unwrap();

        assert_eq!(result.selected_atom_count(), 3);
        for index in 0..3 {
            assert_point_close(
                result.transform().transform_point(moving[index]),
                reference[index],
                3.0e-11,
            );
        }
    }

    #[test]
    fn noisy_fit_has_expected_nonzero_improved_rmsd() {
        let moving = fixture_points();
        let expected = non_axis_transform();
        let mut reference = transform_points(&moving, expected);
        reference[0].x += 0.08;
        reference[1].y -= 0.05;
        reference[4].z += 0.03;
        let (topology, moving_model, reference_model) = models(&moving, &reference);
        let selection = all(&topology);
        let result = kabsch(moving_model.view(), reference_model.view(), &selection).unwrap();
        let unaligned = moving
            .iter()
            .zip(&reference)
            .map(|(moving, reference)| (*moving - *reference).norm_squared())
            .sum::<f64>()
            / moving.len() as f64;
        let unaligned = unaligned.sqrt();

        assert!(result.rmsd().to_value() > 0.01);
        assert!(result.rmsd().to_value() < 0.08);
        assert!(result.rmsd().to_value() < unaligned);
    }

    #[test]
    fn explicit_weights_change_fit_and_common_scale_does_not() {
        let moving = fixture_points();
        let expected = non_axis_transform();
        let mut reference = transform_points(&moving, expected);
        reference[4].x += 2.5;
        reference[4].y -= 1.5;
        let (topology, moving_model, reference_model) = models(&moving, &reference);
        let selection = all(&topology);
        let uniform = kabsch(moving_model.view(), reference_model.view(), &selection).unwrap();
        let weights = [10.0, 10.0, 10.0, 10.0, 0.1];
        let scaled_weights = [30.0, 30.0, 30.0, 30.0, 0.3];
        let weighted = kabsch_with_options(
            moving_model.view(),
            reference_model.view(),
            &selection,
            KabschOptions {
                weighting: AlignmentWeighting::Explicit(&weights),
                ..KabschOptions::default()
            },
        )
        .unwrap();
        let scaled = kabsch_with_options(
            moving_model.view(),
            reference_model.view(),
            &selection,
            KabschOptions {
                weighting: AlignmentWeighting::Explicit(&scaled_weights),
                ..KabschOptions::default()
            },
        )
        .unwrap();

        let uniform_error = (uniform.transform().transform_point(moving[0]) - reference[0]).norm();
        let weighted_error =
            (weighted.transform().transform_point(moving[0]) - reference[0]).norm();
        assert!(weighted_error < uniform_error);
        assert_transform_close(weighted.transform(), scaled.transform(), 2.0e-12);
        assert_close(
            weighted.rmsd().to_value(),
            scaled.rmsd().to_value(),
            2.0e-12,
        );
    }

    #[test]
    fn explicit_weights_follow_selection_indices_order_not_semantic_id_insertion_order() {
        let moving = fixture_points();
        let expected = non_axis_transform();
        let mut reference = transform_points(&moving, expected);
        reference[0].x += 1.7;
        reference[4].x -= 1.9;
        let (topology, moving_model, reference_model) = models(&moving, &reference);
        let atom_ids = topology.atom_ids();
        let selection = AtomSelection::from_atoms(
            &topology,
            [atom_ids[4], atom_ids[1], atom_ids[3], atom_ids[0]],
        )
        .unwrap();

        assert_eq!(
            selection
                .indices()
                .iter()
                .map(|index| index.index())
                .collect::<Vec<_>>(),
            [0, 1, 3, 4]
        );

        // These weights make dense index 0 dominant in selection.indices()
        // order. Reversing the endpoint weights instead makes the first
        // inserted semantic ID (dense index 4) dominant.
        let selection_order_weights = [100.0, 1.0, 1.0, 0.1];
        let insertion_order_interpretation = [0.1, 1.0, 1.0, 100.0];
        let selection_order_fit = kabsch_with_options(
            moving_model.view(),
            reference_model.view(),
            &selection,
            KabschOptions {
                weighting: AlignmentWeighting::Explicit(&selection_order_weights),
                ..KabschOptions::default()
            },
        )
        .unwrap();
        let insertion_order_fit = kabsch_with_options(
            moving_model.view(),
            reference_model.view(),
            &selection,
            KabschOptions {
                weighting: AlignmentWeighting::Explicit(&insertion_order_interpretation),
                ..KabschOptions::default()
            },
        )
        .unwrap();

        let selected_dense_zero_error =
            (selection_order_fit.transform().transform_point(moving[0]) - reference[0]).norm();
        let insertion_dense_zero_error =
            (insertion_order_fit.transform().transform_point(moving[0]) - reference[0]).norm();
        let selected_dense_four_error =
            (selection_order_fit.transform().transform_point(moving[4]) - reference[4]).norm();
        let insertion_dense_four_error =
            (insertion_order_fit.transform().transform_point(moving[4]) - reference[4]).norm();
        assert!(selected_dense_zero_error < insertion_dense_zero_error * 0.1);
        assert!(insertion_dense_four_error < selected_dense_four_error * 0.1);
    }

    #[test]
    fn weighted_noisy_fit_matches_high_precision_svd_golden() {
        let moving = [
            Point3::new(-1.4, 0.2, 2.1),
            Point3::new(0.3, -1.7, 0.5),
            Point3::new(2.2, 0.9, -0.8),
            Point3::new(-0.6, 2.4, 1.3),
            Point3::new(1.5, 1.1, 2.7),
            Point3::new(-2.0, -0.9, -1.2),
        ];
        let reference = [
            Point3::new(1.873307692308, -4.135076923077, -0.201230769231),
            Point3::new(3.672153846154, -2.542538461538, 1.554384615385),
            Point3::new(1.272692307692, -0.149923076923, 2.197230769231),
            Point3::new(-0.144923076923, -2.663769230769, -0.328307692308),
            Point3::new(0.199692307692, -3.570923076923, 2.290230769231),
            Point3::new(4.000538461538, -1.492384615385, -1.161153846154),
        ];
        let weights = [0.7, 3.25, 1.1, 5.5, 0.4, 2.75];
        let (topology, moving_model, reference_model) = models(&moving, &reference);

        // Independent reference: mpmath 1.3.0 at 100 decimal digits, using
        // weighted centroids/cross-covariance and its general SVD, followed by
        // V * diag(1, 1, det(V * U^T)) * U^T. The 3e-12 transform and 5e-14
        // RMSD tolerances cover only f64 rounding relative to that reference.
        let expected = RigidTransform::new(
            Matrix3::from_columns(
                Vector3::new(
                    -0.23369675064066617,
                    0.30951178758237726,
                    0.9217311332962317,
                ),
                Vector3::new(
                    -0.9253543656919893,
                    0.22023658458981056,
                    -0.30856951356702656,
                ),
                Vector3::new(-0.2985048184448125, -0.925039620857278, 0.2349395096816475),
            ),
            Vector3::new(2.333984523603927, -1.795397432011596, 0.6560333319022371),
        )
        .unwrap();
        let expected_rmsd = 0.029397741136029633_f64;

        let result = kabsch_with_options(
            moving_model.view(),
            reference_model.view(),
            &all(&topology),
            KabschOptions {
                weighting: AlignmentWeighting::Explicit(&weights),
                ..KabschOptions::default()
            },
        )
        .unwrap();

        assert_transform_close(result.transform(), expected, 3.0e-12);
        assert_close(result.rmsd().to_value(), expected_rmsd, 5.0e-14);
        assert_close(result.transform().rotation().determinant(), 1.0, 1.0e-12);
    }

    #[test]
    fn mirrored_coordinates_keep_a_proper_nonzero_residual() {
        let moving = fixture_points();
        let reference = moving.map(|point| Point3::new(-point.x, point.y, point.z));
        let (topology, moving, reference) = models(&moving, &reference);

        let result = kabsch(moving.view(), reference.view(), &all(&topology)).unwrap();

        assert_close(result.transform().rotation().determinant(), 1.0, 1.0e-12);
        assert!(result.rmsd().to_value() > 0.1);
        RigidTransform::new(
            result.transform().rotation(),
            result.transform().translation(),
        )
        .unwrap();
    }

    #[test]
    fn planar_non_collinear_selection_succeeds() {
        let moving = [
            Point3::new(-1.0, -0.5, 0.0),
            Point3::new(2.0, -0.5, 0.0),
            Point3::new(0.2, 1.7, 0.0),
            Point3::new(-0.4, 0.8, 0.0),
        ];
        let expected = non_axis_transform();
        let reference = transform_points(&moving, expected);
        let (topology, moving, reference) = models(&moving, &reference);

        let result = kabsch(moving.view(), reference.view(), &all(&topology)).unwrap();

        assert_transform_close(result.transform(), expected, 3.0e-12);
        assert!(result.rmsd().to_value() < 2.0e-12);
    }

    #[test]
    fn rank_two_point_sets_with_rank_one_cross_covariance_are_rejected() {
        let moving = [
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
        ];
        let reference = [
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
        ];
        let (topology, moving, reference) = models(&moving, &reference);

        assert_eq!(
            kabsch(moving.view(), reference.view(), &all(&topology)),
            Err(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::CrossCovariance,
            })
        );
    }

    #[test]
    fn fewer_than_three_selected_atoms_are_structured_failures() {
        let points = fixture_points();
        let (topology, moving, reference) = models(&points, &points);
        for count in 0..3 {
            let selection =
                AtomSelection::from_atoms(&topology, topology.atom_ids()[..count].iter().copied())
                    .unwrap();
            assert_eq!(
                kabsch(moving.view(), reference.view(), &selection),
                Err(AlignmentError::InsufficientSelectedAtoms {
                    selected: count,
                    minimum: 3,
                })
            );
        }
    }

    #[test]
    fn coincident_and_collinear_geometry_are_rejected() {
        let valid = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        let coincident = [Point3::new(2.0, -1.0, 0.5); 4];
        let collinear = [
            Point3::new(-2.0, 0.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ];
        let (topology, moving, reference) = models(&coincident, &valid);
        assert_eq!(
            kabsch(moving.view(), reference.view(), &all(&topology)),
            Err(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Moving,
            })
        );
        let (topology, moving, reference) = models(&valid, &collinear);
        assert_eq!(
            kabsch(moving.view(), reference.view(), &all(&topology)),
            Err(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Reference,
            })
        );
    }

    #[test]
    fn near_collinear_rank_test_is_scale_relative_and_thresholded() {
        let points = |height: f64, scale: f64| {
            [
                Point3::new(-scale, 0.0, 0.0),
                Point3::new(scale, 0.0, 0.0),
                Point3::new(0.0, height * scale, 0.0),
                Point3::new(0.5 * scale, -height * scale, 0.0),
            ]
        };
        for scale in [1.0, 1.0e8] {
            let accepted = points(1.0e-5, scale);
            let (topology, moving, reference) = models(&accepted, &accepted);
            assert!(kabsch(moving.view(), reference.view(), &all(&topology)).is_ok());

            let rejected = points(1.0e-7, scale);
            let (topology, moving, reference) = models(&rejected, &rejected);
            assert!(matches!(
                kabsch(moving.view(), reference.view(), &all(&topology)),
                Err(AlignmentError::DegenerateGeometry {
                    geometry: AlignmentGeometry::Moving
                })
            ));
        }
    }

    #[test]
    fn exact_shared_topology_is_required_for_views_and_selection() {
        let points = fixture_points();
        let topology_a = topology(points.len());
        let topology_b = topology(points.len());
        assert!(topology_a.same_layout(&topology_b));
        assert!(!Arc::ptr_eq(&topology_a, &topology_b));
        let moving_a = model(&topology_a, &points);
        let reference_a = model(&topology_a, &points);
        let reference_b = model(&topology_b, &points);
        assert_eq!(
            kabsch(moving_a.view(), reference_b.view(), &all(&topology_a)),
            Err(AlignmentError::TopologyMismatch)
        );
        assert_eq!(
            kabsch(moving_a.view(), reference_a.view(), &all(&topology_b)),
            Err(AlignmentError::SelectionTopologyMismatch)
        );
    }

    #[test]
    fn invalid_explicit_weights_are_structured_failures() {
        let points = fixture_points();
        let (topology, moving, reference) = models(&points, &points);
        let selection = all(&topology);
        let fit = |weights: &[f64]| {
            kabsch_with_options(
                moving.view(),
                reference.view(),
                &selection,
                KabschOptions {
                    weighting: AlignmentWeighting::Explicit(weights),
                    ..KabschOptions::default()
                },
            )
        };
        assert_eq!(
            fit(&[1.0]),
            Err(AlignmentError::WeightCountMismatch {
                expected: points.len(),
                actual: 1,
            })
        );
        for (weight, expected) in [
            (
                0.0,
                AlignmentError::NonPositiveWeight { selection_index: 2 },
            ),
            (
                -1.0,
                AlignmentError::NonPositiveWeight { selection_index: 2 },
            ),
            (
                f64::NAN,
                AlignmentError::NonFiniteWeight { selection_index: 2 },
            ),
            (
                f64::INFINITY,
                AlignmentError::NonFiniteWeight { selection_index: 2 },
            ),
        ] {
            let weights = [1.0, 1.0, weight, 1.0, 1.0];
            assert_eq!(fit(&weights), Err(expected));
        }
    }

    #[test]
    fn periodic_policy_rejects_by_default_and_can_use_stored_coordinates() {
        let moving_points = fixture_points();
        let expected = non_axis_transform();
        let reference_points = transform_points(&moving_points, expected);
        let (topology, mut moving, mut reference) = models(&moving_points, &reference_points);
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(3.0, 4.0, 5.0), NANOMETER),
            [true; 3],
        )
        .unwrap();
        moving.set_cell(Some(cell));
        reference.set_cell(Some(cell));
        let selection = all(&topology);
        let moving_before = moving.positions().values().value().to_vec();

        assert_eq!(
            kabsch(moving.view(), reference.view(), &selection),
            Err(AlignmentError::PeriodicCoordinates {
                moving: true,
                reference: true,
            })
        );
        let result = kabsch_with_options(
            moving.view(),
            reference.view(),
            &selection,
            KabschOptions {
                periodic_policy: PeriodicAlignmentPolicy::UseStoredCoordinates,
                ..KabschOptions::default()
            },
        )
        .unwrap();
        assert_transform_close(result.transform(), expected, 3.0e-12);
        assert_eq!(moving.cell(), Some(&cell));
        assert_eq!(reference.cell(), Some(&cell));
        assert_eq!(
            moving.positions().values().to_value(),
            moving_before.as_slice()
        );
    }

    #[test]
    fn large_coordinate_offsets_retain_small_internal_geometry() {
        let offset = Vector3::new(1.0e12, -2.0e12, 3.0e12);
        let moving = fixture_points().map(|point| point + offset);
        let rotation = non_axis_transform().rotation();
        let expected =
            RigidTransform::new(rotation, Vector3::new(-4.0e11, 7.0e11, 2.0e11)).unwrap();
        let reference = transform_points(&moving, expected);
        let (topology, moving, reference_model) = models(&moving, &reference);

        let result = kabsch(moving.view(), reference_model.view(), &all(&topology)).unwrap();

        assert_close(result.transform().rotation().determinant(), 1.0, 1.0e-12);
        for (moving, reference) in moving.positions().values().value().iter().zip(reference) {
            assert_point_close(
                result.transform().transform_point(*moving),
                reference,
                2.0e-3,
            );
        }
        assert!(result.rmsd().to_value() < 1.0e-3);
    }
}
