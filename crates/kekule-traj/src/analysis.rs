//! Fixed-topology trajectory superposition and RMSD workflows.
//!
//! [`Trajectory::superpose_to_frame`] is an explicit transactional
//! transformation. [`Trajectory::rmsd_to_frame`] measures coordinates exactly
//! as stored and never centers or fits them. [`Trajectory::aligned_rmsd_to_frame`]
//! is the allocation-light convenience for fitting and measuring in one pass.

use std::fmt;

use kekule::alignment::{kabsch_with_options, AlignmentError, KabschOptions, RigidAlignment};
use kekule::geometry::{PeriodicCell, PeriodicCellError, RigidTransform};
use kekule::structure::Positions;
use kekule::topology::AtomSelection;
use kekule::units::{Quantity, CANONICAL_LENGTH_UNIT};

use crate::{Forces, FrameError, Trajectory, TrajectoryError, TrajectoryFrame, Velocities};

/// Options used to fit every trajectory frame onto one reference frame.
///
/// This is the same fitting contract as Kekule's single-model Kabsch kernel.
/// In particular, periodic frames are rejected by default.
pub type SuperpositionOptions<'a> = KabschOptions<'a>;

/// Per-selected-atom weighting for direct RMSD measurement.
#[derive(Debug, Clone, Copy, Default)]
pub enum RmsdWeighting<'a> {
    /// Give every selected atom equal weight.
    #[default]
    Uniform,
    /// Use positive finite weights in sorted selection order.
    Explicit(&'a [f64]),
}

/// Handling of frames that carry periodic cells during direct RMSD.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum PeriodicRmsdPolicy {
    /// Reject either periodic input because no imaging has been performed.
    #[default]
    RejectPeriodic,
    /// Measure the Cartesian coordinates exactly as stored and ignore cells.
    ///
    /// This performs no imaging, wrapping, unwrapping, minimum-image
    /// correction, or molecule reconstruction.
    UseStoredCoordinates,
}

/// Options for direct RMSD over coordinates exactly as stored.
#[derive(Debug, Clone, Copy, Default)]
pub struct RmsdOptions<'a> {
    /// Per-selected-atom measurement weights.
    pub weighting: RmsdWeighting<'a>,
    /// Handling of frames carrying periodic cells.
    pub periodic_policy: PeriodicRmsdPolicy,
}

/// Options for fused superposition and RMSD measurement.
#[derive(Debug, Clone, Copy, Default)]
pub struct AlignedRmsdOptions<'a> {
    /// Options controlling the fit selection and periodic-cell policy.
    pub superposition: SuperpositionOptions<'a>,
    /// Weights for the independently chosen measurement selection.
    pub measurement_weighting: RmsdWeighting<'a>,
}

/// Complete record of one successful in-place trajectory superposition.
#[derive(Debug, Clone, PartialEq)]
pub struct SuperpositionReport {
    reference_frame: usize,
    alignments: Vec<RigidAlignment>,
}

impl SuperpositionReport {
    /// Returns the zero-based reference-frame index.
    pub const fn reference_frame(&self) -> usize {
        self.reference_frame
    }

    /// Returns one applied alignment per trajectory frame, in frame order.
    pub fn alignments(&self) -> &[RigidAlignment] {
        &self.alignments
    }

    /// Returns the alignment applied to one trajectory frame.
    pub fn alignment(&self, frame: usize) -> Option<&RigidAlignment> {
        self.alignments.get(frame)
    }

    /// Returns the number of transformed frames.
    pub fn len(&self) -> usize {
        self.alignments.len()
    }

    /// Returns whether no frames were transformed.
    pub fn is_empty(&self) -> bool {
        self.alignments.is_empty()
    }
}

impl Trajectory {
    /// Transactionally superposes every frame onto `reference_frame`.
    ///
    /// The fit uses `fit_selection`, then applies the resulting proper rigid
    /// transform to every position in the frame. Velocities, forces, and cell
    /// vectors are rotated without translation. Atom data, time, step, and
    /// frame properties are preserved. The trajectory is unchanged if any fit
    /// or transformed-frame validation fails.
    pub fn superpose_to_frame(
        &mut self,
        reference_frame: usize,
        fit_selection: &AtomSelection,
    ) -> Result<SuperpositionReport, SuperpositionError> {
        self.superpose_to_frame_with_options(
            reference_frame,
            fit_selection,
            SuperpositionOptions::default(),
        )
    }

    /// Transactionally superposes every frame with explicit fitting options.
    pub fn superpose_to_frame_with_options(
        &mut self,
        reference_frame: usize,
        fit_selection: &AtomSelection,
        options: SuperpositionOptions<'_>,
    ) -> Result<SuperpositionReport, SuperpositionError> {
        let alignments = {
            let reference = self.frames().nth(reference_frame).ok_or(
                SuperpositionError::ReferenceFrameOutOfRange {
                    index: reference_frame,
                    frame_count: self.len(),
                },
            )?;
            self.frames()
                .enumerate()
                .map(|(frame, moving)| {
                    kabsch_with_options(
                        moving.model_view(),
                        reference.model_view(),
                        fit_selection,
                        options,
                    )
                    .map_err(|source| SuperpositionError::Alignment { frame, source })
                })
                .collect::<Result<Vec<_>, _>>()?
        };

        let transformed_frames = self
            .frames()
            .zip(&alignments)
            .enumerate()
            .map(|(frame, (source, alignment))| {
                transform_frame(source, alignment.transform()).map_err(|error| match error {
                    TransformFrameError::Frame(source) => {
                        SuperpositionError::FrameTransform { frame, source }
                    }
                    TransformFrameError::Cell(source) => {
                        SuperpositionError::CellTransform { frame, source }
                    }
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let transformed = Trajectory::from_frames(self.shared_topology(), transformed_frames)
            .map_err(|source| SuperpositionError::TrajectoryPublication(Box::new(source)))?;
        *self = transformed;

        Ok(SuperpositionReport {
            reference_frame,
            alignments,
        })
    }

    /// Computes direct RMSD from every frame to `reference_frame`.
    ///
    /// Coordinates are measured exactly as stored. This method never centers,
    /// rotates, translates, or otherwise changes the trajectory.
    pub fn rmsd_to_frame(
        &self,
        reference_frame: usize,
        selection: &AtomSelection,
    ) -> Result<Quantity<Vec<f64>>, RmsdError> {
        self.rmsd_to_frame_with_options(reference_frame, selection, RmsdOptions::default())
    }

    /// Computes direct RMSD with explicit weighting and periodic-cell policy.
    pub fn rmsd_to_frame_with_options(
        &self,
        reference_frame: usize,
        selection: &AtomSelection,
        options: RmsdOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, RmsdError> {
        let reference =
            self.frames()
                .nth(reference_frame)
                .ok_or(RmsdError::ReferenceFrameOutOfRange {
                    index: reference_frame,
                    frame_count: self.len(),
                })?;
        validate_measurement_selection(self, selection)?;
        let weights = NormalizedRmsdWeights::new(options.weighting, selection.indices().len())?;
        let reference_periodic = reference.cell().is_some();
        let mut values = Vec::with_capacity(self.len());
        for (frame, moving) in self.frames().enumerate() {
            let moving_periodic = moving.cell().is_some();
            if options.periodic_policy == PeriodicRmsdPolicy::RejectPeriodic
                && (moving_periodic || reference_periodic)
            {
                return Err(RmsdError::PeriodicCoordinates {
                    frame,
                    moving: moving_periodic,
                    reference: reference_periodic,
                });
            }
            values.push(measure_rmsd(
                moving, reference, selection, weights, None, frame,
            )?);
        }
        Ok(Quantity::new(values, CANONICAL_LENGTH_UNIT))
    }

    /// Fits every frame and measures RMSD without changing the trajectory.
    ///
    /// The transform is determined by `fit_selection`, while RMSD is measured
    /// over `measurement_selection`. This distinction permits workflows such
    /// as fitting a protein backbone and measuring ligand or domain motion.
    pub fn aligned_rmsd_to_frame(
        &self,
        reference_frame: usize,
        fit_selection: &AtomSelection,
        measurement_selection: &AtomSelection,
    ) -> Result<Quantity<Vec<f64>>, RmsdError> {
        self.aligned_rmsd_to_frame_with_options(
            reference_frame,
            fit_selection,
            measurement_selection,
            AlignedRmsdOptions::default(),
        )
    }

    /// Fits and measures RMSD with explicit fit and measurement options.
    pub fn aligned_rmsd_to_frame_with_options(
        &self,
        reference_frame: usize,
        fit_selection: &AtomSelection,
        measurement_selection: &AtomSelection,
        options: AlignedRmsdOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, RmsdError> {
        let reference =
            self.frames()
                .nth(reference_frame)
                .ok_or(RmsdError::ReferenceFrameOutOfRange {
                    index: reference_frame,
                    frame_count: self.len(),
                })?;
        validate_measurement_selection(self, measurement_selection)?;
        let weights = NormalizedRmsdWeights::new(
            options.measurement_weighting,
            measurement_selection.indices().len(),
        )?;
        let mut values = Vec::with_capacity(self.len());
        for (frame, moving) in self.frames().enumerate() {
            let alignment = kabsch_with_options(
                moving.model_view(),
                reference.model_view(),
                fit_selection,
                options.superposition,
            )
            .map_err(|source| RmsdError::Alignment { frame, source })?;
            values.push(measure_rmsd(
                moving,
                reference,
                measurement_selection,
                weights,
                Some(alignment.transform()),
                frame,
            )?);
        }
        Ok(Quantity::new(values, CANONICAL_LENGTH_UNIT))
    }
}

fn validate_measurement_selection(
    trajectory: &Trajectory,
    selection: &AtomSelection,
) -> Result<(), RmsdError> {
    if !std::ptr::eq(selection.topology(), trajectory.topology()) {
        return Err(RmsdError::SelectionTopologyMismatch);
    }
    if selection.indices().is_empty() {
        return Err(RmsdError::EmptySelection);
    }
    Ok(())
}

fn measure_rmsd(
    moving: crate::TrajectoryFrameView<'_>,
    reference: crate::TrajectoryFrameView<'_>,
    selection: &AtomSelection,
    weights: NormalizedRmsdWeights<'_>,
    transform: Option<RigidTransform>,
    frame: usize,
) -> Result<f64, RmsdError> {
    let moving_positions = moving.positions().values();
    let reference_positions = reference.positions().values();
    let moving_positions = moving_positions.value();
    let reference_positions = reference_positions.value();
    let mut squared_residual = CompensatedSum::default();
    let mut weight_sum = CompensatedSum::default();
    for (selection_index, dense_index) in selection.indices().iter().copied().enumerate() {
        let mut point = moving_positions[dense_index.index()];
        if let Some(transform) = transform {
            point = transform.transform_point(point);
        }
        let residual = point - reference_positions[dense_index.index()];
        let weight = weights.at(selection_index);
        squared_residual.add(weight * residual.norm_squared());
        weight_sum.add(weight);
    }
    let mean_squared_residual = squared_residual.value() / weight_sum.value();
    if !mean_squared_residual.is_finite() || mean_squared_residual < 0.0 {
        return Err(RmsdError::NumericalFailure { frame });
    }
    Ok(mean_squared_residual.sqrt())
}

fn transform_frame(
    source: crate::TrajectoryFrameView<'_>,
    transform: RigidTransform,
) -> Result<TrajectoryFrame, TransformFrameError> {
    let topology = source.topology_arc();
    let positions = source.positions().values();
    let positions = positions
        .value()
        .iter()
        .copied()
        .map(|point| transform.transform_point(point))
        .collect::<Vec<_>>();
    let positions = Positions::new(Quantity::new(positions, CANONICAL_LENGTH_UNIT))
        .map_err(FrameError::from)
        .map_err(|source| TransformFrameError::Frame(Box::new(source)))?;
    let cell = source
        .cell()
        .copied()
        .map(|cell| transform_cell(cell, transform))
        .transpose()
        .map_err(|source| TransformFrameError::Cell(Box::new(source)))?;
    let mut transformed = TrajectoryFrame::new(positions, topology.bond_count());
    transformed.set_cell(cell);
    transformed
        .set_properties(source.properties().clone())
        .map_err(|source| TransformFrameError::Frame(Box::new(source)))?;
    if let Some(values) = source.velocities() {
        let rotated = values
            .value()
            .iter()
            .copied()
            .map(|value| transform.transform_vector(value))
            .collect::<Vec<_>>();
        let velocities = Velocities::new(Quantity::new(rotated, values.unit()))
            .map_err(|source| TransformFrameError::Frame(Box::new(source)))?;
        transformed
            .set_velocities(Some(velocities))
            .map_err(|source| TransformFrameError::Frame(Box::new(source)))?;
    }
    if let Some(values) = source.forces() {
        let rotated = values
            .value()
            .iter()
            .copied()
            .map(|value| transform.transform_vector(value))
            .collect::<Vec<_>>();
        let forces = Forces::new(Quantity::new(rotated, values.unit()))
            .map_err(|source| TransformFrameError::Frame(Box::new(source)))?;
        transformed
            .set_forces(Some(forces))
            .map_err(|source| TransformFrameError::Frame(Box::new(source)))?;
    }
    transformed
        .set_time(source.time())
        .map_err(|source| TransformFrameError::Frame(Box::new(source)))?;
    transformed.set_step(source.step());
    Ok(transformed)
}

fn transform_cell(
    cell: PeriodicCell,
    transform: RigidTransform,
) -> Result<PeriodicCell, PeriodicCellError> {
    let vectors = cell
        .vectors()
        .map(|vectors| vectors.map(|vector| transform.transform_vector(vector)));
    PeriodicCell::new(vectors, cell.periodic_axes())
}

#[derive(Debug)]
enum TransformFrameError {
    Frame(Box<FrameError>),
    Cell(Box<PeriodicCellError>),
}

#[derive(Debug, Clone, Copy)]
enum NormalizedRmsdWeights<'a> {
    Uniform,
    Explicit { values: &'a [f64], maximum: f64 },
}

impl<'a> NormalizedRmsdWeights<'a> {
    fn new(weighting: RmsdWeighting<'a>, selected: usize) -> Result<Self, RmsdError> {
        match weighting {
            RmsdWeighting::Uniform => Ok(Self::Uniform),
            RmsdWeighting::Explicit(values) => {
                if values.len() != selected {
                    return Err(RmsdError::WeightCountMismatch {
                        expected: selected,
                        actual: values.len(),
                    });
                }
                let mut maximum: f64 = 0.0;
                for (selection_index, weight) in values.iter().copied().enumerate() {
                    if !weight.is_finite() {
                        return Err(RmsdError::NonFiniteWeight { selection_index });
                    }
                    if weight <= 0.0 {
                        return Err(RmsdError::NonPositiveWeight { selection_index });
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

#[derive(Debug, Default)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    fn add(&mut self, value: f64) {
        let updated = self.sum + value;
        if self.sum.abs() >= value.abs() {
            self.correction += (self.sum - updated) + value;
        } else {
            self.correction += (value - updated) + self.sum;
        }
        self.sum = updated;
    }

    fn value(&self) -> f64 {
        self.sum + self.correction
    }
}

/// Failure to transform a complete in-memory trajectory by rigid superposition.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum SuperpositionError {
    /// The requested reference frame does not exist.
    ReferenceFrameOutOfRange { index: usize, frame_count: usize },
    /// One frame could not be fitted to the reference.
    Alignment {
        frame: usize,
        source: AlignmentError,
    },
    /// One frame's positions or dynamic vectors could not be transformed.
    FrameTransform {
        frame: usize,
        source: Box<FrameError>,
    },
    /// One periodic cell could not be rotated into a valid cell.
    CellTransform {
        frame: usize,
        source: Box<PeriodicCellError>,
    },
    /// Valid transformed frames could not be published as one trajectory.
    TrajectoryPublication(Box<TrajectoryError>),
}

impl fmt::Display for SuperpositionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ReferenceFrameOutOfRange { index, frame_count } => write!(
                formatter,
                "trajectory reference frame {index} is out of range for {frame_count} frames"
            ),
            Self::Alignment { frame, source } => {
                write!(
                    formatter,
                    "cannot superpose trajectory frame {frame}: {source}"
                )
            }
            Self::FrameTransform { frame, source } => write!(
                formatter,
                "cannot publish transformed trajectory frame {frame}: {source}"
            ),
            Self::CellTransform { frame, source } => write!(
                formatter,
                "cannot rotate periodic cell for trajectory frame {frame}: {source}"
            ),
            Self::TrajectoryPublication(source) => {
                write!(formatter, "cannot publish superposed trajectory: {source}")
            }
        }
    }
}

impl std::error::Error for SuperpositionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Alignment { source, .. } => Some(source),
            Self::FrameTransform { source, .. } => Some(source.as_ref()),
            Self::CellTransform { source, .. } => Some(source.as_ref()),
            Self::TrajectoryPublication(source) => Some(source.as_ref()),
            Self::ReferenceFrameOutOfRange { .. } => None,
        }
    }
}

/// Failure to calculate a direct or aligned trajectory RMSD series.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum RmsdError {
    /// The requested reference frame does not exist.
    ReferenceFrameOutOfRange { index: usize, frame_count: usize },
    /// The measurement selection belongs to another exact topology.
    SelectionTopologyMismatch,
    /// Direct RMSD requires at least one selected atom.
    EmptySelection,
    /// Explicit measurement weights do not match the selected atom count.
    WeightCountMismatch { expected: usize, actual: usize },
    /// One explicit measurement weight is NaN or infinite.
    NonFiniteWeight { selection_index: usize },
    /// One explicit measurement weight is zero or negative.
    NonPositiveWeight { selection_index: usize },
    /// Stored-coordinate RMSD was not explicitly permitted for periodic data.
    PeriodicCoordinates {
        frame: usize,
        moving: bool,
        reference: bool,
    },
    /// One frame could not be fitted before RMSD measurement.
    Alignment {
        frame: usize,
        source: AlignmentError,
    },
    /// Finite inputs did not produce a finite RMSD for one frame.
    NumericalFailure { frame: usize },
}

impl fmt::Display for RmsdError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ReferenceFrameOutOfRange { index, frame_count } => write!(
                formatter,
                "trajectory reference frame {index} is out of range for {frame_count} frames"
            ),
            Self::SelectionTopologyMismatch => {
                formatter.write_str("RMSD selection belongs to another topology allocation")
            }
            Self::EmptySelection => formatter.write_str("RMSD requires at least one selected atom"),
            Self::WeightCountMismatch { expected, actual } => write!(
                formatter,
                "RMSD requires {expected} weights, but received {actual}"
            ),
            Self::NonFiniteWeight { selection_index } => write!(
                formatter,
                "RMSD weight at selection index {selection_index} is not finite"
            ),
            Self::NonPositiveWeight { selection_index } => write!(
                formatter,
                "RMSD weight at selection index {selection_index} is not strictly positive"
            ),
            Self::PeriodicCoordinates {
                frame,
                moving,
                reference,
            } => write!(
                formatter,
                "direct RMSD for trajectory frame {frame} rejects periodic coordinates (moving cell: {moving}, reference cell: {reference})"
            ),
            Self::Alignment { frame, source } => {
                write!(formatter, "cannot align trajectory frame {frame} for RMSD: {source}")
            }
            Self::NumericalFailure { frame } => {
                write!(formatter, "RMSD numerical calculation failed for frame {frame}")
            }
        }
    }
}

impl std::error::Error for RmsdError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Alignment { source, .. } => Some(source),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use kekule::alignment::{AlignmentGeometry, PeriodicAlignmentPolicy};
    use kekule::core::{Atom, BondOrder, Element};
    use kekule::geometry::{Matrix3, Point3, Vector3};
    use kekule::properties::{PropertyKey, PropertyValue};
    use kekule::topology::{Topology, TopologyBuilder};
    use kekule::units::{
        Quantity, ANGSTROM, CANONICAL_FORCE_UNIT, CANONICAL_VELOCITY_UNIT, DIMENSIONLESS,
        NANOMETER, PICOSECOND,
    };

    fn make_topology(atom_count: usize) -> Arc<Topology> {
        let mut graph = kekule::core::MoleculeEditor::new();
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

    fn selection(topology: &Arc<Topology>, indices: &[usize]) -> AtomSelection {
        AtomSelection::from_atoms(
            topology,
            indices
                .iter()
                .copied()
                .map(|index| topology.atom_ids()[index]),
        )
        .unwrap()
    }

    fn all(topology: &Arc<Topology>) -> AtomSelection {
        AtomSelection::from_atoms(topology, topology.atom_ids().iter().copied()).unwrap()
    }

    fn frame(topology: &Arc<Topology>, points: &[Point3]) -> TrajectoryFrame {
        TrajectoryFrame::new(
            Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
            topology.bond_count(),
        )
    }

    fn transformed(points: &[Point3], transform: RigidTransform) -> Vec<Point3> {
        points
            .iter()
            .copied()
            .map(|point| transform.transform_point(point))
            .collect()
    }

    fn quarter_turn() -> RigidTransform {
        RigidTransform::new(
            Matrix3::from_columns(
                Vector3::new(0.0, 1.0, 0.0),
                Vector3::new(-1.0, 0.0, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ),
            Vector3::new(4.0, -2.0, 1.0),
        )
        .unwrap()
    }

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected}, received {actual}"
        );
    }

    fn assert_point_close(actual: Point3, expected: Point3, tolerance: f64) {
        assert_close(actual.x, expected.x, tolerance);
        assert_close(actual.y, expected.y, tolerance);
        assert_close(actual.z, expected.z, tolerance);
    }

    fn assert_vector_close(actual: Vector3, expected: Vector3, tolerance: f64) {
        assert_close(actual.x, expected.x, tolerance);
        assert_close(actual.y, expected.y, tolerance);
        assert_close(actual.z, expected.z, tolerance);
    }

    #[test]
    fn direct_rmsd_does_not_fit_and_explicit_weights_follow_selection_order() {
        let topology = make_topology(2);
        let reference = [Point3::origin(), Point3::origin()];
        let moving = [Point3::new(1.0, 0.0, 0.0), Point3::new(3.0, 0.0, 0.0)];
        let trajectory = Trajectory::from_frames(
            Arc::clone(&topology),
            [frame(&topology, &reference), frame(&topology, &moving)],
        )
        .unwrap();
        let selection = all(&topology);

        let uniform = trajectory.rmsd_to_frame(0, &selection).unwrap();
        assert_eq!(uniform.unit(), CANONICAL_LENGTH_UNIT);
        assert_close(uniform.value()[0], 0.0, 1.0e-14);
        assert_close(uniform.value()[1], 5.0_f64.sqrt(), 1.0e-14);

        let weights = [3.0, 1.0];
        let weighted = trajectory
            .rmsd_to_frame_with_options(
                0,
                &selection,
                RmsdOptions {
                    weighting: RmsdWeighting::Explicit(&weights),
                    ..RmsdOptions::default()
                },
            )
            .unwrap();
        assert_close(weighted.value()[1], 3.0_f64.sqrt(), 1.0e-14);
        let scaled_weights = [30.0, 10.0];
        let scaled = trajectory
            .rmsd_to_frame_with_options(
                0,
                &selection,
                RmsdOptions {
                    weighting: RmsdWeighting::Explicit(&scaled_weights),
                    ..RmsdOptions::default()
                },
            )
            .unwrap();
        assert_eq!(scaled, weighted);
    }

    #[test]
    fn split_and_fused_fit_measure_workflows_agree_for_distinct_selections() {
        let topology = make_topology(4);
        let moving = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 1.5, 0.0),
            Point3::new(0.5, 0.4, 1.2),
        ];
        let transform = quarter_turn();
        let mut reference = transformed(&moving, transform);
        reference[3].x += 2.0;
        let trajectory = Trajectory::from_frames(
            Arc::clone(&topology),
            [frame(&topology, &reference), frame(&topology, &moving)],
        )
        .unwrap();
        let fit_selection = selection(&topology, &[0, 1, 2]);
        let measurement_selection = selection(&topology, &[3]);

        let fused = trajectory
            .aligned_rmsd_to_frame(0, &fit_selection, &measurement_selection)
            .unwrap();
        assert_close(fused.value()[0], 0.0, 1.0e-12);
        assert_close(fused.value()[1], 2.0, 2.0e-12);

        let mut split = trajectory.clone();
        let report = split.superpose_to_frame(0, &fit_selection).unwrap();
        assert_eq!(report.reference_frame(), 0);
        assert_eq!(report.len(), 2);
        let measured = split.rmsd_to_frame(0, &measurement_selection).unwrap();
        assert_close(measured.value()[0], fused.value()[0], 1.0e-12);
        assert_close(measured.value()[1], fused.value()[1], 2.0e-12);
    }

    #[test]
    fn superposition_rotates_complete_geometric_state_and_preserves_metadata() {
        let topology = make_topology(4);
        let moving = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 1.5, 0.0),
            Point3::new(0.5, 0.4, 1.2),
        ];
        let transform = quarter_turn();
        let reference = transformed(&moving, transform);
        let moving_cell = PeriodicCell::new(
            Quantity::new(
                [
                    Vector3::new(10.0, 0.0, 0.0),
                    Vector3::new(0.0, 20.0, 0.0),
                    Vector3::new(0.0, 0.0, 30.0),
                ],
                ANGSTROM,
            ),
            [true; 3],
        )
        .unwrap();
        let reference_cell = transform_cell(moving_cell, transform).unwrap();
        let mut reference_frame = TrajectoryFrame::new(
            Positions::new(Quantity::new(&reference, NANOMETER)).unwrap(),
            topology.bond_count(),
        );
        reference_frame.set_cell(Some(reference_cell));
        reference_frame
            .set_velocities(Some(
                Velocities::new(Quantity::new(
                    vec![Vector3::new(0.0, 1.0, 0.0); 4],
                    CANONICAL_VELOCITY_UNIT,
                ))
                .unwrap(),
            ))
            .unwrap();
        reference_frame
            .set_forces(Some(
                Forces::new(Quantity::new(
                    vec![Vector3::new(-1.0, 0.0, 0.0); 4],
                    CANONICAL_FORCE_UNIT,
                ))
                .unwrap(),
            ))
            .unwrap();

        let mut moving_frame = TrajectoryFrame::new(
            Positions::new(Quantity::new(&moving, NANOMETER)).unwrap(),
            topology.bond_count(),
        );
        moving_frame.set_cell(Some(moving_cell));
        moving_frame
            .set_velocities(Some(
                Velocities::new(Quantity::new(
                    vec![Vector3::new(1.0, 0.0, 0.0); 4],
                    CANONICAL_VELOCITY_UNIT,
                ))
                .unwrap(),
            ))
            .unwrap();
        moving_frame
            .set_forces(Some(
                Forces::new(Quantity::new(
                    vec![Vector3::new(0.0, 1.0, 0.0); 4],
                    CANONICAL_FORCE_UNIT,
                ))
                .unwrap(),
            ))
            .unwrap();
        moving_frame
            .set_time(Some(Quantity::new(2.5, PICOSECOND)))
            .unwrap();
        moving_frame.set_step(Some(25));
        let score_key = PropertyKey::new("score").unwrap();
        moving_frame
            .properties_mut()
            .atoms_mut()
            .set_value(
                score_key.clone(),
                0,
                Some(PropertyValue::Real {
                    value: 0.7,
                    unit: DIMENSIONLESS,
                }),
            )
            .unwrap();
        let label_key = PropertyKey::new("label").unwrap();
        moving_frame
            .properties_mut()
            .insert(
                label_key.clone(),
                PropertyValue::String("moving".to_owned()),
            )
            .unwrap();
        let expected_properties = moving_frame.properties().clone();
        let mut trajectory =
            Trajectory::from_frames(Arc::clone(&topology), [reference_frame, moving_frame])
                .unwrap();

        let report = trajectory
            .superpose_to_frame_with_options(
                0,
                &all(&topology),
                SuperpositionOptions {
                    periodic_policy: PeriodicAlignmentPolicy::UseStoredCoordinates,
                    ..SuperpositionOptions::default()
                },
            )
            .unwrap();
        assert_eq!(report.alignments().len(), 2);
        let transformed = trajectory.frame(1).unwrap();
        for (actual, expected) in transformed
            .positions()
            .values()
            .value()
            .iter()
            .copied()
            .zip(reference)
        {
            assert_point_close(actual, expected, 3.0e-12);
        }
        let transformed_cell = transformed.cell().copied().unwrap();
        assert_eq!(
            transformed_cell.periodic_axes(),
            reference_cell.periodic_axes()
        );
        for (actual, expected) in transformed_cell
            .vectors()
            .to_value()
            .into_iter()
            .zip(reference_cell.vectors().to_value())
        {
            assert_vector_close(actual, expected, 1.0e-12);
        }
        for velocity in transformed.velocities().unwrap().values().to_value() {
            assert_vector_close(*velocity, Vector3::new(0.0, 1.0, 0.0), 2.0e-12);
        }
        for force in transformed.forces().unwrap().values().to_value() {
            assert_vector_close(*force, Vector3::new(-1.0, 0.0, 0.0), 2.0e-12);
        }
        assert_eq!(transformed.time(), Some(Quantity::new(2.5, PICOSECOND)));
        assert_eq!(transformed.step(), Some(25));
        assert_eq!(transformed.properties(), &expected_properties);
        assert_eq!(
            transformed.properties().get(&label_key),
            Some(&PropertyValue::String("moving".to_owned()))
        );
    }

    #[test]
    fn late_superposition_failure_leaves_complete_trajectory_unchanged() {
        let topology = make_topology(4);
        let reference = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 1.5, 0.0),
            Point3::new(0.5, 0.4, 1.2),
        ];
        let moving = transformed(&reference, quarter_turn());
        let collinear = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ];
        let mut trajectory = Trajectory::from_frames(
            Arc::clone(&topology),
            [
                frame(&topology, &reference),
                frame(&topology, &moving),
                frame(&topology, &collinear),
            ],
        )
        .unwrap();
        let before = (0..trajectory.len())
            .map(|frame| trajectory.frame(frame).unwrap().clone())
            .collect::<Vec<_>>();

        assert_eq!(
            trajectory.superpose_to_frame(3, &all(&topology)),
            Err(SuperpositionError::ReferenceFrameOutOfRange {
                index: 3,
                frame_count: 3,
            })
        );
        assert_eq!(
            trajectory.superpose_to_frame(0, &all(&topology)),
            Err(SuperpositionError::Alignment {
                frame: 2,
                source: AlignmentError::DegenerateGeometry {
                    geometry: AlignmentGeometry::Moving,
                },
            })
        );
        let after = (0..trajectory.len())
            .map(|frame| trajectory.frame(frame).unwrap().clone())
            .collect::<Vec<_>>();
        assert_eq!(after, before);
    }

    #[test]
    fn rmsd_reports_reference_selection_weight_and_periodic_failures() {
        let topology = make_topology(4);
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
            [true; 3],
        )
        .unwrap();
        let mut periodic_frame = TrajectoryFrame::new(
            Positions::new(Quantity::new(&points, ANGSTROM)).unwrap(),
            topology.bond_count(),
        );
        periodic_frame.set_cell(Some(cell));
        let trajectory = Trajectory::from_frames(Arc::clone(&topology), [periodic_frame]).unwrap();

        assert_eq!(
            trajectory.rmsd_to_frame(1, &all(&topology)),
            Err(RmsdError::ReferenceFrameOutOfRange {
                index: 1,
                frame_count: 1,
            })
        );
        let empty = selection(&topology, &[]);
        assert_eq!(
            trajectory.rmsd_to_frame(0, &empty),
            Err(RmsdError::EmptySelection)
        );
        let independent = make_topology(4);
        assert_eq!(
            trajectory.rmsd_to_frame(0, &all(&independent)),
            Err(RmsdError::SelectionTopologyMismatch)
        );
        let weights = [1.0];
        assert_eq!(
            trajectory.rmsd_to_frame_with_options(
                0,
                &all(&topology),
                RmsdOptions {
                    weighting: RmsdWeighting::Explicit(&weights),
                    periodic_policy: PeriodicRmsdPolicy::UseStoredCoordinates,
                },
            ),
            Err(RmsdError::WeightCountMismatch {
                expected: 4,
                actual: 1,
            })
        );
        for (weight, expected) in [
            (f64::NAN, RmsdError::NonFiniteWeight { selection_index: 2 }),
            (0.0, RmsdError::NonPositiveWeight { selection_index: 2 }),
        ] {
            let weights = [1.0, 1.0, weight, 1.0];
            assert_eq!(
                trajectory.rmsd_to_frame_with_options(
                    0,
                    &all(&topology),
                    RmsdOptions {
                        weighting: RmsdWeighting::Explicit(&weights),
                        periodic_policy: PeriodicRmsdPolicy::UseStoredCoordinates,
                    },
                ),
                Err(expected)
            );
        }
        assert_eq!(
            trajectory.rmsd_to_frame(0, &all(&topology)),
            Err(RmsdError::PeriodicCoordinates {
                frame: 0,
                moving: true,
                reference: true,
            })
        );
        let allowed = trajectory
            .rmsd_to_frame_with_options(
                0,
                &all(&topology),
                RmsdOptions {
                    periodic_policy: PeriodicRmsdPolicy::UseStoredCoordinates,
                    ..RmsdOptions::default()
                },
            )
            .unwrap();
        assert_eq!(allowed.value(), &[0.0]);
    }

    #[test]
    fn aligned_rmsd_reports_the_failing_frame_without_mutating_input() {
        let topology = make_topology(4);
        let reference = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 1.5, 0.0),
            Point3::new(0.5, 0.4, 1.2),
        ];
        let collinear = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ];
        let trajectory = Trajectory::from_frames(
            Arc::clone(&topology),
            [frame(&topology, &reference), frame(&topology, &collinear)],
        )
        .unwrap();
        let before = trajectory
            .frame(1)
            .unwrap()
            .positions()
            .values()
            .value()
            .to_vec();

        assert_eq!(
            trajectory.aligned_rmsd_to_frame(0, &all(&topology), &all(&topology)),
            Err(RmsdError::Alignment {
                frame: 1,
                source: AlignmentError::DegenerateGeometry {
                    geometry: AlignmentGeometry::Moving,
                },
            })
        );
        assert_eq!(
            trajectory.frame(1).unwrap().positions().values().to_value(),
            before.as_slice()
        );
    }
}
