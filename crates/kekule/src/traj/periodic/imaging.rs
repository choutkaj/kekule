//! Whole-molecule imaging and trajectory unwrapping over periodic cells, and
//! the periodic error type.

use std::collections::VecDeque;
use std::fmt;

use crate::math::{PeriodicCell, Point3, Vector3};
use crate::structure::{
    AsStructureView, DenseArrayError, Place, Positions, StructureSet, StructureSetError,
};
use crate::topology::{AtomSelection, InstanceBondId, Topology, TopologyAtomIndex};
use crate::traj::FrameView;
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};

use super::{MoleculeImager, TrajectoryUnwrapper};

/// Failure to reconstruct or publish periodic coordinates.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PeriodicError {
    /// A frame has no periodic cell.
    MissingCell {
        /// The frame index.
        frame: usize,
    },
    /// The anchor selection belongs to another topology.
    SelectionTopologyMismatch,
    /// The reference structure belongs to another topology.
    ReferenceTopologyMismatch,
    /// A frame belongs to another topology.
    TopologyMismatch {
        /// The frame index.
        frame: usize,
    },
    /// Stateful unwrapping requires consecutive source-frame indices.
    NonSequentialFrame {
        /// The previously unwrapped frame index.
        previous: usize,
        /// The frame index that does not follow it.
        frame: usize,
    },
    /// Available time values must not decrease, including across missing values.
    NonMonotonicTime {
        /// The frame index whose time decreases.
        frame: usize,
    },
    /// Imaging was requested with no anchor atom.
    EmptyAnchors,
    /// Selected shortest bond images do not form a consistent finite molecule (e.g. a winding ring).
    InconsistentBondImages {
        /// The frame index.
        frame: usize,
        /// The bond whose images are inconsistent.
        bond: InstanceBondId,
    },
    /// The periodic axes of the cell change between unwrapped frames.
    PeriodicAxesChanged {
        /// The frame index whose axes differ.
        frame: usize,
    },
    /// A displacement lies on a half-cell boundary, with no unique temporal image.
    AmbiguousDisplacement {
        /// The frame index.
        frame: usize,
        /// The atom whose displacement is ambiguous.
        atom: TopologyAtomIndex,
        /// The cell axis (0, 1 or 2) of the ambiguity.
        axis: usize,
    },
    /// An ill-conditioned cell requires more lattice images than the bounded
    /// nearest-image search examines.
    LimitExceeded {
        /// The frame whose cell is ill-conditioned.
        frame: usize,
        /// The exceeded bound of lattice images.
        limit: crate::diagnostics::LimitExceeded,
    },
    /// A periodic coordinate calculation exceeds numerical precision or
    /// range.
    NumericalFailure {
        /// The frame index.
        frame: usize,
    },
    /// A frame's positions are invalid; see the source.
    Position {
        /// The frame index.
        frame: usize,
        /// The position failure.
        source: Box<DenseArrayError>,
    },
    /// Valid positions could not be published into the destination.
    Publication(Box<StructureSetError>),
}

impl fmt::Display for PeriodicError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingCell { frame } => write!(f, "Frame {frame} has no periodic cell."),
            Self::SelectionTopologyMismatch => {
                f.write_str("The anchor selection belongs to another topology.")
            }
            Self::ReferenceTopologyMismatch => {
                f.write_str("The reference structure belongs to another topology.")
            }
            Self::TopologyMismatch { frame } => {
                write!(f, "Frame {frame} belongs to another topology.")
            }
            Self::NonSequentialFrame { previous, frame } => write!(
                f,
                "Frame {frame} does not follow frame {previous}."
            ),
            Self::NonMonotonicTime { frame } => {
                write!(f, "Time decreases at frame {frame}.")
            }
            Self::EmptyAnchors => f.write_str("No anchor atom is selected."),
            Self::InconsistentBondImages { frame, bond } => write!(
                f,
                "The periodic images around {bond} in frame {frame} are inconsistent."
            ),
            Self::PeriodicAxesChanged { frame } => {
                write!(f, "The periodic axes change at frame {frame}.")
            }
            Self::AmbiguousDisplacement { frame, atom, axis } => write!(
                f,
                "The displacement of {atom} along cell axis {} in frame {frame} is half a cell. Its image is ambiguous.",
                ["a", "b", "c"].get(*axis).copied().unwrap_or("?")
            ),
            Self::LimitExceeded { frame, limit } => write!(
                f,
                "The nearest-image search in frame {frame} exceeded its {}.",
                limit.describe()
            ),
            Self::NumericalFailure { frame } => write!(
                f,
                "A periodic coordinate in frame {frame} exceeds the precision or range of double-precision numbers."
            ),
            Self::Position { frame, .. } => write!(f, "Frame {frame} has invalid positions."),
            Self::Publication(_) => f.write_str("Cannot publish the periodic coordinates."),
        }
    }
}

impl std::error::Error for PeriodicError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Position { source, .. } => Some(source.as_ref()),
            Self::Publication(source) => Some(source.as_ref()),
            _ => None,
        }
    }
}

/// Makes each bonded molecule whole independently in every frame.
///
/// The first atom of each molecule stays at its original position; other atoms
/// are placed using shortest Cartesian bond images. All bonds, including ring
/// closures, must agree with the reconstruction. Molecules are neither centered
/// nor joined to other molecules. This does not establish continuity across time.
///
/// # Errors
///
/// Returns [`PeriodicError::MissingCell`] for a frame without a cell,
/// [`PeriodicError::InconsistentBondImages`] when bond images disagree,
/// [`PeriodicError::LimitExceeded`] or [`PeriodicError::NumericalFailure`] for
/// an ill-conditioned cell or out-of-range coordinates, and
/// [`PeriodicError::Position`] or [`PeriodicError::Publication`] when the
/// result cannot be stored. On error the set is unchanged.
pub fn make_molecules_whole(trajectory: &mut StructureSet) -> Result<(), PeriodicError> {
    let positions = whole_positions(trajectory, None)?;
    publish(trajectory, positions)
}

/// Images whole molecules around explicitly selected anchors in every frame.
///
/// Every molecule containing a selected atom becomes an anchor. Other anchors
/// are placed at their nearest centroid image to the first anchor in topology
/// order. The combined geometric center of all atoms in the anchor molecules
/// is moved to the cell center along periodic axes. Other whole molecules are
/// placed at their nearest centroid image to that combined center. No anchor
/// heuristics are used.
///
/// This includes making molecules whole. It is a per-frame centering operation,
/// not temporal unwrapping. Coordinates along nonperiodic fractional axes remain
/// unchanged; a molecule may extend outside the primary cell.
///
/// # Errors
///
/// Returns [`PeriodicError::SelectionTopologyMismatch`] when `anchors`
/// belongs to another layout, [`PeriodicError::EmptyAnchors`] when it selects
/// no atom, and the errors of [`make_molecules_whole`]. On error the set is
/// unchanged.
pub fn image_molecules(
    trajectory: &mut StructureSet,
    anchors: &AtomSelection,
) -> Result<(), PeriodicError> {
    let positions = whole_positions(trajectory, Some(anchors))?;
    publish(trajectory, positions)
}

/// Moves each molecule of every frame to its periodic image nearest the same
/// molecule in `reference`.
///
/// Each bonded molecule of a frame is translated as a whole by the lattice
/// vector of the frame's own cell that brings its geometric center nearest,
/// by the exact Cartesian nearest-image search, to the center of the same
/// molecule in `reference`. Positions change only by lattice vectors, so each
/// frame stays the same periodic system; the frame's cell, velocities, forces,
/// time, step, and properties are retained.
///
/// Molecules rather than atoms are imaged, as in [`image_molecules`], so that
/// molecules made whole by earlier steps stay whole: per-atom images would
/// split a molecule whose atoms lie on either side of a half-cell distance from
/// the reference. Make molecules whole in both structures first (see
/// [`make_molecules_whole`]), or the centers are not meaningful; this function
/// does not repair split molecules.
///
/// The frame's cell defines the images. The reference's cell is not used, so
/// the reference may have a different cell or none; with a changing cell, the
/// nearest image is the nearest under the frame's cell. A frame without a cell
/// is an error. Ties between equally near images are resolved
/// deterministically. This is a per-frame operation, independent of frame
/// order and unlike [`unwrap_trajectory`] free of temporal state. To image onto
/// a frame of the same set, retain a copy of it first, for example with
/// `trajectory.get(0).unwrap().to_structure()`.
///
/// # Errors
///
/// Returns [`PeriodicError::ReferenceTopologyMismatch`] when `reference` does
/// not share the set's layout, [`PeriodicError::MissingCell`] for a
/// frame without a cell, [`PeriodicError::LimitExceeded`] or
/// [`PeriodicError::NumericalFailure`] for an ill-conditioned cell or
/// out-of-range coordinates, and [`PeriodicError::Position`] or
/// [`PeriodicError::Publication`] when the result cannot be stored. On error
/// the set is unchanged.
pub fn image_onto(
    trajectory: &mut StructureSet,
    reference: &(impl AsStructureView + ?Sized),
) -> Result<(), PeriodicError> {
    let imager = MoleculeImager::new(trajectory.shared_topology());
    let reference = reference.as_structure_view();
    let positions = trajectory
        .iter()
        .enumerate()
        .map(|(index, item)| imager.reference_positions(index, item.as_structure_view(), reference))
        .collect::<Result<Vec<_>, _>>()?;
    publish(trajectory, positions)
}

/// Unwraps a trajectory through time using fractional-coordinate continuity.
///
/// The geometries are the frames, in set order, and the times of a time order
/// are their times; a set of another order or none has no times. Frame zero
/// is unchanged. For each subsequent frame, integer periodic images
/// are chosen so each atom moves less than half a cell in each periodic
/// fractional coordinate relative to the previous unwrapped frame. For changing
/// cells, previous fractional coordinates are retained and the current cell maps
/// the new fractional coordinates back to Cartesian space. This is the lattice
/// convention in Kulke & Vermaas (2022), equation B6, DOI: 10.1021/acs.jctc.2c00327.
///
/// Available times must not decrease. Requires fixed periodic-axis flags and
/// sequential, sufficiently closely sampled frames. Multiple crossings between
/// saved frames cannot be inferred; exact half-cell displacements are rejected
/// as ambiguous. Initial split molecules are not repaired: call
/// [`make_molecules_whole`] first when necessary. Apply before frame-dependent
/// alignment. This convention is not a universal prescription for diffusion
/// analysis under a fluctuating simulation cell.
///
/// # Errors
///
/// Returns [`PeriodicError::MissingCell`] for a frame without a cell,
/// [`PeriodicError::PeriodicAxesChanged`] when the periodic axes change,
/// [`PeriodicError::NonMonotonicTime`] when an available time decreases,
/// [`PeriodicError::AmbiguousDisplacement`] for an exact half-cell move,
/// [`PeriodicError::NumericalFailure`] for out-of-range coordinates, and
/// [`PeriodicError::Position`] or [`PeriodicError::Publication`] when the
/// result cannot be stored. On error the set is unchanged.
pub fn unwrap_trajectory(trajectory: &mut StructureSet) -> Result<(), PeriodicError> {
    let mut unwrapper = TrajectoryUnwrapper::new(trajectory.shared_topology());
    let positions = trajectory
        .iter()
        .enumerate()
        .map(|(index, geometry)| {
            let time = match geometry.place() {
                Some(Place::Time { time, .. }) => time,
                _ => None,
            };
            let frame =
                FrameView::bound(geometry.topology_handle(), geometry.payload(), time, None);
            unwrapper.next_positions(index, frame)
        })
        .collect::<Result<Vec<_>, _>>()?;
    publish(trajectory, positions)
}

fn whole_positions(
    trajectory: &StructureSet,
    anchors: Option<&AtomSelection>,
) -> Result<Vec<Positions>, PeriodicError> {
    let imager = MoleculeImager::new(trajectory.shared_topology());
    let anchors = anchors
        .map(|selection| imager.anchor_groups(selection))
        .transpose()?;
    trajectory
        .iter()
        .enumerate()
        .map(|(index, item)| {
            imager.frame_positions(index, item.as_structure_view(), anchors.as_deref())
        })
        .collect()
}

fn publish(trajectory: &mut StructureSet, positions: Vec<Positions>) -> Result<(), PeriodicError> {
    trajectory
        .replace_positions(positions)
        .map_err(|error| PeriodicError::Publication(Box::new(error)))
}

pub(super) fn positions(points: Vec<Point3>, frame: usize) -> Result<Positions, PeriodicError> {
    Positions::new(Quantity::new(points, CANONICAL_LENGTH_UNIT)).map_err(|source| {
        PeriodicError::Position {
            frame,
            source: Box::new(source),
        }
    })
}

pub(super) struct MoleculePlan {
    pub(super) groups: Vec<Vec<usize>>,
    pub(super) atom_group: Vec<usize>,
    tree: Vec<(usize, usize, usize)>,
    bonds: Vec<(InstanceBondId, usize, usize)>,
}

impl MoleculePlan {
    pub(super) fn new(topology: &Topology) -> Self {
        let mut adjacency = vec![Vec::new(); topology.atom_count()];
        let mut bonds = Vec::with_capacity(topology.bond_count());
        for bond in topology.bonds() {
            let id = bond.id();
            let [a, b] = bond.atoms().map(|atom| atom.index().index());
            adjacency[a].push((b, bonds.len()));
            adjacency[b].push((a, bonds.len()));
            bonds.push((id, a, b));
        }
        let mut groups = Vec::with_capacity(topology.instance_count());
        let mut atom_group = vec![0; topology.atom_count()];
        let mut tree = Vec::new();
        let mut visited = vec![false; topology.atom_count()];
        for molecule in topology.instances() {
            let atoms: Vec<_> = molecule.atoms().map(|atom| atom.index().index()).collect();
            for &atom in &atoms {
                atom_group[atom] = groups.len();
            }
            let root = atoms[0]; // Published molecules are nonempty and connected.
            visited[root] = true;
            let mut queue = VecDeque::from([root]);
            while let Some(parent) = queue.pop_front() {
                for &(child, bond) in &adjacency[parent] {
                    if !visited[child] {
                        visited[child] = true;
                        tree.push((parent, child, bond));
                        queue.push_back(child);
                    }
                }
            }
            groups.push(atoms);
        }
        Self {
            groups,
            atom_group,
            tree,
            bonds,
        }
    }

    pub(super) fn make_whole(
        &self,
        source: &[Point3],
        lattice: &Lattice,
    ) -> Result<Vec<Point3>, PeriodicError> {
        let mut points = source.to_vec();
        self.make_whole_in_place(&mut points, lattice)?;
        Ok(points)
    }

    /// Makes every molecule of `points` whole in place. Every image is
    /// chosen and checked before any point moves, so failure changes
    /// nothing.
    pub(super) fn make_whole_in_place(
        &self,
        points: &mut [Point3],
        lattice: &Lattice,
    ) -> Result<(), PeriodicError> {
        let mut atom_images = vec![[0_i64; 3]; points.len()];
        // Choose each bond image once, in its authoritative endpoint order. Tree
        // traversal may reverse a bond; reusing its image also makes ties consistent.
        let bond_images = self
            .bonds
            .iter()
            .map(|&(_, a, b)| {
                lattice
                    .nearest_image(points[b] - points[a])
                    .map(|(_, images)| images)
            })
            .collect::<Result<Vec<_>, _>>()?;
        for &(parent, child, bond) in &self.tree {
            let direction = if parent == self.bonds[bond].1 { 1 } else { -1 };
            let images = bond_images[bond];
            for axis in 0..3 {
                let image = atom_images[parent][axis] + direction * images[axis];
                checked_image(image as f64, lattice.frame)?;
                atom_images[child][axis] = image;
            }
        }
        for (&(bond, a, b), images) in self.bonds.iter().zip(bond_images) {
            // Check closure in integer lattice space, without a length tolerance
            // that could conceal a winding cycle in a short cell direction.
            if (0..3).any(|axis| atom_images[b][axis] - atom_images[a][axis] != images[axis]) {
                return Err(PeriodicError::InconsistentBondImages {
                    frame: lattice.frame,
                    bond,
                });
            }
        }
        for &(_, child, _) in &self.tree {
            points[child] = points[child] - lattice.cartesian(atom_images[child].map(|n| n as f64));
        }
        Ok(())
    }

    /// Adds each atom of `points`, its molecule moved to the periodic image
    /// nearest the same molecule in `reference` as [`Self::image_onto`]
    /// moves it, to the atom's entry of `sum`. An imaged atom that is not
    /// finite fails as the frame's position, and `sum` is then partly added.
    pub(super) fn add_imaged_onto(
        &self,
        points: &[Point3],
        reference: &[Point3],
        lattice: &Lattice,
        sum: &mut [Vector3],
    ) -> Result<(), PeriodicError> {
        for group in &self.groups {
            let delta =
                center(points, group, lattice.frame)? - center(reference, group, lattice.frame)?;
            let (_, images) = lattice.nearest_image(delta)?;
            let shift = lattice.cartesian(images.map(|n| n as f64));
            for &atom in group {
                let imaged = points[atom] - shift;
                if !imaged.is_finite() {
                    return Err(PeriodicError::Position {
                        frame: lattice.frame,
                        source: Box::new(DenseArrayError::NonFinite { index: atom }),
                    });
                }
                sum[atom] += imaged - Point3::origin();
            }
        }
        Ok(())
    }

    pub(super) fn image_onto(
        &self,
        points: &mut [Point3],
        reference: &[Point3],
        lattice: &Lattice,
    ) -> Result<(), PeriodicError> {
        for group in &self.groups {
            let delta =
                center(points, group, lattice.frame)? - center(reference, group, lattice.frame)?;
            let (_, images) = lattice.nearest_image(delta)?;
            let shift = lattice.cartesian(images.map(|n| n as f64));
            for &atom in group {
                points[atom] = points[atom] - shift;
            }
        }
        Ok(())
    }

    pub(super) fn image(
        &self,
        points: &mut [Point3],
        lattice: &Lattice,
        anchors: &[bool],
    ) -> Result<(), PeriodicError> {
        let first = anchors
            .iter()
            .position(|selected| *selected)
            .expect("nonempty checked anchors");
        let first_center = center(points, &self.groups[first], lattice.frame)?;
        for (group, &anchor) in self.groups.iter().zip(anchors) {
            if anchor {
                let delta = center(points, group, lattice.frame)? - first_center;
                let shift = lattice.minimum_image(delta)? - delta;
                for &atom in group {
                    points[atom] = points[atom] + shift;
                }
            }
        }
        let anchor_atoms: Vec<_> = self
            .groups
            .iter()
            .zip(anchors)
            .filter(|(_, anchor)| **anchor)
            .flat_map(|(atoms, _)| atoms.iter().copied())
            .collect();
        let anchor_center = center(points, &anchor_atoms, lattice.frame)?;
        let fractional = lattice.fractional(anchor_center - Point3::origin())?;
        let shift = lattice.cartesian(std::array::from_fn(|axis| {
            if lattice.periodic[axis] {
                0.5 - fractional[axis]
            } else {
                0.0
            }
        }));
        for point in points.iter_mut() {
            *point = *point + shift;
        }
        let anchor_center = anchor_center + shift;
        for (group, &anchor) in self.groups.iter().zip(anchors) {
            if !anchor {
                let delta = center(points, group, lattice.frame)? - anchor_center;
                let shift = lattice.minimum_image(delta)? - delta;
                for &atom in group {
                    points[atom] = points[atom] + shift;
                }
            }
        }
        Ok(())
    }
}

fn center(points: &[Point3], atoms: &[usize], frame: usize) -> Result<Point3, PeriodicError> {
    let scale = 1.0 / atoms.len() as f64;
    let mut mean = Vector3::zero();
    for &atom in atoms {
        mean += (points[atom] - Point3::origin()) * scale;
    }
    if !mean.is_finite() {
        return Err(PeriodicError::NumericalFailure { frame });
    }
    Ok(Point3::origin() + mean)
}

pub(super) struct Lattice {
    geometry: crate::math::Lattice,
    pub(super) periodic: [bool; 3],
    frame: usize,
}
impl Lattice {
    pub(super) fn new(cell: Option<PeriodicCell>, frame: usize) -> Result<Self, PeriodicError> {
        let cell = cell.ok_or(PeriodicError::MissingCell { frame })?;
        Ok(Self {
            geometry: crate::math::Lattice::new(cell)
                .map_err(|error| geometry_error(error, frame))?,
            periodic: cell.periodic_axes(),
            frame,
        })
    }
    pub(super) fn fractional(&self, vector: Vector3) -> Result<[f64; 3], PeriodicError> {
        self.geometry
            .fractional(Quantity::new(vector, CANONICAL_LENGTH_UNIT))
            .map_err(|error| geometry_error(error, self.frame))
    }
    pub(super) fn cartesian(&self, fractional: [f64; 3]) -> Vector3 {
        self.geometry.cartesian(fractional).into_value()
    }
    fn minimum_image(&self, delta: Vector3) -> Result<Vector3, PeriodicError> {
        self.geometry
            .minimum_image(Quantity::new(delta, CANONICAL_LENGTH_UNIT))
            .map(Quantity::into_value)
            .map_err(|error| geometry_error(error, self.frame))
    }
    fn nearest_image(&self, delta: Vector3) -> Result<(Vector3, [i64; 3]), PeriodicError> {
        self.geometry
            .nearest_image(Quantity::new(delta, CANONICAL_LENGTH_UNIT))
            .map(|(vector, images)| (vector.into_value(), images))
            .map_err(|error| geometry_error(error, self.frame))
    }
}
fn geometry_error(error: crate::math::LatticeError, frame: usize) -> PeriodicError {
    match error {
        crate::math::LatticeError::NumericalFailure => PeriodicError::NumericalFailure { frame },
        crate::math::LatticeError::LimitExceeded(limit) => {
            PeriodicError::LimitExceeded { frame, limit }
        }
        _ => PeriodicError::NumericalFailure { frame },
    }
}
pub(super) fn checked_image(value: f64, frame: usize) -> Result<f64, PeriodicError> {
    if !value.is_finite() || value.abs() >= 2.0_f64.powi(52) {
        return Err(PeriodicError::NumericalFailure { frame });
    }
    Ok(value.round())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::NANOMETER;

    #[test]
    fn nearest_images_match_exhaustive_cartesian_search_for_all_axis_flags() {
        let base = [
            Vector3::new(1.0, 0.2, -0.1),
            Vector3::new(0.4, 1.2, 0.3),
            Vector3::new(0.2, -0.1, 0.8),
        ];
        for basis in [
            base,
            base.map(|v| Vector3::new(-v.y, v.x, v.z)),
            [base[1], base[0], base[2]],
        ] {
            for mask in 1..8 {
                let axes = std::array::from_fn(|axis| mask & (1 << axis) != 0);
                let cell = PeriodicCell::new(Quantity::new(basis, NANOMETER), axes).unwrap();
                let lattice = Lattice::new(Some(cell), 7).unwrap();
                for geometry in 0..19 {
                    let delta = Vector3::new(
                        f64::from((geometry * 17) % 31 - 15) / 11.0,
                        f64::from((geometry * 11) % 29 - 14) / 9.0,
                        f64::from((geometry * 7) % 23 - 11) / 8.0,
                    );
                    let (nearest, images) = lattice.nearest_image(delta).unwrap();
                    let mut exhaustive = f64::INFINITY;
                    for a in -4..=4 {
                        for b in -4..=4 {
                            for c in -4..=4 {
                                if (!axes[0] && a != 0)
                                    || (!axes[1] && b != 0)
                                    || (!axes[2] && c != 0)
                                {
                                    continue;
                                }
                                let residual = delta
                                    - basis[0] * f64::from(a)
                                    - basis[1] * f64::from(b)
                                    - basis[2] * f64::from(c);
                                exhaustive = exhaustive.min(residual.norm_squared());
                            }
                        }
                    }
                    assert!(
                        (nearest.norm_squared() - exhaustive).abs() < 1.0e-12,
                        "basis={basis:?}, axes={axes:?}, delta={delta:?}, images={images:?}"
                    );
                    for axis in 0..3 {
                        assert!(axes[axis] || images[axis] == 0);
                    }
                    assert!(
                        (nearest - (delta - lattice.cartesian(images.map(|n| n as f64)))).norm()
                            < 1.0e-12
                    );
                }
            }
        }
    }

    #[test]
    fn rectangular_cells_use_independent_images_even_with_extreme_aspect_ratios() {
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(1.0, 1.0e-6, 1.0e6), NANOMETER),
            [true; 3],
        )
        .unwrap();
        let lattice = Lattice::new(Some(cell), 0).unwrap();
        let (nearest, images) = lattice
            .nearest_image(Vector3::new(1.25, 1.25e-6, 1.25e6))
            .unwrap();
        assert_eq!(images, [1; 3]);
        assert!((nearest.x - 0.25).abs() < 1.0e-14);
        assert!((nearest.y - 0.25e-6).abs() < 1.0e-20);
        assert!((nearest.z - 0.25e6).abs() < 1.0e-8);
    }

    #[test]
    fn nearest_image_fails_explicitly_on_excessive_search_or_lost_integer_precision() {
        let cell = PeriodicCell::new(
            Quantity::new(
                [
                    Vector3::new(1.0, 0.0, 0.0),
                    Vector3::new(0.99999, 1.0e-6, 0.0),
                    Vector3::new(0.0, 0.0, 1.0),
                ],
                NANOMETER,
            ),
            [true; 3],
        )
        .unwrap();
        let lattice = Lattice::new(Some(cell), 7).unwrap();
        assert!(matches!(
            lattice.nearest_image(Vector3::new(0.49, 0.49, 0.0)),
            Err(PeriodicError::LimitExceeded { frame: 7, limit })
                if limit.resource() == crate::diagnostics::Resource::LatticeImages
        ));
        assert!(matches!(
            lattice.nearest_image(Vector3::new(1.0e100, 0.0, 0.0)),
            Err(PeriodicError::NumericalFailure { frame: 7 })
        ));
    }
}
