//! Periodic coordinate reconstruction, imaging, and temporal unwrapping.
//!
//! These operations change positions only. The exact shared topology, cells,
//! velocities, forces, time, step, and all properties are retained. Copy-returning
//! methods leave their source unchanged; `_in_place` methods publish only after
//! every frame succeeds. Every processed frame must have a periodic cell.
//!
//! Molecules and bonds come from the authoritative topology. No bonds are guessed
//! from distances or hierarchy. Orthorhombic, triclinic, rotated, and partially
//! periodic cells are supported. Molecular reconstruction uses shortest Cartesian
//! bond images; temporal unwrapping uses continuity in fractional cell coordinates.
//! Molecular image ties are resolved deterministically. The exact nearest-image
//! search permits up to one million candidates per displacement and fails explicitly
//! for cells requiring more; it never substitutes an approximate image silently.
//!
//! ```no_run
//! use kekule::{mmcif, topology::AtomSelection};
//! use kekule_traj::io::read_trajectory;
//!
//! let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?)?;
//! let topology = document.interpret()?.to_topology();
//! let mut trajectory = read_trajectory("trajectory.xtc", topology.clone())?;
//! trajectory.make_molecules_whole_in_place()?;
//!
//! // A continuous path through time, retaining the initial molecular images.
//! let continuous = trajectory.unwrap()?;
//!
//! // A separate per-frame view with molecules centered around selected anchors.
//! let anchors = AtomSelection::all(&topology); // Or select the solute's atoms.
//! let imaged = trajectory.image_molecules(&anchors)?;
//! let aligned = imaged.superpose_to_frame(0, &anchors)?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! For a large trajectory, retain the bond plan and temporal state across reads:
//!
//! ```no_run
//! use kekule::mmcif;
//! use kekule_traj::{
//!     io::open_trajectory, periodic::{MoleculeImager, TrajectoryUnwrapper},
//!     TrajectoryReader,
//! };
//! let topology = mmcif::parse_str(&std::fs::read_to_string("system.cif")?)?
//!     .interpret()?.to_topology();
//! let mut reader = open_trajectory("trajectory.xtc", topology.clone())?;
//! let mut frame = reader.frame_buffer();
//! let imager = MoleculeImager::new(topology.clone());
//! let mut unwrapper = TrajectoryUnwrapper::new(topology);
//! let mut index = 0;
//! while reader.read_next(&mut frame)? {
//!     imager.make_whole_in_place(index, &mut frame)?;
//!     unwrapper.unwrap_in_place(index, &mut frame)?;
//!     // Analyze or write this frame here. Downsample only after unwrapping.
//!     index += 1;
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::collections::VecDeque;
use std::fmt;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::structure::{PositionError, Positions};
use kekule::topology::{
    AtomSelection, InstanceAtomId, InstanceBondId, Topology, TopologyAtomIndex,
};
use kekule::units::{Quantity, CANONICAL_LENGTH_UNIT};

use crate::{Trajectory, TrajectoryError};

mod stream;
pub use stream::{MoleculeImager, TrajectoryUnwrapper};

/// Failure to reconstruct or publish periodic coordinates.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PeriodicError {
    MissingCell {
        frame: usize,
    },
    SelectionTopologyMismatch,
    TopologyMismatch {
        frame: usize,
    },
    /// Stateful unwrapping requires consecutive source-frame indices.
    NonSequentialFrame {
        previous: usize,
        frame: usize,
    },
    /// Available time values must not decrease, including across missing values.
    NonMonotonicTime {
        frame: usize,
    },
    EmptyAnchors,
    /// Selected shortest bond images do not form a consistent finite molecule (e.g. a winding ring).
    InconsistentBondImages {
        frame: usize,
        bond: InstanceBondId,
    },
    PeriodicAxesChanged {
        frame: usize,
    },
    /// A displacement lies on a half-cell boundary, with no unique temporal image.
    AmbiguousDisplacement {
        frame: usize,
        atom: TopologyAtomIndex,
        axis: usize,
    },
    /// An ill-conditioned cell requires more than the bounded nearest-image search permits.
    ImageSearchLimit {
        frame: usize,
    },
    NumericalFailure {
        frame: usize,
    },
    Position {
        frame: usize,
        source: Box<PositionError>,
    },
    Publication(Box<TrajectoryError>),
}

impl fmt::Display for PeriodicError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingCell { frame } => write!(f, "frame {frame} has no periodic cell"),
            Self::SelectionTopologyMismatch => f.write_str("anchor selection belongs to another topology"),
            Self::TopologyMismatch { frame } => write!(f, "frame {frame} belongs to another topology"),
            Self::NonSequentialFrame { previous, frame } => write!(f, "frame {frame} must immediately follow frame {previous} for temporal unwrapping"),
            Self::NonMonotonicTime { frame } => write!(f, "time decreases at frame {frame} during temporal unwrapping"),
            Self::EmptyAnchors => f.write_str("imaging requires at least one anchor atom"),
            Self::InconsistentBondImages { frame, bond } => write!(f, "frame {frame} has inconsistent periodic images around bond {bond}"),
            Self::PeriodicAxesChanged { frame } => write!(f, "periodic axes change at frame {frame}"),
            Self::AmbiguousDisplacement { frame, atom, axis } => write!(f, "frame {frame}, atom {atom}, cell axis {axis}: ambiguous half-cell displacement"),
            Self::ImageSearchLimit { frame } => write!(f, "frame {frame}: cell geometry exceeds the nearest-image search limit"),
            Self::NumericalFailure { frame } => write!(f, "frame {frame}: periodic coordinate calculation exceeds numerical precision or range"),
            Self::Position { frame, source } => write!(f, "frame {frame}: {source}"),
            Self::Publication(source) => write!(f, "cannot publish periodic coordinates: {source}"),
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

impl Trajectory {
    /// Returns a copy with each bonded molecule made whole independently in each frame.
    ///
    /// The first atom of each molecule stays at its original position; other atoms
    /// are placed using shortest Cartesian bond images. All bonds, including ring
    /// closures, must agree with the reconstruction. Molecules are neither centered
    /// nor joined to other molecules. This does not establish continuity across time.
    pub fn make_molecules_whole(&self) -> Result<Self, PeriodicError> {
        let positions = self.whole_positions(None)?;
        self.with_positions(positions)
            .map_err(|e| PeriodicError::Publication(Box::new(e)))
    }

    /// Makes molecules whole transactionally, retaining all non-position state in place.
    pub fn make_molecules_whole_in_place(&mut self) -> Result<(), PeriodicError> {
        let positions = self.whole_positions(None)?;
        self.replace_positions(positions)
            .map_err(|e| PeriodicError::Publication(Box::new(e)))
    }

    /// Returns a copy with whole molecules imaged around explicitly selected anchors.
    ///
    /// Every molecule containing a selected atom becomes an anchor. Other anchors
    /// are placed at their nearest centroid image to the first anchor in topology
    /// order. The combined geometric center of all atoms in the anchor molecules
    /// is moved to the cell center along periodic axes. Other whole molecules are placed at their
    /// nearest centroid image to that combined center. No anchor heuristics are used.
    ///
    /// This includes making molecules whole. It is a per-frame centering operation,
    /// not temporal unwrapping. Coordinates along nonperiodic fractional axes remain
    /// unchanged; a molecule may extend outside the primary cell.
    pub fn image_molecules(&self, anchors: &AtomSelection) -> Result<Self, PeriodicError> {
        let positions = self.whole_positions(Some(anchors))?;
        self.with_positions(positions)
            .map_err(|e| PeriodicError::Publication(Box::new(e)))
    }

    /// Images whole molecules transactionally around the selected anchors.
    pub fn image_molecules_in_place(
        &mut self,
        anchors: &AtomSelection,
    ) -> Result<(), PeriodicError> {
        let positions = self.whole_positions(Some(anchors))?;
        self.replace_positions(positions)
            .map_err(|e| PeriodicError::Publication(Box::new(e)))
    }

    /// Returns a temporally unwrapped copy using fractional-coordinate continuity.
    ///
    /// Frame zero is unchanged. For each subsequent frame, integer periodic images
    /// are chosen so each atom moves less than half a cell in each periodic
    /// fractional coordinate relative to the previous unwrapped frame. For changing
    /// cells, previous fractional coordinates are retained and the current cell maps
    /// the new fractional coordinates back to Cartesian space. This is the lattice
    /// convention in Kulke & Vermaas (2022), equation B6, DOI: 10.1021/acs.jctc.2c00327.
    ///
    /// Available times must not decrease. Requires fixed periodic-axis flags and sequential, sufficiently closely
    /// sampled frames. Multiple crossings between saved frames cannot be inferred;
    /// exact half-cell displacements are rejected as ambiguous. Initial split
    /// molecules are not repaired: use `make_molecules_whole` first when necessary.
    /// Apply before frame-dependent alignment. This convention is not a universal
    /// prescription for diffusion analysis under a fluctuating simulation cell.
    pub fn unwrap(&self) -> Result<Self, PeriodicError> {
        let positions = self.unwrapped_positions()?;
        self.with_positions(positions)
            .map_err(|e| PeriodicError::Publication(Box::new(e)))
    }

    /// Applies fractional-coordinate temporal unwrapping transactionally in place.
    pub fn unwrap_in_place(&mut self) -> Result<(), PeriodicError> {
        let positions = self.unwrapped_positions()?;
        self.replace_positions(positions)
            .map_err(|e| PeriodicError::Publication(Box::new(e)))
    }

    fn whole_positions(
        &self,
        anchors: Option<&AtomSelection>,
    ) -> Result<Vec<Positions>, PeriodicError> {
        let imager = MoleculeImager::new(self.shared_topology());
        let anchors = anchors
            .map(|selection| imager.anchor_groups(selection))
            .transpose()?;
        self.frames()
            .enumerate()
            .map(|(index, frame)| imager.frame_positions(index, frame, anchors.as_deref()))
            .collect()
    }

    fn unwrapped_positions(&self) -> Result<Vec<Positions>, PeriodicError> {
        let mut unwrapper = TrajectoryUnwrapper::new(self.shared_topology());
        self.frames()
            .enumerate()
            .map(|(index, frame)| unwrapper.next_positions(index, frame))
            .collect()
    }
}

fn positions(points: Vec<Point3>, frame: usize) -> Result<Positions, PeriodicError> {
    Positions::new(Quantity::new(points, CANONICAL_LENGTH_UNIT)).map_err(|source| {
        PeriodicError::Position {
            frame,
            source: Box::new(source),
        }
    })
}

struct MoleculePlan {
    groups: Vec<Vec<usize>>,
    atom_group: Vec<usize>,
    tree: Vec<(usize, usize, usize)>,
    bonds: Vec<(InstanceBondId, usize, usize)>,
}

impl MoleculePlan {
    fn new(topology: &Topology) -> Self {
        let mut adjacency = vec![Vec::new(); topology.atom_count()];
        let mut bonds = Vec::with_capacity(topology.bond_count());
        for (id, bond) in topology.bonds() {
            let (a, b) = bond.endpoints();
            let a = topology
                .atom_index(InstanceAtomId::new(id.molecule(), a))
                .expect("validated bond atom")
                .index();
            let b = topology
                .atom_index(InstanceAtomId::new(id.molecule(), b))
                .expect("validated bond atom")
                .index();
            adjacency[a].push((b, bonds.len()));
            adjacency[b].push((a, bonds.len()));
            bonds.push((id, a, b));
        }
        let mut groups = Vec::with_capacity(topology.instance_count());
        let mut atom_group = vec![0; topology.atom_count()];
        let mut tree = Vec::new();
        let mut visited = vec![false; topology.atom_count()];
        for molecule in topology.molecules() {
            let atoms: Vec<_> = molecule
                .atoms()
                .map(|(id, _)| topology.atom_index(id).expect("validated atom").index())
                .collect();
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

    fn make_whole(
        &self,
        source: &[Point3],
        lattice: &Lattice,
    ) -> Result<Vec<Point3>, PeriodicError> {
        let mut points = source.to_vec();
        let mut atom_images = vec![[0_i64; 3]; source.len()];
        // Choose each bond image once, in its authoritative endpoint order. Tree
        // traversal may reverse a bond; reusing its image also makes ties consistent.
        let bond_images = self
            .bonds
            .iter()
            .map(|&(_, a, b)| {
                lattice
                    .nearest_image(source[b] - source[a])
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
            points[child] = source[child] - lattice.cartesian(atom_images[child].map(|n| n as f64));
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
        Ok(points)
    }

    fn image(
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

struct Lattice {
    basis: [Vector3; 3],
    reciprocal: [Vector3; 3],
    periodic: [bool; 3],
    frame: usize,
}

impl Lattice {
    fn new(cell: Option<PeriodicCell>, frame: usize) -> Result<Self, PeriodicError> {
        let cell = cell.ok_or(PeriodicError::MissingCell { frame })?;
        let [a, b, c] = cell.vectors().to_value();
        let inverse_volume = 1.0 / a.dot(b.cross(c));
        let reciprocal = [
            b.cross(c) * inverse_volume,
            c.cross(a) * inverse_volume,
            a.cross(b) * inverse_volume,
        ];
        if !reciprocal.iter().all(|v| v.is_finite()) {
            return Err(PeriodicError::NumericalFailure { frame });
        }
        Ok(Self {
            basis: [a, b, c],
            reciprocal,
            periodic: cell.periodic_axes(),
            frame,
        })
    }

    fn fractional(&self, vector: Vector3) -> Result<[f64; 3], PeriodicError> {
        let fractional = self.reciprocal.map(|dual| dual.dot(vector));
        if !fractional.iter().all(|v| v.is_finite()) {
            return Err(PeriodicError::NumericalFailure { frame: self.frame });
        }
        Ok(fractional)
    }

    fn cartesian(&self, fractional: [f64; 3]) -> Vector3 {
        self.basis[0] * fractional[0]
            + self.basis[1] * fractional[1]
            + self.basis[2] * fractional[2]
    }

    fn minimum_image(&self, delta: Vector3) -> Result<Vector3, PeriodicError> {
        self.nearest_image(delta).map(|(vector, _)| vector)
    }

    fn nearest_image(&self, delta: Vector3) -> Result<(Vector3, [i64; 3]), PeriodicError> {
        let fractional = self.fractional(delta)?;
        let mut images = [0.0; 3];
        for axis in 0..3 {
            if self.periodic[axis] {
                images[axis] = checked_image(fractional[axis], self.frame)?;
            }
        }
        let mut best = delta - self.cartesian(images);
        let mut best_images = images.map(|n| n as i64);
        let mut best_squared = best.norm_squared();
        if !best_squared.is_finite() {
            return Err(PeriodicError::NumericalFailure { frame: self.frame });
        }
        // In an exactly orthogonal basis the three image choices are independent.
        // Besides avoiding enumeration, this handles very elongated rectangular
        // cells without an unnecessarily large Cartesian search radius.
        if self.basis[0].dot(self.basis[1]) == 0.0
            && self.basis[0].dot(self.basis[2]) == 0.0
            && self.basis[1].dot(self.basis[2]) == 0.0
        {
            return Ok((best, best_images));
        }
        // If r bounds the best Cartesian residual, reciprocal-vector norms bound
        // every fractional residual. Enumerating this finite box therefore includes
        // the true closest lattice image, even for skewed or partially periodic cells.
        let radius = best_squared.sqrt();
        let mut bounds = [(0_i64, 0_i64); 3];
        let mut candidates = 1_u64;
        for axis in 0..3 {
            if self.periodic[axis] {
                let span = radius * self.reciprocal[axis].norm();
                let padding = 64.0 * f64::EPSILON * (1.0 + fractional[axis].abs() + span);
                let low = (fractional[axis] - span - padding).ceil();
                let high = (fractional[axis] + span + padding).floor();
                if !low.is_finite()
                    || !high.is_finite()
                    || low.abs().max(high.abs()) >= 2.0_f64.powi(52)
                {
                    return Err(PeriodicError::NumericalFailure { frame: self.frame });
                }
                let (low, high) = (low as i64, high as i64);
                let count = (high - low + 1).max(0) as u64;
                candidates = candidates
                    .checked_mul(count)
                    .filter(|n| *n <= 1_000_000)
                    .ok_or(PeriodicError::ImageSearchLimit { frame: self.frame })?;
                bounds[axis] = (low, high);
            }
        }
        for a in bounds[0].0..=bounds[0].1 {
            for b in bounds[1].0..=bounds[1].1 {
                for c in bounds[2].0..=bounds[2].1 {
                    let residual = delta - self.cartesian([a as f64, b as f64, c as f64]);
                    let squared = residual.norm_squared();
                    if squared < best_squared {
                        best = residual;
                        best_squared = squared;
                        best_images = [a, b, c];
                    }
                }
            }
        }
        Ok((best, best_images))
    }
}

fn checked_image(value: f64, frame: usize) -> Result<f64, PeriodicError> {
    if !value.is_finite() || value.abs() >= 2.0_f64.powi(52) {
        return Err(PeriodicError::NumericalFailure { frame });
    }
    Ok(value.round())
}

#[cfg(test)]
mod tests {
    use super::*;
    use kekule::units::NANOMETER;

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
                for sample in 0..19 {
                    let delta = Vector3::new(
                        f64::from((sample * 17) % 31 - 15) / 11.0,
                        f64::from((sample * 11) % 29 - 14) / 9.0,
                        f64::from((sample * 7) % 23 - 11) / 8.0,
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
            Err(PeriodicError::ImageSearchLimit { frame: 7 })
        ));
        assert!(matches!(
            lattice.nearest_image(Vector3::new(1.0e100, 0.0, 0.0)),
            Err(PeriodicError::NumericalFailure { frame: 7 })
        ));
    }
}
