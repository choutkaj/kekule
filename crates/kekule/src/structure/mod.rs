//! Geometry-bearing structures of coordinate-free topologies.
//!
//! - Dense per-atom arrays: [`Positions`] (canonical length unit, nm),
//!   [`Velocities`] (nm/ps) and [`Forces`] (kJ/(mol·nm)). They carry no atom
//!   identity and are validated against a topology by their owner.
//! - [`Model`]: one shared, immutable [`Topology`](crate::topology::Topology)
//!   plus one complete structure of it: positions, an optional periodic cell,
//!   optional occupancies and B-factors, and structure-scoped properties.
//! - Publication: [`ModelBuilder`] assembles a new topology and its positions
//!   together; [`ModelEditor`] edits a detached draft of a model whose
//!   positions follow atoms through structural edits, and
//!   [`ModelEditor::append_model`] imports another structure
//!   ([`ModelAppend`], [`ModelAppendReport`]). [`ModelEditor::topology_mut`]
//!   ([`ModelEditorTopologyMut`]) and [`ModelEditor::properties_mut`]
//!   ([`ModelEditorPropertiesMut`]) reach annotations without geometry.
//! - Collections of one shared topology: [`Ensemble`] is a distribution of
//!   weighted samples ([`EnsembleSample`]) grouped into states
//!   ([`EnsembleState`]) by an [`EnsemblePartition`]; [`Trajectory`] is a time-ordered
//!   in-memory sequence of [`TrajectoryFrame`]s, which may carry time, step,
//!   velocities and forces. File codecs and streaming live in `kekule-traj`.
//!   [`Trajectory::to_ensemble`] samples chosen frames into an ensemble of
//!   one state ([`TrajectorySampling`]); [`TrajectorySampler`] does so for
//!   frames that arrive one at a time. [`Series`] is an ordered sequence of
//!   [`SeriesItem`]s that is neither temporal nor weighted, such as the
//!   steps of a geometry optimization or the points of a relaxed scan, each
//!   at a step or a coordinate value of its [`SeriesAxis`].
//! - Borrowed access: [`ModelView`] is the common input of coordinate-dependent
//!   kernels. Models, ensemble samples ([`EnsembleSampleView`]), trajectory
//!   frames ([`TrajectoryFrameView`]), and series items ([`SeriesItemView`])
//!   all borrow as one without copying ([`AsModelView`]); [`ModelAtomView`] pairs one topology atom with its
//!   structure state. [`ModelView::kind`] says which structure a view borrows
//!   ([`StructureKind`]), so sample and frame state is never silently dropped.
//! - [`measure`]: checked distances, angles, dihedrals, bond dihedrals and
//!   distance-based selections.
//! - [`placement`]: builder geometry on a [`ModelEditor`] draft: hydrogens
//!   completed at ideal positions or collapsed back into counts, an atom's
//!   ideal shape and free directions, atoms attached to existing ones or
//!   placed as new molecules, and names for the atoms placement adds.
//! - Geometry edits: [`DistanceEdit`], [`AngleEdit`] and [`DihedralEdit`] (see
//!   below), and rigid motions with [`Model::apply_rigid_transform`] and
//!   [`Model::transform_atoms`].
//!
//! Owners accept semantic atom and bond IDs (and, in `*_at` methods, dense
//! indexes), translate IDs to dense indexes, and check every array's length.
//! Each owner has its own error type ([`ModelError`], [`ModelBuildError`],
//! [`ModelEditError`], [`EnsembleError`], [`TrajectoryError`],
//! [`SeriesError`], [`GeometryEditError`]); [`StructureError`] and [`DenseArrayError`] report
//! invalid structure state shared by all of them.
//!
//! ```
//! use kekule::{geometry::Point3, smiles, units::{Quantity, NANOMETER}};
//! use kekule::structure::{Model, Positions};
//!
//! let positions = Positions::new(Quantity::new(
//!     [Point3::origin(), Point3::new(0.125, 0.0, 0.0)], NANOMETER,
//! ))?;
//! let model = Model::new(smiles::to_topology("CO")?, positions)?;
//! assert_eq!(model.atom_count(), 2);
//! let oxygen = model.topology().atom_ids()[1];
//! let position = model.position(oxygen).unwrap().into_unit(NANOMETER)?;
//! assert_eq!(position.into_value(), Point3::new(0.125, 0.0, 0.0));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Geometry editing and scans
//!
//! [`Model::set_distance`], [`Model::set_angle`], and [`Model::set_dihedral`]
//! move connected fragments. Prepare a [`DihedralEdit`] (or [`DistanceEdit`] /
//! [`AngleEdit`]) once for repeated absolute targets on the same topology:
//! ```
//! use kekule::{smiles, geometry::Point3, structure::{Model, Positions, DihedralEdit},
//!     units::{Quantity, ANGSTROM, DEGREE}};
//! let topology = smiles::to_topology("CCCC")?;
//! let positions = Positions::new(Quantity::new([
//!     Point3::new(0.0, 1.0, 0.0), Point3::origin(),
//!     Point3::new(1.5, 0.0, 0.0), Point3::new(1.5, 1.0, 1.0),
//! ], ANGSTROM))?;
//! let mut model = Model::new(topology, positions)?;
//! let [a, b, c, d] = model.topology().atom_ids().try_into().unwrap();
//! model.set_distance(b, c, Quantity::new(1.6, ANGSTROM))?;
//! let edit = DihedralEdit::new(model.shared_topology(), [a, b, c, d])?;
//! for degrees in [-180.0, -60.0, 60.0, 180.0] {
//!     let mut sample = model.clone();
//!     edit.apply(&mut sample, Quantity::new(degrees, DEGREE))?;
//!     // Evaluate or collect sample; model remains the scan baseline.
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Explicit selections can intentionally deform rings or boundary bonds:
//! ```
//! use kekule::{structure::{AngleEdit, Model}, topology::AtomSelection,
//!     units::{Quantity, DEGREE}};
//! # fn example(model: &mut Model) -> Result<(), Box<dyn std::error::Error>> {
//! let [a, b, c] = model.topology().atom_ids()[..3].try_into().unwrap();
//! let topology = model.shared_topology();
//! let moving = AtomSelection::from_atoms(topology.clone(), [c])?;
//! let edit = AngleEdit::with_moving_atoms(topology.clone(), [a, b, c], &moving)?;
//! edit.apply(model, Quantity::new(120.0, DEGREE))?;
//! // Only C moves; bonds to other unselected neighbors can deform.
//! # Ok(())
//! # }
//! ```
//!
//! Geometry edits are atomic, use stored Cartesian coordinates, and preserve
//! topology, properties, and cell. Apply periodic preprocessing explicitly.
//! [`Model::translate`], [`Model::rotate`], and [`Model::transform_atoms`] also
//! support rigid selection manipulation and keep the cell;
//! [`Model::apply_rigid_transform`] moves the whole model and rotates its cell.
//! No relaxation is performed.

mod dynamics;
mod ensemble;
mod geometry_edit;
pub mod measure;
mod model;
mod model_editor;
mod model_state;
mod partition;
pub mod placement;
mod positions;
mod sampling;
mod scalars;
mod series;
mod state;
mod trajectory;

pub use dynamics::{Forces, Velocities};
pub use ensemble::{
    Ensemble, EnsembleError, EnsembleIter, EnsembleSample, EnsembleSampleMut, EnsembleSampleView,
    EnsembleState, EnsembleStateIter, EnsembleStateView,
};
pub use geometry_edit::{AngleEdit, DihedralEdit, DistanceEdit, GeometryEditError};
pub use model::{
    AsModelView, Model, ModelAtomView, ModelBuildError, ModelBuilder, ModelError, ModelTopologyMut,
    ModelView, StructureKind,
};
pub use model_editor::{
    ModelAppend, ModelAppendReport, ModelEditError, ModelEditor, ModelEditorProperties,
    ModelEditorPropertiesMut, ModelEditorTopologyMut,
};
pub use partition::{EnsemblePartition, PartitionError};
pub use positions::{DenseArrayError, Positions};
pub use sampling::{TrajectorySampler, TrajectorySampling};
pub use series::{
    Series, SeriesAxis, SeriesError, SeriesItem, SeriesItemMut, SeriesItemView, SeriesIter,
};
pub use state::StructureError;
pub use trajectory::{
    MissingTimePolicy, Trajectory, TrajectoryError, TrajectoryFrame, TrajectoryFrameMut,
    TrajectoryFrameView, TrajectoryIter,
};

#[cfg(test)]
mod model_tests;
#[cfg(test)]
mod tests;
