//! Streaming superposition and trajectory reductions.
//!
//! In-memory superposition and RMSD are trajectory and ensemble methods in
//! [`kekule::alignment`] (`trajectory.superpose(0, &selection, options)`,
//! `trajectory.rmsd(&reference, &correspondence, options)`, ...). This module
//! adds the streaming counterpart, [`FrameSuperposer`], and per-atom trajectory
//! reductions, in which every frame counts once, that run either over a loaded
//! trajectory ([`rmsf`], [`contact_occupancy`]) or frame by frame through
//! accumulators.

use kekule::alignment::{self, AlignmentError, AlignmentOptions, FitAtoms, RigidAlignment};
use kekule::structure::{AsModelView, ModelView, TrajectoryFrameView};

use crate::FrameBuffer;

mod sealed {
    /// Restricts [`super::AsFrameView`] to trajectory frames and buffers.
    pub trait Sealed {}
}

/// A trajectory frame that reductions can observe: a borrowed
/// [`TrajectoryFrameView`] or a [`FrameBuffer`].
///
/// Reductions count every frame once, so they accept trajectory frames only.
/// Ensemble members are weighted samples whose statistics must use their
/// weights, and a [`ModelView`] may borrow either kind of structure, so
/// neither implements this sealed trait:
///
/// ```compile_fail,E0277
/// # use kekule::{structure::Ensemble, topology::AtomSelection};
/// # use kekule_traj::analysis::RmsfAccumulator;
/// fn unweighted(ensemble: &Ensemble, atoms: &AtomSelection) {
///     let mut rmsf = RmsfAccumulator::new(atoms).unwrap();
///     let _ = rmsf.observe(0, &ensemble.get(0).unwrap());
/// }
/// ```
///
/// ```compile_fail,E0277
/// # use kekule::{structure::Ensemble, topology::InstanceAtomId, units::{Quantity, ANGSTROM}};
/// # use kekule_traj::analysis::ContactOccupancyAccumulator;
/// fn unweighted(ensemble: &Ensemble, pairs: Vec<(InstanceAtomId, InstanceAtomId)>) {
///     let topology = ensemble.shared_topology();
///     let cutoff = Quantity::new(4.0, ANGSTROM);
///     let mut contacts = ContactOccupancyAccumulator::new(&topology, pairs, cutoff).unwrap();
///     let _ = contacts.observe(0, &ensemble.get(0).unwrap().as_model_view());
/// }
/// ```
pub trait AsFrameView: sealed::Sealed {
    /// Borrows the frame.
    fn as_frame_view(&self) -> TrajectoryFrameView<'_>;
}

impl sealed::Sealed for TrajectoryFrameView<'_> {}
impl sealed::Sealed for FrameBuffer {}

impl AsFrameView for TrajectoryFrameView<'_> {
    fn as_frame_view(&self) -> TrajectoryFrameView<'_> {
        *self
    }
}

impl AsFrameView for FrameBuffer {
    fn as_frame_view(&self) -> TrajectoryFrameView<'_> {
        self.frame_view()
    }
}

mod reductions;
pub use reductions::{
    contact_occupancy, rmsf, ContactOccupancyAccumulator, ContactOccupancyResult, ReductionError,
    RmsfAccumulator, RmsfResult,
};

/// A borrowed reference and fitting pairs for superposing streamed frames.
///
/// Uses the same kernel as in-memory collection superposition. Each call fits
/// one buffered frame and transforms its positions, cell, velocities, and
/// forces together. To keep a reference while reusing its input buffer,
/// retain it first, for example with `buffer.frame_view().to_model()`.
///
/// ```no_run
/// use kekule::topology::AtomSelection;
/// use kekule_traj::{analysis::FrameSuperposer, TrajectoryReader};
/// # fn fit(reader: &mut impl TrajectoryReader) -> Result<(), Box<dyn std::error::Error>> {
/// let mut buffer = reader.frame_buffer();
/// reader.read_next(&mut buffer)?;
/// let reference = buffer.frame_view().to_model();
/// let atoms = AtomSelection::all(reference.shared_topology());
/// let superposer = FrameSuperposer::new(&reference, &atoms, Default::default());
/// let mut index = 1;
/// while reader.read_next(&mut buffer)? {
///     let fit = superposer.superpose(index, &mut buffer)?;
///     # let _ = fit;
///     index += 1;
/// }
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy)]
pub struct FrameSuperposer<'a> {
    reference: ModelView<'a>,
    atoms: FitAtoms<'a>,
    options: AlignmentOptions<'a>,
}

impl<'a> FrameSuperposer<'a> {
    /// Fits onto `reference` through a selection (shared layout) or an
    /// [`kekule::alignment::AtomCorrespondence`] (independent reference), with
    /// the fit weights and periodic policy of `options`.
    pub fn new(
        reference: &'a (impl AsModelView + ?Sized),
        atoms: impl Into<FitAtoms<'a>>,
        options: AlignmentOptions<'a>,
    ) -> Self {
        Self {
            reference: reference.as_model_view(),
            atoms: atoms.into(),
            options,
        }
    }

    /// Fits the buffered frame and transforms it in place, returning the
    /// applied transform and post-fit RMSD. The buffer changes only when the
    /// complete transformed frame is valid.
    ///
    /// `frame_index` is the caller's zero-based source index; failures are
    /// reported as [`AlignmentError::Item`] with that index, like the items of
    /// an in-memory collection.
    pub fn superpose(
        &self,
        frame_index: usize,
        buffer: &mut FrameBuffer,
    ) -> Result<RigidAlignment, AlignmentError> {
        self.fit(buffer).map_err(|source| AlignmentError::Item {
            index: frame_index,
            source: Box::new(source),
        })
    }

    fn fit(&self, buffer: &mut FrameBuffer) -> Result<RigidAlignment, AlignmentError> {
        let alignment = alignment::fit(
            buffer.as_model_view(),
            self.reference,
            self.atoms,
            self.options,
        )?;
        buffer
            .frame_mut()
            .apply_rigid_transform(&alignment.transform())
            .map_err(AlignmentError::Transform)?;
        Ok(alignment)
    }
}
