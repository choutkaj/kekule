//! Streaming superposition and trajectory reductions.
//!
//! In-memory superposition and RMSD live on every realization collection in
//! [`kekule::alignment`] (`trajectory.superpose(0, &selection)`,
//! `trajectory.rmsd(&reference, &correspondence)`, ...). This module adds the
//! streaming counterpart, [`FrameSuperposer`], and per-atom reductions that
//! run either over a loaded collection ([`rmsf`], [`contact_occupancy`]) or
//! frame by frame through accumulators.

use kekule::alignment::{
    kabsch_with_options, AlignmentError, AlignmentOptions, FitAtoms, RigidAlignment,
};
use kekule::structure::{AsModelView, ModelView};

use crate::FrameBuffer;

mod reductions;
pub use reductions::*;

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
/// let atoms = AtomSelection::all(&reference.shared_topology());
/// let superposer = FrameSuperposer::new(&reference, &atoms);
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
    /// [`kekule::alignment::AtomCorrespondence`] (independent reference).
    pub fn new(reference: &'a (impl AsModelView + ?Sized), atoms: impl Into<FitAtoms<'a>>) -> Self {
        Self::with_options(reference, atoms, AlignmentOptions::default())
    }

    /// [`Self::new`] with explicit fit weights and periodic policy.
    pub fn with_options(
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
        let alignment = kabsch_with_options(
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
