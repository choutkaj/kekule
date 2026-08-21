//! Fixed-topology molecular trajectory storage, streaming, file I/O, and
//! trajectory-oriented workflows built on [`kekule`].
//!
//! `kekule-traj` owns ordered frame state and trajectory-specific operations.
//! It reuses Kekule's immutable [`kekule::topology::Topology`], structural
//! topology-bound positions, cells, atom and bond data, units, geometry, selections,
//! and the borrowed [`kekule::structure::ModelView`] contract rather than
//! defining a second molecular model.
//!
//! The crate provides owned [`TrajectoryFrame`] values, reusable
//! [`FrameBuffer`] storage, finite in-memory [`Trajectory`] collections, and
//! streaming reader/writer traits. Production file codecs and format-agnostic
//! path factories live under [`io`].
//!
//! In-memory trajectory superposition and RMSD workflows live under
//! [`analysis`]. Direct RMSD never fits coordinates implicitly; the explicitly
//! named aligned RMSD convenience performs fitting without materializing a
//! transformed trajectory.
//!
//! Coordinate-dependent kernels in `kekule` and external potential adapters
//! can consume [`TrajectoryFrameView::model_view`] or
//! [`FrameBuffer::model_view`] without copying coordinates.
#![forbid(unsafe_code)]
#![warn(rustdoc::broken_intra_doc_links)]
// Kekule consistently names owned conversions `to_*`, including consuming ones.
#![allow(clippy::wrong_self_convention)]

mod trajectory;

pub mod analysis;
pub mod io;

pub use trajectory::*;
