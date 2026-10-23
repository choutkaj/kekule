//! Opening, loading, saving and streaming trajectories through the public API.

#[path = "../support/mod.rs"]
mod support;

mod dense_vectors;
mod frame_buffers;
mod loading;
mod readers;
mod reference_imaging;
mod saving;
mod streaming;
mod streaming_transforms;
