//! Opening, loading, saving and streaming trajectories through the public API.

#[path = "../support/mod.rs"]
mod support;

mod dense_vectors;
mod frame_buffers;
mod loading;
mod readers;
mod saving;
mod streaming;
mod streaming_transforms;
