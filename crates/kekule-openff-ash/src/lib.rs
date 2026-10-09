//! The OpenFF Ash NAGL charge model, `openff-gnn-am1bcc-1.0.0`, as data.
//!
//! This crate only carries the converted model; load it with
//! `kekule_openff::NaglModel::ash()`. The weights are unchanged from the
//! checksum-pinned upstream checkpoint (see `ATTRIBUTION.md`) and are licensed
//! under CC BY 4.0, which requires retaining that attribution when
//! redistributing them. The Rust code is MIT OR Apache-2.0.
#![forbid(unsafe_code)]
#![no_std]

/// Upstream checkpoint file name, as declared by the Rosemary force field.
pub const MODEL_FILE: &str = "openff-gnn-am1bcc-1.0.0.pt";

/// SHA-256 of the upstream PyTorch checkpoint the bundle was exported from.
pub const CHECKPOINT_SHA256: &str =
    "7981e7f5b0b1e424c9e10a40d9e7606d96dcd3dd2b095cb4eeff6829f92238ee";

/// The schema-2 bundle manifest (`model.json`): configuration, tensor layout,
/// chemical domain, lookup table, and the SHA-256 of the decoded weights.
pub const MANIFEST: &[u8] = include_bytes!("../data/model.json");

/// The little-endian float32 weights, stored losslessly as four byte planes
/// (every value's first byte, then every second byte, and so on) compressed
/// with zlib. Decoding must reproduce the manifest's `weights_sha256`.
pub const WEIGHT_PLANES_ZLIB: &[u8] = include_bytes!("../data/weights.planes.zlib");
