//! Reading, interpreting, and writing SMILES, Molfile, SDF, mmCIF, PDB, XYZ,
//! and the native format.

#[path = "../support/mod.rs"]
mod support;

mod assemblies;
mod cancellation;
mod canonical_smiles;
mod components;
mod documents;
mod export;
mod mmcif;
mod native;
mod pdb;
mod pdb_write;
mod sample_sets;
mod sequences;
mod smiles_stereo;
mod writers;
mod xyz;

/// The FNV-1a 64 fingerprint of `bytes`, which pins written output too long
/// to state as text.
pub(crate) fn fingerprint(bytes: impl AsRef<[u8]>) -> u64 {
    bytes
        .as_ref()
        .iter()
        .fold(0xcbf2_9ce4_8422_2325, |hash: u64, &byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3)
        })
}

/// The [`fingerprint`] of the native encoding of `source`.
pub(crate) fn native_fingerprint<'a>(source: impl Into<kekule::native::NativeSource<'a>>) -> u64 {
    fingerprint(
        kekule::native::write(source, kekule::native::NativeWriteOptions::default()).unwrap(),
    )
}
