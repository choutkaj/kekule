use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::error::Error;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process;

use flate2::{read::GzDecoder, Compression, GzBuilder};
use kekule::{
    canon,
    core::{
        Atom, AtomId, AtomRadical, AxisOrientation, Bond, BondId, BondOrder, DoubleBondOrientation,
        Molecule, StereoCarrier, StereoDescriptor, StereoElement, StereoElementId,
        StereoElementKind, StereoGroup, StereoGroupKind, TetrahedralOrientation,
    },
    dssp, hydrogens,
    mmcif::{self, MmcifInterpretOptions, MmcifModelSelection},
    molfile,
    perception::{
        rings,
        valence::{self, ValenceModel, ValenceOptions},
    },
    query,
    sdf::{self, SdfDataField, SdfRecordInterpretation},
    smiles::{self},
    stereo::{
        self, CoordinateStereoError, CoordinateStereoMaterializationReport, StereoCandidate,
        StereoValidationIssue,
    },
    substructure,
};
use serde::Deserialize;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

const GOLDEN_SCHEMA_VERSION: u32 = 1;
const COMPARISON_MODE_IMPLEMENTATION_GOLDEN: &str = "implementation-golden";

mod benchmark;
mod cli;
mod corpus;
mod support;

pub(crate) use benchmark::*;
pub(crate) use cli::*;
pub(crate) use corpus::*;
pub(crate) use support::*;

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests;
