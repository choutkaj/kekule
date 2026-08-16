use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::error::Error;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::{read::GzDecoder, Compression, GzBuilder};
use kekule::{
    canon,
    core::{
        Atom, AtomId, AtomRadical, AxisOrientation, Bond, BondId, BondOrder, DoubleBondOrientation,
        Molecule, StereoBondMark, StereoBondMarkKind, StereoCarrier, StereoDescriptor,
        StereoElement, StereoElementKind, StereoGroup, StereoGroupKind, StereoSource,
        StereoSpecifiedness, TetrahedralOrientation,
    },
    dssp, hydrogens,
    mmcif::{self, MmcifInterpretOptions, MmcifModelSelection, MmcifParseOptions},
    molfile,
    normalization::{NormalizationReport, NormalizationWarning},
    perception::{
        rings,
        valence::{self, ValenceModel, ValenceOptions},
    },
    query,
    sdf::{self, SdfDataField, SdfParseOptions, SdfRecord},
    small::SmallMolecule,
    smiles::{self, CanonicalSmilesWriteOptions, SmilesWriteOptions},
    stereo::{
        self, CoordinateStereoError, CoordinateStereoMaterializationReport, StereoCandidate,
        StereoValidationIssue,
    },
    substructure,
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BenchmarkCorpus {
    pub(crate) id: &'static str,
    pub(crate) label: &'static str,
    pub(crate) local_only: bool,
    pub(crate) default: bool,
}

const BENCHMARK_CORPORA: &[BenchmarkCorpus] = &[
    BenchmarkCorpus {
        id: "pubchem-1k",
        label: "PubChem 1k",
        local_only: true,
        default: true,
    },
    BenchmarkCorpus {
        id: "pubchem-100k",
        label: "PubChem 100k",
        local_only: true,
        default: false,
    },
    BenchmarkCorpus {
        id: "pl-rex",
        label: "PL-REX",
        local_only: true,
        default: false,
    },
    BenchmarkCorpus {
        id: "enamine-diversity",
        label: "Enamine diversity",
        local_only: true,
        default: false,
    },
    BenchmarkCorpus {
        id: "pdb-100",
        label: "PDB 100",
        local_only: true,
        default: true,
    },
    BenchmarkCorpus {
        id: "pdb-1000",
        label: "PDB 1000",
        local_only: true,
        default: false,
    },
];
const DASHBOARD_PATH: &str = "features/DASHBOARD.html";
const BENCHMARK_INPUT_DIGEST_SCHEMA_VERSION: u32 = 1;
const GOLDEN_SCHEMA_VERSION: u32 = 1;
const COMPARISON_MODE_IMPLEMENTATION_GOLDEN: &str = "implementation-golden";

mod benchmark;
mod cli;
mod corpus;
mod dashboard;
mod features;
mod skills;
mod support;

pub(crate) use benchmark::*;
pub(crate) use cli::*;
pub(crate) use corpus::*;
pub(crate) use dashboard::*;
pub(crate) use features::*;
pub(crate) use skills::*;
pub(crate) use support::*;

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests;
