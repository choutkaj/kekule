use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::error::Error;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::{read::GzDecoder, Compression, GzBuilder};
use molecules::{
    core::{
        Atom, AtomId, AtomRadical, AxisOrientation, Bond, BondId, BondOrder, DoubleBondOrientation,
        Molecule, StereoBondMark, StereoBondMarkKind, StereoCarrier, StereoDescriptor,
        StereoElement, StereoElementKind, StereoGroup, StereoGroupKind, StereoSource,
        StereoSpecifiedness, TetrahedralOrientation,
    },
    molfile,
    perception::{
        self,
        aromaticity::{self, AromaticityModel},
        rings,
        stereo::{self, StereoCandidate, StereoPerceptionIssue, StereoPerceptionReport},
        valence::{self, ValenceModel, ValenceOptions},
        SanitizeError, SanitizeOptions,
    },
    sdf::{self, SdfDataField, SdfParseOptions, SdfRecord},
    small::SmallMolecule,
    smiles::{self, CanonicalSmilesWriteOptions, SmilesWriteOptions},
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

const VALIDATION_CORPORA: &[(&str, &str)] = &[
    ("smoke", "Smoke"),
    ("pubchem-100", "PubChem 100"),
    ("pubchem-1k", "PubChem 1k"),
    ("pubchem-100k", "PubChem 100k"),
    ("pl-rex", "PL-REX"),
    ("enamine-diversity", "Enamine diversity"),
    ("pdb-10", "PDB 10"),
    ("pdb-100", "PDB 100"),
];
const DASHBOARD_PATH: &str = "features/DASHBOARD.html";
const VALIDATION_EVIDENCE_SCHEMA_VERSION: u32 = 2;
const GOLDEN_SCHEMA_VERSION: u32 = 1;
const COMPARISON_MODE_IMPLEMENTATION_GOLDEN: &str = "implementation-golden";

mod cli;
mod corpus;
mod dashboard;
mod features;
mod skills;
mod support;
mod validation;

pub(crate) use cli::*;
pub(crate) use corpus::*;
pub(crate) use dashboard::*;
pub(crate) use features::*;
pub(crate) use skills::*;
pub(crate) use support::*;
pub(crate) use validation::*;

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests;
