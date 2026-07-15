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
    perception::{
        self,
        aromaticity::{self, AromaticityModel},
        rings,
        stereo::{self, StereoCandidate, StereoPerceptionIssue, StereoPerceptionReport},
        valence::{self, ValenceModel, ValenceOptions},
        SanitizeError, SanitizeOptions,
    },
    query,
    sdf::{self, SdfDataField, SdfParseOptions, SdfRecord},
    small::SmallMolecule,
    smiles::{self, CanonicalSmilesWriteOptions, SmilesWriteOptions},
    substructure,
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ValidationCorpus {
    pub(crate) id: &'static str,
    pub(crate) label: &'static str,
    pub(crate) local_only: bool,
}

const VALIDATION_CORPORA: &[ValidationCorpus] = &[
    ValidationCorpus {
        id: "smoke",
        label: "Smoke",
        local_only: false,
    },
    ValidationCorpus {
        id: "pubchem-100",
        label: "PubChem 100",
        local_only: true,
    },
    ValidationCorpus {
        id: "pubchem-1k",
        label: "PubChem 1k",
        local_only: true,
    },
    ValidationCorpus {
        id: "pubchem-100k",
        label: "PubChem 100k",
        local_only: true,
    },
    ValidationCorpus {
        id: "pl-rex",
        label: "PL-REX",
        local_only: true,
    },
    ValidationCorpus {
        id: "enamine-diversity",
        label: "Enamine diversity",
        local_only: true,
    },
    ValidationCorpus {
        id: "pdb-10",
        label: "PDB 10",
        local_only: true,
    },
    ValidationCorpus {
        id: "pdb-100",
        label: "PDB 100",
        local_only: true,
    },
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
