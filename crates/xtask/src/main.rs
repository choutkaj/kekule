use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::error::Error;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process;
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::GzDecoder;
use molecules::{
    bio::{self, MacroMolecule, MmcifParseOptions, Residue},
    core::{Atom, AtomId, AtomRadical, Bond, BondOrder, BondStereo, Molecule, Point3, PropValue},
    molfile,
    perception::{self, AromaticityModel, SanitizeError, SanitizeOptions, ValenceModel},
    sdf::{self, SdfParseOptions, SdfRecord},
    small::SmallMolecule,
    smiles::{self, CanonicalSmilesWriteOptions, SmilesParseOptions, SmilesWriteOptions},
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

const VALIDATION_CORPORA: &[(&str, &str)] = &[
    ("tiny", "Tiny"),
    ("pubchem-100", "PubChem 100"),
    ("pubchem-1000", "PubChem 1000"),
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
