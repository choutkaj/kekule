use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::error::Error;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process;

use flate2::read::GzDecoder;
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
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

mod compare;
mod dataset;
mod features;
mod runner;
use compare::*;
use dataset::*;
use features::*;

fn boxed_error(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::other(message.into()).into()
}

fn main() {
    if let Err(error) = runner::run() {
        eprintln!("error: {error}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests;
