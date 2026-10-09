//! Data-only model import; no checkpoint deserialization or network access.
use super::{
    config::Config,
    network::{Network, Tensor},
    ModelIdentity,
};
use crate::{Error, ErrorKind, Result};
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::{collections::BTreeMap, io::Read, path::Path};

pub(super) const PREPARATION: &str = "openff-nagl-0.6.1";
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Entry {
    pub inchi: String,
    pub mapped_smiles: String,
    pub charges: Vec<f32>,
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Domain {
    pub allowed_elements: Vec<u8>,
    pub forbidden_patterns: Vec<String>,
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {
    schema: u32,
    model: String,
    checkpoint_sha256: String,
    weights_sha256: String,
    tensors: BTreeMap<String, Tensor>,
    config: Config,
    domain: Domain,
    lookup_tables: BTreeMap<String, Vec<Entry>>,
    preparation: Option<String>,
}
pub(super) struct Bundle {
    pub identity: ModelIdentity,
    pub config: Config,
    pub network: Network,
    pub domain: Domain,
    pub lookup: BTreeMap<String, Entry>,
}

const METADATA_LIMIT: u64 = 64 * 1024 * 1024;
pub(super) const WEIGHTS_LIMIT: u64 = 256 * 1024 * 1024;

fn error(detail: impl std::fmt::Display) -> Error {
    Error::new(ErrorKind::Model, detail)
}

fn read_bounded(path: &Path, limit: u64) -> Result<Vec<u8>> {
    let file = std::fs::File::open(path).map_err(|e| Error::io(path, e))?;
    let metadata = file.metadata().map_err(|e| Error::io(path, e))?;
    if !metadata.is_file() || metadata.len() > limit {
        return Err(
            error("NAGL bundle file is not a regular file within the size limit").at_path(path),
        );
    }
    let mut bytes = Vec::new();
    file.take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| Error::io(path, e))?;
    if bytes.len() as u64 > limit {
        return Err(error("NAGL bundle file size exceeds limit").at_path(path));
    }
    Ok(bytes)
}

/// Loads `model.json` and `weights.bin` from one bundle directory.
pub(super) fn load(directory: &Path) -> Result<Bundle> {
    let metadata = read_bounded(&directory.join("model.json"), METADATA_LIMIT)?;
    let weights = read_bounded(&directory.join("weights.bin"), WEIGHTS_LIMIT)?;
    parse(&metadata, &weights).map_err(|e| e.at_path(directory))
}

/// Validates one complete schema-2 bundle held in memory.
pub(super) fn parse(metadata: &[u8], weights: &[u8]) -> Result<Bundle> {
    if metadata.len() as u64 > METADATA_LIMIT || weights.len() as u64 > WEIGHTS_LIMIT {
        return Err(error("NAGL bundle size exceeds limit"));
    }
    let manifest: Manifest =
        serde_json::from_slice(metadata).map_err(|e| Error::wrap(ErrorKind::Model, e))?;
    if manifest.schema != 2 || manifest.preparation.as_deref() != Some(PREPARATION) {
        return Err(error(format!(
            "unsupported NAGL bundle: expected schema 2 with preparation profile {PREPARATION}"
        )));
    }
    let identity =
        ModelIdentity::new(manifest.model, manifest.checkpoint_sha256).ok_or_else(|| {
            error("NAGL bundle requires a nonempty model and a 64-digit checkpoint SHA-256")
        })?;
    manifest.config.validate()?;
    if manifest
        .domain
        .allowed_elements
        .iter()
        .any(|&z| !(1..=118).contains(&z))
        || manifest.domain.allowed_elements.len() > 118
        || manifest.domain.forbidden_patterns.len() > 1024
    {
        return Err(error("invalid NAGL chemical domain"));
    }
    if format!("{:x}", Sha256::digest(weights)) != manifest.weights_sha256 {
        return Err(error("NAGL weights checksum mismatch"));
    }
    let network = Network::load(&manifest.config, manifest.tensors, weights)?;
    let readout = manifest.config.readouts.keys().next().unwrap();
    let mut lookup = BTreeMap::new();
    for (name, entries) in manifest.lookup_tables {
        if &name != readout {
            return Err(error("NAGL lookup table has no matching charge readout"));
        }
        for entry in entries {
            // Upstream Ash includes one entry with an empty InChI. Retain it
            // verbatim; the identity adapter never returns an empty success.
            if entry.mapped_smiles.is_empty()
                || entry.charges.is_empty()
                || entry.charges.iter().any(|q| !q.is_finite())
            {
                return Err(error("invalid NAGL lookup entry"));
            }
            if lookup.insert(entry.inchi.clone(), entry).is_some() {
                return Err(error("duplicate NAGL lookup identifier"));
            }
        }
    }
    Ok(Bundle {
        identity,
        config: manifest.config,
        network,
        domain: manifest.domain,
        lookup,
    })
}
