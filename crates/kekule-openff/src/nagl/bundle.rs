//! Data-only model import; no checkpoint deserialization or network access.
use super::{
    config::Config,
    network::{Network, Tensor},
    ModelIdentity,
};
use crate::{error, Result};
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

fn read_bounded(path: &Path, limit: u64) -> Result<Vec<u8>> {
    let file = std::fs::File::open(path).map_err(error)?;
    let metadata = file.metadata().map_err(error)?;
    if !metadata.is_file() || metadata.len() > limit {
        return Err(error(format!(
            "NAGL bundle file size exceeds limit: {}",
            path.display()
        )));
    }
    let mut bytes = Vec::new();
    file.take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(error)?;
    if bytes.len() as u64 > limit {
        return Err(error("NAGL bundle file size exceeds limit"));
    }
    Ok(bytes)
}
pub(super) fn load(directory: &Path) -> Result<Bundle> {
    let metadata = read_bounded(&directory.join("model.json"), 64 * 1024 * 1024)?;
    let manifest: Manifest = serde_json::from_slice(&metadata).map_err(error)?;
    match manifest.schema {
        // Preserve the original immutable Ash export as a supported preset.
        1 if format!("{:x}", Sha256::digest(&metadata)) == "f359ed50ada12a120195464ef0259837884e73826829e562652cb5e46ddb9425" => {},
        2 if manifest.preparation.as_deref() == Some(PREPARATION) => {},
        _ => return Err(error("unsupported NAGL bundle schema or preparation profile (legacy schema 1 requires the pinned Ash manifest)")),
    }
    let identity = ModelIdentity::new(manifest.model, manifest.checkpoint_sha256)?;
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
    let raw = read_bounded(&directory.join("weights.bin"), 256 * 1024 * 1024)?;
    if format!("{:x}", Sha256::digest(&raw)) != manifest.weights_sha256 {
        return Err(error("NAGL weights checksum mismatch"));
    }
    let network = Network::load(&manifest.config, manifest.tensors, &raw)?;
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
                || entry.charges.len() > 4096
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
