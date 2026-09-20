use super::{Case, Outcome, Reference};
use crate::{boxed_error, dataset::sha256};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::io::Read;
use std::{
    error::Error,
    fs,
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
    process,
    time::{SystemTime, UNIX_EPOCH},
};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Metadata {
    pub schema: u32,
    pub contract_sha256: String,
    pub dataset: String,
    pub feature: String,
    pub input_lock_sha256: String,
    pub sha256: String,
    pub reference: Option<Reference>,
    pub reference_code_sha256: Option<String>,
    pub origin: String,
    pub cases: usize,
}

pub(super) fn contract_hash() -> String {
    text_hash(include_str!("../../contract.json"))
}
pub(super) fn feature_contract_hash(feature: &str) -> String {
    if feature == "query.smarts" {
        text_hash(&format!(
            "{}\n{}",
            contract_hash(),
            crate::features::query::CONTRACT
        ))
    } else {
        contract_hash()
    }
}
pub(super) fn text_hash(text: &str) -> String {
    sha256(text.replace("\r\n", "\n").as_bytes())
}
pub(super) fn file_hash(path: &Path) -> Result<String, Box<dyn Error>> {
    let mut reader = fs::File::open(path)?;
    let mut hash = Sha256::new();
    let mut buffer = [0; 65536];
    loop {
        let read = reader.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        hash.update(&buffer[..read]);
    }
    Ok(hash
        .finalize()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect())
}
pub(super) fn reference_code_hash(root: &Path) -> Result<String, Box<dyn Error>> {
    let mut hash = Sha256::new();
    for path in [
        "reference/run.py",
        "reference/rdkit/run_feature.py",
        "reference/rdkit/molecule.py",
        "reference/rdkit/source_radicals.py",
        "reference/biopython/run_feature.py",
        "queries.smarts",
        "query-smarts.json",
    ] {
        hash.update(path.as_bytes());
        hash.update(
            fs::read_to_string(root.join(path))?
                .replace("\r\n", "\n")
                .as_bytes(),
        );
    }
    Ok(hash
        .finalize()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect())
}
fn metadata_path(path: &Path) -> PathBuf {
    path.with_extension("meta.json")
}

/// The manifest is the commit marker. Readers require both files and their
/// exact digest. Interrupted generation and replacement are never accepted.
pub(super) struct GeneratedGoldens {
    path: PathBuf,
    temporary: PathBuf,
    writer: Option<flate2::write::GzEncoder<fs::File>>,
}
impl GeneratedGoldens {
    pub(super) fn create(path: &Path) -> Result<Self, Box<dyn Error>> {
        if path.exists() || metadata_path(path).exists() {
            return Err(boxed_error(format!(
                "stored goldens already exist: {}; choose a new --goldens directory",
                path.display()
            )));
        }
        let temporary = path.with_extension(format!(
            "{}-{}.tmp",
            process::id(),
            SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
        ));
        let file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)?;
        Ok(Self {
            path: path.into(),
            temporary,
            writer: Some(flate2::write::GzEncoder::new(
                file,
                flate2::Compression::best(),
            )),
        })
    }
    pub(super) fn finish(mut self, mut metadata: Metadata) -> Result<(), Box<dyn Error>> {
        let file = self.writer.take().ok_or("golden writer closed")?.finish()?;
        file.sync_all()?;
        drop(file);
        metadata.sha256 = file_hash(&self.temporary)?;
        let manifest = metadata_path(&self.path);
        let temporary_manifest = metadata_path(&self.temporary);
        let publish = (|| -> Result<(), Box<dyn Error>> {
            let mut file = fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary_manifest)?;
            serde_json::to_writer_pretty(&mut file, &metadata)?;
            file.sync_all()?;
            drop(file);
            fs::hard_link(&self.temporary, &self.path)?;
            if let Err(error) = fs::hard_link(&temporary_manifest, &manifest) {
                fs::remove_file(&self.path)?;
                return Err(error.into());
            }
            Ok(())
        })();
        let _ = fs::remove_file(temporary_manifest);
        publish
    }
}
impl Write for GeneratedGoldens {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.writer.as_mut().expect("open writer").write(buf)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.as_mut().expect("open writer").flush()
    }
}
impl Drop for GeneratedGoldens {
    fn drop(&mut self) {
        drop(self.writer.take());
        let _ = fs::remove_file(&self.temporary);
    }
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Golden {
    dataset: String,
    feature: String,
    id: String,
    fixture: Option<String>,
    record_index: Option<usize>,
    input_sha256: Option<String>,
    reference: Option<Reference>,
    expected: Outcome,
}
type Key = (String, usize, String);
fn key(id: &str, fixture: Option<&str>, index: Option<usize>) -> Key {
    (
        fixture.unwrap_or("\u{10ffff}").to_owned(),
        index.unwrap_or(usize::MAX),
        id.to_owned(),
    )
}

/// One lookahead record: memory is independent of corpus size. File order is
/// fixture, record index, source ID; absent-format entries follow real inputs.
pub(super) struct StoredGoldens {
    reader: Box<dyn BufRead>,
    pub metadata: Metadata,
    pending: Option<Golden>,
    previous: Option<Key>,
    records: usize,
}
impl StoredGoldens {
    pub(super) fn load(
        path: &Path,
        dataset: &str,
        feature: &str,
        lock_hash: &str,
    ) -> Result<Self, Box<dyn Error>> {
        let read = || -> Result<Self, Box<dyn Error>> {
            let metadata: Metadata = serde_json::from_reader(fs::File::open(metadata_path(path))?)?;
            if metadata.schema != 2
                || metadata.contract_sha256 != feature_contract_hash(feature)
                || metadata.dataset != dataset
                || metadata.feature != feature
                || metadata.input_lock_sha256 != lock_hash
            {
                return Err(boxed_error(
                    "golden manifest has stale schema, contract or input identity",
                ));
            }
            if file_hash(path)? != metadata.sha256 {
                return Err(boxed_error("golden checksum differs from manifest"));
            }
            let mut stored = Self {
                reader: Box::new(BufReader::new(flate2::read::GzDecoder::new(
                    fs::File::open(path)?,
                ))),
                metadata,
                pending: None,
                previous: None,
                records: 0,
            };
            stored.advance()?;
            Ok(stored)
        };
        read().map_err(|error| boxed_error(format!("cannot load stored goldens {}: {error}\nGenerate them explicitly with cargo benchmark generate --feature {feature} --dataset {dataset}",path.display())))
    }
    fn advance(&mut self) -> Result<(), Box<dyn Error>> {
        let mut line = String::new();
        if self.reader.read_line(&mut line)? == 0 {
            self.pending = None;
            if self.records != self.metadata.cases {
                return Err(boxed_error("golden case count differs from manifest"));
            }
            return Ok(());
        }
        let entry: Golden = serde_json::from_str(&line)?;
        if entry.dataset != self.metadata.dataset
            || entry.feature != self.metadata.feature
            || (entry.reference.is_some() && entry.reference != self.metadata.reference)
        {
            return Err(boxed_error(
                "golden dataset, feature or reference identity differs",
            ));
        }
        if let Some(reference) = &entry.reference {
            let tool = if matches!(
                entry.feature.as_str(),
                "io.mmcif.parse" | "bio.secondary-structure.dssp"
            ) {
                "biopython"
            } else {
                "rdkit"
            };
            if reference.tool != tool || reference.version.is_empty() {
                return Err(boxed_error("invalid independent reference"));
            }
        }
        if let Outcome::Ok { value } = &entry.expected {
            if entry.reference.is_none()
                || entry.fixture.is_none()
                || entry.record_index.is_none()
                || entry.input_sha256.as_ref().is_none_or(|hash| {
                    hash.len() != 64 || !hash.bytes().all(|b| b.is_ascii_hexdigit())
                })
            {
                return Err(boxed_error(
                    "successful golden missing reference or input identity",
                ));
            }
            crate::observation::validate(&entry.feature, value)?;
        }
        let key = key(&entry.id, entry.fixture.as_deref(), entry.record_index);
        if self
            .previous
            .as_ref()
            .is_some_and(|previous| previous >= &key)
        {
            return Err(boxed_error("duplicate or unordered golden case"));
        }
        self.previous = Some(key);
        self.records += 1;
        self.pending = Some(entry);
        Ok(())
    }
    pub(super) fn expected(&mut self, case: &Case) -> Result<Outcome, Box<dyn Error>> {
        let wanted = key(&case.id, Some(&case.fixture), Some(case.index));
        while self.pending.as_ref().is_some_and(|entry| {
            key(&entry.id, entry.fixture.as_deref(), entry.record_index) < wanted
        }) {
            self.advance()?;
        }
        let Some(entry) = self
            .pending
            .as_ref()
            .filter(|entry| key(&entry.id, entry.fixture.as_deref(), entry.record_index) == wanted)
        else {
            return Ok(Outcome::Error {
                message: "missing stored golden for selected case".into(),
            });
        };
        if entry.input_sha256.as_deref() != Some(sha256(case.input.text.as_bytes()).as_str()) {
            return Ok(Outcome::Error {
                message: "stored golden input checksum differs".into(),
            });
        }
        let outcome = self.pending.take().unwrap().expected;
        self.advance()?;
        Ok(outcome)
    }
    pub(super) fn finish(&mut self) -> Result<(), Box<dyn Error>> {
        while self.pending.is_some() {
            self.advance()?;
        }
        Ok(())
    }
}
