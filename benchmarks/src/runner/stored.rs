use super::*;

/// Publish a complete file without replacing an existing set of goldens.
pub(super) struct GeneratedGoldens {
    path: PathBuf,
    temporary: PathBuf,
    writer: Option<flate2::write::GzEncoder<fs::File>>,
}
impl GeneratedGoldens {
    pub(super) fn create(path: &Path) -> Result<Self, Box<dyn Error>> {
        if path.exists() {
            return Err(boxed_error(format!(
                "stored goldens already exist: {}; choose a new --goldens directory to regenerate",
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
                flate2::Compression::default(),
            )),
        })
    }
    pub(super) fn finish(mut self) -> Result<(), Box<dyn Error>> {
        let file = self
            .writer
            .take()
            .ok_or("golden writer already closed")?
            .finish()?;
        file.sync_all()?;
        drop(file);
        // A hard link publishes atomically and fails if the destination exists.
        fs::hard_link(&self.temporary, &self.path)?;
        Ok(())
    }
}
impl Write for GeneratedGoldens {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.writer.as_mut().expect("open golden writer").write(buf)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.as_mut().expect("open golden writer").flush()
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

type Key = (String, Option<String>, Option<usize>);
pub(super) struct StoredGoldens(BTreeMap<Key, Golden>);

impl StoredGoldens {
    pub(super) fn read(
        reader: impl BufRead,
        dataset: &str,
        feature: &str,
        selected: &BTreeSet<String>,
    ) -> Result<Self, Box<dyn Error>> {
        let mut entries = BTreeMap::new();
        let mut seen = BTreeSet::new();
        let mut versions = BTreeSet::new();
        for line in reader.lines() {
            let golden: Golden = serde_json::from_str(&line?)?;
            if golden.dataset != dataset || golden.feature != feature {
                return Err(boxed_error("golden dataset/feature mismatch"));
            }
            if let Some(reference) = &golden.reference {
                let tool = if matches!(feature, "io.mmcif.parse" | "bio.secondary-structure.dssp") {
                    "biopython"
                } else {
                    "rdkit"
                };
                if reference.tool != tool || reference.version.is_empty() {
                    return Err(boxed_error("golden has no valid independent reference"));
                }
                versions.insert((reference.tool.clone(), reference.version.clone()));
            }
            if matches!(golden.expected, Outcome::Ok { .. })
                && (golden.reference.is_none()
                    || golden.fixture.is_none()
                    || golden.record_index.is_none()
                    || golden.input_sha256.as_ref().is_none_or(|hash| {
                        hash.len() != 64 || !hash.bytes().all(|b| b.is_ascii_hexdigit())
                    }))
            {
                return Err(boxed_error(
                    "successful golden is missing reference/input identity",
                ));
            }
            let key = (
                golden.id.clone(),
                golden.fixture.clone(),
                golden.record_index,
            );
            if !seen.insert(key.clone()) {
                return Err(boxed_error("duplicate golden case"));
            }
            if selected.contains(&golden.id) {
                entries.insert(key, golden);
            }
        }
        if versions.len() > 1 {
            return Err(boxed_error("mixed reference versions in one golden file"));
        }
        Ok(Self(entries))
    }

    pub(super) fn load(
        path: &Path,
        dataset: &str,
        feature: &str,
        selected: &BTreeSet<String>,
    ) -> Result<Self, Box<dyn Error>> {
        let input = fs::File::open(path).map_err(|error| {
            let mut message = format!("cannot read stored goldens {}: {error}", path.display());
            if error.kind() == std::io::ErrorKind::NotFound {
                message.push_str(&format!(
                    "\nGenerate them once with cargo benchmark generate --feature {feature} --dataset {dataset} (use --python PATH to select the reference environment and --goldens DIR for a custom storage directory)."
                ));
            }
            boxed_error(message)
        })?;
        Self::read(
            BufReader::new(flate2::read::GzDecoder::new(input)),
            dataset,
            feature,
            selected,
        )
        .map_err(|error| {
            boxed_error(format!(
                "cannot load stored goldens {}: {error}",
                path.display()
            ))
        })
    }

    pub(super) fn expected(&self, case: &Case) -> Outcome {
        let key = (
            case.id.clone(),
            Some(case.fixture.clone()),
            Some(case.index),
        );
        match self.0.get(&key) {
            None => Outcome::Error {
                message: "missing stored golden for this case; generate it explicitly".into(),
            },
            Some(golden)
                if golden.input_sha256.as_deref()
                    != Some(sha256(case.input.text.as_bytes()).as_str()) =>
            {
                Outcome::Error {
                    message:
                        "stored golden input checksum differs; it cannot be used for this input"
                            .into(),
                }
            }
            Some(golden) => golden.expected.clone(),
        }
    }
}
