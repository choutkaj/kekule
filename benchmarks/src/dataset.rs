use crate::*;

pub(crate) const DATASETS: &[&str] = &[
    "pubchem-100k",
    "enamine-diversity",
    "pl-rex",
    "pdb-1000",
    "smoke",
];

/// Preloaded source bytes. Timing never rereads the input file.
pub(crate) struct Input {
    pub(crate) path: PathBuf,
    pub(crate) text: String,
}

impl std::ops::Deref for Input {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.path
    }
}

impl Input {
    pub(crate) fn read(path: &Path) -> Result<Self, Box<dyn Error>> {
        Ok(Self {
            path: path.to_owned(),
            text: fs::read_to_string(path)?,
        })
    }
}

#[derive(Deserialize)]
pub(crate) struct Reference {
    pub(crate) tool: String,
    pub(crate) version: String,
    pub(crate) fixtures: Vec<String>,
    pub(crate) notes: Vec<String>,
}

pub(crate) struct Dataset {
    pub(crate) root: PathBuf,
    pub(crate) lock: Value,
    pub(crate) references: BTreeMap<String, Reference>,
    pub(crate) lock_sha256: String,
}

pub(crate) fn sha256(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

pub(crate) fn safe_join(root: &Path, relative: &str) -> Result<PathBuf, Box<dyn Error>> {
    if relative.is_empty()
        || Path::new(relative)
            .components()
            .any(|c| !matches!(c, std::path::Component::Normal(_)))
    {
        return Err(boxed_error(format!("unsafe relative path: {relative}")));
    }
    Ok(root.join(relative))
}

impl Dataset {
    pub(crate) fn open(id: &str) -> Result<Self, Box<dyn Error>> {
        if !DATASETS.contains(&id) {
            return Err(boxed_error(format!("unknown dataset: {id}")));
        }
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("corpora")
            .join(id);
        let bytes = fs::read(root.join("sources.lock.json"))?;
        let lock: Value = serde_json::from_slice(&bytes)?;
        if lock["corpus_id"] != id {
            return Err(boxed_error("source lock dataset identity mismatch"));
        }
        Ok(Self {
            references: serde_json::from_slice(&fs::read(root.join("references.json"))?)?,
            lock,
            root,
            lock_sha256: sha256(&bytes),
        })
    }

    /// The PDB source order preserves its existing 10/100/1000 nesting.
    /// Molecules use one hash order independent of feature and input format.
    pub(crate) fn selection(&self, limit: usize) -> Result<BTreeSet<String>, Box<dyn Error>> {
        let entries = self.lock["entries"]
            .as_array()
            .ok_or("missing source entries")?;
        let mut ids = entries
            .iter()
            .map(|entry| {
                entry["id"]
                    .as_str()
                    .map(str::to_owned)
                    .ok_or("missing source ID")
            })
            .collect::<Result<Vec<_>, _>>()?;
        if self.lock["corpus_id"] != "pdb-1000" {
            ids.sort_by_cached_key(|id| {
                (
                    sha256(format!("kekule-benchmark-v1:{id}").as_bytes()),
                    id.clone(),
                )
            });
        }
        if ids.iter().collect::<BTreeSet<_>>().len() != ids.len() {
            return Err(boxed_error("duplicate source IDs"));
        }
        Ok(ids.into_iter().take(limit).collect())
    }

    pub(crate) fn members(&self, fixture: &str) -> Result<Vec<String>, Box<dyn Error>> {
        if let Some(packs) = self.lock["packs"].as_array() {
            if let Some(pack) = packs.iter().find(|pack| pack["path"] == fixture) {
                return pack["members"]
                    .as_array()
                    .ok_or("missing pack members")?
                    .iter()
                    .map(|id| {
                        id.as_str()
                            .map(str::to_owned)
                            .ok_or_else(|| boxed_error("invalid pack member"))
                    })
                    .collect();
            }
        }
        let entries = self.lock["entries"]
            .as_array()
            .ok_or("missing source entries")?;
        let entry = entries
            .iter()
            .find(|entry| {
                entry["files"]
                    .as_array()
                    .is_some_and(|files| files.iter().any(|f| f["path"] == fixture))
            })
            .ok_or_else(|| boxed_error(format!("fixture absent from source lock: {fixture}")))?;
        Ok(vec![entry["id"]
            .as_str()
            .ok_or("missing source ID")?
            .to_owned()])
    }

    pub(crate) fn reference(
        &self,
        feature: &str,
        fixture: &str,
        input: &Input,
        reference: &Reference,
    ) -> Result<Value, Box<dyn Error>> {
        let path = self
            .root
            .join("golden")
            .join(feature)
            .join(format!("{}.json.gz", slugify_fixture(fixture)));
        let mut bytes = Vec::new();
        GzDecoder::new(fs::File::open(&path)?).read_to_end(&mut bytes)?;
        let golden: Value = serde_json::from_slice(&bytes)?;
        let pinned = self.lock["packs"]
            .as_array()
            .into_iter()
            .flatten()
            .chain(
                self.lock["entries"]
                    .as_array()
                    .into_iter()
                    .flatten()
                    .flat_map(|entry| entry["files"].as_array().into_iter().flatten()),
            )
            .find(|file| file["path"] == fixture)
            .and_then(|file| file["sha256"].as_str())
            .ok_or("fixture checksum missing from source lock")?;
        if pinned != sha256(input.text.as_bytes()) {
            return Err(boxed_error(format!(
                "input differs from source lock: {fixture}"
            )));
        }
        let schema = if matches!(
            feature,
            "stereo.representation"
                | "stereo.perception"
                | "io.smiles.isomeric"
                | "io.smiles.canonical"
        ) {
            2
        } else {
            1
        };
        let version = golden["reference"]["version"]
            .as_str()
            .ok_or("missing reference version")?;
        let version_matches = reference.version == version
            || reference.version == format!("RDKit {version}")
            || reference.version == format!("Biopython {version}");
        if golden["schema_version"] != schema
            || golden["feature_id"] != feature
            || golden["corpus_id"] != self.lock["corpus_id"]
            || golden["fixture_path"] != fixture
            || golden["input_sha256"] != sha256(input.text.as_bytes())
            || golden["reference"]["tool"] != reference.tool
            || !version_matches
            || golden["reference"]["runtime_dependency"] != false
        {
            return Err(boxed_error(format!(
                "reference provenance mismatch: {}",
                path.display()
            )));
        }
        Ok(golden
            .get("expected")
            .ok_or("missing reference output")?
            .clone())
    }
}

pub(crate) fn split_records(input: &Input, count: usize) -> Result<Vec<String>, Box<dyn Error>> {
    let records = if count == 1 {
        vec![input.text.clone()]
    } else if input.extension().is_some_and(|ext| ext == "sdf") {
        input
            .text
            .split_inclusive("$$$$")
            .filter(|s| !s.trim().is_empty())
            .enumerate()
            .map(|(i, s)| {
                if i == 0 {
                    s.to_owned()
                } else {
                    s.strip_prefix("\r\n")
                        .or_else(|| s.strip_prefix('\n'))
                        .unwrap_or(s)
                        .to_owned()
                }
            })
            .collect()
    } else {
        input
            .text
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(str::to_owned)
            .collect()
    };
    if records.len() != count {
        return Err(boxed_error(format!(
            "source membership count {count} differs from input record count {}",
            records.len()
        )));
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn changed_input_is_rejected_even_when_reference_metadata_is_consistent() {
        let root = env::temp_dir().join(format!(
            "kekule-reference-integrity-{}-{}",
            process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let directory = root.join("golden/io.smiles.parse");
        fs::create_dir_all(&directory).unwrap();
        let path = directory.join("data_input.smi.json.gz");
        let bytes = b"CC\n";
        let golden = json!({"schema_version":1,"feature_id":"io.smiles.parse","corpus_id":"test","fixture_path":"data/input.smi","input_sha256":sha256(bytes),"reference":{"tool":"rdkit","version":"test","runtime_dependency":false},"expected":{"records":[]}});
        let mut encoder = flate2::write::GzEncoder::new(
            fs::File::create(&path).unwrap(),
            flate2::Compression::default(),
        );
        encoder.write_all(golden.to_string().as_bytes()).unwrap();
        encoder.finish().unwrap();
        let dataset = Dataset {
            root: root.clone(),
            references: BTreeMap::new(),
            lock_sha256: String::new(),
            lock: json!({"corpus_id":"test","entries":[{"id":"1","files":[{"path":"data/input.smi","sha256":sha256(bytes)}]}]}),
        };
        let reference = Reference {
            tool: "rdkit".into(),
            version: "RDKit test".into(),
            fixtures: vec!["data/input.smi".into()],
            notes: vec![],
        };
        let mut input = Input {
            path: "data/input.smi".into(),
            text: "CC\n".into(),
        };
        assert!(dataset
            .reference("io.smiles.parse", "data/input.smi", &input, &reference)
            .is_ok());
        input.text = "CCC\n".into();
        assert!(dataset
            .reference("io.smiles.parse", "data/input.smi", &input, &reference)
            .is_err());
        fs::remove_file(path).unwrap();
        fs::remove_dir(directory).unwrap();
        fs::remove_dir(root.join("golden")).unwrap();
        fs::remove_dir(root).unwrap();
    }
    #[test]
    fn molecular_selections_are_nested_and_independent_of_file_order() {
        let mut dataset = Dataset {
            root: PathBuf::new(),
            references: BTreeMap::new(),
            lock_sha256: String::new(),
            lock: json!({"corpus_id":"pubchem-100k","entries":[{"id":"1"},{"id":"2"},{"id":"3"},{"id":"4"}]}),
        };
        let small = dataset.selection(2).unwrap();
        assert!(small.is_subset(&dataset.selection(3).unwrap()));
        dataset.lock["entries"].as_array_mut().unwrap().reverse();
        assert_eq!(small, dataset.selection(2).unwrap());
        dataset.lock["entries"][1]["id"] = dataset.lock["entries"][0]["id"].clone();
        assert!(dataset.selection(2).is_err());
    }
    #[test]
    fn pdb_selection_retains_locked_prefix() {
        let dataset = Dataset {
            root: PathBuf::new(),
            references: BTreeMap::new(),
            lock_sha256: String::new(),
            lock: json!({"corpus_id":"pdb-1000","entries":[{"id":"Z"},{"id":"A"},{"id":"B"}]}),
        };
        assert_eq!(
            dataset.selection(2).unwrap(),
            BTreeSet::from(["Z".into(), "A".into()])
        );
    }
    #[test]
    fn rejects_escaping_paths() {
        for path in ["../secret", "/absolute", "C:\\secret", "a/../../b", ""] {
            assert!(safe_join(Path::new("data"), path).is_err());
        }
    }
    #[test]
    fn sdf_split_preserves_blank_titles_and_rejects_missing_members() {
        let input = Input {
            path: "pack.sdf".into(),
            text: "\nprogram\n\nfirst\n$$$$\n\nprogram\n\nsecond\n$$$$\n".into(),
        };
        let records = split_records(&input, 2).unwrap();
        assert!(records[0].starts_with("\nprogram"));
        assert!(records[1].starts_with("\nprogram"));
        assert!(split_records(&input, 3).is_err());
    }
}
