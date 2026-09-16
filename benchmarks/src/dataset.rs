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

pub(crate) struct Dataset {
    pub(crate) root: PathBuf,
    pub(crate) lock: Value,
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
        Ok(Self { lock, root })
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

    /// Source membership and format alone select inputs. Every available input
    /// in that format is tested, even if an ID has multiple source files.
    pub(crate) fn fixtures(&self, feature: &str) -> Vec<String> {
        let accepts = |path: &str| {
            let ext = Path::new(path)
                .extension()
                .and_then(|s| s.to_str())
                .unwrap_or("");
            if feature == "query.smarts" {
                matches!(ext, "smarts" | "sma" | "smi" | "smiles" | "txt")
            } else if feature.starts_with("io.smiles.") {
                matches!(ext, "smi" | "smiles" | "txt")
            } else if feature.starts_with("io.mmcif.") || feature.starts_with("bio.") {
                matches!(ext, "cif" | "mmcif")
            } else if feature.starts_with("io.sdf.") || feature.starts_with("io.mol.") {
                matches!(ext, "sdf" | "mol" | "mdl")
            } else {
                matches!(ext, "sdf" | "mol" | "mdl" | "smi" | "smiles" | "txt")
            }
        };
        let mut paths = self.lock["packs"]
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
            .filter_map(|file| file["path"].as_str())
            .filter(|path| accepts(path))
            .map(str::to_owned)
            .collect::<Vec<_>>();
        paths.sort();
        paths.dedup();
        paths
    }

    pub(crate) fn verify(&self, fixture: &str, input: &Input) -> Result<(), Box<dyn Error>> {
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
            .ok_or("input checksum missing from source lock")?;
        if pinned != sha256(input.text.as_bytes()) {
            return Err(boxed_error(format!(
                "input differs from source lock: {fixture}"
            )));
        }
        Ok(())
    }
}

pub(crate) fn split_records(input: &Input, count: usize) -> Result<Vec<String>, Box<dyn Error>> {
    let records = match input.extension().and_then(|e| e.to_str()) {
        Some("mol" | "mdl" | "cif" | "mmcif") => vec![input.text.clone()],
        Some("sdf") => {
            let mut records = Vec::new();
            let mut record = String::new();
            for line in input.text.split_inclusive('\n') {
                record.push_str(line);
                if line.trim_end_matches(['\r', '\n']) == "$$$$" {
                    records.push(std::mem::take(&mut record));
                }
            }
            if !record.trim().is_empty() {
                records.push(record);
            }
            records
        }
        _ => input
            .text
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(str::to_owned)
            .collect(),
    };
    if records.is_empty() || (count != 1 && records.len() != count) {
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
    fn selection_does_not_exclude_inputs_based_on_their_names() {
        let dataset = Dataset {
            root: PathBuf::new(),
            lock: json!({"entries":[{"id":"external", "files":[{"path":"data/smarts-example.smi"}]}]}),
        };
        assert_eq!(
            dataset.fixtures("io.smiles.parse"),
            vec!["data/smarts-example.smi"]
        );
        assert_eq!(
            dataset.fixtures("algo.rings.fast"),
            vec!["data/smarts-example.smi"]
        );
    }
    #[test]
    fn molecular_selections_are_nested_and_independent_of_file_order() {
        let mut dataset = Dataset {
            root: PathBuf::new(),
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
