use std::collections::BTreeMap;
use std::env;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use std::process;

use molecules::prelude::{
    perceive_aromaticity, perceive_ring_membership, perceive_ring_set, read_mmcif_str,
    read_mol_v2000_str, read_smiles_str, sanitize_small_molecule, write_mol_v2000, write_sdf_v2000,
    write_smiles, AromaticityModel, Atom, AtomId, Bond, BondOrder, BondStereo, MacroMolecule,
    MmcifParseOptions, Molecule, PropValue, SanitizeOptions, SdfParseOptions, SdfRecord,
    SmallMolecule, SmilesParseOptions, SmilesWriteOptions,
};
use molecules::read_sdf_v2000_records;
use serde_json::{json, Value};

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let mut args = env::args().skip(1);
    match args.next().as_deref() {
        Some("dashboard") => dashboard(args.collect()),
        Some("validate") => validate(args.collect()),
        Some("features") => list_features(),
        Some("skills") => skills(args.collect()),
        _ => {
            print_help();
            Ok(())
        }
    }
}

fn print_help() {
    eprintln!(
        "usage:\n  cargo xtask dashboard [--check]\n  cargo xtask validate --feature FEATURE_ID\n  cargo xtask skills --check\n  cargo xtask features"
    );
}

fn dashboard(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    let check = args.iter().any(|arg| arg == "--check");
    let features = read_features()?;
    let rendered = render_dashboard(&features);
    let path = Path::new("features/DASHBOARD.md");

    if check {
        let existing = fs::read_to_string(path)?;
        if existing != rendered {
            return Err(boxed_error(
                "features/DASHBOARD.md is out of date; run `cargo xtask dashboard`",
            ));
        }
    } else {
        fs::write(path, rendered)?;
    }
    Ok(())
}

fn validate(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    let feature = value_after_flag(&args, "--feature")
        .ok_or_else(|| boxed_error("missing required flag: --feature FEATURE_ID"))?;
    let features = read_features()?;
    if !features.iter().any(|candidate| candidate.id == feature) {
        return Err(boxed_error(format!("unknown feature: {feature}")));
    }

    println!("validation harness found feature `{feature}`");
    let manifest_path = validation_manifest_path(feature);
    if manifest_path.exists() {
        let manifest = read_validation_manifest(&manifest_path)?;
        if manifest.feature_id != feature {
            return Err(boxed_error(format!(
                "{} declares feature_id `{}`, expected `{feature}`",
                manifest_path.display(),
                manifest.feature_id
            )));
        }
        println!(
            "validation manifest uses {} {}",
            manifest.reference_tool, manifest.reference_version
        );
        validate_manifest_paths(&manifest_path, &manifest)?;
        println!(
            "validation manifest lists {} fixture(s)",
            manifest.fixtures.len()
        );
        let compared = validate_golden_outputs(&manifest_path, &manifest)?;
        if compared > 0 {
            println!("validation compared {compared} golden file(s)");
        }
    } else {
        println!("no reference validation manifest configured for `{feature}`");
    }
    Ok(())
}

fn list_features() -> Result<(), Box<dyn Error>> {
    for feature in read_features()? {
        println!(
            "{}\t{}\tv{}\timplemented={}\tvalidated={}",
            feature.id, feature.area, feature.version, feature.implemented, feature.validated
        );
    }
    Ok(())
}

fn skills(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    if args.iter().any(|arg| arg != "--check") {
        return Err(boxed_error("usage: cargo xtask skills --check"));
    }
    check_skills(Path::new(".codex/skills"))?;
    println!("repo-local feature skills are in sync");
    Ok(())
}

fn value_after_flag<'a>(args: &'a [String], flag: &str) -> Option<&'a str> {
    args.windows(2)
        .find(|window| window[0] == flag)
        .map(|window| window[1].as_str())
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Feature {
    id: String,
    title: String,
    area: String,
    version: u32,
    implemented: bool,
    validated: bool,
    description: String,
    depends_on: Vec<String>,
}

fn read_features() -> Result<Vec<Feature>, Box<dyn Error>> {
    read_features_from(Path::new("features"))
}

fn read_features_from(root: &Path) -> Result<Vec<Feature>, Box<dyn Error>> {
    let mut features = Vec::new();
    for entry in fs::read_dir(root)? {
        let entry = entry?;
        let path = entry.path();
        if !path.is_dir() || is_hidden_or_template(&path) {
            continue;
        }
        let metadata_path = path.join("feature.toml");
        if metadata_path.exists() {
            features.push(read_feature(&metadata_path)?);
        }
    }
    features.sort_by(|a, b| a.id.cmp(&b.id));
    validate_feature_set(&features)?;
    Ok(features)
}

fn is_hidden_or_template(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(|name| name.starts_with('_'))
        .unwrap_or(false)
}

fn read_feature(path: &Path) -> Result<Feature, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let map = parse_simple_toml(&text);
    reject_deprecated_feature_keys(&map, path)?;
    let feature = Feature {
        id: required(&map, "id", path)?,
        title: required(&map, "title", path)?,
        area: required(&map, "area", path)?,
        version: required_u32(&map, "version", path)?,
        implemented: required_bool(&map, "implemented", path)?,
        validated: required_bool(&map, "validated", path)?,
        description: required(&map, "description", path)?,
        depends_on: required_string_array(&map, "depends_on", path)?,
    };
    validate_feature(&feature, path)?;
    Ok(feature)
}

fn required(
    map: &BTreeMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<String, Box<dyn Error>> {
    map.get(key)
        .cloned()
        .ok_or_else(|| boxed_error(format!("{} is missing `{key}`", path.display())))
}

fn reject_deprecated_feature_keys(
    map: &BTreeMap<String, String>,
    path: &Path,
) -> Result<(), Box<dyn Error>> {
    for key in ["priority", "status", "last_ai_review"] {
        if map.contains_key(key) {
            return Err(boxed_error(format!(
                "{} uses deprecated feature metadata key `{key}`",
                path.display()
            )));
        }
    }
    Ok(())
}

fn parse_simple_toml(text: &str) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    let mut pending: Option<(String, String)> = None;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some((key, mut value)) = pending.take() {
            value.push(' ');
            value.push_str(line);
            if line.ends_with(']') {
                map.insert(key, normalize_value(&value));
            } else {
                pending = Some((key, value));
            }
            continue;
        }
        if let Some((key, value)) = line.split_once('=') {
            let key = key.trim().to_owned();
            let value = value.trim();
            if value.starts_with('[') && !value.ends_with(']') {
                pending = Some((key, value.to_owned()));
            } else {
                map.insert(key, normalize_value(value));
            }
        }
    }
    map
}

fn required_bool(
    map: &BTreeMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<bool, Box<dyn Error>> {
    match required(map, key, path)?.as_str() {
        "true" => Ok(true),
        "false" => Ok(false),
        value => Err(boxed_error(format!(
            "{} has invalid boolean `{key}` value `{value}`",
            path.display()
        ))),
    }
}

fn required_u32(
    map: &BTreeMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<u32, Box<dyn Error>> {
    let value = required(map, key, path)?;
    let parsed = value.parse::<u32>().map_err(|_| {
        boxed_error(format!(
            "{} has invalid integer `{key}` value `{value}`",
            path.display()
        ))
    })?;
    if parsed == 0 {
        return Err(boxed_error(format!(
            "{} has invalid zero `{key}` value",
            path.display()
        )));
    }
    Ok(parsed)
}

fn required_string_array(
    map: &BTreeMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<Vec<String>, Box<dyn Error>> {
    let value = required(map, key, path)?;
    parse_string_array(&value).ok_or_else(|| {
        boxed_error(format!(
            "{} has invalid string array `{key}` value `{value}`",
            path.display()
        ))
    })
}

fn parse_string_array(value: &str) -> Option<Vec<String>> {
    let value = value.trim();
    let inner = value.strip_prefix('[')?.strip_suffix(']')?.trim();
    if inner.is_empty() {
        return Some(Vec::new());
    }
    inner
        .split(',')
        .filter(|item| !item.trim().is_empty())
        .map(|item| {
            let item = item.trim();
            if item.starts_with('"') && item.ends_with('"') && item.len() >= 2 {
                Some(item[1..item.len() - 1].to_owned())
            } else {
                None
            }
        })
        .collect()
}

fn validate_feature(feature: &Feature, path: &Path) -> Result<(), Box<dyn Error>> {
    let expected_id = path
        .parent()
        .and_then(Path::file_name)
        .and_then(|name| name.to_str())
        .ok_or_else(|| {
            boxed_error(format!(
                "cannot determine feature directory for {}",
                path.display()
            ))
        })?;
    if feature.id != expected_id {
        return Err(boxed_error(format!(
            "{} declares id `{}`, expected `{expected_id}`",
            path.display(),
            feature.id
        )));
    }
    for (key, value) in [
        ("title", feature.title.as_str()),
        ("area", feature.area.as_str()),
        ("description", feature.description.as_str()),
    ] {
        if value.trim().is_empty() {
            return Err(boxed_error(format!(
                "{} has empty required field `{key}`",
                path.display()
            )));
        }
    }
    let feature_doc = path.with_file_name("feature.md");
    if !feature_doc.exists() {
        return Err(boxed_error(format!(
            "{} is missing required feature.md",
            path.parent().unwrap_or_else(|| Path::new(".")).display()
        )));
    }
    Ok(())
}

fn validate_feature_set(features: &[Feature]) -> Result<(), Box<dyn Error>> {
    let mut seen = BTreeMap::<&str, ()>::new();
    for feature in features {
        if seen.insert(feature.id.as_str(), ()).is_some() {
            return Err(boxed_error(format!(
                "duplicate feature id `{}`",
                feature.id
            )));
        }
    }
    for feature in features {
        for dependency in &feature.depends_on {
            if !seen.contains_key(dependency.as_str()) {
                return Err(boxed_error(format!(
                    "feature `{}` depends on unknown feature `{dependency}`",
                    feature.id
                )));
            }
        }
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ValidationManifest {
    feature_id: String,
    reference_tool: String,
    reference_version: String,
    fixtures: Vec<String>,
    fixture_sources: Vec<String>,
}

fn validation_manifest_path(feature: &str) -> PathBuf {
    Path::new("validation")
        .join("features")
        .join(feature)
        .join("validation.toml")
}

fn read_validation_manifest(path: &Path) -> Result<ValidationManifest, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let map = parse_simple_toml(&text);
    Ok(ValidationManifest {
        feature_id: required(&map, "feature_id", path)?,
        reference_tool: required(&map, "reference_tool", path)?,
        reference_version: required(&map, "reference_version", path)?,
        fixtures: optional_string_array(&map, "fixtures", path)?,
        fixture_sources: optional_string_array(&map, "fixture_sources", path)?,
    })
}

fn optional_string_array(
    map: &BTreeMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<Vec<String>, Box<dyn Error>> {
    match map.get(key) {
        Some(value) => parse_string_array(value).ok_or_else(|| {
            boxed_error(format!(
                "{} has invalid string array `{key}` value `{value}`",
                path.display()
            ))
        }),
        None => Ok(Vec::new()),
    }
}

fn validate_manifest_paths(
    manifest_path: &Path,
    manifest: &ValidationManifest,
) -> Result<(), Box<dyn Error>> {
    let base = manifest_path.parent().ok_or_else(|| {
        boxed_error(format!(
            "{} has no parent directory",
            manifest_path.display()
        ))
    })?;
    for fixture in &manifest.fixtures {
        let path = base.join(fixture);
        if !path.exists() {
            return Err(boxed_error(format!(
                "{} references missing fixture `{fixture}`",
                manifest_path.display()
            )));
        }
    }
    validate_fixture_sources(manifest_path, manifest)?;
    Ok(())
}

fn validate_golden_outputs(
    manifest_path: &Path,
    manifest: &ValidationManifest,
) -> Result<usize, Box<dyn Error>> {
    if manifest.fixtures.is_empty() {
        return Ok(0);
    }
    let base = manifest_path.parent().ok_or_else(|| {
        boxed_error(format!(
            "{} has no parent directory",
            manifest_path.display()
        ))
    })?;
    let mut compared = 0;
    for fixture in &manifest.fixtures {
        let fixture_path = base.join(fixture);
        let golden_path = base
            .join("golden")
            .join(format!("{}.json", slugify_fixture(fixture)));
        if !golden_path.exists() {
            return Err(boxed_error(format!(
                "{} is missing golden file for fixture `{fixture}`",
                manifest_path.display()
            )));
        }
        let golden: Value = serde_json::from_str(&fs::read_to_string(&golden_path)?)?;
        validate_golden_metadata(&golden_path, &golden, manifest, fixture)?;
        let expected = golden.get("expected").ok_or_else(|| {
            boxed_error(format!("{} is missing `expected`", golden_path.display()))
        })?;
        let actual = implementation_expected(&manifest.feature_id, &fixture_path)?;
        let expected = normalize_for_comparison(expected);
        let actual = normalize_for_comparison(&actual);
        if let Some(diff) = first_json_diff("$", &expected, &actual) {
            return Err(boxed_error(format!(
                "{} differs from implementation output for fixture `{fixture}`: {diff}",
                golden_path.display()
            )));
        }
        compared += 1;
    }
    Ok(compared)
}

fn validate_golden_metadata(
    golden_path: &Path,
    golden: &Value,
    manifest: &ValidationManifest,
    fixture: &str,
) -> Result<(), Box<dyn Error>> {
    if golden.get("schema_version") != Some(&json!(1)) {
        return Err(boxed_error(format!(
            "{} has unsupported schema_version",
            golden_path.display()
        )));
    }
    if golden.get("feature_id").and_then(Value::as_str) != Some(manifest.feature_id.as_str()) {
        return Err(boxed_error(format!(
            "{} feature_id does not match manifest",
            golden_path.display()
        )));
    }
    if golden.get("fixture_path").and_then(Value::as_str) != Some(fixture) {
        return Err(boxed_error(format!(
            "{} fixture_path does not match manifest",
            golden_path.display()
        )));
    }
    Ok(())
}

fn implementation_expected(feature: &str, fixture_path: &Path) -> Result<Value, Box<dyn Error>> {
    match feature {
        "io.sdf.v2000.parse" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(sdf_record_json).collect::<Vec<_>>() }))
        }
        "io.sdf.v2000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let molecules = records
                .into_iter()
                .map(|record| record.molecule)
                .collect::<Vec<_>>();
            let written = write_sdf_v2000(&molecules)?;
            let records = read_sdf_v2000_records(&written, SdfParseOptions::default())?
                .into_iter()
                .enumerate()
                .map(|(index, record)| IndexedSmallRecord {
                    record_index: index,
                    title: record.title,
                    molecule: record.molecule,
                })
                .collect::<Vec<_>>();
            Ok(json!({ "records": records.iter().map(sdf_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v2000.parse" | "core.conformers" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(conformer_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v2000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let records = records
                .into_iter()
                .enumerate()
                .map(|(index, record)| {
                    let written = write_mol_v2000(&record.molecule)?;
                    let molecule = read_mol_v2000_str(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title: molecule_title(&molecule.mol),
                        molecule,
                    })
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
            Ok(json!({ "records": records.iter().map(sdf_record_json).collect::<Vec<_>>() }))
        }
        "io.smiles.parse" => {
            let records = read_smiles_records(fixture_path)?;
            Ok(
                json!({ "records": records.iter().map(smiles_parse_record_json).collect::<Vec<_>>() }),
            )
        }
        "io.smiles.write" => {
            let records = read_smiles_records(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter()
                    .map(smiles_write_record_json)
                    .collect::<Result<Vec<_>, Box<dyn Error>>>()?
            }))
        }
        "algo.rings.fast" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(ring_membership_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.rings.sssr" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(ring_set_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.valence.rdkit-like" | "chem.sanitize.rdkit-like" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(sanitized_atom_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.aromaticity.rdkit-like" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(aromaticity_record_json).collect::<Vec<_>>() }),
            )
        }
        "io.mmcif.parse" | "bio.hierarchy.smcra" => {
            let input = fs::read_to_string(fixture_path)?;
            let molecule = read_mmcif_str(&input, MmcifParseOptions::default())?;
            Ok(mmcif_expected_json(&molecule))
        }
        _ => Err(boxed_error(format!(
            "no implementation comparison configured for feature `{feature}`"
        ))),
    }
}

#[derive(Debug, Clone)]
struct IndexedSmallRecord {
    record_index: usize,
    title: String,
    molecule: SmallMolecule,
}

#[derive(Debug, Clone)]
struct IndexedSmilesRecord {
    record_index: usize,
    title: String,
    input_smiles: String,
    molecule: SmallMolecule,
}

fn read_small_records_by_suffix(path: &Path) -> Result<Vec<IndexedSmallRecord>, Box<dyn Error>> {
    let input = fs::read_to_string(path)?;
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("mol" | "mdl")
    ) {
        let molecule = read_mol_v2000_str(&input)?;
        return Ok(vec![IndexedSmallRecord {
            record_index: 0,
            title: molecule_title(&molecule.mol),
            molecule,
        }]);
    }
    Ok(read_sdf_v2000_records(&input, SdfParseOptions::default())?
        .into_iter()
        .enumerate()
        .map(|(index, record)| small_record(index, record))
        .collect())
}

fn small_record(index: usize, record: SdfRecord) -> IndexedSmallRecord {
    IndexedSmallRecord {
        record_index: index,
        title: record.title,
        molecule: record.molecule,
    }
}

fn read_smiles_records(path: &Path) -> Result<Vec<IndexedSmilesRecord>, Box<dyn Error>> {
    let mut records = Vec::new();
    for (index, raw_line) in fs::read_to_string(path)?.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(2, char::is_whitespace);
        let smiles = parts.next().unwrap_or_default().to_owned();
        let title = parts.next().unwrap_or_default().trim().to_owned();
        let molecule = read_smiles_str(&smiles, SmilesParseOptions)?;
        records.push(IndexedSmilesRecord {
            record_index: index,
            title,
            input_smiles: smiles,
            molecule,
        });
    }
    Ok(records)
}

fn sdf_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": atoms_json(mol),
        "bonds": bonds_json(mol),
        "properties": sdf_properties_json(mol),
    })
}

fn conformer_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "conformers": mol.conformers().map(|(_, conformer)| {
            mol.atom_ids()
                .filter_map(|atom_id| {
                    conformer.position(atom_id).map(|point| json!({
                        "atom_index": atom_id.raw(),
                        "x": point.x,
                        "y": point.y,
                        "z": point.z,
                    }))
                })
                .collect::<Vec<_>>()
        }).collect::<Vec<_>>(),
        "atoms": atoms_json(mol),
    })
}

fn ring_membership_record_json(record: &mut IndexedSmallRecord) -> Value {
    let membership = perceive_ring_membership(&mut record.molecule.mol);
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_in_ring": mol.atom_ids().map(|id| membership.atom_in_ring(id)).collect::<Vec<_>>(),
        "bond_in_ring": mol.bond_ids().map(|id| membership.bond_in_ring(id)).collect::<Vec<_>>(),
    })
}

fn ring_set_record_json(record: &mut IndexedSmallRecord) -> Value {
    let ring_set = perceive_ring_set(&mut record.molecule.mol);
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "rings": ring_set
            .rings()
            .iter()
            .map(|ring| ring.atoms.iter().map(|atom| atom.raw()).collect::<Vec<_>>())
            .collect::<Vec<_>>(),
    })
}

fn sanitized_atom_record_json(record: &mut IndexedSmallRecord) -> Value {
    match sanitize_small_molecule(&mut record.molecule, SanitizeOptions::default()) {
        Ok(_) => json!({
            "record_index": record.record_index,
            "status": "ok",
            "title": record.title,
            "atoms": atoms_json(&record.molecule.mol),
        }),
        Err(_) => json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
        }),
    }
}

fn aromaticity_record_json(record: &mut IndexedSmallRecord) -> Value {
    let status = sanitize_small_molecule(
        &mut record.molecule,
        SanitizeOptions {
            perceive_valence: true,
            perceive_rings: true,
            perceive_aromaticity: false,
        },
    )
    .and_then(|_| {
        perceive_aromaticity(&mut record.molecule.mol, AromaticityModel::RdkitLike)
            .map_err(molecules::prelude::SanitizeError::Aromaticity)
    });
    if status.is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
        });
    }
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_aromatic": mol.atoms().map(|(_, atom)| atom.aromatic).collect::<Vec<_>>(),
        "bond_aromatic": mol.bonds().map(|(_, bond)| bond.aromatic).collect::<Vec<_>>(),
    })
}

fn smiles_write_record_json(record: &IndexedSmilesRecord) -> Result<Value, Box<dyn Error>> {
    Ok(json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "input_smiles": record.input_smiles,
        "output_smiles": write_smiles(&record.molecule, SmilesWriteOptions)?,
    }))
}

fn smiles_parse_record_json(record: &IndexedSmilesRecord) -> Value {
    sdf_record_json(&IndexedSmallRecord {
        record_index: record.record_index,
        title: record.title.clone(),
        molecule: record.molecule.clone(),
    })
}

fn atoms_json(mol: &Molecule) -> Vec<Value> {
    mol.atoms()
        .map(|(id, atom)| atom_json(id, atom))
        .collect::<Vec<_>>()
}

fn atom_json(id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "atom_map": atom.atom_map,
        "aromatic": atom.aromatic,
    })
}

fn bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| bond_json(id.raw(), bond))
        .collect::<Vec<_>>()
}

fn bond_json(index: u32, bond: &Bond) -> Value {
    json!({
        "index": index,
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": bond.aromatic,
        "stereo": bond_stereo_json(bond.stereo),
    })
}

fn bond_order_json(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Zero => "ZERO",
        BondOrder::Single => "SINGLE",
        BondOrder::Double => "DOUBLE",
        BondOrder::Triple => "TRIPLE",
        BondOrder::Quadruple => "QUADRUPLE",
        BondOrder::Aromatic => "AROMATIC",
        BondOrder::Dative => "DATIVE",
    }
}

fn bond_stereo_json(stereo: Option<BondStereo>) -> &'static str {
    match stereo {
        None | Some(BondStereo::Unspecified) => "STEREONONE",
        Some(BondStereo::E) => "STEREOE",
        Some(BondStereo::Z) => "STEREOZ",
        Some(BondStereo::Up) => "STEREOANY",
        Some(BondStereo::Down) => "STEREOANY",
    }
}

fn sdf_properties_json(mol: &Molecule) -> Value {
    let mut props = serde_json::Map::new();
    for (key, value) in mol.props() {
        let Some(name) = key.strip_prefix("sdf.field.") else {
            continue;
        };
        if let PropValue::String(text) = value {
            props.insert(name.to_owned(), json!(text));
        }
    }
    Value::Object(props)
}

fn molecule_title(mol: &Molecule) -> String {
    match mol.props().get("sdf.title") {
        Some(PropValue::String(title)) => title.clone(),
        _ => String::new(),
    }
}

fn mmcif_expected_json(molecule: &MacroMolecule) -> Value {
    json!({
        "atom_site_rows": atom_site_rows_json(molecule),
        "structure": structure_json(molecule),
    })
}

fn atom_site_rows_json(molecule: &MacroMolecule) -> Value {
    let rows = molecule
        .hierarchy
        .atom_sites()
        .map(|(site_id, site)| {
            let residue = molecule
                .hierarchy
                .residue(site.residue)
                .expect("residue exists");
            let chain = molecule
                .hierarchy
                .chain(residue.chain)
                .expect("chain exists");
            let model = molecule.hierarchy.model(chain.model).expect("model exists");
            let atom = molecule.mol.atom(site.atom).expect("atom exists");
            json!({
                "group_PDB": Value::Null,
                "id": (site_id.raw() + 1).to_string(),
                "type_symbol": atom.element.symbol(),
                "label_atom_id": site.metadata.label_atom_id,
                "auth_atom_id": site.metadata.auth_atom_id,
                "label_alt_id": site.metadata.label_alt_id,
                "label_comp_id": residue.name,
                "auth_comp_id": residue.name,
                "label_asym_id": chain.label_id,
                "auth_asym_id": chain.author_id,
                "label_seq_id": residue.label_seq_id.map(|value| value.to_string()),
                "auth_seq_id": residue.author_seq_id,
                "pdbx_PDB_ins_code": residue.insertion_code,
                "occupancy": site.metadata.occupancy.map(|value| format!("{value:.2}")),
                "B_iso_or_equiv": site.metadata.b_factor.map(|value| format!("{value:.2}")),
                "Cartn_x": Value::Null,
                "Cartn_y": Value::Null,
                "Cartn_z": Value::Null,
                "pdbx_PDB_model_num": model.model_id,
            })
        })
        .collect::<Vec<_>>();
    json!({
        "status": "ok",
        "row_count": rows.len(),
        "rows": rows,
    })
}

fn structure_json(molecule: &MacroMolecule) -> Value {
    json!({
        "status": "ok",
        "models": molecule.hierarchy.models().map(|(model_id, model)| {
            json!({
                "id": model_id.raw(),
                "chains": model.chains.iter().map(|chain_id| {
                    let chain = molecule.hierarchy.chain(*chain_id).expect("chain exists");
                    json!({
                        "id": chain.label_id,
                        "residues": chain.residues.iter().map(|residue_id| {
                            let residue = molecule.hierarchy.residue(*residue_id).expect("residue exists");
                            json!({
                                "name": residue.name,
                                "hetflag": Value::Null,
                                "sequence_id": residue.label_seq_id,
                                "insertion_code": residue.insertion_code,
                                "atoms": residue.atom_sites.iter().map(|site_id| {
                                    let site = molecule.hierarchy.atom_site(*site_id).expect("site exists");
                                    let atom = molecule.mol.atom(site.atom).expect("atom exists");
                                    let name = site
                                        .metadata
                                        .label_atom_id
                                        .clone()
                                        .unwrap_or_else(|| atom.element.symbol().to_owned());
                                    json!({
                                        "name": name,
                                        "full_name": name,
                                        "altloc": site.metadata.label_alt_id,
                                        "element": atom.element.symbol(),
                                        "occupancy": site.metadata.occupancy,
                                        "bfactor": site.metadata.b_factor,
                                        "coord": Value::Null,
                                    })
                                }).collect::<Vec<_>>(),
                            })
                        }).collect::<Vec<_>>(),
                    })
                }).collect::<Vec<_>>(),
            })
        }).collect::<Vec<_>>(),
    })
}

fn first_json_diff(path: &str, expected: &Value, actual: &Value) -> Option<String> {
    match (expected, actual) {
        (Value::Object(expected), Value::Object(actual)) => {
            for key in expected.keys() {
                let next = format!("{path}.{key}");
                let Some(actual_value) = actual.get(key) else {
                    return Some(format!("{next} missing from actual output"));
                };
                if let Some(diff) = first_json_diff(&next, &expected[key], actual_value) {
                    return Some(diff);
                }
            }
            for key in actual.keys() {
                if !expected.contains_key(key) {
                    return Some(format!("{path}.{key} present only in actual output"));
                }
            }
            None
        }
        (Value::Array(expected), Value::Array(actual)) => {
            if expected.len() != actual.len() {
                return Some(format!(
                    "{path} length differs: expected {}, actual {}",
                    expected.len(),
                    actual.len()
                ));
            }
            for (index, (expected_value, actual_value)) in expected.iter().zip(actual).enumerate() {
                if let Some(diff) =
                    first_json_diff(&format!("{path}[{index}]"), expected_value, actual_value)
                {
                    return Some(diff);
                }
            }
            None
        }
        _ if expected == actual => None,
        _ => Some(format!(
            "{path} differs: expected {}, actual {}",
            expected, actual
        )),
    }
}

fn normalize_for_comparison(value: &Value) -> Value {
    match value {
        Value::Array(items) => Value::Array(
            items
                .iter()
                .map(normalize_for_comparison)
                .collect::<Vec<_>>(),
        ),
        Value::Object(object) => {
            let mut normalized = serde_json::Map::new();
            for (key, value) in object {
                normalized.insert(key.clone(), normalize_for_comparison(value));
            }
            normalize_undirected_bond_object(&mut normalized);
            normalize_bond_array_object(&mut normalized);
            normalize_ring_set_object(&mut normalized);
            Value::Object(normalized)
        }
        _ => value.clone(),
    }
}

fn normalize_undirected_bond_object(object: &mut serde_json::Map<String, Value>) {
    let Some(begin) = object.get("begin_atom_index").and_then(Value::as_u64) else {
        return;
    };
    let Some(end) = object.get("end_atom_index").and_then(Value::as_u64) else {
        return;
    };
    if begin > end {
        object.insert("begin_atom_index".to_owned(), json!(end));
        object.insert("end_atom_index".to_owned(), json!(begin));
    }
}

fn normalize_bond_array_object(object: &mut serde_json::Map<String, Value>) {
    let Some(Value::Array(bonds)) = object.get_mut("bonds") else {
        return;
    };
    for bond in bonds.iter_mut() {
        if let Value::Object(bond) = bond {
            bond.remove("index");
        }
    }
    bonds.sort_by_key(bond_sort_key);
    for (index, bond) in bonds.iter_mut().enumerate() {
        if let Value::Object(bond) = bond {
            bond.insert("index".to_owned(), json!(index));
        }
    }
}

fn bond_sort_key(value: &Value) -> (u64, u64, String, String) {
    let Some(object) = value.as_object() else {
        return (u64::MAX, u64::MAX, String::new(), String::new());
    };
    (
        object
            .get("begin_atom_index")
            .and_then(Value::as_u64)
            .unwrap_or(u64::MAX),
        object
            .get("end_atom_index")
            .and_then(Value::as_u64)
            .unwrap_or(u64::MAX),
        object
            .get("bond_type")
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_owned(),
        object
            .get("stereo")
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_owned(),
    )
}

fn normalize_ring_set_object(object: &mut serde_json::Map<String, Value>) {
    let Some(Value::Array(rings)) = object.get_mut("rings") else {
        return;
    };
    for ring in rings.iter_mut() {
        let Value::Array(atoms) = ring else {
            continue;
        };
        atoms.sort_by_key(|value| value.as_u64().unwrap_or(u64::MAX));
    }
    rings.sort_by(|left, right| {
        let left = left
            .as_array()
            .map(|items| items.iter().filter_map(Value::as_u64).collect::<Vec<_>>())
            .unwrap_or_default();
        let right = right
            .as_array()
            .map(|items| items.iter().filter_map(Value::as_u64).collect::<Vec<_>>())
            .unwrap_or_default();
        left.cmp(&right)
    });
}

fn slugify_fixture(fixture: &str) -> String {
    let mut slug = String::new();
    let mut previous_was_separator = false;
    for ch in fixture.chars() {
        if ch.is_ascii_alphanumeric() || matches!(ch, '.' | '-' | '_') {
            slug.push(ch);
            previous_was_separator = false;
        } else if !previous_was_separator {
            slug.push('_');
            previous_was_separator = true;
        }
    }
    slug.trim_matches(['.', '_', '-']).to_owned()
}

fn validate_fixture_sources(
    manifest_path: &Path,
    manifest: &ValidationManifest,
) -> Result<(), Box<dyn Error>> {
    if manifest.fixtures.is_empty() {
        if !manifest.fixture_sources.is_empty() {
            return Err(boxed_error(format!(
                "{} declares fixture_sources without fixtures",
                manifest_path.display()
            )));
        }
        return Ok(());
    }
    if manifest.fixture_sources.len() != manifest.fixtures.len() {
        return Err(boxed_error(format!(
            "{} must declare one external fixture_sources entry for each fixture",
            manifest_path.display()
        )));
    }
    for fixture in &manifest.fixtures {
        let source = manifest
            .fixture_sources
            .iter()
            .find(|entry| fixture_source_path(entry).as_deref() == Some(fixture.as_str()))
            .ok_or_else(|| {
                boxed_error(format!(
                    "{} is missing external provenance for fixture `{fixture}`",
                    manifest_path.display()
                ))
            })?;
        validate_fixture_source_entry(manifest_path, fixture, source)?;
    }
    Ok(())
}

fn fixture_source_path(entry: &str) -> Option<String> {
    entry.split('|').next().map(|part| part.trim().to_owned())
}

fn validate_fixture_source_entry(
    manifest_path: &Path,
    fixture: &str,
    source: &str,
) -> Result<(), Box<dyn Error>> {
    let mut has_source = false;
    let mut has_url = false;
    let mut has_sha256 = false;
    for part in source.split('|').skip(1) {
        let Some((key, value)) = part.trim().split_once('=') else {
            continue;
        };
        let key = key.trim();
        let value = value.trim();
        match key {
            "source" if !value.is_empty() && !value.eq_ignore_ascii_case("manual") => {
                has_source = true;
            }
            "url" if value.starts_with("https://") => has_url = true,
            "sha256" if is_sha256(value) => has_sha256 = true,
            _ => {}
        }
    }
    if has_source && has_url && has_sha256 {
        Ok(())
    } else {
        Err(boxed_error(format!(
            "{} provenance for fixture `{fixture}` must include non-manual source, https url, and sha256",
            manifest_path.display()
        )))
    }
}

fn is_sha256(value: &str) -> bool {
    value.len() == 64 && value.bytes().all(|byte| byte.is_ascii_hexdigit())
}

fn normalize_value(value: &str) -> String {
    let value = value.trim().trim_end_matches(',').trim();
    if value.starts_with('"') && value.ends_with('"') && value.len() >= 2 {
        value[1..value.len() - 1].to_owned()
    } else {
        value.to_owned()
    }
}

fn render_dashboard(features: &[Feature]) -> String {
    let mut out = String::new();
    out.push_str("# Feature Dashboard\n\n");
    out.push_str("Generated from `features/*/feature.toml`. Do not hand-edit this file.\n\n");
    out.push_str("| Feature | Title | Area | Version | Implemented | Validated |\n");
    out.push_str("|---|---|---|---:|:---:|:---:|\n");
    for feature in features {
        out.push_str(&format!(
            "| `{}` | {} | {} | {} | {} | {} |\n",
            feature.id,
            feature.title,
            feature.area,
            feature.version,
            checkmark(feature.implemented),
            checkmark(feature.validated)
        ));
    }
    out
}

fn checkmark(value: bool) -> &'static str {
    if value {
        "✅"
    } else {
        "❌"
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SkillMetadata {
    name: String,
    description: String,
}

fn check_skills(root: &Path) -> Result<(), Box<dyn Error>> {
    for expected in expected_skills() {
        let path = root.join(expected.name).join("SKILL.md");
        if !path.exists() {
            return Err(boxed_error(format!(
                "missing repo-local skill `{}` at {}",
                expected.name,
                path.display()
            )));
        }
        let text = fs::read_to_string(&path)?;
        let metadata = parse_skill_metadata(&text, &path)?;
        if metadata.name != expected.name {
            return Err(boxed_error(format!(
                "{} declares skill name `{}`, expected `{}`",
                path.display(),
                metadata.name,
                expected.name
            )));
        }
        let lower = text.to_lowercase();
        for required in expected.required_phrases {
            if !lower.contains(&required.to_lowercase()) {
                return Err(boxed_error(format!(
                    "{} is missing required phrase `{required}`",
                    path.display()
                )));
            }
        }
    }
    Ok(())
}

struct ExpectedSkill {
    name: &'static str,
    required_phrases: &'static [&'static str],
}

fn expected_skills() -> &'static [ExpectedSkill] {
    &[
        ExpectedSkill {
            name: "feature-work",
            required_phrases: &[
                "add -> optional research -> plan -> implement",
                "feature.md",
                "implemented = true",
                "validated = true",
                "externally supplied",
                "cargo xtask dashboard --check",
                "cargo xtask validate --feature",
            ],
        },
        ExpectedSkill {
            name: "feature-review",
            required_phrases: &[
                "independent audit",
                "architecture",
                "validation claims",
                "feature.md",
                "cargo test --workspace",
                "cargo xtask validate --feature",
            ],
        },
    ]
}

fn parse_skill_metadata(text: &str, path: &Path) -> Result<SkillMetadata, Box<dyn Error>> {
    let mut lines = text.lines();
    if lines.next() != Some("---") {
        return Err(boxed_error(format!(
            "{} is missing YAML frontmatter",
            path.display()
        )));
    }
    let mut fields = BTreeMap::<String, String>::new();
    for line in lines.by_ref() {
        if line == "---" {
            let name = fields
                .get("name")
                .cloned()
                .ok_or_else(|| boxed_error(format!("{} is missing `name`", path.display())))?;
            let description = fields.get("description").cloned().ok_or_else(|| {
                boxed_error(format!("{} is missing `description`", path.display()))
            })?;
            if name.trim().is_empty() || description.trim().is_empty() {
                return Err(boxed_error(format!(
                    "{} has empty skill frontmatter",
                    path.display()
                )));
            }
            return Ok(SkillMetadata { name, description });
        }
        if let Some((key, value)) = line.split_once(':') {
            fields.insert(key.trim().to_owned(), value.trim().to_owned());
        }
    }
    Err(boxed_error(format!(
        "{} has unterminated YAML frontmatter",
        path.display()
    )))
}

fn boxed_error(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::other(message.into()).into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn value_after_flag_finds_following_value() {
        let args = vec![
            "validate".to_owned(),
            "--feature".to_owned(),
            "core.graph".to_owned(),
        ];

        assert_eq!(value_after_flag(&args, "--feature"), Some("core.graph"));
    }

    #[test]
    fn simple_toml_parser_normalizes_strings_booleans_and_arrays() {
        let parsed = parse_simple_toml(
            r#"
            id = "core.graph"
            implemented = true
            depends_on = []
            "#,
        );

        assert_eq!(parsed.get("id"), Some(&"core.graph".to_owned()));
        assert_eq!(parsed.get("implemented"), Some(&"true".to_owned()));
        assert_eq!(parsed.get("depends_on"), Some(&"[]".to_owned()));
    }

    #[test]
    fn string_array_parser_accepts_empty_and_string_lists() {
        assert_eq!(parse_string_array("[]"), Some(Vec::new()));
        assert_eq!(
            parse_string_array(r#"["core.graph", "validation.harness"]"#),
            Some(vec![
                "core.graph".to_owned(),
                "validation.harness".to_owned()
            ])
        );
        assert_eq!(parse_string_array("[core.graph]"), None);
    }

    #[test]
    fn simple_toml_parser_accepts_multiline_string_arrays() {
        let parsed = parse_simple_toml(
            r#"
            feature_id = "io.sdf.v2000.parse"
            fixtures = [
              "fixtures/a.sdf",
              "fixtures/b.sdf",
            ]
            "#,
        );

        assert_eq!(
            parse_string_array(parsed.get("fixtures").expect("fixtures should parse")),
            Some(vec![
                "fixtures/a.sdf".to_owned(),
                "fixtures/b.sdf".to_owned()
            ])
        );
    }

    #[test]
    fn read_feature_parses_typed_metadata() {
        let root = temp_feature_root("read-feature");
        write_feature(
            &root,
            "example.feature",
            r#"id = "example.feature"
title = "Example"
area = "infrastructure"
version = 2
implemented = false
validated = true
description = "Example feature."
depends_on = ["core.graph"]
"#,
        );

        let feature = read_feature(&root.join("example.feature").join("feature.toml"))
            .expect("feature should parse");

        assert_eq!(feature.id, "example.feature");
        assert_eq!(feature.version, 2);
        assert!(!feature.implemented);
        assert!(feature.validated);
        assert_eq!(feature.depends_on, vec!["core.graph"]);
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn read_feature_rejects_bad_boolean_deprecated_keys_missing_docs_and_directory_mismatch() {
        let root = temp_feature_root("bad-feature");
        write_feature(
            &root,
            "bad.bool",
            r#"id = "bad.bool"
title = "Bad"
area = "infrastructure"
version = 1
implemented = maybe
validated = false
description = "Bad feature."
depends_on = []
"#,
        );
        assert!(read_feature(&root.join("bad.bool").join("feature.toml")).is_err());

        write_feature(
            &root,
            "bad.deprecated",
            r#"id = "bad.deprecated"
title = "Bad"
area = "infrastructure"
version = 1
priority = "P0"
implemented = false
validated = false
description = "Bad feature."
depends_on = []
"#,
        );
        assert!(read_feature(&root.join("bad.deprecated").join("feature.toml")).is_err());

        write_feature(
            &root,
            "bad.version",
            r#"id = "bad.version"
title = "Bad"
area = "infrastructure"
version = 0
implemented = false
validated = false
description = "Bad feature."
depends_on = []
"#,
        );
        assert!(read_feature(&root.join("bad.version").join("feature.toml")).is_err());

        write_feature_without_doc(
            &root,
            "missing.doc",
            r#"id = "missing.doc"
title = "Bad"
area = "infrastructure"
version = 1
implemented = false
validated = false
description = "Bad feature."
depends_on = []
"#,
        );
        assert!(read_feature(&root.join("missing.doc").join("feature.toml")).is_err());

        write_feature(
            &root,
            "real.id",
            r#"id = "wrong.id"
title = "Bad"
area = "infrastructure"
version = 1
implemented = false
validated = false
description = "Bad feature."
depends_on = []
"#,
        );
        assert!(read_feature(&root.join("real.id").join("feature.toml")).is_err());
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn read_features_sorts_skips_templates_and_validates_dependencies() {
        let root = temp_feature_root("feature-set");
        write_feature(
            &root,
            "z.feature",
            r#"id = "z.feature"
title = "Zed"
area = "core"
version = 1
implemented = true
validated = false
description = "Z feature."
depends_on = ["a.feature"]
"#,
        );
        write_feature(
            &root,
            "a.feature",
            r#"id = "a.feature"
title = "Aye"
area = "core"
version = 1
implemented = false
validated = false
description = "A feature."
depends_on = []
"#,
        );
        fs::create_dir_all(root.join("_template")).expect("template dir should create");
        fs::write(root.join("_template").join("feature.toml"), "not = valid")
            .expect("template metadata should write");

        let features = read_features_from(&root).expect("feature set should parse");

        assert_eq!(
            features
                .iter()
                .map(|feature| feature.id.as_str())
                .collect::<Vec<_>>(),
            vec!["a.feature", "z.feature"]
        );

        write_feature(
            &root,
            "bad.dependency",
            r#"id = "bad.dependency"
title = "Bad"
area = "core"
version = 1
implemented = false
validated = false
description = "Bad dependency."
depends_on = ["missing.feature"]
"#,
        );
        assert!(read_features_from(&root).is_err());
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn render_dashboard_is_stable_and_uses_boolean_labels() {
        let features = vec![
            Feature {
                id: "a.feature".to_owned(),
                title: "Aye".to_owned(),
                area: "core".to_owned(),
                version: 1,
                implemented: false,
                validated: false,
                description: "A feature.".to_owned(),
                depends_on: Vec::new(),
            },
            Feature {
                id: "z.feature".to_owned(),
                title: "Zed".to_owned(),
                area: "io".to_owned(),
                version: 3,
                implemented: true,
                validated: false,
                description: "Z feature.".to_owned(),
                depends_on: vec!["a.feature".to_owned()],
            },
        ];

        let dashboard = render_dashboard(&features);

        assert!(
            dashboard.contains("| Feature | Title | Area | Version | Implemented | Validated |")
        );
        assert!(dashboard.contains("| `a.feature` | Aye | core | 1 | ❌ | ❌ |"));
        assert!(dashboard.contains("| `z.feature` | Zed | io | 3 | ✅ | ❌ |"));
        assert!(dashboard.ends_with('\n'));
    }

    #[test]
    fn skill_metadata_parser_and_check_validate_repo_skill_contract() {
        let root = temp_feature_root("skills-check");
        write_skill(
            &root,
            "feature-work",
            r#"---
name: feature-work
description: Builder skill.
---
# Feature Work
add -> optional research -> plan -> implement
Use feature.md. Only set implemented = true and validated = true with evidence.
Molecular validation fixtures must be externally supplied.
Run cargo xtask dashboard --check and cargo xtask validate --feature <feature-id>.
"#,
        );
        write_skill(
            &root,
            "feature-review",
            r#"---
name: feature-review
description: Review skill.
---
# Feature Review
Independent audit for architecture and validation claims.
Read feature.md. Run cargo test --workspace and cargo xtask validate --feature <feature-id>.
"#,
        );

        check_skills(&root).expect("skills should pass");

        fs::write(
            root.join("feature-review").join("SKILL.md"),
            "# Missing frontmatter",
        )
        .expect("skill should rewrite");
        assert!(check_skills(&root).is_err());
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn validation_manifest_path_is_feature_scoped() {
        assert_eq!(
            validation_manifest_path("core.graph"),
            PathBuf::from("validation/features/core.graph/validation.toml")
        );
    }

    #[test]
    fn validation_manifest_reads_and_checks_fixture_paths() {
        let root = temp_feature_root("validation-manifest");
        let feature_dir = root.join("validation").join("features").join("example");
        let fixture_dir = feature_dir.join("fixtures");
        fs::create_dir_all(&fixture_dir).expect("fixture dir should create");
        fs::write(fixture_dir.join("ok.txt"), "{}").expect("fixture should write");
        let manifest_path = feature_dir.join("validation.toml");
        fs::write(
            &manifest_path,
            r#"feature_id = "example"
reference_tool = "manual-fixtures"
reference_version = "test"
fixtures = [
  "fixtures/ok.txt",
]
fixture_sources = [
  "fixtures/ok.txt | source=Example External Source | url=https://example.org/ok.txt | sha256=0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
]
"#,
        )
        .expect("manifest should write");

        let manifest = read_validation_manifest(&manifest_path).expect("manifest should parse");
        assert_eq!(manifest.fixtures, vec!["fixtures/ok.txt"]);
        assert_eq!(manifest.fixture_sources.len(), 1);
        validate_manifest_paths(&manifest_path, &manifest).expect("fixture should exist");

        fs::write(
            &manifest_path,
            r#"feature_id = "example"
reference_tool = "manual-fixtures"
reference_version = "test"
fixtures = [
  "fixtures/missing.txt",
]
fixture_sources = [
  "fixtures/missing.txt | source=Example External Source | url=https://example.org/missing.txt | sha256=0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
]
"#,
        )
        .expect("manifest should rewrite");
        let manifest = read_validation_manifest(&manifest_path).expect("manifest should parse");
        assert!(validate_manifest_paths(&manifest_path, &manifest).is_err());

        fs::write(
            &manifest_path,
            r#"feature_id = "example"
reference_tool = "manual-fixtures"
reference_version = "test"
fixtures = [
  "fixtures/ok.txt",
]
"#,
        )
        .expect("manifest should rewrite");
        let manifest = read_validation_manifest(&manifest_path).expect("manifest should parse");
        assert!(validate_manifest_paths(&manifest_path, &manifest).is_err());
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn validation_requires_matching_golden_files() {
        let root = temp_feature_root("validation-goldens");
        let feature_dir = root
            .join("validation")
            .join("features")
            .join("io.smiles.parse");
        let fixture_dir = feature_dir.join("fixtures");
        fs::create_dir_all(&fixture_dir).expect("fixture dir should create");
        fs::write(fixture_dir.join("case.smi"), "CCO\n").expect("fixture should write");
        let manifest_path = feature_dir.join("validation.toml");
        fs::write(
            &manifest_path,
            r#"feature_id = "io.smiles.parse"
reference_tool = "rdkit"
reference_version = "test"
fixtures = [
  "fixtures/case.smi",
]
fixture_sources = [
  "fixtures/case.smi | source=PubChem | url=https://example.org/case.smi | sha256=0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
]
"#,
        )
        .expect("manifest should write");

        let manifest = read_validation_manifest(&manifest_path).expect("manifest should parse");
        validate_manifest_paths(&manifest_path, &manifest).expect("manifest paths should pass");
        assert!(validate_golden_outputs(&manifest_path, &manifest).is_err());

        fs::create_dir_all(feature_dir.join("golden")).expect("golden dir should create");
        fs::write(
            feature_dir.join("golden").join("fixtures_case.smi.json"),
            r#"{
  "schema_version": 1,
  "feature_id": "wrong.feature",
  "fixture_path": "fixtures/case.smi",
  "expected": {}
}
"#,
        )
        .expect("golden should write");
        assert!(validate_golden_outputs(&manifest_path, &manifest).is_err());
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn comparison_normalizes_undirected_bonds_and_ring_order() {
        let expected = json!({
            "records": [{
                "bonds": [
                    {"index": 0, "begin_atom_index": 5, "end_atom_index": 0, "bond_type": "SINGLE", "stereo": "STEREONONE"}
                ],
                "rings": [[5, 3, 1]]
            }]
        });
        let actual = json!({
            "records": [{
                "bonds": [
                    {"index": 7, "begin_atom_index": 0, "end_atom_index": 5, "bond_type": "SINGLE", "stereo": "STEREONONE"}
                ],
                "rings": [[1, 3, 5]]
            }]
        });

        assert_eq!(
            normalize_for_comparison(&expected),
            normalize_for_comparison(&actual)
        );
    }

    fn temp_feature_root(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time should be available")
            .as_nanos();
        let root =
            env::temp_dir().join(format!("molecules-xtask-{label}-{}-{nonce}", process::id()));
        fs::create_dir_all(&root).expect("temp feature root should create");
        root
    }

    fn write_feature(root: &Path, id: &str, metadata: &str) {
        let dir = root.join(id);
        fs::create_dir_all(&dir).expect("feature dir should create");
        fs::write(dir.join("feature.toml"), metadata).expect("feature metadata should write");
        fs::write(dir.join("feature.md"), "# Feature\n").expect("feature doc should write");
    }

    fn write_feature_without_doc(root: &Path, id: &str, metadata: &str) {
        let dir = root.join(id);
        fs::create_dir_all(&dir).expect("feature dir should create");
        fs::write(dir.join("feature.toml"), metadata).expect("feature metadata should write");
    }

    fn write_skill(root: &Path, name: &str, text: &str) {
        let dir = root.join(name);
        fs::create_dir_all(&dir).expect("skill dir should create");
        fs::write(dir.join("SKILL.md"), text).expect("skill should write");
    }
}
