use crate::*;

pub(crate) fn first_json_diff(
    feature_id: &str,
    path: &str,
    expected: &Value,
    actual: &Value,
) -> Option<String> {
    match (expected, actual) {
        (Value::Object(expected), Value::Object(actual)) => {
            for key in expected.keys() {
                let next = format!("{path}.{key}");
                let Some(actual_value) = actual.get(key) else {
                    return Some(format!("{next} missing from actual output"));
                };
                if let Some(diff) = first_json_diff(feature_id, &next, &expected[key], actual_value)
                {
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
                if let Some(diff) = first_json_diff(
                    feature_id,
                    &format!("{path}[{index}]"),
                    expected_value,
                    actual_value,
                ) {
                    return Some(diff);
                }
            }
            None
        }
        (Value::Number(expected), Value::Number(actual))
            if path.contains(".coord[")
                && expected
                    .as_f64()
                    .zip(actual.as_f64())
                    .map(|(expected, actual)| (expected - actual).abs() <= 0.0015)
                    .unwrap_or(false) =>
        {
            None
        }
        (Value::Number(expected), Value::Number(actual))
            if feature_id == "bio.secondary-structure.dssp"
                && (path.ends_with(".phi_degrees")
                    || path.ends_with(".psi_degrees")
                    || path.ends_with(".kappa_degrees")
                    || path.ends_with(".alpha_degrees"))
                && expected
                    .as_f64()
                    .zip(actual.as_f64())
                    .map(|(expected, actual)| (expected - actual).abs() <= 0.15)
                    .unwrap_or(false) =>
        {
            None
        }
        (Value::Number(expected), Value::Number(actual))
            if feature_id == "descriptor.molecular"
                && path.ends_with(".average_mass_da")
                && expected
                    .as_f64()
                    .zip(actual.as_f64())
                    .map(|(expected, actual)| (expected - actual).abs() <= 0.05)
                    .unwrap_or(false) =>
        {
            None
        }
        (Value::Number(expected), Value::Number(actual))
            if feature_id == "descriptor.molecular"
                && path.ends_with(".monoisotopic_mass_da")
                && expected
                    .as_f64()
                    .zip(actual.as_f64())
                    .map(|(expected, actual)| (expected - actual).abs() <= 5.0e-5)
                    .unwrap_or(false) =>
        {
            None
        }
        (Value::Number(expected), Value::Number(actual))
            if feature_id == "bio.secondary-structure.dssp"
                && path.ends_with(".tco")
                && expected
                    .as_f64()
                    .zip(actual.as_f64())
                    .map(|(expected, actual)| (expected - actual).abs() <= 0.0015)
                    .unwrap_or(false) =>
        {
            None
        }
        (Value::Number(expected), Value::Number(actual))
            if feature_id == "bio.secondary-structure.dssp"
                && path.ends_with(".energy_kcal_per_mol")
                && expected
                    .as_f64()
                    .zip(actual.as_f64())
                    .map(|(expected, actual)| (expected - actual).abs() <= 0.051)
                    .unwrap_or(false) =>
        {
            None
        }
        _ if expected == actual => None,
        _ => Some(format!(
            "{path} differs: expected {}, actual {}",
            expected, actual
        )),
    }
}

pub(crate) fn normalize_feature_for_comparison_in_place(feature_id: &str, value: &mut Value) {
    normalize_for_comparison_in_place(value);
    if feature_id == "bio.secondary-structure.dssp" {
        normalize_dssp_residue_order(value);
    }
}

fn normalize_dssp_residue_order(value: &mut Value) {
    let Some(residues) = value
        .as_object_mut()
        .and_then(|object| object.get_mut("residues"))
        .and_then(Value::as_array_mut)
    else {
        return;
    };
    residues.sort_by_key(dssp_residue_sort_key);

    let mut sheets = BTreeMap::new();
    let mut strands = BTreeMap::new();
    let mut ladders = BTreeMap::new();
    for residue in residues {
        let Some(residue) = residue.as_object_mut() else {
            continue;
        };
        canonicalize_dssp_id(residue.get_mut("sheet"), &mut sheets);
        canonicalize_dssp_id(residue.get_mut("strand"), &mut strands);
        if let Some(values) = residue.get_mut("ladders").and_then(Value::as_array_mut) {
            for value in values {
                canonicalize_dssp_id(Some(value), &mut ladders);
            }
        }
    }
}

fn canonicalize_dssp_id(value: Option<&mut Value>, ids: &mut BTreeMap<i64, i64>) {
    let Some(value) = value else {
        return;
    };
    let Some(id) = value.as_i64() else {
        return;
    };
    let next = ids.len() as i64;
    let canonical = *ids.entry(id).or_insert(next);
    *value = json!(canonical);
}

fn dssp_residue_sort_key(value: &Value) -> (String, i64, String, String, i64, String) {
    let field = |name: &str| value.get(name);
    (
        field("chain_id")
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_owned(),
        field("sequence_id")
            .and_then(Value::as_i64)
            .unwrap_or(i64::MIN),
        field("insertion_code")
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_owned(),
        field("label_chain_id")
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_owned(),
        field("label_sequence_id")
            .and_then(Value::as_i64)
            .unwrap_or(i64::MIN),
        field("residue_name")
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_owned(),
    )
}

pub(crate) fn normalize_for_comparison_in_place(value: &mut Value) {
    match value {
        Value::Array(items) => {
            for item in items {
                normalize_for_comparison_in_place(item);
            }
        }
        Value::Object(object) => {
            for value in object.values_mut() {
                normalize_for_comparison_in_place(value);
            }
            normalize_undirected_bond_object(object);
            normalize_bond_array_object(object);
            normalize_ring_set_object(object);
            normalize_coord_object(object);
        }
        _ => {}
    }
}

#[cfg(test)]
pub(crate) fn normalize_for_comparison(value: &Value) -> Value {
    let mut normalized = value.clone();
    normalize_for_comparison_in_place(&mut normalized);
    normalized
}

pub(crate) fn normalize_coord_object(object: &mut serde_json::Map<String, Value>) {
    let Some(Value::Array(coord)) = object.get_mut("coord") else {
        return;
    };
    for value in coord.iter_mut() {
        if let Some(number) = value.as_f64() {
            *value = json!((number * 1000.0).round() / 1000.0);
        }
    }
}

pub(crate) fn normalize_undirected_bond_object(object: &mut serde_json::Map<String, Value>) {
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

pub(crate) fn normalize_bond_array_object(object: &mut serde_json::Map<String, Value>) {
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

pub(crate) fn bond_sort_key(value: &Value) -> (u64, u64, String, String) {
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

pub(crate) fn normalize_ring_set_object(object: &mut serde_json::Map<String, Value>) {
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

pub(crate) fn slugify_fixture(fixture: &str) -> String {
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

pub(crate) fn is_sha256(value: &str) -> bool {
    value.len() == 64 && value.bytes().all(|byte| byte.is_ascii_hexdigit())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dssp_numeric_tolerances_do_not_relax_other_features() {
        let expected = json!({ "phi_degrees": -75.0 });
        let actual = json!({ "phi_degrees": -74.9 });

        assert!(
            first_json_diff("bio.secondary-structure.dssp", "$", &expected, &actual,).is_none()
        );
        assert!(first_json_diff("unrelated.feature", "$", &expected, &actual).is_some());
    }

    #[test]
    fn molecular_descriptor_mass_tolerances_are_field_specific() {
        let expected = json!({
            "average_mass_da": 100.0,
            "monoisotopic_mass_da": 99.0,
            "formal_charge": 0,
        });
        let actual = json!({
            "average_mass_da": 100.049,
            "monoisotopic_mass_da": 99.000049,
            "formal_charge": 0,
        });
        assert!(first_json_diff("descriptor.molecular", "$", &expected, &actual).is_none());

        let wrong_charge = json!({
            "average_mass_da": 100.0,
            "monoisotopic_mass_da": 99.0,
            "formal_charge": 1,
        });
        assert!(first_json_diff("descriptor.molecular", "$", &expected, &wrong_charge).is_some());
        assert!(first_json_diff("unrelated.feature", "$", &expected, &actual).is_some());
    }
}
