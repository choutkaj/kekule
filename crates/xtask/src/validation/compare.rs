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

pub(crate) fn normalize_for_comparison(value: &Value) -> Value {
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
            normalize_coord_object(&mut normalized);
            Value::Object(normalized)
        }
        _ => value.clone(),
    }
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
}
