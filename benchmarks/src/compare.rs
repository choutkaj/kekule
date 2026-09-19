use serde::Serialize;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};

/// Raw differences remain in the report, including differences within the
/// documented reference precision. No values are rounded or rewritten.
#[derive(Default, Serialize)]
pub(crate) struct Differences {
    pub(crate) structural: usize,
    pub(crate) numerical: usize,
    pub(crate) within_precision: usize,
    pub(crate) details: Vec<Value>,
}

impl Differences {
    pub(crate) fn agrees(&self) -> bool {
        self.structural == 0 && self.numerical == self.within_precision
    }
    pub(crate) fn exact(&self) -> bool {
        self.structural == 0 && self.numerical == 0
    }
}

pub(crate) fn differences(feature: &str, expected: &Value, actual: &Value) -> Differences {
    fn visit(feature: &str, path: &str, a: &Value, b: &Value, result: &mut Differences) {
        match (a, b) {
            (Value::Object(a), Value::Object(b)) => {
                for key in a.keys().chain(b.keys()).collect::<BTreeSet<_>>() {
                    match (a.get(key), b.get(key)) {
                        (Some(a), Some(b)) => {
                            visit(feature, &format!("{path}.{key}"), a, b, result)
                        }
                        (a, b) => {
                            result.structural += 1;
                            result.details.push(json!({"path":format!("{path}.{key}"),"kind":"field_presence","expected":a,"actual":b}));
                        }
                    }
                }
            }
            (Value::Array(a), Value::Array(b)) if a.len() == b.len() => {
                for (index, (a, b)) in a.iter().zip(b).enumerate() {
                    visit(feature, &format!("{path}[{index}]"), a, b, result);
                }
            }
            _ if a == b => (),
            (Value::Number(a), Value::Number(b)) if a.is_f64() && b.is_f64() => {
                let (a, b) = (a.as_f64().unwrap(), b.as_f64().unwrap());
                let angle = path.ends_with("_degrees");
                let absolute_error = if angle {
                    ((a - b + 180.0).rem_euclid(360.0) - 180.0).abs()
                } else {
                    (a - b).abs()
                };
                let precision = reference_precision(feature, path, a, b);
                let accepted = absolute_error <= precision;
                result.numerical += 1;
                result.within_precision += usize::from(accepted);
                result.details.push(json!({"path":path,"kind":"number","expected":a,"actual":b,
                    "absolute_error":absolute_error,"allowed_error":precision,"within_precision":accepted}));
            }
            _ => {
                result.structural += 1;
                result
                    .details
                    .push(json!({"path":path,"kind":"structure","expected":a,"actual":b}));
            }
        }
    }
    let mut result = Differences::default();
    visit(feature, "$", expected, actual, &mut result);
    result
}

fn reference_precision(feature: &str, path: &str, a: f64, b: f64) -> f64 {
    // Sixteen machine epsilons cover the small fixed chain of unit conversions
    // and decimal parsing. This is not an empirical tolerance fit to a corpus.
    let roundoff = 16.0 * f64::EPSILON * a.abs().max(b.abs());
    let quantization = if feature == "bio.secondary-structure.dssp" {
        if path.ends_with("energy_kcal_per_mol")
            || [
                "phi_degrees",
                "psi_degrees",
                "alpha_degrees",
                "kappa_degrees",
            ]
            .iter()
            .any(|key| path.ends_with(key))
        {
            0.05 // mkdssp's documented legacy/mmCIF output: one decimal place.
        } else if path.ends_with("omega_degrees") {
            // DSSP evaluates angles in f32; Biopython evaluates omega in f64
            // from f32 coordinates. Allow arithmetic roundoff on the angular
            // domain, including near zero. This is not decimal quantization.
            16.0 * f64::from(f32::EPSILON) * 180.0
        } else if path.ends_with(".tco") {
            0.0005
        } else {
            0.0
        }
    } else if feature.contains(".v2000.write") && path.contains(".coord[") {
        0.00005 // V2000 coordinate fields carry four decimal places in angstroms.
    } else {
        0.0
    };
    quantization + roundoff
}

pub(crate) fn normalize_benchmark_for_comparison_in_place(benchmark_id: &str, value: &mut Value) {
    normalize_for_comparison_in_place(value);
    if benchmark_id == "bio.secondary-structure.dssp" {
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
            if let Some(Value::Array(descriptors)) = object.get_mut("bond_descriptors") {
                descriptors.sort_by_key(bond_sort_key);
            }
            normalize_ring_set_object(object);
        }
        _ => {}
    }
}

pub(crate) fn normalize_undirected_bond_object(object: &mut serde_json::Map<String, Value>) {
    let Some(begin) = object.get("begin_atom_index").and_then(Value::as_u64) else {
        return;
    };
    let Some(end) = object.get("end_atom_index").and_then(Value::as_u64) else {
        return;
    };
    if object.get("bond_type").and_then(Value::as_str) == Some("DATIVE") {
        return;
    }
    if begin > end {
        object.insert("begin_atom_index".to_owned(), json!(end));
        object.insert("end_atom_index".to_owned(), json!(begin));
    }
}

pub(crate) fn normalize_bond_array_object(object: &mut serde_json::Map<String, Value>) {
    let Some(Value::Array(bonds)) = object.get_mut("bonds") else {
        return;
    };
    bonds.sort_by_key(bond_sort_key);
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
        // A cycle's vertices are unique. Fix its starting vertex and direction
        // without losing the edges implied by consecutive atom indices.
        let Some(start) = atoms
            .iter()
            .enumerate()
            .min_by_key(|(_, value)| value.as_u64())
            .map(|(index, _)| index)
        else {
            continue;
        };
        atoms.rotate_left(start);
        if atoms
            .iter()
            .skip(1)
            .map(Value::as_u64)
            .cmp(atoms.iter().skip(1).rev().map(Value::as_u64))
            .is_gt()
        {
            atoms[1..].reverse();
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn dssp_omega_uses_single_precision_arithmetic_without_hiding_geometry_changes() {
        let expected = json!({"omega_degrees":179.0});
        let actual = json!({"omega_degrees":179.0 + f64::from(f32::EPSILON) * 180.0});
        let result = differences("bio.secondary-structure.dssp", &expected, &actual);
        assert!(result.agrees());
        assert!(!result.exact());
        assert_eq!(result.details[0]["actual"], actual["omega_degrees"]);
        let changed = json!({"omega_degrees":179.01});
        assert!(!differences("bio.secondary-structure.dssp", &expected, &changed).agrees());
        assert!(!differences("other.feature", &expected, &actual).agrees());
    }

    #[test]
    fn precision_never_hides_raw_values_or_structural_changes() {
        let expected =
            json!({"energy_kcal_per_mol":-2.8,"phi_degrees":179.99,"partner_sequence_id":12});
        let actual =
            json!({"energy_kcal_per_mol":-2.811,"phi_degrees":-179.99,"partner_sequence_id":12});
        let result = differences("bio.secondary-structure.dssp", &expected, &actual);
        assert!(result.agrees());
        assert!(!result.exact());
        assert_eq!(result.numerical, 2);
        assert_eq!(result.details.len(), 2);
        assert_eq!(result.details[0]["expected"], -2.8);
        assert_eq!(result.details[0]["actual"], -2.811);
        let mut changed = actual;
        changed["partner_sequence_id"] = json!(13);
        assert!(!differences("bio.secondary-structure.dssp", &expected, &changed).agrees());
        changed["energy_kcal_per_mol"] = json!(-2.86);
        assert_eq!(
            differences("bio.secondary-structure.dssp", &expected, &changed).within_precision,
            1
        );
    }
    #[test]
    fn roundoff_is_distinct_from_scientific_disagreement() {
        let a = json!({"coord":[1.0,0.0,0.0]});
        let b = json!({"coord":[1.0000000000000002,0.0,0.0]});
        assert!(differences("io.mol.parse", &a, &b).agrees());
        let c = json!({"coord":[1.00001,0.0,0.0]});
        assert!(!differences("io.mol.parse", &a, &c).agrees());
        assert!(differences("io.mol.v2000.write", &a, &c).agrees());
        assert!(!differences(
            "descriptor.molecular",
            &json!({"average_mass_da":12.011}),
            &json!({"average_mass_da":12.0107})
        )
        .agrees());
    }
    #[test]
    fn descriptor_order_follows_normalized_bond_endpoints() {
        let mut a = json!({"bond_descriptors":[{"begin_atom_index":0,"end_atom_index":3,"descriptor":"E"},{"begin_atom_index":1,"end_atom_index":2,"descriptor":"Z"}]});
        let mut b = json!({"bond_descriptors":[{"begin_atom_index":2,"end_atom_index":1,"descriptor":"Z"},{"begin_atom_index":3,"end_atom_index":0,"descriptor":"E"}]});
        normalize_for_comparison_in_place(&mut a);
        normalize_for_comparison_in_place(&mut b);
        assert_eq!(a, b);
        b["bond_descriptors"][0]["descriptor"] = json!("Z");
        assert!(!differences("", &a, &b).exact());
    }
    #[test]
    fn no_numeric_rounding_or_tolerance_hides_differences() {
        for _feature in [
            "io.sdf.parse",
            "descriptor.molecular",
            "bio.secondary-structure.dssp",
        ] {
            for key in [
                "coord",
                "average_mass_da",
                "phi_degrees",
                "energy_kcal_per_mol",
            ] {
                let expected = json!({key: 1.0});
                let actual = json!({key: 1.000000001});
                assert!(!differences(_feature, &expected, &actual).exact());
            }
        }
    }
}
