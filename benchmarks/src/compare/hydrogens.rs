//! Identify hydrogen representation differences only when all other graph
//! observations and every atom's current total hydrogen count agree.
use serde::Serialize;
use serde_json::{json, Value};

#[derive(PartialEq, Serialize)]
struct HydrogenState {
    declared: u64,
    inferred: u64,
    no_implicit_hydrogens: bool,
}

impl HydrogenState {
    fn read(atom: &Value) -> Option<Self> {
        let state = Self {
            declared: atom.get("explicit_hydrogens")?.as_u64()?,
            inferred: atom.get("implicit_hydrogens")?.as_u64()?,
            no_implicit_hydrogens: atom.get("no_implicit_hydrogens")?.as_bool()?,
        };
        if state.no_implicit_hydrogens && state.inferred != 0 {
            return None;
        }
        Some(state)
    }

    fn total(&self) -> Option<u64> {
        self.declared.checked_add(self.inferred)
    }

    fn bond_valence(&self, atom: &Value) -> Option<u64> {
        // Both observations include declared non-graph H in explicit valence.
        atom.get("explicit_valence")?
            .as_u64()?
            .checked_sub(self.declared)
    }
}

fn equal_except(expected: &Value, actual: &Value, fields: &[&str]) -> bool {
    let (Some(expected), Some(actual)) = (expected.as_object(), actual.as_object()) else {
        return false;
    };
    expected
        .iter()
        .filter(|(key, _)| !fields.contains(&key.as_str()))
        .eq(actual
            .iter()
            .filter(|(key, _)| !fields.contains(&key.as_str())))
}

fn graph_diagnostic(expected: &Value, actual: &Value, path: &str) -> Option<Value> {
    if expected == actual || !equal_except(expected, actual, &["atoms"]) {
        return None;
    }
    let (expected, actual) = (
        expected.get("atoms")?.as_array()?,
        actual.get("atoms")?.as_array()?,
    );
    if expected.len() != actual.len() {
        return None;
    }
    let mut changes = Vec::new();
    for (expected, actual) in expected.iter().zip(actual) {
        if !equal_except(
            expected,
            actual,
            &[
                "explicit_hydrogens",
                "implicit_hydrogens",
                "no_implicit_hydrogens",
                "explicit_valence",
            ],
        ) {
            return None;
        }
        let e = HydrogenState::read(expected)?;
        let a = HydrogenState::read(actual)?;
        let total = e.total()?;
        if total != a.total()? || e.bond_valence(expected)? != a.bond_valence(actual)? {
            return None;
        }
        if e != a {
            changes.push(json!({"atom_index":expected.get("index")?.as_u64()?,
                "non_graph_hydrogens":total,"expected":e,"actual":a}));
        }
    }
    if changes.is_empty() {
        return None;
    }
    Some(json!({"kind":"hydrogen_representation_only","path":path,"atoms":changes}))
}

pub(super) fn diagnostics(expected: &Value, actual: &Value) -> Vec<Value> {
    let (Some(expected), Some(actual)) = (
        expected.get("records").and_then(Value::as_array),
        actual.get("records").and_then(Value::as_array),
    ) else {
        return Vec::new();
    };
    if expected.len() != actual.len() {
        return Vec::new();
    }
    let mut result = Vec::new();
    for (index, (expected, actual)) in expected.iter().zip(actual).enumerate() {
        if expected.get("status").and_then(Value::as_str) != Some("ok")
            || actual.get("status").and_then(Value::as_str) != Some("ok")
            || expected.get("record_index").is_none()
            || expected.get("record_index") != actual.get("record_index")
            || expected.get("title").is_none()
            || expected.get("title") != actual.get("title")
        {
            continue;
        }
        for stage in ["added_graph", "round_trip"] {
            if let (Some(expected), Some(actual)) = (expected.get(stage), actual.get(stage)) {
                if let Some(diagnostic) =
                    graph_diagnostic(expected, actual, &format!("$.records[{index}].{stage}"))
                {
                    result.push(diagnostic);
                }
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare::differences;

    fn atom(declared: u64, inferred: u64, fixed: bool) -> Value {
        json!({"index":0,"atomic_number":7,"symbol":"N","isotope":null,
            "formal_charge":0,"radical_electrons":0,"spin_multiplicity":null,
            "aromatic":true,"atom_map":null,"coord":null,
            "explicit_hydrogens":declared,"implicit_hydrogens":inferred,
            "no_implicit_hydrogens":fixed,"explicit_valence":2+declared})
    }

    fn observation(atom: Value) -> Value {
        json!({"records":[{"record_index":0,"status":"ok","title":"",
            "added_graph":{"atoms":[atom],"bonds":[],"stereo":[],"groups":[]},
            "round_trip":{"atoms":[atom],"bonds":[],"stereo":[],"groups":[]}}]})
    }

    #[test]
    fn declared_and_inferred_hydrogens_are_diagnosed_without_accepting_them() {
        let expected = observation(atom(1, 0, true));
        let actual = observation(atom(0, 1, false));
        let originals = (expected.clone(), actual.clone());
        let diff = differences("chem.hydrogen-transforms", &expected, &actual);
        assert!(!diff.agrees());
        assert!(!diff.exact());
        assert_eq!(diff.structural, 8);
        assert_eq!(diff.details.len(), 8);
        assert_eq!(diff.diagnostics.len(), 2);
        assert_eq!(diff.diagnostics[1]["path"], "$.records[0].round_trip");
        assert_eq!(diff.diagnostics[1]["atoms"][0]["non_graph_hydrogens"], 1);
        assert_eq!(diff.diagnostics[1]["atoms"][0]["expected"]["declared"], 1);
        assert_eq!(diff.diagnostics[1]["atoms"][0]["actual"]["inferred"], 1);
        assert_eq!((expected, actual), originals);
    }

    #[test]
    fn inference_policy_with_zero_hydrogens_is_still_a_raw_difference() {
        let expected = observation(atom(0, 0, false));
        let actual = observation(atom(0, 0, true));
        let diff = differences("chem.hydrogen-transforms", &expected, &actual);
        assert!(!diff.agrees());
        assert_eq!(diff.diagnostics.len(), 2);
        assert_eq!(diff.diagnostics[0]["atoms"][0]["non_graph_hydrogens"], 0);
    }

    #[test]
    fn other_atom_and_graph_changes_prevent_a_storage_only_diagnostic() {
        let expected = observation(atom(1, 0, true));
        let actual = observation(atom(0, 1, false));
        for (key, value) in [
            ("implicit_hydrogens", json!(0)),
            ("formal_charge", json!(1)),
            ("radical_electrons", json!(1)),
            ("spin_multiplicity", json!(2)),
            ("isotope", json!(15)),
            ("index", json!(1)),
            ("explicit_valence", json!(3)),
        ] {
            let mut changed = actual["records"][0]["round_trip"].clone();
            changed["atoms"][0][key] = value;
            assert!(
                graph_diagnostic(&expected["records"][0]["round_trip"], &changed, "test").is_none()
            );
        }
        for key in ["bonds", "stereo", "groups", "atoms"] {
            let mut changed = actual["records"][0]["round_trip"].clone();
            changed[key]
                .as_array_mut()
                .unwrap()
                .push(json!({"different":true}));
            assert!(
                graph_diagnostic(&expected["records"][0]["round_trip"], &changed, "test").is_none()
            );
        }
    }

    #[test]
    fn invalid_counts_and_record_correspondence_never_prove_equivalence() {
        let expected = observation(atom(1, 0, true));
        let mut actual = observation(atom(0, 1, false));
        actual["records"][0]["title"] = json!("other");
        assert!(diagnostics(&expected, &actual).is_empty());
        assert!(diagnostics(&expected, &json!({"records":[]})).is_empty());
        assert!(diagnostics(&expected, &observation(atom(0, 1, true))).is_empty());
        let mut invalid = atom(1, 0, false);
        invalid["explicit_valence"] = json!(0);
        assert!(diagnostics(&expected, &observation(invalid)).is_empty());
        let mut overflow = atom(0, 0, false);
        overflow["explicit_hydrogens"] = json!(u64::MAX);
        overflow["implicit_hydrogens"] = json!(1);
        assert!(diagnostics(&expected, &observation(overflow)).is_empty());
        assert!(
            differences("another.feature", &expected, &observation(atom(0, 1, true)))
                .diagnostics
                .is_empty()
        );
    }
}
