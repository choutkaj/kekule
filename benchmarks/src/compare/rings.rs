//! Describe selected-cycle differences without changing the raw comparison.
//! These measurements concern reported paths, not their validity in the input
//! graph or whether either selection is a minimum or complete cycle basis.
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};

type Edge = (u64, u64);
type Cycle = BTreeSet<Edge>;

#[derive(Default)]
struct CycleSpan(BTreeMap<Edge, Cycle>);

impl CycleSpan {
    // Sparse Gaussian elimination over GF(2): adding paths toggles their edges.
    fn reduce(&self, mut cycle: Cycle) -> Cycle {
        while let Some(edge) = cycle.first() {
            let Some(pivot) = self.0.get(edge) else {
                break;
            };
            cycle = cycle.symmetric_difference(pivot).copied().collect();
        }
        cycle
    }

    fn insert(&mut self, cycle: Cycle) {
        let reduced = self.reduce(cycle);
        if let Some(&edge) = reduced.first() {
            self.0.insert(edge, reduced);
        }
    }

    fn same_span(&self, other: &Self) -> bool {
        self.0.len() == other.0.len()
            && other
                .0
                .values()
                .all(|cycle| self.reduce(cycle.clone()).is_empty())
    }
}

struct Selection {
    count: usize,
    edges: BTreeSet<Edge>,
    span: CycleSpan,
}

impl Selection {
    fn from_paths(value: &Value) -> Option<Self> {
        let paths = value.as_array()?;
        let mut selection = Self {
            count: paths.len(),
            edges: BTreeSet::new(),
            span: CycleSpan::default(),
        };
        for path in paths {
            let vertices: Vec<_> = path
                .as_array()?
                .iter()
                .map(Value::as_u64)
                .collect::<Option<_>>()?;
            if vertices.len() < 3
                || vertices.iter().collect::<BTreeSet<_>>().len() != vertices.len()
            {
                return None;
            }
            let cycle: Cycle = vertices
                .iter()
                .zip(vertices.iter().cycle().skip(1))
                .map(|(&a, &b)| (a.min(b), a.max(b)))
                .collect();
            selection.edges.extend(&cycle);
            selection.span.insert(cycle);
        }
        Some(selection)
    }

    fn summary(&self) -> Value {
        json!({"selected_cycles":self.count,"covered_edges":self.edges.len(),
            "cycle_span_rank":self.span.0.len()})
    }
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
    expected
        .iter()
        .zip(actual)
        .enumerate()
        .filter_map(|(index, (expected, actual))| {
            if expected.get("status")?.as_str()? != "ok"
                || actual.get("status")?.as_str()? != "ok"
                || expected.get("record_index")? != actual.get("record_index")?
                || expected.get("title")? != actual.get("title")?
                || expected.get("rings")? == actual.get("rings")?
            {
                return None;
            }
            let expected = Selection::from_paths(expected.get("rings")?)?;
            let actual = Selection::from_paths(actual.get("rings")?)?;
            Some(
                json!({"kind":"selected_cycle_span","path":format!("$.records[{index}].rings"),
                "expected":expected.summary(),"actual":actual.summary(),
                "same_covered_edges":expected.edges == actual.edges,
                "same_cycle_span":expected.span.same_span(&actual.span)}),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare::{differences, normalize_benchmark_for_comparison_in_place};

    fn observation(rings: Value) -> Value {
        json!({"records":[{"record_index":0,"status":"ok","title":"","rings":rings}]})
    }

    #[test]
    fn alternative_bases_and_redundant_cycles_remain_disagreements() {
        let expected = observation(json!([[0, 1, 2], [0, 2, 3]]));
        for paths in [
            json!([[0, 1, 2], [0, 1, 2, 3]]),
            json!([[0, 1, 2], [0, 2, 3], [0, 1, 2, 3]]),
        ] {
            let actual = observation(paths);
            let original = (expected.clone(), actual.clone());
            let diff = differences("algo.rings.sssr", &expected, &actual);
            assert!(!diff.agrees());
            assert!(!diff.exact());
            assert!(!diff.details.is_empty());
            assert_eq!(diff.diagnostics.len(), 1);
            let diagnostic = &diff.diagnostics[0];
            assert_eq!(diagnostic["same_covered_edges"], true);
            assert_eq!(diagnostic["same_cycle_span"], true);
            assert_eq!(diagnostic["expected"]["cycle_span_rank"], 2);
            assert_eq!(diagnostic["actual"]["cycle_span_rank"], 2);
            assert_eq!((&expected, &actual), (&original.0, &original.1));
        }
    }

    #[test]
    fn edge_coverage_does_not_imply_equal_cycle_span() {
        // Three independent triangles span K4's cycle space. Its three
        // Hamilton cycles cover the same edges but span only its even cycles.
        let expected = observation(json!([[0, 1, 2], [0, 1, 3], [0, 2, 3]]));
        let actual = observation(json!([[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3]]));
        let diff = differences("algo.rings.sssr", &expected, &actual);
        assert_eq!(diff.diagnostics[0]["same_covered_edges"], true);
        assert_eq!(diff.diagnostics[0]["same_cycle_span"], false);
        assert_eq!(diff.diagnostics[0]["expected"]["cycle_span_rank"], 3);
        assert_eq!(diff.diagnostics[0]["actual"]["cycle_span_rank"], 2);

        let missing = observation(json!([[0, 1, 2]]));
        let diff = differences("algo.rings.sssr", &expected, &missing);
        assert_eq!(diff.diagnostics[0]["same_covered_edges"], false);
        assert_eq!(diff.diagnostics[0]["same_cycle_span"], false);
    }

    #[test]
    fn equivalent_paths_need_no_diagnostic() {
        let mut expected = observation(json!([[0, 1, 2], [0, 2, 3]]));
        let mut actual = observation(json!([[3, 2, 0], [1, 2, 0]]));
        for value in [&mut expected, &mut actual] {
            normalize_benchmark_for_comparison_in_place("algo.rings.sssr", value);
        }
        let diff = differences("algo.rings.sssr", &expected, &actual);
        assert!(diff.exact());
        assert!(diff.diagnostics.is_empty());
        assert!(serde_json::to_value(diff)
            .unwrap()
            .get("diagnostics")
            .is_none());
    }

    #[test]
    fn diagnostics_require_matching_records_and_simple_paths() {
        let expected = observation(json!([[0, 1, 2]]));
        let actual = observation(json!([[0, 1, 2, 3]]));
        for (field, value) in [
            ("record_index", json!(1)),
            ("title", json!("different")),
            ("status", json!("error")),
            ("rings", json!([[0, 1]])),
            ("rings", json!([[0, 1, 2, 1]])),
            ("rings", json!([[0, 1, -2]])),
        ] {
            let mut invalid = actual.clone();
            invalid["records"][0][field] = value;
            assert!(diagnostics(&expected, &invalid).is_empty());
        }
        assert!(diagnostics(&expected, &json!({"records":[]})).is_empty());
        assert!(differences("another.feature", &expected, &actual)
            .diagnostics
            .is_empty());
    }
}
