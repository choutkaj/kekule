//! Explain decoded multiline-value differences without changing exact comparison.
//! CIF 1.1 paragraph 17 permits eliding line-end whitespace in text fields, but
//! requires preserving leading whitespace. Single-line values are excluded:
//! decoded observations cannot distinguish their original quoting syntax.
use serde_json::{json, Value};

fn trailing_whitespace_only(expected: &str, actual: &str) -> bool {
    expected != actual
        && expected.contains('\n')
        && expected
            .split('\n')
            .map(|line| line.trim_end_matches([' ', '\t']))
            .eq(actual
                .split('\n')
                .map(|line| line.trim_end_matches([' ', '\t'])))
}

pub(super) fn diagnostics(expected: &Value, actual: &Value) -> Vec<Value> {
    let (Some(expected), Some(actual)) = (
        expected.get("blocks").and_then(Value::as_array),
        actual.get("blocks").and_then(Value::as_array),
    ) else {
        return Vec::new();
    };
    if expected.len() != actual.len() {
        return Vec::new();
    }
    let mut result = Vec::new();
    for (block, (expected, actual)) in expected.iter().zip(actual).enumerate() {
        if expected.get("name").is_none() || expected.get("name") != actual.get("name") {
            continue;
        }
        let (Some(expected), Some(actual)) = (
            expected.get("values").and_then(Value::as_object),
            actual.get("values").and_then(Value::as_object),
        ) else {
            continue;
        };
        for (tag, expected) in expected {
            let (Some(expected), Some(actual)) = (
                expected.as_array(),
                actual.get(tag).and_then(Value::as_array),
            ) else {
                continue;
            };
            if expected.len() != actual.len() {
                continue;
            }
            for (index, (expected, actual)) in expected.iter().zip(actual).enumerate() {
                let (Some(expected), Some(actual)) = (expected.as_str(), actual.as_str()) else {
                    continue;
                };
                if trailing_whitespace_only(expected, actual) {
                    result.push(json!({"kind":"multiline_trailing_whitespace",
                        "path":format!("$.blocks[{block}].values.{tag}[{index}]")}));
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

    fn observation(text: &str) -> Value {
        json!({"blocks":[{"name":"sample","values":{"_details.text":[text]}}]})
    }

    #[test]
    fn whitespace_diagnostics_preserve_values_and_disagreement() {
        let expected = observation("  first\n\n\tlast");
        let actual = observation("  first \t\n \t\n\tlast  ");
        let original = (expected.clone(), actual.clone());
        for (expected, actual) in [(&expected, &actual), (&actual, &expected)] {
            let diff = differences("io.mmcif.parse", expected, actual);
            assert!(!diff.agrees());
            assert!(!diff.exact());
            assert_eq!(diff.structural, 1);
            assert_eq!(
                diff.details[0]["expected"],
                expected["blocks"][0]["values"]["_details.text"][0]
            );
            assert_eq!(
                diff.details[0]["actual"],
                actual["blocks"][0]["values"]["_details.text"][0]
            );
            assert_eq!(
                diff.diagnostics,
                vec![json!({
                    "kind":"multiline_trailing_whitespace",
                    "path":"$.blocks[0].values._details.text[0]"
                })]
            );
        }
        assert_eq!((expected, actual), original);
    }

    #[test]
    fn other_text_differences_and_single_line_values_are_not_explained_away() {
        for (expected, actual) in [
            ("first\nsecond", "first\n second"),
            ("first\nsecond", "first\nchanged"),
            ("first\nsecond", "first\nsecond\n"),
            ("first\nsecond", "first\n\nsecond"),
            ("first\nsecond", "first\r\nsecond"),
            ("first\nsecond", "first\nsecond\u{a0}"),
            ("first\nsecond", "first\nsecond\u{b}"),
            ("one line", "one line "),
            ("?", "."),
        ] {
            let diff = differences(
                "io.mmcif.parse",
                &observation(expected),
                &observation(actual),
            );
            assert!(!diff.agrees());
            assert!(diff.diagnostics.is_empty(), "{expected:?} / {actual:?}");
        }
        let value = observation("first\nsecond ");
        let diff = differences("io.mmcif.parse", &value, &value);
        assert!(diff.exact());
        assert!(diff.diagnostics.is_empty());
    }

    #[test]
    fn a_whitespace_diagnostic_does_not_cover_other_changed_values() {
        let mut expected = observation("first\nsecond");
        let mut actual = observation("first \nsecond");
        expected["blocks"][0]["values"]["_details.other"] = json!(["original"]);
        actual["blocks"][0]["values"]["_details.other"] = json!(["changed"]);
        let diff = differences("io.mmcif.parse", &expected, &actual);
        assert!(!diff.agrees());
        assert_eq!(diff.structural, 2);
        assert_eq!(diff.details.len(), 2);
        assert_eq!(diff.diagnostics.len(), 1);
        assert_eq!(
            diff.diagnostics[0]["path"],
            "$.blocks[0].values._details.text[0]"
        );
    }

    #[test]
    fn diagnostic_is_scoped_to_matching_blocks_and_value_positions() {
        let expected = observation("first\nsecond");
        let actual = observation("first \nsecond");
        let mut different_block = actual.clone();
        different_block["blocks"][0]["name"] = json!("other");
        let mut different_column = actual.clone();
        different_column["blocks"][0]["values"]["_details.text"] =
            json!(["first \nsecond", "extra"]);
        let mut missing_tag = actual.clone();
        missing_tag["blocks"][0]["values"] = json!({"_other.text":["first \nsecond"]});
        for other in [
            different_block,
            different_column,
            missing_tag,
            json!({"blocks":[]}),
        ] {
            assert!(diagnostics(&expected, &other).is_empty());
        }
        assert!(differences("another.feature", &expected, &actual)
            .diagnostics
            .is_empty());
    }
}
