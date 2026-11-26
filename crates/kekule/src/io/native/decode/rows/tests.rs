use std::collections::HashSet;

use super::{distinct_rows, StringColumn};
use crate::io::native::encode::ABSENT_STRING;

fn column(entries: &[Option<String>]) -> StringColumn {
    StringColumn {
        text: entries.iter().flatten().map(String::as_str).collect(),
        lengths: entries
            .iter()
            .map(|entry| {
                entry
                    .as_ref()
                    .map_or(ABSENT_STRING, |entry| entry.len() as u32)
            })
            .collect(),
    }
}

#[test]
fn every_node_gets_the_row_of_its_own_entries() {
    // Thousands of distinct rows, more than the table of recent rows has
    // slots, met again out of order; absent and empty entries differ.
    let count = 30_000;
    let names = (0..count)
        .map(|node| (node % 7 != 0).then(|| format!("é{}", (node * 7919) % 5_000)))
        .collect::<Vec<_>>();
    let flags = (0..count)
        .map(|node| (node % 3 != 0).then(String::new))
        .collect::<Vec<_>>();
    let extras = (0..count)
        .map(|node| (node % 5 == 0).then_some((node % 11) as i32))
        .collect::<Vec<_>>();
    let columns = [column(&names), column(&flags)];
    let (keys, rows) = distinct_rows(
        &columns,
        extras.iter().copied(),
        count,
        |[name, flag], extra| (name.map(str::to_owned), flag.map(str::to_owned), extra),
    );
    assert_eq!(keys.len(), count);
    for (node, &key) in keys.iter().enumerate() {
        assert_eq!(
            rows[key as usize],
            (names[node].clone(), flags[node].clone(), extras[node]),
            "node {node}"
        );
    }
    assert_eq!(rows.iter().collect::<HashSet<_>>().len(), rows.len());
    assert!(rows.len() > 4_096 && rows.len() < count);
    // Keys number the rows in order of first appearance.
    let mut next = 0;
    for &key in &keys {
        assert!(key <= next);
        if key == next {
            next += 1;
        }
    }
    assert_eq!(next as usize, rows.len());
}
