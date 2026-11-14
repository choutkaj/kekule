use super::*;

#[test]
fn shared_rows_append_and_adopt_tables() {
    let mut rows = ["a", "b"]
        .into_iter()
        .map(str::to_owned)
        .collect::<SharedRows<_>>();
    let copy = rows.share();
    assert_eq!(rows.push("c".to_owned()), 2);
    assert_eq!(copy.len(), 2);
    assert_eq!(rows.get(2), "c");
    let row = SharedRow::new("d".to_owned());
    let key = rows.key(&row);
    assert_eq!((key, rows.get(key).as_str()), (3, "d"));
    // A small foreign table's row is copied, not its table adopted.
    assert_eq!(rows.key(&row), 4);
    let large = (0..100).map(|i| i.to_string()).collect::<SharedRows<_>>();
    let adopted = rows.key(&large.row(7));
    assert_eq!((adopted, rows.get(adopted).as_str()), (12, "7"));
    assert_eq!(rows.key(&large.row(9)), 14);
    assert_eq!(&*rows.row(1), "b");
    let shared = copy.row(0);
    assert_eq!(rows.key(&shared), 0);
}

#[test]
fn compaction_merges_small_tables_and_rewrites_keys() {
    let mut rows = (0..2000).collect::<SharedRows<u32>>();
    let mut keys = vec![5, 99];
    for value in 1000..1020 {
        let table = (value..value + 64).collect::<SharedRows<u32>>();
        keys.push(rows.key(&table.row(0)));
    }
    // An unreferenced row of a small table is dropped by compaction.
    rows.key(&(7..71).collect::<SharedRows<u32>>().row(0));
    rows.compact(&mut keys);
    assert!(rows.tables.len() <= 2);
    let values = keys.iter().map(|&key| *rows.get(key)).collect::<Vec<_>>();
    let mut expected = vec![5, 99];
    expected.extend(1000..1020);
    assert_eq!(values, expected);
    assert_eq!(rows.len(), 2020);
}
