//! Row tables shared between topology snapshots and drafts.
//!
//! Hierarchy labels are large and rarely change. A hierarchy stores each
//! residue's labels and each atom site's metadata as a key into
//! [`SharedRows`]: immutable tables shared by every snapshot and draft that
//! keeps them. Publishing an edit copies keys, which are plain integers, and
//! adds only the rows the edit changed, as one small new table.

use std::fmt;
use std::ops::Deref;
use std::sync::Arc;

/// One row of a shared table: the table and the row's position in it.
pub(crate) struct SharedRow<T> {
    table: Arc<Vec<T>>,
    index: u32,
}

impl<T> SharedRow<T> {
    /// A row of its own one-row table.
    pub(crate) fn new(value: T) -> Self {
        Self {
            table: Arc::new(vec![value]),
            index: 0,
        }
    }
}

impl<T> Clone for SharedRow<T> {
    fn clone(&self) -> Self {
        Self {
            table: Arc::clone(&self.table),
            index: self.index,
        }
    }
}

impl<T> Deref for SharedRow<T> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.table[self.index as usize]
    }
}

impl<T: fmt::Debug> fmt::Debug for SharedRow<T> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        T::fmt(self, formatter)
    }
}

impl<T: PartialEq> PartialEq for SharedRow<T> {
    fn eq(&self, other: &Self) -> bool {
        **self == **other
    }
}

/// Rows addressed by key: the concatenation of shared tables.
#[derive(Clone)]
pub(crate) struct SharedRows<T> {
    tables: Vec<Arc<Vec<T>>>,
    /// The key of each table's first row, then the row count.
    starts: Vec<u32>,
}

impl<T> Default for SharedRows<T> {
    fn default() -> Self {
        Self {
            tables: Vec::new(),
            starts: vec![0],
        }
    }
}

impl<T: fmt::Debug> fmt::Debug for SharedRows<T> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.debug_list().entries(self.rows()).finish()
    }
}

/// More tables than this are merged, keeping lookups short.
const MAX_TABLES: usize = 16;

/// Rows of smaller foreign tables are copied rather than their tables
/// adopted, so one-row tables of edited nodes do not accumulate.
const ADOPTED_ROWS: usize = 64;

impl<T> SharedRows<T> {
    /// The number of rows.
    pub(crate) fn len(&self) -> usize {
        *self.starts.last().expect("starts end with the row count") as usize
    }

    fn locate(&self, key: u32) -> (usize, usize) {
        // The last table holding a start at or below the key holds it.
        let table = self.starts.partition_point(|&start| start <= key) - 1;
        (table, (key - self.starts[table]) as usize)
    }

    /// The row with `key`.
    pub(crate) fn get(&self, key: u32) -> &T {
        let (table, index) = self.locate(key);
        &self.tables[table][index]
    }

    /// The row with `key` as a shared row.
    pub(crate) fn row(&self, key: u32) -> SharedRow<T> {
        let (table, index) = self.locate(key);
        SharedRow {
            table: Arc::clone(&self.tables[table]),
            index: index as u32,
        }
    }

    /// Every row in key order.
    fn rows(&self) -> impl Iterator<Item = &T> {
        self.tables.iter().flat_map(|table| table.iter())
    }

    /// Appends one row and returns its key, growing the last table when no
    /// other snapshot shares it.
    pub(crate) fn push(&mut self, value: T) -> u32 {
        let key = self.len() as u32;
        if let Some(table) = self.tables.last_mut().and_then(Arc::get_mut) {
            table.push(value);
        } else {
            self.tables.push(Arc::new(Vec::new()));
            self.starts.push(key);
            Arc::get_mut(self.tables.last_mut().expect("just pushed"))
                .expect("a new table is unshared")
                .push(value);
        }
        *self
            .starts
            .last_mut()
            .expect("starts end with the row count") = key + 1;
        key
    }

    /// Rows with the same tables, for keys copied from `self`.
    pub(crate) fn share(&self) -> Self {
        Self {
            tables: self.tables.clone(),
            starts: self.starts.clone(),
        }
    }
}

impl<T: Clone> SharedRows<T> {
    /// The key of a shared row: its key here when its table is one of these
    /// tables, else the key of a copy of a row from a small table, or of the
    /// row after adopting its large table.
    pub(crate) fn key(&mut self, row: &SharedRow<T>) -> u32 {
        // Rows mostly come from the first table, a published source's, or
        // the last one.
        let found = match self.tables.len() {
            0 => None,
            count => [0, count - 1]
                .into_iter()
                .chain(1..count - 1)
                .find(|&table| Arc::ptr_eq(&self.tables[table], &row.table)),
        };
        if let Some(table) = found {
            return self.starts[table] + row.index;
        }
        if row.table.len() < ADOPTED_ROWS {
            return self.push((**row).clone());
        }
        let start = self.len() as u32;
        self.tables.push(Arc::clone(&row.table));
        self.starts
            .push(start + u32::try_from(row.table.len()).expect("rows fit keys"));
        start + row.index
    }

    /// Merges the small tables into one when there are many, keeping the
    /// rows `keys` name and rewriting `keys` to match.
    pub(crate) fn compact(&mut self, keys: &mut [u32]) {
        if self.tables.len() <= MAX_TABLES {
            return;
        }
        let largest = self
            .tables
            .iter()
            .map(|table| table.len())
            .max()
            .unwrap_or(0);
        let small = |table: &Arc<Vec<T>>| table.len() * MAX_TABLES < largest.max(1);
        let mut moved = Vec::with_capacity(self.tables.len());
        let (mut tables, mut starts) = (Vec::new(), vec![0u32]);
        for table in &self.tables {
            if small(table) {
                moved.push(None);
            } else {
                let start = *starts.last().expect("starts end with the row count");
                moved.push(Some(start));
                tables.push(Arc::clone(table));
                starts.push(start + table.len() as u32);
            }
        }
        // Each kept row of a small table moves into one merged table.
        let mut merged = Vec::new();
        let merged_start = *starts.last().expect("nonempty");
        let mut remap = std::collections::HashMap::new();
        for key in keys.iter_mut() {
            let (table, index) = self.locate(*key);
            *key = match moved[table] {
                Some(start) => start + index as u32,
                None => *remap.entry(*key).or_insert_with(|| {
                    merged.push(self.tables[table][index].clone());
                    merged_start + merged.len() as u32 - 1
                }),
            };
        }
        if !merged.is_empty() {
            let end = merged_start + merged.len() as u32;
            tables.push(Arc::new(merged));
            starts.push(end);
        }
        self.tables = tables;
        self.starts = starts;
    }
}

impl<T> FromIterator<T> for SharedRows<T> {
    fn from_iter<I: IntoIterator<Item = T>>(rows: I) -> Self {
        let table = rows.into_iter().collect::<Vec<_>>();
        if table.is_empty() {
            return Self::default();
        }
        let len = u32::try_from(table.len()).expect("rows fit keys");
        Self {
            tables: vec![Arc::new(table)],
            starts: vec![0, len],
        }
    }
}

#[cfg(test)]
mod tests {
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
}
