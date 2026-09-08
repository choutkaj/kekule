//! Borrowed semantic category access for interpretation. The public document
//! retains scalar items and loops exactly as parsed; no coordinate table is copied.
use super::mmcif_interpret::MmcifInterpretError;
use super::{MmcifBlock, MmcifEntry, MmcifItem, MmcifLoopTable, MmcifValue};

pub(crate) enum MmcifCategory<'a> {
    Scalars(Vec<&'a MmcifItem>),
    Loop(&'a MmcifLoopTable),
}

pub(super) fn in_category(tag: &str, category: &str) -> bool {
    tag.split_once('.')
        .is_some_and(|(prefix, _)| prefix.eq_ignore_ascii_case(category))
}

impl MmcifBlock {
    /// Presence is independent of a particular required column. Malformed
    /// structural data must reach validation, not disappear during block selection.
    pub(crate) fn has_category(&self, category: &str) -> bool {
        self.entries().iter().any(|entry| match entry {
            MmcifEntry::Item(item) => in_category(item.tag(), category),
            MmcifEntry::Loop(table) => table.tags().iter().any(|tag| in_category(tag, category)),
        })
    }

    pub(crate) fn category(
        &self,
        category: &str,
    ) -> Result<Option<MmcifCategory<'_>>, MmcifInterpretError> {
        let mut scalars = Vec::new();
        let mut table = None;
        for entry in self.entries() {
            match entry {
                MmcifEntry::Item(item) if in_category(item.tag(), category) => scalars.push(item),
                MmcifEntry::Loop(candidate)
                    if candidate
                        .tags()
                        .iter()
                        .any(|tag| in_category(tag, category)) =>
                {
                    let line = candidate
                        .row(0)
                        .and_then(|row| row.first())
                        .map(MmcifValue::line);
                    if candidate
                        .tags()
                        .iter()
                        .any(|tag| !in_category(tag, category))
                    {
                        return Err(MmcifInterpretError::new(
                            line,
                            format!("category {category} shares a loop with another category"),
                        ));
                    }
                    if table.replace(candidate).is_some() {
                        return Err(MmcifInterpretError::new(
                            line,
                            format!("category {category} is split across multiple loops"),
                        ));
                    }
                }
                _ => {}
            }
        }
        match (table, scalars.is_empty()) {
            (Some(_), false) => Err(MmcifInterpretError::new(
                Some(scalars[0].value().line()),
                format!("category {category} mixes scalar items and a loop"),
            )),
            (Some(table), true) => Ok(Some(MmcifCategory::Loop(table))),
            (None, false) => Ok(Some(MmcifCategory::Scalars(scalars))),
            (None, true) => Ok(None),
        }
    }
}

impl<'a> MmcifCategory<'a> {
    pub(crate) fn row_count(&self) -> usize {
        match self {
            Self::Scalars(_) => 1,
            Self::Loop(table) => table.row_count(),
        }
    }

    pub(crate) fn value(&self, row: usize, tag: &str) -> Option<&'a MmcifValue> {
        match self {
            Self::Scalars(items) if row == 0 => items
                .iter()
                .find(|item| item.tag().eq_ignore_ascii_case(tag))
                .map(|item| item.value()),
            Self::Scalars(_) => None,
            Self::Loop(table) => table.value(row, tag),
        }
    }

    pub(crate) fn row_line(&self, row: usize) -> Option<usize> {
        match self {
            Self::Scalars(items) if row == 0 => items.first().map(|item| item.value().line()),
            Self::Scalars(_) => None,
            Self::Loop(table) => table
                .row(row)
                .and_then(|values| values.first())
                .map(MmcifValue::line),
        }
    }
}
