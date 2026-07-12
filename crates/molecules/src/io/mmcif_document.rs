use std::collections::BTreeSet;

use super::{tokenize_mmcif, MmcifParseError, MmcifParseOptions, MmcifToken};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifDocument {
    blocks: Vec<MmcifDataBlock>,
}

impl MmcifDocument {
    pub fn blocks(&self) -> &[MmcifDataBlock] {
        &self.blocks
    }

    pub fn block(&self, name: &str) -> Option<&MmcifDataBlock> {
        self.blocks
            .iter()
            .find(|block| block.name.eq_ignore_ascii_case(name))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifDataBlock {
    name: String,
    entries: Vec<MmcifEntry>,
}

impl MmcifDataBlock {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn entries(&self) -> &[MmcifEntry] {
        &self.entries
    }

    pub fn item(&self, tag: &str) -> Option<&MmcifValue> {
        self.entries.iter().find_map(|entry| match entry {
            MmcifEntry::Item(item) if item.tag.eq_ignore_ascii_case(tag) => Some(&item.value),
            _ => None,
        })
    }

    pub fn loop_with_tag(&self, tag: &str) -> Option<&MmcifLoopTable> {
        self.entries.iter().find_map(|entry| match entry {
            MmcifEntry::Loop(table) if table.column_index(tag).is_some() => Some(table),
            _ => None,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifEntry {
    Item(MmcifItem),
    Loop(MmcifLoopTable),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifItem {
    tag: String,
    value: MmcifValue,
}

impl MmcifItem {
    pub fn tag(&self) -> &str {
        &self.tag
    }

    pub fn value(&self) -> &MmcifValue {
        &self.value
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifLoopTable {
    tags: Vec<String>,
    values: Vec<MmcifValue>,
}

impl MmcifLoopTable {
    pub fn tags(&self) -> &[String] {
        &self.tags
    }

    pub fn row_count(&self) -> usize {
        self.values.len() / self.tags.len()
    }

    pub fn column_index(&self, tag: &str) -> Option<usize> {
        self.tags
            .iter()
            .position(|candidate| candidate.eq_ignore_ascii_case(tag))
    }

    pub fn value(&self, row: usize, tag: &str) -> Option<&MmcifValue> {
        let column = self.column_index(tag)?;
        self.values.get(row.checked_mul(self.tags.len())? + column)
    }

    pub fn row(&self, row: usize) -> Option<&[MmcifValue]> {
        let start = row.checked_mul(self.tags.len())?;
        self.values.get(start..start + self.tags.len())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifValue {
    text: String,
    line: usize,
    bare: bool,
}

impl MmcifValue {
    pub fn text(&self) -> &str {
        &self.text
    }

    pub const fn line(&self) -> usize {
        self.line
    }

    pub fn is_missing(&self) -> bool {
        self.bare && matches!(self.text.as_str(), "." | "?")
    }

    pub fn optional_text(&self) -> Option<&str> {
        (!self.is_missing()).then_some(self.text.as_str())
    }
}

pub fn parse_mmcif_str(
    input: &str,
    options: MmcifParseOptions,
) -> Result<MmcifDocument, MmcifParseError> {
    let tokens = tokenize_mmcif(input, options)?;
    parse_tokens(&tokens, options.max_atom_site_rows)
}

fn parse_tokens(
    tokens: &[MmcifToken],
    max_atom_site_rows: usize,
) -> Result<MmcifDocument, MmcifParseError> {
    let mut blocks = Vec::new();
    let mut block_names = BTreeSet::new();
    let mut index = 0usize;
    while index < tokens.len() {
        let token = &tokens[index];
        if !has_bare_prefix(token, "data_") {
            return Err(MmcifParseError::new(
                token.line,
                "content appears before the first data block",
            ));
        }
        let name = token.text[5..].to_owned();
        if name.is_empty() {
            return Err(MmcifParseError::new(token.line, "data block has no name"));
        }
        if !block_names.insert(name.to_ascii_lowercase()) {
            return Err(MmcifParseError::new(
                token.line,
                format!("duplicate data block {name}"),
            ));
        }
        index += 1;
        let mut entries = Vec::new();
        let mut data_names = BTreeSet::new();
        while index < tokens.len() && !has_bare_prefix(&tokens[index], "data_") {
            if is_bare_word(&tokens[index], "loop_") {
                let (table, next) = parse_loop(tokens, index)?;
                if table
                    .tags()
                    .iter()
                    .any(|tag| tag.to_ascii_lowercase().starts_with("_atom_site."))
                    && table.row_count() > max_atom_site_rows
                {
                    return Err(MmcifParseError::new(
                        tokens[index].line,
                        "atom-site row count exceeds configured limit",
                    ));
                }
                for tag in table.tags() {
                    if !data_names.insert(tag.to_ascii_lowercase()) {
                        return Err(MmcifParseError::new(
                            tokens[index].line,
                            format!("duplicate data name {tag}"),
                        ));
                    }
                }
                entries.push(MmcifEntry::Loop(table));
                index = next;
            } else if tokens[index].bare && tokens[index].text.starts_with('_') {
                let tag = tokens[index].text.clone();
                let line = tokens[index].line;
                if !data_names.insert(tag.to_ascii_lowercase()) {
                    return Err(MmcifParseError::new(
                        line,
                        format!("duplicate data name {tag}"),
                    ));
                }
                index += 1;
                let value = tokens.get(index).ok_or_else(|| {
                    MmcifParseError::new(line, format!("missing value for {tag}"))
                })?;
                if is_control(value) {
                    return Err(MmcifParseError::new(
                        value.line,
                        format!("missing value for {tag}"),
                    ));
                }
                entries.push(MmcifEntry::Item(MmcifItem {
                    tag,
                    value: owned_value(value),
                }));
                index += 1;
            } else {
                return Err(MmcifParseError::new(
                    tokens[index].line,
                    format!("unexpected token `{}` in data block", tokens[index].text),
                ));
            }
        }
        blocks.push(MmcifDataBlock { name, entries });
    }
    if blocks.is_empty() {
        return Err(MmcifParseError::new(1, "missing data block"));
    }
    Ok(MmcifDocument { blocks })
}

fn parse_loop(
    tokens: &[MmcifToken],
    start: usize,
) -> Result<(MmcifLoopTable, usize), MmcifParseError> {
    let line = tokens[start].line;
    let mut index = start + 1;
    let mut tags = Vec::new();
    let mut seen = BTreeSet::new();
    while index < tokens.len() && tokens[index].bare && tokens[index].text.starts_with('_') {
        if !seen.insert(tokens[index].text.to_ascii_lowercase()) {
            return Err(MmcifParseError::new(
                tokens[index].line,
                format!("duplicate loop tag {}", tokens[index].text),
            ));
        }
        tags.push(tokens[index].text.clone());
        index += 1;
    }
    if tags.is_empty() {
        return Err(MmcifParseError::new(line, "loop has no tags"));
    }
    let mut values = Vec::new();
    while index < tokens.len() && !is_control(&tokens[index]) {
        values.push(owned_value(&tokens[index]));
        index += 1;
    }
    if values.len() % tags.len() != 0 {
        let line = values.last().map(MmcifValue::line).unwrap_or(line);
        return Err(MmcifParseError::new(line, "loop has ragged rows"));
    }
    if index < tokens.len() && is_bare_word(&tokens[index], "stop_") {
        index += 1;
    }
    Ok((MmcifLoopTable { tags, values }, index))
}

fn is_control(token: &MmcifToken) -> bool {
    token.bare
        && (token.text.eq_ignore_ascii_case("loop_")
            || token.text.eq_ignore_ascii_case("stop_")
            || token.text.eq_ignore_ascii_case("global_")
            || token
                .text
                .get(..5)
                .is_some_and(|prefix| prefix.eq_ignore_ascii_case("data_"))
            || token.text.starts_with('_')
            || token
                .text
                .get(..5)
                .is_some_and(|prefix| prefix.eq_ignore_ascii_case("save_")))
}

fn is_bare_word(token: &MmcifToken, word: &str) -> bool {
    token.bare && token.text.eq_ignore_ascii_case(word)
}

fn has_bare_prefix(token: &MmcifToken, prefix: &str) -> bool {
    token.bare
        && token
            .text
            .get(..prefix.len())
            .is_some_and(|value| value.eq_ignore_ascii_case(prefix))
}

fn owned_value(token: &MmcifToken) -> MmcifValue {
    MmcifValue {
        text: token.text.clone(),
        line: token.line,
        bare: token.bare,
    }
}
