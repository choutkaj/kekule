use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

use crate::diagnostics::{normalize_line_endings, read_bounded, DiagnosticKind, SourceLocation};

use super::structure_interpret::{
    interpret_mmcif, interpret_mmcif_block, interpret_mmcif_conformations,
    interpret_mmcif_conformations_block, interpret_mmcif_ensemble, interpret_mmcif_ensemble_block,
    StructureEnsembleInterpretError, StructureEnsembleInterpretOptions,
    StructureEnsembleInterpretation, StructureInterpretError, StructureInterpretOptions,
    StructureInterpretation,
};

/// Number of `_atom_site` columns in a typical wwPDB entry, used to keep the
/// default token budget consistent with the default atom-site row limit.
const TYPICAL_ATOM_SITE_COLUMNS: usize = 21;

/// Syntax and resource limits of CIF/mmCIF parsing.
///
/// Parsed values are stored as compact spans into one shared copy of the
/// source, about 16 bytes per tag or value. The defaults are mutually
/// consistent: `max_tokens` admits `max_atom_site_rows` rows of a typical
/// 21-column wwPDB atom-site table plus the other categories of the entry, and
/// `max_input_bytes` admits the largest current wwPDB entries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MmcifParseOptions {
    /// Maximum UTF-8 byte length of the input; at most 4 GiB.
    pub max_input_bytes: usize,
    /// Maximum number of tags, values, and control words.
    pub max_tokens: usize,
    /// Maximum byte length of one value.
    pub max_token_bytes: usize,
    /// Maximum number of `_atom_site` rows in one data block.
    pub max_atom_site_rows: usize,
}

impl Default for MmcifParseOptions {
    fn default() -> Self {
        let max_atom_site_rows = 5_000_000;
        Self {
            max_input_bytes: 1024 * 1024 * 1024,
            max_tokens: max_atom_site_rows * TYPICAL_ATOM_SITE_COLUMNS + 20_000_000,
            max_token_bytes: 16 * 1024 * 1024,
            max_atom_site_rows,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifParseError {
    kind: DiagnosticKind,
    location: SourceLocation,
    message: String,
}

impl MmcifParseError {
    fn new(kind: DiagnosticKind, line: usize, message: impl Into<String>) -> Self {
        Self {
            kind,
            location: SourceLocation::at_line(line),
            message: message.into(),
        }
    }

    pub const fn kind(&self) -> DiagnosticKind {
        self.kind
    }

    pub const fn location(&self) -> SourceLocation {
        self.location
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for MmcifParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "mmCIF parse error at {}: {}",
            self.location, self.message
        )
    }
}

impl std::error::Error for MmcifParseError {}

/// One parsed tag or value: a byte span of the shared source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Span {
    start: u32,
    end: u32,
    line: u32,
    bare: bool,
}

impl Span {
    fn text(self, source: &str) -> &str {
        &source[self.start as usize..self.end as usize]
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifDocument {
    blocks: Vec<MmcifBlock>,
}

impl MmcifDocument {
    pub fn blocks(&self) -> &[MmcifBlock] {
        &self.blocks
    }

    /// Interprets caller-supplied selections as an unweighted ensemble, in the
    /// supplied order. Repeated selections are retained. Exactly one atom-site
    /// block is required; atom identities, chemistry and topology layout must
    /// agree between selections. No combinations or statistical weights are inferred.
    pub fn interpret_conformations(
        &self,
        selections: &[StructureInterpretOptions],
    ) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
        interpret_mmcif_conformations(self, selections)
    }

    pub fn block(&self, name: &str) -> Option<&MmcifBlock> {
        self.blocks
            .iter()
            .find(|block| block.name.eq_ignore_ascii_case(name))
    }

    /// Interprets the exactly one structural block as one selected coordinate
    /// model under the model, alternate-location, and connectivity policy of
    /// `options`.
    pub fn interpret(
        &self,
        options: StructureInterpretOptions,
    ) -> Result<StructureInterpretation, StructureInterpretError> {
        interpret_mmcif(self, options)
    }

    /// Interprets the coordinate models that `options` selects, by default
    /// all of them, in the exactly one structural block as an ensemble of
    /// equally weighted members.
    pub fn interpret_ensemble(
        &self,
        options: StructureEnsembleInterpretOptions,
    ) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
        interpret_mmcif_ensemble(self, options)
    }
}

/// One independently interpretable CIF/mmCIF data block.
///
/// The block shares one copy of the source text with its document; items and
/// loops are read through borrowed views.
#[derive(Clone, PartialEq, Eq)]
pub struct MmcifBlock {
    name: String,
    source: Arc<str>,
    entries: Vec<EntryData>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum EntryData {
    Item(ItemData),
    Loop(LoopData),
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ItemData {
    tag: String,
    value: Span,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LoopData {
    tags: Vec<String>,
    values: Vec<Span>,
}

impl fmt::Debug for MmcifBlock {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MmcifBlock")
            .field("name", &self.name)
            .field("entries", &self.entries().collect::<Vec<_>>())
            .finish()
    }
}

impl MmcifBlock {
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Interprets explicit selections from this block as an unweighted ensemble.
    /// See [`MmcifDocument::interpret_conformations`] for compatibility requirements.
    pub fn interpret_conformations(
        &self,
        selections: &[StructureInterpretOptions],
    ) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
        interpret_mmcif_conformations_block(self, selections)
    }

    /// Scalar items and loops in source order.
    pub fn entries(&self) -> impl ExactSizeIterator<Item = MmcifEntry<'_>> {
        self.entries.iter().map(|entry| match entry {
            EntryData::Item(data) => MmcifEntry::Item(MmcifItem {
                source: &self.source,
                data,
            }),
            EntryData::Loop(data) => MmcifEntry::Loop(MmcifLoopTable {
                source: &self.source,
                data,
            }),
        })
    }

    /// Inventories alternate labels and stable residue identities without
    /// selecting coordinates or changing the source document. Identities can be
    /// used in [`super::AltLocSelection`] overrides and residue groups.
    pub fn alternate_locations(
        &self,
    ) -> Result<Vec<super::AltLocResidue>, StructureInterpretError> {
        super::structure_interpret::alternate_location_inventory(self)
    }

    /// The value of one scalar item, matching the tag case-insensitively.
    pub fn item(&self, tag: &str) -> Option<MmcifValue<'_>> {
        self.entries().find_map(|entry| match entry {
            MmcifEntry::Item(item) if item.tag().eq_ignore_ascii_case(tag) => Some(item.value()),
            _ => None,
        })
    }

    /// The loop containing `tag`, matching it case-insensitively.
    pub fn loop_with_tag(&self, tag: &str) -> Option<MmcifLoopTable<'_>> {
        self.entries().find_map(|entry| match entry {
            MmcifEntry::Loop(table) if table.column_index(tag).is_some() => Some(table),
            _ => None,
        })
    }

    /// Interprets one selected coordinate model of this block under the
    /// policy of `options`.
    pub fn interpret(
        &self,
        options: StructureInterpretOptions,
    ) -> Result<StructureInterpretation, StructureInterpretError> {
        interpret_mmcif_block(self, options)
    }

    /// Interprets the coordinate models that `options` selects, by default
    /// all of them, in this block as an ensemble of equally weighted members.
    pub fn interpret_ensemble(
        &self,
        options: StructureEnsembleInterpretOptions,
    ) -> Result<StructureEnsembleInterpretation, StructureEnsembleInterpretError> {
        interpret_mmcif_ensemble_block(self, options)
    }
}

/// One scalar item or loop of a data block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MmcifEntry<'a> {
    Item(MmcifItem<'a>),
    Loop(MmcifLoopTable<'a>),
}

/// One scalar `_tag value` item.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct MmcifItem<'a> {
    source: &'a str,
    data: &'a ItemData,
}

impl<'a> MmcifItem<'a> {
    pub fn tag(&self) -> &'a str {
        &self.data.tag
    }

    pub fn value(&self) -> MmcifValue<'a> {
        MmcifValue::new(self.source, self.data.value)
    }
}

impl fmt::Debug for MmcifItem<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MmcifItem")
            .field("tag", &self.tag())
            .field("value", &self.value())
            .finish()
    }
}

/// One `loop_` table with complete rows.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct MmcifLoopTable<'a> {
    source: &'a str,
    data: &'a LoopData,
}

impl<'a> MmcifLoopTable<'a> {
    pub fn tags(&self) -> &'a [String] {
        &self.data.tags
    }

    pub fn row_count(&self) -> usize {
        self.data.values.len() / self.data.tags.len()
    }

    /// The column of `tag`, matching it case-insensitively. Resolve columns
    /// once and read cells with [`Self::value_at`] when scanning many rows.
    pub fn column_index(&self, tag: &str) -> Option<usize> {
        self.data
            .tags
            .iter()
            .position(|candidate| candidate.eq_ignore_ascii_case(tag))
    }

    pub fn value(&self, row: usize, tag: &str) -> Option<MmcifValue<'a>> {
        self.value_at(row, self.column_index(tag)?)
    }

    /// The cell at a row and a column from [`Self::column_index`].
    pub fn value_at(&self, row: usize, column: usize) -> Option<MmcifValue<'a>> {
        if row >= self.row_count() || column >= self.data.tags.len() {
            return None;
        }
        // Parsed loops contain complete rows, so the offset is in bounds.
        let span = self.data.values[row * self.data.tags.len() + column];
        Some(MmcifValue::new(self.source, span))
    }

    /// The cells of one row in tag order.
    pub fn row(&self, row: usize) -> Option<impl ExactSizeIterator<Item = MmcifValue<'a>> + 'a> {
        if row >= self.row_count() {
            return None;
        }
        let width = self.data.tags.len();
        let source = self.source;
        let cells = &self.data.values[row * width..(row + 1) * width];
        Some(cells.iter().map(move |span| MmcifValue::new(source, *span)))
    }
}

impl fmt::Debug for MmcifLoopTable<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MmcifLoopTable")
            .field("tags", &self.tags())
            .field("rows", &self.row_count())
            .finish()
    }
}

/// One parsed value with its source line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MmcifValue<'a> {
    text: &'a str,
    line: usize,
    bare: bool,
}

impl<'a> MmcifValue<'a> {
    fn new(source: &'a str, span: Span) -> Self {
        Self {
            text: span.text(source),
            line: span.line as usize,
            bare: span.bare,
        }
    }

    pub fn text(&self) -> &'a str {
        self.text
    }

    pub const fn line(&self) -> usize {
        self.line
    }

    pub const fn location(&self) -> SourceLocation {
        SourceLocation::at_line(self.line)
    }

    /// Whether this is an unquoted `.` (inapplicable) or `?` (unknown).
    pub fn is_missing(&self) -> bool {
        self.bare && matches!(self.text, "." | "?")
    }

    pub fn optional_text(&self) -> Option<&'a str> {
        (!self.is_missing()).then_some(self.text)
    }
}

/// Reads and parses CIF/mmCIF text, enforcing the byte limit while reading.
pub(crate) fn parse_mmcif_reader(
    reader: impl std::io::Read,
    options: MmcifParseOptions,
) -> Result<MmcifDocument, MmcifParseError> {
    let input = read_bounded(reader, options.max_input_bytes)
        .map_err(|(kind, message)| MmcifParseError::new(kind, 1, message))?;
    parse_mmcif_str(&input, options)
}

pub(crate) fn parse_mmcif_str(
    input: &str,
    options: MmcifParseOptions,
) -> Result<MmcifDocument, MmcifParseError> {
    if input.len() > options.max_input_bytes {
        return Err(MmcifParseError::new(
            DiagnosticKind::Limit,
            1,
            "input exceeds configured byte limit",
        ));
    }
    if u32::try_from(input.len()).is_err() {
        return Err(MmcifParseError::new(
            DiagnosticKind::Limit,
            1,
            "input exceeds the 4 GiB addressable mmCIF limit",
        ));
    }
    let source: Arc<str> = normalize_line_endings(input).into();
    let mut tokens = Tokens {
        tokenizer: Tokenizer::new(&source, options),
        peeked: None,
    };
    parse_blocks(&source, &mut tokens, options.max_atom_site_rows)
}

fn parse_blocks(
    source: &Arc<str>,
    tokens: &mut Tokens<'_>,
    max_atom_site_rows: usize,
) -> Result<MmcifDocument, MmcifParseError> {
    let text = |span: Span| span.text(source);
    let mut blocks = Vec::new();
    let mut block_names = BTreeSet::new();
    while let Some(token) = tokens.next()? {
        if !has_bare_prefix(token, text(token), "data_") {
            return Err(syntax_error(
                token,
                "content appears before the first data block",
            ));
        }
        let name = text(token)[5..].to_owned();
        if name.is_empty() {
            return Err(syntax_error(token, "data block has no name"));
        }
        if !block_names.insert(name.to_ascii_lowercase()) {
            return Err(MmcifParseError::new(
                DiagnosticKind::Inconsistent,
                token.line as usize,
                format!("duplicate data block {name}"),
            ));
        }
        let mut entries = Vec::new();
        let mut data_names = BTreeSet::new();
        while let Some(token) = tokens.peek()? {
            if has_bare_prefix(token, text(token), "data_") {
                break;
            }
            tokens.next()?;
            if is_bare_word(token, text(token), "loop_") {
                let table = parse_loop(source, tokens, token, max_atom_site_rows)?;
                for tag in &table.tags {
                    if !data_names.insert(tag.to_ascii_lowercase()) {
                        return Err(MmcifParseError::new(
                            DiagnosticKind::Inconsistent,
                            token.line as usize,
                            format!("duplicate data name {tag}"),
                        ));
                    }
                }
                entries.push(EntryData::Loop(table));
            } else if token.bare && text(token).starts_with('_') {
                let tag = text(token).to_owned();
                if !data_names.insert(tag.to_ascii_lowercase()) {
                    return Err(MmcifParseError::new(
                        DiagnosticKind::Inconsistent,
                        token.line as usize,
                        format!("duplicate data name {tag}"),
                    ));
                }
                let value = tokens
                    .next()?
                    .ok_or_else(|| syntax_error(token, format!("missing value for {tag}")))?;
                if is_control(value, text(value)) {
                    return Err(syntax_error(value, format!("missing value for {tag}")));
                }
                entries.push(EntryData::Item(ItemData { tag, value }));
            } else {
                return Err(syntax_error(
                    token,
                    format!("unexpected token `{}` in data block", text(token)),
                ));
            }
        }
        // Scalar atom-site fields describe one row, subject to the same bound
        // as a one-row loop. Other category semantics belong to interpretation.
        if max_atom_site_rows == 0 {
            if let Some(item) = entries.iter().find_map(|entry| match entry {
                EntryData::Item(item)
                    if super::mmcif_category::in_category(&item.tag, "_atom_site") =>
                {
                    Some(item)
                }
                _ => None,
            }) {
                return Err(MmcifParseError::new(
                    DiagnosticKind::Limit,
                    item.value.line as usize,
                    "atom-site row count exceeds configured limit",
                ));
            }
        }
        blocks.push(MmcifBlock {
            name,
            source: Arc::clone(source),
            entries,
        });
    }
    if blocks.is_empty() {
        return Err(MmcifParseError::new(
            DiagnosticKind::Syntax,
            1,
            "missing data block",
        ));
    }
    Ok(MmcifDocument { blocks })
}

/// Parses one loop after its `loop_` keyword. An atom-site loop fails at the
/// first value of its first row beyond `max_atom_site_rows`, before buffering
/// the rest of the loop.
fn parse_loop(
    source: &str,
    tokens: &mut Tokens<'_>,
    keyword: Span,
    max_atom_site_rows: usize,
) -> Result<LoopData, MmcifParseError> {
    let text = |span: Span| span.text(source);
    let mut tags = Vec::new();
    let mut seen = BTreeSet::new();
    while let Some(token) = tokens.peek()? {
        if !(token.bare && text(token).starts_with('_')) {
            break;
        }
        tokens.next()?;
        if !seen.insert(text(token).to_ascii_lowercase()) {
            return Err(MmcifParseError::new(
                DiagnosticKind::Inconsistent,
                token.line as usize,
                format!("duplicate loop tag {}", text(token)),
            ));
        }
        tags.push(text(token).to_owned());
    }
    if tags.is_empty() {
        return Err(syntax_error(keyword, "loop has no tags"));
    }
    let max_values = tags
        .iter()
        .any(|tag| tag.to_ascii_lowercase().starts_with("_atom_site."))
        .then(|| max_atom_site_rows.saturating_mul(tags.len()));
    let mut values = Vec::new();
    while let Some(token) = tokens.peek()? {
        if is_control(token, text(token)) {
            break;
        }
        if max_values == Some(values.len()) {
            return Err(MmcifParseError::new(
                DiagnosticKind::Limit,
                token.line as usize,
                "atom-site row count exceeds configured limit",
            ));
        }
        tokens.next()?;
        values.push(token);
    }
    if values.len() % tags.len() != 0 {
        let line = values.last().map_or(keyword.line, |value| value.line);
        return Err(MmcifParseError::new(
            DiagnosticKind::Syntax,
            line as usize,
            "loop has ragged rows",
        ));
    }
    if let Some(token) = tokens.peek()? {
        if is_bare_word(token, text(token), "stop_") {
            tokens.next()?;
        }
    }
    Ok(LoopData { tags, values })
}

fn syntax_error(token: Span, message: impl Into<String>) -> MmcifParseError {
    MmcifParseError::new(DiagnosticKind::Syntax, token.line as usize, message)
}

fn is_control(token: Span, text: &str) -> bool {
    token.bare
        && (text.eq_ignore_ascii_case("loop_")
            || text.eq_ignore_ascii_case("stop_")
            || text.eq_ignore_ascii_case("global_")
            || text
                .get(..5)
                .is_some_and(|prefix| prefix.eq_ignore_ascii_case("data_"))
            || text.starts_with('_')
            || text
                .get(..5)
                .is_some_and(|prefix| prefix.eq_ignore_ascii_case("save_")))
}

fn is_bare_word(token: Span, text: &str, word: &str) -> bool {
    token.bare && text.eq_ignore_ascii_case(word)
}

fn has_bare_prefix(token: Span, text: &str, prefix: &str) -> bool {
    token.bare
        && text
            .get(..prefix.len())
            .is_some_and(|value| value.eq_ignore_ascii_case(prefix))
}

/// One-token lookahead over the lazy tokenizer.
struct Tokens<'a> {
    tokenizer: Tokenizer<'a>,
    peeked: Option<Option<Span>>,
}

impl Tokens<'_> {
    fn peek(&mut self) -> Result<Option<Span>, MmcifParseError> {
        if self.peeked.is_none() {
            self.peeked = Some(self.tokenizer.next_token()?);
        }
        Ok(self.peeked.expect("peeked token is present"))
    }

    fn next(&mut self) -> Result<Option<Span>, MmcifParseError> {
        match self.peeked.take() {
            Some(token) => Ok(token),
            None => self.tokenizer.next_token(),
        }
    }
}

/// Splits line-normalized CIF text into value spans without copying it.
struct Tokenizer<'a> {
    source: &'a str,
    line_start: usize,
    /// End of the physical line at `line_start`, found once per line so that
    /// tokenizing a line is linear in its length.
    line_end: usize,
    line: usize,
    column: usize,
    count: usize,
    options: MmcifParseOptions,
    /// Number of line-end searches, which pins that each line is searched once.
    #[cfg(test)]
    line_end_searches: usize,
}

impl<'a> Tokenizer<'a> {
    fn new(source: &'a str, options: MmcifParseOptions) -> Self {
        let mut tokenizer = Self {
            source,
            line_start: 0,
            line_end: 0,
            line: 1,
            column: 0,
            count: 0,
            options,
            #[cfg(test)]
            line_end_searches: 0,
        };
        tokenizer.enter_line(0);
        tokenizer
    }

    /// Makes `start`, the start of a physical line or one past the end of the
    /// source, the current line.
    fn enter_line(&mut self, start: usize) {
        self.line_start = start;
        if start <= self.source.len() {
            self.line_end = self.line_end(start);
        }
    }

    fn next_token(&mut self) -> Result<Option<Span>, MmcifParseError> {
        loop {
            if self.line_start > self.source.len() {
                return Ok(None);
            }
            let line_end = self.line_end;
            if self.column == 0 && self.source[self.line_start..].starts_with(';') {
                return self.semicolon_text(line_end).map(Some);
            }
            let line = &self.source[self.line_start..line_end];
            let bytes = line.as_bytes();
            while self.column < bytes.len() && bytes[self.column].is_ascii_whitespace() {
                self.column += 1;
            }
            if self.column == bytes.len() || bytes[self.column] == b'#' {
                self.enter_line(line_end + 1);
                self.line += 1;
                self.column = 0;
                continue;
            }
            let start = self.column;
            let (value_start, value_end, bare) =
                if bytes[self.column] == b'\'' || bytes[self.column] == b'"' {
                    let quote = bytes[self.column];
                    self.column += 1;
                    let value_start = self.column;
                    // In CIF 1.1 an embedded matching quote belongs to the value
                    // unless it is followed by whitespace or the end of the line.
                    while self.column < bytes.len()
                        && !(bytes[self.column] == quote
                            && bytes
                                .get(self.column + 1)
                                .is_none_or(u8::is_ascii_whitespace))
                    {
                        self.column += 1;
                    }
                    if self.column == bytes.len() {
                        return Err(MmcifParseError::new(
                            DiagnosticKind::Syntax,
                            self.line,
                            "unterminated quoted value",
                        ));
                    }
                    let value_end = self.column;
                    self.column += 1;
                    (value_start, value_end, false)
                } else {
                    while self.column < bytes.len() && !bytes[self.column].is_ascii_whitespace() {
                        self.column += 1;
                    }
                    (start, self.column, true)
                };
            return self
                .span(
                    self.line_start + value_start,
                    self.line_start + value_end,
                    self.line,
                    bare,
                )
                .map(Some);
        }
    }

    /// A `;`-delimited text field starting at the current line. CIF text starts
    /// immediately after the opening semicolon; only the newline introducing
    /// the closing semicolon is excluded.
    fn semicolon_text(&mut self, opening_end: usize) -> Result<Span, MmcifParseError> {
        let opening_line = self.line;
        let text_start = self.line_start + 1;
        let mut next = opening_end + 1;
        let mut line = self.line + 1;
        let closing = loop {
            if next > self.source.len() {
                return Err(MmcifParseError::new(
                    DiagnosticKind::Syntax,
                    opening_line,
                    "unterminated semicolon text",
                ));
            }
            if self.source[next..].starts_with(';') {
                break next;
            }
            let end = self.line_end(next);
            if end - text_start > self.options.max_token_bytes {
                return Err(MmcifParseError::new(
                    DiagnosticKind::Limit,
                    opening_line,
                    "semicolon value exceeds configured token limit",
                ));
            }
            next = end + 1;
            line += 1;
        };
        let text_end = closing - 1;
        if text_end - text_start > self.options.max_token_bytes {
            return Err(MmcifParseError::new(
                DiagnosticKind::Limit,
                opening_line,
                "semicolon value exceeds configured token limit",
            ));
        }
        // A closing delimiter can be followed by more tokens on its line.
        if self
            .source
            .as_bytes()
            .get(closing + 1)
            .is_some_and(|byte| !byte.is_ascii_whitespace())
        {
            return Err(MmcifParseError::new(
                DiagnosticKind::Syntax,
                line,
                "closing semicolon must be followed by whitespace or end of line",
            ));
        }
        self.enter_line(closing);
        self.line = line;
        self.column = 1;
        self.span(text_start, text_end, opening_line, false)
    }

    fn line_end(&mut self, start: usize) -> usize {
        #[cfg(test)]
        {
            self.line_end_searches += 1;
        }
        self.source[start..]
            .find('\n')
            .map_or(self.source.len(), |offset| start + offset)
    }

    fn span(
        &mut self,
        start: usize,
        end: usize,
        line: usize,
        bare: bool,
    ) -> Result<Span, MmcifParseError> {
        if end - start > self.options.max_token_bytes {
            return Err(MmcifParseError::new(
                DiagnosticKind::Limit,
                line,
                "value exceeds configured token limit",
            ));
        }
        if self.count >= self.options.max_tokens {
            return Err(MmcifParseError::new(
                DiagnosticKind::Limit,
                line,
                "input exceeds configured token count",
            ));
        }
        self.count += 1;
        let narrow = |value: usize| u32::try_from(value).expect("input length fits u32");
        Ok(Span {
            start: narrow(start),
            end: narrow(end),
            line: u32::try_from(line).unwrap_or(u32::MAX),
            bare,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_limits_admit_the_default_atom_site_row_limit() {
        let options = MmcifParseOptions::default();
        assert!(
            options.max_tokens >= options.max_atom_site_rows * (TYPICAL_ATOM_SITE_COLUMNS + 1),
            "token budget must cover a full default-sized atom-site table"
        );
        assert!(u32::try_from(options.max_input_bytes).is_ok());
    }

    #[test]
    fn tokenizer_searches_each_line_end_once() {
        let repeats = 1_000;
        let source = format!(
            "data_x\nloop_\n_x.a\n_x.b\n_x.c\n{}# trailing comment\n;text\n; 2 3\n",
            "1 'a b' \"c\" ".repeat(repeats)
        );
        let mut tokenizer = Tokenizer::new(&source, MmcifParseOptions::default());
        let mut tokens = Vec::new();
        while let Some(span) = tokenizer.next_token().unwrap() {
            tokens.push((span.text(&source), span.line, span.bare));
        }
        let mut expected = vec![
            ("data_x", 1, true),
            ("loop_", 2, true),
            ("_x.a", 3, true),
            ("_x.b", 4, true),
            ("_x.c", 5, true),
        ];
        for _ in 0..repeats {
            expected.extend([("1", 6, true), ("a b", 6, false), ("c", 6, false)]);
        }
        expected.extend([("text", 7, false), ("2", 8, true), ("3", 8, true)]);
        assert_eq!(tokens, expected);
        // Eight lines plus the empty remainder after the final newline.
        assert_eq!(tokenizer.line_end_searches, 9);
    }
}
