//! The public SMILES entry points that `crate::smiles` re-exports: record
//! parsing, the parse-and-interpret shortcuts, and record writing with its
//! options and sources.

use std::sync::Arc;

use super::canonical::write_canonical_smiles;
use super::parse::{parse_smiles_document, parse_smiles_reader};
use super::write::{write_isomeric_smiles, write_smiles};
use super::{SmilesDocument, SmilesInterpretOptions, SmilesParseOptions, SmilesRecords};
use crate::diagnostics::{DiagnosticKind, Format, FormatError};
use crate::molecule::Molecule;
use crate::topology::Topology;

/// SMILES output policy of [`write()`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum SmilesWriteMode {
    /// Ordinary non-canonical SMILES without stereo. Declared and inferred
    /// hydrogen counts may be combined in a bracket atom; the original
    /// inference policy is not serialized.
    #[default]
    Ordinary,
    /// Non-canonical SMILES preserving represented stereo.
    ///
    /// SMILES has no notation for an explicitly unknown configuration, so
    /// this and canonical mode write such centers, bonds and axes
    /// unmarked, as RDKit does; they read back as unspecified.
    ///
    /// Tetrahedral absolute, AND, OR and relative groups are emitted as
    /// CXSMILES (`a`, `&`, `o`, `r`). Unsupported member geometries,
    /// quantitative racemic groups and multiple independent relative groups
    /// return an error rather than losing relationships. The global `r`
    /// flag shields other specified centers with explicit absolute
    /// membership. Directional encoding is bounded by 4,096 carrier
    /// combinations and 50,000,000 graph visits; exceeding either returns a
    /// resource-limit error.
    Isomeric,
    /// Deterministic canonical isomeric SMILES.
    ///
    /// Atom priorities, stereo refinement, branch traversal, ring closures
    /// and AND/OR representative selection follow RDKit-style conventions;
    /// byte-for-byte RDKit compatibility is not guaranteed, and canonical
    /// spelling can change when these ordering rules are improved. Output is
    /// invariant under atom numbering and preserves supported stereo,
    /// isotopes, formal charges, and atom maps. Neutral unmapped
    /// nonisotopic terminal hydrogen vertices may collapse into hydrogen
    /// counts. Supplied local stereo assertions are retained even when CIP
    /// perception finds no stereogenic unit. Group numbering and member
    /// order are canonical, and inverting every member of a non-absolute
    /// group produces the same output.
    ///
    /// Canonical output of a molecule with a ring requires installed
    /// aromaticity perception, such as [`Molecule::perceive`] installs, and
    /// fails with
    /// [`DiagnosticKind::MissingPerception`]
    /// without it.
    /// Perceived
    /// aromatic atoms and bonds are written aromatic, so every Kekule form
    /// of a molecule gives one string under the installed aromaticity
    /// model. The exception is a molecule with a stereo center or E/Z
    /// endpoint on an aromatic atom: SMILES cannot mark that configuration
    /// on aromatic syntax, so the whole molecule is written in its stored
    /// Kekule form.
    ///
    /// Complete canonical labeling is bounded by 100,000 search states,
    /// 50,000,000 atom/edge/twin visits, and 2,000,000 pending atom labels.
    /// Serialization additionally bounds the input to 2,000,000 combined
    /// atom/bond slots and a complexity score `2*n*(n+2*m)` of at most
    /// 50,000,000 for `n` live atoms and `m` live bonds. Exceeding a bound
    /// returns
    /// [`DiagnosticKind::Limit`]
    /// without a partial result.
    Canonical,
}

/// Options of [`write()`] and [`write_to()`]; the default writes
/// [`SmilesWriteMode::Ordinary`] SMILES.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SmilesWriteOptions {
    mode: SmilesWriteMode,
}

impl SmilesWriteOptions {
    /// The output policy.
    pub const fn mode(&self) -> SmilesWriteMode {
        self.mode
    }

    /// Sets [`Self::mode`].
    #[must_use]
    pub const fn with_mode(mut self, mode: SmilesWriteMode) -> Self {
        self.mode = mode;
        self
    }
}

impl SmilesWriteOptions {
    /// Options for [`SmilesWriteMode::Ordinary`] output.
    pub const fn ordinary() -> Self {
        Self {
            mode: SmilesWriteMode::Ordinary,
        }
    }

    /// Options for [`SmilesWriteMode::Isomeric`] output.
    pub const fn isomeric() -> Self {
        Self {
            mode: SmilesWriteMode::Isomeric,
        }
    }

    /// Options for [`SmilesWriteMode::Canonical`] output.
    pub const fn canonical() -> Self {
        Self {
            mode: SmilesWriteMode::Canonical,
        }
    }
}

/// Parses one record: base SMILES, an optional CX extension with its fields
/// typed, and an optional name. `options` bounds the input and the parsed
/// graph.
///
/// # Errors
///
/// Returns a [`FormatError`] of kind [`DiagnosticKind::Limit`] when the
/// input exceeds `options`, and [`DiagnosticKind::Syntax`],
/// [`DiagnosticKind::InvalidValue`] or [`DiagnosticKind::Unsupported`]
/// when the record is malformed or uses syntax Kekule does not read.
pub fn parse_str(input: &str, options: SmilesParseOptions) -> Result<SmilesDocument, FormatError> {
    parse_smiles_document(input, options)
}

/// Reads and parses one record from `reader`, as [`parse_str`] parses the
/// same text. Use [`records`] for line-delimited multi-record input.
///
/// # Errors
///
/// As [`parse_str`]. The configured byte limit is enforced while reading;
/// a read failure or non-UTF-8 input fails with [`DiagnosticKind::Io`].
pub fn parse_reader(
    reader: impl std::io::Read,
    options: SmilesParseOptions,
) -> Result<SmilesDocument, FormatError> {
    parse_smiles_reader(reader, options)
}

/// Streams the records of a multi-record SMILES input such as a `.smi`
/// file, one record per line, each bounded by `options`.
///
/// Blank lines are skipped. Each item is one record or the failure of one
/// line, located at its input line; reading continues after a malformed
/// line. Memory is bounded by one line.
pub fn records<R: std::io::Read>(reader: R, options: SmilesParseOptions) -> SmilesRecords<R> {
    SmilesRecords::new(reader, options)
}

/// Parses and interprets one SMILES record into source-ordered connected components.
///
/// This is the concise form of [`parse_str`] followed by
/// [`SmilesDocument::interpret`], both with default options.
/// Dot-delimited components remain separate, and no perception is run implicitly.
///
/// # Errors
///
/// The errors of [`parse_str`], and of interpretation: chemistry that
/// cannot be published fails with [`DiagnosticKind::Chemistry`].
pub fn to_molecules(input: &str) -> Result<Vec<Molecule>, FormatError> {
    let document = parse_str(input, SmilesParseOptions::default())?;
    Ok(document
        .interpret(SmilesInterpretOptions::default())?
        .into_molecules())
}

/// Parses and interprets one SMILES record as a shared coordinate-free
/// topology.
///
/// Every connected component becomes one explicit molecule
/// occurrence in source order. No hierarchy or perception is fabricated.
///
/// # Errors
///
/// The errors of [`to_molecules`]. A topology that cannot be assembled
/// fails with [`DiagnosticKind::Chemistry`] and the
/// [`TopologyBuildError`](crate::topology::TopologyBuildError) as its
/// source.
pub fn to_topology(input: &str) -> Result<Arc<Topology>, FormatError> {
    let document = parse_str(input, SmilesParseOptions::default())?;
    let topology = document
        .interpret(SmilesInterpretOptions::default())?
        .into_topology()
        .map_err(|error| {
            FormatError::interpret(
                Format::Smiles,
                DiagnosticKind::Chemistry,
                "Cannot build the topology.",
            )
            .with_source(error)
        })?;
    Ok(Arc::new(topology))
}

/// What one SMILES record describes: one connected molecule, or every
/// molecule occurrence of a topology joined with `.`, with an optional
/// record name written after the SMILES as `.smi` files carry it.
///
/// Convert from `&Molecule`, `&Topology`, or `&Arc<Topology>`, or build
/// one with [`Self::molecule`] or [`Self::topology`] to add a name.
///
/// ```
/// use kekule::smiles::{self, SmilesSource};
/// let ethanol = smiles::to_molecules("CCO")?.remove(0);
/// let line = smiles::write(SmilesSource::molecule(&ethanol).with_name("ethanol"), Default::default())?;
/// assert_eq!(line, "CCO ethanol");
/// let document = smiles::parse_str(&line, Default::default())?;
/// assert_eq!(document.name(), Some("ethanol"));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct SmilesSource<'a> {
    content: SmilesContent<'a>,
    name: Option<&'a str>,
}

#[derive(Debug, Clone, Copy)]
enum SmilesContent<'a> {
    Molecule(&'a Molecule),
    /// Every molecule occurrence of a topology.
    Topology(&'a Topology),
}

impl<'a> SmilesSource<'a> {
    /// One connected molecule.
    pub const fn molecule(molecule: &'a Molecule) -> Self {
        Self {
            content: SmilesContent::Molecule(molecule),
            name: None,
        }
    }

    /// Every molecule occurrence of a topology, joined with `.`.
    pub const fn topology(topology: &'a Topology) -> Self {
        Self {
            content: SmilesContent::Topology(topology),
            name: None,
        }
    }

    /// Writes `name` after the SMILES, separated by one space. Writing
    /// fails for an empty name, one with surrounding whitespace, or one
    /// spanning lines, which would not read back as the same name.
    #[must_use]
    pub const fn with_name(mut self, name: &'a str) -> Self {
        self.name = Some(name);
        self
    }
}

impl<'a> From<&'a Molecule> for SmilesSource<'a> {
    fn from(molecule: &'a Molecule) -> Self {
        Self::molecule(molecule)
    }
}

impl<'a> From<&'a Topology> for SmilesSource<'a> {
    fn from(topology: &'a Topology) -> Self {
        Self::topology(topology)
    }
}

impl<'a> From<&'a Arc<Topology>> for SmilesSource<'a> {
    fn from(topology: &'a Arc<Topology>) -> Self {
        Self::topology(topology)
    }
}

/// Writes one SMILES record.
///
/// [`SmilesWriteMode`] selects ordinary, isomeric, or canonical output.
/// Canonical output of a molecule with a ring requires installed
/// aromaticity perception.
/// Every mode preserves isotope labels and total hydrogen counts. An atom
/// that permits inferred hydrogens and requires bracket syntax must have
/// hydrogen perception installed; otherwise writing returns an error. Call
/// [`Molecule::perceive`] explicitly before exporting such atoms. Radical
/// electrons are encoded by bracket valence, with redundant CX annotations
/// for one, two or three electrons (`^1`, `^2`, `^5`); explicit spin
/// multiplicities remain unsupported.
///
/// A topology emits each occurrence of a reused definition once, joined
/// with `.`; canonical mode sorts components, other modes keep instance
/// order. Enhanced groups share one CX extension with record-global atom
/// indices.
///
/// # Errors
///
/// Returns a [`FormatError`] of kind [`DiagnosticKind::MissingPerception`]
/// when required perception is not installed,
/// [`DiagnosticKind::Unsupported`] when the mode cannot represent the
/// chemistry without loss, and [`DiagnosticKind::Limit`] when a
/// canonicalization or stereo-encoding bound of [`SmilesWriteMode`] is
/// exceeded.
///
/// ```
/// use kekule::smiles::{self, SmilesWriteOptions};
/// let canonical = |source: &str| -> Result<String, Box<dyn std::error::Error>> {
///     let mut molecule = smiles::to_molecules(source)?.remove(0);
///     molecule.perceive()?;
///     Ok(smiles::write(&molecule, SmilesWriteOptions::canonical())?)
/// };
/// assert_eq!(canonical("OCC")?, "CCO");
/// assert_eq!(canonical("CC1=C(C)C=CC=C1")?, canonical("CC1=CC=CC=C1C")?);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn write<'a>(
    source: impl Into<SmilesSource<'a>>,
    options: SmilesWriteOptions,
) -> Result<String, FormatError> {
    let source = source.into();
    if let Some(name) = source.name {
        if name.is_empty() || name.trim() != name || name.contains(['\r', '\n']) {
            return Err(FormatError::write(
                Format::Smiles,
                DiagnosticKind::Unsupported,
                "The record name is empty, has surrounding spaces or spans lines.",
            )
            .with_detail(crate::diagnostics::FormatErrorDetail::Title));
        }
    }
    let mut record = match source.content {
        SmilesContent::Molecule(molecule) => match options.mode {
            SmilesWriteMode::Ordinary => write_smiles(molecule),
            SmilesWriteMode::Isomeric => write_isomeric_smiles(molecule),
            SmilesWriteMode::Canonical => write_canonical_smiles(molecule),
        },
        SmilesContent::Topology(topology) => write_topology(topology, options),
    }?;
    if let Some(name) = source.name {
        record.push(' ');
        record.push_str(name);
    }
    Ok(record)
}

/// Writes one SMILES record as a line terminated by `\n`, so repeated
/// calls produce a line-delimited SMILES file. See [`write()`].
///
/// # Errors
///
/// The errors of [`write()`], before anything reaches `writer`; a failing
/// `writer` gives [`DiagnosticKind::Io`].
pub fn write_to<'a>(
    writer: &mut impl std::io::Write,
    source: impl Into<SmilesSource<'a>>,
    options: SmilesWriteOptions,
) -> Result<(), FormatError> {
    let mut line = write(source, options)?;
    line.push('\n');
    writer
        .write_all(line.as_bytes())
        .map_err(|error| FormatError::sink(Format::Smiles, error))
}

fn write_topology(
    topology: &crate::topology::Topology,
    options: SmilesWriteOptions,
) -> Result<String, crate::diagnostics::FormatError> {
    let mut parts = topology
        .instances()
        .map(|m| match options.mode() {
            SmilesWriteMode::Canonical => super::canonical::write_canonical_emission(m.molecule()),
            SmilesWriteMode::Isomeric => super::write::write_source_order_emission(
                m.molecule(),
                super::write::StereoWriteMode::Encode,
                super::write::CanonicalAtomStyle::StoredKekule,
            ),
            SmilesWriteMode::Ordinary => super::write::write_source_order_emission(
                m.molecule(),
                super::write::StereoWriteMode::Reject,
                super::write::CanonicalAtomStyle::Aromatic,
            ),
        })
        .collect::<Result<Vec<_>, _>>()?;
    if options.mode() == SmilesWriteMode::Canonical {
        let mut keyed = parts
            .into_iter()
            .map(|p| Ok((p.render()?, p)))
            .collect::<Result<Vec<_>, crate::diagnostics::FormatError>>()?;
        keyed.sort_by(|a, b| a.0.cmp(&b.0));
        parts = keyed.into_iter().map(|(_, p)| p).collect();
    }
    super::emit::Emission::join(parts).render()
}
