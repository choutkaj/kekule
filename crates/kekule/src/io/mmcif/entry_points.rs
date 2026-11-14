//! The public mmCIF entry points that `crate::mmcif` re-exports.

use crate::diagnostics::FormatError;
use crate::io::{MmcifBlockSource, MmcifDocument, MmcifParseOptions, MmcifWriteOptions};

/// Parses a structural mmCIF data document without assigning molecular
/// meaning. `options` sets the syntax and resource limits.
///
/// # Errors
///
/// Returns a [`FormatError`] of kind
/// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
/// when the input exceeds `options`,
/// [`DiagnosticKind::Syntax`](crate::diagnostics::DiagnosticKind::Syntax)
/// for malformed CIF, and
/// [`DiagnosticKind::Inconsistent`](crate::diagnostics::DiagnosticKind::Inconsistent)
/// for contradictory structure, such as a repeated block or item.
pub fn parse_str(input: &str, options: MmcifParseOptions) -> Result<MmcifDocument, FormatError> {
    crate::io::parse_mmcif_str(input, options)
}

/// Reads and parses CIF/mmCIF text from `reader`.
///
/// # Errors
///
/// As [`parse_str`]. The configured byte limit is enforced while reading,
/// before the input is held in memory; a read failure or non-UTF-8 input
/// fails with [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
pub fn parse_reader(
    reader: impl std::io::Read,
    options: MmcifParseOptions,
) -> Result<MmcifDocument, FormatError> {
    crate::io::parse_mmcif_reader(reader, options)
}

/// Writes one data block per source and returns the text with a
/// [`WriteReport`](crate::diagnostics::WriteReport) of what the options
/// chose to omit, summed over the blocks.
///
/// Sources convert from structures (one block), structure sets (one
/// multi-model block), and interpretations (keeping their source reports).
/// Entity kinds derive from canonical topology classification unless a
/// source carries explicit classifications or reports
/// ([`MmcifBlockSource`]). Blocks are named by
/// [`MmcifWriteOptions::block_name`], numbered when there are several,
/// unless [`MmcifBlockSource::with_block_name`] names them.
///
/// As the PDB's own files do, a block states as `_struct_conn` rows only
/// the bonds its reader cannot derive. A residue whose bundled component
/// template, or the protonation variant its hydrogens select, gives it
/// exactly its bonds is left to the template, and the peptide or
/// phosphodiester link between residues adjacent in sequence is left to
/// polymer linkage when it lies within bonding distance in the written
/// coordinates. Every other bond is stated: links between residues,
/// dative bonds, and all bonds of a residue without such a template,
/// such as a ligand or a structure read from XYZ. Reading the block back
/// gives every bond with its order, though bonds may come back in another
/// order within a molecule. Each partner is selected by label and author
/// sequence numbers, so that waters of one asym are told apart.
///
/// Metal coordination between residues is no bond of the structure: an
/// interpretation report keeps it
/// ([`StructureInterpretationReport::coordination`](crate::macromolecular::StructureInterpretationReport::coordination)).
/// A block written with its report states it as `metalc` rows; a block
/// written from a structure alone has none. A record whose partners the
/// structure no longer has as the report recorded them is left out and
/// counted as
/// [`WriteOmission::Coordination`](crate::diagnostics::WriteOmission::Coordination).
///
/// # Errors
///
/// Returns a [`FormatError`] of kind
/// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
/// when mmCIF cannot hold the content, and `Inconsistent`,
/// `InvalidReference` or `InvalidValue` when a source's classifications or
/// reports do not fit its topology, and `InvalidValue` when two blocks
/// would share a name.
///
/// ```no_run
/// use kekule::mmcif::{self, MmcifBlockSource, MmcifWriteOptions};
/// # fn example(structure: &kekule::structure::Structure, other: &kekule::structure::Structure,
/// #     interpretation: &kekule::macromolecular::StructureInterpretation)
/// #     -> Result<(), kekule::diagnostics::FormatError> {
/// let (one, _) = mmcif::write([structure], MmcifWriteOptions::default())?;
/// let (two, _) = mmcif::write([structure, other], MmcifWriteOptions::default())?;
/// let (faithful, _) = mmcif::write([interpretation], MmcifWriteOptions::default())?;
/// # let _ = (one, two, faithful);
/// # Ok(())
/// # }
/// ```
pub fn write<'a, B: Into<MmcifBlockSource<'a>>>(
    blocks: impl IntoIterator<Item = B>,
    options: MmcifWriteOptions,
) -> Result<(String, crate::diagnostics::WriteReport), FormatError> {
    crate::io::write_mmcif(blocks, options)
}

/// Streams the blocks [`write()`] would return and returns its report.
/// Values reach `writer` as
/// they are formatted, so memory stays bounded however large the
/// structure; wrap an unbuffered sink such as a `File` in a
/// [`BufWriter`](std::io::BufWriter).
///
/// # Errors
///
/// The errors of [`write()`]; output before the failure has already
/// reached `writer`. A failing `writer` gives
/// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
pub fn write_to<'a, B: Into<MmcifBlockSource<'a>>>(
    writer: &mut impl std::io::Write,
    blocks: impl IntoIterator<Item = B>,
    options: MmcifWriteOptions,
) -> Result<crate::diagnostics::WriteReport, FormatError> {
    crate::io::write_mmcif_to(writer, blocks, options)
}
