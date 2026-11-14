//! The public PDB entry points that `crate::pdb` re-exports.

use crate::diagnostics::FormatError;
use crate::io::{PdbDocument, PdbParseOptions, PdbSource, PdbWriteOptions};

/// Parses a PDB entry into typed records without assigning molecular
/// meaning. `options` sets the resource limits.
///
/// # Errors
///
/// Returns a [`FormatError`] of kind
/// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
/// when the input exceeds `options`, and `Syntax`, `InvalidValue`,
/// `Inconsistent` or `Unsupported` when a typed record is malformed,
/// contradicts another, or uses content Kekule does not read.
pub fn parse_str(input: &str, options: PdbParseOptions) -> Result<PdbDocument, FormatError> {
    crate::io::parse_pdb_str(input, options)
}

/// Reads and parses a PDB entry from `reader`.
///
/// # Errors
///
/// As [`parse_str`]. The configured byte limit is enforced while reading,
/// before the input is held in memory; a read failure or non-UTF-8 input
/// fails with [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
pub fn parse_reader(
    reader: impl std::io::Read,
    options: PdbParseOptions,
) -> Result<PdbDocument, FormatError> {
    crate::io::parse_pdb_reader(reader, options)
}

/// Writes one entry: a structure or geometry without `MODEL` records, or a
/// structure set with one `MODEL` per geometry. Returns the text with a
/// [`WriteReport`](crate::diagnostics::WriteReport) of what the options
/// chose to omit.
///
/// Records follow dense atom order, which interpretation keeps, so
/// interpreting the output gives the same atoms in the same order.
/// Identity comes from the hierarchy:
/// - author chain, residue, and atom identifiers, falling back to label
///   identifiers;
/// - residues of protein and nucleic-acid instances with an amino-acid or
///   nucleotide class are `ATOM` records, everything else `HETATM`;
/// - a `TER` record, which takes the next serial number, follows each
///   chain's last polymer atom;
/// - an instance without atom sites becomes one `MOL` (or `HOH`) residue
///   without a chain identifier, and atoms without a name are named by
///   their element and 1-based molecule-local number.
///
/// Serial numbers beyond 99,999 and residue numbers beyond 9,999 use
/// hybrid-36. Coordinates are written to 0.001 Å, occupancies and B
/// factors to 0.01, and absent occupancies and B factors are blank. A
/// representable periodic cell is a `CRYST1` record in space group P 1.
///
/// `CONECT` lists every bond a reader could not derive, each from both
/// atoms and once per step of bond order:
/// - bonds within `HETATM` residues and residues the bundled
///   [`ComponentLibrary`](crate::components::ComponentLibrary) does not
///   cover;
/// - bonds between residues other than the linkage of standard residues
///   adjacent in one chain, such as disulfides and ligand attachments.
///
/// Bonds within standard residues the dictionary covers, polymer
/// linkage, and bonds within water are derived by readers and not listed.
///
/// Metal coordination between residues is no bond of the structure: an
/// interpretation report keeps it
/// ([`StructureInterpretationReport::coordination`](crate::macromolecular::StructureInterpretationReport::coordination)).
/// An entry written with its reports states it as `LINK` records before
/// the coordinates ([`PdbSource::with_reports`]); an entry written from
/// a structure alone has none. Records the entry cannot state are
/// counted as
/// [`WriteOmission::Coordination`](crate::diagnostics::WriteOmission::Coordination).
///
/// # Errors
///
/// A write-stage [`FormatError`] of kind
/// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
/// with a [`PdbWriteIssue`](crate::pdb::PdbWriteIssue) when PDB cannot hold the content, and
/// `Inconsistent` for an empty structure set or a report count that
/// does not match the structures. Identifiers
/// that do not fit their columns fail: chain identifiers longer
/// than one character (unless [`PdbChainIdPolicy::SegmentId`](crate::pdb::PdbChainIdPolicy::SegmentId) writes one
/// of up to four as a segment identifier), residue names longer than three, atom names longer
/// than four, insertion codes longer than one, residue numbers outside
/// -999..=2,436,111, more than 87,440,031 atom and `TER` records, and
/// residues whose atoms are not consecutive. So do values outside their
/// columns and formal charges outside -9..=9. Fixed implicit-hydrogen
/// counts, other represented chemistry, ensemble or frame state, and
/// cells PDB cannot carry fail unless [`PdbWriteOptions`] omits them,
/// each under its own policy, and every omission is in the report.
/// Nothing is written when validation fails.
///
/// ```
/// use kekule::pdb;
///
/// let text = "\
/// HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
/// HETATM    2  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H
/// HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H
/// END
/// ";
/// let structure = pdb::parse_str(text, Default::default())?
///     .interpret(Default::default())?
///     .into_structure();
/// let (written, report) = pdb::write(&structure, Default::default())?;
/// assert_eq!(written, text);
/// assert!(report.omissions().is_empty());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// An interpretation keeps its report, and so its coordination; a structure
/// or a structure set takes reports from [`PdbSource::with_reports`]:
///
/// ```no_run
/// use kekule::pdb::{self, PdbSource};
/// # fn example(interpretation: &kekule::macromolecular::StructureInterpretation,
/// #     set: &kekule::macromolecular::StructureSetInterpretation)
/// #     -> Result<(), kekule::diagnostics::FormatError> {
/// let (entry, _) = pdb::write(interpretation, Default::default())?;
/// let (set_entry, _) = pdb::write(
///     PdbSource::structure_set(set.structure_set()).with_reports(set.reports()),
///     Default::default(),
/// )?;
/// # let _ = (entry, set_entry);
/// # Ok(())
/// # }
/// ```
pub fn write<'a>(
    source: impl Into<PdbSource<'a>>,
    options: PdbWriteOptions,
) -> Result<(String, crate::diagnostics::WriteReport), FormatError> {
    crate::io::write_pdb(source, options)
}

/// Streams the entry [`write()`] would return and returns its report.
/// The entry is validated
/// before the first record reaches `writer`; wrap an unbuffered sink
/// such as a `File` in a [`BufWriter`](std::io::BufWriter).
///
/// # Errors
///
/// The errors of [`write()`]; a failing `writer` gives
/// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
pub fn write_to<'a>(
    writer: &mut impl std::io::Write,
    source: impl Into<PdbSource<'a>>,
    options: PdbWriteOptions,
) -> Result<crate::diagnostics::WriteReport, FormatError> {
    crate::io::write_pdb_to(writer, source, options)
}
