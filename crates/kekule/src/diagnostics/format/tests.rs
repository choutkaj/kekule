use super::*;

#[test]
fn format_errors_display_format_stage_record_and_location_but_not_their_source() {
    let cause = std::io::Error::new(std::io::ErrorKind::BrokenPipe, "closed");
    let error = FormatError::write(Format::Sdf, DiagnosticKind::Io, "")
        .in_record(3)
        .with_source(cause);
    assert_eq!(
        error.to_string(),
        "Cannot write the SDF output in record 3."
    );
    assert_eq!(
        std::error::Error::source(&error)
            .and_then(|source| source.downcast_ref::<std::io::Error>())
            .map(std::io::Error::kind),
        Some(std::io::ErrorKind::BrokenPipe)
    );

    let located = FormatError::parse(
        Format::Smiles,
        DiagnosticKind::Syntax,
        "Ring bond 1 is not closed.",
    )
    .at(SourceLocation::at_column(1, 4))
    .on_line(7)
    .shift_column(2);
    assert_eq!(located.location(), Some(SourceLocation::at_column(7, 6)));
    assert_eq!(located.record(), None);
    assert_eq!(
        located.to_string(),
        "Cannot parse the SMILES input at line 7, column 6. Ring bond 1 is not closed."
    );
    assert!(std::error::Error::source(&located).is_none());

    let lined = FormatError::interpret(
        Format::Pdb,
        DiagnosticKind::Inconsistent,
        "Atom 3 appears twice.",
    )
    .on_line(2);
    assert_eq!(lined.location(), Some(SourceLocation::at_line(2)));
    assert_eq!(
        lined.to_string(),
        "Cannot interpret the PDB input at line 2. Atom 3 appears twice."
    );
}

#[test]
fn locations_display_line_and_optional_column() {
    assert_eq!(SourceLocation::at_line(4).to_string(), "line 4");
    assert_eq!(
        SourceLocation::at_column(4, 7).to_string(),
        "line 4, column 7"
    );
}
