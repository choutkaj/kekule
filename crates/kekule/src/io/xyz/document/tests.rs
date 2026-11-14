use super::*;

#[test]
fn frames_keep_comments_extra_columns_and_lines() {
    let text = "2\r\n water? \r\nO 0 0 0\r\nH 0.96 0 0 0.41 charge\r\n";
    let document = parse_xyz_str(text, XyzParseOptions::default()).unwrap();
    let [frame] = document.frames() else {
        panic!("one frame")
    };
    assert_eq!(frame.comment(), " water? ");
    assert_eq!(frame.atoms()[1].element(), "H");
    assert_eq!(frame.atoms()[1].position(), Point3::new(0.96, 0.0, 0.0));
    assert_eq!(frame.atoms()[1].extra(), Some("0.41 charge"));
    assert_eq!(frame.atoms()[0].extra(), None);
    assert_eq!(frame.atoms()[1].location().line(), 4);
}

#[test]
fn malformed_frames_fail_with_typed_kinds_and_lines() {
    for (text, kind, line) in [
        ("x\n\n", DiagnosticKind::Syntax, 1),
        ("0\n\n", DiagnosticKind::InvalidValue, 1),
        ("2\n\nHe 0 0 0\n", DiagnosticKind::Syntax, 1),
        ("1\n\nHe 0 0\n", DiagnosticKind::Syntax, 3),
        ("1\n\nHe 0 nan 0\n", DiagnosticKind::InvalidValue, 3),
        ("1\n\nHe 0 0 0\n\n1\n", DiagnosticKind::Syntax, 5),
        ("1\n", DiagnosticKind::Syntax, 1),
    ] {
        let error = parse_xyz_str(text, XyzParseOptions::default()).unwrap_err();
        assert_eq!(error.kind(), kind, "{text:?}");
        assert_eq!(error.location().map(|at| at.line()), Some(line), "{text:?}");
    }
    assert!(parse_xyz_str("", XyzParseOptions::default())
        .unwrap()
        .frames()
        .is_empty());
    assert_eq!(
        parse_xyz_str("1\n\nHe 0 0 0\n\n \n", XyzParseOptions::default())
            .unwrap()
            .frames()
            .len(),
        1
    );
}

#[test]
fn limits_bound_input_and_atom_lines() {
    let options = XyzParseOptions::default().with_max_atom_records(1);
    let error = parse_xyz_str("2\n\nHe 0 0 0\nHe 3 0 0\n", options).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Limit);
    assert_eq!(error.limit_exceeded().unwrap().resource(), Resource::Atoms);
    // A huge declared count fails on the limit before allocating.
    let error = parse_xyz_str("18446744073709551615\n\n", XyzParseOptions::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Limit);
    let options = XyzParseOptions::default().with_max_input_bytes(3);
    let error = parse_xyz_str("1\n\nHe 0 0 0\n", options).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Limit);
}
