//! Focused regressions for the CIF source-syntax boundary.
use crate::macromolecular;
use crate::mmcif::{self, MmcifParseOptions};

// Test-only re-encoding of focused source fixtures. No chemistry is generated.
fn token(value: mmcif::MmcifValue<'_>) -> String {
    let text = value.text();
    if value.is_missing() {
        text.to_owned()
    } else if text.is_empty()
        || text.chars().any(char::is_whitespace)
        || text.starts_with(['_', '#', ';', '\'', '"', '$', '[', ']'])
        || matches!(text, "." | "?")
    {
        assert!(
            !text.contains('\'') && !text.contains('\n'),
            "focused fixture token"
        );
        format!("'{text}'")
    } else {
        text.to_owned()
    }
}

pub(crate) fn singleton_loops_as_scalars(input: &str) -> String {
    let doc = mmcif::parse_str(input).unwrap();
    let mut output = String::new();
    for block in doc.blocks() {
        output.push_str(&format!("data_{}\n", block.name()));
        for entry in block.entries() {
            match entry {
                mmcif::MmcifEntry::Item(item) => {
                    output.push_str(&format!("{} {}\n", item.tag(), token(item.value())))
                }
                mmcif::MmcifEntry::Loop(table) if table.row_count() == 1 => {
                    for tag in table.tags() {
                        output
                            .push_str(&format!("{tag} {}\n", token(table.value(0, tag).unwrap())));
                    }
                }
                mmcif::MmcifEntry::Loop(table) => {
                    output.push_str(&format!("loop_\n{}\n", table.tags().join("\n")));
                    for row in 0..table.row_count() {
                        for value in table.row(row).unwrap() {
                            output.push_str(&format!("{}\n", token(value)));
                        }
                    }
                }
            }
        }
    }
    output
}

const SINGLE_SITE: &str = "data_x\nloop_\n_entity.id\n_entity.type\n1 water\nloop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nO O HOH A 1 1.0 2.0 3.0\n";

#[test]
fn mmcif_scalar_categories_support_document_block_and_ensemble_paths() {
    let looped = mmcif::parse_str(SINGLE_SITE).unwrap();
    let scalar = mmcif::parse_str(&singleton_loops_as_scalars(SINGLE_SITE)).unwrap();
    assert!(scalar.blocks()[0]
        .loop_with_tag("_atom_site.type_symbol")
        .is_none());
    assert!(scalar.blocks()[0].item("_atom_site.type_symbol").is_some());
    let expected = looped.interpret().unwrap();
    for actual in [
        scalar.interpret().unwrap(),
        scalar.blocks()[0].interpret().unwrap(),
    ] {
        assert!(expected.topology().same_layout(actual.topology()));
        assert_eq!(expected.model().positions(), actual.model().positions());
        assert_eq!(actual.model().atom_count(), 1);
        assert_eq!(actual.report().entity_definitions(), 1);
        assert_eq!(actual.report().solvent_molecules(), 1);
        assert_eq!(actual.model().topology().residues().count(), 1);
    }
    let ensemble = scalar.interpret_ensemble().unwrap();
    assert_eq!(ensemble.ensemble().len(), 1);
    assert!(ensemble.topology().same_layout(expected.topology()));
    assert_eq!(
        ensemble
            .ensemble()
            .get(0)
            .unwrap()
            .as_model_view()
            .positions(),
        expected.model().positions()
    );
    let options = MmcifParseOptions {
        max_atom_site_rows: 0,
        ..Default::default()
    };
    for input in [
        SINGLE_SITE.to_owned(),
        singleton_loops_as_scalars(SINGLE_SITE),
    ] {
        assert!(mmcif::parse_str_with_options(&input, options)
            .unwrap_err()
            .message()
            .contains("row count"));
    }
}

#[test]
fn mmcif_present_categories_missing_required_fields_reject_in_both_encodings() {
    // Every currently consumed category is discovered by category, even when
    // the old identifying column is absent. No malformed category vanishes.
    for (category, field, value) in [
        ("_entity", "type", "water"),
        ("_struct_asym", "entity_id", "1"),
        ("_atom_site", "label_atom_id", "O"),
        ("_entity_poly", "type", "polypeptide(L)"),
        ("_chem_comp_bond", "atom_id_1", "N"),
        ("_pdbx_branch_scheme", "asym_id", "A"),
        ("_pdbx_entity_branch_link", "atom_id_1", "N"),
        ("_pdbx_poly_seq_scheme", "seq_id", "1"),
        ("_struct_conn", "id", "link"),
    ] {
        let source = if matches!(category, "_entity" | "_struct_asym" | "_atom_site") {
            // Remove this category from the valid fixture through document entries.
            let doc = mmcif::parse_str(SINGLE_SITE).unwrap();
            let mut s = String::from("data_x\n");
            for entry in doc.blocks()[0].entries() {
                if let mmcif::MmcifEntry::Loop(table) = entry {
                    if !table.tags()[0].starts_with(&format!("{category}.")) {
                        s.push_str(&format!(
                            "loop_\n{}\n{}\n",
                            table.tags().join("\n"),
                            table
                                .row(0)
                                .unwrap()
                                .map(token)
                                .collect::<Vec<_>>()
                                .join(" ")
                        ));
                    }
                }
            }
            s
        } else {
            SINGLE_SITE.to_owned()
        };
        let input = format!("{source}loop_\n{category}.{field}\n{value}\n");
        let loop_error = mmcif::parse_str(&input).unwrap().interpret().unwrap_err();
        let scalar_doc = mmcif::parse_str(&singleton_loops_as_scalars(&input)).unwrap();
        let scalar_error = scalar_doc.interpret().unwrap_err();
        assert_eq!(loop_error.message(), scalar_error.message(), "{category}");
        assert!(
            scalar_error.message().contains("missing required"),
            "{category}: {scalar_error}"
        );
        assert!(loop_error.location().is_some() && scalar_error.location().is_some());
        assert!(scalar_doc.interpret_ensemble().is_err());
    }
}

#[test]
fn mmcif_category_ambiguity_is_rejected_without_normalizing_the_document() {
    for (extra, message) in [
        (
            "_entity_poly.entity_id 1\nloop_\n_entity_poly.type\npolypeptide(L)\n",
            "mixes scalar",
        ),
        (
            "loop_\n_entity_poly.entity_id\n1\nloop_\n_entity_poly.type\npolypeptide(L)\n",
            "multiple loops",
        ),
        (
            "loop_\n_entity_poly.entity_id\n_other.type\n1 polypeptide(L)\n",
            "another category",
        ),
    ] {
        let doc = mmcif::parse_str(&format!("{SINGLE_SITE}{extra}")).unwrap();
        let error = doc.interpret().unwrap_err();
        assert!(error.message().contains(message), "{error}");
        assert!(error.location().is_some());
        assert!(doc.interpret_ensemble().is_err());
    }
    let doc = mmcif::parse_str(&format!(
        "{SINGLE_SITE}data_broken\n_atom_site.label_atom_id O\n"
    ))
    .unwrap();
    assert!(doc
        .interpret()
        .unwrap_err()
        .message()
        .contains("more than one data block"));
    assert!(matches!(
        doc.interpret_ensemble(),
        Err(macromolecular::StructureEnsembleInterpretError::MultipleAtomSiteBlocks)
    ));
}

#[test]
fn mmcif_1ake_preserves_declared_sequence_and_polymer_connectivity() {
    const SOURCE: &str = include_str!("../../tests/fixtures/mmcif/1AKE.cif");
    let doc = mmcif::parse_str(SOURCE).unwrap();
    let sequence = doc.blocks()[0]
        .item("_entity_poly.pdbx_seq_one_letter_code_can")
        .unwrap();
    assert!(sequence.text().starts_with("MRIILLGAPGAGKGT"));
    assert_eq!(
        sequence
            .text()
            .chars()
            .filter(|c| !c.is_whitespace())
            .count(),
        214
    );
    // Re-encode the complete scalar category in place, preserving raw values.
    let start = SOURCE.find("_entity_poly.entity_id").unwrap();
    let end = start + SOURCE[start..].find("\n#").unwrap();
    let category = &SOURCE[start..end];
    let tags = doc.blocks()[0]
        .entries()
        .filter_map(|entry| match entry {
            mmcif::MmcifEntry::Item(item) if item.tag().starts_with("_entity_poly.") => {
                Some(item.tag())
            }
            _ => None,
        })
        .collect::<Vec<_>>();
    let values = tags
        .iter()
        .enumerate()
        .map(|(i, tag)| {
            let a = category.find(tag).unwrap() + tag.len();
            let b = tags
                .get(i + 1)
                .map_or(category.len(), |next| category.find(next).unwrap());
            category[a..b].trim()
        })
        .collect::<Vec<_>>();
    let looped = format!(
        "{}loop_\n{}\n{}{}",
        &SOURCE[..start],
        tags.join("\n"),
        values.join("\n"),
        &SOURCE[end..]
    );
    let a = doc.interpret().unwrap();
    let b = mmcif::parse_str(&looped).unwrap().interpret().unwrap();
    assert_eq!(a.model().atom_count(), 3804);
    assert_eq!(a.topology().bond_count(), 3484);
    assert_eq!(a.topology().instance_count(), 382);
    assert!(a.topology().same_layout(b.topology()));
    assert_eq!(a.model().positions(), b.model().positions());
    assert_eq!(a.model().properties(), b.model().properties());
    let protein = crate::topology::AtomSelection::for_molecule_classes(
        &a.model().shared_topology(),
        [crate::topology::MoleculeClass::Protein],
    )
    .unwrap();
    assert_eq!(protein.indices().len(), 3312);
    assert_eq!(a.model().occupancies(), b.model().occupancies());
    assert_eq!(a.model().b_factors(), b.model().b_factors());
}

#[test]
fn mmcif_quoted_values_allow_embedded_delimiters() {
    for (literal, expected) in [
        ("'don't stop'", "don't stop"),
        ("\"a\"quoted\"word\"", "a\"quoted\"word"),
        ("'a'#b'", "a'#b"),
        ("'a''", "a'"),
        ("''", ""),
    ] {
        for separator in [" ", "\t", "\n"] {
            let source = format!("data_x\n_x.text {literal}{separator}_x.after yes\n");
            let document = mmcif::parse_str(&source).unwrap();
            let block = &document.blocks()[0];
            assert_eq!(block.item("_x.text").unwrap().text(), expected);
            assert_eq!(block.item("_x.after").unwrap().text(), "yes");
        }
    }
    let document =
        mmcif::parse_str("data_x\nloop_\n_x.text\n_x.id\n'don't stop' 1\n\"a\"b\" 2\n").unwrap();
    let table = document.blocks()[0].loop_with_tag("_x.text").unwrap();
    assert_eq!(table.value(0, "_x.text").unwrap().text(), "don't stop");
    assert_eq!(table.value(1, "_x.text").unwrap().text(), "a\"b");
}

#[test]
fn mmcif_loop_values_on_one_long_line_parse_in_order() {
    let rows = 20_000;
    let values = (0..rows).map(|row| row.to_string()).collect::<Vec<_>>();
    let source = format!("data_x\nloop_\n_test.value\n{} \n", values.join(" "));
    let document = mmcif::parse_str(&source).unwrap();
    let table = document.blocks()[0].loop_with_tag("_test.value").unwrap();
    assert_eq!(table.row_count(), rows);
    let parsed = (0..rows)
        .map(|row| {
            let value = table.value_at(row, 0).unwrap();
            assert_eq!(value.line(), 4);
            value.text()
        })
        .collect::<Vec<_>>();
    assert_eq!(parsed, values);
}

#[test]
fn mmcif_quotes_must_close_before_whitespace_or_end_of_line() {
    for literal in ["'value'junk", "\"value\"junk", "'value'#comment"] {
        let error = mmcif::parse_str(&format!("data_x\n_x.text {literal}\n")).unwrap_err();
        assert_eq!(error.location().line(), 2);
        assert!(error.message().contains("unterminated quoted value"));
    }
}

#[test]
fn mmcif_semicolon_text_preserves_opening_content_and_whitespace() {
    // CIF 1.1 syntax paragraphs 17–20: the closing delimiter's preceding
    // newline is excluded; other text (including initial blank lines) remains.
    for (field, expected) in [
        (";first\n;", "first"),
        (";first\nsecond\n;", "first\nsecond"),
        (";\n;", ""),
        (";\n\n  body  \n\n;", "\n\n  body  \n"),
        ("; foo\n  bar\n;", " foo\n  bar"),
        (";?\n;", "?"),
        (";.\n;", "."),
        (
            ";# literal\n ; not a delimiter\n;",
            "# literal\n ; not a delimiter",
        ),
    ] {
        let input = format!("data_x\n_x.text\n{field}\n_x.after yes\n");
        for newline in ["\n", "\r\n", "\r"] {
            let document = mmcif::parse_str(&input.replace('\n', newline)).unwrap();
            let value = document.blocks()[0].item("_x.text").unwrap();
            assert_eq!(value.text(), expected);
            assert_eq!(value.line(), 3);
            assert!(
                !value.is_missing(),
                "delimited placeholders are literal text"
            );
            assert_eq!(document.blocks()[0].item("_x.after").unwrap().text(), "yes");
        }
    }
}

#[test]
fn mmcif_semicolon_closing_line_preserves_following_tokens() {
    let doc = mmcif::parse_str("data_x\n_x.text\n;value\n; _x.after yes\n").unwrap();
    assert_eq!(doc.blocks()[0].item("_x.text").unwrap().text(), "value");
    assert_eq!(doc.blocks()[0].item("_x.after").unwrap().line(), 4);
    let error = mmcif::parse_str("data_x\n_x.text\n;value\n;not_delimited\n").unwrap_err();
    assert_eq!(error.location().line(), 4);
    assert!(error.message().contains("closing semicolon"));
    let error = mmcif::parse_str("data_x\n_x.text\n;value\nno end").unwrap_err();
    assert_eq!(error.location().line(), 3);
    assert!(error.message().contains("unterminated semicolon"));
}

#[test]
fn mmcif_semicolon_size_limit_counts_opening_line_and_newlines() {
    let options = MmcifParseOptions {
        max_token_bytes: 8,
        ..Default::default()
    };
    for field in [";12345678\n;", ";1234\n678\n;"] {
        let doc =
            mmcif::parse_str_with_options(&format!("data_x\n_x.v\n{field}\n"), options).unwrap();
        assert_eq!(doc.blocks()[0].item("_x.v").unwrap().text().len(), 8);
    }
    for field in [";123456789\n;", ";1234\n6789\n;", ";12345678\n\n;"] {
        let error = mmcif::parse_str_with_options(&format!("data_x\n_x.v\n{field}\n"), options)
            .unwrap_err();
        assert_eq!(error.location().line(), 3);
        assert!(error.message().contains("token limit"));
    }
}

#[test]
fn mmcif_failures_carry_structured_kinds_and_original_lines() {
    use crate::diagnostics::{DiagnosticKind, SourceLocation};
    for (input, kind, line) in [
        ("loop_\n", DiagnosticKind::Syntax, 1),
        ("data_x\n_a.b 'open\n", DiagnosticKind::Syntax, 2),
        ("data_x\r\n_a.b 'open\r\n", DiagnosticKind::Syntax, 2),
        ("data_x\r_a.b 1\r_a.b 2\r", DiagnosticKind::Inconsistent, 3),
        ("data_x\nloop_\n_a.b\n_a.c\n1\n", DiagnosticKind::Syntax, 5),
        ("data_x\n_a.b\n;text\n", DiagnosticKind::Syntax, 3),
    ] {
        let error = mmcif::parse_str(input).unwrap_err();
        assert_eq!(
            (error.kind(), error.location()),
            (kind, SourceLocation::at_line(line)),
            "{input:?}: {error}"
        );
    }
    let limited = MmcifParseOptions {
        max_tokens: 2,
        ..MmcifParseOptions::default()
    };
    let error = mmcif::parse_str_with_options("data_x\n_a.b 1\n_a.c 2\n", limited).unwrap_err();
    assert_eq!(
        (error.kind(), error.location()),
        (DiagnosticKind::Limit, SourceLocation::at_line(2))
    );

    let unknown_element = SINGLE_SITE.replace("O O HOH", "Xx O HOH");
    let error = mmcif::parse_str(&unknown_element)
        .unwrap()
        .interpret()
        .unwrap_err();
    assert_eq!(
        (error.kind(), error.location()),
        (
            DiagnosticKind::InvalidValue,
            Some(SourceLocation::at_line(19))
        )
    );
    let error = mmcif::parse_str("data_x\n_cell.length_a 1\n")
        .unwrap()
        .interpret()
        .unwrap_err();
    assert_eq!(
        (error.kind(), error.location()),
        (DiagnosticKind::Inconsistent, None)
    );
    let error = mmcif::parse_str(SINGLE_SITE)
        .unwrap()
        .interpret_ensemble_with_options(macromolecular::StructureEnsembleInterpretOptions {
            model_ids: Some(vec!["9".into()]),
            ..Default::default()
        })
        .unwrap_err();
    assert_eq!(
        (error.kind(), error.location()),
        (DiagnosticKind::Inconsistent, None)
    );
}

#[test]
fn mmcif_atom_site_row_limit_stops_at_the_first_excess_row() {
    use crate::diagnostics::{DiagnosticKind, SourceLocation};
    let options = |max_atom_site_rows| MmcifParseOptions {
        max_atom_site_rows,
        ..MmcifParseOptions::default()
    };
    // Each unterminated quote follows the excess row and is never tokenized.
    for (input, max_rows, line) in [
        ("data_x\nloop_\n_atom_site.id\n1\n2\n'open\n", 1, 5),
        (
            "data_x\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1 C\n2 'open\n",
            1,
            6,
        ),
        (
            "data_x\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1 C\n2\n",
            1,
            6,
        ),
        (
            "data_x\nloop_\n_x.id\n_atom_site.id\n1 1\n2 2 'open\n",
            1,
            6,
        ),
        ("data_x\nloop_\n_atom_site.id\n1\n'open\n", 0, 4),
    ] {
        let error = mmcif::parse_str_with_options(input, options(max_rows)).unwrap_err();
        assert_eq!(
            (error.kind(), error.location(), error.message()),
            (
                DiagnosticKind::Limit,
                SourceLocation::at_line(line),
                "atom-site row count exceeds configured limit"
            ),
            "{input:?}"
        );
    }
    for (input, max_rows, tag, rows) in [
        (
            "data_x\nloop_\n_atom_site.id\n1\nloop_\n_x.id\n1\n",
            1,
            "_atom_site.id",
            1,
        ),
        ("data_x\nloop_\n_atom_site.id\n", 0, "_atom_site.id", 0),
        ("data_x\nloop_\n_x.id\n1\n2\n", 1, "_x.id", 2),
    ] {
        let document = mmcif::parse_str_with_options(input, options(max_rows)).unwrap();
        let table = document.blocks()[0].loop_with_tag(tag).unwrap();
        assert_eq!(table.row_count(), rows, "{input:?}");
    }
}

#[test]
fn mmcif_readers_enforce_limits_while_reading_and_match_string_parsing() {
    use crate::diagnostics::{DiagnosticKind, SourceLocation};
    assert_eq!(
        mmcif::parse_reader(SINGLE_SITE.as_bytes()).unwrap(),
        mmcif::parse_str(SINGLE_SITE).unwrap()
    );
    let limited = MmcifParseOptions {
        max_input_bytes: 8,
        ..MmcifParseOptions::default()
    };
    let error = mmcif::parse_reader_with_options(SINGLE_SITE.as_bytes(), limited).unwrap_err();
    assert_eq!(
        (error.kind(), error.location()),
        (DiagnosticKind::Limit, SourceLocation::at_line(1))
    );
    let error = mmcif::parse_reader([b'd', b'a', 0xff].as_slice()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Io);
}

#[test]
fn mmcif_values_borrow_the_source_and_columns_resolve_once() {
    let doc = mmcif::parse_str(SINGLE_SITE).unwrap();
    let table = doc.blocks()[0].loop_with_tag("_ATOM_SITE.CARTN_Y").unwrap();
    let column = table.column_index("_atom_site.Cartn_y").unwrap();
    assert_eq!(column, 6);
    let value = table.value_at(0, column).unwrap();
    assert_eq!((value.text(), value.line()), ("2.0", 19));
    assert_eq!(table.value(0, "_atom_site.cartn_y"), Some(value));
    assert_eq!(table.value_at(1, column), None);
    assert_eq!(table.value_at(0, 8), None);
    assert_eq!(
        table
            .row(0)
            .unwrap()
            .map(|value| value.text())
            .collect::<Vec<_>>(),
        ["O", "O", "HOH", "A", "1", "1.0", "2.0", "3.0"]
    );
}
