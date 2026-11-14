use super::*;

#[test]
fn only_whole_p1_unit_cells_are_stated_cells() {
    let box_ = [
        Some(20.0),
        Some(30.0),
        Some(40.0),
        Some(90.0),
        Some(90.0),
        Some(90.0),
    ];
    let placeholder = PLACEHOLDER_CELL.map(Some);
    let mut incomplete = box_;
    incomplete[2] = None;
    let space_group = |name: &str| CellOmissionReason::SpaceGroup {
        name: name.to_owned(),
    };
    for (parameters, space_groups, expected) in [
        (box_, &[][..], Ok(())),
        (box_, &["P 1"][..], Ok(())),
        (box_, &["p1", "P 1"][..], Ok(())),
        (
            box_,
            &["P 1", "P 21 21 21"][..],
            Err(space_group("P 21 21 21")),
        ),
        (box_, &["C 1 2 1"][..], Err(space_group("C 1 2 1"))),
        (
            placeholder,
            &["P 1"][..],
            Err(CellOmissionReason::Placeholder),
        ),
        (placeholder, &[][..], Err(CellOmissionReason::Placeholder)),
        (
            incomplete,
            &["P 1"][..],
            Err(CellOmissionReason::Incomplete),
        ),
    ] {
        let actual = stated_cell(parameters, space_groups)
            .unwrap()
            // Lengths pass through the canonical nanometre unit.
            .map(|cell| {
                cell.to_parameters()
                    .unwrap()
                    .lengths
                    .into_value()
                    .map(|value| (value * 1e9).round() / 1e9)
            });
        assert_eq!(
            actual,
            expected.map(|()| [20.0, 30.0, 40.0]),
            "{space_groups:?}"
        );
    }
    let flat = [
        Some(20.0),
        Some(30.0),
        Some(0.0),
        Some(90.0),
        Some(90.0),
        Some(90.0),
    ];
    assert!(stated_cell(flat, &[]).is_err());
}
