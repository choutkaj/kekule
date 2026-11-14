use super::*;

fn entries(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> OperatorEntries {
    Some((matrix, vector))
}

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
/// A quarter turn about z.
const TURN: [[f64; 3]; 3] = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];

#[test]
fn oversized_expressions_are_refused_before_they_expand() {
    let operators = BTreeMap::from([("1".to_owned(), entries(IDENTITY, [0.0; 3]))]);
    let refused = |expression: &str, budget: usize| {
        let mut left = budget;
        let result = expression_operations(expression, &operators, &mut left);
        (result.err(), left)
    };
    // A product of ranges that would be 10^10 operations, and a list of a
    // hundred million identifiers, refused by their counts alone.
    let too_many = Some(AssemblyOmissionReason::TooManyOperations);
    assert_eq!(
        refused("(1-99999)(1-99999)", MAX_OPERATOR_USES),
        (too_many.clone(), MAX_OPERATOR_USES)
    );
    let list = "1-99999,".repeat(1_000);
    assert_eq!(
        refused(list.trim_end_matches(','), MAX_OPERATOR_USES).0,
        too_many
    );
    // One operation of 30,000 operators uses 30,000.
    let long = "(1)".repeat(30_000);
    assert_eq!(refused(&long, 29_999).0, too_many);
    assert_eq!(refused(&long, 30_000), (None, 0));
    // The budget is shared: what one expression takes, the next lacks.
    let mut budget = 5;
    assert!(expression_operations("(1,1,1)", &operators, &mut budget).is_ok());
    assert_eq!(budget, 2);
    assert_eq!(
        expression_operations("(1,1,1)", &operators, &mut budget).err(),
        too_many
    );
}

#[test]
fn expressions_name_lists_ranges_and_products() {
    let names = |expression| {
        expression_groups(expression).map(|groups| {
            groups
                .into_iter()
                .map(|group| group.ids().join(" "))
                .collect::<Vec<_>>()
        })
    };
    assert_eq!(names("1"), Some(vec!["1".to_owned()]));
    assert_eq!(names("1,2,5"), Some(vec!["1 2 5".to_owned()]));
    assert_eq!(names("(1-4)"), Some(vec!["1 2 3 4".to_owned()]));
    assert_eq!(
        names("(X0)(1-2)"),
        Some(vec!["X0".to_owned(), "1 2".to_owned()])
    );
    assert_eq!(
        names("(1,2)(3, 4)"),
        Some(vec!["1 2".to_owned(), "3 4".to_owned()])
    );
    for unreadable in ["", "(1", "1)", "(1)x", "1,,2", "(4-1)", "1 2"] {
        assert_eq!(names(unreadable), None, "{unreadable:?}");
    }
}

#[test]
fn products_apply_the_rightmost_operator_first() {
    let operators = BTreeMap::from([
        ("1".to_owned(), entries(IDENTITY, [0.0; 3])),
        ("T".to_owned(), entries(IDENTITY, [10.0, 0.0, 0.0])),
        ("R".to_owned(), entries(TURN, [0.0; 3])),
    ]);
    let operations = expression_operations("(R)(T)", &operators, &mut MAX_OPERATOR_USES.clone())
        .expect("readable");
    assert_eq!(operations.len(), 1);
    let operation = &operations[0];
    assert_eq!(operation.operators(), ["R", "T"]);
    assert_eq!(operation.name(), "R-T");
    // Translated along x first, then turned onto y.
    let moved = operation
        .translation()
        .into_unit(ANGSTROM)
        .expect("a length")
        .into_value();
    assert!((moved.x - 0.0).abs() < 1e-12 && (moved.y - 10.0).abs() < 1e-12);
    let product = expression_operations("(1,R)(1,T)", &operators, &mut MAX_OPERATOR_USES.clone())
        .expect("readable");
    let names = product
        .iter()
        .map(AssemblyOperation::name)
        .collect::<Vec<_>>();
    assert_eq!(names, ["1-1", "1-T", "R-1", "R-T"]);
    assert!(product[0].is_identity() && !product[1].is_identity());
}

#[test]
fn rounded_rotations_compose_without_drift() {
    // Within the tolerance, yet the product of 1,000 as stated would
    // stretch x by half.
    let stretched = [[1.0004, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let operators = BTreeMap::from([("S".to_owned(), entries(stretched, [0.0; 3]))]);
    let product = expression_operations(
        &"(S)".repeat(1000),
        &operators,
        &mut MAX_OPERATOR_USES.clone(),
    )
    .unwrap();
    assert_eq!(product.len(), 1);
    assert!(is_rotation(product[0].rotation, 1.0e-12));
    assert!(product[0].is_rigid());
    assert!(operator("I", entries(IDENTITY, [0.0; 3]))
        .unwrap()
        .is_identity());
}

#[test]
fn operators_that_are_not_rotations_are_refused() {
    let mirror = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let scaled = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for matrix in [mirror, scaled] {
        assert_eq!(
            operator("2", entries(matrix, [0.0; 3])),
            Err(AssemblyOmissionReason::NotRotation {
                operator: "2".to_owned()
            })
        );
    }
    assert_eq!(
        operator("3", None),
        Err(AssemblyOmissionReason::UnreadableOperator {
            operator: "3".to_owned()
        })
    );
    // Six decimals, as files state them.
    let rounded = [
        [0.309017, -0.809017, 0.5],
        [0.809017, 0.5, 0.309017],
        [-0.5, 0.309017, 0.809017],
    ];
    let operation = operator("4", entries(rounded, [1.0, 2.0, 3.0])).unwrap();
    assert!(is_rotation(operation.rotation, 1.0e-12));
    let operators = BTreeMap::from([("1".to_owned(), entries(IDENTITY, [0.0; 3]))]);
    assert_eq!(
        expression_operations("(1-2)", &operators, &mut MAX_OPERATOR_USES.clone()),
        Err(AssemblyOmissionReason::UnknownOperator {
            operator: "2".to_owned()
        })
    );
}
