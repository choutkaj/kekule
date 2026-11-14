use super::*;

#[test]
fn permutation_materialization_and_sorting_require_budget() {
    // Twenty independently resonating groups would produce 2^20 forms.
    let mut work = Work {
        left: 10_000,
        limit: 10_000,
    };
    assert!(matches!(
        build_permutations(&[2; 20], &vec![vec![2]; 20], 1_000_000, &mut work),
        Err(ResonanceError::LimitExceeded(limit))
            if limit.resource() == Resource::Work && limit.limit() == 10_000
    ));
    // Construction alone fits, but construction plus sorting does not.
    let mut work = Work {
        left: 60,
        limit: 60,
    };
    assert!(matches!(
        build_permutations(&[2, 2], &[vec![1, 1], vec![1, 1]], 4, &mut work),
        Err(ResonanceError::LimitExceeded(limit))
            if limit.resource() == Resource::Work && limit.limit() == 60
    ));
}

#[test]
fn budgeted_permutations_preserve_reference_priority_order() {
    let mut work = Work {
        left: 1_000,
        limit: 1_000,
    };
    let permutations =
        build_permutations(&[2, 2], &[vec![1, 1], vec![1, 1]], 4, &mut work).unwrap();
    assert_eq!(
        permutations
            .into_iter()
            .map(|(_, chosen)| chosen)
            .collect::<Vec<_>>(),
        vec![vec![0, 0], vec![0, 1], vec![1, 0], vec![1, 1]]
    );
    assert!(work.left < work.limit);
    assert!(matches!(
        build_permutations(&[2], &[vec![2]], usize::MAX, &mut work),
        Err(ResonanceError::LimitExceeded(_))
    ));
}
