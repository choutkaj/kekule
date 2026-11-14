use super::*;

#[test]
fn incomplete_search_never_publishes_a_candidate() {
    let molecule = crate::tests::read_smiles("C1CCCCC1").expect("cyclohexane");
    let error = CanonicalOrder::with_limits(
        &molecule,
        CanonicalAtomStyle::Aromatic,
        1,
        MAX_CANONICAL_REFINEMENT_WORK,
    )
    .err()
    .expect("search is incomplete");
    assert_eq!(error.kind(), crate::diagnostics::DiagnosticKind::Limit);
    assert_eq!(
        error.limit_exceeded(),
        Some(LimitExceeded::new(Resource::SearchStates, 1, 2))
    );
    let error = CanonicalOrder::with_limits(
        &molecule,
        CanonicalAtomStyle::Aromatic,
        MAX_CANONICAL_SEARCH_STATES,
        1,
    )
    .err()
    .expect("refinement is incomplete");
    assert_eq!(error.kind(), crate::diagnostics::DiagnosticKind::Limit);
    assert_eq!(
        error
            .limit_exceeded()
            .map(|limit| (limit.resource(), limit.limit())),
        Some((Resource::Work, 1))
    );
}
