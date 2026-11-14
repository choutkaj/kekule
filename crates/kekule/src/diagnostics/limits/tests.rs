use super::*;
use crate::diagnostics::{DiagnosticKind, Format, FormatError};

#[test]
fn limits_report_resource_bound_and_first_excess() {
    assert_eq!(LimitExceeded::check(Resource::Atoms, 3, 3), Ok(()));
    let exceeded = LimitExceeded::check(Resource::Atoms, 3, 5).unwrap_err();
    assert_eq!(
        (exceeded.resource(), exceeded.limit(), exceeded.observed()),
        (Resource::Atoms, 3, 5)
    );
    assert_eq!(
        exceeded.to_string(),
        "Exceeded the atom limit: 5 atoms of 3."
    );
    let work = LimitExceeded::new(Resource::Work, 69_248_000, 69_400_898);
    assert_eq!(
        format!("Ring perception exceeded its {}.", work.describe()),
        "Ring perception exceeded its work limit: 69,400,898 steps of 69,248,000."
    );
    assert_eq!(
        LimitExceeded::first_beyond(Resource::Records, usize::MAX).observed(),
        usize::MAX
    );
    let error = FormatError::parse(Format::Sdf, DiagnosticKind::Limit, "too many records")
        .exceeding(LimitExceeded::first_beyond(Resource::Records, 2));
    assert_eq!(
        error.limit_exceeded(),
        Some(LimitExceeded::new(Resource::Records, 2, 3))
    );
}
