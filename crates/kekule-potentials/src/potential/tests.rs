use super::*;

#[test]
fn backend_failures_keep_their_cause_as_the_source() {
    let cause = std::io::Error::other("model server unavailable");
    let error = EvaluationError::backend(cause);
    assert_eq!(error.to_string(), "The potential backend failed.");
    let source = std::error::Error::source(&error).unwrap();
    assert_eq!(source.to_string(), "model server unavailable");
    assert!(source.downcast_ref::<std::io::Error>().is_some());
    assert_eq!(
        error,
        EvaluationError::backend(std::io::Error::other("model server unavailable"))
    );
    assert_ne!(
        error,
        EvaluationError::backend(std::io::Error::other("timeout"))
    );
}
