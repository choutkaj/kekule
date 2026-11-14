use super::*;

#[test]
fn drafts_return_with_their_error_displayed_transparently() {
    let error = DraftError::new(std::fmt::Error, vec![1, 2]);
    assert_eq!(error.draft(), &[1, 2]);
    assert_eq!(error.to_string(), std::fmt::Error.to_string());
    assert!(std::error::Error::source(&error).is_none());
    let (cause, draft) = error.into_parts();
    assert_eq!((cause, draft), (std::fmt::Error, vec![1, 2]));
}
