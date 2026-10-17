//! The recoverable-draft failure shared by builders and editors.
use std::fmt;

/// A failed recoverable publication that returns the unmodified draft with
/// the error.
///
/// `try_build` and `try_finish` of [`MoleculeEditor`](crate::molecule::MoleculeEditor),
/// [`TopologyBuilder`](crate::topology::TopologyBuilder),
/// [`TopologyEditor`](crate::topology::edit::TopologyEditor),
/// [`ModelBuilder`](crate::structure::ModelBuilder) and
/// [`ModelEditor`](crate::structure::ModelEditor) fail with it, so a caller can
/// repair the draft and retry. [`error`](Self::error) is also the
/// [`source`](std::error::Error::source); `Display` does not repeat it.
#[derive(Debug)]
pub struct DraftError<D, E> {
    parts: Box<(E, D)>,
}

impl<D, E> DraftError<D, E> {
    pub(crate) fn new(error: E, draft: D) -> Self {
        Self {
            parts: Box::new((error, draft)),
        }
    }

    /// Why the draft could not be published.
    pub fn error(&self) -> &E {
        &self.parts.0
    }

    /// The draft, exactly as it was before the attempt.
    pub fn draft(&self) -> &D {
        &self.parts.1
    }

    /// Returns the draft for repair, discarding the error.
    pub fn into_draft(self) -> D {
        self.parts.1
    }

    /// Returns the error and the draft.
    pub fn into_parts(self) -> (E, D) {
        *self.parts
    }
}

impl<D, E> fmt::Display for DraftError<D, E> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("draft could not be published; it is returned for repair")
    }
}

impl<D: fmt::Debug, E: std::error::Error + 'static> std::error::Error for DraftError<D, E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(self.error())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drafts_return_with_their_error_as_the_only_source() {
        let error = DraftError::new(std::fmt::Error, vec![1, 2]);
        assert_eq!(error.draft(), &[1, 2]);
        assert_eq!(
            error.to_string(),
            "draft could not be published; it is returned for repair"
        );
        assert_eq!(
            std::error::Error::source(&error).map(ToString::to_string),
            Some(std::fmt::Error.to_string())
        );
        let (cause, draft) = error.into_parts();
        assert_eq!((cause, draft), (std::fmt::Error, vec![1, 2]));
    }
}
