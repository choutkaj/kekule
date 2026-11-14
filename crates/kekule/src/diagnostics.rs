//! Shared diagnostics of format parsing, interpretation and writing.
//!
//! Every format reports failures as one [`FormatError`]: its [`FormatStage`],
//! a structured [`DiagnosticKind`] that callers can match on, and a
//! [`SourceLocation`] in the original input. Messages remain human-readable
//! detail, not an API contract. Writers share [`UnrepresentablePolicy`] for
//! state a format cannot represent, and structure writers return a
//! [`WriteReport`] of what that policy left out; readers share [`UnsupportedContent`]
//! for valid content Kekule does not support.
//!
//! Bounded operations, in formats and elsewhere, report an exhausted resource
//! limit as a [`LimitExceeded`] naming the counted [`Resource`].
//!
//! Long-running operations that accept a [`Cancellation`] stop when its check
//! turns true and report [`Cancelled`]; in a [`FormatError`] the kind is
//! [`DiagnosticKind::Cancelled`].
//!
//! # Example
//!
//! ```
//! use kekule::diagnostics::{DiagnosticKind, Format, FormatStage};
//!
//! let error = kekule::smiles::to_molecules("C1CC").unwrap_err();
//! assert_eq!(error.format(), Format::Smiles);
//! assert_eq!(error.stage(), FormatStage::Parse);
//! assert_eq!(error.kind(), DiagnosticKind::Syntax);
//! ```

mod cancellation;
mod chain;
mod format;
mod input;
mod limits;
mod wording;

pub use self::cancellation::{Cancellation, Cancelled};
pub use self::chain::ErrorChain;
pub use self::format::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, FormatStage, GeometryField,
    SourceLocation, UnsupportedContent,
};
pub use self::limits::{LimitExceeded, Resource};
pub use crate::io::{UnrepresentablePolicy, WriteOmission, WriteReport};

pub(crate) use self::input::{normalize_line_endings, read_bounded};
pub(crate) use self::wording::{counted, listed, Capitalized, Count};
