use std::fmt;

use super::Capitalized;

/// An error followed by every error of its
/// [`source`](std::error::Error::source) chain, written as consecutive
/// sentences.
///
/// Every Kekule error displays as one or two sentences and does not repeat
/// its source, so the chain reads as one message, such as `Cannot build the
/// solvated structure. Atom 12 overlaps atom 3.` Errors from outside Kekule in
/// the chain, such as I/O errors, get an upper-case first letter and a
/// closing full stop when they lack them. A cause whose text equals the
/// previous one, as a transparent wrapper writes it, appears once. The
/// alternate form, `{:#}`, writes one error per line.
///
/// ```
/// use kekule::diagnostics::ErrorChain;
///
/// let error = kekule::smiles::to_molecules("C1CC").unwrap_err();
/// let message = ErrorChain::new(&error).to_string();
/// assert!(message.starts_with("Cannot parse the SMILES input"));
/// assert!(message.ends_with('.'));
/// ```
#[derive(Clone, Copy)]
pub struct ErrorChain<'a> {
    error: &'a (dyn std::error::Error + 'static),
}

impl<'a> ErrorChain<'a> {
    /// The chain that starts at `error`.
    pub fn new(error: &'a (dyn std::error::Error + 'static)) -> Self {
        Self { error }
    }
}

impl fmt::Debug for ErrorChain<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ErrorChain")
            .field("error", &self.error)
            .finish()
    }
}

impl fmt::Display for ErrorChain<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let separator = if formatter.alternate() { "\n" } else { " " };
        let mut previous: Option<String> = None;
        let mut next = Some(self.error);
        while let Some(error) = next {
            next = error.source();
            let text = error.to_string();
            let text = text.trim();
            if text.is_empty() || previous.as_deref() == Some(text) {
                continue;
            }
            if previous.is_some() {
                formatter.write_str(separator)?;
            }
            write!(formatter, "{}", Capitalized(text))?;
            if !text.ends_with(['.', '!', '?']) {
                formatter.write_str(".")?;
            }
            previous = Some(text.to_owned());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
