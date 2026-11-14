use std::fmt;

/// A count written with thousands separators, such as `69,248,000`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Count(pub(crate) usize);

impl fmt::Display for Count {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let digits = self.0.to_string();
        let lead = digits.len() % 3;
        for (index, digit) in digits.char_indices() {
            if index != 0 && index % 3 == lead {
                formatter.write_str(",")?;
            }
            write!(formatter, "{digit}")?;
        }
        Ok(())
    }
}

/// A count with its noun, such as `1 atom` or `1,024 atoms`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Counted<'a> {
    count: usize,
    singular: &'a str,
    plural: &'a str,
}

/// `count` with `singular` or `plural`, such as `1 atom` or `1,024 atoms`.
pub(crate) const fn counted<'a>(count: usize, singular: &'a str, plural: &'a str) -> Counted<'a> {
    Counted {
        count,
        singular,
        plural,
    }
}

impl fmt::Display for Counted<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let noun = if self.count == 1 {
            self.singular
        } else {
            self.plural
        };
        write!(formatter, "{} {noun}", Count(self.count))
    }
}

/// Items joined as a list in a sentence, such as `"CG", "CD" and "NE"`.
pub(crate) fn listed<T: fmt::Display>(items: &[T]) -> String {
    match items {
        [] => String::new(),
        [only] => only.to_string(),
        [rest @ .., last] => {
            let rest = rest.iter().map(ToString::to_string).collect::<Vec<_>>();
            format!("{} and {last}", rest.join(", "))
        }
    }
}

/// Writes a value's display text with its first letter in upper case, for a
/// noun phrase such as an identifier (`atom 3`) that starts a sentence.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Capitalized<T>(pub(crate) T);

impl<T: fmt::Display> fmt::Display for Capitalized<T> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = self.0.to_string();
        let mut characters = text.chars();
        if let Some(first) = characters.next() {
            write!(formatter, "{}", first.to_uppercase())?;
            formatter.write_str(characters.as_str())?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::{LimitExceeded, Resource};

    #[test]
    fn counted_nouns_agree_with_their_count() {
        assert_eq!(counted(1, "atom", "atoms").to_string(), "1 atom");
        assert_eq!(counted(0, "atom", "atoms").to_string(), "0 atoms");
        assert_eq!(
            counted(1_024, "entry", "entries").to_string(),
            "1,024 entries"
        );
        assert_eq!(
            LimitExceeded::first_beyond(Resource::Rows, 0)
                .describe()
                .to_string(),
            "row limit: 1 row of 0"
        );
        assert_eq!(
            LimitExceeded::first_beyond(Resource::MatrixEntries, 0)
                .describe()
                .to_string(),
            "matrix entry limit: 1 entry of 0"
        );
    }

    #[test]
    fn lists_join_with_and() {
        assert_eq!(listed::<&str>(&[]), "");
        assert_eq!(listed(&["a"]), "a");
        assert_eq!(listed(&["a", "b"]), "a and b");
        assert_eq!(listed(&["a", "b", "c"]), "a, b and c");
    }

    #[test]
    fn counts_group_thousands() {
        for (count, text) in [
            (0, "0"),
            (999, "999"),
            (1_000, "1,000"),
            (12_345, "12,345"),
            (123_456, "123,456"),
            (1_234_567, "1,234,567"),
        ] {
            assert_eq!(Count(count).to_string(), text);
        }
    }
}
