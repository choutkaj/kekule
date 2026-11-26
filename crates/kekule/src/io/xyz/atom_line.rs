//! The atom lines of XYZ text, as the document parser and the trajectory
//! reader of `traj::io::xyz` both read them.
//!
//! An atom line holds an element token, three Cartesian coordinates, and
//! optionally further whitespace-separated columns, as extended XYZ writes
//! them. The token is an element symbol in any letter case or an atomic
//! number.

use crate::math::Point3;
use crate::molecule::Element;

/// An atom line split into its element token, its position, and the columns
/// after the coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct AtomLine<'a> {
    /// The element token as written.
    pub(crate) element: &'a str,
    /// The coordinates, each a finite number.
    pub(crate) position: Point3,
    /// The columns after the coordinates, verbatim from the first to the
    /// last non-whitespace character, or `None` when there are none.
    pub(crate) extra: Option<&'a str>,
}

/// Why a line is not an atom line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AtomLineError<'a> {
    /// The line lacks an element or one of three coordinates.
    MissingField,
    /// The coordinate is not a number.
    NotNumber(&'a str),
    /// The coordinate is a number that is not finite.
    NotFinite(&'a str),
}

/// Splits `text` into an atom line at Unicode whitespace.
pub(crate) fn parse_atom_line(text: &str) -> Result<AtomLine<'_>, AtomLineError<'_>> {
    let mut rest = text.trim_start();
    let mut field = || {
        let end = rest.find(char::is_whitespace).unwrap_or(rest.len());
        let (value, tail) = rest.split_at(end);
        rest = tail.trim_start();
        (!value.is_empty()).then_some(value)
    };
    let (Some(element), Some(x), Some(y), Some(z)) = (field(), field(), field(), field()) else {
        return Err(AtomLineError::MissingField);
    };
    let position = Point3::new(coordinate(x)?, coordinate(y)?, coordinate(z)?);
    let extra = rest.trim_end();
    Ok(AtomLine {
        element,
        position,
        extra: (!extra.is_empty()).then_some(extra),
    })
}

/// The finite number `value` holds.
fn coordinate(value: &str) -> Result<f64, AtomLineError<'_>> {
    let number = value
        .parse::<f64>()
        .map_err(|_| AtomLineError::NotNumber(value))?;
    if number.is_finite() {
        Ok(number)
    } else {
        Err(AtomLineError::NotFinite(value))
    }
}

/// The element `token` names: a symbol in any letter case, or an atomic
/// number.
pub(crate) fn element_of_token(token: &str) -> Option<Element> {
    if token.bytes().all(|byte| byte.is_ascii_digit()) {
        return token
            .parse::<u8>()
            .ok()
            .and_then(Element::from_atomic_number);
    }
    let mut characters = token.chars();
    let first = characters.next()?;
    let symbol = first
        .to_uppercase()
        .chain(characters.flat_map(char::to_lowercase))
        .collect::<String>();
    Element::from_symbol(&symbol)
}
