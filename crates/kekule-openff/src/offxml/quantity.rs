//! Bounded, non-evaluating parser for OFFXML multiplicative unit expressions.
use super::{error, invalid};
use crate::Result;
use kekule::units::*;

struct Value {
    magnitude: f64,
    unit: Unit,
}
struct Parser<'a> {
    source: &'a str,
    position: usize,
}
impl Parser<'_> {
    fn whitespace(&mut self) {
        while self.peek().is_some_and(|c| c.is_ascii_whitespace()) {
            self.position += 1;
        }
    }
    fn peek(&self) -> Option<u8> {
        self.source.as_bytes().get(self.position).copied()
    }
    fn take(&mut self, token: &str) -> bool {
        self.whitespace();
        if self.source[self.position..].starts_with(token) {
            self.position += token.len();
            true
        } else {
            false
        }
    }
    fn expression(&mut self, depth: usize) -> Result<Value> {
        let mut value = self.factor(depth)?;
        loop {
            let divide = if self.take("*") {
                false
            } else if self.take("/") {
                true
            } else {
                break;
            };
            let right = self.factor(depth)?;
            if divide {
                if right.magnitude == 0.0 {
                    return Err(error("division by zero in quantity"));
                }
                value.magnitude /= right.magnitude;
                value.unit = value.unit.try_div(right.unit).map_err(invalid)?;
            } else {
                value.magnitude *= right.magnitude;
                value.unit = value.unit.try_mul(right.unit).map_err(invalid)?;
            }
            finite(value.magnitude)?;
        }
        Ok(value)
    }
    fn factor(&mut self, depth: usize) -> Result<Value> {
        if depth > 32 {
            return Err(error("quantity nesting limit exceeded (32)"));
        }
        self.whitespace();
        // Unary signs bind less tightly than powers, as in the reference parser.
        if self.take("+") {
            return self.factor(depth + 1);
        }
        if self.take("-") {
            let mut value = self.factor(depth + 1)?;
            value.magnitude = -value.magnitude;
            return Ok(value);
        }
        let mut value = if self.take("(") {
            let value = self.expression(depth + 1)?;
            if !self.take(")") {
                return Err(error("missing closing parenthesis in quantity"));
            }
            value
        } else if self.peek().is_some_and(|c| c.is_ascii_digit() || c == b'.') {
            let start = self.position;
            while self.peek().is_some_and(|c| c.is_ascii_digit() || c == b'.') {
                self.position += 1;
            }
            if self.peek().is_some_and(|c| c == b'e' || c == b'E') {
                self.position += 1;
                if self.peek().is_some_and(|c| c == b'+' || c == b'-') {
                    self.position += 1;
                }
                while self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    self.position += 1;
                }
            }
            Value {
                magnitude: super::number(&self.source[start..self.position])?,
                unit: DIMENSIONLESS,
            }
        } else {
            let start = self.position;
            while self
                .peek()
                .is_some_and(|c| c.is_ascii_alphabetic() || c == b'_')
            {
                self.position += 1;
            }
            Value {
                magnitude: 1.0,
                unit: named_unit(&self.source[start..self.position])?,
            }
        };
        if self.take("**") {
            let parenthesized = self.take("(");
            self.whitespace();
            let start = self.position;
            if self.peek().is_some_and(|c| c == b'+' || c == b'-') {
                self.position += 1;
            }
            while self.peek().is_some_and(|c| c.is_ascii_digit()) {
                self.position += 1;
            }
            let power = self.source[start..self.position]
                .parse::<i32>()
                .map_err(invalid)?;
            if !(-32..=32).contains(&power) {
                return Err(error("quantity exponent limit exceeded (32)"));
            }
            if parenthesized && !self.take(")") {
                return Err(error("invalid parenthesized quantity exponent"));
            }
            value.magnitude = value.magnitude.powi(power);
            value.unit = value.unit.try_powi(power).map_err(invalid)?;
            finite(value.magnitude)?;
        }
        Ok(value)
    }
}

fn finite(value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(error("nonfinite quantity"))
    }
}

fn named_unit(name: &str) -> Result<Unit> {
    Ok(match name {
        "angstrom" | "angstroms" => ANGSTROM,
        "nanometer" | "nanometers" | "nm" => NANOMETER,
        "degree" | "degrees" => DEGREE,
        "radian" | "radians" => RADIAN,
        "kilocalorie_per_mole" | "kilocalories_per_mole" => KILOCALORIE_PER_MOLE,
        "kilojoule_per_mole" | "kilojoules_per_mole" => KILOJOULE_PER_MOLE,
        "kilocalorie" | "kilocalories" | "kcal" => KILOCALORIE,
        "kilojoule" | "kilojoules" | "kJ" => KILOJOULE,
        "mole" | "moles" | "mol" => MOLE,
        "elementary_charge" | "elementary_charges" => ELEMENTARY_CHARGE,
        "dimensionless" => DIMENSIONLESS,
        _ => {
            return Err(error(format!(
                "unsupported unit or quantity token {name:?}"
            )))
        }
    })
}

pub(crate) fn parse_quantity(source: &str, target: Unit) -> Result<Quantity<f64>> {
    if source.len() > 4096 || !source.is_ascii() {
        return Err(error("quantity must be ASCII and at most 4096 bytes"));
    }
    let mut parser = Parser {
        source,
        position: 0,
    };
    let value = parser.expression(0)?;
    parser.whitespace();
    if parser.position != source.len() {
        return Err(error(format!(
            "unexpected quantity text at byte {}",
            parser.position
        )));
    }
    let converted = Quantity::new(value.magnitude, value.unit)
        .to_unit(target)
        .map_err(invalid)?;
    finite(*converted.value())?;
    Ok(converted)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equivalent_unit_spellings_and_parentheses() {
        for text in [
            "2*kilocalorie_per_mole*angstrom**-2",
            "2 * kilocalorie / (mole * angstrom**2)",
            "2e0*kcal/mol/angstroms**2",
            "2.*(kilocalories/moles)*angstrom**(-2)",
            "+2 * kilocalorie_per_mole / angstrom ** 2",
        ] {
            let value = parse_quantity(text, CANONICAL_FORCE_CONSTANT_UNIT).unwrap();
            assert!((*value.value() - 836.8).abs() < 1e-10, "{text}");
        }
        assert_eq!(
            *parse_quantity("-2**2 * elementary_charge", ELEMENTARY_CHARGE)
                .unwrap()
                .value(),
            -4.0
        );
        assert_eq!(
            *parse_quantity("(-2)**2 * elementary_charge", ELEMENTARY_CHARGE)
                .unwrap()
                .value(),
            4.0
        );
    }

    #[test]
    fn invalid_or_unbounded_expressions_are_rejected() {
        for text in [
            "",
            "1 nm",
            "1*nm+2*nm",
            "1*nm/0",
            "NaN*nm",
            "1e999*nm",
            "1*nm**.5",
            "1*nm**999999999999",
            "1*(nm",
            "1*nm)",
            "1*nm**2**3",
            "1*kelvin",
            "1*__import__('os')",
            "1*nm # comment",
            "1*nm***2",
            "1*nm**33",
            "1*nm/1e-999",
            "1*nm*",
            "1*Å",
        ] {
            assert!(parse_quantity(text, NANOMETER).is_err(), "{text}");
        }
        assert!(parse_quantity(
            &format!("{}nm{}", "(".repeat(34), ")".repeat(34)),
            NANOMETER
        )
        .is_err());
        assert!(parse_quantity(&" ".repeat(4097), NANOMETER).is_err());
        assert!(parse_quantity("1*kJ/mol", NANOMETER).is_err());
    }
}
