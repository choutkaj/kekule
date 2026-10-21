//! Canonical printing of selection expressions.

use std::fmt::{self, Write};

use super::ast::{
    DistanceUnit, Field, Flag, HostTermKind, Level, Literal, Node, NodeKind, ValueMatch,
};
use super::syntax::KEYWORDS;
use crate::topology::{MoleculeClass, ResidueClass};

/// Binding strength: an operand weaker than its position requires is
/// parenthesized.
fn strength(node: &Node) -> u8 {
    match node.kind {
        NodeKind::Or(_) => 1,
        NodeKind::And(_) => 2,
        _ => 3,
    }
}

pub(crate) fn node(formatter: &mut fmt::Formatter<'_>, node: &Node, required: u8) -> fmt::Result {
    let parenthesize = strength(node) < required;
    if parenthesize {
        formatter.write_char('(')?;
    }
    match &node.kind {
        NodeKind::All => formatter.write_str("all")?,
        NodeKind::None => formatter.write_str("none")?,
        NodeKind::Or(operands) => operator(formatter, operands, " or ", 2)?,
        NodeKind::And(operands) => operator(formatter, operands, " and ", 3)?,
        NodeKind::Not(inner) => {
            formatter.write_str("not ")?;
            self::node(formatter, inner, 3)?;
        }
        NodeKind::Within {
            distance,
            unit,
            reference,
        } => {
            write!(formatter, "within {distance}")?;
            match unit {
                DistanceUnit::Angstrom => {}
                DistanceUnit::Nanometer => formatter.write_str(" nm")?,
                DistanceUnit::Picometer => formatter.write_str(" pm")?,
                DistanceUnit::Bohr => formatter.write_str(" bohr")?,
            }
            formatter.write_str(" of ")?;
            self::node(formatter, &reference.root, 3)?;
        }
        NodeKind::Expand { level, inner } => {
            formatter.write_str(match level {
                Level::Residue => "byres ",
                Level::Chain => "bychain ",
                Level::Molecule => "bymolecule ",
            })?;
            self::node(formatter, inner, 3)?;
        }
        NodeKind::Bonded(inner) => {
            formatter.write_str("bonded to ")?;
            self::node(formatter, inner, 3)?;
        }
        NodeKind::Flag(flag) => formatter.write_str(flag_name(*flag))?,
        NodeKind::Match { field, values } => {
            field_name(formatter, field)?;
            for value in values {
                formatter.write_char(' ')?;
                match value {
                    ValueMatch::Exact(literal) => self::literal(formatter, literal)?,
                    ValueMatch::Range(low, high) => {
                        self::literal(formatter, low)?;
                        if is_nonnegative(low) && is_nonnegative(high) {
                            formatter.write_char('-')?;
                        } else {
                            formatter.write_str(" to ")?;
                        }
                        self::literal(formatter, high)?;
                    }
                }
            }
        }
        NodeKind::Compare { field, op, value } => {
            field_name(formatter, field)?;
            write!(formatter, " {} ", op.symbol())?;
            literal(formatter, value)?;
        }
        NodeKind::PropertyTrue(key) => write!(formatter, "property {key}")?,
        NodeKind::Smarts(pattern) => {
            formatter.write_str("smarts ")?;
            quoted(formatter, pattern)?;
        }
        NodeKind::Host(kind, text) => match kind {
            HostTermKind::Reference => write!(formatter, "#{text}")?,
            HostTermKind::Named => write!(formatter, "${text}")?,
            HostTermKind::Keyword => formatter.write_str(text)?,
        },
    }
    if parenthesize {
        formatter.write_char(')')?;
    }
    Ok(())
}

fn operator(
    formatter: &mut fmt::Formatter<'_>,
    operands: &[Node],
    separator: &str,
    required: u8,
) -> fmt::Result {
    for (index, operand) in operands.iter().enumerate() {
        if index > 0 {
            formatter.write_str(separator)?;
        }
        node(formatter, operand, required)?;
    }
    Ok(())
}

pub(crate) const fn flag_name(flag: Flag) -> &'static str {
    match flag {
        Flag::Aromatic => "aromatic",
        Flag::Ring => "ring",
        Flag::Hydrogen => "hydrogen",
        Flag::Heavy => "heavy",
        Flag::PolarHydrogen => "polar_hydrogen",
        Flag::Protein => "protein",
        Flag::Nucleic => "nucleic",
        Flag::Water => "water",
        Flag::Ion => "ion",
        Flag::Ligand => "ligand",
        Flag::Backbone => "backbone",
        Flag::Sidechain => "sidechain",
    }
}

pub(crate) fn field_name(formatter: &mut fmt::Formatter<'_>, field: &Field) -> fmt::Result {
    formatter.write_str(match field {
        Field::Index => "index",
        Field::Element => "element",
        Field::Name => "name",
        Field::LabelName => "label_name",
        Field::Charge => "charge",
        Field::Isotope => "isotope",
        Field::BFactor => "bfactor",
        Field::Occupancy => "occupancy",
        Field::Property(key) => return write!(formatter, "property {key}"),
        Field::ResName => "resname",
        Field::ResId => "resid",
        Field::LabelResId => "label_resid",
        Field::ResClass => "resclass",
        Field::Chain => "chain",
        Field::LabelChain => "label_chain",
        Field::MolClass => "molclass",
        Field::Instance => "instance",
        Field::Definition => "definition",
    })
}

fn is_nonnegative(literal: &Literal) -> bool {
    match literal {
        Literal::Integer(value) => *value >= 0,
        Literal::ResId { number, .. } => *number >= 0,
        Literal::Real { value, unit } => {
            unit.is_none() && *value >= 0.0 && !value.is_sign_negative()
        }
        _ => false,
    }
}

fn literal(formatter: &mut fmt::Formatter<'_>, literal: &Literal) -> fmt::Result {
    match literal {
        Literal::Integer(value) => write!(formatter, "{value}"),
        Literal::Real { value, unit } => {
            write!(formatter, "{value}")?;
            if let Some((_, spelling)) = unit {
                write!(formatter, " {spelling}")?;
            }
            Ok(())
        }
        Literal::Text(text) => word(formatter, text),
        Literal::Element(element) => formatter.write_str(element.symbol()),
        Literal::ResId {
            number,
            insertion_code,
        } => {
            write!(formatter, "{number}")?;
            if let Some(code) = insertion_code {
                formatter.write_str(code)?;
            }
            Ok(())
        }
        Literal::ResidueClass(class) => formatter.write_str(residue_class_name(*class)),
        Literal::MoleculeClass(class) => formatter.write_str(molecule_class_name(*class)),
        Literal::Property { text, number, unit } => {
            // An unquoted number reads back as a number; other text that
            // looks like one must stay quoted.
            if number.is_none() && text.parse::<f64>().is_ok_and(f64::is_finite) {
                quoted(formatter, text)?;
            } else {
                word(formatter, text)?;
            }
            if let Some((_, spelling)) = unit {
                write!(formatter, " {spelling}")?;
            }
            Ok(())
        }
    }
}

/// Writes `text` bare when it reads back as the same value, else quoted.
///
/// Lower-case identifiers are always quoted, so a value can never read back
/// as a keyword, including host keywords and keywords added later.
fn word(formatter: &mut fmt::Formatter<'_>, text: &str) -> fmt::Result {
    let identifier_like = text
        .chars()
        .all(|character| character.is_ascii_lowercase() || character == '_');
    let safe = !text.is_empty()
        && !identifier_like
        && !KEYWORDS.contains(&text)
        && !text.starts_with(['\'', '#', '$'])
        && text.chars().all(|character| {
            !character.is_whitespace()
                && !matches!(character, '(' | ')' | '"' | '<' | '>' | '=' | '!' | '\\')
        });
    if safe {
        formatter.write_str(text)
    } else {
        quoted(formatter, text)
    }
}

fn quoted(formatter: &mut fmt::Formatter<'_>, text: &str) -> fmt::Result {
    formatter.write_char('"')?;
    for character in text.chars() {
        if matches!(character, '"' | '\\') {
            formatter.write_char('\\')?;
        }
        formatter.write_char(character)?;
    }
    formatter.write_char('"')
}

pub(crate) const fn residue_class_name(class: ResidueClass) -> &'static str {
    match class {
        ResidueClass::AminoAcid => "amino_acid",
        ResidueClass::DnaNucleotide => "dna",
        ResidueClass::RnaNucleotide => "rna",
        ResidueClass::Carbohydrate => "carbohydrate",
        ResidueClass::Water => "water",
        ResidueClass::Ion => "ion",
        ResidueClass::Other => "other",
    }
}

pub(crate) const fn molecule_class_name(class: MoleculeClass) -> &'static str {
    match class {
        MoleculeClass::Protein => "protein",
        MoleculeClass::Dna => "dna",
        MoleculeClass::Rna => "rna",
        MoleculeClass::Carbohydrate => "carbohydrate",
        MoleculeClass::Water => "water",
        MoleculeClass::Ion => "ion",
        MoleculeClass::SmallMolecule => "small_molecule",
        MoleculeClass::Other => "other",
    }
}
