use std::fmt;
use std::ops::Range;

use crate::core::{AtomId, BondId};
use crate::small::model::SmallMolecule;

use super::smiles::{self, SmilesDocument, SmilesParseOptions};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesInterpretError {
    offset: usize,
    message: String,
}

impl SmilesInterpretError {
    pub const fn offset(&self) -> usize {
        self.offset
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for SmilesInterpretError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "SMILES interpretation error at {}: {}",
            self.offset, self.message
        )
    }
}

impl std::error::Error for SmilesInterpretError {}

/// A single-molecule accessor was used for a component-aware interpretation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmilesComponentCountError {
    actual: usize,
}

impl SmilesComponentCountError {
    pub const fn actual(self) -> usize {
        self.actual
    }
}

impl fmt::Display for SmilesComponentCountError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "single-molecule SMILES access requires exactly one component, found {}",
            self.actual
        )
    }
}

impl std::error::Error for SmilesComponentCountError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesAtomMapping {
    atom: AtomId,
    source_span: Range<usize>,
}

impl SmilesAtomMapping {
    pub const fn atom(&self) -> AtomId {
        self.atom
    }

    pub fn source_span(&self) -> Range<usize> {
        self.source_span.clone()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmilesBondMapping {
    bond: BondId,
    source_offset: usize,
}

impl SmilesBondMapping {
    pub const fn bond(&self) -> BondId {
        self.bond
    }

    pub const fn source_offset(&self) -> usize {
        self.source_offset
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SmilesInterpretationReport {
    atom_mappings: Vec<SmilesAtomMapping>,
    bond_mappings: Vec<SmilesBondMapping>,
}

impl SmilesInterpretationReport {
    pub fn atom_mappings(&self) -> &[SmilesAtomMapping] {
        &self.atom_mappings
    }

    pub fn bond_mappings(&self) -> &[SmilesBondMapping] {
        &self.bond_mappings
    }
}

/// Interpretation of one dot-delimited connected SMILES component.
#[derive(Debug, Clone, PartialEq)]
pub struct SmilesComponentInterpretation {
    source_span: Range<usize>,
    molecule: SmallMolecule,
    report: SmilesInterpretationReport,
}

impl SmilesComponentInterpretation {
    pub fn source_span(&self) -> Range<usize> {
        self.source_span.clone()
    }

    pub fn molecule(&self) -> &SmallMolecule {
        &self.molecule
    }

    pub fn report(&self) -> &SmilesInterpretationReport {
        &self.report
    }

    pub fn into_molecule(self) -> SmallMolecule {
        self.molecule
    }

    pub fn into_parts(self) -> (SmallMolecule, SmilesInterpretationReport) {
        (self.molecule, self.report)
    }
}

/// Canonical interpretation of one SMILES document as connected molecules.
#[derive(Debug, Clone, PartialEq)]
pub struct SmilesInterpretation {
    components: Vec<SmilesComponentInterpretation>,
}

impl SmilesInterpretation {
    pub fn components(&self) -> &[SmilesComponentInterpretation] {
        &self.components
    }

    pub fn molecules(&self) -> impl ExactSizeIterator<Item = &SmallMolecule> + DoubleEndedIterator {
        self.components
            .iter()
            .map(SmilesComponentInterpretation::molecule)
    }

    pub fn into_molecules(self) -> Vec<SmallMolecule> {
        self.components
            .into_iter()
            .map(SmilesComponentInterpretation::into_molecule)
            .collect()
    }

    /// Convenience access for callers that require exactly one component.
    ///
    /// Prefer [`Self::components`] for general SMILES input. This method keeps
    /// the pre-component API convenient for callers whose input contract is
    /// already single-molecule and fails loudly rather than discarding data.
    pub fn molecule(&self) -> Result<&SmallMolecule, SmilesComponentCountError> {
        Ok(self.single_component()?.molecule())
    }

    /// Convenience report access for an interpretation known to contain one component.
    pub fn report(&self) -> Result<&SmilesInterpretationReport, SmilesComponentCountError> {
        Ok(self.single_component()?.report())
    }

    /// Consumes an interpretation known to contain exactly one component.
    pub fn into_molecule(self) -> Result<SmallMolecule, SmilesComponentCountError> {
        Ok(self.into_single_component()?.into_molecule())
    }

    /// Consumes an interpretation known to contain exactly one component and its report.
    pub fn into_parts(
        self,
    ) -> Result<(SmallMolecule, SmilesInterpretationReport), SmilesComponentCountError> {
        Ok(self.into_single_component()?.into_parts())
    }

    fn single_component(
        &self,
    ) -> Result<&SmilesComponentInterpretation, SmilesComponentCountError> {
        match self.components.as_slice() {
            [component] => Ok(component),
            components => Err(SmilesComponentCountError {
                actual: components.len(),
            }),
        }
    }

    fn into_single_component(
        mut self,
    ) -> Result<SmilesComponentInterpretation, SmilesComponentCountError> {
        if self.components.len() != 1 {
            return Err(SmilesComponentCountError {
                actual: self.components.len(),
            });
        }
        Ok(self
            .components
            .pop()
            .expect("length was checked to contain one SMILES component"))
    }
}

/// Interprets each dot-delimited SMILES component independently.
///
/// Parsing remains record-level: [`SmilesDocument`] preserves the complete
/// source and component separators. Interpretation turns each syntactic
/// component into one connected [`SmallMolecule`] with component-local atom and
/// bond identifiers while retaining mappings to the original source offsets.
pub fn interpret_smiles_document(
    document: &SmilesDocument,
) -> Result<SmilesInterpretation, SmilesInterpretError> {
    let mut components = Vec::with_capacity(document.component_token_ranges().len());
    for token_range in document.component_token_ranges() {
        let source_span = component_source_span(document, token_range.clone())?;
        let source =
            document
                .source()
                .get(source_span.clone())
                .ok_or_else(|| SmilesInterpretError {
                    offset: source_span.start,
                    message: "component source span is outside the SMILES document".to_owned(),
                })?;

        // The complete document has already passed the caller's resource policy.
        // Reparse the isolated component only to reuse the mature single-component
        // semantic interpreter without imposing a second, stricter default limit.
        let local_document = smiles::parse_smiles_document_with_options(
            source,
            SmilesParseOptions {
                max_input_bytes: source.len(),
                max_atoms: u32::MAX as usize,
                max_bonds: u32::MAX as usize,
            },
        )
        .map_err(|error| SmilesInterpretError {
            offset: source_span.start.saturating_add(error.offset()),
            message: error.message().to_owned(),
        })?;
        let local = smiles::interpret_smiles_document(&local_document).map_err(|error| {
            SmilesInterpretError {
                offset: source_span.start.saturating_add(error.offset()),
                message: error.message().to_owned(),
            }
        })?;
        let (molecule, report) = local.into_parts();
        molecule
            .graph()
            .validate_connected()
            .map_err(|error| SmilesInterpretError {
                offset: source_span.start,
                message: error.to_string(),
            })?;

        let atom_mappings = report
            .atom_mappings()
            .iter()
            .map(|mapping| {
                let local_span = mapping.source_span();
                SmilesAtomMapping {
                    atom: mapping.atom(),
                    source_span: source_span.start.saturating_add(local_span.start)
                        ..source_span.start.saturating_add(local_span.end),
                }
            })
            .collect();
        let bond_mappings = report
            .bond_mappings()
            .iter()
            .map(|mapping| SmilesBondMapping {
                bond: mapping.bond(),
                source_offset: source_span.start.saturating_add(mapping.source_offset()),
            })
            .collect();

        components.push(SmilesComponentInterpretation {
            source_span,
            molecule,
            report: SmilesInterpretationReport {
                atom_mappings,
                bond_mappings,
            },
        });
    }
    Ok(SmilesInterpretation { components })
}

fn component_source_span(
    document: &SmilesDocument,
    token_range: Range<usize>,
) -> Result<Range<usize>, SmilesInterpretError> {
    if token_range.is_empty() {
        return Err(SmilesInterpretError {
            offset: document.source().len(),
            message: "empty SMILES component".to_owned(),
        });
    }
    let first = document
        .tokens()
        .get(token_range.start)
        .ok_or_else(|| SmilesInterpretError {
            offset: document.source().len(),
            message: "component token range starts outside the SMILES document".to_owned(),
        })?;
    let last = document
        .tokens()
        .get(token_range.end.saturating_sub(1))
        .ok_or_else(|| SmilesInterpretError {
            offset: document.source().len(),
            message: "component token range ends outside the SMILES document".to_owned(),
        })?;
    Ok(first.span().start..last.span().end)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot_smiles_interprets_as_connected_molecules() {
        let document = smiles::parse_smiles_document("CC(=O)[O-].[Na+]").expect("valid salt");
        let interpretation = interpret_smiles_document(&document).expect("interpret salt");
        assert_eq!(interpretation.components().len(), 2);
        assert!(interpretation
            .molecules()
            .all(|molecule| molecule.graph().is_connected()));
        assert_eq!(interpretation.components()[0].source_span(), 0..10);
        assert_eq!(interpretation.components()[1].source_span(), 11..16);
    }

    #[test]
    fn component_mappings_retain_document_offsets() {
        let document = smiles::parse_smiles_document("C.[Na+]").expect("valid components");
        let interpretation = interpret_smiles_document(&document).expect("interpret components");
        assert_eq!(
            interpretation.components()[1].report().atom_mappings()[0].source_span(),
            2..7
        );
    }

    #[test]
    fn single_component_convenience_rejects_dot_smiles_without_panicking() {
        let document = smiles::parse_smiles_document("C.O").expect("valid components");
        let error = interpret_smiles_document(&document)
            .expect("interpret components")
            .into_molecule()
            .expect_err("single-component access must reject two components");
        assert_eq!(error.actual(), 2);
    }
}
