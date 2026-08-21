use super::structure_documents as raw;
use super::{MolfileDocument, MolfileInterpretError, MolfileInterpretation};

/// Interprets one Molfile record as exactly one connected small molecule.
///
/// The Molfile document remains a faithful syntax container. A CTAB containing
/// multiple disconnected graph components is rejected at the semantic
/// `SmallMolecule` boundary rather than leaking a disconnected `Molecule`.
pub fn interpret_molfile_document(
    document: &MolfileDocument,
) -> Result<MolfileInterpretation, MolfileInterpretError> {
    let interpretation = raw::interpret_molfile_document(document)?;
    interpretation
        .molecule()
        .as_molecule()
        .validate_connected()
        .map_err(|error| MolfileInterpretError {
            line: document
                .atom_records()
                .first()
                .map(|record| record.number())
                .unwrap_or(1),
            message: error.to_string(),
        })?;
    Ok(interpretation)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disconnected_molfile_is_not_published_as_one_small_molecule() {
        let input = "disconnected\nkekule\n\n  2  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
        let document = raw::parse_molfile_document(input).expect("parse CTAB");
        let error = interpret_molfile_document(&document).expect_err("reject disconnected CTAB");
        assert!(error.message().contains("molecule must be connected"));
    }
}
