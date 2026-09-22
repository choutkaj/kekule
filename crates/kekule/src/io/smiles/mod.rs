mod canonical;
mod cx;
mod emit;
mod interpret;
mod parse;
mod write;

pub use canonical::write_canonical_smiles;
pub use interpret::{
    SmilesAtomMapping, SmilesBondMapping, SmilesComponentCountError, SmilesComponentInterpretation,
    SmilesInterpretError, SmilesInterpretation, SmilesInterpretationReport,
};
pub use parse::{
    parse_smiles_document, parse_smiles_document_with_options, SmilesDocument, SmilesDocumentToken,
    SmilesDocumentTokenKind, SmilesParseError, SmilesParseOptions,
};
pub use write::{write_isomeric_smiles, write_smiles};

pub(crate) fn write_topology(
    topology: &crate::topology::Topology,
    options: crate::smiles::SmilesWriteOptions,
) -> Result<String, crate::io::MolWriteError> {
    use crate::smiles::SmilesWriteMode;
    let mut parts = topology
        .molecules()
        .map(|m| match options.mode {
            SmilesWriteMode::Canonical => canonical::write_canonical_emission(m.molecule()),
            SmilesWriteMode::Isomeric => write::write_source_order_emission(
                m.molecule(),
                write::StereoWriteMode::Encode,
                write::CanonicalAtomStyle::StoredKekule,
            ),
            SmilesWriteMode::Ordinary => write::write_source_order_emission(
                m.molecule(),
                write::StereoWriteMode::Reject,
                write::CanonicalAtomStyle::Aromatic,
            ),
        })
        .collect::<Result<Vec<_>, _>>()?;
    if options.mode == SmilesWriteMode::Canonical {
        let mut keyed = parts
            .into_iter()
            .map(|p| Ok((p.render()?, p)))
            .collect::<Result<Vec<_>, crate::io::MolWriteError>>()?;
        keyed.sort_by(|a, b| a.0.cmp(&b.0));
        parts = keyed.into_iter().map(|(_, p)| p).collect();
    }
    emit::Emission::join(parts).render()
}
