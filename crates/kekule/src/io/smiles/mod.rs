mod canonical;
mod cx;
mod emit;
mod interpret;
mod parse;
mod write;

pub(crate) use canonical::write_canonical_smiles;
pub use cx::{SmilesCxField, SmilesCxFieldKind, SmilesCxStereoGroupKind};
pub use interpret::{
    SmilesAtomMapping, SmilesBondMapping, SmilesComponentCountError, SmilesComponentInterpretation,
    SmilesInterpretError, SmilesInterpretOptions, SmilesInterpretation, SmilesInterpretationReport,
};
pub(crate) use parse::{parse_smiles_document, parse_smiles_reader};
pub use parse::{
    SmilesDocument, SmilesDocumentToken, SmilesDocumentTokenKind, SmilesParseError,
    SmilesParseOptions, SmilesRecords,
};
pub(crate) use write::{write_isomeric_smiles, write_smiles};

pub(crate) fn write_topology(
    topology: &crate::topology::Topology,
    options: crate::smiles::SmilesWriteOptions,
) -> Result<String, crate::io::MoleculeWriteError> {
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
            .collect::<Result<Vec<_>, crate::io::MoleculeWriteError>>()?;
        keyed.sort_by(|a, b| a.0.cmp(&b.0));
        parts = keyed.into_iter().map(|(_, p)| p).collect();
    }
    emit::Emission::join(parts).render()
}
