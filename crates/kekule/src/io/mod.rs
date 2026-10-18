//! Format parsing, interpretation, and writing. Items are published through
//! the format modules of the crate root (`smiles`, `molfile`, `sdf`, `mmcif`,
//! `pdb`, `components`, `macromolecular`); this module only organizes them.

mod components;
mod element_symbol;
mod lines;
mod mmcif_category;
mod mmcif_cell;
mod mmcif_document;
mod mmcif_write;
mod molfile_document;
mod molfile_write;
mod native;
mod pdb_document;
mod sdf_document;
mod sdf_ensemble;
mod sdf_write;
pub(crate) mod smiles;
mod staged_coordinates;
mod structure_connectivity;
mod structure_documents;
mod structure_interpret;
mod v2000;
mod v3000;
pub(crate) mod write_error;

pub use components::{
    Component, ComponentAtom, ComponentBond, ComponentBondOrder, ComponentLibrary, ComponentLinking,
};
pub(crate) use mmcif_document::{parse_mmcif_reader, parse_mmcif_str};
pub use mmcif_document::{
    MmcifBlock, MmcifDocument, MmcifEntry, MmcifItem, MmcifLoopTable, MmcifParseOptions, MmcifValue,
};
pub(crate) use mmcif_write::{write_mmcif, write_mmcif_to};
pub use mmcif_write::{
    MmcifAtomField, MmcifBlockSource, MmcifCellReason, MmcifEntityClassifications,
    MmcifMemberField, MmcifWriteIssue, MmcifWriteOptions,
};
pub(crate) use molfile_document::{parse_molfile_document, parse_molfile_reader};
pub use molfile_document::{
    MolfileAtom, MolfileAtomParity, MolfileBond, MolfileBondOrder, MolfileBondStereo,
    MolfileDocument, MolfileHeader, MolfileParseOptions, MolfileRadical, MolfileStereoGroup,
    MolfileStereoGroupKind, MolfileUnsupported, MolfileUnsupportedContent, MolfileVersion,
};
pub(crate) use molfile_write::{write_molfile, write_molfile_to};
pub use molfile_write::{MolfileSource, MolfileWriteOptions, MolfileWriteVersion};
pub(crate) use native::{read_native, write_native, write_native_to};
pub use native::{
    NativeKind, NativeObject, NativeReadOptions, NativeSource, NativeWriteOptions,
    DEFAULT_NATIVE_MAX_ATOMS, DEFAULT_NATIVE_MAX_BONDS, DEFAULT_NATIVE_MAX_INPUT_BYTES,
    DEFAULT_NATIVE_MAX_ITEMS, DEFAULT_NATIVE_MAX_STRING_BYTES, NATIVE_FORMAT_VERSION, NATIVE_MAGIC,
};
pub(crate) use pdb_document::{parse_pdb_reader, parse_pdb_str};
pub use pdb_document::{
    PdbAtom, PdbAtomKind, PdbConect, PdbCryst1, PdbDocument, PdbLink, PdbModel, PdbParseOptions,
    PdbPartner, PdbRecord, PdbSsbond, PdbTer,
};
pub(crate) use sdf_document::{parse_sdf_document, parse_sdf_reader};
pub use sdf_document::{
    SdfDataField, SdfDocument, SdfParseOptions, SdfRecord, SdfRecordInterpretation,
    SdfRecordInterpretationReport, SdfRecords,
};
pub(crate) use sdf_ensemble::sdf_ensemble;
pub(crate) use sdf_write::{write_sdf, write_sdf_to};
pub use sdf_write::{SdfRecordSource, SdfWriteOptions, SDF_ENERGY_FIELD, SDF_WEIGHT_FIELD};
pub(crate) use smiles::{
    parse_smiles_document, parse_smiles_reader, write_canonical_smiles, write_isomeric_smiles,
    write_smiles,
};
pub use smiles::{
    SmilesAtomMapping, SmilesBondMapping, SmilesComponentInterpretation, SmilesCxField,
    SmilesCxFieldKind, SmilesCxStereoGroupKind, SmilesDocument, SmilesDocumentToken,
    SmilesDocumentTokenKind, SmilesInterpretOptions, SmilesInterpretation,
    SmilesInterpretationReport, SmilesParseOptions, SmilesRecords, MAX_CANONICAL_GRAPH_SLOTS,
    MAX_CANONICAL_INPUT_COMPLEXITY, MAX_CANONICAL_PENDING_LABELS, MAX_CANONICAL_REFINEMENT_WORK,
    MAX_CANONICAL_SEARCH_STATES, MAX_DIRECTIONAL_COMBINATIONS, MAX_DIRECTIONAL_VISITS,
};
pub use structure_documents::{
    MolfileAtomMapping, MolfileBondMapping, MolfileInterpretOptions, MolfileInterpretation,
    MolfileInterpretationReport, MolfileInterpretationWarning,
};
pub use structure_interpret::{
    AltLocDecision, AltLocPolicy, AltLocPreference, AltLocResidue, AltLocSelection,
    AltLocSelectionReason, AtomProvenance, BondEvidence, BondProvenance, CellOmissionReason,
    ConflictingSelector, ConnectionResolutionReason, ConnectivityOptions, EntityKind,
    EntityMetadataPolicy, InstanceProvenance, ModelSelection, SourceResidueId,
    SourceResiduePosition, StructureEnsembleInterpretOptions, StructureEnsembleInterpretation,
    StructureInterpretIssue, StructureInterpretOptions, StructureInterpretation,
    StructureInterpretationReport,
};
pub use write_error::UnrepresentablePolicy;
