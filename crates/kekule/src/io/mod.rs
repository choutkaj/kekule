//! Format parsing, interpretation, and writing. Items are published through
//! the format modules of the crate root (`smiles`, `molfile`, `sdf`, `mmcif`,
//! `pdb`, `xyz`, `components`, `macromolecular`); this module only organizes them.
//!
//! Each format has one directory holding its parser, document types,
//! interpretation, writer, and entry points. Plumbing shared by several
//! formats stays at this level.

mod components;
mod coordinate_source;
mod definition_reuse;
mod element_symbol;
mod lines;
mod macromolecular;
pub(crate) mod mmcif;
pub(crate) mod molfile;
pub(crate) mod native;
mod normalization;
pub(crate) mod pdb;
pub(crate) mod sdf;
pub(crate) mod smiles;
mod source_stereo;
mod staged_coordinates;
pub(crate) mod write_policy;
pub(crate) mod xyz;

pub(crate) use components::dictionary_code;
pub use components::{
    Component, ComponentAtom, ComponentBond, ComponentBondOrder, ComponentLibrary, ComponentLinking,
};
pub(crate) use coordinate_source::AtomPositionSource;
pub(crate) use macromolecular::is_metal;
pub use macromolecular::{
    AltLocDecision, AltLocPolicy, AltLocPreference, AltLocResidue, AltLocSelection,
    AltLocSelectionReason, Assembly, AssemblyError, AssemblyGenerator, AssemblyOmissionReason,
    AssemblyOperation, AtomProvenance, BondEvidence, BondProvenance, CellOmissionReason,
    ConflictingSelector, ConnectionResolutionReason, ConnectivityOptions, CoordinationProvenance,
    EntityKind, EntityMetadataPolicy, EntityMonomer, InstanceProvenance, ModelSelection,
    PolymerEntity, PolymerEntityInstance, PolymerType, SequencePosition, SequenceSource,
    SourceResidueId, SourceResiduePosition, StructureInterpretIssue, StructureInterpretOptions,
    StructureInterpretation, StructureInterpretationReport, StructureSetInterpretOptions,
    StructureSetInterpretation,
};
pub(crate) use mmcif::{parse_mmcif_reader, parse_mmcif_str, write_mmcif, write_mmcif_to};
pub use mmcif::{
    MmcifAtomField, MmcifBlock, MmcifBlockSource, MmcifCellReason, MmcifDocument,
    MmcifEntityClassifications, MmcifEntry, MmcifGeometryField, MmcifItem, MmcifLoopTable,
    MmcifParseOptions, MmcifValue, MmcifWriteIssue, MmcifWriteOptions,
};
pub(crate) use molfile::{
    parse_molfile_document, parse_molfile_reader, write_molfile, write_molfile_to,
};
pub use molfile::{
    MolfileAtom, MolfileAtomMapping, MolfileAtomParity, MolfileBond, MolfileBondMapping,
    MolfileBondOrder, MolfileBondStereo, MolfileDocument, MolfileHeader, MolfileInterpretOptions,
    MolfileInterpretation, MolfileInterpretationReport, MolfileInterpretationWarning,
    MolfileParseOptions, MolfileRadical, MolfileSource, MolfileStereoGroup, MolfileStereoGroupKind,
    MolfileUnsupported, MolfileUnsupportedContent, MolfileVersion, MolfileWriteOptions,
    MolfileWriteVersion,
};
pub(crate) use native::{read_native, write_native, write_native_to};
pub use native::{
    EnsembleDescription, NativeKind, NativeObject, NativeReadOptions, NativeSource,
    NativeWriteOptions, StructureSetDescription, DEFAULT_NATIVE_MAX_ATOMS,
    DEFAULT_NATIVE_MAX_BONDS, DEFAULT_NATIVE_MAX_INPUT_BYTES, DEFAULT_NATIVE_MAX_ITEMS,
    DEFAULT_NATIVE_MAX_STRING_BYTES, NATIVE_FORMAT_VERSION, NATIVE_MAGIC,
};
pub(crate) use normalization::{
    canonicalize_molecule_for_publication, localize_source_aromatic_bonds, NormalizationError,
    NormalizationWarning,
};
#[cfg(test)]
pub(crate) use normalization::{
    NormalizationReport, SourceStereoNormalizationError, SourceStereoNormalizationIssue,
};
pub(crate) use pdb::{parse_pdb_reader, parse_pdb_str, write_pdb, write_pdb_to};
pub use pdb::{
    PdbAtom, PdbAtomField, PdbAtomKind, PdbCellReason, PdbChainIdPolicy, PdbConect, PdbCryst1,
    PdbDocument, PdbLink, PdbModel, PdbParseOptions, PdbPartner, PdbRecord, PdbSource, PdbSsbond,
    PdbTer, PdbWriteIssue, PdbWriteOptions,
};
pub(crate) use sdf::{
    parse_sdf_document, parse_sdf_reader, sdf_structure_set, write_sdf, write_sdf_to,
};
pub use sdf::{
    SdfDataField, SdfDocument, SdfParseOptions, SdfRecord, SdfRecordInterpretation,
    SdfRecordInterpretationReport, SdfRecordSource, SdfRecords, SdfWriteOptions,
    SDF_COORDINATE_FIELD, SDF_COORDINATE_UNIT_FIELD, SDF_STEP_FIELD,
};
pub use smiles::{
    SmilesAtomMapping, SmilesBondMapping, SmilesComponentInterpretation, SmilesCxField,
    SmilesCxFieldKind, SmilesCxStereoGroupKind, SmilesDocument, SmilesDocumentToken,
    SmilesDocumentTokenKind, SmilesInterpretOptions, SmilesInterpretation,
    SmilesInterpretationReport, SmilesParseOptions, SmilesRecords, MAX_CANONICAL_GRAPH_SLOTS,
    MAX_CANONICAL_INPUT_COMPLEXITY, MAX_CANONICAL_PENDING_LABELS, MAX_CANONICAL_REFINEMENT_WORK,
    MAX_CANONICAL_SEARCH_STATES, MAX_DIRECTIONAL_COMBINATIONS, MAX_DIRECTIONAL_VISITS,
};
pub(crate) use source_stereo::{
    has_specified_double_bond_stereo, molfile_double_bond_geometry_marks, normalize_source_stereo,
    project_molfile_stereo_bond_marks, remove_unkept_source_stereo, source_tetrahedral_carriers,
    SourceStereoBondMark, SourceStereoBondMarkKind,
};
pub use write_policy::{UnrepresentablePolicy, WriteOmission, WriteReport};
pub(crate) use xyz::{parse_xyz_reader, parse_xyz_str};
pub use xyz::{
    XyzAtom, XyzDocument, XyzFrame, XyzFrameSelection, XyzInferredBond, XyzInterpretIssue,
    XyzInterpretOptions, XyzInterpretation, XyzInterpretationReport, XyzParseOptions,
    XyzStructureSetInterpretOptions, XyzStructureSetInterpretation, DEFAULT_XYZ_MAX_ATOM_CONTACTS,
};
