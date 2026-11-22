//! Interpretation types shared by the macromolecular structure formats.
//!
//! Macromolecular formats interpret atom sites, alternate locations, coordinate
//! models, and connectivity through one pipeline. Its options, reports, issues,
//! and source provenance live here, and only here; the format modules [`mmcif`](crate::mmcif)
//! and [`pdb`](crate::pdb) take and return them.

pub use crate::io::{
    AltLocDecision, AltLocPolicy, AltLocPreference, AltLocResidue, AltLocSelection,
    AltLocSelectionReason, Assembly, AssemblyBond, AssemblyCut, AssemblyError, AssemblyGenerator,
    AssemblyOmissionReason, AssemblyOperation, AssemblyReport, AssemblyStereoLoss,
    AssemblyStructure, AtomProvenance, BondEvidence, BondProvenance, CellOmissionReason,
    ConflictingSelector, ConnectionResolutionReason, ConnectivityOptions, CoordinationProvenance,
    EntityKind, EntityMetadataPolicy, EntityMonomer, InstanceProvenance, ModelSelection,
    PolymerEntity, PolymerEntityInstance, PolymerType, SequencePosition, SequenceSource,
    SourceResidueId, SourceResiduePosition, StructureInterpretIssue, StructureInterpretOptions,
    StructureInterpretation, StructureInterpretationReport, StructureSetInterpretOptions,
    StructureSetInterpretation,
};
