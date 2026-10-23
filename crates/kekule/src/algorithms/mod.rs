mod aromaticity;
mod canonical;
mod cip;
mod conjugation;
mod hydrogens;
mod resonance;
mod rings;
mod rotatable_bonds;
mod stereo;
mod substructure;
mod valence;

pub(crate) use crate::molecule::RingMembership;
pub(crate) use aromaticity::perceive_aromaticity_in_place;
pub use aromaticity::{perceive_aromaticity, AromaticityError, AromaticityOptions};
pub use canonical::{canonical_atom_ranking, CanonicalAtomRanking};
pub(crate) use cip::bond_reference_atoms;
pub use cip::{
    assign_cip_descriptors, CipAssignment, CipAssignmentError, CipAssignmentIssue,
    CipAssignmentOptions, CipAssignmentReport, CipRankingError, CipSkipped, CipSkippedReason,
};
pub use conjugation::{perceive_conjugation, ConjugationError};
pub(crate) use hydrogens::{
    add_hydrogens_to_molecule, remove_hydrogens_from_molecule, remove_hydrogens_keeping,
};
pub use hydrogens::{
    AddHydrogensOptions, AddHydrogensReport, AddHydrogensScope, AddedHydrogen, AddedHydrogenOrigin,
    HydrogenCountAdjustment, HydrogenTransformError, RemoveHydrogensReport, RemovedHydrogen,
    RetainedHydrogen, RetainedHydrogenReason,
};
pub use resonance::{
    enumerate_resonance, perceive_resonance, ResonanceContributor, ResonanceError, ResonanceFlags,
    ResonanceOptions, ResonanceStructures,
};
pub(crate) use rings::{compute_ring_membership, ordered_atom_pair, perceive_rings_cancellable};
pub use rings::{perceive_rings, RingPerceptionError, RingPerceptionOptions, RingPerceptionScope};
pub use rotatable_bonds::{rotatable_bonds, RotatableBondOptions, RotatableBondSet};
pub(crate) use stereo::{
    atom_axis_carriers, atom_hydrogen_count, atom_is_atropisomeric_sp2_endpoint,
    coordinates_are_planar, coordinates_contradict, double_bond_endpoint_carriers,
    double_bond_geometry_is_supported, double_bond_orientation_from_points,
    locally_eligible_double_bonds, nitrogen_lone_pair_is_stereogenic, tetrahedral_orientation_at,
    tetrahedral_orientation_from_points,
};
pub use stereo::{
    cleanup_stereo, detect_stereo_candidates, infer_coordinate_stereo,
    materialize_coordinate_stereo, CoordinateStereoError, CoordinateStereoMaterializationReport,
    CoordinateStereoOptions, CoordinateStereoResult, StereoCandidate, StereoCleanupReport,
    StereoSearchError, StereoSearchOptions, StereoValidationError, StereoValidationIssue,
};
pub use substructure::{
    MatchCompletion, PreparedTarget, PreparedTopologyTarget, QueryMatch, QueryPerception,
    SubstructureMatchError, SubstructureMatchOptions, SubstructureMatchWork, TaggedMatchError,
    TaggedQuery, TopologyQueryMatch, TupleMode, MAX_SUBSTRUCTURE_QUERY_ATOMS,
};
pub(crate) use valence::{
    allowed_valences, explicit_valence, rdkit_bracket_radical_electrons,
    rdkit_charge_adjusted_default_valence, rdkit_default_valence, rdkit_inferred_hydrogen_count,
};
pub use valence::{
    perceive_valence, represented_valence, ValenceError, ValenceIssue, ValenceOptions,
};
