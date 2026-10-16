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

pub(crate) use crate::core::RingMembership;
pub(crate) use aromaticity::perceive_aromaticity_in_place;
pub use aromaticity::{
    perceive_aromaticity, perceive_aromaticity_with_options,
    perceive_aromaticity_with_ring_options, AromaticityError, AromaticityOptions,
};
pub(crate) use canonical::canonical_atom_ranking;
pub use canonical::CanonicalAtomRanking;
pub(crate) use cip::bond_reference_atoms;
pub use cip::{
    assign_cip_descriptors, assign_cip_descriptors_with_options, CipAssignment, CipAssignmentError,
    CipAssignmentIssue, CipAssignmentOptions, CipAssignmentReport, CipRankingError, CipSkipped,
    CipSkippedReason,
};
pub use conjugation::{perceive_conjugation, ConjugationError};
pub(crate) use hydrogens::{
    add_hydrogens_to_molecule, remove_hydrogens_from_molecule, remove_hydrogens_keeping,
};
pub use hydrogens::{
    AddHydrogensOptions, AddHydrogensReport, AddedHydrogen, AddedHydrogenOrigin,
    HydrogenCountAdjustment, HydrogenTransformError, RemoveHydrogensReport, RemovedHydrogen,
    RetainedHydrogen, RetainedHydrogenReason,
};
pub use resonance::{
    enumerate_resonance, perceive_resonance, ResonanceContributor, ResonanceError, ResonanceFlags,
    ResonanceOptions, ResonanceStructures,
};
pub(crate) use rings::{compute_ring_membership, ordered_atom_pair};
pub use rings::{
    perceive_ring_membership, perceive_ring_set, perceive_ring_set_with_options,
    RingPerceptionError, RingPerceptionOptions,
};
pub(crate) use rotatable_bonds::detect_rotatable_bonds;
pub use rotatable_bonds::{RotatableBondOptions, RotatableBondSet};
pub(crate) use stereo::{
    atom_axis_carriers, atom_hydrogen_count, atom_is_atropisomeric_sp2_endpoint,
    coordinates_are_planar, coordinates_contradict, double_bond_endpoint_carriers,
    double_bond_geometry_is_supported, double_bond_orientation_from_points,
    locally_eligible_double_bonds, nitrogen_lone_pair_is_stereogenic, tetrahedral_orientation_at,
    tetrahedral_orientation_from_points,
};
pub use stereo::{
    cleanup_stereo, detect_stereo_candidates, detect_stereo_candidates_with_options,
    infer_coordinate_stereo, infer_coordinate_stereo_with_options, materialize_coordinate_stereo,
    materialize_coordinate_stereo_with_options, validate_stereo, CoordinateStereoError,
    CoordinateStereoMaterializationReport, CoordinateStereoOptions, CoordinateStereoResult,
    StereoCandidate, StereoCleanupReport, StereoPerceptionError, StereoPerceptionOptions,
    StereoValidationError, StereoValidationIssue,
};
pub use substructure::{
    find_match, find_matches, find_matches_with_options, find_topology_matches,
    find_topology_matches_with_options, visit_matches, visit_matches_with_options, MatchCompletion,
    PreparedTarget, PreparedTopologyTarget, QueryMatch, QueryPerception, SubstructureMatchError,
    SubstructureMatchOptions, SubstructureMatchWork, TaggedMatchError, TaggedQuery,
    TopologyQueryMatch, MAX_SUBSTRUCTURE_QUERY_ATOMS,
};
pub(crate) use valence::{
    allowed_valences, explicit_valence, rdkit_bracket_radical_electrons,
    rdkit_charge_adjusted_default_valence, rdkit_default_valence, rdkit_inferred_hydrogen_count,
};
pub use valence::{
    perceive_valence, perceive_valence_with_options, represented_valence, ValenceError,
    ValenceIssue, ValenceOptions,
};
