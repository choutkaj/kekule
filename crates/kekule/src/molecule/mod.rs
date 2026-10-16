//! Connected molecular graphs and their represented and perceived chemistry.
//!
//! [`Molecule`] is the foundational chemical object. Its [`Graph`] stores
//! authoritative atoms, bonds, connectivity, and represented stereochemistry
//! (whose types live in [`crate::stereo`]); [`Perception`] stores
//! reconstructible derived chemistry. Use
//! [`MoleculeEditor`] for construction and structural mutation so the
//! non-empty, connected publication invariant is checked transactionally.
//!
//! Coordinates and residue/chain hierarchy deliberately do not live here. See
//! [`crate::structure`] and [`crate::topology`] for those layers.

mod atom_bond;
pub(crate) mod delocalization;
mod element;
mod element_reference;
mod graph;
pub(crate) mod ids;
mod molecule_edit;
mod molecule_editor_ops;
mod molecule_workflows;
mod owner;
pub(crate) mod perception;
pub(crate) mod stereo;

pub use atom_bond::{Atom, AtomRadical, Bond, BondOrder, ImplicitHydrogens};
pub(crate) use delocalization::{
    ConjugationModel, ConjugationPerception, ResonanceGroup, ResonancePerception,
};
pub use element::Element;
pub use element_reference::{ElementReferenceData, VanDerWaalsRadiusSource};
pub use graph::Graph;
pub(crate) use ids::{checked_fixed_id_collection_len, checked_raw_id};
pub use ids::{AtomId, BondId};
pub(crate) use ids::{StereoElementId, StereoGroupId};
pub(crate) use molecule_edit::canonicalize_represented_chemistry;
pub use molecule_edit::{
    GraphValidationError, MoleculeConnectivityError, MoleculeEditor, MoleculeFinishError,
    MoleculePublicationError, StereoPublicationError,
};
pub use molecule_editor_ops::MoleculeCorrespondence;
pub use owner::{AtomMut, BondMut, Molecule, MoleculeError, MoleculeIdKind};
pub(crate) use owner::{Result, SlotCompaction};
use perception::validate_perception;
pub(crate) use perception::{
    AromaticityModel, AromaticityPerception, Ring, RingBasisModel, RingBasisState, RingMembership,
    RingPerception, RingSet, ValenceModel, ValencePerception,
};
pub(crate) use perception::{CipPerception, Perception, PerceptionInstallError};
pub(crate) use stereo::{has_axis_stereo, has_double_bond_stereo, has_tetrahedral_stereo};
pub(crate) use stereo::{
    AxisOrientation, AxisStereo, DoubleBondOrientation, DoubleBondStereo, StereoCarrier,
    StereoDescriptor, StereoElement, StereoElementKind, StereoGroup, StereoGroupKind,
    TetrahedralOrientation, TetrahedralStereo,
};
