//! Connected molecular graphs and their represented and perceived chemistry.
//!
//! [`Molecule`] is the foundational chemical object: one published, non-empty,
//! connected, geometry-independent molecular entity. Its [`Graph`] stores
//! authoritative represented chemistry: [`Atom`]s (element, isotope, formal
//! charge, [`AtomRadical`], [`ImplicitHydrogens`], atom map) and [`Bond`]s
//! with a localized [`BondOrder`], addressed by molecule-local [`AtomId`] and
//! [`BondId`], plus represented stereochemistry (whose types live in
//! [`crate::stereo`]). Derived chemistry, such as valence, rings, and
//! aromaticity, is reconstructible [`Perception`] installed by
//! [`Molecule::perceive`].
//!
//! Use [`MoleculeEditor`] for construction and structural mutation so the
//! non-empty, connected publication invariant is checked transactionally;
//! [`AtomMut`] and [`BondMut`] are its checked mutable views, and
//! [`MoleculeCorrespondence`] maps IDs across appends and publication.
//! Operations fail with [`MoleculeError`] (naming a [`MoleculeIdKind`] when an
//! ID space is exhausted) or, at publication, [`MoleculePublicationError`].
//! [`Element`] provides element identity and reference data, including
//! van der Waals radii selected by [`VanDerWaalsRadiusSource`]; unknown
//! symbols and atomic numbers are reported as [`ElementError`].
//!
//! Coordinates and residue/chain hierarchy deliberately do not live here. See
//! [`crate::structure`] and [`crate::topology`] for those layers.
//!
//! ```
//! use kekule::molecule::{AtomId, BondId, BondOrder, Element};
//!
//! let mut phenol = kekule::smiles::to_molecules("Oc1ccccc1")?.pop().unwrap();
//! assert_eq!(phenol.atom_count(), 7);
//! assert_eq!(phenol.atom(AtomId::new(0)).unwrap().element(), Element::O);
//! assert_eq!(phenol.bond(BondId::new(0)).unwrap().order, BondOrder::Single);
//! // Aromaticity is perceived state, not a represented bond order.
//! assert!(phenol.perception().aromaticity().is_none());
//! phenol.perceive()?;
//! let aromaticity = phenol.perception().aromaticity().expect("perceived");
//! assert!(aromaticity.atom_is_aromatic(AtomId::new(1)));
//! assert_eq!(phenol.implicit_hydrogens(AtomId::new(0)), Some(1));
//! let rings = phenol.perception().rings().and_then(|rings| rings.ring_set());
//! assert_eq!(rings.map(|rings| rings.len()), Some(1));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

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
pub use element::{Element, ElementError};
pub use element_reference::VanDerWaalsRadiusSource;
pub use graph::Graph;
pub(crate) use ids::{checked_fixed_id_collection_len, checked_raw_id};
pub use ids::{AtomId, BondId};
pub(crate) use ids::{StereoElementId, StereoGroupId};
pub(crate) use molecule_edit::{canonicalize_represented_chemistry, validate_represented_stereo};
pub use molecule_edit::{MoleculeEditor, MoleculePublicationError};
pub use molecule_editor_ops::MoleculeCorrespondence;
pub use owner::{AtomMut, BondMut, Molecule, MoleculeError, MoleculeIdKind};
pub(crate) use owner::{Result, SlotCompaction};
use perception::validate_perception;
pub(crate) use perception::{
    AromaticityModel, AromaticityPerception, Ring, RingBasisModel, RingBasisState, RingMembership,
    RingPerception, RingSet, ValenceModel, ValencePerception, ValenceValidation,
};
pub(crate) use perception::{CipPerception, Perception, PerceptionInstallError};
use std::num::{NonZeroU16, NonZeroU32};
pub(crate) use stereo::{has_axis_stereo, has_double_bond_stereo, has_tetrahedral_stereo};
pub(crate) use stereo::{
    AxisOrientation, AxisStereo, DoubleBondOrientation, DoubleBondStereo, StereoCarrier,
    StereoDescriptor, StereoElement, StereoElementKind, StereoGroup, StereoGroupKind,
    TetrahedralOrientation, TetrahedralStereo,
};
