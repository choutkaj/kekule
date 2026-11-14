//! Derived valence, rings, aromaticity, conjugation, resonance and explicit
//! resonance enumeration.
//!
//! [`crate::molecule::Molecule::perceive`] installs Kekule's default transactional perception
//! profile. The nested modules expose the individual expert algorithms.
//! Perception does not alter represented graph chemistry and is invalidated by
//! relevant graph edits.
//!
//! A [`Perception`] holds one optional section per
//! step (valence, rings, aromaticity, conjugation, resonance and CIP). Each
//! section records the model or policy that produced it, and an absent
//! section has not been perceived.
//!
//! ```
//! use kekule::smiles;
//!
//! let mut phenol = smiles::to_molecules("Oc1ccccc1")?.remove(0);
//! assert!(phenol.perception().aromaticity().is_none());
//! phenol.perceive()?;
//!
//! let perception = phenol.perception();
//! let rings = perception.rings().and_then(|rings| rings.ring_set());
//! assert_eq!(rings.map(|rings| rings.len()), Some(1));
//! let aromaticity = perception.aromaticity().expect("perceived");
//! assert_eq!(aromaticity.atoms().len(), 6);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

/// Resonance delocalization perception, installed by default perception,
/// and explicit RDKit-like contributor enumeration.
pub mod resonance {
    pub use crate::algorithms::{
        enumerate_resonance, perceive_resonance, ResonanceContributor, ResonanceError,
        ResonanceFlags, ResonanceOptions, ResonancePerceptionError, ResonanceStructures,
    };
    pub use crate::molecule::delocalization::{
        ResonanceGroup, ResonanceModel, ResonancePerception,
    };
}
/// Bond conjugation, including aromatic bonds, without contributor enumeration.
pub mod conjugation {
    pub use crate::algorithms::{perceive_conjugation, ConjugationError};
    pub use crate::molecule::delocalization::{ConjugationModel, ConjugationPerception};
}
pub use crate::chemistry::PerceptionError;
pub use crate::molecule::perception::{
    DelocalizationIssue, Perception, PerceptionBuildError, PerceptionBuilder, PerceptionComponent,
    PerceptionInstallError,
};

/// Expert valence perception for canonical represented chemistry.
///
/// The RDKit-like model derives the inferred contribution to implicit H
/// from ordinary localized bond orders and represented atom state. It does
/// not require installed ring or aromaticity perception. Source aromatic
/// bonds are localized during format interpretation before this layer is
/// reached.
pub mod valence {
    pub use crate::algorithms::{
        perceive_valence, represented_valence, ValenceError, ValenceIssue, ValenceOptions,
    };
    pub use crate::molecule::perception::{ValenceModel, ValencePerception, ValenceValidation};
}

/// Ring membership and deterministic ring-basis perception.
///
/// Ring perception installs derived state and never changes represented
/// bond orders.
pub mod rings {
    pub use crate::algorithms::{
        perceive_rings, RingPerceptionError, RingPerceptionOptions, RingPerceptionScope,
    };
    pub use crate::molecule::perception::{
        MalformedRingReason, Ring, RingBasisModel, RingBasisState, RingMembership, RingPerception,
        RingSet,
    };
}

/// Aromaticity perception over canonical localized bond orders.
///
/// Aromatic membership is derived state; Kekule does not store an aromatic
/// bond order in represented graph chemistry.
pub mod aromaticity {
    pub use crate::algorithms::{perceive_aromaticity, AromaticityError, AromaticityOptions};
    pub use crate::molecule::perception::{AromaticityModel, AromaticityPerception};
}
