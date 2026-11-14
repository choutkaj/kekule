use std::fmt;

use crate::algorithms::{
    perceive_aromaticity_in_place, perceive_rings_cancellable, perceive_valence, AromaticityError,
    RingPerceptionError, RingPerceptionOptions, ValenceError, ValenceOptions,
};
use crate::algorithms::{
    perceive_conjugation, perceive_resonance_cancellable, ConjugationError,
    ResonancePerceptionError,
};
use crate::diagnostics::Cancelled;
use crate::molecule::{AromaticityModel, Molecule, ValenceModel};
use crate::molecule::{ConjugationModel, ResonanceModel};
use crate::parallel::{Progress, Run};

/// Failure from the default discrete chemical perception profile.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PerceptionError {
    /// Conjugation perception failed; see the source.
    Conjugation(ConjugationError),
    /// Resonance perception failed; see the source.
    Resonance(ResonancePerceptionError),
    /// Valence perception failed; see the source.
    Valence(ValenceError),
    /// Ring perception failed; see the source.
    Rings(RingPerceptionError),
    /// Aromaticity perception failed; see the source.
    Aromaticity(AromaticityError),
    /// The run's [`Cancellation`](crate::diagnostics::Cancellation) stopped
    /// perception, as a `perceive_with` form such as
    /// [`Molecule::perceive_with`] reports.
    Cancelled(Cancelled),
}

impl fmt::Display for PerceptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Conjugation(_) => f.write_str("Cannot perceive conjugation."),
            Self::Resonance(_) => f.write_str("Cannot perceive resonance."),
            Self::Valence(_) => f.write_str("Cannot perceive valence."),
            Self::Rings(_) => f.write_str("Cannot perceive rings."),
            Self::Aromaticity(_) => f.write_str("Cannot perceive aromaticity."),
            Self::Cancelled(_) => f.write_str("Perception stopped."),
        }
    }
}

impl std::error::Error for PerceptionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Conjugation(error) => Some(error),
            Self::Resonance(error) => Some(error),
            Self::Valence(error) => Some(error),
            Self::Rings(error) => Some(error),
            Self::Aromaticity(error) => Some(error),
            Self::Cancelled(cancelled) => Some(cancelled),
        }
    }
}

/// Install the default discrete perception state on canonical chemistry.
///
/// The fixed profile is RDKit-like valence and implicit hydrogens, followed by
/// the default deterministic ring set, RDKit-like aromaticity and conjugation,
/// and charge-move resonance delocalization. The complete
/// operation is transactional and never rewrites represented chemistry.
pub(crate) fn perceive_molecule(molecule: &mut Molecule) -> Result<(), PerceptionError> {
    perceive_molecule_with(molecule, Run::default())
}

/// The stages of [`perceive_molecule_with`]: valence, rings, aromaticity,
/// conjugation, and resonance.
const STAGES: usize = 5;

/// [`perceive_molecule`] under `run`, polling its cancellation between
/// stages and during ring and resonance perception, which dominate on large
/// molecules, and reporting the stages done.
pub(crate) fn perceive_molecule_with(
    molecule: &mut Molecule,
    run: Run<'_>,
) -> Result<(), PerceptionError> {
    let previous = molecule.perception().clone();
    if let Err(error) = perceive_molecule_in_place(molecule, run) {
        molecule
            .install_perception(previous)
            .expect("previous perception state must remain valid");
        return Err(error);
    }
    run.report(Progress::of(STAGES, STAGES));
    Ok(())
}

fn perceive_molecule_in_place(
    molecule: &mut Molecule,
    run: Run<'_>,
) -> Result<(), PerceptionError> {
    let cancellation = run.cancellation();
    let stage = |done: usize| {
        run.report(Progress::of(done, STAGES));
        run.check().map_err(PerceptionError::Cancelled)
    };
    stage(0)?;
    perceive_valence(molecule, ValenceModel::RdkitLike, ValenceOptions::default())
        .map_err(PerceptionError::Valence)?;
    stage(1)?;
    let ring_options = RingPerceptionOptions::for_molecule(molecule);
    perceive_rings_cancellable(molecule, ring_options, cancellation).map_err(
        |error| match error {
            RingPerceptionError::Cancelled(cancelled) => PerceptionError::Cancelled(cancelled),
            error => PerceptionError::Rings(error),
        },
    )?;
    stage(2)?;
    perceive_aromaticity_in_place(molecule, AromaticityModel::RdkitLike)
        .map_err(PerceptionError::Aromaticity)?;
    stage(3)?;
    perceive_conjugation(molecule, ConjugationModel::RdkitLike)
        .map_err(PerceptionError::Conjugation)?;
    stage(4)?;
    perceive_resonance_cancellable(molecule, ResonanceModel::ChargeMoves, cancellation).map_err(
        |error| match error {
            ResonancePerceptionError::Cancelled(cancelled) => PerceptionError::Cancelled(cancelled),
            error => PerceptionError::Resonance(error),
        },
    )
}
