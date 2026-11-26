use std::fmt;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, PoisonError};

use crate::algorithms::{perceive_molecule_with, PerceptionError};
use crate::parallel::{map_weighted_ranges, Progress, Run};

use crate::molecule::perception::PerceptionWarning;
use crate::molecule::Molecule;

use super::{MoleculeDefinition, MoleculeDefinitionId, Topology};

/// Default perception failed for one reusable molecule definition.
///
/// The definition ID belongs to the source topology. Every source owner and
/// its installed perception remain unchanged on failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologyPerceptionError {
    /// The failed definition, in the source topology.
    pub definition: MoleculeDefinitionId,
    /// The perception failure.
    pub source: PerceptionError,
}

impl fmt::Display for TopologyPerceptionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "Cannot perceive molecule {}.", self.definition)
    }
}

impl std::error::Error for TopologyPerceptionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.source)
    }
}

impl Topology {
    /// Returns a new snapshot with the default molecular perception profile.
    ///
    /// Each reusable definition is perceived once, in definition order, using
    /// [`crate::molecule::Molecule::perceive`]: RDKit-like valence and implicit
    /// hydrogens, the default deterministic ring set, and RDKit-like aromaticity.
    /// Existing perception is recomputed and dependent CIP state is cleared,
    /// except on definitions whose perception this profile installed and
    /// nothing changed since, such as the definitions a structural edit left
    /// untouched: perception is deterministic, so they are shared with the
    /// result as they are. Perceiving after an edit costs the definitions the
    /// edit changed, whatever the size of the system.
    /// Failure identifies the first failing definition and leaves `self` intact.
    ///
    /// Represented graphs, definition reuse, all semantic IDs, dense ordering,
    /// hierarchy, classifications, and properties are preserved exactly, so the
    /// result shares this topology's layout and its static storage without
    /// copying it ([`Self::shares_layout`]). Selections, structures, buffers,
    /// and prepared potentials bound to this topology remain usable with the
    /// result. Values that read perception, such as prepared substructure
    /// targets, stay bound to the snapshot they were prepared from.
    ///
    /// ```
    /// use kekule::{smiles, topology::Topology};
    ///
    /// let molecules = smiles::to_molecules("c1ccccc1.[Na+]")?;
    /// let source = Topology::from_molecules(molecules)?;
    /// let perceived = source.perceived()?;
    /// assert!(perceived.shares_layout(&source));
    /// assert!(perceived.instances().all(|m| m.molecule().perception().aromaticity().is_some()));
    /// assert!(source.instances().all(|m| m.molecule().perception().aromaticity().is_none()));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn perceived(&self) -> Result<Self, TopologyPerceptionError> {
        self.perceived_with(Run::default())
    }

    /// [`Self::perceived`] under a [`Run`]. Definitions are independent, so
    /// the run's executor perceives them in parallel, in work items of about
    /// equal atom counts; the result, and the error, are those of the calling
    /// thread alone. Its cancellation is polled between definitions and
    /// inside each, as [`Molecule::perceive_with`] polls it. Progress counts
    /// the definitions done, of the definition count; the shared ones count
    /// first.
    ///
    /// # Errors
    ///
    /// The errors of [`Self::perceived`], and, once `run`'s cancellation is
    /// requested, an error whose source is [`PerceptionError::Cancelled`],
    /// naming a definition being perceived.
    pub fn perceived_with(&self, run: Run<'_>) -> Result<Self, TopologyPerceptionError> {
        let total = self.definitions.len();
        run.report(Progress::of(0, total));
        // Perception is deterministic, so a definition already holding
        // exactly this profile's result is shared, not perceived again.
        let pending = self
            .definitions
            .iter()
            .filter(|definition| !definition.default_perception)
            .collect::<Vec<_>>();
        let shared = total - pending.len();
        if shared > 0 && !pending.is_empty() {
            run.report(Progress::of(shared, total));
        }
        let mut perceived = perceive_definitions(&pending, run, shared)?.into_iter();
        let definitions = self
            .definitions
            .iter()
            .map(|definition| {
                if definition.default_perception {
                    return definition.clone();
                }
                MoleculeDefinition {
                    molecule: perceived
                        .next()
                        .expect("one molecule per pending definition"),
                    default_perception: true,
                    ..definition.clone()
                }
            })
            .collect();
        run.report(Progress::of(total, total));

        // Perception never requires reclassification, reindexing, or rebuilding
        // hierarchy correspondence, so the snapshot keeps the published layout.
        Ok(Self {
            definitions,
            instances: self.instances.clone(),
            atom_order: Arc::clone(&self.atom_order),
            bond_order: Arc::clone(&self.bond_order),
            hierarchy: Arc::clone(&self.hierarchy),
            hierarchy_index: Arc::clone(&self.hierarchy_index),
            molecule_class_overrides: self.molecule_class_overrides.clone(),
            properties: Arc::clone(&self.properties),
            layout: self.layout,
        })
    }

    /// The warnings of each definition's installed perception
    /// ([`crate::perception::Perception::warnings`]), in definition order.
    /// Empty when perception is complete or not installed.
    pub fn perception_warnings(
        &self,
    ) -> impl Iterator<Item = (MoleculeDefinitionId, PerceptionWarning)> + '_ {
        self.definitions.iter().flat_map(|definition| {
            definition
                .molecule
                .perception()
                .warnings()
                .into_iter()
                .map(move |warning| (definition.id, warning))
        })
    }
}

/// Perceives `pending` on `run`'s executor and returns their molecules in
/// order. Each work item, a contiguous range of about equal atom count,
/// stops at its first failure, or once an earlier definition has failed, so
/// the first failure in range order is the first in definition order, as on
/// one thread. Progress counts on from `done` definitions.
fn perceive_definitions(
    pending: &[&MoleculeDefinition],
    run: Run<'_>,
    done: usize,
) -> Result<Vec<Arc<Molecule>>, TopologyPerceptionError> {
    // Each definition polls the cancellation; only the topology reports
    // progress.
    let inner = match run.cancellation() {
        Some(cancellation) => Run::new().with_cancellation(cancellation),
        None => Run::new(),
    };
    let total = done + pending.len();
    // Counting and reporting under one lock keeps calls one at a time and
    // their counts nondecreasing.
    let completed = Mutex::new(done);
    let first_failure = AtomicUsize::new(usize::MAX);
    let weights = pending
        .iter()
        .map(|definition| definition.molecule.atom_count())
        .collect::<Vec<_>>();
    let ranges = map_weighted_ranges(run.executor(), &weights, |range| {
        let mut molecules = Vec::with_capacity(range.len());
        for position in range {
            if first_failure.load(Ordering::Relaxed) < position {
                // An earlier failure is the error, whatever this range holds.
                return Err(None);
            }
            let definition = pending[position];
            let mut molecule = Molecule::clone(&definition.molecule);
            if let Err(source) = perceive_molecule_with(&mut molecule, inner) {
                first_failure.fetch_min(position, Ordering::Relaxed);
                return Err(Some(TopologyPerceptionError {
                    definition: definition.id,
                    source,
                }));
            }
            molecules.push(Arc::new(molecule));
            if run.reports_progress() {
                let mut completed = completed.lock().unwrap_or_else(PoisonError::into_inner);
                *completed += 1;
                // The topology reports the last count once all are done.
                if *completed < total {
                    run.report(Progress::of(*completed, total));
                }
            }
        }
        Ok(molecules)
    });
    let mut molecules = Vec::with_capacity(pending.len());
    for range in ranges {
        match range {
            Ok(range) => molecules.extend(range),
            Err(error) => return Err(error.expect("a stopped range follows a failed one")),
        }
    }
    Ok(molecules)
}
