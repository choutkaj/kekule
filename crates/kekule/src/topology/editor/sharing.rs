//! Definitions an appended molecule shares when it is identical to one the
//! draft already holds.
use super::*;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::hash::Hasher;

/// What a draft names one atom: its residue's component and its atom name,
/// empty without an atom site.
pub(super) type AtomLabel = (Box<str>, Box<str>);

/// The definitions a draft holds that an appended molecule may share, keyed
/// by the hash of their graph.
///
/// Two molecules are identical by the identity interpretation uses: equal
/// graphs in local atom order, equal properties, and equal per-atom labels.
/// Sharing also needs equal perception, class, and explicit class, so that a
/// shared definition describes every instance as its own would.
#[derive(Debug, Clone, Default)]
pub(super) struct SharedDefinitions {
    /// Whether the source snapshot's definitions are indexed.
    source_indexed: bool,
    by_graph: HashMap<u64, Vec<Shareable>>,
}

#[derive(Debug, Clone)]
struct Shareable {
    chemistry: Chemistry,
    /// Filled when a candidate first needs them.
    labels: Option<Arc<[AtomLabel]>>,
    class: MoleculeClass,
    class_explicit: bool,
}

#[derive(Debug, Clone)]
enum Chemistry {
    /// A source definition and its first instance, which names its labels.
    Source(MoleculeDefinitionId, MoleculeInstanceId),
    /// A definition appended earlier, and whether its perception is the
    /// default profile's.
    Added(Arc<Molecule>, bool),
}

/// The definitions an edit has left to be classified again: some draft
/// occurrence of each lost its class, so publication infers the class of the
/// whole definition anew. An appended molecule does not share them, so its
/// occurrences keep their own class.
#[derive(Default)]
pub(super) struct Reclassified {
    source: BTreeSet<MoleculeDefinitionId>,
    added: BTreeSet<*const Molecule>,
}

impl Reclassified {
    /// The definitions of the draft's `groups` that publication classifies
    /// again.
    pub(super) fn of(groups: &[GroupSlot]) -> Self {
        let mut reclassified = Self::default();
        for slot in groups {
            let GroupSlot::Live(group) = slot else {
                continue;
            };
            if group.class.is_some() {
                continue;
            }
            match &group.chemistry {
                GroupChemistry::Source(id) => {
                    reclassified.source.insert(*id);
                }
                GroupChemistry::Added(molecule, _) => {
                    reclassified.added.insert(Arc::as_ptr(molecule));
                }
                GroupChemistry::Draft(_) => {}
            }
        }
        reclassified
    }

    fn contains(&self, chemistry: &Chemistry) -> bool {
        match chemistry {
            Chemistry::Source(id, _) => self.source.contains(id),
            Chemistry::Added(molecule, _) => self.added.contains(&Arc::as_ptr(molecule)),
        }
    }
}

/// One definition of an appended topology and an instance of it.
pub(super) struct Appended<'a> {
    pub(super) molecule: &'a Molecule,
    pub(super) default_perception: bool,
    pub(super) labels: &'a Arc<[AtomLabel]>,
    pub(super) class: MoleculeClass,
    pub(super) class_explicit: bool,
}

fn graph_hash(molecule: &Molecule) -> u64 {
    let mut hasher = DefaultHasher::new();
    molecule.hash_graph(&mut hasher);
    hasher.finish()
}

/// The labels of one instance's atoms in local atom order.
pub(super) fn instance_labels(
    topology: &Topology,
    instance: MoleculeInstanceId,
) -> Arc<[AtomLabel]> {
    topology
        .instance(instance)
        .expect("labelled instances exist")
        .atoms()
        .map(|atom| match atom.atom_site() {
            Some(site) => {
                let residue = site.residue();
                let component = residue
                    .label_comp_id()
                    .or(residue.author_comp_id())
                    .unwrap_or_else(|| residue.name());
                let metadata = site.metadata();
                let name = metadata
                    .label_atom_id()
                    .or(metadata.auth_atom_id())
                    .unwrap_or_default();
                (component.into(), name.into())
            }
            None => (Box::default(), Box::default()),
        })
        .collect()
}

impl SharedDefinitions {
    /// The chemistry of an identical definition the draft holds, from its
    /// `source` snapshot or an earlier append, if any, unless it is
    /// `reclassified`.
    pub(super) fn find(
        &mut self,
        source: Option<&Topology>,
        reclassified: &Reclassified,
        appended: &Appended<'_>,
    ) -> Option<GroupChemistry> {
        if !self.source_indexed {
            self.source_indexed = true;
            if let Some(source) = source {
                self.index_source(source);
            }
        }
        let candidates = self.by_graph.get_mut(&graph_hash(appended.molecule))?;
        for candidate in candidates {
            if candidate.class != appended.class
                || candidate.class_explicit != appended.class_explicit
                || reclassified.contains(&candidate.chemistry)
            {
                continue;
            }
            let (molecule, default_perception) = match &candidate.chemistry {
                Chemistry::Source(id, _) => {
                    let definition = source
                        .and_then(|source| source.definition(*id))
                        .expect("indexed source definitions are the source's");
                    (definition.molecule(), definition.default_perception)
                }
                Chemistry::Added(molecule, default_perception) => {
                    (&**molecule, *default_perception)
                }
            };
            let identical = default_perception == appended.default_perception
                && *molecule == *appended.molecule
                && molecule.properties == appended.molecule.properties
                && molecule.perception == appended.molecule.perception;
            if !identical {
                continue;
            }
            let labels = candidate
                .labels
                .get_or_insert_with(|| match &candidate.chemistry {
                    Chemistry::Source(_, instance) => {
                        instance_labels(source.expect("source candidates have a source"), *instance)
                    }
                    Chemistry::Added(..) => unreachable!("appended definitions are labelled"),
                });
            if labels != appended.labels {
                continue;
            }
            return Some(match &candidate.chemistry {
                Chemistry::Source(id, _) => GroupChemistry::Source(*id),
                Chemistry::Added(molecule, default_perception) => {
                    GroupChemistry::Added(Arc::clone(molecule), *default_perception)
                }
            });
        }
        None
    }

    /// Offers an appended definition to later appends.
    pub(super) fn insert(&mut self, molecule: Arc<Molecule>, appended: &Appended<'_>) {
        self.by_graph
            .entry(graph_hash(&molecule))
            .or_default()
            .push(Shareable {
                chemistry: Chemistry::Added(molecule, appended.default_perception),
                labels: Some(Arc::clone(appended.labels)),
                class: appended.class,
                class_explicit: appended.class_explicit,
            });
    }

    /// Forgets every definition, as clearing the draft does.
    pub(super) fn clear(&mut self) {
        self.by_graph.clear();
        self.source_indexed = true;
    }

    fn index_source(&mut self, source: &Topology) {
        let mut first = vec![None; source.definition_count()];
        for (instance, record) in source.instance_records() {
            first[record.definition().index()].get_or_insert(instance);
        }
        for definition in source.definitions() {
            let Some(instance) = first[definition.id().index()] else {
                continue;
            };
            self.by_graph
                .entry(graph_hash(definition.molecule()))
                .or_default()
                .push(Shareable {
                    chemistry: Chemistry::Source(definition.id(), instance),
                    labels: None,
                    class: definition.class(),
                    class_explicit: source
                        .molecule_class_overrides
                        .contains_key(&definition.id()),
                });
        }
    }
}
