//! Entity classification and source-provenance planning before row preparation.
use super::one_based_serial;
use super::{MmcifEntityClassifications, MmcifWriteError};
use crate::io::mmcif_interpret::{MmcifEntityKind, MmcifInterpretationReport};
use crate::structure::ModelView;
use crate::topology::{
    InstanceAtomId, MoleculeClass, MoleculeInstance, MoleculeInstanceId, ResidueClass, Topology,
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum EntityKind {
    Polymer,
    Branched,
    NonPolymer,
    Water,
}

impl EntityKind {
    pub(super) const fn as_mmcif(self) -> &'static str {
        match self {
            Self::Polymer => "polymer",
            Self::Branched => "branched",
            Self::NonPolymer => "non-polymer",
            Self::Water => "water",
        }
    }

    pub(super) const fn default_group_pdb(self) -> &'static str {
        match self {
            Self::Polymer => "ATOM",
            Self::Branched | Self::NonPolymer | Self::Water => "HETATM",
        }
    }

    fn as_source(self) -> MmcifEntityKind {
        match self {
            Self::Polymer => MmcifEntityKind::Polymer,
            Self::Branched => MmcifEntityKind::Branched,
            Self::NonPolymer => MmcifEntityKind::NonPolymer,
            Self::Water => MmcifEntityKind::Water,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct EntityRow {
    pub(super) id: String,
    pub(super) kind: EntityKind,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct AsymRow {
    pub(super) id: String,
    pub(super) entity_id: String,
}

#[derive(Debug, Clone)]
pub(super) struct AtomEntityAssignment {
    pub(super) entity_id: String,
    pub(super) asym_id: String,
    pub(super) kind: EntityKind,
}

#[derive(Debug, Clone)]
pub(super) struct EntityPlan {
    pub(super) entities: Vec<EntityRow>,
    pub(super) asyms: Vec<AsymRow>,
    pub(super) atoms: BTreeMap<InstanceAtomId, AtomEntityAssignment>,
}

pub(super) fn ensemble_entity_plan(
    model: ModelView<'_>,
    classifications: Option<&MmcifEntityClassifications>,
    report: Option<&MmcifInterpretationReport>,
) -> Result<EntityPlan, MmcifWriteError> {
    if let Some(report) = report {
        entity_plan_from_report(model, report)
    } else {
        let normalized = normalize_entity_classifications(
            model,
            classifications
                .into_iter()
                .flat_map(|set| set.iter().map(|(id, kind)| (id, vec![kind.clone()]))),
        )?;
        generic_entity_plan(model, &normalized)
    }
}

pub(super) fn entity_plan_from_report(
    model: ModelView<'_>,
    report: &MmcifInterpretationReport,
) -> Result<EntityPlan, MmcifWriteError> {
    if report
        .instances()
        .iter()
        .all(|instance| instance.atoms().is_empty())
    {
        let classifications = normalize_entity_classifications(
            model,
            report
                .instances()
                .iter()
                .map(|instance| (instance.molecule(), instance.entity_kinds().to_vec())),
        )?;
        generic_entity_plan(model, &classifications)
    } else {
        report_entity_plan(model, report)
    }
}

fn hierarchy_asym_ids(model: ModelView<'_>) -> Result<BTreeSet<String>, MmcifWriteError> {
    let mut reserved_asym_ids = BTreeSet::new();
    for (_, chain) in model.topology().hierarchy().chains() {
        if chain.label_id().is_empty() {
            return Err(MmcifWriteError::InvalidModel(
                "topology hierarchy chain label ID is empty".to_owned(),
            ));
        }
        if !reserved_asym_ids.insert(chain.label_id().to_owned()) {
            return Err(MmcifWriteError::DuplicateAsymId(
                chain.label_id().to_owned(),
            ));
        }
    }
    Ok(reserved_asym_ids)
}

pub(super) fn generic_entity_plan(
    model: ModelView<'_>,
    classifications: &BTreeMap<MoleculeInstanceId, EntityKind>,
) -> Result<EntityPlan, MmcifWriteError> {
    let hierarchy = model.topology().hierarchy();
    let mut reserved_asym_ids = hierarchy_asym_ids(model)?;
    let mut entities = Vec::new();
    let mut asyms = Vec::new();
    let mut atoms = BTreeMap::new();

    for (_, chain) in hierarchy.chains() {
        let sites = chain
            .residues()
            .iter()
            .flat_map(|residue| {
                hierarchy
                    .residue(*residue)
                    .expect("published hierarchy chain references a live residue")
                    .atom_sites()
                    .iter()
                    .copied()
            })
            .map(|site| {
                hierarchy
                    .atom_site(site)
                    .expect("published hierarchy residue references a live atom site")
            })
            .collect::<Vec<_>>();
        if sites.is_empty() {
            continue;
        }
        let source_kinds = sites
            .iter()
            .filter_map(|site| classifications.get(&site.atom().molecule()).copied())
            .map(EntityKind::as_source)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let kind = match source_kinds.as_slice() {
            [kind] => entity_kind_from_source(sites[0].atom().molecule(), kind)?,
            _ => {
                return Err(MmcifWriteError::ConflictingAsymEntityClassifications {
                    asym_id: chain.label_id().to_owned(),
                    classifications: source_kinds,
                });
            }
        };
        let entity_id = (entities.len() + 1).to_string();
        entities.push(EntityRow {
            id: entity_id.clone(),
            kind,
        });
        asyms.push(AsymRow {
            id: chain.label_id().to_owned(),
            entity_id: entity_id.clone(),
        });
        for site in sites {
            let atom = site.atom();
            if atoms
                .insert(
                    atom,
                    AtomEntityAssignment {
                        entity_id: entity_id.clone(),
                        asym_id: chain.label_id().to_owned(),
                        kind,
                    },
                )
                .is_some()
            {
                return Err(MmcifWriteError::DuplicateAtomSite(atom));
            }
        }
    }

    for (id, molecule) in model.topology().instances() {
        let definition = model
            .topology()
            .definition_for_instance(id)
            .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?;
        let qualified_atoms = definition
            .molecule()
            .atoms()
            .map(|(atom, _)| molecule.qualify_atom(atom))
            .collect::<Vec<_>>();
        let assigned = qualified_atoms
            .iter()
            .filter(|atom| atoms.contains_key(atom))
            .count();
        if assigned != 0 {
            if assigned != qualified_atoms.len() {
                let missing = qualified_atoms
                    .into_iter()
                    .find(|atom| !atoms.contains_key(atom))
                    .expect("partially assigned molecule has a missing atom site");
                return Err(MmcifWriteError::MissingAtomSite(missing));
            }
            continue;
        }
        let base = format!("M{}", one_based_serial(id.raw()));
        let mut candidate = base.clone();
        let mut suffix = 2usize;
        while reserved_asym_ids.contains(&candidate) {
            candidate = format!("{base}_{suffix}");
            suffix += 1;
        }
        reserved_asym_ids.insert(candidate.clone());
        let kind = entity_kind(molecule, classifications)?;
        let entity_id = (entities.len() + 1).to_string();
        entities.push(EntityRow {
            id: entity_id.clone(),
            kind,
        });
        asyms.push(AsymRow {
            id: candidate.clone(),
            entity_id: entity_id.clone(),
        });
        for atom in qualified_atoms {
            atoms.insert(
                atom,
                AtomEntityAssignment {
                    entity_id: entity_id.clone(),
                    asym_id: candidate.clone(),
                    kind,
                },
            );
        }
    }

    Ok(EntityPlan {
        entities,
        asyms,
        atoms,
    })
}

fn report_entity_plan(
    model: ModelView<'_>,
    report: &MmcifInterpretationReport,
) -> Result<EntityPlan, MmcifWriteError> {
    let _ = hierarchy_asym_ids(model)?;
    let mut seen_instances = BTreeSet::new();
    let mut provenance = BTreeMap::new();
    let mut reserved_entity_ids = BTreeSet::new();
    for instance in report.instances() {
        if model.topology().instance(instance.molecule()).is_err() {
            return Err(MmcifWriteError::UnknownClassifiedMolecule(
                instance.molecule(),
            ));
        }
        if !seen_instances.insert(instance.molecule()) {
            return Err(MmcifWriteError::DuplicateEntityClassification(
                instance.molecule(),
            ));
        }
        for atom in instance.atoms() {
            if atom.atom().molecule() != instance.molecule()
                || model.topology().atom(atom.atom()).is_err()
            {
                return Err(MmcifWriteError::UnknownAtomProvenance(atom.atom()));
            }
            if let Some(entity_id) = atom.entity_id() {
                reserved_entity_ids.insert(entity_id.to_owned());
            }
            if provenance.insert(atom.atom(), atom).is_some() {
                return Err(MmcifWriteError::DuplicateAtomProvenance(atom.atom()));
            }
        }
    }

    let hierarchy = model.topology().hierarchy();
    let mut entities = Vec::new();
    let mut entity_kinds = BTreeMap::<String, EntityKind>::new();
    let mut asyms = Vec::new();
    let mut atoms = BTreeMap::new();
    let mut generated_serial = 1usize;
    for (_, chain) in hierarchy.chains() {
        let sites = chain
            .residues()
            .iter()
            .flat_map(|residue| {
                hierarchy
                    .residue(*residue)
                    .expect("published hierarchy chain references a live residue")
                    .atom_sites()
                    .iter()
                    .copied()
            })
            .map(|site| {
                hierarchy
                    .atom_site(site)
                    .expect("published hierarchy residue references a live atom site")
            })
            .collect::<Vec<_>>();
        if sites.is_empty() {
            continue;
        }
        let source_atoms = sites
            .iter()
            .map(|site| {
                provenance
                    .get(&site.atom())
                    .copied()
                    .ok_or(MmcifWriteError::MissingAtomProvenance(site.atom()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        for atom in &source_atoms {
            if atom.asym_id() != chain.label_id() {
                return Err(MmcifWriteError::InconsistentAtomSite {
                    atom: atom.atom(),
                    field: "source asym_id",
                });
            }
        }
        let source_entity_ids = source_atoms
            .iter()
            .filter_map(|atom| atom.entity_id().map(str::to_owned))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        if source_entity_ids.len() > 1 {
            return Err(MmcifWriteError::ConflictingAsymEntityIds {
                asym_id: chain.label_id().to_owned(),
                entity_ids: source_entity_ids,
            });
        }
        let source_kinds = source_atoms
            .iter()
            .map(|atom| atom.entity_kind.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let kind = match source_kinds.as_slice() {
            [kind] => entity_kind_from_source(source_atoms[0].atom().molecule(), kind)?,
            _ => {
                return Err(MmcifWriteError::ConflictingAsymEntityClassifications {
                    asym_id: chain.label_id().to_owned(),
                    classifications: source_kinds,
                });
            }
        };
        let entity_id = if let Some(entity_id) = source_entity_ids.first() {
            entity_id.clone()
        } else {
            loop {
                let candidate = format!("K{generated_serial}");
                generated_serial += 1;
                if reserved_entity_ids.insert(candidate.clone()) {
                    break candidate;
                }
            }
        };
        if let Some(existing) = entity_kinds.get(&entity_id) {
            if *existing != kind {
                return Err(MmcifWriteError::ConflictingSourceEntityClassifications {
                    entity_id,
                    classifications: vec![existing.as_source(), kind.as_source()],
                });
            }
        } else {
            entity_kinds.insert(entity_id.clone(), kind);
            entities.push(EntityRow {
                id: entity_id.clone(),
                kind,
            });
        }
        asyms.push(AsymRow {
            id: chain.label_id().to_owned(),
            entity_id: entity_id.clone(),
        });
        for site in sites {
            let atom = site.atom();
            if atoms
                .insert(
                    atom,
                    AtomEntityAssignment {
                        entity_id: entity_id.clone(),
                        asym_id: chain.label_id().to_owned(),
                        kind,
                    },
                )
                .is_some()
            {
                return Err(MmcifWriteError::DuplicateAtomSite(atom));
            }
        }
    }

    for atom in model.topology().atom_ids() {
        if !atoms.contains_key(atom) {
            return Err(MmcifWriteError::MissingAtomProvenance(*atom));
        }
    }
    for atom in provenance.keys() {
        if !atoms.contains_key(atom) {
            return Err(MmcifWriteError::UnknownAtomProvenance(*atom));
        }
    }

    Ok(EntityPlan {
        entities,
        asyms,
        atoms,
    })
}

fn entity_kind(
    molecule: &MoleculeInstance,
    classifications: &BTreeMap<MoleculeInstanceId, EntityKind>,
) -> Result<EntityKind, MmcifWriteError> {
    classifications
        .get(&molecule.id())
        .copied()
        .ok_or(MmcifWriteError::MissingEntityClassification(molecule.id()))
}

fn entity_kind_from_source(
    molecule: MoleculeInstanceId,
    kind: &MmcifEntityKind,
) -> Result<EntityKind, MmcifWriteError> {
    match kind {
        MmcifEntityKind::Polymer => Ok(EntityKind::Polymer),
        MmcifEntityKind::Branched => Ok(EntityKind::Branched),
        MmcifEntityKind::NonPolymer => Ok(EntityKind::NonPolymer),
        MmcifEntityKind::Water => Ok(EntityKind::Water),
        MmcifEntityKind::Other(classification) => {
            Err(MmcifWriteError::UnsupportedEntityClassification {
                molecule,
                classification: classification.clone(),
            })
        }
    }
}

pub(super) fn normalize_entity_classifications(
    model: ModelView<'_>,
    entries: impl IntoIterator<Item = (MoleculeInstanceId, Vec<MmcifEntityKind>)>,
) -> Result<BTreeMap<MoleculeInstanceId, EntityKind>, MmcifWriteError> {
    let mut explicit = BTreeMap::new();
    for (molecule, source_kinds) in entries {
        if model.topology().instance(molecule).is_err() {
            return Err(MmcifWriteError::UnknownClassifiedMolecule(molecule));
        }
        if explicit.contains_key(&molecule) {
            return Err(MmcifWriteError::DuplicateEntityClassification(molecule));
        }
        let kinds = source_kinds
            .into_iter()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let kind = match kinds.as_slice() {
            [kind] => entity_kind_from_source(molecule, kind)?,
            [] => return Err(MmcifWriteError::MissingEntityClassification(molecule)),
            _ => {
                return Err(MmcifWriteError::ConflictingEntityClassifications {
                    molecule,
                    classifications: kinds,
                });
            }
        };
        explicit.insert(molecule, kind);
    }
    model
        .topology()
        .instances()
        .map(|(molecule, _)| {
            explicit
                .get(&molecule)
                .copied()
                .map(Ok)
                .unwrap_or_else(|| canonical_entity_kind(model.topology(), molecule))
                .map(|kind| (molecule, kind))
        })
        .collect()
}

fn canonical_entity_kind(
    topology: &Topology,
    molecule: MoleculeInstanceId,
) -> Result<EntityKind, MmcifWriteError> {
    let class = topology
        .molecule_class(molecule)
        .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?;
    match class {
        MoleculeClass::Protein | MoleculeClass::Dna | MoleculeClass::Rna => Ok(EntityKind::Polymer),
        MoleculeClass::Water => Ok(EntityKind::Water),
        MoleculeClass::Ion | MoleculeClass::SmallMolecule => Ok(EntityKind::NonPolymer),
        MoleculeClass::Carbohydrate => carbohydrate_entity_kind(topology, molecule),
        MoleculeClass::Other => Err(MmcifWriteError::UnresolvedCanonicalEntityClassification {
            molecule,
            classification: class,
        }),
    }
}

fn carbohydrate_entity_kind(
    topology: &Topology,
    molecule: MoleculeInstanceId,
) -> Result<EntityKind, MmcifWriteError> {
    let view = topology
        .molecule(molecule)
        .map_err(|error| MmcifWriteError::InvalidModel(error.to_string()))?;
    let mut residues = view.residues();
    if residues
        .next()
        .is_some_and(|residue| residue.class() == ResidueClass::Carbohydrate)
        && residues.next().is_none()
    {
        return Ok(EntityKind::NonPolymer);
    }
    Err(MmcifWriteError::UnresolvedCanonicalEntityClassification {
        molecule,
        classification: MoleculeClass::Carbohydrate,
    })
}
