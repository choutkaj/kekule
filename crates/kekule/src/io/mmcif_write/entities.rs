//! Entity classification and source-provenance planning before row preparation.
use super::one_based_serial;
use super::{fail, invalid_model, MmcifAtomField, MmcifEntityClassifications, MmcifWriteIssue};
use crate::diagnostics::FormatError;
use crate::io::structure_interpret::{EntityKind, StructureInterpretationReport};
use crate::structure::ModelView;
use crate::topology::{
    AtomSite, Chain, Hierarchy, InstanceAtomId, MoleculeClass, MoleculeInstance,
    MoleculeInstanceId, Residue, ResidueClass, Topology,
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum AssignedEntityKind {
    Polymer,
    Branched,
    NonPolymer,
    Water,
}

impl AssignedEntityKind {
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

    fn as_source(self) -> EntityKind {
        match self {
            Self::Polymer => EntityKind::Polymer,
            Self::Branched => EntityKind::Branched,
            Self::NonPolymer => EntityKind::NonPolymer,
            Self::Water => EntityKind::Water,
        }
    }
}

/// A molecule's entity kind, and whether canonical classification derived it
/// rather than the caller or the source stating it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct MoleculeEntityKind {
    kind: AssignedEntityKind,
    canonical: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct EntityRow {
    pub(super) id: String,
    pub(super) kind: AssignedEntityKind,
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
    pub(super) kind: AssignedEntityKind,
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
    report: Option<&StructureInterpretationReport>,
) -> Result<EntityPlan, FormatError> {
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
    report: &StructureInterpretationReport,
) -> Result<EntityPlan, FormatError> {
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
                .map(|instance| (instance.instance(), instance.entity_kinds().to_vec())),
        )?;
        generic_entity_plan(model, &classifications)
    } else {
        report_entity_plan(model, report)
    }
}

fn hierarchy_asym_ids(model: ModelView<'_>) -> Result<BTreeSet<String>, FormatError> {
    let mut reserved_asym_ids = BTreeSet::new();
    for (_, chain) in model.topology().hierarchy().chains() {
        if chain.label_id().is_empty() {
            return Err(invalid_model(
                "topology hierarchy chain label ID is empty".to_owned(),
            ));
        }
        if !reserved_asym_ids.insert(chain.label_id().to_owned()) {
            return Err(fail(MmcifWriteIssue::DuplicateAsymId(
                chain.label_id().to_owned(),
            )));
        }
    }
    Ok(reserved_asym_ids)
}

pub(super) fn generic_entity_plan(
    model: ModelView<'_>,
    classifications: &BTreeMap<MoleculeInstanceId, MoleculeEntityKind>,
) -> Result<EntityPlan, FormatError> {
    let hierarchy = model.topology().hierarchy();
    let mut reserved_asym_ids = hierarchy_asym_ids(model)?;
    let mut entities = Vec::new();
    let mut asyms = Vec::new();
    let mut atoms = BTreeMap::new();

    for (_, chain) in hierarchy.chains() {
        let segments = chain_segments(hierarchy, chain, |residue, sites| {
            // Canonical classes describe molecules; a residue that a polymer
            // holds by a bond may still be an entity of another kind. Stated
            // classifications hold as given.
            let molecule_kinds = sites
                .iter()
                .filter_map(|site| classifications.get(&site.atom().instance()).copied())
                .map(|molecule| {
                    if molecule.canonical {
                        residue_entity_kind(molecule.kind, residue.class())
                    } else {
                        molecule.kind
                    }
                })
                .map(AssignedEntityKind::as_source)
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            let kind = match molecule_kinds.as_slice() {
                [kind] => entity_kind_from_source(sites[0].atom().instance(), kind)?,
                _ => {
                    return Err(fail(
                        MmcifWriteIssue::ConflictingAsymEntityClassifications {
                            asym_id: chain.label_id().to_owned(),
                            classifications: molecule_kinds,
                        },
                    ));
                }
            };
            Ok((kind, None))
        })?;
        for (asym_id, segment) in
            segment_asym_ids(chain.label_id(), segments, &mut reserved_asym_ids)
        {
            let entity_id = (entities.len() + 1).to_string();
            entities.push(EntityRow {
                id: entity_id.clone(),
                kind: segment.kind,
            });
            asyms.push(AsymRow {
                id: asym_id.clone(),
                entity_id: entity_id.clone(),
            });
            assign_sites(&mut atoms, &segment, &entity_id, &asym_id)?;
        }
    }

    for (id, molecule) in model.topology().instance_records() {
        let definition = model
            .topology()
            .instance(id)
            .expect("listed instances are topology molecules")
            .definition();
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
                return Err(fail(MmcifWriteIssue::MissingAtomSite(missing)));
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
    report: &StructureInterpretationReport,
) -> Result<EntityPlan, FormatError> {
    let mut reserved_asym_ids = hierarchy_asym_ids(model)?;
    let mut seen_instances = BTreeSet::new();
    let mut provenance = BTreeMap::new();
    let mut reserved_entity_ids = BTreeSet::new();
    for instance in report.instances() {
        if model.topology().instance(instance.instance()).is_none() {
            return Err(fail(MmcifWriteIssue::UnknownClassifiedMolecule(
                instance.instance(),
            )));
        }
        if !seen_instances.insert(instance.instance()) {
            return Err(fail(MmcifWriteIssue::DuplicateEntityClassification(
                instance.instance(),
            )));
        }
        for atom in instance.atoms() {
            if atom.atom().instance() != instance.instance()
                || model.topology().atom(atom.atom()).is_none()
            {
                return Err(fail(MmcifWriteIssue::UnknownAtomProvenance(atom.atom())));
            }
            if let Some(entity_id) = atom.entity_id() {
                reserved_entity_ids.insert(entity_id.to_owned());
            }
            if provenance.insert(atom.atom(), atom).is_some() {
                return Err(fail(MmcifWriteIssue::DuplicateAtomProvenance(atom.atom())));
            }
        }
    }

    let hierarchy = model.topology().hierarchy();
    let mut entities = Vec::new();
    let mut entity_kinds = BTreeMap::<String, AssignedEntityKind>::new();
    let mut asyms = Vec::new();
    let mut atoms = BTreeMap::new();
    let mut generated_serial = 1usize;
    for (_, chain) in hierarchy.chains() {
        let segments = chain_segments(hierarchy, chain, |_, sites| {
            let source_atoms = sites
                .iter()
                .map(|site| {
                    provenance
                        .get(&site.atom())
                        .copied()
                        .ok_or_else(|| fail(MmcifWriteIssue::MissingAtomProvenance(site.atom())))
                })
                .collect::<Result<Vec<_>, _>>()?;
            for atom in &source_atoms {
                if atom.asym_id() != chain.label_id() {
                    return Err(fail(MmcifWriteIssue::InconsistentAtomSite {
                        atom: atom.atom(),
                        field: MmcifAtomField::SourceAsymId,
                    }));
                }
            }
            let source_entity_ids = source_atoms
                .iter()
                .filter_map(|atom| atom.entity_id().map(str::to_owned))
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            if source_entity_ids.len() > 1 {
                return Err(fail(MmcifWriteIssue::ConflictingAsymEntityIds {
                    asym_id: chain.label_id().to_owned(),
                    entity_ids: source_entity_ids,
                }));
            }
            let source_kinds = source_atoms
                .iter()
                .map(|atom| atom.entity_kind.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            let kind = match source_kinds.as_slice() {
                [kind] => entity_kind_from_source(source_atoms[0].atom().instance(), kind)?,
                _ => {
                    return Err(fail(
                        MmcifWriteIssue::ConflictingAsymEntityClassifications {
                            asym_id: chain.label_id().to_owned(),
                            classifications: source_kinds,
                        },
                    ));
                }
            };
            Ok((kind, source_entity_ids.into_iter().next()))
        })?;
        for (asym_id, segment) in
            segment_asym_ids(chain.label_id(), segments, &mut reserved_asym_ids)
        {
            let entity_id = if let Some(entity_id) = &segment.source_entity_id {
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
                if *existing != segment.kind {
                    return Err(fail(
                        MmcifWriteIssue::ConflictingSourceEntityClassifications {
                            entity_id,
                            classifications: vec![existing.as_source(), segment.kind.as_source()],
                        },
                    ));
                }
            } else {
                entity_kinds.insert(entity_id.clone(), segment.kind);
                entities.push(EntityRow {
                    id: entity_id.clone(),
                    kind: segment.kind,
                });
            }
            asyms.push(AsymRow {
                id: asym_id.clone(),
                entity_id: entity_id.clone(),
            });
            assign_sites(&mut atoms, &segment, &entity_id, &asym_id)?;
        }
    }

    for atom in model.topology().atom_ids() {
        if !atoms.contains_key(atom) {
            return Err(fail(MmcifWriteIssue::MissingAtomProvenance(*atom)));
        }
    }
    for atom in provenance.keys() {
        if !atoms.contains_key(atom) {
            return Err(fail(MmcifWriteIssue::UnknownAtomProvenance(*atom)));
        }
    }

    Ok(EntityPlan {
        entities,
        asyms,
        atoms,
    })
}

/// The atom sites of one hierarchy chain that become one label asym.
struct ChainSegment<'h> {
    kind: AssignedEntityKind,
    source_entity_id: Option<String>,
    sites: Vec<&'h AtomSite>,
}

/// Splits a chain into label asyms, one per entity instance, in residue
/// order. mmCIF gives every entity instance an asym of its own, while a
/// hierarchy chain without a source `label_asym_id` may be a PDB author chain
/// holding a polymer, ligands, and water together. Residues of one entity
/// kind and source entity share an asym, except that every non-polymer
/// residue is an instance of its own. A chain the source names as one asym
/// stays one, so residues it would split are conflicting classifications.
fn chain_segments<'h>(
    hierarchy: &'h Hierarchy,
    chain: &'h Chain,
    mut classify: impl FnMut(
        &'h Residue,
        &[&'h AtomSite],
    ) -> Result<(AssignedEntityKind, Option<String>), FormatError>,
) -> Result<Vec<ChainSegment<'h>>, FormatError> {
    let mut segments = Vec::<ChainSegment<'h>>::new();
    let source_asym = chain.residues().iter().any(|residue| {
        hierarchy
            .residue(*residue)
            .expect("published hierarchy chain references a live residue")
            .atom_sites()
            .iter()
            .any(|site| {
                hierarchy
                    .atom_site(*site)
                    .expect("published hierarchy residue references a live atom site")
                    .metadata()
                    .label_asym_id
                    .is_some()
            })
    });
    for residue in chain.residues() {
        let residue = hierarchy
            .residue(*residue)
            .expect("published hierarchy chain references a live residue");
        let sites = residue
            .atom_sites()
            .iter()
            .map(|site| {
                hierarchy
                    .atom_site(*site)
                    .expect("published hierarchy residue references a live atom site")
            })
            .collect::<Vec<_>>();
        if sites.is_empty() {
            continue;
        }
        let (kind, source_entity_id) = classify(residue, &sites)?;
        let shared = (source_asym || kind != AssignedEntityKind::NonPolymer)
            .then(|| {
                segments.iter_mut().find(|segment| {
                    segment.kind == kind && segment.source_entity_id == source_entity_id
                })
            })
            .flatten();
        match shared {
            Some(segment) => segment.sites.extend(sites),
            None => segments.push(ChainSegment {
                kind,
                source_entity_id,
                sites,
            }),
        }
    }
    if source_asym && segments.len() > 1 {
        let classifications = segments
            .iter()
            .map(|segment| segment.kind.as_source())
            .collect::<BTreeSet<_>>();
        return Err(if classifications.len() > 1 {
            fail(MmcifWriteIssue::ConflictingAsymEntityClassifications {
                asym_id: chain.label_id().to_owned(),
                classifications: classifications.into_iter().collect(),
            })
        } else {
            fail(MmcifWriteIssue::ConflictingAsymEntityIds {
                asym_id: chain.label_id().to_owned(),
                entity_ids: segments
                    .iter()
                    .filter_map(|segment| segment.source_entity_id.clone())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect(),
            })
        });
    }
    Ok(segments)
}

/// Pairs each segment with its label asym ID. The first keeps the chain's
/// label; the others extend it with a serial no other asym uses.
fn segment_asym_ids<'h>(
    label: &str,
    segments: Vec<ChainSegment<'h>>,
    reserved_asym_ids: &mut BTreeSet<String>,
) -> Vec<(String, ChainSegment<'h>)> {
    let mut serial = 2usize;
    segments
        .into_iter()
        .enumerate()
        .map(|(index, segment)| {
            if index == 0 {
                return (label.to_owned(), segment);
            }
            loop {
                let candidate = format!("{label}_{serial}");
                serial += 1;
                if reserved_asym_ids.insert(candidate.clone()) {
                    return (candidate, segment);
                }
            }
        })
        .collect()
}

fn assign_sites(
    atoms: &mut BTreeMap<InstanceAtomId, AtomEntityAssignment>,
    segment: &ChainSegment<'_>,
    entity_id: &str,
    asym_id: &str,
) -> Result<(), FormatError> {
    for site in &segment.sites {
        let atom = site.atom();
        if atoms
            .insert(
                atom,
                AtomEntityAssignment {
                    entity_id: entity_id.to_owned(),
                    asym_id: asym_id.to_owned(),
                    kind: segment.kind,
                },
            )
            .is_some()
        {
            return Err(fail(MmcifWriteIssue::DuplicateAtomSite(atom)));
        }
    }
    Ok(())
}

/// The entity kind of one residue of a molecule of `molecule_kind`. Water,
/// ions, and saccharides bonded to a polymer, such as a coordinated metal or
/// an N-linked glycan, are entities of their own rather than polymer
/// residues; each saccharide residue is a non-polymer instance.
fn residue_entity_kind(
    molecule_kind: AssignedEntityKind,
    class: ResidueClass,
) -> AssignedEntityKind {
    match (molecule_kind, class) {
        (AssignedEntityKind::Polymer | AssignedEntityKind::Branched, ResidueClass::Water) => {
            AssignedEntityKind::Water
        }
        (AssignedEntityKind::Polymer, ResidueClass::Ion | ResidueClass::Carbohydrate)
        | (AssignedEntityKind::Branched, ResidueClass::Ion) => AssignedEntityKind::NonPolymer,
        (kind, _) => kind,
    }
}

fn entity_kind(
    molecule: &MoleculeInstance,
    classifications: &BTreeMap<MoleculeInstanceId, MoleculeEntityKind>,
) -> Result<AssignedEntityKind, FormatError> {
    classifications
        .get(&molecule.id())
        .map(|classification| classification.kind)
        .ok_or_else(|| fail(MmcifWriteIssue::MissingEntityClassification(molecule.id())))
}

fn entity_kind_from_source(
    molecule: MoleculeInstanceId,
    kind: &EntityKind,
) -> Result<AssignedEntityKind, FormatError> {
    match kind {
        EntityKind::Polymer => Ok(AssignedEntityKind::Polymer),
        EntityKind::Branched => Ok(AssignedEntityKind::Branched),
        EntityKind::NonPolymer => Ok(AssignedEntityKind::NonPolymer),
        EntityKind::Water => Ok(AssignedEntityKind::Water),
        EntityKind::Other(classification) => {
            Err(fail(MmcifWriteIssue::UnsupportedEntityClassification {
                molecule,
                classification: classification.clone(),
            }))
        }
    }
}

pub(super) fn normalize_entity_classifications(
    model: ModelView<'_>,
    entries: impl IntoIterator<Item = (MoleculeInstanceId, Vec<EntityKind>)>,
) -> Result<BTreeMap<MoleculeInstanceId, MoleculeEntityKind>, FormatError> {
    let mut explicit = BTreeMap::new();
    for (molecule, source_kinds) in entries {
        if model.topology().instance(molecule).is_none() {
            return Err(fail(MmcifWriteIssue::UnknownClassifiedMolecule(molecule)));
        }
        if explicit.contains_key(&molecule) {
            return Err(fail(MmcifWriteIssue::DuplicateEntityClassification(
                molecule,
            )));
        }
        let kinds = source_kinds
            .into_iter()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let kind = match kinds.as_slice() {
            [kind] => entity_kind_from_source(molecule, kind)?,
            [] => return Err(fail(MmcifWriteIssue::MissingEntityClassification(molecule))),
            _ => {
                return Err(fail(MmcifWriteIssue::ConflictingEntityClassifications {
                    molecule,
                    classifications: kinds,
                }));
            }
        };
        explicit.insert(molecule, kind);
    }
    model
        .topology()
        .instance_records()
        .map(|(molecule, _)| {
            let classification = match explicit.get(&molecule) {
                Some(&kind) => MoleculeEntityKind {
                    kind,
                    canonical: false,
                },
                None => MoleculeEntityKind {
                    kind: canonical_entity_kind(model.topology(), molecule)?,
                    canonical: true,
                },
            };
            Ok((molecule, classification))
        })
        .collect()
}

fn canonical_entity_kind(
    topology: &Topology,
    molecule: MoleculeInstanceId,
) -> Result<AssignedEntityKind, FormatError> {
    let class = topology
        .instance(molecule)
        .ok_or_else(|| invalid_molecule(molecule))?
        .class();
    match class {
        MoleculeClass::Protein | MoleculeClass::Dna | MoleculeClass::Rna => {
            Ok(AssignedEntityKind::Polymer)
        }
        MoleculeClass::Water => Ok(AssignedEntityKind::Water),
        MoleculeClass::Ion | MoleculeClass::SmallMolecule => Ok(AssignedEntityKind::NonPolymer),
        MoleculeClass::Carbohydrate => carbohydrate_entity_kind(topology, molecule),
        MoleculeClass::Other => Err(fail(
            MmcifWriteIssue::UnresolvedCanonicalEntityClassification {
                molecule,
                classification: class,
            },
        )),
    }
}

fn carbohydrate_entity_kind(
    topology: &Topology,
    molecule: MoleculeInstanceId,
) -> Result<AssignedEntityKind, FormatError> {
    let view = topology
        .instance(molecule)
        .ok_or_else(|| invalid_molecule(molecule))?;
    let mut residues = view.residues();
    if residues
        .next()
        .is_some_and(|residue| residue.class() == ResidueClass::Carbohydrate)
        && residues.next().is_none()
    {
        return Ok(AssignedEntityKind::NonPolymer);
    }
    Err(fail(
        MmcifWriteIssue::UnresolvedCanonicalEntityClassification {
            molecule,
            classification: MoleculeClass::Carbohydrate,
        },
    ))
}

fn invalid_molecule(molecule: MoleculeInstanceId) -> FormatError {
    invalid_model(format!("invalid molecule instance: {molecule}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residues_coordinated_or_linked_to_a_polymer_are_entities_of_their_own() {
        use AssignedEntityKind::{Branched, NonPolymer, Polymer, Water};
        for (molecule, class, expected) in [
            (Polymer, ResidueClass::AminoAcid, Polymer),
            (Polymer, ResidueClass::RnaNucleotide, Polymer),
            // A modified residue cannot be told from a residue's class alone.
            (Polymer, ResidueClass::Other, Polymer),
            (Polymer, ResidueClass::Ion, NonPolymer),
            (Polymer, ResidueClass::Carbohydrate, NonPolymer),
            (Polymer, ResidueClass::Water, Water),
            (Branched, ResidueClass::Carbohydrate, Branched),
            (Branched, ResidueClass::Ion, NonPolymer),
            (Branched, ResidueClass::Water, Water),
            // A metal complex keeps its coordinated water.
            (NonPolymer, ResidueClass::Water, NonPolymer),
            (NonPolymer, ResidueClass::Ion, NonPolymer),
            (Water, ResidueClass::Water, Water),
        ] {
            assert_eq!(
                residue_entity_kind(molecule, class),
                expected,
                "{molecule:?} {class:?}"
            );
        }
    }
}
