//! Entity classification and source-provenance planning before row preparation.
use super::one_based_serial;
use super::{fail, invalid_structure, MmcifEntityClassifications, MmcifWriteErrorKind};
use crate::diagnostics::{AtomField, FormatError};
use crate::io::macromolecular::{EntityKind, StructureInterpretationReport};
use crate::structure::StructureView;
use crate::topology::{
    AtomSiteView, ChainView, InstanceAtomId, MoleculeClass, MoleculeInstanceId, ResidueClass,
    ResidueId, ResidueView, Topology,
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
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
    /// The `auth_asym_id` of an atom whose site and chain state none; see
    /// [`AuthorIds`].
    pub(super) author_asym_id: String,
}

#[derive(Debug, Clone)]
pub(super) struct EntityPlan {
    pub(super) entities: Vec<EntityRow>,
    pub(super) asyms: Vec<AsymRow>,
    pub(super) atoms: BTreeMap<InstanceAtomId, AtomEntityAssignment>,
    /// The report's metal coordination between residues, written as
    /// `metalc` rows: the records that still match the structure.
    pub(super) coordination: Vec<[InstanceAtomId; 2]>,
    /// The report's coordination records left out because they no longer
    /// match the structure.
    pub(super) omitted_coordination: usize,
    /// The report's polymer sequences by entity ID: each monomer's
    /// component by its number.
    pub(super) sequences: BTreeMap<String, BTreeMap<i32, String>>,
}

pub(super) fn multi_model_entity_plan(
    structure: StructureView<'_>,
    classifications: Option<&MmcifEntityClassifications>,
    report: Option<&StructureInterpretationReport>,
) -> Result<EntityPlan, FormatError> {
    if let Some(report) = report {
        entity_plan_from_report(structure, report)
    } else {
        let normalized = normalize_entity_classifications(
            structure,
            classifications
                .into_iter()
                .flat_map(|set| set.iter().map(|(id, kind)| (id, vec![kind.clone()]))),
        )?;
        generic_entity_plan(structure, &normalized)
    }
}

pub(super) fn entity_plan_from_report(
    structure: StructureView<'_>,
    report: &StructureInterpretationReport,
) -> Result<EntityPlan, FormatError> {
    if report
        .instances()
        .iter()
        .all(|instance| instance.atoms().is_empty())
    {
        let classifications = normalize_entity_classifications(
            structure,
            report
                .instances()
                .iter()
                .map(|instance| (instance.instance(), instance.entity_kinds().to_vec())),
        )?;
        generic_entity_plan(structure, &classifications)
    } else {
        report_entity_plan(structure, report)
    }
}

/// The label asym IDs of a structure's hierarchy chains. Chains may share a
/// label, as copies
/// [`append_structure`](crate::structure::StructureEditor::append_structure)
/// appends do, while mmCIF requires every structural instance to have an ID of
/// its own. The first chain with a label keeps it; each later one is the label
/// with a copy serial, `A-2`, `A-3`, and so on, that no chain label or other
/// written asym uses. A chain without a label, as a PDB chain with a blank
/// chain column, takes the first of `A` to `Z`, `AA`, `AB`, and so on that no
/// chain label or other written asym uses. Author chain identifiers are
/// written as they are; [`AuthorIds`] names chains without one.
struct AsymIds {
    reserved: BTreeSet<String>,
    named: BTreeSet<String>,
}

impl AsymIds {
    fn new(structure: StructureView<'_>) -> Self {
        let reserved = structure
            .topology()
            .chains()
            .map(|chain| chain.label_id().to_owned())
            .filter(|label| !label.is_empty())
            .collect();
        Self {
            reserved,
            named: BTreeSet::new(),
        }
    }

    /// The asym ID of the next chain, in hierarchy order, with this label.
    fn chain(&mut self, label: &str) -> String {
        if label.is_empty() {
            return (1..)
                .map(letters)
                .find(|candidate| self.reserved.insert(candidate.clone()))
                .expect("letter IDs are unbounded");
        }
        if self.named.insert(label.to_owned()) {
            return label.to_owned();
        }
        let mut serial = 2usize;
        loop {
            let candidate = format!("{label}-{serial}");
            serial += 1;
            if self.reserved.insert(candidate.clone()) {
                return candidate;
            }
        }
    }
}

/// The author chain ID an atom site states, or else its chain: none where
/// both are absent or empty.
pub(super) fn stated_author_id<'h>(
    site: AtomSiteView<'h>,
    chain: ChainView<'h>,
) -> Option<&'h str> {
    site.metadata()
        .auth_asym_id()
        .filter(|id| !id.is_empty())
        .or_else(|| chain.author_id().filter(|id| !id.is_empty()))
}

/// The author chain IDs of chains that state none, and of molecules outside
/// the hierarchy. mmCIF readers group asyms into author chains by
/// `auth_asym_id`, so all asyms of one chain take one. A chain whose other
/// sites state one name takes it; otherwise the chain takes its own asym
/// ID, that of its first asym. Where another chain states that ID as its
/// author ID, or already takes it, the chain takes the ID with a
/// serial, `A-2`, `A-3`, and so on, that none takes, so that chains stay
/// apart.
struct AuthorIds {
    taken: BTreeSet<String>,
}

impl AuthorIds {
    fn new(structure: StructureView<'_>) -> Self {
        let taken = structure
            .topology()
            .atom_sites()
            .filter_map(|site| stated_author_id(site, site.residue().chain()))
            .map(str::to_owned)
            .collect();
        Self { taken }
    }

    /// The author ID of the chain written as asym `asym_id` whose sites
    /// state none.
    fn fallback(&mut self, asym_id: &str) -> String {
        if self.taken.insert(asym_id.to_owned()) {
            return asym_id.to_owned();
        }
        (2usize..)
            .map(|serial| format!("{asym_id}-{serial}"))
            .find(|candidate| self.taken.insert(candidate.clone()))
            .expect("serials are unbounded")
    }

    /// The author ID of the sites of `chain`, written as asym `asym_id`, that
    /// state none: the one its other sites state, else a fallback; `asym_id`
    /// itself, unclaimed, when every site states one.
    fn chain(&mut self, chain: ChainView<'_>, asym_id: &str) -> String {
        let mut unstated = false;
        let mut stated = BTreeSet::new();
        for site in chain.residues().flat_map(|residue| residue.atom_sites()) {
            match stated_author_id(site, chain) {
                Some(author) => {
                    stated.insert(author);
                }
                None => unstated = true,
            }
        }
        match (unstated, stated.first()) {
            (false, _) => asym_id.to_owned(),
            (true, Some(&author)) if stated.len() == 1 => author.to_owned(),
            (true, _) => self.fallback(asym_id),
        }
    }
}

/// What makes asyms instances of one entity: the same components, or
/// water. Each polymer and branched asym is an entity of its own.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum EntityKey {
    NonPolymer(Vec<ResidueChemistry>),
    Water,
}

impl EntityKey {
    fn of(segment: &ChainSegment) -> Option<Self> {
        match segment.kind {
            AssignedEntityKind::NonPolymer => Some(Self::NonPolymer(segment.chemistry.clone())),
            AssignedEntityKind::Water => Some(Self::Water),
            AssignedEntityKind::Polymer | AssignedEntityKind::Branched => None,
        }
    }
}

/// The `number`th letter ID, from 1: `A` to `Z`, then `AA`, `AB`, and so on.
fn letters(mut number: usize) -> String {
    let mut id = Vec::new();
    while number > 0 {
        number -= 1;
        id.push(b'A' + u8::try_from(number % 26).expect("a letter offset"));
        number /= 26;
    }
    id.reverse();
    String::from_utf8(id).expect("ASCII letters")
}

pub(super) fn generic_entity_plan(
    structure: StructureView<'_>,
    classifications: &BTreeMap<MoleculeInstanceId, MoleculeEntityKind>,
) -> Result<EntityPlan, FormatError> {
    let mut asym_ids = AsymIds::new(structure);
    let mut author_ids = AuthorIds::new(structure);
    let mut entities = Vec::<EntityRow>::new();
    let mut grouped = BTreeMap::<EntityKey, String>::new();
    let mut asyms = Vec::new();
    let mut atoms = BTreeMap::new();
    // The entity of an asym: its group's, or a new one.
    let mut entity_of = |kind: AssignedEntityKind, key: Option<EntityKey>| {
        if let Some(entity_id) = key.as_ref().and_then(|key| grouped.get(key)) {
            return entity_id.clone();
        }
        let entity_id = (entities.len() + 1).to_string();
        entities.push(EntityRow {
            id: entity_id.clone(),
            kind,
        });
        if let Some(key) = key {
            grouped.insert(key, entity_id.clone());
        }
        entity_id
    };

    for chain in structure.topology().chains() {
        let roles = residue_roles(chain);
        let segments = chain_segments(chain, |residue, sites| {
            // Canonical classes describe molecules; a residue that a polymer
            // holds by a bond may still be an entity of another kind. Stated
            // classifications hold as given.
            let role = roles[&residue.id()];
            let molecule_kinds = sites
                .iter()
                .filter_map(|atom| classifications.get(&atom.instance()).copied())
                .map(|molecule| {
                    if molecule.canonical {
                        residue_entity_kind(molecule.kind, residue.class(), role)
                    } else {
                        molecule.kind
                    }
                })
                .map(AssignedEntityKind::as_source)
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            let kind = match molecule_kinds.as_slice() {
                [kind] => entity_kind_from_source(sites[0].instance(), kind)?,
                _ => {
                    return Err(fail(
                        MmcifWriteErrorKind::ConflictingAsymEntityClassifications {
                            asym_id: chain.label_id().to_owned(),
                            classifications: molecule_kinds,
                        },
                    ));
                }
            };
            let canonical = sites.iter().all(|atom| {
                classifications
                    .get(&atom.instance())
                    .is_some_and(|molecule| molecule.canonical)
            });
            Ok(ResidueKind {
                kind,
                source_entity_id: None,
                canonical,
                glycan: classifications
                    .get(&sites[0].instance())
                    .filter(|_| canonical)
                    .and_then(|molecule| role.glycan_of(molecule.kind, kind)),
            })
        })?;
        let label = asym_ids.chain(chain.label_id());
        let author = author_ids.chain(chain, &label);
        for (asym_id, segment) in segment_asym_ids(&label, segments, &mut asym_ids.reserved) {
            let entity_id = entity_of(segment.kind, EntityKey::of(&segment));
            asyms.push(AsymRow {
                id: asym_id.clone(),
                entity_id: entity_id.clone(),
            });
            assign_sites(&mut atoms, &segment, &entity_id, &asym_id, &author)?;
        }
    }

    for (id, molecule_atoms) in unplaced_molecules(structure, &atoms)? {
        let asym_id = molecule_asym_id(id, &mut asym_ids.reserved);
        let kind = entity_kind(id, classifications)?;
        // Molecules outside the hierarchy are written as the generic
        // component `MOL`, which names no chemistry, so only water groups.
        let entity_id = entity_of(
            kind,
            (kind == AssignedEntityKind::Water).then_some(EntityKey::Water),
        );
        let author = author_ids.fallback(&asym_id);
        asyms.push(AsymRow {
            id: asym_id.clone(),
            entity_id: entity_id.clone(),
        });
        assign_atoms(
            &mut atoms,
            molecule_atoms,
            &entity_id,
            &asym_id,
            kind,
            &author,
        );
    }

    Ok(EntityPlan {
        entities,
        asyms,
        atoms,
        coordination: Vec::new(),
        omitted_coordination: 0,
        sequences: BTreeMap::new(),
    })
}

fn report_entity_plan(
    structure: StructureView<'_>,
    report: &StructureInterpretationReport,
) -> Result<EntityPlan, FormatError> {
    let mut asym_ids = AsymIds::new(structure);
    let mut author_ids = AuthorIds::new(structure);
    let mut seen_instances = BTreeSet::new();
    let mut provenance = BTreeMap::new();
    let mut reserved_entity_ids = BTreeSet::new();
    for instance in report.instances() {
        if structure.topology().instance(instance.instance()).is_none() {
            return Err(fail(MmcifWriteErrorKind::UnknownClassifiedMolecule(
                instance.instance(),
            )));
        }
        if !seen_instances.insert(instance.instance()) {
            return Err(fail(MmcifWriteErrorKind::DuplicateEntityClassification(
                instance.instance(),
            )));
        }
        for atom in instance.atoms() {
            if atom.atom().instance() != instance.instance()
                || structure.topology().atom(atom.atom()).is_none()
            {
                return Err(fail(MmcifWriteErrorKind::UnknownAtomProvenance(
                    atom.atom(),
                )));
            }
            if let Some(entity_id) = atom.entity_id() {
                reserved_entity_ids.insert(entity_id.to_owned());
            }
            if provenance.insert(atom.atom(), atom).is_some() {
                return Err(fail(MmcifWriteErrorKind::DuplicateAtomProvenance(
                    atom.atom(),
                )));
            }
        }
    }

    let mut chains = Vec::new();
    for chain in structure.topology().chains() {
        let roles = residue_roles(chain);
        let segments = chain_segments(chain, |residue, sites| {
            let source_atoms = sites
                .iter()
                .filter_map(|atom| provenance.get(atom).copied())
                .collect::<Vec<_>>();
            // Atoms an edit added have no provenance. A residue of such
            // atoms alone takes the kind the structure gives it, as without
            // a report; another takes its source atoms' kind and entity.
            if source_atoms.is_empty() {
                let role = roles[&residue.id()];
                let kind = canonical_residue_kind(structure, chain, residue, sites, role)?;
                return Ok(ResidueKind {
                    kind,
                    source_entity_id: None,
                    canonical: true,
                    glycan: role.glycan_of(
                        canonical_entity_kind(structure.topology(), sites[0].instance())?,
                        kind,
                    ),
                });
            }
            for atom in &source_atoms {
                if atom.asym_id() != chain.label_id() {
                    return Err(fail(MmcifWriteErrorKind::InconsistentAtomSite {
                        atom: atom.atom(),
                        field: AtomField::SourceAsymId,
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
                return Err(fail(MmcifWriteErrorKind::ConflictingAsymEntityIds {
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
                        MmcifWriteErrorKind::ConflictingAsymEntityClassifications {
                            asym_id: chain.label_id().to_owned(),
                            classifications: source_kinds,
                        },
                    ));
                }
            };
            Ok(ResidueKind {
                kind,
                source_entity_id: source_entity_ids.into_iter().next(),
                canonical: false,
                glycan: None,
            })
        })?;
        let label = asym_ids.chain(chain.label_id());
        let author = author_ids.chain(chain, &label);
        // The report's sequence of the polymer on this chain.
        let sequence = report
            .polymer_entities()
            .iter()
            .find(|entity| {
                entity
                    .instances()
                    .iter()
                    .any(|instance| instance.chain() == Some(chain.id()))
            })
            .map(|entity| {
                entity
                    .monomers()
                    .iter()
                    .map(|monomer| (monomer.number(), monomer.component_id().to_owned()))
                    .collect::<BTreeMap<_, _>>()
            });
        let segments = segment_asym_ids(&label, segments, &mut asym_ids.reserved);
        chains.push((author, sequence, segments));
    }

    // Entities by what their asyms share: a source entity takes the asyms an
    // edit added of its chemistry.
    let mut grouped = BTreeMap::<EntityKey, String>::new();
    for (_, segment) in chains.iter().flat_map(|(_, _, segments)| segments) {
        if let (Some(entity_id), Some(key)) = (&segment.source_entity_id, EntityKey::of(segment)) {
            grouped.entry(key).or_insert_with(|| entity_id.clone());
        }
    }
    let mut sequences = BTreeMap::new();
    let mut entities = Vec::new();
    let mut entity_kinds = BTreeMap::<String, AssignedEntityKind>::new();
    let mut generated_serial = 1usize;
    // The entity of an asym: its source entity, its group's, or a new one.
    let mut entity_of = |kind: AssignedEntityKind,
                         source_entity_id: Option<&String>,
                         key: Option<EntityKey>|
     -> Result<String, FormatError> {
        let entity_id = if let Some(entity_id) = source_entity_id {
            entity_id.clone()
        } else if let Some(entity_id) = key.as_ref().and_then(|key| grouped.get(key)) {
            entity_id.clone()
        } else {
            let entity_id = loop {
                let candidate = format!("K{generated_serial}");
                generated_serial += 1;
                if reserved_entity_ids.insert(candidate.clone()) {
                    break candidate;
                }
            };
            if let Some(key) = key {
                grouped.insert(key, entity_id.clone());
            }
            entity_id
        };
        if let Some(existing) = entity_kinds.get(&entity_id) {
            if *existing != kind {
                return Err(fail(
                    MmcifWriteErrorKind::ConflictingSourceEntityClassifications {
                        entity_id,
                        classifications: vec![existing.as_source(), kind.as_source()],
                    },
                ));
            }
        } else {
            entity_kinds.insert(entity_id.clone(), kind);
            entities.push(EntityRow {
                id: entity_id.clone(),
                kind,
            });
        }
        Ok(entity_id)
    };
    let mut asyms = Vec::new();
    let mut atoms = BTreeMap::new();
    for (author, sequence, segments) in chains {
        for (asym_id, segment) in segments {
            let entity_id = entity_of(
                segment.kind,
                segment.source_entity_id.as_ref(),
                EntityKey::of(&segment),
            )?;
            asyms.push(AsymRow {
                id: asym_id.clone(),
                entity_id: entity_id.clone(),
            });
            // The sequence is the source polymer's. A polymer an edit added
            // to the chain has no provenance, and its own residues number it.
            let source_polymer = segment.kind == AssignedEntityKind::Polymer
                && segment
                    .sites
                    .iter()
                    .any(|atom| provenance.contains_key(atom));
            if let (true, Some(sequence)) = (source_polymer, &sequence) {
                sequences
                    .entry(entity_id.clone())
                    .or_insert_with(|| sequence.clone());
            }
            assign_sites(&mut atoms, &segment, &entity_id, &asym_id, &author)?;
        }
    }

    // Molecules an edit added outside the hierarchy are written as without
    // a report.
    for (id, molecule_atoms) in unplaced_molecules(structure, &atoms)? {
        if let Some(&atom) = molecule_atoms
            .iter()
            .find(|atom| provenance.contains_key(atom))
        {
            return Err(fail(MmcifWriteErrorKind::MissingAtomSite(atom)));
        }
        let asym_id = molecule_asym_id(id, &mut asym_ids.reserved);
        let kind = canonical_entity_kind(structure.topology(), id)?;
        let key = (kind == AssignedEntityKind::Water).then_some(EntityKey::Water);
        let entity_id = entity_of(kind, None, key)?;
        let author = author_ids.fallback(&asym_id);
        asyms.push(AsymRow {
            id: asym_id.clone(),
            entity_id: entity_id.clone(),
        });
        assign_atoms(
            &mut atoms,
            molecule_atoms,
            &entity_id,
            &asym_id,
            kind,
            &author,
        );
    }
    for atom in provenance.keys() {
        if !atoms.contains_key(atom) {
            return Err(fail(MmcifWriteErrorKind::UnknownAtomProvenance(*atom)));
        }
    }

    let coordination = report.coordination_in(structure.topology());

    Ok(EntityPlan {
        entities,
        asyms,
        atoms,
        coordination: coordination.records,
        omitted_coordination: coordination.unmatched,
        sequences,
    })
}

/// The atoms of the atom sites of one hierarchy chain that become one label
/// asym.
struct ChainSegment {
    kind: AssignedEntityKind,
    source_entity_id: Option<String>,
    /// The molecules its atoms lie in.
    molecules: BTreeSet<MoleculeInstanceId>,
    /// The places in a sequence its residues state.
    places: BTreeSet<i32>,
    /// The first residue of the glycan it is, when it is one that a
    /// polymer holds.
    glycan: Option<ResidueId>,
    /// The chemistry of its non-polymer residues, in residue order.
    chemistry: Vec<ResidueChemistry>,
    /// Its residues' author numbers and insertion codes, which identify
    /// residues that have no label number.
    author_keys: BTreeSet<(Option<String>, Option<String>)>,
    sites: Vec<InstanceAtomId>,
}

impl ChainSegment {
    fn new(entry: ResidueEntry) -> Self {
        let mut segment = Self {
            kind: entry.kind,
            source_entity_id: entry.source_entity_id.clone(),
            molecules: BTreeSet::new(),
            places: BTreeSet::new(),
            glycan: entry.glycan,
            chemistry: Vec::new(),
            author_keys: BTreeSet::new(),
            sites: Vec::new(),
        };
        segment.push(entry);
        segment
    }

    /// Whether it is a polymer whose residues state their places in its
    /// sequence, and the place `entry` states is free.
    fn places_free(&self, entry: &ResidueEntry) -> bool {
        self.kind == AssignedEntityKind::Polymer
            && !self.places.is_empty()
            && entry
                .place
                .is_some_and(|place| !self.places.contains(&place))
    }

    fn push(&mut self, entry: ResidueEntry) {
        self.molecules
            .extend(entry.sites.iter().map(|atom| atom.instance()));
        self.places.extend(entry.place);
        self.chemistry.extend(entry.chemistry);
        self.author_keys.insert(entry.author_key);
        self.sites.extend(entry.sites);
    }
}

/// A component and, for a placeholder name, its atoms' chemistry, which
/// decide whether two non-polymer residues are instances of one entity.
type ResidueChemistry = (String, Vec<AtomChemistry>);

/// An atom's name, atomic number, mass number, charge, and hydrogens, and
/// the names of its bonded atoms in the residue with the bond orders, or
/// the residue's other atoms when it states no name.
type AtomChemistry = (
    String,
    u8,
    Option<u16>,
    i8,
    Option<usize>,
    Vec<(String, String)>,
);

/// A residue's component and, when its name is a placeholder, the
/// chemistry of its atoms by name. A named component is one chemistry, as
/// the chemical component dictionary defines it, however much of it a
/// residue models.
fn residue_chemistry(residue: ResidueView<'_>) -> ResidueChemistry {
    let component = residue.label_comp_id().unwrap_or(residue.name());
    if !placeholder_component(component) {
        return (component.to_owned(), Vec::new());
    }
    let name = |site: AtomSiteView<'_>| site_name(site).to_owned();
    let mut atoms = residue
        .atom_sites()
        .map(|site| {
            let atom = site.atom();
            let mut bonds = atom
                .bonds()
                .map(|bond| {
                    let [left, right] = bond.atoms();
                    let other = if left.id() == atom.id() { right } else { left };
                    let other = other
                        .atom_site()
                        .filter(|other| other.residue().id() == residue.id())
                        .map_or_else(|| other.atom().element.symbol().to_owned(), name);
                    (other, format!("{:?}", bond.bond().order))
                })
                .collect::<Vec<_>>();
            bonds.sort();
            (
                name(site),
                atom.atom().element.atomic_number(),
                atom.atom().isotope.map(std::num::NonZero::get),
                atom.atom().formal_charge,
                atom.total_hydrogens(),
                bonds,
            )
        })
        .collect::<Vec<_>>();
    atoms.sort();
    (component.to_owned(), atoms)
}

/// Whether a component name names no one chemistry: no name, the codes the
/// wwPDB reserves for unknown and in-house ligands (`UNL`, `UNX`, `UNK`,
/// `LIG`, `DRG`, `INH`, and `01` to `99`), or `MOL`, the component of a
/// molecule outside the hierarchy.
fn placeholder_component(name: &str) -> bool {
    matches!(
        name,
        "" | "UNL" | "UNX" | "UNK" | "LIG" | "DRG" | "INH" | "MOL"
    ) || (name.len() == 2 && name.bytes().all(|byte| byte.is_ascii_digit()) && name != "00")
}

/// An atom site's name: its label, or else its author name.
fn site_name(site: AtomSiteView<'_>) -> &str {
    site.metadata()
        .label_atom_id()
        .or(site.metadata().auth_atom_id())
        .unwrap_or_default()
}

/// The atoms of a polymer's monomers that its backbone links join.
const BACKBONE_ATOMS: [&str; 4] = ["N", "C", "P", "O3'"];

/// Where canonical classification places a residue of a chain.
#[derive(Debug, Clone, Copy)]
struct ResidueRole {
    /// Whether it is a monomer of the chain's sequence.
    monomer: bool,
    /// The saccharides bonded to it and to each other, with it: their
    /// first residue and their number.
    glycan: Option<(ResidueId, usize)>,
}

impl ResidueRole {
    /// The glycan a residue of `kind` in a molecule of `molecule_kind` is
    /// an instance of: one a polymer holds.
    fn glycan_of(
        self,
        molecule_kind: AssignedEntityKind,
        kind: AssignedEntityKind,
    ) -> Option<ResidueId> {
        self.glycan
            .filter(|_| {
                molecule_kind == AssignedEntityKind::Polymer && kind == AssignedEntityKind::Branched
            })
            .map(|(first, _)| first)
    }
}

/// The roles of a chain's residues. A residue is a monomer of the chain's
/// sequence when it states its place (`label_seq_id`) or is an amino acid
/// or a nucleotide. Where no residue of the chain states a place, as in a
/// PDB entry without `SEQRES` records, a residue of another class is a
/// monomer too when a backbone atom of a monomer (`N`, `C`, `P`, `O3'`)
/// bonds it, as a modified residue's neighbours do. Saccharides bonded to
/// each other are one glycan.
fn residue_roles(chain: ChainView<'_>) -> BTreeMap<ResidueId, ResidueRole> {
    let residues = chain.residues().collect::<Vec<_>>();
    let classes = residues
        .iter()
        .map(|residue| (residue.id(), residue.class()))
        .collect::<BTreeMap<_, _>>();
    // Each residue's bonds to the chain's other residues: the other residue
    // and the name of its atom.
    let links = residues
        .iter()
        .map(|residue| {
            let mut links = Vec::new();
            for site in residue.atom_sites() {
                let atom = site.atom();
                for bond in atom.bonds() {
                    let [left, right] = bond.atoms();
                    let other = if left.id() == atom.id() { right } else { left };
                    if let Some(other) = other.atom_site().filter(|other| {
                        other.residue().id() != residue.id()
                            && classes.contains_key(&other.residue().id())
                    }) {
                        links.push((other.residue().id(), site_name(other)));
                    }
                }
            }
            (residue.id(), links)
        })
        .collect::<BTreeMap<_, _>>();
    let mut monomers = residues
        .iter()
        .filter(|residue| residue.label_seq_id().is_some() || polymer_class(residue.class()))
        .map(|residue| residue.id())
        .collect::<BTreeSet<_>>();
    if residues
        .iter()
        .all(|residue| residue.label_seq_id().is_none())
    {
        loop {
            let joined = residues
                .iter()
                .filter(|residue| {
                    residue.class() == ResidueClass::Other && !monomers.contains(&residue.id())
                })
                .filter(|residue| {
                    links[&residue.id()].iter().any(|(other, name)| {
                        monomers.contains(other) && BACKBONE_ATOMS.contains(name)
                    })
                })
                .map(|residue| residue.id())
                .collect::<Vec<_>>();
            if joined.is_empty() {
                break;
            }
            monomers.extend(joined);
        }
    }
    let saccharide = |residue: &ResidueId| classes[residue] == ResidueClass::Carbohydrate;
    let mut glycans = BTreeMap::new();
    for residue in residues.iter().map(|residue| residue.id()) {
        if !saccharide(&residue) || glycans.contains_key(&residue) {
            continue;
        }
        let mut members = vec![residue];
        let mut next = 0;
        while let Some(&at) = members.get(next) {
            next += 1;
            for (other, _) in &links[&at] {
                if saccharide(other) && !members.contains(other) {
                    members.push(*other);
                }
            }
        }
        for &member in &members {
            glycans.insert(member, (residue, members.len()));
        }
    }
    residues
        .iter()
        .map(|residue| {
            let id = residue.id();
            let role = ResidueRole {
                monomer: monomers.contains(&id),
                glycan: glycans.get(&id).copied(),
            };
            (id, role)
        })
        .collect()
}

/// Whether residues of `class` are monomers of a polymer.
fn polymer_class(class: ResidueClass) -> bool {
    matches!(
        class,
        ResidueClass::AminoAcid | ResidueClass::DnaNucleotide | ResidueClass::RnaNucleotide
    )
}

/// A residue's entity kind and source entity.
struct ResidueKind {
    kind: AssignedEntityKind,
    source_entity_id: Option<String>,
    /// Whether canonical classification gave the kind, rather than the
    /// caller or the source.
    canonical: bool,
    /// The first residue of the glycan it is part of, when canonical
    /// classification makes it branched.
    glycan: Option<ResidueId>,
}

/// One residue's entity kind and atoms, before it joins a segment.
struct ResidueEntry {
    kind: AssignedEntityKind,
    source_entity_id: Option<String>,
    class: ResidueClass,
    /// Whether canonical classification gave the kind, rather than the
    /// caller or the source.
    canonical: bool,
    /// Its stated place in a sequence (`label_seq_id`).
    place: Option<i32>,
    glycan: Option<ResidueId>,
    chemistry: Option<ResidueChemistry>,
    author_key: (Option<String>, Option<String>),
    sites: Vec<InstanceAtomId>,
}

impl ResidueEntry {
    /// Whether it is a monomer of a polymer that canonical classification
    /// makes no polymer, as one a chain break leaves alone: a residue of a
    /// polymer's residue class that states its place in the sequence.
    fn fragment(&self) -> bool {
        self.kind == AssignedEntityKind::NonPolymer
            && self.canonical
            && self.place.is_some()
            && polymer_class(self.class)
    }
}

/// Splits a chain into label asyms, one per entity instance, with the asym
/// the source names first. mmCIF gives every entity instance an asym of its
/// own, while a hierarchy chain may be a PDB author chain holding a polymer,
/// ligands, and water together, or an asym an edit added residues to.
///
/// The residues whose sites state a `label_asym_id` are the source's asym,
/// which stays one; residues of other kinds or source entities in it are
/// conflicting classifications, except that a residue of a polymer's
/// residue class that states its place in the sequence and that canonical
/// classification makes no polymer, as one a chain break leaves alone, is
/// a polymer residue. Of the other residues, each non-polymer residue is an
/// instance of its own; a polymer residue joins the asym of its molecule,
/// so that the sequence of one molecule does not continue into the next,
/// or the polymer asym whose residues state their places when it states a
/// free one, as do the pieces of a chain a PDB reader numbers by its
/// `SEQRES` sequence and such a fragment; a branched residue joins the asym
/// of its glycan; and water joins the chain's water asym, or the next one,
/// so that no two waters of an asym share an author number.
fn chain_segments<'h>(
    chain: ChainView<'h>,
    mut classify: impl FnMut(ResidueView<'h>, &[InstanceAtomId]) -> Result<ResidueKind, FormatError>,
) -> Result<Vec<ChainSegment>, FormatError> {
    let mut sourced = Vec::new();
    let mut added = Vec::new();
    for residue in chain.residues() {
        let sites = residue
            .atom_sites()
            .map(|site| site.atom().id())
            .collect::<Vec<_>>();
        if sites.is_empty() {
            continue;
        }
        let ResidueKind {
            kind,
            source_entity_id,
            canonical,
            glycan,
        } = classify(residue, &sites)?;
        let entry = ResidueEntry {
            kind,
            source_entity_id,
            class: residue.class(),
            canonical,
            place: residue.label_seq_id(),
            glycan,
            chemistry: (kind == AssignedEntityKind::NonPolymer).then(|| residue_chemistry(residue)),
            author_key: (
                residue.author_seq_id().map(str::to_owned),
                residue.insertion_code().map(str::to_owned),
            ),
            sites,
        };
        if residue
            .atom_sites()
            .any(|site| site.metadata().label_asym_id.is_some())
        {
            sourced.push(entry);
        } else {
            added.push(entry);
        }
    }
    let mut segments = source_segment(chain, sourced)?
        .into_iter()
        .collect::<Vec<_>>();
    // Fragments join a polymer once the chain's polymers have their asyms.
    let (fragments, others): (Vec<_>, Vec<_>) = added.into_iter().partition(ResidueEntry::fragment);
    for entry in others.into_iter().chain(fragments) {
        let shared = match entry.kind {
            AssignedEntityKind::NonPolymer if entry.fragment() => segments
                .iter_mut()
                .find(|segment| segment.places_free(&entry)),
            AssignedEntityKind::NonPolymer => None,
            AssignedEntityKind::Water => segments.iter_mut().find(|segment| {
                segment.kind == AssignedEntityKind::Water
                    && !segment.author_keys.contains(&entry.author_key)
            }),
            AssignedEntityKind::Polymer => segments.iter_mut().find(|segment| {
                segment.kind == AssignedEntityKind::Polymer
                    && (segment.molecules.contains(&entry.sites[0].instance())
                        || segment.places_free(&entry))
            }),
            AssignedEntityKind::Branched => segments.iter_mut().find(|segment| {
                segment.kind == AssignedEntityKind::Branched
                    && segment.glycan == entry.glycan
                    && segment.molecules.contains(&entry.sites[0].instance())
            }),
        };
        match shared {
            Some(segment) => segment.push(entry),
            None => segments.push(ChainSegment::new(entry)),
        }
    }
    Ok(segments)
}

/// The segment of the residues of a chain whose sites state a
/// `label_asym_id`: the source's asym.
fn source_segment(
    chain: ChainView<'_>,
    entries: Vec<ResidueEntry>,
) -> Result<Option<ChainSegment>, FormatError> {
    let kinds = entries
        .iter()
        .map(|entry| entry.kind)
        .collect::<BTreeSet<_>>();
    let kind = match kinds.iter().copied().collect::<Vec<_>>().as_slice() {
        [] => return Ok(None),
        [kind] => *kind,
        [AssignedEntityKind::Polymer, AssignedEntityKind::NonPolymer]
            if entries
                .iter()
                .all(|entry| entry.kind == AssignedEntityKind::Polymer || entry.fragment()) =>
        {
            AssignedEntityKind::Polymer
        }
        _ => {
            return Err(fail(
                MmcifWriteErrorKind::ConflictingAsymEntityClassifications {
                    asym_id: chain.label_id().to_owned(),
                    classifications: kinds
                        .into_iter()
                        .map(AssignedEntityKind::as_source)
                        .collect(),
                },
            ));
        }
    };
    let entity_ids = entries
        .iter()
        .filter_map(|entry| entry.source_entity_id.clone())
        .collect::<BTreeSet<_>>();
    if entity_ids.len() > 1 {
        return Err(fail(MmcifWriteErrorKind::ConflictingAsymEntityIds {
            asym_id: chain.label_id().to_owned(),
            entity_ids: entity_ids.into_iter().collect(),
        }));
    }
    let mut entries = entries.into_iter();
    let mut segment = ChainSegment::new(entries.next().expect("one entry per kind"));
    segment.kind = kind;
    segment.source_entity_id = entity_ids.into_iter().next();
    for entry in entries {
        segment.push(entry);
    }
    Ok(Some(segment))
}

/// Pairs each segment with its label asym ID. The first keeps the chain's
/// asym ID; the others extend it with a serial no other asym uses.
fn segment_asym_ids(
    label: &str,
    segments: Vec<ChainSegment>,
    reserved_asym_ids: &mut BTreeSet<String>,
) -> Vec<(String, ChainSegment)> {
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
    segment: &ChainSegment,
    entity_id: &str,
    asym_id: &str,
    author_asym_id: &str,
) -> Result<(), FormatError> {
    for &atom in &segment.sites {
        if atoms
            .insert(
                atom,
                AtomEntityAssignment {
                    entity_id: entity_id.to_owned(),
                    asym_id: asym_id.to_owned(),
                    kind: segment.kind,
                    author_asym_id: author_asym_id.to_owned(),
                },
            )
            .is_some()
        {
            return Err(fail(MmcifWriteErrorKind::DuplicateAtomSite(atom)));
        }
    }
    Ok(())
}

/// Assigns every atom of a molecule outside the hierarchy.
fn assign_atoms(
    atoms: &mut BTreeMap<InstanceAtomId, AtomEntityAssignment>,
    molecule_atoms: Vec<InstanceAtomId>,
    entity_id: &str,
    asym_id: &str,
    kind: AssignedEntityKind,
    author_asym_id: &str,
) {
    for atom in molecule_atoms {
        atoms.insert(
            atom,
            AtomEntityAssignment {
                entity_id: entity_id.to_owned(),
                asym_id: asym_id.to_owned(),
                kind,
                author_asym_id: author_asym_id.to_owned(),
            },
        );
    }
}

/// The molecules none of whose atoms the chains assign, with their atoms:
/// molecules outside the hierarchy. A molecule the chains assign in part
/// fails with the first atom without an atom site.
fn unplaced_molecules(
    structure: StructureView<'_>,
    atoms: &BTreeMap<InstanceAtomId, AtomEntityAssignment>,
) -> Result<Vec<(MoleculeInstanceId, Vec<InstanceAtomId>)>, FormatError> {
    let mut unplaced = Vec::new();
    for (id, molecule) in structure.topology().instance_records() {
        let molecule_atoms = structure
            .topology()
            .instance(id)
            .expect("listed instances are topology molecules")
            .definition()
            .molecule()
            .atoms()
            .map(|(atom, _)| molecule.qualify_atom(atom))
            .collect::<Vec<_>>();
        let assigned = molecule_atoms
            .iter()
            .filter(|atom| atoms.contains_key(atom))
            .count();
        if assigned == 0 {
            unplaced.push((id, molecule_atoms));
        } else if assigned != molecule_atoms.len() {
            let missing = molecule_atoms
                .into_iter()
                .find(|atom| !atoms.contains_key(atom))
                .expect("partially assigned molecule has a missing atom site");
            return Err(fail(MmcifWriteErrorKind::MissingAtomSite(missing)));
        }
    }
    Ok(unplaced)
}

/// The asym ID of a molecule outside the hierarchy: `M` and its 1-based
/// instance number, with a serial where a chain or asym uses that.
fn molecule_asym_id(id: MoleculeInstanceId, reserved: &mut BTreeSet<String>) -> String {
    let base = format!("M{}", one_based_serial(id.raw()));
    let mut candidate = base.clone();
    let mut suffix = 2usize;
    while reserved.contains(&candidate) {
        candidate = format!("{base}_{suffix}");
        suffix += 1;
    }
    reserved.insert(candidate.clone());
    candidate
}

/// The entity kind of a residue from its molecules' canonical classes, as a
/// write without a report gives it.
fn canonical_residue_kind(
    structure: StructureView<'_>,
    chain: ChainView<'_>,
    residue: ResidueView<'_>,
    sites: &[InstanceAtomId],
    role: ResidueRole,
) -> Result<AssignedEntityKind, FormatError> {
    let kinds = sites
        .iter()
        .map(|atom| {
            canonical_entity_kind(structure.topology(), atom.instance())
                .map(|kind| residue_entity_kind(kind, residue.class(), role).as_source())
        })
        .collect::<Result<BTreeSet<_>, _>>()?
        .into_iter()
        .collect::<Vec<_>>();
    match kinds.as_slice() {
        [kind] => entity_kind_from_source(sites[0].instance(), kind),
        _ => Err(fail(
            MmcifWriteErrorKind::ConflictingAsymEntityClassifications {
                asym_id: chain.label_id().to_owned(),
                classifications: kinds,
            },
        )),
    }
}

/// The entity kind of one residue of a molecule of `molecule_kind`, in
/// `role`. Only the monomers of a polymer's sequence are polymer residues.
/// Water, ions, and other residues a polymer holds by a bond, such as a
/// coordinated metal or a covalently bound ligand, are entities of their
/// own; a glycan bonded to a polymer is a branched entity, or a non-polymer
/// one when it is a single saccharide.
fn residue_entity_kind(
    molecule_kind: AssignedEntityKind,
    class: ResidueClass,
    role: ResidueRole,
) -> AssignedEntityKind {
    match (molecule_kind, class) {
        (AssignedEntityKind::Polymer | AssignedEntityKind::Branched, ResidueClass::Water) => {
            AssignedEntityKind::Water
        }
        (AssignedEntityKind::Polymer | AssignedEntityKind::Branched, ResidueClass::Ion) => {
            AssignedEntityKind::NonPolymer
        }
        (AssignedEntityKind::Polymer, ResidueClass::Carbohydrate) => match role.glycan {
            Some((_, residues)) if residues > 1 => AssignedEntityKind::Branched,
            _ => AssignedEntityKind::NonPolymer,
        },
        (AssignedEntityKind::Polymer, _) if !role.monomer => AssignedEntityKind::NonPolymer,
        (kind, _) => kind,
    }
}

fn entity_kind(
    molecule: MoleculeInstanceId,
    classifications: &BTreeMap<MoleculeInstanceId, MoleculeEntityKind>,
) -> Result<AssignedEntityKind, FormatError> {
    classifications
        .get(&molecule)
        .map(|classification| classification.kind)
        .ok_or_else(|| fail(MmcifWriteErrorKind::MissingEntityClassification(molecule)))
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
            Err(fail(MmcifWriteErrorKind::UnsupportedEntityClassification {
                molecule,
                classification: classification.clone(),
            }))
        }
    }
}

pub(super) fn normalize_entity_classifications(
    structure: StructureView<'_>,
    entries: impl IntoIterator<Item = (MoleculeInstanceId, Vec<EntityKind>)>,
) -> Result<BTreeMap<MoleculeInstanceId, MoleculeEntityKind>, FormatError> {
    let mut explicit = BTreeMap::new();
    for (molecule, source_kinds) in entries {
        if structure.topology().instance(molecule).is_none() {
            return Err(fail(MmcifWriteErrorKind::UnknownClassifiedMolecule(
                molecule,
            )));
        }
        if explicit.contains_key(&molecule) {
            return Err(fail(MmcifWriteErrorKind::DuplicateEntityClassification(
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
            [] => {
                return Err(fail(MmcifWriteErrorKind::MissingEntityClassification(
                    molecule,
                )))
            }
            _ => {
                return Err(fail(
                    MmcifWriteErrorKind::ConflictingEntityClassifications {
                        molecule,
                        classifications: kinds,
                    },
                ));
            }
        };
        explicit.insert(molecule, kind);
    }
    structure
        .topology()
        .instance_records()
        .map(|(molecule, _)| {
            let classification = match explicit.get(&molecule) {
                Some(&kind) => MoleculeEntityKind {
                    kind,
                    canonical: false,
                },
                None => MoleculeEntityKind {
                    kind: canonical_entity_kind(structure.topology(), molecule)?,
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
            MmcifWriteErrorKind::UnresolvedCanonicalEntityClassification {
                molecule,
                classification: class,
            },
        )),
    }
}

/// A saccharide of one residue is a non-polymer entity, and one of several
/// residues a branched entity.
fn carbohydrate_entity_kind(
    topology: &Topology,
    molecule: MoleculeInstanceId,
) -> Result<AssignedEntityKind, FormatError> {
    let view = topology
        .instance(molecule)
        .ok_or_else(|| invalid_molecule(molecule))?;
    Ok(match view.residues().count() {
        0 | 1 => AssignedEntityKind::NonPolymer,
        _ => AssignedEntityKind::Branched,
    })
}

fn invalid_molecule(molecule: MoleculeInstanceId) -> FormatError {
    invalid_structure(format!("The structure has no {molecule}."))
}

#[cfg(test)]
mod tests;
