//! Polymer entity sequences an mmCIF block states: `_entity_poly`,
//! `_entity_poly_seq`, and `_pdbx_poly_seq_scheme`, with the block's own
//! `_chem_comp` codes.

use std::collections::{BTreeMap, HashMap, HashSet};

use crate::diagnostics::{DiagnosticKind, FormatError};
use crate::io::{ComponentLibrary, ComponentLinking};
use crate::sequence::{CodeOrigin, CodeSources, ResidueCode, StatedComponent};
use crate::topology::{ChainId, ResidueId, Topology};

use super::super::components::parent_ids;
use super::super::mmcif_category::MmcifCategory;
use super::super::MmcifBlock;
use super::atom_site::{optional, optional_i32, required, row_error, AtomRow, Field, FieldKey};

/// The type of a polymer entity, its `_entity_poly.type`.
///
/// Matching of source values is case-insensitive.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum PolymerType {
    /// `polypeptide(L)`.
    PolypeptideL,
    /// `polypeptide(D)`.
    PolypeptideD,
    /// `polydeoxyribonucleotide`: DNA.
    Polydeoxyribonucleotide,
    /// `polyribonucleotide`: RNA.
    Polyribonucleotide,
    /// `polydeoxyribonucleotide/polyribonucleotide hybrid`.
    NucleotideHybrid,
    /// `peptide nucleic acid`.
    PeptideNucleicAcid,
    /// `cyclic-pseudo-peptide`.
    CyclicPseudoPeptide,
    /// `polysaccharide(D)`.
    PolysaccharideD,
    /// `polysaccharide(L)`.
    PolysaccharideL,
    /// Any other value, such as `other`, as written in the source.
    Other(String),
}

impl PolymerType {
    fn from_mmcif(value: &str) -> Self {
        match value.to_ascii_lowercase().as_str() {
            "polypeptide(l)" => Self::PolypeptideL,
            "polypeptide(d)" => Self::PolypeptideD,
            "polydeoxyribonucleotide" => Self::Polydeoxyribonucleotide,
            "polyribonucleotide" => Self::Polyribonucleotide,
            "polydeoxyribonucleotide/polyribonucleotide hybrid" => Self::NucleotideHybrid,
            "peptide nucleic acid" => Self::PeptideNucleicAcid,
            "cyclic-pseudo-peptide" => Self::CyclicPseudoPeptide,
            "polysaccharide(d)" => Self::PolysaccharideD,
            "polysaccharide(l)" => Self::PolysaccharideL,
            _ => Self::Other(value.to_owned()),
        }
    }

    /// Whether the entity is DNA, RNA, or a hybrid of the two.
    pub const fn is_nucleic_acid(&self) -> bool {
        matches!(
            self,
            Self::Polydeoxyribonucleotide | Self::Polyribonucleotide | Self::NucleotideHybrid
        )
    }
}

/// Whether the file states a part of a polymer entity or Kekule derives it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SequenceSource {
    /// The file states it.
    Stated,
    /// Derived: monomers from the atom sites, a sequence from the monomers'
    /// codes.
    Derived,
}

/// One polymer entity of an mmCIF block: its monomers, its canonical
/// one-letter sequence, and its instances in the interpreted structure.
///
/// Monomers come from `_entity_poly_seq`, which lists the whole polymer,
/// residues the structure does not observe included. Without it, they are the
/// residues the atom sites place on the entity, by `label_seq_id`. The
/// canonical sequence is `_entity_poly.pdbx_seq_one_letter_code_can` without
/// its line breaks, or else the monomers' codes joined.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolymerEntity {
    pub(crate) id: String,
    pub(crate) polymer_type: Option<PolymerType>,
    pub(crate) monomers: Vec<EntityMonomer>,
    pub(crate) monomers_source: SequenceSource,
    pub(crate) sequence: String,
    pub(crate) sequence_source: SequenceSource,
    pub(crate) instances: Vec<PolymerEntityInstance>,
}

impl PolymerEntity {
    /// The entity identifier, `_entity.id`.
    pub fn id(&self) -> &str {
        &self.id
    }

    /// The polymer type, `None` when the block states none.
    pub const fn polymer_type(&self) -> Option<&PolymerType> {
        self.polymer_type.as_ref()
    }

    /// The monomers in sequence order.
    pub fn monomers(&self) -> &[EntityMonomer] {
        &self.monomers
    }

    /// Whether `_entity_poly_seq` lists the monomers.
    pub const fn monomers_source(&self) -> SequenceSource {
        self.monomers_source
    }

    /// The canonical one-letter sequence, where each modified monomer has its
    /// parent's code, such as `MSKGEELF…`.
    pub fn sequence(&self) -> &str {
        &self.sequence
    }

    /// Whether the block states the canonical sequence.
    pub const fn sequence_source(&self) -> SequenceSource {
        self.sequence_source
    }

    /// The entity's instances (`label_asym_id`), in the block's order.
    pub fn instances(&self) -> &[PolymerEntityInstance] {
        &self.instances
    }
}

/// One monomer of a polymer entity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EntityMonomer {
    pub(crate) number: i32,
    pub(crate) component_id: String,
    pub(crate) alternatives: Vec<String>,
    pub(crate) code: ResidueCode,
}

impl EntityMonomer {
    /// Its position in the entity, `_entity_poly_seq.num`, which an atom
    /// site's `label_seq_id` names.
    pub const fn number(&self) -> i32 {
        self.number
    }

    /// The component identifier, such as `MSE`.
    pub fn component_id(&self) -> &str {
        &self.component_id
    }

    /// Further components the block lists at the same position, where the
    /// geometry holds either (microheterogeneity).
    pub fn alternatives(&self) -> &[String] {
        &self.alternatives
    }

    /// Its one-letter code. Where the canonical sequence the block states
    /// gives the monomer other letters than its component, the block's
    /// letters apply, with [`CodeOrigin::Sequence`].
    pub const fn code(&self) -> &ResidueCode {
        &self.code
    }
}

/// One instance of a polymer entity, a `label_asym_id`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolymerEntityInstance {
    pub(crate) asym_id: String,
    pub(crate) author_chain_id: Option<String>,
    pub(crate) chain: Option<ChainId>,
    pub(crate) positions: Vec<SequencePosition>,
}

impl PolymerEntityInstance {
    /// The label chain identifier, `label_asym_id`.
    pub fn asym_id(&self) -> &str {
        &self.asym_id
    }

    /// The author chain identifier, from `_pdbx_poly_seq_scheme` or the atom
    /// sites.
    pub fn author_chain_id(&self) -> Option<&str> {
        self.author_chain_id.as_deref()
    }

    /// The structure's chain of this instance; `None` when the structure has no
    /// atoms of it.
    pub const fn chain(&self) -> Option<ChainId> {
        self.chain
    }

    /// One position per [entity monomer](PolymerEntity::monomers), in the
    /// same order.
    pub fn positions(&self) -> &[SequencePosition] {
        &self.positions
    }
}

/// One monomer of one entity instance: the structure's residue there and its
/// author numbering.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SequencePosition {
    pub(crate) residue: Option<ResidueId>,
    pub(crate) author_seq_id: Option<String>,
    pub(crate) insertion_code: Option<String>,
}

impl SequencePosition {
    /// The structure's residue at this position; `None` where the structure
    /// does not observe the monomer.
    pub const fn residue(&self) -> Option<ResidueId> {
        self.residue
    }

    /// The author residue number, from `_pdbx_poly_seq_scheme.pdb_seq_num`
    /// for observed and unobserved monomers alike, or from the observed
    /// residue.
    pub fn author_seq_id(&self) -> Option<&str> {
        self.author_seq_id.as_deref()
    }

    /// The author insertion code, from the same source.
    pub fn insertion_code(&self) -> Option<&str> {
        self.insertion_code.as_deref()
    }
}

/// The polymer entities a block states, with every instance's positions but
/// without the structure's chains and residues, which
/// [`place_polymer_entities`] adds.
pub(super) fn read_polymer_entities(
    block: &MmcifBlock,
    rows: &[AtomRow],
    library: Option<&ComponentLibrary>,
) -> Result<Vec<PolymerEntity>, FormatError> {
    let mut entities = Entities::default();
    if let Some(table) = block.category("_entity")? {
        for row in 0..table.row_count() {
            if required(&table, row, "_entity.type")?.eq_ignore_ascii_case("polymer") {
                entities.add(required(&table, row, "_entity.id")?);
            }
        }
    }
    let mut stated_sequences = BTreeMap::new();
    if let Some(table) = block.category("_entity_poly")? {
        for row in 0..table.row_count() {
            let id = required(&table, row, "_entity_poly.entity_id")?;
            let index = entities.add(id);
            if stated_sequences.contains_key(&index) {
                return Err(row_error(
                    &table,
                    row,
                    DiagnosticKind::Inconsistent,
                    format!("Polymer entity {id:?} appears twice."),
                ));
            }
            entities.list[index].polymer_type =
                optional(&table, row, "_entity_poly.type").map(PolymerType::from_mmcif);
            let sequence = optional(&table, row, "_entity_poly.pdbx_seq_one_letter_code_can")
                .map(|sequence| sequence.split_whitespace().collect::<String>());
            stated_sequences.insert(index, sequence);
        }
    }
    if let Some(table) = block.category("_entity_poly_seq")? {
        let entity = Field::resolve(&table, "_entity_poly_seq.entity_id");
        let number = Field::resolve(&table, "_entity_poly_seq.num");
        let component = Field::resolve(&table, "_entity_poly_seq.mon_id");
        for row in 0..table.row_count() {
            let index = entities.add(required(&table, row, entity)?);
            let number = required_i32(&table, row, number)?;
            let entity = &mut entities.list[index];
            add_monomer(entity, number, required(&table, row, component)?);
            entity.monomers_source = SequenceSource::Stated;
        }
    }
    let Entities {
        list: mut entities,
        indices,
    } = entities;
    if entities.is_empty() {
        return Ok(entities);
    }

    let derived = entities
        .iter()
        .filter(|entity| entity.monomers_source == SequenceSource::Derived)
        .map(|entity| entity.id.as_str())
        .collect::<HashSet<_>>();

    // Instances in `_struct_asym` order, then any others the scheme or the
    // atom sites name, and monomers of entities `_entity_poly_seq` omits.
    let mut instances = Vec::<(usize, String)>::new();
    let mut listed = HashSet::<String>::new();
    if let Some(table) = block.category("_struct_asym")? {
        for row in 0..table.row_count() {
            let entity = required(&table, row, "_struct_asym.entity_id")?;
            if let Some(&index) = indices.get(entity) {
                let asym = required(&table, row, "_struct_asym.id")?;
                if listed.insert(asym.to_owned()) {
                    instances.push((index, asym.to_owned()));
                }
            }
        }
    }
    let mut schemes = HashMap::<String, Scheme>::new();
    if let Some(table) = block.category("_pdbx_poly_seq_scheme")? {
        read_schemes(&table, &indices, &mut instances, &mut listed, &mut schemes)?;
    }
    let mut derived_monomers = BTreeMap::<String, BTreeMap<i32, Vec<String>>>::new();
    let mut previous = None::<(&str, Option<i32>, &str)>;
    for row in rows {
        // Rows run residue by residue, so most repeat the previous residue.
        if previous.is_some_and(|(asym, number, component)| {
            asym == row.asym_id
                && (derived.is_empty() || (number == row.label_seq_id && component == row.comp_id))
        }) {
            continue;
        }
        previous = Some((&row.asym_id, row.label_seq_id, &row.comp_id));
        let Some(entity) = row.entity_id.as_deref() else {
            continue;
        };
        let Some(&index) = indices.get(entity) else {
            continue;
        };
        if !listed.contains(&row.asym_id) {
            listed.insert(row.asym_id.clone());
            instances.push((index, row.asym_id.clone()));
        }
        if let (true, Some(number)) = (derived.contains(entity), row.label_seq_id) {
            let components = derived_monomers
                .entry(entity.to_owned())
                .or_default()
                .entry(number)
                .or_default();
            if !components.contains(&row.comp_id) {
                components.push(row.comp_id.clone());
            }
        }
    }
    for (entity, monomers) in derived_monomers {
        let entity = &mut entities[indices[&entity]];
        for (number, components) in monomers {
            for component in components {
                add_monomer(entity, number, &component);
            }
        }
    }

    let stated = read_stated_components(block)?;
    let sources = CodeSources::new(Some(&stated), library);
    let mut codes = HashMap::<(String, &str), ResidueCode>::new();
    for (index, entity) in entities.iter_mut().enumerate() {
        sort_monomers(&mut entity.monomers);
        let unknown = match &entity.polymer_type {
            Some(kind) if kind.is_nucleic_acid() => "N",
            Some(_) => "X",
            None if entity.monomers.iter().any(|monomer| {
                sources.linking(&monomer.component_id) == Some(ComponentLinking::Nucleotide)
            }) =>
            {
                "N"
            }
            None => "X",
        };
        for monomer in &mut entity.monomers {
            monomer.code = codes
                .entry((monomer.component_id.clone(), unknown))
                .or_insert_with(|| {
                    sources
                        .code(&monomer.component_id, Some(unknown))
                        .expect("an unknown is supplied")
                })
                .clone();
        }
        match stated_sequences.remove(&index).flatten() {
            Some(sequence) => {
                adopt_stated_letters(&mut entity.monomers, &sequence);
                entity.sequence = sequence;
                entity.sequence_source = SequenceSource::Stated;
            }
            None => {
                entity.sequence = entity
                    .monomers
                    .iter()
                    .map(|monomer| monomer.code.as_str())
                    .collect();
            }
        }
    }
    for (index, asym) in instances {
        let entity = &mut entities[index];
        let scheme = schemes.remove(&asym).unwrap_or_default();
        let mut positions = vec![SequencePosition::default(); entity.monomers.len()];
        for (number, author_seq_id, insertion_code) in scheme.positions {
            if let Some(position) = position_of(&entity.monomers, number) {
                let position = &mut positions[position];
                if position.author_seq_id.is_none() {
                    position.author_seq_id = Some(author_seq_id);
                    position.insertion_code = insertion_code;
                }
            }
        }
        entity.instances.push(PolymerEntityInstance {
            asym_id: asym,
            author_chain_id: scheme.author_chain_id,
            chain: None,
            positions,
        });
    }
    entities.retain(|entity| !entity.monomers.is_empty());
    Ok(entities)
}

/// The positions one asym's `_pdbx_poly_seq_scheme` rows number.
#[derive(Default)]
struct Scheme {
    author_chain_id: Option<String>,
    /// Sequence number, author number, and insertion code.
    positions: Vec<(i32, String, Option<String>)>,
}

fn read_schemes(
    table: &MmcifCategory<'_>,
    indices: &HashMap<String, usize>,
    instances: &mut Vec<(usize, String)>,
    listed: &mut HashSet<String>,
    schemes: &mut HashMap<String, Scheme>,
) -> Result<(), FormatError> {
    let asym = Field::resolve(table, "_pdbx_poly_seq_scheme.asym_id");
    let entity = Field::resolve(table, "_pdbx_poly_seq_scheme.entity_id");
    let number = Field::resolve(table, "_pdbx_poly_seq_scheme.seq_id");
    let pdb_number = Field::resolve(table, "_pdbx_poly_seq_scheme.pdb_seq_num");
    let author_number = Field::resolve(table, "_pdbx_poly_seq_scheme.auth_seq_num");
    let insertion = Field::resolve(table, "_pdbx_poly_seq_scheme.pdb_ins_code");
    let strand = Field::resolve(table, "_pdbx_poly_seq_scheme.pdb_strand_id");
    for row in 0..table.row_count() {
        let asym_id = required(table, row, asym)?;
        let entity = optional(table, row, entity).and_then(|entity| indices.get(entity));
        if let Some(&index) = entity {
            if listed.insert(asym_id.to_owned()) {
                instances.push((index, asym_id.to_owned()));
            }
        }
        let scheme = schemes.entry(asym_id.to_owned()).or_default();
        if scheme.author_chain_id.is_none() {
            scheme.author_chain_id = optional(table, row, strand).map(str::to_owned);
        }
        let author =
            optional(table, row, pdb_number).or_else(|| optional(table, row, author_number));
        if let (Some(number), Some(author)) = (optional_i32(table, row, number)?, author) {
            scheme.positions.push((
                number,
                author.to_owned(),
                optional(table, row, insertion).map(str::to_owned),
            ));
        }
    }
    Ok(())
}

/// The block's `_chem_comp` rows, keyed by upper-case identifier.
fn read_stated_components(
    block: &MmcifBlock,
) -> Result<BTreeMap<String, StatedComponent>, FormatError> {
    let mut components = BTreeMap::new();
    let Some(table) = block.category("_chem_comp")? else {
        return Ok(components);
    };
    for row in 0..table.row_count() {
        let id = required(&table, row, "_chem_comp.id")?;
        let component = StatedComponent {
            linking: optional(&table, row, "_chem_comp.type").map(ComponentLinking::from_type),
            code: optional(&table, row, "_chem_comp.one_letter_code").map(str::to_owned),
            parents: optional(&table, row, "_chem_comp.mon_nstd_parent_comp_id")
                .map(parent_ids)
                .unwrap_or_default(),
        };
        if components
            .insert(id.to_ascii_uppercase(), component)
            .is_some()
        {
            return Err(row_error(
                &table,
                row,
                DiagnosticKind::Inconsistent,
                format!("Component {id:?} appears twice."),
            ));
        }
    }
    Ok(components)
}

/// Polymer entities in the order the block first names them.
#[derive(Default)]
struct Entities {
    list: Vec<PolymerEntity>,
    indices: HashMap<String, usize>,
}

impl Entities {
    fn add(&mut self, id: &str) -> usize {
        if let Some(&index) = self.indices.get(id) {
            return index;
        }
        self.list.push(PolymerEntity {
            id: id.to_owned(),
            polymer_type: None,
            monomers: Vec::new(),
            monomers_source: SequenceSource::Derived,
            sequence: String::new(),
            sequence_source: SequenceSource::Derived,
            instances: Vec::new(),
        });
        self.indices.insert(id.to_owned(), self.list.len() - 1);
        self.list.len() - 1
    }
}

/// Adds a monomer, or an alternative when the last one has its number.
fn add_monomer(entity: &mut PolymerEntity, number: i32, component: &str) {
    if let Some(monomer) = entity
        .monomers
        .last_mut()
        .filter(|monomer| monomer.number == number)
    {
        if monomer.component_id != component
            && !monomer.alternatives.iter().any(|other| other == component)
        {
            monomer.alternatives.push(component.to_owned());
        }
        return;
    }
    entity.monomers.push(EntityMonomer {
        number,
        component_id: component.to_owned(),
        alternatives: Vec::new(),
        code: ResidueCode::new("X", CodeOrigin::Unknown),
    });
}

/// Sorts monomers by number, merging those of one number into alternatives.
fn sort_monomers(monomers: &mut Vec<EntityMonomer>) {
    if monomers
        .windows(2)
        .all(|pair| pair[0].number < pair[1].number)
    {
        return;
    }
    monomers.sort_by_key(|monomer| monomer.number);
    let mut merged = Vec::<EntityMonomer>::with_capacity(monomers.len());
    for monomer in monomers.drain(..) {
        match merged.last_mut() {
            Some(last) if last.number == monomer.number => {
                for component in std::iter::once(monomer.component_id).chain(monomer.alternatives) {
                    if last.component_id != component && !last.alternatives.contains(&component) {
                        last.alternatives.push(component);
                    }
                }
            }
            _ => merged.push(monomer),
        }
    }
    *monomers = merged;
}

/// Gives each monomer its letters of the stated canonical sequence where they
/// can be told apart: when the monomers' codes are as long as the sequence
/// in all, or each monomer has one letter.
fn adopt_stated_letters(monomers: &mut [EntityMonomer], sequence: &str) {
    if !sequence.is_ascii() {
        return;
    }
    let total = monomers
        .iter()
        .map(|monomer| monomer.code.as_str().len())
        .sum::<usize>();
    let lengths = if total == sequence.len() {
        monomers
            .iter()
            .map(|monomer| monomer.code.as_str().len())
            .collect::<Vec<_>>()
    } else if monomers.len() == sequence.len() {
        vec![1; monomers.len()]
    } else {
        return;
    };
    let mut start = 0;
    for (monomer, length) in monomers.iter_mut().zip(lengths) {
        let letters = &sequence[start..start + length];
        start += length;
        if letters != monomer.code.as_str() {
            monomer.code = monomer.code.restated(letters, CodeOrigin::Sequence);
        }
    }
}

fn position_of(monomers: &[EntityMonomer], number: i32) -> Option<usize> {
    monomers
        .binary_search_by_key(&number, |monomer| monomer.number)
        .ok()
}

fn required_i32(table: &MmcifCategory<'_>, row: usize, field: Field) -> Result<i32, FormatError> {
    let text = required(table, row, field)?;
    text.parse().map_err(|_| {
        row_error(
            table,
            row,
            DiagnosticKind::InvalidValue,
            format!("The {} value {text:?} is not a whole number.", field.tag()),
        )
    })
}

/// Adds the structure's chains and residues to the entities' instances.
pub(super) fn place_polymer_entities(
    entities: &[PolymerEntity],
    topology: &Topology,
) -> Vec<PolymerEntity> {
    if entities.is_empty() {
        return Vec::new();
    }
    let chains = topology
        .chains()
        .map(|chain| (chain.label_id(), chain))
        .collect::<HashMap<_, _>>();
    let mut placed = entities.to_vec();
    for entity in &mut placed {
        for instance in &mut entity.instances {
            let Some(&chain) = chains.get(instance.asym_id.as_str()) else {
                continue;
            };
            instance.chain = Some(chain.id());
            if instance.author_chain_id.is_none() {
                instance.author_chain_id = chain.author_id().map(str::to_owned);
            }
            let mut by_author = None::<HashMap<(String, Option<String>), usize>>;
            for residue in chain.residues() {
                let position = match residue.label_seq_id() {
                    Some(number) => position_of(&entity.monomers, number),
                    None => residue.author_seq_id().and_then(|author| {
                        by_author
                            .get_or_insert_with(|| {
                                instance
                                    .positions
                                    .iter()
                                    .enumerate()
                                    .filter_map(|(index, position)| {
                                        Some((
                                            (
                                                position.author_seq_id.clone()?,
                                                position.insertion_code.clone(),
                                            ),
                                            index,
                                        ))
                                    })
                                    .collect()
                            })
                            .get(&(
                                author.to_owned(),
                                residue.insertion_code().map(str::to_owned),
                            ))
                            .copied()
                    }),
                };
                let Some(position) = position.map(|index| &mut instance.positions[index]) else {
                    continue;
                };
                if position.residue.is_some() {
                    continue;
                }
                position.residue = Some(residue.id());
                if position.author_seq_id.is_none() {
                    position.author_seq_id = residue.author_seq_id().map(str::to_owned);
                    position.insertion_code = residue.insertion_code().map(str::to_owned);
                }
            }
        }
    }
    placed
}
