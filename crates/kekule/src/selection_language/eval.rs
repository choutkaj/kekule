//! Evaluation of a parsed selection over dense atom bit sets.

use std::cmp::Ordering;
use std::ops::{ControlFlow, Range};
use std::sync::Arc;

use super::ast::{
    DistanceUnit, Field, Flag, HostTermKind, Level, Literal, Node, NodeKind, ValueMatch,
};
use super::{
    HostRequest, HostRole, HostTerm, HostValue, SelectEvaluateError, SelectEvaluateOptions,
    SelectionExpr,
};
use crate::geometry::Point3;
use crate::molecule::{Element, Molecule};
use crate::properties::{PropertyColumn, PropertyKey};
use crate::structure::ModelView;
use crate::topology::{
    AtomSelection, InstanceAtomId, MoleculeClass, ResidueClass, ResidueView, Topology,
    TopologyAtomIndex,
};
use crate::units::{self, Quantity, Unit, CANONICAL_LENGTH_UNIT};

type Result<T> = std::result::Result<T, SelectEvaluateError>;

/// Orders an entity's value against a literal; `None` when unordered.
type Compare<'c> = &'c dyn Fn(&Literal) -> Option<Ordering>;

/// A set of dense atom indices.
#[derive(Debug, Clone)]
struct Bits {
    words: Vec<u64>,
    len: usize,
}

impl Bits {
    fn empty(len: usize) -> Self {
        Self {
            words: vec![0; len.div_ceil(64)],
            len,
        }
    }

    fn full(len: usize) -> Self {
        let mut bits = Self {
            words: vec![!0; len.div_ceil(64)],
            len,
        };
        bits.trim();
        bits
    }

    fn trim(&mut self) {
        if !self.len.is_multiple_of(64) {
            if let Some(last) = self.words.last_mut() {
                *last &= (1u64 << (self.len % 64)) - 1;
            }
        }
    }

    fn set(&mut self, index: usize) {
        if let Some(word) = self.words.get_mut(index / 64) {
            *word |= 1 << (index % 64);
        }
    }

    fn get(&self, index: usize) -> bool {
        self.words
            .get(index / 64)
            .is_some_and(|word| word & (1 << (index % 64)) != 0)
    }

    fn and(&mut self, other: &Self) {
        for (word, other) in self.words.iter_mut().zip(&other.words) {
            *word &= other;
        }
    }

    fn or(&mut self, other: &Self) {
        for (word, other) in self.words.iter_mut().zip(&other.words) {
            *word |= other;
        }
    }

    fn not(&mut self) {
        for word in &mut self.words {
            *word = !*word;
        }
        self.trim();
    }

    fn ones(&self) -> impl Iterator<Item = usize> + '_ {
        self.words.iter().enumerate().flat_map(|(index, &word)| {
            let mut word = word;
            std::iter::from_fn(move || {
                (word != 0).then(|| {
                    let bit = word.trailing_zeros() as usize;
                    word &= word - 1;
                    index * 64 + bit
                })
            })
        })
    }
}

/// A field value of one atom, residue, chain, or molecule.
#[derive(Debug, Clone, Copy)]
enum Value<'a> {
    Integer(i64),
    Real(f64),
    Text(&'a str),
    Element(Element),
    ResId(i64, Option<&'a str>),
    ResidueClass(ResidueClass),
    MoleculeClass(MoleculeClass),
}

const PROTEIN_BACKBONE: [&str; 4] = ["N", "CA", "C", "O"];
const PROTEIN_MAIN_CHAIN: [&str; 5] = ["N", "CA", "C", "O", "OXT"];
const NUCLEIC_BACKBONE: [&str; 11] = [
    "P", "OP1", "OP2", "OP3", "O1P", "O2P", "O5'", "C5'", "C4'", "C3'", "O3'",
];

pub(crate) struct Evaluator<'a> {
    expression: &'a SelectionExpr,
    topology: &'a Arc<Topology>,
    structure: Option<ModelView<'a>>,
    options: &'a SelectEvaluateOptions<'a>,
    /// The atom site of each dense atom, computed on first use.
    sites: Option<Vec<Option<crate::topology::AtomSiteId>>>,
}

impl<'a> Evaluator<'a> {
    pub(crate) fn new(
        expression: &'a SelectionExpr,
        topology: &'a Arc<Topology>,
        structure: Option<ModelView<'a>>,
        options: &'a SelectEvaluateOptions<'a>,
    ) -> Self {
        Self {
            expression,
            topology,
            structure,
            options,
            sites: None,
        }
    }

    pub(crate) fn run(mut self) -> Result<AtomSelection> {
        let bits = self.node(&self.expression.root)?;
        let indices = bits
            .ones()
            .map(|index| TopologyAtomIndex::new(index as u32))
            .collect();
        Ok(AtomSelection::from_sorted_indices(
            Arc::clone(self.topology),
            indices,
        ))
    }

    fn len(&self) -> usize {
        self.topology.atom_count()
    }

    fn term(&self, node: &Node) -> (String, Range<usize>) {
        (super::NodeDisplay(node).to_string(), node.span.0.clone())
    }

    fn requires_structure(&self, node: &Node) -> Result<ModelView<'a>> {
        self.structure.ok_or_else(|| {
            let (term, span) = self.term(node);
            SelectEvaluateError::RequiresStructure { term, span }
        })
    }

    fn missing_perception(&self, node: &Node) -> SelectEvaluateError {
        let (term, span) = self.term(node);
        SelectEvaluateError::MissingPerception { term, span }
    }

    fn atom(&self, index: usize) -> (InstanceAtomId, &'a Molecule) {
        let topology: &'a Topology = self.topology;
        let id = topology.atom_ids()[index];
        (id, topology.definition_molecule(id.instance()))
    }

    fn element(&self, index: usize) -> Option<Element> {
        let (id, molecule) = self.atom(index);
        molecule.atom(id.atom()).map(|atom| atom.element())
    }

    fn dense(&self, atom: InstanceAtomId) -> Option<usize> {
        self.topology.atom_index(atom).map(|index| index.index())
    }

    fn sites(&mut self) -> &[Option<crate::topology::AtomSiteId>] {
        if self.sites.is_none() {
            let topology: &Topology = self.topology;
            let mut sites = vec![None; topology.atom_count()];
            for site in topology.atom_sites() {
                if let Some(index) = topology.atom_index(site.atom().id()) {
                    sites[index.index()] = Some(site.id());
                }
            }
            self.sites = Some(sites);
        }
        self.sites.as_deref().unwrap_or_default()
    }

    /// The atom name: the author atom ID, else the label atom ID.
    fn name(&mut self, index: usize, label: bool) -> Option<&'a str> {
        let site = (*self.sites().get(index)?)?;
        let topology: &'a Topology = self.topology;
        let metadata = topology.atom_site(site)?.metadata();
        if label {
            metadata.label_atom_id()
        } else {
            metadata.auth_atom_id().or(metadata.label_atom_id())
        }
    }

    fn residue_atoms(&self, residue: ResidueView<'_>, bits: &mut Bits) {
        for site in residue.atom_sites() {
            if let Some(index) = self.dense(site.atom().id()) {
                bits.set(index);
            }
        }
    }

    fn node(&mut self, node: &Node) -> Result<Bits> {
        let len = self.len();
        Ok(match &node.kind {
            NodeKind::All => Bits::full(len),
            NodeKind::None => Bits::empty(len),
            NodeKind::And(operands) => {
                let mut bits = Bits::full(len);
                for operand in operands {
                    let other = self.node(operand)?;
                    bits.and(&other);
                }
                bits
            }
            NodeKind::Or(operands) => {
                let mut bits = Bits::empty(len);
                for operand in operands {
                    let other = self.node(operand)?;
                    bits.or(&other);
                }
                bits
            }
            NodeKind::Not(inner) => {
                let mut bits = self.node(inner)?;
                bits.not();
                bits
            }
            NodeKind::Within {
                distance,
                unit,
                reference,
            } => self.within(node, *distance, *unit, reference)?,
            NodeKind::Expand { level, inner } => {
                let inner = self.node(inner)?;
                self.expand(*level, &inner)
            }
            NodeKind::Bonded(inner) => {
                let inner = self.node(inner)?;
                self.bonded(&inner)
            }
            NodeKind::Flag(flag) => self.flag(node, *flag)?,
            NodeKind::Match { field, values } => self.field(node, field, |value, literals| {
                values.iter().any(|candidate| match candidate {
                    ValueMatch::Exact(literal) => {
                        literals(literal).is_some_and(|ordering| ordering == Ordering::Equal)
                    }
                    ValueMatch::Range(low, high) => {
                        literals(low).is_some_and(|ordering| ordering != Ordering::Less)
                            && literals(high).is_some_and(|ordering| ordering != Ordering::Greater)
                    }
                }) && value
            })?,
            NodeKind::Compare { field, op, value } => {
                self.field(node, field, |present, literals| {
                    present && literals(value).is_some_and(|ordering| op.holds(ordering))
                })?
            }
            NodeKind::PropertyTrue(key) => {
                let literal = Literal::Property {
                    text: "true".into(),
                    number: None,
                    unit: None,
                };
                self.property(node, key, &|ordering| {
                    ordering(&literal) == Some(Ordering::Equal)
                })?
            }
            NodeKind::Smarts(pattern) => self.smarts(node, pattern)?,
            NodeKind::Host(kind, text) => match self.host(node, *kind, text, HostRole::Atoms)? {
                HostValue::Atoms(selection) => {
                    let mut bits = Bits::empty(len);
                    for index in selection.indices() {
                        bits.set(index.index());
                    }
                    bits
                }
                HostValue::Points(_) => {
                    let (term, span) = self.term(node);
                    return Err(SelectEvaluateError::HostPoints { term, span });
                }
            },
        })
    }

    fn host(
        &self,
        node: &Node,
        kind: HostTermKind,
        text: &str,
        role: HostRole,
    ) -> Result<HostValue> {
        let (term, span) = self.term(node);
        let Some(resolver) = self.options.resolver() else {
            return Err(SelectEvaluateError::NoResolver { term, span });
        };
        let request = HostRequest {
            topology: self.topology,
            structure: self.structure,
            role,
        };
        let value = resolver
            .resolve(&HostTerm { kind, text }, request)
            .map_err(|source| SelectEvaluateError::Host {
                term: term.clone(),
                span: span.clone(),
                source,
            })?;
        if let HostValue::Atoms(selection) = &value {
            if selection.ensure_compatible(self.topology).is_err() {
                return Err(SelectEvaluateError::HostTopologyMismatch { term, span });
            }
        }
        Ok(value)
    }

    fn within(
        &mut self,
        node: &Node,
        distance: f64,
        unit: DistanceUnit,
        reference: &Node,
    ) -> Result<Bits> {
        let view = self.requires_structure(node)?;
        let unit = match unit {
            DistanceUnit::Angstrom => units::ANGSTROM,
            DistanceUnit::Nanometer => units::NANOMETER,
            DistanceUnit::Picometer => Unit::new(units::Dimension::LENGTH, 1.0e-12, Some("pm"))
                .map_err(|source| {
                    let (term, span) = self.term(node);
                    SelectEvaluateError::Unit {
                        term,
                        span,
                        source: Box::new(source),
                    }
                })?,
            DistanceUnit::Bohr => units::BOHR,
        };
        let cutoff = Quantity::new(distance, unit)
            .value_in(CANONICAL_LENGTH_UNIT)
            .map_err(|source| {
                let (term, span) = self.term(node);
                SelectEvaluateError::Unit {
                    term,
                    span,
                    source: Box::new(source),
                }
            })?;
        let positions: &[Point3] = view.positions().values().value();
        let references: Vec<Point3> = match &reference.kind {
            NodeKind::Host(kind, text) => {
                match self.host(reference, *kind, text, HostRole::Reference)? {
                    HostValue::Atoms(selection) => selection
                        .indices()
                        .iter()
                        .filter_map(|index| positions.get(index.index()).copied())
                        .collect(),
                    HostValue::Points(points) => points.values().value().to_vec(),
                }
            }
            _ => {
                let bits = self.node(reference)?;
                bits.ones()
                    .filter_map(|index| positions.get(index).copied())
                    .collect()
            }
        };
        let mut bits = Bits::empty(self.len());
        super::grid::within(positions, &references, cutoff, |index| bits.set(index));
        Ok(bits)
    }

    fn expand(&mut self, level: Level, inner: &Bits) -> Bits {
        let topology: &'a Topology = self.topology;
        let mut bits = inner.clone();
        match level {
            Level::Molecule => {
                let mut instances = vec![false; topology.instance_count()];
                for index in inner.ones() {
                    let (id, _) = self.atom(index);
                    instances[id.instance().index()] = true;
                }
                for instance in topology.instances() {
                    if instances[instance.id().index()] {
                        for atom in instance.atoms() {
                            bits.set(atom.index().index());
                        }
                    }
                }
            }
            Level::Residue | Level::Chain => {
                let mut residues = vec![false; topology.residues().len()];
                let mut chains = vec![false; topology.chains().len()];
                for index in inner.ones() {
                    if let Some(site) = self.sites().get(index).copied().flatten() {
                        if let Some(site) = topology.atom_site(site) {
                            let residue = site.residue();
                            residues[residue.id().index()] = true;
                            chains[residue.chain().id().index()] = true;
                        }
                    }
                }
                for residue in topology.residues() {
                    let selected = match level {
                        Level::Residue => residues[residue.id().index()],
                        _ => chains[residue.chain().id().index()],
                    };
                    if selected {
                        self.residue_atoms(residue, &mut bits);
                    }
                }
            }
        }
        bits
    }

    fn bonded(&self, inner: &Bits) -> Bits {
        let mut bits = Bits::empty(self.len());
        for index in inner.ones() {
            let (id, molecule) = self.atom(index);
            if let Some(neighbors) = molecule.neighbors(id.atom()) {
                for neighbor in neighbors {
                    if let Some(dense) = self.dense(InstanceAtomId::new(id.instance(), neighbor)) {
                        bits.set(dense);
                    }
                }
            }
        }
        bits
    }

    fn molecules_of(&self, classes: &[MoleculeClass]) -> Bits {
        let mut bits = Bits::empty(self.len());
        for instance in self.topology.instances() {
            if classes.contains(&instance.class()) {
                for atom in instance.atoms() {
                    bits.set(atom.index().index());
                }
            }
        }
        bits
    }

    fn flag(&mut self, node: &Node, flag: Flag) -> Result<Bits> {
        let len = self.len();
        let mut bits = Bits::empty(len);
        match flag {
            Flag::Hydrogen | Flag::Heavy => {
                for index in 0..len {
                    let hydrogen = self.element(index) == Some(Element::H);
                    if hydrogen == (flag == Flag::Hydrogen) {
                        bits.set(index);
                    }
                }
            }
            Flag::PolarHydrogen => {
                for index in 0..len {
                    let (id, molecule) = self.atom(index);
                    if self.element(index) != Some(Element::H) {
                        continue;
                    }
                    let polar = molecule.neighbors(id.atom()).is_some_and(|mut neighbors| {
                        neighbors.any(|neighbor| {
                            molecule.atom(neighbor).is_some_and(|atom| {
                                matches!(atom.element(), Element::N | Element::O)
                            })
                        })
                    });
                    if polar {
                        bits.set(index);
                    }
                }
            }
            Flag::Aromatic | Flag::Ring => {
                for index in 0..len {
                    let (id, molecule) = self.atom(index);
                    let perception = molecule.perception();
                    let member = if flag == Flag::Aromatic {
                        perception
                            .aromaticity()
                            .map(|aromaticity| aromaticity.atom_is_aromatic(id.atom()))
                    } else {
                        perception
                            .rings()
                            .map(|rings| rings.membership().atom_in_ring(id.atom()))
                    };
                    match member {
                        Some(true) => bits.set(index),
                        Some(false) => {}
                        None => return Err(self.missing_perception(node)),
                    }
                }
            }
            Flag::Protein => bits = self.molecules_of(&[MoleculeClass::Protein]),
            Flag::Nucleic => bits = self.molecules_of(&[MoleculeClass::Dna, MoleculeClass::Rna]),
            Flag::Water => bits = self.molecules_of(&[MoleculeClass::Water]),
            Flag::Ion => bits = self.molecules_of(&[MoleculeClass::Ion]),
            Flag::Ligand => bits = self.molecules_of(&[MoleculeClass::SmallMolecule]),
            Flag::Backbone | Flag::Sidechain => {
                let protein = self.molecules_of(&[MoleculeClass::Protein]);
                let nucleic = self.molecules_of(&[MoleculeClass::Dna, MoleculeClass::Rna]);
                let topology: &'a Topology = self.topology;
                let mut main_chain = Bits::empty(len);
                let mut amino_acid = Bits::empty(len);
                for residue in topology.residues() {
                    let class = residue.class();
                    for site in residue.atom_sites() {
                        let Some(index) = self.dense(site.atom().id()) else {
                            continue;
                        };
                        let metadata = site.metadata();
                        let name = metadata.auth_atom_id().or(metadata.label_atom_id());
                        let name = name.unwrap_or_default();
                        if class == ResidueClass::AminoAcid && protein.get(index) {
                            amino_acid.set(index);
                            if PROTEIN_MAIN_CHAIN.contains(&name) {
                                main_chain.set(index);
                            }
                            if flag == Flag::Backbone && PROTEIN_BACKBONE.contains(&name) {
                                bits.set(index);
                            }
                        }
                        if flag == Flag::Backbone
                            && matches!(
                                class,
                                ResidueClass::DnaNucleotide | ResidueClass::RnaNucleotide
                            )
                            && nucleic.get(index)
                            && NUCLEIC_BACKBONE.contains(&name)
                        {
                            bits.set(index);
                        }
                    }
                }
                if flag == Flag::Sidechain {
                    for index in amino_acid.ones() {
                        if main_chain.get(index) {
                            continue;
                        }
                        let (id, molecule) = self.atom(index);
                        let on_main_chain = self.element(index) == Some(Element::H)
                            && molecule.neighbors(id.atom()).is_some_and(|mut neighbors| {
                                neighbors.any(|neighbor| {
                                    self.dense(InstanceAtomId::new(id.instance(), neighbor))
                                        .is_some_and(|dense| main_chain.get(dense))
                                })
                            });
                        if !on_main_chain {
                            bits.set(index);
                        }
                    }
                }
            }
        }
        Ok(bits)
    }

    fn smarts(&self, node: &Node, pattern: &str) -> Result<Bits> {
        let error = |source| {
            let (term, span) = self.term(node);
            SelectEvaluateError::Smarts {
                term,
                span,
                source: Box::new(source),
            }
        };
        let query = crate::smarts::parse_str(pattern, Default::default()).map_err(|source| {
            let (term, span) = self.term(node);
            SelectEvaluateError::InvalidSmarts { term, span, source }
        })?;
        let target = crate::substructure::PreparedTopologyTarget::new(self.topology);
        let mut bits = Bits::empty(self.len());
        target
            .visit(&query, self.options.smarts(), |matched| {
                for &atom in matched.atoms() {
                    if let Some(index) = self.dense(atom) {
                        bits.set(index);
                    }
                }
                ControlFlow::Continue(())
            })
            .map_err(error)?;
        Ok(bits)
    }

    /// Evaluates a field term: `test(present, compare)` decides one entity,
    /// where `compare` orders the entity's value against a literal.
    fn field(
        &mut self,
        node: &Node,
        field: &Field,
        test: impl Fn(bool, Compare<'_>) -> bool,
    ) -> Result<Bits> {
        let len = self.len();
        let mut bits = Bits::empty(len);
        let topology: &'a Topology = self.topology;
        match field {
            Field::Property(key) => {
                return self.property(node, key, &|compare| test(true, compare));
            }
            Field::Index => {
                for index in 0..len {
                    if decide(&test, Some(Value::Integer(index as i64))) {
                        bits.set(index);
                    }
                }
            }
            Field::Element | Field::Charge | Field::Isotope => {
                for index in 0..len {
                    let (id, molecule) = self.atom(index);
                    let Some(atom) = molecule.atom(id.atom()) else {
                        continue;
                    };
                    let value = match field {
                        Field::Element => Some(Value::Element(atom.element())),
                        Field::Charge => Some(Value::Integer(i64::from(atom.formal_charge()))),
                        _ => atom
                            .isotope()
                            .map(|isotope| Value::Integer(i64::from(isotope.get()))),
                    };
                    if decide(&test, value) {
                        bits.set(index);
                    }
                }
            }
            Field::Name | Field::LabelName => {
                let label = matches!(field, Field::LabelName);
                for index in 0..len {
                    let value = self.name(index, label).map(Value::Text);
                    if decide(&test, value) {
                        bits.set(index);
                    }
                }
            }
            Field::BFactor | Field::Occupancy => {
                let view = self.requires_structure(node)?;
                let values = if matches!(field, Field::BFactor) {
                    view.b_factors().map(|values| *values.value())
                } else {
                    view.occupancies()
                };
                // Literals of B factors are in square ångströms; stored values
                // are canonical.
                let factor = if matches!(field, Field::BFactor) {
                    units::SQUARE_ANGSTROM
                        .conversion_factor_to(units::CANONICAL_AREA_UNIT)
                        .map_err(|source| {
                            let (term, span) = self.term(node);
                            SelectEvaluateError::Unit {
                                term,
                                span,
                                source: Box::new(source),
                            }
                        })?
                } else {
                    1.0
                };
                if let Some(values) = values {
                    for (index, value) in values.iter().enumerate() {
                        let value = value.map(Value::Real);
                        if decide_scaled(&test, value, factor) {
                            bits.set(index);
                        }
                    }
                }
            }
            Field::ResName
            | Field::ResId
            | Field::LabelResId
            | Field::ResClass
            | Field::Chain
            | Field::LabelChain => {
                for residue in topology.residues() {
                    let value = match field {
                        Field::ResName => Some(Value::Text(residue.name())),
                        Field::ResId => match residue.author_seq_id() {
                            Some(text) => text
                                .trim()
                                .parse::<i64>()
                                .ok()
                                .map(|number| Value::ResId(number, residue.insertion_code())),
                            None => residue
                                .label_seq_id()
                                .map(|number| Value::ResId(i64::from(number), None)),
                        },
                        Field::LabelResId => residue
                            .label_seq_id()
                            .map(|number| Value::Integer(i64::from(number))),
                        Field::ResClass => Some(Value::ResidueClass(residue.class())),
                        Field::Chain => {
                            let chain = residue.chain();
                            Some(Value::Text(chain.author_id().unwrap_or(chain.label_id())))
                        }
                        _ => Some(Value::Text(residue.chain().label_id())),
                    };
                    if decide(&test, value) {
                        self.residue_atoms(residue, &mut bits);
                    }
                }
            }
            Field::MolClass | Field::Instance | Field::Definition => {
                for instance in topology.instances() {
                    let value = match field {
                        Field::MolClass => Value::MoleculeClass(instance.class()),
                        Field::Instance => Value::Integer(i64::from(instance.id().raw())),
                        _ => Value::Integer(i64::from(instance.definition_id().raw())),
                    };
                    if decide(&test, Some(value)) {
                        for atom in instance.atoms() {
                            bits.set(atom.index().index());
                        }
                    }
                }
            }
        }
        Ok(bits)
    }

    /// Evaluates a property term against the one scope holding `key`.
    fn property(
        &self,
        node: &Node,
        key: &PropertyKey,
        test: &dyn Fn(Compare<'_>) -> bool,
    ) -> Result<Bits> {
        let len = self.len();
        let topology: &'a Topology = self.topology;
        let structure_column = self
            .structure
            .and_then(|view| view.properties().atoms().get(key));
        let topology_column = topology.properties().atoms().get(key);
        let molecule_scope = topology.definitions().any(|definition| {
            definition
                .molecule()
                .properties()
                .atoms()
                .get(key)
                .is_some()
        });
        let scopes = usize::from(structure_column.is_some())
            + usize::from(topology_column.is_some())
            + usize::from(molecule_scope);
        let (term, span) = self.term(node);
        match scopes {
            0 => {
                return Err(SelectEvaluateError::UnknownProperty {
                    key: key.clone(),
                    span,
                })
            }
            1 => {}
            _ => {
                return Err(SelectEvaluateError::AmbiguousProperty {
                    key: key.clone(),
                    span,
                })
            }
        }
        let mut bits = Bits::empty(len);
        let mismatch = std::cell::Cell::new(false);
        for index in 0..len {
            let cell = if let Some(column) = structure_column.or(topology_column) {
                cell(column, index)
            } else {
                let (id, molecule) = self.atom(index);
                molecule
                    .properties()
                    .atoms()
                    .get(key)
                    .and_then(|column| cell(column, id.atom().index()))
            };
            let Some(cell) = cell else {
                continue;
            };
            let compare = |literal: &Literal| -> Option<Ordering> {
                let ordering = compare_property(cell, literal);
                if ordering.is_none() {
                    mismatch.set(true);
                }
                ordering.flatten()
            };
            if test(&compare) {
                bits.set(index);
            }
        }
        if mismatch.get() {
            return Err(SelectEvaluateError::PropertyType {
                key: key.clone(),
                term,
                span,
            });
        }
        Ok(bits)
    }
}

/// One property cell: its value and, for a real column, its unit.
#[derive(Debug, Clone, Copy)]
enum Cell<'a> {
    Bool(bool),
    Int(i64),
    Real(f64, Unit),
    Text(&'a str),
}

fn cell(column: &PropertyColumn, index: usize) -> Option<Cell<'_>> {
    match column {
        PropertyColumn::Bool(values) => values.get(index).copied().flatten().map(Cell::Bool),
        PropertyColumn::Int(values) => values.get(index).copied().flatten().map(Cell::Int),
        PropertyColumn::Real { unit, values } => values
            .get(index)
            .copied()
            .flatten()
            .map(|value| Cell::Real(value, *unit)),
        PropertyColumn::String(values) => values
            .get(index)
            .and_then(|value| value.as_deref())
            .map(Cell::Text),
    }
}

/// Orders a property cell against a literal: `None` when the literal cannot
/// apply to the column's type, `Some(None)` when the values are unordered.
fn compare_property(cell: Cell<'_>, literal: &Literal) -> Option<Option<Ordering>> {
    let Literal::Property { text, number, unit } = literal else {
        return None;
    };
    match cell {
        Cell::Bool(value) => {
            if unit.is_some() {
                return None;
            }
            let literal = match text.as_str() {
                "true" => true,
                "false" => false,
                _ => return None,
            };
            Some(Some(value.cmp(&literal)))
        }
        Cell::Int(value) => {
            let number = (*number)?;
            if unit.is_some() {
                return None;
            }
            Some((value as f64).partial_cmp(&number))
        }
        Cell::Real(value, column_unit) => {
            let number = (*number)?;
            let number = match unit {
                Some((unit, _)) => Quantity::new(number, *unit).value_in(column_unit).ok()?,
                None => number,
            };
            Some(value.partial_cmp(&number))
        }
        Cell::Text(value) => {
            if unit.is_some() {
                return None;
            }
            Some(Some(value.cmp(text.as_str())))
        }
    }
}

fn decide(test: &impl Fn(bool, Compare<'_>) -> bool, value: Option<Value<'_>>) -> bool {
    decide_scaled(test, value, 1.0)
}

/// Decides one entity whose real literals are scaled by `factor` into the
/// stored unit first.
fn decide_scaled(
    test: &impl Fn(bool, Compare<'_>) -> bool,
    value: Option<Value<'_>>,
    factor: f64,
) -> bool {
    let Some(value) = value else {
        return false;
    };
    test(true, &|literal| compare(value, literal, factor))
}

fn compare(value: Value<'_>, literal: &Literal, factor: f64) -> Option<Ordering> {
    match (value, literal) {
        (Value::Integer(value), Literal::Integer(literal)) => Some(value.cmp(literal)),
        (Value::Real(value), Literal::Real { value: literal, .. }) => {
            value.partial_cmp(&(literal * factor))
        }
        (Value::Text(value), Literal::Text(literal)) => Some(value.cmp(literal.as_str())),
        (Value::Element(value), Literal::Element(literal)) => Some(value.cmp(literal)),
        (
            Value::ResId(number, code),
            Literal::ResId {
                number: literal,
                insertion_code,
            },
        ) => Some(
            number
                .cmp(literal)
                .then_with(|| code.cmp(&insertion_code.as_deref())),
        ),
        (Value::ResidueClass(value), Literal::ResidueClass(literal)) => {
            Some(if value == *literal {
                Ordering::Equal
            } else {
                Ordering::Less
            })
        }
        (Value::MoleculeClass(value), Literal::MoleculeClass(literal)) => {
            Some(if value == *literal {
                Ordering::Equal
            } else {
                Ordering::Less
            })
        }
        _ => None,
    }
}
