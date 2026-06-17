#![forbid(unsafe_code)]

use std::collections::BTreeMap;
use std::fmt;

pub mod prelude {
    pub use crate::{
        perceive_aromaticity, perceive_ring_membership, read_sdf_v2000_str, AromaticityError,
        AromaticityModel, Atom, AtomId, AtomStereo, BioHierarchy, Bond, BondId, BondOrder,
        BondStereo, ComputedState, Element, MacroMolecule, Molecule, MoleculeError, PropMap,
        PropValue, Result, RingMembership, SdfParseError, SdfParseOptions, SdfRecord,
        SmallMolecule,
    };
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomId(u32);

impl AtomId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl fmt::Display for AtomId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "a{}", self.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondId(u32);

impl BondId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl fmt::Display for BondId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "b{}", self.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Element {
    atomic_number: u8,
}

const ELEMENT_SYMBOLS: [&str; 119] = [
    "?", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

impl Element {
    pub fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        if (1..=118).contains(&atomic_number) {
            Some(Self { atomic_number })
        } else {
            None
        }
    }

    pub fn from_symbol(symbol: &str) -> Option<Self> {
        let atomic_number = ELEMENT_SYMBOLS
            .iter()
            .position(|candidate| *candidate == symbol)?;
        if atomic_number == 0 {
            return None;
        }
        Some(Self {
            atomic_number: atomic_number as u8,
        })
    }

    pub const fn atomic_number(self) -> u8 {
        self.atomic_number
    }

    pub fn symbol(self) -> &'static str {
        ELEMENT_SYMBOLS
            .get(self.atomic_number as usize)
            .copied()
            .unwrap_or("?")
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

pub type PropMap = BTreeMap<String, PropValue>;

#[derive(Debug, Clone, PartialEq)]
pub enum PropValue {
    String(String),
    Int(i64),
    Float(f64),
    Bool(bool),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub element: Element,
    pub isotope: Option<u16>,
    pub formal_charge: i8,
    pub radical_electrons: u8,
    pub explicit_hydrogens: u8,
    pub implicit_hydrogens: Option<u8>,
    pub aromatic: bool,
    pub chiral: Option<AtomStereo>,
    pub atom_map: Option<u32>,
    pub props: PropMap,
}

impl Atom {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            isotope: None,
            formal_charge: 0,
            radical_electrons: 0,
            explicit_hydrogens: 0,
            implicit_hydrogens: None,
            aromatic: false,
            chiral: None,
            atom_map: None,
            props: PropMap::new(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomStereo {
    TetrahedralClockwise,
    TetrahedralCounterClockwise,
    Unspecified,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    a: AtomId,
    b: AtomId,
    pub order: BondOrder,
    pub aromatic: bool,
    pub stereo: Option<BondStereo>,
    pub props: PropMap,
}

impl Bond {
    pub fn new(a: AtomId, b: AtomId, order: BondOrder) -> Self {
        Self {
            a,
            b,
            order,
            aromatic: matches!(order, BondOrder::Aromatic),
            stereo: None,
            props: PropMap::new(),
        }
    }

    pub const fn a(&self) -> AtomId {
        self.a
    }

    pub const fn b(&self) -> AtomId {
        self.b
    }

    pub const fn endpoints(&self) -> (AtomId, AtomId) {
        (self.a, self.b)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Zero,
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    Dative,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondStereo {
    E,
    Z,
    Up,
    Down,
    Unspecified,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub enum ComputedState {
    #[default]
    Absent,
    Stale,
    Fresh,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PerceptionState {
    pub valence: ComputedState,
    pub rings: ComputedState,
    pub aromaticity: ComputedState,
    pub stereo: ComputedState,
}

impl PerceptionState {
    pub fn invalidate_all(&mut self) {
        self.valence = invalidate(self.valence);
        self.rings = invalidate(self.rings);
        self.aromaticity = invalidate(self.aromaticity);
        self.stereo = invalidate(self.stereo);
    }
}

fn invalidate(state: ComputedState) -> ComputedState {
    match state {
        ComputedState::Fresh => ComputedState::Stale,
        ComputedState::Stale | ComputedState::Absent => state,
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Molecule {
    atoms: Vec<Option<Atom>>,
    bonds: Vec<Option<Bond>>,
    adjacency: Vec<Vec<BondId>>,
    props: PropMap,
    perception: PerceptionState,
    ring_membership: Option<RingMembership>,
}

impl Molecule {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.iter().flatten().count()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.iter().flatten().count()
    }

    pub fn add_atom(&mut self, atom: Atom) -> AtomId {
        let id = AtomId::new(self.atoms.len() as u32);
        self.atoms.push(Some(atom));
        self.adjacency.push(Vec::new());
        self.invalidate_topology();
        id
    }

    pub fn delete_atom(&mut self, id: AtomId) -> Result<Atom> {
        self.atom(id)?;
        let incident = self.adjacency[id.index()].clone();
        for bond_id in incident {
            if self
                .bonds
                .get(bond_id.index())
                .and_then(Option::as_ref)
                .is_some()
            {
                self.delete_bond(bond_id)?;
            }
        }
        self.adjacency[id.index()].clear();
        let atom = self.atoms[id.index()]
            .take()
            .ok_or(MoleculeError::InvalidAtomId(id))?;
        self.invalidate_topology();
        Ok(atom)
    }

    pub fn atom(&self, id: AtomId) -> Result<&Atom> {
        self.atoms
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidAtomId(id))
    }

    pub fn atom_mut(&mut self, id: AtomId) -> Result<&mut Atom> {
        self.atom(id)?;
        self.invalidate_topology();
        self.atoms[id.index()]
            .as_mut()
            .ok_or(MoleculeError::InvalidAtomId(id))
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        self.atoms
            .iter()
            .enumerate()
            .filter_map(|(index, atom)| atom.as_ref().map(|atom| (AtomId::new(index as u32), atom)))
    }

    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms().map(|(id, _)| id)
    }

    pub fn add_bond(&mut self, a: AtomId, b: AtomId, order: BondOrder) -> Result<BondId> {
        self.atom(a)?;
        self.atom(b)?;
        if a == b {
            return Err(MoleculeError::SelfBond(a));
        }
        if self.bond_between(a, b)?.is_some() {
            return Err(MoleculeError::DuplicateBond { a, b });
        }
        let id = BondId::new(self.bonds.len() as u32);
        self.bonds.push(Some(Bond::new(a, b, order)));
        self.adjacency[a.index()].push(id);
        self.adjacency[b.index()].push(id);
        self.invalidate_topology();
        Ok(id)
    }

    pub fn delete_bond(&mut self, id: BondId) -> Result<Bond> {
        let bond = self
            .bonds
            .get_mut(id.index())
            .and_then(Option::take)
            .ok_or(MoleculeError::InvalidBondId(id))?;
        self.remove_incident_bond(bond.a, id);
        self.remove_incident_bond(bond.b, id);
        self.invalidate_topology();
        Ok(bond)
    }

    pub fn bond(&self, id: BondId) -> Result<&Bond> {
        self.bonds
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidBondId(id))
    }

    pub fn bond_mut(&mut self, id: BondId) -> Result<&mut Bond> {
        self.bond(id)?;
        self.invalidate_topology();
        self.bonds[id.index()]
            .as_mut()
            .ok_or(MoleculeError::InvalidBondId(id))
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &Bond)> {
        self.bonds
            .iter()
            .enumerate()
            .filter_map(|(index, bond)| bond.as_ref().map(|bond| (BondId::new(index as u32), bond)))
    }

    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bonds().map(|(id, _)| id)
    }

    pub fn neighbors(&self, id: AtomId) -> Result<impl Iterator<Item = AtomId> + '_> {
        self.atom(id)?;
        Ok(self.adjacency[id.index()]
            .iter()
            .filter_map(|bond_id| self.bond(*bond_id).ok())
            .map(move |bond| bond.other_atom(id)))
    }

    pub fn incident_bonds(&self, id: AtomId) -> Result<impl Iterator<Item = (BondId, &Bond)> + '_> {
        self.atom(id)?;
        Ok(self.adjacency[id.index()]
            .iter()
            .filter_map(|bond_id| self.bond(*bond_id).ok().map(|bond| (*bond_id, bond))))
    }

    pub fn bond_between(&self, a: AtomId, b: AtomId) -> Result<Option<BondId>> {
        self.atom(a)?;
        self.atom(b)?;
        Ok(self.adjacency[a.index()].iter().copied().find(|bond_id| {
            self.bond(*bond_id)
                .map(|bond| bond.connects(a, b))
                .unwrap_or(false)
        }))
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }

    pub fn perception(&self) -> &PerceptionState {
        &self.perception
    }

    pub fn perception_mut(&mut self) -> &mut PerceptionState {
        &mut self.perception
    }

    pub fn ring_membership(&self) -> Option<&RingMembership> {
        self.ring_membership.as_ref()
    }

    pub fn invalidate_topology(&mut self) {
        self.perception.invalidate_all();
    }

    fn remove_incident_bond(&mut self, atom: AtomId, bond: BondId) {
        if let Some(incident) = self.adjacency.get_mut(atom.index()) {
            incident.retain(|id| *id != bond);
        }
    }
}

impl Bond {
    fn connects(&self, a: AtomId, b: AtomId) -> bool {
        (self.a == a && self.b == b) || (self.a == b && self.b == a)
    }

    fn other_atom(&self, atom: AtomId) -> AtomId {
        if self.a == atom {
            self.b
        } else {
            self.a
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingMembership {
    atom_flags: Vec<bool>,
    bond_flags: Vec<bool>,
}

impl RingMembership {
    pub fn atom_in_ring(&self, atom: AtomId) -> bool {
        self.atom_flags.get(atom.index()).copied().unwrap_or(false)
    }

    pub fn bond_in_ring(&self, bond: BondId) -> bool {
        self.bond_flags.get(bond.index()).copied().unwrap_or(false)
    }

    pub fn ring_atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atom_flags
            .iter()
            .enumerate()
            .filter_map(|(index, in_ring)| in_ring.then_some(AtomId::new(index as u32)))
    }

    pub fn ring_bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bond_flags
            .iter()
            .enumerate()
            .filter_map(|(index, in_ring)| in_ring.then_some(BondId::new(index as u32)))
    }
}

pub fn perceive_ring_membership(mol: &mut Molecule) -> RingMembership {
    let membership = compute_ring_membership(mol);
    mol.ring_membership = Some(membership.clone());
    mol.perception.rings = ComputedState::Fresh;
    membership
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityModel {
    RdkitLikeBasic,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AromaticityError {
    UnsupportedElement(AtomId),
}

impl fmt::Display for AromaticityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedElement(id) => {
                write!(f, "unsupported aromaticity element at atom {id}")
            }
        }
    }
}

impl std::error::Error for AromaticityError {}

pub fn perceive_aromaticity(
    mol: &mut Molecule,
    model: AromaticityModel,
) -> std::result::Result<(), AromaticityError> {
    match model {
        AromaticityModel::RdkitLikeBasic => perceive_rdkit_like_basic_aromaticity(mol),
    }
}

fn perceive_rdkit_like_basic_aromaticity(
    mol: &mut Molecule,
) -> std::result::Result<(), AromaticityError> {
    for atom in mol.atoms.iter_mut().flatten() {
        atom.aromatic = false;
    }
    for bond in mol.bonds.iter_mut().flatten() {
        bond.aromatic = false;
    }

    let membership = if mol.perception.rings == ComputedState::Fresh {
        mol.ring_membership
            .clone()
            .unwrap_or_else(|| compute_ring_membership(mol))
    } else {
        perceive_ring_membership(mol)
    };

    for component in aromatic_ring_components(mol, &membership) {
        let electrons = aromatic_pi_electrons(mol, &component)?;
        if electrons >= 2 && (electrons - 2) % 4 == 0 {
            for atom_id in component.atoms {
                if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
                    atom.aromatic = true;
                }
            }
            for bond_id in component.bonds {
                if let Some(bond) = mol.bonds[bond_id.index()].as_mut() {
                    bond.aromatic = true;
                }
            }
        }
    }

    mol.perception.aromaticity = ComputedState::Fresh;
    Ok(())
}

#[derive(Debug, Clone)]
struct AromaticComponent {
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
}

fn aromatic_ring_components(mol: &Molecule, membership: &RingMembership) -> Vec<AromaticComponent> {
    let mut visited_atoms = vec![false; mol.atoms.len()];
    let mut components = Vec::new();

    for start in membership.ring_atom_ids() {
        if visited_atoms[start.index()] {
            continue;
        }
        let mut stack = vec![start];
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        visited_atoms[start.index()] = true;

        while let Some(atom_id) = stack.pop() {
            atoms.push(atom_id);
            if let Ok(incident) = mol.incident_bonds(atom_id) {
                for (bond_id, bond) in incident {
                    if !membership.bond_in_ring(bond_id) {
                        continue;
                    }
                    bonds.push(bond_id);
                    let neighbor = bond.other_atom(atom_id);
                    if !visited_atoms[neighbor.index()] {
                        visited_atoms[neighbor.index()] = true;
                        stack.push(neighbor);
                    }
                }
            }
        }
        bonds.sort();
        bonds.dedup();
        components.push(AromaticComponent { atoms, bonds });
    }

    components
}

fn aromatic_pi_electrons(
    mol: &Molecule,
    component: &AromaticComponent,
) -> std::result::Result<u8, AromaticityError> {
    let mut electrons = 0u8;
    for bond_id in &component.bonds {
        let bond = mol.bond(*bond_id).expect("component bond should be live");
        if matches!(bond.order, BondOrder::Double | BondOrder::Aromatic) {
            electrons += 2;
        }
    }

    for atom_id in &component.atoms {
        let atom = mol.atom(*atom_id).expect("component atom should be live");
        match atom.element.symbol() {
            "C" | "N" => {}
            "O" | "S" | "P" => {
                if !component_atom_has_pi_bond(mol, component, *atom_id) {
                    electrons += 2;
                }
            }
            _ => return Err(AromaticityError::UnsupportedElement(*atom_id)),
        }
    }

    Ok(electrons)
}

fn component_atom_has_pi_bond(
    mol: &Molecule,
    component: &AromaticComponent,
    atom_id: AtomId,
) -> bool {
    component.bonds.iter().any(|bond_id| {
        mol.bond(*bond_id)
            .map(|bond| {
                (bond.a == atom_id || bond.b == atom_id)
                    && matches!(bond.order, BondOrder::Double | BondOrder::Aromatic)
            })
            .unwrap_or(false)
    })
}

fn compute_ring_membership(mol: &Molecule) -> RingMembership {
    let mut graph = vec![Vec::<(AtomId, BondId)>::new(); mol.atoms.len()];
    let mut live_bonds = Vec::new();
    for (bond_id, bond) in mol.bonds() {
        graph[bond.a.index()].push((bond.b, bond_id));
        graph[bond.b.index()].push((bond.a, bond_id));
        live_bonds.push(bond_id);
    }

    let mut discovery = vec![None; mol.atoms.len()];
    let mut low = vec![0usize; mol.atoms.len()];
    let mut bridge = vec![false; mol.bonds.len()];
    let mut time = 0usize;

    for atom_id in mol.atom_ids().collect::<Vec<_>>() {
        if discovery[atom_id.index()].is_none() {
            ring_dfs(
                atom_id,
                None,
                &graph,
                &mut discovery,
                &mut low,
                &mut bridge,
                &mut time,
            );
        }
    }

    let mut membership = RingMembership {
        atom_flags: vec![false; mol.atoms.len()],
        bond_flags: vec![false; mol.bonds.len()],
    };
    for bond_id in live_bonds {
        if !bridge[bond_id.index()] {
            let bond = mol.bond(bond_id).expect("live bond should be readable");
            membership.bond_flags[bond_id.index()] = true;
            membership.atom_flags[bond.a.index()] = true;
            membership.atom_flags[bond.b.index()] = true;
        }
    }
    membership
}

fn ring_dfs(
    atom: AtomId,
    parent_bond: Option<BondId>,
    graph: &[Vec<(AtomId, BondId)>],
    discovery: &mut [Option<usize>],
    low: &mut [usize],
    bridge: &mut [bool],
    time: &mut usize,
) {
    discovery[atom.index()] = Some(*time);
    low[atom.index()] = *time;
    *time += 1;

    for (neighbor, bond_id) in &graph[atom.index()] {
        if Some(*bond_id) == parent_bond {
            continue;
        }
        if discovery[neighbor.index()].is_none() {
            ring_dfs(
                *neighbor,
                Some(*bond_id),
                graph,
                discovery,
                low,
                bridge,
                time,
            );
            low[atom.index()] = low[atom.index()].min(low[neighbor.index()]);
            if low[neighbor.index()] > discovery[atom.index()].expect("atom is discovered") {
                bridge[bond_id.index()] = true;
            }
        } else {
            low[atom.index()] =
                low[atom.index()].min(discovery[neighbor.index()].expect("neighbor discovered"));
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SmallMolecule {
    pub mol: Molecule,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SdfParseOptions {
    pub allow_missing_final_delimiter: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecord {
    pub title: String,
    pub molecule: SmallMolecule,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfParseError {
    pub record: usize,
    pub line: usize,
    pub message: String,
}

impl SdfParseError {
    fn new(record: usize, line: usize, message: impl Into<String>) -> Self {
        Self {
            record,
            line,
            message: message.into(),
        }
    }
}

impl fmt::Display for SdfParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SDF parse error in record {} at line {}: {}",
            self.record, self.line, self.message
        )
    }
}

impl std::error::Error for SdfParseError {}

pub fn read_sdf_v2000_str(
    input: &str,
    options: SdfParseOptions,
) -> std::result::Result<Vec<SmallMolecule>, SdfParseError> {
    read_sdf_v2000_records(input, options)
        .map(|records| records.into_iter().map(|record| record.molecule).collect())
}

pub fn read_sdf_v2000_records(
    input: &str,
    options: SdfParseOptions,
) -> std::result::Result<Vec<SdfRecord>, SdfParseError> {
    let normalized = input.replace("\r\n", "\n").replace('\r', "\n");
    let mut records = Vec::new();
    let mut current = Vec::new();
    let mut start_line = 1usize;
    let mut saw_delimiter = false;

    for (offset, line) in normalized.lines().enumerate() {
        let line_number = offset + 1;
        if line.trim() == "$$$$" {
            saw_delimiter = true;
            if current.iter().any(|line: &&str| !line.trim().is_empty()) {
                records.push(parse_sdf_record(records.len() + 1, start_line, &current)?);
            }
            current.clear();
            start_line = line_number + 1;
        } else {
            current.push(line);
        }
    }

    if current.iter().any(|line| !line.trim().is_empty()) {
        if saw_delimiter || options.allow_missing_final_delimiter {
            records.push(parse_sdf_record(records.len() + 1, start_line, &current)?);
        } else {
            return Err(SdfParseError::new(
                records.len() + 1,
                start_line + current.len().saturating_sub(1),
                "missing final $$$$ record delimiter",
            ));
        }
    }

    Ok(records)
}

fn parse_sdf_record(
    record: usize,
    start_line: usize,
    lines: &[&str],
) -> std::result::Result<SdfRecord, SdfParseError> {
    if lines.len() < 4 {
        return Err(SdfParseError::new(
            record,
            start_line,
            "record must contain three header lines and a counts line",
        ));
    }
    let title = lines[0].to_owned();
    let counts = lines[3];
    if counts.contains("V3000") {
        return Err(SdfParseError::new(
            record,
            start_line + 3,
            "V3000 records are not supported by the V2000 parser",
        ));
    }
    if !counts.contains("V2000") {
        return Err(SdfParseError::new(
            record,
            start_line + 3,
            "counts line must declare V2000",
        ));
    }
    let (atom_count, bond_count) = parse_counts_line(counts)
        .ok_or_else(|| SdfParseError::new(record, start_line + 3, "invalid V2000 counts line"))?;

    let mut mol = Molecule::new();
    mol.props_mut()
        .insert("sdf.title".to_owned(), PropValue::String(title.clone()));
    mol.props_mut().insert(
        "sdf.program".to_owned(),
        PropValue::String(lines[1].to_owned()),
    );
    mol.props_mut().insert(
        "sdf.comment".to_owned(),
        PropValue::String(lines[2].to_owned()),
    );

    let atom_start = 4;
    let bond_start = atom_start + atom_count;
    let property_start = bond_start + bond_count;
    if lines.len() < property_start {
        return Err(SdfParseError::new(
            record,
            start_line + lines.len(),
            "record ended before declared atom and bond blocks",
        ));
    }

    let mut atom_ids = Vec::with_capacity(atom_count);
    for atom_index in 0..atom_count {
        let line_number = start_line + atom_start + atom_index;
        let symbol = atom_symbol_from_v2000_line(lines[atom_start + atom_index])
            .ok_or_else(|| SdfParseError::new(record, line_number, "invalid atom line"))?;
        let element = Element::from_symbol(symbol).ok_or_else(|| {
            SdfParseError::new(
                record,
                line_number,
                format!("unknown element symbol `{symbol}`"),
            )
        })?;
        atom_ids.push(mol.add_atom(Atom::new(element)));
    }

    for bond_index in 0..bond_count {
        let line_number = start_line + bond_start + bond_index;
        let (a, b, order) = parse_v2000_bond_line(lines[bond_start + bond_index])
            .ok_or_else(|| SdfParseError::new(record, line_number, "invalid bond line"))?;
        let a = atom_ids.get(a - 1).copied().ok_or_else(|| {
            SdfParseError::new(record, line_number, "bond endpoint outside atom block")
        })?;
        let b = atom_ids.get(b - 1).copied().ok_or_else(|| {
            SdfParseError::new(record, line_number, "bond endpoint outside atom block")
        })?;
        mol.add_bond(a, b, order).map_err(|error| {
            SdfParseError::new(record, line_number, format!("invalid graph bond: {error}"))
        })?;
    }

    let end_index = lines[property_start..]
        .iter()
        .position(|line| line.trim() == "M  END")
        .map(|index| property_start + index)
        .ok_or_else(|| SdfParseError::new(record, start_line + property_start, "missing M  END"))?;
    parse_sdf_data_fields(record, start_line, &mut mol, &lines[end_index + 1..])?;

    Ok(SdfRecord {
        title,
        molecule: SmallMolecule { mol },
    })
}

fn parse_counts_line(line: &str) -> Option<(usize, usize)> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    let atoms = fields.first()?.parse().ok()?;
    let bonds = fields.get(1)?.parse().ok()?;
    Some((atoms, bonds))
}

fn atom_symbol_from_v2000_line(line: &str) -> Option<&str> {
    line.get(31..34)
        .map(str::trim)
        .filter(|symbol| !symbol.is_empty())
        .or_else(|| line.split_whitespace().nth(3))
}

fn parse_v2000_bond_line(line: &str) -> Option<(usize, usize, BondOrder)> {
    let mut fields = line.split_whitespace();
    let a = fields.next()?.parse().ok()?;
    let b = fields.next()?.parse().ok()?;
    let order_code: u8 = fields.next()?.parse().ok()?;
    let order = match order_code {
        0 => BondOrder::Zero,
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        9 => BondOrder::Dative,
        _ => return None,
    };
    Some((a, b, order))
}

fn parse_sdf_data_fields(
    record: usize,
    start_line: usize,
    mol: &mut Molecule,
    lines: &[&str],
) -> std::result::Result<(), SdfParseError> {
    let mut index = 0;
    while index < lines.len() {
        let line = lines[index];
        if !line.trim_start().starts_with('>') {
            index += 1;
            continue;
        }
        let field_name = sdf_field_name(line).ok_or_else(|| {
            SdfParseError::new(record, start_line + index, "invalid SDF data field header")
        })?;
        index += 1;
        let mut values = Vec::new();
        while index < lines.len() && !lines[index].trim_start().starts_with('>') {
            if lines[index].is_empty() {
                index += 1;
                break;
            }
            values.push(lines[index]);
            index += 1;
        }
        mol.props_mut().insert(
            format!("sdf.field.{field_name}"),
            PropValue::String(values.join("\n")),
        );
    }
    Ok(())
}

fn sdf_field_name(line: &str) -> Option<String> {
    let start = line.find('<')?;
    let end = line[start + 1..].find('>')? + start + 1;
    let name = line[start + 1..end].trim();
    if name.is_empty() {
        None
    } else {
        Some(name.to_owned())
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct BioHierarchy {
    pub props: PropMap,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MacroMolecule {
    pub mol: Molecule,
    pub hierarchy: BioHierarchy,
}

pub type Result<T> = std::result::Result<T, MoleculeError>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoleculeError {
    InvalidAtomId(AtomId),
    InvalidBondId(BondId),
    SelfBond(AtomId),
    DuplicateBond { a: AtomId, b: AtomId },
    UnsupportedFeature(&'static str),
}

impl fmt::Display for MoleculeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(id) => write!(f, "invalid atom id: {id}"),
            Self::InvalidBondId(id) => write!(f, "invalid bond id: {id}"),
            Self::SelfBond(id) => write!(f, "cannot create a bond from atom {id} to itself"),
            Self::DuplicateBond { a, b } => write!(f, "duplicate bond between {a} and {b}"),
            Self::UnsupportedFeature(name) => write!(f, "unsupported feature: {name}"),
        }
    }
}

impl std::error::Error for MoleculeError {}

#[cfg(test)]
mod tests {
    use super::*;

    fn carbon() -> Atom {
        Atom::new(Element::from_symbol("C").expect("carbon should be available"))
    }

    fn oxygen() -> Atom {
        Atom::new(Element::from_symbol("O").expect("oxygen should be available"))
    }

    fn ring_molecule(
        symbols: &[&str],
        orders: &[BondOrder],
    ) -> (Molecule, Vec<AtomId>, Vec<BondId>) {
        assert_eq!(symbols.len(), orders.len());
        let mut mol = Molecule::new();
        let atoms = symbols
            .iter()
            .map(|symbol| {
                mol.add_atom(Atom::new(
                    Element::from_symbol(symbol).expect("test element should be available"),
                ))
            })
            .collect::<Vec<_>>();
        let mut bonds = Vec::new();
        for index in 0..atoms.len() {
            let next = (index + 1) % atoms.len();
            bonds.push(
                mol.add_bond(atoms[index], atoms[next], orders[index])
                    .expect("ring bond should be valid"),
            );
        }
        (mol, atoms, bonds)
    }

    fn sorted_atom_ids(ids: impl IntoIterator<Item = AtomId>) -> Vec<AtomId> {
        let mut ids = ids.into_iter().collect::<Vec<_>>();
        ids.sort();
        ids
    }

    fn sorted_bond_ids(ids: impl IntoIterator<Item = BondId>) -> Vec<BondId> {
        let mut ids = ids.into_iter().collect::<Vec<_>>();
        ids.sort();
        ids
    }

    fn mark_all_fresh(mol: &mut Molecule) {
        mol.perception_mut().valence = ComputedState::Fresh;
        mol.perception_mut().rings = ComputedState::Fresh;
        mol.perception_mut().aromaticity = ComputedState::Fresh;
        mol.perception_mut().stereo = ComputedState::Fresh;
    }

    fn assert_all_stale(mol: &Molecule) {
        assert_eq!(mol.perception().valence, ComputedState::Stale);
        assert_eq!(mol.perception().rings, ComputedState::Stale);
        assert_eq!(mol.perception().aromaticity, ComputedState::Stale);
        assert_eq!(mol.perception().stereo, ComputedState::Stale);
    }

    #[test]
    fn element_from_atomic_number_accepts_periodic_table_bounds() {
        assert_eq!(
            Element::from_atomic_number(1)
                .expect("hydrogen exists")
                .symbol(),
            "H"
        );
        assert_eq!(
            Element::from_atomic_number(118)
                .expect("oganesson exists")
                .symbol(),
            "Og"
        );
    }

    #[test]
    fn element_from_atomic_number_rejects_out_of_range_values() {
        assert_eq!(Element::from_atomic_number(0), None);
        assert_eq!(Element::from_atomic_number(119), None);
    }

    #[test]
    fn element_from_symbol_is_canonical_and_case_sensitive() {
        assert_eq!(
            Element::from_symbol("C")
                .expect("carbon exists")
                .atomic_number(),
            6
        );
        assert_eq!(
            Element::from_symbol("Cl")
                .expect("chlorine exists")
                .atomic_number(),
            17
        );
        assert_eq!(
            Element::from_symbol("Og")
                .expect("oganesson exists")
                .atomic_number(),
            118
        );
        assert_eq!(Element::from_symbol("CL"), None);
        assert_eq!(Element::from_symbol("Xx"), None);
        assert_eq!(Element::from_symbol("?"), None);
    }

    #[test]
    fn element_symbol_and_display_are_canonical() {
        let iron = Element::from_atomic_number(26).expect("iron exists");

        assert_eq!(iron.symbol(), "Fe");
        assert_eq!(iron.to_string(), "Fe");
    }

    #[test]
    fn atom_new_sets_chemically_general_defaults() {
        let atom = carbon();

        assert_eq!(atom.element.symbol(), "C");
        assert_eq!(atom.isotope, None);
        assert_eq!(atom.formal_charge, 0);
        assert_eq!(atom.radical_electrons, 0);
        assert_eq!(atom.explicit_hydrogens, 0);
        assert_eq!(atom.implicit_hydrogens, None);
        assert!(!atom.aromatic);
        assert_eq!(atom.chiral, None);
        assert_eq!(atom.atom_map, None);
        assert!(atom.props.is_empty());
    }

    #[test]
    fn atom_payload_fields_can_be_set_and_read() {
        let mut atom = carbon();
        atom.isotope = Some(13);
        atom.formal_charge = -1;
        atom.radical_electrons = 1;
        atom.explicit_hydrogens = 3;
        atom.implicit_hydrogens = Some(1);
        atom.aromatic = true;
        atom.chiral = Some(AtomStereo::TetrahedralClockwise);
        atom.atom_map = Some(7);
        atom.props
            .insert("label".to_owned(), PropValue::String("alpha".to_owned()));

        assert_eq!(atom.isotope, Some(13));
        assert_eq!(atom.formal_charge, -1);
        assert_eq!(atom.radical_electrons, 1);
        assert_eq!(atom.explicit_hydrogens, 3);
        assert_eq!(atom.implicit_hydrogens, Some(1));
        assert!(atom.aromatic);
        assert_eq!(atom.chiral, Some(AtomStereo::TetrahedralClockwise));
        assert_eq!(atom.atom_map, Some(7));
        assert_eq!(
            atom.props.get("label"),
            Some(&PropValue::String("alpha".to_owned()))
        );
    }

    #[test]
    fn bond_new_sets_endpoints_order_and_aromatic_default() {
        let a = AtomId::new(3);
        let b = AtomId::new(4);
        let single = Bond::new(a, b, BondOrder::Single);
        let aromatic = Bond::new(a, b, BondOrder::Aromatic);

        assert_eq!(single.a(), a);
        assert_eq!(single.b(), b);
        assert_eq!(single.endpoints(), (a, b));
        assert_eq!(single.order, BondOrder::Single);
        assert!(!single.aromatic);
        assert_eq!(single.stereo, None);
        assert!(single.props.is_empty());
        assert!(aromatic.aromatic);
    }

    #[test]
    fn bond_payload_fields_can_be_set_and_read() {
        let mut bond = Bond::new(AtomId::new(1), AtomId::new(2), BondOrder::Dative);
        bond.stereo = Some(BondStereo::Up);
        bond.props
            .insert("score".to_owned(), PropValue::Float(1.25));

        assert_eq!(bond.order, BondOrder::Dative);
        assert_eq!(bond.stereo, Some(BondStereo::Up));
        assert_eq!(bond.props.get("score"), Some(&PropValue::Float(1.25)));
    }

    #[test]
    fn prop_value_equality_covers_all_initial_variants() {
        assert_eq!(
            PropValue::String("value".to_owned()),
            PropValue::String("value".to_owned())
        );
        assert_eq!(PropValue::Int(42), PropValue::Int(42));
        assert_eq!(PropValue::Float(2.5), PropValue::Float(2.5));
        assert_eq!(PropValue::Bool(true), PropValue::Bool(true));
    }

    #[test]
    fn mutable_payload_access_invalidates_fresh_perception() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");

        mark_all_fresh(&mut mol);
        mol.atom_mut(a).expect("atom exists").formal_charge = 1;
        assert_all_stale(&mol);

        mark_all_fresh(&mut mol);
        mol.bond_mut(bond).expect("bond exists").order = BondOrder::Double;
        assert_all_stale(&mol);
    }

    #[test]
    fn sdf_v2000_parses_single_record_atoms_bonds_and_fields() {
        let input = "\
Water
  molecules
comment
  2  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0
    1.0000    0.0000    0.0000 H   0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
>  <NAME>
water

$$$$
";

        let molecules =
            read_sdf_v2000_str(input, SdfParseOptions::default()).expect("record should parse");
        let mol = &molecules[0].mol;

        assert_eq!(molecules.len(), 1);
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(
            mol.atom(AtomId::new(0))
                .expect("atom exists")
                .element
                .symbol(),
            "O"
        );
        assert_eq!(
            mol.bond(BondId::new(0)).expect("bond exists").order,
            BondOrder::Single
        );
        assert_eq!(
            mol.props().get("sdf.field.NAME"),
            Some(&PropValue::String("water".to_owned()))
        );
    }

    #[test]
    fn sdf_v2000_parses_multiple_records_in_order() {
        let input = "\
One
  molecules

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
$$$$
Two
  molecules

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0
M  END
$$$$
";

        let records = read_sdf_v2000_records(input, SdfParseOptions::default())
            .expect("records should parse");

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].title, "One");
        assert_eq!(records[1].title, "Two");
        assert_eq!(
            records[1]
                .molecule
                .mol
                .atom(AtomId::new(0))
                .expect("atom exists")
                .element
                .symbol(),
            "O"
        );
    }

    #[test]
    fn sdf_v2000_can_allow_missing_final_delimiter() {
        let input = "\
Methane
  molecules

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
M  END
";

        let molecules = read_sdf_v2000_str(
            input,
            SdfParseOptions {
                allow_missing_final_delimiter: true,
            },
        )
        .expect("record should parse");

        assert_eq!(molecules.len(), 1);
        assert_eq!(molecules[0].mol.atom_count(), 1);
    }

    #[test]
    fn sdf_v2000_rejects_v3000_and_bad_endpoints() {
        let v3000 = "\
V3000
  molecules

  0  0  0  0  0  0            999 V3000
M  END
$$$$
";
        let err =
            read_sdf_v2000_str(v3000, SdfParseOptions::default()).expect_err("V3000 should fail");
        assert!(err.message.contains("V3000"));

        let bad_endpoint = "\
Bad
  molecules

  1  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
$$$$
";
        let err = read_sdf_v2000_str(bad_endpoint, SdfParseOptions::default())
            .expect_err("bad endpoint should fail");
        assert!(err.message.contains("outside atom block"));
    }

    #[test]
    fn sdf_v2000_parse_does_not_perceive_chemistry() {
        let input = "\
Benzene-ish
  molecules

  2  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0
  1  2  4  0  0  0  0
M  END
$$$$
";

        let molecules =
            read_sdf_v2000_str(input, SdfParseOptions::default()).expect("record should parse");
        let mol = &molecules[0].mol;

        assert_eq!(mol.perception().rings, ComputedState::Absent);
        assert_eq!(mol.perception().aromaticity, ComputedState::Absent);
        assert_eq!(
            mol.bond(BondId::new(0)).expect("bond exists").order,
            BondOrder::Aromatic
        );
    }

    #[test]
    fn ring_membership_empty_and_linear_molecules_have_no_rings() {
        let mut empty = Molecule::new();
        let empty_membership = perceive_ring_membership(&mut empty);
        assert!(empty_membership.ring_atom_ids().next().is_none());
        assert!(empty_membership.ring_bond_ids().next().is_none());

        let mut chain = Molecule::new();
        let a = chain.add_atom(carbon());
        let b = chain.add_atom(carbon());
        let c = chain.add_atom(carbon());
        let ab = chain
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");
        let bc = chain
            .add_bond(b, c, BondOrder::Single)
            .expect("bond should be valid");
        let chain_membership = perceive_ring_membership(&mut chain);

        assert!(!chain_membership.atom_in_ring(a));
        assert!(!chain_membership.atom_in_ring(b));
        assert!(!chain_membership.bond_in_ring(ab));
        assert!(!chain_membership.bond_in_ring(bc));
        assert_eq!(chain.perception().rings, ComputedState::Fresh);
    }

    #[test]
    fn ring_membership_marks_triangle_atoms_and_bonds() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let c = mol.add_atom(carbon());
        let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
        let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
        let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");

        let membership = perceive_ring_membership(&mut mol);

        assert_eq!(sorted_atom_ids(membership.ring_atom_ids()), vec![a, b, c]);
        assert_eq!(
            sorted_bond_ids(membership.ring_bond_ids()),
            vec![ab, bc, ca]
        );
    }

    #[test]
    fn ring_membership_excludes_tail_from_ring() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let c = mol.add_atom(carbon());
        let tail = mol.add_atom(oxygen());
        let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
        let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
        let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");
        let tail_bond = mol.add_bond(c, tail, BondOrder::Single).expect("bond");

        let membership = perceive_ring_membership(&mut mol);

        assert_eq!(sorted_atom_ids(membership.ring_atom_ids()), vec![a, b, c]);
        assert_eq!(
            sorted_bond_ids(membership.ring_bond_ids()),
            vec![ab, bc, ca]
        );
        assert!(!membership.atom_in_ring(tail));
        assert!(!membership.bond_in_ring(tail_bond));
    }

    #[test]
    fn ring_membership_handles_fused_and_disconnected_components() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let c = mol.add_atom(carbon());
        let d = mol.add_atom(carbon());
        let isolated_a = mol.add_atom(oxygen());
        let isolated_b = mol.add_atom(oxygen());
        let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
        let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
        let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");
        let cd = mol.add_bond(c, d, BondOrder::Single).expect("bond");
        let da = mol.add_bond(d, a, BondOrder::Single).expect("bond");
        let bridge = mol
            .add_bond(isolated_a, isolated_b, BondOrder::Single)
            .expect("bond");

        let membership = perceive_ring_membership(&mut mol);

        assert_eq!(
            sorted_atom_ids(membership.ring_atom_ids()),
            vec![a, b, c, d]
        );
        assert_eq!(
            sorted_bond_ids(membership.ring_bond_ids()),
            vec![ab, bc, ca, cd, da]
        );
        assert!(!membership.bond_in_ring(bridge));
    }

    #[test]
    fn ring_membership_ignores_deleted_bonds_and_becomes_stale_after_mutation() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let c = mol.add_atom(carbon());
        let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
        let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
        let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");
        mol.delete_bond(ca).expect("bond should delete");

        let membership = perceive_ring_membership(&mut mol);
        assert!(!membership.bond_in_ring(ab));
        assert!(!membership.bond_in_ring(bc));
        assert!(!membership.bond_in_ring(ca));

        mol.add_bond(c, a, BondOrder::Single).expect("bond");
        assert_eq!(mol.perception().rings, ComputedState::Stale);
        assert!(mol.ring_membership().is_some());
    }

    #[test]
    fn aromaticity_marks_benzene_like_ring() {
        let (mut mol, atoms, bonds) = ring_molecule(
            &["C", "C", "C", "C", "C", "C"],
            &[
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
            ],
        );

        perceive_aromaticity(&mut mol, AromaticityModel::RdkitLikeBasic)
            .expect("benzene should be supported");

        assert_eq!(mol.perception().aromaticity, ComputedState::Fresh);
        assert!(atoms
            .iter()
            .all(|atom| mol.atom(*atom).expect("atom exists").aromatic));
        assert!(bonds
            .iter()
            .all(|bond| mol.bond(*bond).expect("bond exists").aromatic));
    }

    #[test]
    fn aromaticity_leaves_cyclohexane_and_cyclobutadiene_non_aromatic() {
        let (mut cyclohexane, atoms, bonds) =
            ring_molecule(&["C", "C", "C", "C", "C", "C"], &[BondOrder::Single; 6]);
        perceive_aromaticity(&mut cyclohexane, AromaticityModel::RdkitLikeBasic)
            .expect("cyclohexane should be supported");
        assert!(atoms
            .iter()
            .all(|atom| !cyclohexane.atom(*atom).expect("atom exists").aromatic));
        assert!(bonds
            .iter()
            .all(|bond| !cyclohexane.bond(*bond).expect("bond exists").aromatic));

        let (mut cyclobutadiene, atoms, bonds) = ring_molecule(
            &["C", "C", "C", "C"],
            &[
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
            ],
        );
        perceive_aromaticity(&mut cyclobutadiene, AromaticityModel::RdkitLikeBasic)
            .expect("cyclobutadiene should be supported");
        assert!(atoms
            .iter()
            .all(|atom| !cyclobutadiene.atom(*atom).expect("atom exists").aromatic));
        assert!(bonds
            .iter()
            .all(|bond| !cyclobutadiene.bond(*bond).expect("bond exists").aromatic));
    }

    #[test]
    fn aromaticity_supports_basic_heteroaromatic_ring() {
        let (mut furan_like, atoms, bonds) = ring_molecule(
            &["O", "C", "C", "C", "C"],
            &[
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
            ],
        );

        perceive_aromaticity(&mut furan_like, AromaticityModel::RdkitLikeBasic)
            .expect("furan-like ring should be supported");

        assert!(atoms
            .iter()
            .all(|atom| furan_like.atom(*atom).expect("atom exists").aromatic));
        assert!(bonds
            .iter()
            .all(|bond| furan_like.bond(*bond).expect("bond exists").aromatic));
    }

    #[test]
    fn aromaticity_uses_ring_membership_not_acyclic_double_bonds() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let c = mol.add_atom(carbon());
        mol.add_bond(a, b, BondOrder::Double).expect("bond");
        mol.add_bond(b, c, BondOrder::Single).expect("bond");

        perceive_aromaticity(&mut mol, AromaticityModel::RdkitLikeBasic)
            .expect("acyclic molecule should be supported");

        assert!(!mol.atom(a).expect("atom exists").aromatic);
        assert!(!mol.bond(BondId::new(0)).expect("bond exists").aromatic);
    }

    #[test]
    fn aromaticity_clears_existing_flags_before_assignment() {
        let (mut mol, atoms, bonds) =
            ring_molecule(&["C", "C", "C", "C", "C", "C"], &[BondOrder::Single; 6]);
        for atom in &atoms {
            mol.atom_mut(*atom).expect("atom exists").aromatic = true;
        }
        for bond in &bonds {
            mol.bond_mut(*bond).expect("bond exists").aromatic = true;
        }

        perceive_aromaticity(&mut mol, AromaticityModel::RdkitLikeBasic)
            .expect("cyclohexane should be supported");

        assert!(atoms
            .iter()
            .all(|atom| !mol.atom(*atom).expect("atom exists").aromatic));
        assert!(bonds
            .iter()
            .all(|bond| !mol.bond(*bond).expect("bond exists").aromatic));
    }

    #[test]
    fn aromaticity_becomes_stale_after_topology_mutation() {
        let (mut mol, atoms, _) = ring_molecule(
            &["C", "C", "C", "C", "C", "C"],
            &[
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
            ],
        );
        perceive_aromaticity(&mut mol, AromaticityModel::RdkitLikeBasic)
            .expect("benzene should be supported");

        mol.add_atom(oxygen());
        assert_eq!(mol.perception().aromaticity, ComputedState::Stale);
        assert!(atoms
            .iter()
            .all(|atom| mol.atom(*atom).expect("atom exists").aromatic));
    }

    #[test]
    fn empty_molecule_has_no_atoms_or_bonds() {
        let mol = Molecule::new();

        assert_eq!(mol.atom_count(), 0);
        assert_eq!(mol.bond_count(), 0);
        assert!(mol.atoms().next().is_none());
        assert!(mol.bonds().next().is_none());
    }

    #[test]
    fn atom_insertion_assigns_stable_typed_ids() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());

        assert_eq!(a.raw(), 0);
        assert_eq!(b.raw(), 1);
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(
            mol.atom(a).expect("first atom exists").element.symbol(),
            "C"
        );
        assert_eq!(
            mol.atom(b).expect("second atom exists").element.symbol(),
            "O"
        );
        assert_eq!(mol.atom_ids().collect::<Vec<_>>(), vec![a, b]);
    }

    #[test]
    fn bond_insertion_assigns_stable_typed_ids() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");

        assert_eq!(bond.raw(), 0);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(
            mol.bond(bond).expect("bond should exist").endpoints(),
            (a, b)
        );
        assert_eq!(mol.bond_ids().collect::<Vec<_>>(), vec![bond]);
    }

    #[test]
    fn invalid_atom_ids_are_rejected() {
        let mut mol = Molecule::new();
        let atom = mol.add_atom(carbon());

        assert_eq!(
            mol.atom(AtomId::new(99))
                .expect_err("missing atom should fail"),
            MoleculeError::InvalidAtomId(AtomId::new(99))
        );
        mol.delete_atom(atom).expect("atom should delete");
        assert_eq!(
            mol.atom(atom).expect_err("deleted atom should fail"),
            MoleculeError::InvalidAtomId(atom)
        );
        assert_eq!(
            mol.add_bond(atom, AtomId::new(99), BondOrder::Single)
                .expect_err("deleted endpoint should fail"),
            MoleculeError::InvalidAtomId(atom)
        );
    }

    #[test]
    fn invalid_bond_ids_are_rejected() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");

        assert_eq!(
            mol.bond(BondId::new(99))
                .expect_err("missing bond should fail"),
            MoleculeError::InvalidBondId(BondId::new(99))
        );
        mol.delete_bond(bond).expect("bond should delete");
        assert_eq!(
            mol.bond(bond).expect_err("deleted bond should fail"),
            MoleculeError::InvalidBondId(bond)
        );
        assert_eq!(
            mol.delete_bond(bond)
                .expect_err("deleting bond twice should fail"),
            MoleculeError::InvalidBondId(bond)
        );
    }

    #[test]
    fn self_bonds_are_rejected() {
        let mut mol = Molecule::new();
        let atom = mol.add_atom(carbon());

        let err = mol
            .add_bond(atom, atom, BondOrder::Single)
            .expect_err("self-bond should fail");
        assert_eq!(err, MoleculeError::SelfBond(atom));
    }

    #[test]
    fn duplicate_bond_is_rejected() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        mol.add_bond(a, b, BondOrder::Single)
            .expect("first bond should be valid");

        let err = mol
            .add_bond(a, b, BondOrder::Double)
            .expect_err("duplicate should fail");
        assert_eq!(err, MoleculeError::DuplicateBond { a, b });

        let reverse_err = mol
            .add_bond(b, a, BondOrder::Double)
            .expect_err("reverse duplicate should fail");
        assert_eq!(reverse_err, MoleculeError::DuplicateBond { a: b, b: a });
    }

    #[test]
    fn neighbor_iteration_reports_live_adjacent_atoms() {
        let mut mol = Molecule::new();
        let center = mol.add_atom(carbon());
        let left = mol.add_atom(carbon());
        let right = mol.add_atom(oxygen());
        let isolated = mol.add_atom(carbon());
        mol.add_bond(center, left, BondOrder::Single)
            .expect("left bond should be valid");
        mol.add_bond(center, right, BondOrder::Double)
            .expect("right bond should be valid");

        assert_eq!(
            sorted_atom_ids(mol.neighbors(center).expect("center exists")),
            vec![left, right]
        );
        assert_eq!(
            mol.neighbors(isolated)
                .expect("isolated atom exists")
                .collect::<Vec<_>>(),
            Vec::<AtomId>::new()
        );
        match mol.neighbors(AtomId::new(99)) {
            Ok(_) => panic!("missing atom should fail"),
            Err(err) => assert_eq!(err, MoleculeError::InvalidAtomId(AtomId::new(99))),
        };
    }

    #[test]
    fn incident_bond_iteration_reports_live_bonds() {
        let mut mol = Molecule::new();
        let center = mol.add_atom(carbon());
        let left = mol.add_atom(carbon());
        let right = mol.add_atom(oxygen());
        let left_bond = mol
            .add_bond(center, left, BondOrder::Single)
            .expect("left bond should be valid");
        let right_bond = mol
            .add_bond(center, right, BondOrder::Double)
            .expect("right bond should be valid");

        assert_eq!(
            sorted_bond_ids(
                mol.incident_bonds(center)
                    .expect("center exists")
                    .map(|(id, _)| id)
            ),
            vec![left_bond, right_bond]
        );

        mol.delete_bond(left_bond).expect("left bond should delete");
        assert_eq!(
            mol.incident_bonds(center)
                .expect("center still exists")
                .map(|(id, _)| id)
                .collect::<Vec<_>>(),
            vec![right_bond]
        );
        match mol.incident_bonds(AtomId::new(99)) {
            Ok(_) => panic!("missing atom should fail"),
            Err(err) => assert_eq!(err, MoleculeError::InvalidAtomId(AtomId::new(99))),
        };
    }

    #[test]
    fn bond_between_finds_live_undirected_bonds() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());
        let c = mol.add_atom(carbon());
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");

        assert_eq!(mol.bond_between(a, b).expect("atoms exist"), Some(bond));
        assert_eq!(mol.bond_between(b, a).expect("atoms exist"), Some(bond));
        assert_eq!(mol.bond_between(a, c).expect("atoms exist"), None);
    }

    #[test]
    fn bond_deletion_preserves_remaining_ids_and_counts() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let c = mol.add_atom(oxygen());
        let first = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("first bond should be valid");
        let second = mol
            .add_bond(b, c, BondOrder::Double)
            .expect("second bond should be valid");

        let removed = mol.delete_bond(first).expect("first bond should delete");

        assert_eq!(removed.a(), a);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(mol.bond(first), Err(MoleculeError::InvalidBondId(first)));
        assert_eq!(
            mol.bond(second).expect("second bond remains").order,
            BondOrder::Double
        );
        assert_eq!(mol.bond_ids().collect::<Vec<_>>(), vec![second]);
        assert_eq!(
            mol.neighbors(b)
                .expect("middle atom exists")
                .collect::<Vec<_>>(),
            vec![c]
        );
    }

    #[test]
    fn atom_deletion_removes_incident_bonds_and_preserves_remaining_ids() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let c = mol.add_atom(oxygen());
        let first = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("first bond should be valid");
        let second = mol
            .add_bond(b, c, BondOrder::Double)
            .expect("second bond should be valid");

        let removed = mol.delete_atom(b).expect("middle atom should delete");

        assert_eq!(removed.element.symbol(), "C");
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 0);
        assert_eq!(mol.atom(b), Err(MoleculeError::InvalidAtomId(b)));
        assert_eq!(
            mol.atom(a).expect("first atom remains").element.symbol(),
            "C"
        );
        assert_eq!(
            mol.atom(c).expect("third atom remains").element.symbol(),
            "O"
        );
        assert_eq!(mol.bond(first), Err(MoleculeError::InvalidBondId(first)));
        assert_eq!(mol.bond(second), Err(MoleculeError::InvalidBondId(second)));
        assert_eq!(mol.atom_ids().collect::<Vec<_>>(), vec![a, c]);
        assert_eq!(
            mol.neighbors(a)
                .expect("first atom exists")
                .collect::<Vec<_>>(),
            Vec::<AtomId>::new()
        );
    }

    #[test]
    fn adding_after_deletion_allocates_new_ids() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        let first_bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");
        mol.delete_bond(first_bond).expect("bond should delete");
        mol.delete_atom(a).expect("atom should delete");

        let c = mol.add_atom(oxygen());
        let second_bond = mol
            .add_bond(b, c, BondOrder::Double)
            .expect("new bond should be valid");

        assert_eq!(c.raw(), 2);
        assert_eq!(second_bond.raw(), 1);
        assert_eq!(mol.atom_ids().collect::<Vec<_>>(), vec![b, c]);
        assert_eq!(mol.bond_ids().collect::<Vec<_>>(), vec![second_bond]);
    }

    #[test]
    fn every_topology_mutation_invalidates_fresh_perception() {
        let mut mol = Molecule::new();
        mark_all_fresh(&mut mol);
        let a = mol.add_atom(carbon());
        assert_all_stale(&mol);

        mark_all_fresh(&mut mol);
        let b = mol.add_atom(oxygen());
        assert_all_stale(&mol);

        mark_all_fresh(&mut mol);
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");
        assert_all_stale(&mol);

        mark_all_fresh(&mut mol);
        mol.delete_bond(bond).expect("bond should delete");
        assert_all_stale(&mol);

        mark_all_fresh(&mut mol);
        mol.delete_atom(a).expect("atom should delete");
        assert_all_stale(&mol);
    }

    #[test]
    fn absent_perception_remains_absent_after_topology_mutation() {
        let mut mol = Molecule::new();

        mol.add_atom(carbon());

        assert_eq!(mol.perception().valence, ComputedState::Absent);
        assert_eq!(mol.perception().rings, ComputedState::Absent);
        assert_eq!(mol.perception().aromaticity, ComputedState::Absent);
        assert_eq!(mol.perception().stereo, ComputedState::Absent);
    }

    #[test]
    fn property_maps_can_be_mutated_without_topology_changes() {
        let mut mol = Molecule::new();
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(oxygen());
        let bond = mol
            .add_bond(a, b, BondOrder::Single)
            .expect("bond should be valid");
        mol.props_mut().insert(
            "name".to_owned(),
            PropValue::String("carbon monoxide".to_owned()),
        );
        mol.atom_mut(a)
            .expect("atom exists")
            .props
            .insert("role".to_owned(), PropValue::String("donor".to_owned()));
        mol.bond_mut(bond)
            .expect("bond exists")
            .props
            .insert("source".to_owned(), PropValue::Bool(true));

        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(
            mol.props().get("name"),
            Some(&PropValue::String("carbon monoxide".to_owned()))
        );
        assert_eq!(
            mol.atom(a).expect("atom exists").props.get("role"),
            Some(&PropValue::String("donor".to_owned()))
        );
        assert_eq!(
            mol.bond(bond).expect("bond exists").props.get("source"),
            Some(&PropValue::Bool(true))
        );
    }

    #[test]
    fn wrappers_share_the_core_molecule_graph() {
        let mut small = SmallMolecule::default();
        let a = small.mol.add_atom(carbon());
        let b = small.mol.add_atom(oxygen());
        small
            .mol
            .add_bond(a, b, BondOrder::Single)
            .expect("small molecule graph should accept bonds");

        let mut macro_mol = MacroMolecule::default();
        let c = macro_mol.mol.add_atom(carbon());

        assert_eq!(small.mol.atom_count(), 2);
        assert_eq!(small.mol.bond_count(), 1);
        assert_eq!(
            macro_mol
                .mol
                .atom(c)
                .expect("macro atom exists")
                .element
                .symbol(),
            "C"
        );
    }
}
