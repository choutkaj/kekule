#![forbid(unsafe_code)]

use std::collections::BTreeMap;
use std::fmt;

pub mod prelude {
    pub use crate::{
        perceive_aromaticity, perceive_ring_membership, perceive_ring_set, perceive_valence,
        read_mmcif_str, read_mol_v2000_str, read_sdf_v2000_str, read_smiles_str,
        sanitize_small_molecule, write_mol_v2000, write_sdf_v2000, write_smiles, AromaticityError,
        AromaticityModel, Atom, AtomId, AtomSite, AtomSiteId, AtomSiteMetadata, AtomStereo,
        BioHierarchy, BioHierarchyError, Bond, BondId, BondOrder, BondStereo, Chain, ChainId,
        ComputedState, Conformer, ConformerId, Element, MacroMolecule, MmcifParseError,
        MmcifParseOptions, Model, ModelId, MolWriteError, Molecule, MoleculeError, Point3, PropMap,
        PropValue, Residue, ResidueId, Result, Ring, RingMembership, RingSet, SanitizeError,
        SanitizeOptions, SanitizeReport, SdfParseError, SdfParseOptions, SdfRecord, SmallMolecule,
        SmilesParseError, SmilesParseOptions, SmilesWriteOptions, ValenceIssue, ValenceModel,
        ValenceReport,
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
pub struct ConformerId(u32);

impl ConformerId {
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

impl fmt::Display for ConformerId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "c{}", self.0)
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Conformer {
    positions: Vec<Option<Point3>>,
}

impl Conformer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_atom_capacity(atom_capacity: usize) -> Self {
        Self {
            positions: vec![None; atom_capacity],
        }
    }

    pub fn set_position(&mut self, atom: AtomId, point: Point3) {
        if self.positions.len() <= atom.index() {
            self.positions.resize(atom.index() + 1, None);
        }
        self.positions[atom.index()] = Some(point);
    }

    pub fn clear_position(&mut self, atom: AtomId) {
        if let Some(position) = self.positions.get_mut(atom.index()) {
            *position = None;
        }
    }

    pub fn position(&self, atom: AtomId) -> Option<Point3> {
        self.positions.get(atom.index()).copied().flatten()
    }

    pub fn positions(&self) -> impl Iterator<Item = (AtomId, Point3)> + '_ {
        self.positions
            .iter()
            .enumerate()
            .filter_map(|(index, point)| point.map(|point| (AtomId::new(index as u32), point)))
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
    conformers: Vec<Option<Conformer>>,
    props: PropMap,
    perception: PerceptionState,
    ring_membership: Option<RingMembership>,
    ring_set: Option<RingSet>,
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
        for conformer in self.conformers.iter_mut().flatten() {
            conformer.clear_position(id);
        }
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

    pub fn ring_set(&self) -> Option<&RingSet> {
        self.ring_set.as_ref()
    }

    pub fn add_conformer(&mut self, mut conformer: Conformer) -> ConformerId {
        if conformer.positions.len() < self.atoms.len() {
            conformer.positions.resize(self.atoms.len(), None);
        }
        let id = ConformerId::new(self.conformers.len() as u32);
        self.conformers.push(Some(conformer));
        id
    }

    pub fn conformer(&self, id: ConformerId) -> Result<&Conformer> {
        self.conformers
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidConformerId(id))
    }

    pub fn conformer_mut(&mut self, id: ConformerId) -> Result<&mut Conformer> {
        self.conformers
            .get_mut(id.index())
            .and_then(Option::as_mut)
            .ok_or(MoleculeError::InvalidConformerId(id))
    }

    pub fn conformers(&self) -> impl Iterator<Item = (ConformerId, &Conformer)> {
        self.conformers
            .iter()
            .enumerate()
            .filter_map(|(index, conformer)| {
                conformer
                    .as_ref()
                    .map(|conformer| (ConformerId::new(index as u32), conformer))
            })
    }

    pub fn first_conformer(&self) -> Option<(ConformerId, &Conformer)> {
        self.conformers().next()
    }

    pub fn invalidate_topology(&mut self) {
        self.perception.invalidate_all();
        self.ring_set = None;
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
    RdkitLike,
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
        AromaticityModel::RdkitLike => perceive_rdkit_like_aromaticity(mol),
    }
}

fn perceive_rdkit_like_aromaticity(
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingSet {
    rings: Vec<Ring>,
}

impl RingSet {
    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }

    pub fn len(&self) -> usize {
        self.rings.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rings.is_empty()
    }
}

pub fn perceive_ring_set(mol: &mut Molecule) -> RingSet {
    let membership = perceive_ring_membership(mol);
    let mut graph = BTreeMap::<AtomId, Vec<(AtomId, BondId)>>::new();
    for (bond_id, bond) in mol.bonds() {
        if membership.bond_in_ring(bond_id) {
            graph.entry(bond.a()).or_default().push((bond.b(), bond_id));
            graph.entry(bond.b()).or_default().push((bond.a(), bond_id));
        }
    }
    for edges in graph.values_mut() {
        edges.sort_by_key(|(atom, bond)| (*atom, *bond));
    }

    let mut parent = BTreeMap::<AtomId, (AtomId, BondId)>::new();
    let mut depth = BTreeMap::<AtomId, usize>::new();
    let mut visited_edges = BTreeMap::<(AtomId, AtomId), ()>::new();
    let mut rings = Vec::new();

    for start in mol.atom_ids() {
        if !graph.contains_key(&start) || depth.contains_key(&start) {
            continue;
        }
        depth.insert(start, 0);
        let mut stack = vec![start];
        while let Some(atom) = stack.pop() {
            for (neighbor, bond_id) in graph.get(&atom).into_iter().flatten().copied() {
                let edge_key = ordered_atom_pair(atom, neighbor);
                if visited_edges.insert(edge_key, ()).is_some() {
                    continue;
                }
                if !depth.contains_key(&neighbor) {
                    parent.insert(neighbor, (atom, bond_id));
                    depth.insert(neighbor, depth[&atom] + 1);
                    stack.push(neighbor);
                } else if parent.get(&atom).map(|(p, _)| *p) != Some(neighbor)
                    && parent.get(&neighbor).map(|(p, _)| *p) != Some(atom)
                {
                    if let Some(ring) = fundamental_cycle(atom, neighbor, bond_id, &parent, &depth)
                    {
                        if !rings
                            .iter()
                            .any(|existing: &Ring| same_ring(existing, &ring))
                        {
                            rings.push(ring);
                        }
                    }
                }
            }
        }
    }

    rings.sort_by_key(|ring| (ring.atoms.len(), ring.atoms.clone()));
    let cyclomatic = mol.bond_count().saturating_add(connected_components(mol)) - mol.atom_count();
    rings.truncate(cyclomatic);
    let ring_set = RingSet { rings };
    mol.ring_set = Some(ring_set.clone());
    ring_set
}

fn ordered_atom_pair(a: AtomId, b: AtomId) -> (AtomId, AtomId) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

fn fundamental_cycle(
    a: AtomId,
    b: AtomId,
    closing_bond: BondId,
    parent: &BTreeMap<AtomId, (AtomId, BondId)>,
    depth: &BTreeMap<AtomId, usize>,
) -> Option<Ring> {
    let mut left = a;
    let mut right = b;
    let mut left_atoms = vec![left];
    let mut right_atoms = vec![right];
    let mut left_bonds = Vec::new();
    let mut right_bonds = Vec::new();

    while left != right {
        let left_depth = *depth.get(&left)?;
        let right_depth = *depth.get(&right)?;
        if left_depth >= right_depth {
            let (next, bond) = *parent.get(&left)?;
            left_bonds.push(bond);
            left = next;
            left_atoms.push(left);
        } else {
            let (next, bond) = *parent.get(&right)?;
            right_bonds.push(bond);
            right = next;
            right_atoms.push(right);
        }
    }

    right_atoms.reverse();
    right_bonds.reverse();
    let mut atoms = left_atoms;
    atoms.extend(right_atoms.into_iter().skip(1));
    atoms.sort();
    atoms.dedup();
    let mut bonds = left_bonds;
    bonds.extend(right_bonds);
    bonds.push(closing_bond);
    bonds.sort();
    bonds.dedup();
    Some(Ring { atoms, bonds })
}

fn same_ring(a: &Ring, b: &Ring) -> bool {
    a.atoms == b.atoms && a.bonds == b.bonds
}

fn connected_components(mol: &Molecule) -> usize {
    let mut seen = BTreeMap::<AtomId, ()>::new();
    let mut count = 0;
    for start in mol.atom_ids() {
        if seen.contains_key(&start) {
            continue;
        }
        count += 1;
        let mut stack = vec![start];
        seen.insert(start, ());
        while let Some(atom) = stack.pop() {
            if let Ok(neighbors) = mol.neighbors(atom) {
                for neighbor in neighbors {
                    if seen.insert(neighbor, ()).is_none() {
                        stack.push(neighbor);
                    }
                }
            }
        }
    }
    count
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValenceModel {
    RdkitLike,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValenceIssue {
    UnsupportedElement(AtomId),
    ValenceExceeded {
        atom: AtomId,
        explicit_valence: u8,
        max_allowed: u8,
    },
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ValenceReport {
    pub issues: Vec<ValenceIssue>,
}

impl ValenceReport {
    pub fn is_ok(&self) -> bool {
        self.issues.is_empty()
    }
}

pub fn perceive_valence(mol: &mut Molecule, model: ValenceModel) -> ValenceReport {
    match model {
        ValenceModel::RdkitLike => perceive_rdkit_like_valence(mol),
    }
}

fn perceive_rdkit_like_valence(mol: &mut Molecule) -> ValenceReport {
    let mut assignments = Vec::<(AtomId, u8)>::new();
    let mut issues = Vec::new();
    for (atom_id, atom) in mol.atoms() {
        let explicit = explicit_valence(mol, atom_id).saturating_add(atom.explicit_hydrogens);
        match allowed_valences(atom) {
            Some(allowed) => {
                if let Some(max_allowed) = allowed.iter().copied().max() {
                    if explicit > max_allowed {
                        issues.push(ValenceIssue::ValenceExceeded {
                            atom: atom_id,
                            explicit_valence: explicit,
                            max_allowed,
                        });
                        assignments.push((atom_id, 0));
                    } else {
                        let target = allowed
                            .iter()
                            .copied()
                            .find(|allowed| *allowed >= explicit)
                            .unwrap_or(explicit);
                        assignments.push((atom_id, target - explicit));
                    }
                }
            }
            None => issues.push(ValenceIssue::UnsupportedElement(atom_id)),
        }
    }
    for (atom_id, hydrogens) in assignments {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.implicit_hydrogens = Some(hydrogens);
        }
    }
    mol.perception.valence = ComputedState::Fresh;
    ValenceReport { issues }
}

fn explicit_valence(mol: &Molecule, atom: AtomId) -> u8 {
    mol.incident_bonds(atom)
        .ok()
        .into_iter()
        .flatten()
        .map(|(_, bond)| bond_order_valence(bond.order))
        .sum()
}

fn bond_order_valence(order: BondOrder) -> u8 {
    match order {
        BondOrder::Zero | BondOrder::Dative => 0,
        BondOrder::Single | BondOrder::Aromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
    }
}

fn allowed_valences(atom: &Atom) -> Option<&'static [u8]> {
    match (atom.element.symbol(), atom.formal_charge) {
        ("H", 0) => Some(&[1]),
        ("B", _) => Some(&[3]),
        ("C", 0) => Some(&[4]),
        ("C", 1 | -1) => Some(&[3]),
        ("N", 1) => Some(&[4]),
        ("N", -1) => Some(&[2]),
        ("N", 0) => Some(&[3, 5]),
        ("O", 0) => Some(&[2]),
        ("O", -1 | 1) => Some(&[1]),
        ("F" | "Cl" | "Br" | "I", 0) => Some(&[1]),
        ("P", 0) => Some(&[3, 5]),
        ("S", 0) => Some(&[2, 4, 6]),
        ("S", -1 | 1) => Some(&[1, 3, 5]),
        _ => None,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SanitizeOptions {
    pub perceive_valence: bool,
    pub perceive_rings: bool,
    pub perceive_aromaticity: bool,
}

impl Default for SanitizeOptions {
    fn default() -> Self {
        Self {
            perceive_valence: true,
            perceive_rings: true,
            perceive_aromaticity: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SanitizeReport {
    pub valence: Option<ValenceReport>,
    pub ring_count: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SanitizeError {
    Valence(ValenceReport),
    Aromaticity(AromaticityError),
}

impl fmt::Display for SanitizeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Valence(report) => write!(
                f,
                "valence perception reported {} issue(s)",
                report.issues.len()
            ),
            Self::Aromaticity(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for SanitizeError {}

pub fn sanitize_small_molecule(
    molecule: &mut SmallMolecule,
    options: SanitizeOptions,
) -> std::result::Result<SanitizeReport, SanitizeError> {
    let valence = if options.perceive_valence {
        let report = perceive_valence(&mut molecule.mol, ValenceModel::RdkitLike);
        if !report.is_ok() {
            return Err(SanitizeError::Valence(report));
        }
        Some(report)
    } else {
        None
    };
    let ring_count = if options.perceive_rings {
        Some(perceive_ring_set(&mut molecule.mol).len())
    } else {
        None
    };
    if options.perceive_aromaticity {
        perceive_aromaticity(&mut molecule.mol, AromaticityModel::RdkitLike)
            .map_err(SanitizeError::Aromaticity)?;
    }
    Ok(SanitizeReport {
        valence,
        ring_count,
    })
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SmallMolecule {
    pub mol: Molecule,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MolParseOptions;

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

pub fn read_mol_v2000_str(input: &str) -> std::result::Result<SmallMolecule, SdfParseError> {
    let normalized = input.replace("\r\n", "\n").replace('\r', "\n");
    let lines = normalized.lines().collect::<Vec<_>>();
    parse_mol_v2000_lines(1, 1, &lines)
}

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
    let title = lines.first().copied().unwrap_or_default().to_owned();
    let end_index = lines
        .iter()
        .position(|line| line.trim() == "M  END")
        .ok_or_else(|| SdfParseError::new(record, start_line, "missing M  END"))?;
    let mut molecule = parse_mol_v2000_lines(record, start_line, &lines[..=end_index])?;
    parse_sdf_data_fields(
        record,
        start_line + end_index + 1,
        &mut molecule.mol,
        &lines[end_index + 1..],
    )?;
    Ok(SdfRecord { title, molecule })
}

fn parse_mol_v2000_lines(
    record: usize,
    start_line: usize,
    lines: &[&str],
) -> std::result::Result<SmallMolecule, SdfParseError> {
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
    let mut conformer = Conformer::with_atom_capacity(atom_count);
    for atom_index in 0..atom_count {
        let line_number = start_line + atom_start + atom_index;
        let atom_line = lines[atom_start + atom_index];
        let symbol = atom_symbol_from_v2000_line(atom_line)
            .ok_or_else(|| SdfParseError::new(record, line_number, "invalid atom line"))?;
        let element = Element::from_symbol(symbol).ok_or_else(|| {
            SdfParseError::new(
                record,
                line_number,
                format!("unknown element symbol `{symbol}`"),
            )
        })?;
        let mut atom = Atom::new(element);
        apply_atom_v2000_fields(&mut atom, atom_line);
        let atom_id = mol.add_atom(atom);
        if let Some(point) = atom_coordinates_from_v2000_line(atom_line) {
            conformer.set_position(atom_id, point);
        }
        atom_ids.push(atom_id);
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
    parse_m_records(
        record,
        start_line + property_start,
        &mut mol,
        &atom_ids,
        &lines[property_start..end_index],
    )?;
    if conformer.positions().next().is_some() {
        mol.add_conformer(conformer);
    }

    Ok(SmallMolecule { mol })
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

fn atom_coordinates_from_v2000_line(line: &str) -> Option<Point3> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    Some(Point3::new(
        fields.first()?.parse().ok()?,
        fields.get(1)?.parse().ok()?,
        fields.get(2)?.parse().ok()?,
    ))
}

fn apply_atom_v2000_fields(atom: &mut Atom, line: &str) {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    if let Some(charge_code) = fields.get(5).and_then(|value| value.parse::<i8>().ok()) {
        atom.formal_charge = match charge_code {
            1 => 3,
            2 => 2,
            3 => 1,
            5 => -1,
            6 => -2,
            7 => -3,
            _ => 0,
        };
    }
    if let Some(atom_map) = fields
        .get(13)
        .or_else(|| fields.get(12))
        .and_then(|value| value.parse::<u32>().ok())
    {
        if atom_map != 0 {
            atom.atom_map = Some(atom_map);
        }
    }
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

fn parse_m_records(
    record: usize,
    start_line: usize,
    mol: &mut Molecule,
    atom_ids: &[AtomId],
    lines: &[&str],
) -> std::result::Result<(), SdfParseError> {
    for (offset, line) in lines.iter().enumerate() {
        let fields = line.split_whitespace().collect::<Vec<_>>();
        match fields.as_slice() {
            ["M", "CHG", count, rest @ ..] => {
                parse_atom_value_pairs(
                    record,
                    start_line + offset,
                    count,
                    rest,
                    mol,
                    atom_ids,
                    |atom, value| {
                        atom.formal_charge = value as i8;
                    },
                )?;
            }
            ["M", "ISO", count, rest @ ..] => {
                parse_atom_value_pairs(
                    record,
                    start_line + offset,
                    count,
                    rest,
                    mol,
                    atom_ids,
                    |atom, value| {
                        atom.isotope = (value > 0).then_some(value as u16);
                    },
                )?;
            }
            ["M", "RAD", count, rest @ ..] => {
                parse_atom_value_pairs(
                    record,
                    start_line + offset,
                    count,
                    rest,
                    mol,
                    atom_ids,
                    |atom, value| {
                        atom.radical_electrons = match value {
                            1 => 2,
                            2 => 1,
                            3 => 2,
                            _ => 0,
                        };
                    },
                )?;
            }
            _ => {}
        }
    }
    Ok(())
}

fn parse_atom_value_pairs<F>(
    record: usize,
    line: usize,
    count: &str,
    rest: &[&str],
    mol: &mut Molecule,
    atom_ids: &[AtomId],
    mut apply: F,
) -> std::result::Result<(), SdfParseError>
where
    F: FnMut(&mut Atom, i32),
{
    let count = count
        .parse::<usize>()
        .map_err(|_| SdfParseError::new(record, line, "invalid M record count"))?;
    if rest.len() < count * 2 {
        return Err(SdfParseError::new(record, line, "truncated M record"));
    }
    for pair in rest.chunks(2).take(count) {
        let atom_index = pair[0]
            .parse::<usize>()
            .map_err(|_| SdfParseError::new(record, line, "invalid M record atom index"))?;
        let value = pair[1]
            .parse::<i32>()
            .map_err(|_| SdfParseError::new(record, line, "invalid M record value"))?;
        let atom_id = atom_ids
            .get(atom_index.saturating_sub(1))
            .copied()
            .ok_or_else(|| SdfParseError::new(record, line, "M record atom outside atom block"))?;
        let atom = mol
            .atom_mut(atom_id)
            .map_err(|error| SdfParseError::new(record, line, error.to_string()))?;
        apply(atom, value);
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MolWriteError {
    pub message: String,
}

impl MolWriteError {
    fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl fmt::Display for MolWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for MolWriteError {}

pub fn write_mol_v2000(molecule: &SmallMolecule) -> std::result::Result<String, MolWriteError> {
    let mol = &molecule.mol;
    if mol.atom_count() > 999 || mol.bond_count() > 999 {
        return Err(MolWriteError::new(
            "V2000 writer supports at most 999 atoms and 999 bonds",
        ));
    }
    let atoms = mol.atom_ids().collect::<Vec<_>>();
    let bonds = mol.bond_ids().collect::<Vec<_>>();
    let mut atom_index = BTreeMap::new();
    for (index, atom_id) in atoms.iter().enumerate() {
        atom_index.insert(*atom_id, index + 1);
    }

    let title = prop_string(mol, "sdf.title").unwrap_or_default();
    let program = prop_string(mol, "sdf.program").unwrap_or_else(|| "molecules".to_owned());
    let comment = prop_string(mol, "sdf.comment").unwrap_or_default();
    let conformer = mol.first_conformer().map(|(_, conformer)| conformer);
    let mut out = String::new();
    out.push_str(&format!("{title}\n{program}\n{comment}\n"));
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0            999 V2000\n",
        atoms.len(),
        bonds.len()
    ));

    for atom_id in &atoms {
        let atom = mol
            .atom(*atom_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        let point = conformer
            .and_then(|conformer| conformer.position(*atom_id))
            .unwrap_or_default();
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3}{:>2}{:>3}  0  0  0  0  0  0  0{:>3}  0  0\n",
            point.x,
            point.y,
            point.z,
            atom.element.symbol(),
            0,
            v2000_charge_code(atom.formal_charge),
            atom.atom_map.unwrap_or(0)
        ));
    }

    for bond_id in &bonds {
        let bond = mol
            .bond(*bond_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        let a = atom_index
            .get(&bond.a())
            .ok_or_else(|| MolWriteError::new("bond endpoint missing from atom table"))?;
        let b = atom_index
            .get(&bond.b())
            .ok_or_else(|| MolWriteError::new("bond endpoint missing from atom table"))?;
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}  0  0  0  0\n",
            a,
            b,
            v2000_bond_code(bond.order)?
        ));
    }

    push_m_record(
        &mut out,
        "CHG",
        atoms
            .iter()
            .filter_map(|id| {
                let atom = mol.atom(*id).ok()?;
                (atom.formal_charge != 0)
                    .then_some((*atom_index.get(id)? as i32, atom.formal_charge as i32))
            })
            .collect(),
    );
    push_m_record(
        &mut out,
        "ISO",
        atoms
            .iter()
            .filter_map(|id| {
                let atom = mol.atom(*id).ok()?;
                atom.isotope.map(|isotope| {
                    (
                        *atom_index.get(id).expect("atom indexed") as i32,
                        isotope as i32,
                    )
                })
            })
            .collect(),
    );
    push_m_record(
        &mut out,
        "RAD",
        atoms
            .iter()
            .filter_map(|id| {
                let atom = mol.atom(*id).ok()?;
                (atom.radical_electrons != 0)
                    .then_some((*atom_index.get(id)? as i32, atom.radical_electrons as i32))
            })
            .collect(),
    );
    out.push_str("M  END\n");
    Ok(out)
}

pub fn write_sdf_v2000(molecules: &[SmallMolecule]) -> std::result::Result<String, MolWriteError> {
    let mut out = String::new();
    for molecule in molecules {
        out.push_str(&write_mol_v2000(molecule)?);
        for (key, value) in molecule.mol.props() {
            if let (Some(name), PropValue::String(text)) = (key.strip_prefix("sdf.field."), value) {
                out.push_str(&format!(">  <{name}>\n{text}\n\n"));
            }
        }
        out.push_str("$$$$\n");
    }
    Ok(out)
}

fn prop_string(mol: &Molecule, key: &str) -> Option<String> {
    match mol.props().get(key) {
        Some(PropValue::String(value)) => Some(value.clone()),
        _ => None,
    }
}

fn v2000_charge_code(charge: i8) -> i8 {
    match charge {
        3 => 1,
        2 => 2,
        1 => 3,
        -1 => 5,
        -2 => 6,
        -3 => 7,
        _ => 0,
    }
}

fn v2000_bond_code(order: BondOrder) -> std::result::Result<u8, MolWriteError> {
    match order {
        BondOrder::Zero => Ok(0),
        BondOrder::Single => Ok(1),
        BondOrder::Double => Ok(2),
        BondOrder::Triple => Ok(3),
        BondOrder::Aromatic => Ok(4),
        BondOrder::Dative => Ok(9),
        BondOrder::Quadruple => Err(MolWriteError::new(
            "V2000 writer does not support quadruple bonds",
        )),
    }
}

fn push_m_record(out: &mut String, code: &str, pairs: Vec<(i32, i32)>) {
    for chunk in pairs.chunks(8) {
        if chunk.is_empty() {
            continue;
        }
        out.push_str(&format!("M  {code} {:>2}", chunk.len()));
        for (atom, value) in chunk {
            out.push_str(&format!("{atom:>4}{value:>4}"));
        }
        out.push('\n');
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SmilesParseOptions;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SmilesWriteOptions;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesParseError {
    pub offset: usize,
    pub message: String,
}

impl SmilesParseError {
    fn new(offset: usize, message: impl Into<String>) -> Self {
        Self {
            offset,
            message: message.into(),
        }
    }
}

impl fmt::Display for SmilesParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SMILES parse error at {}: {}", self.offset, self.message)
    }
}

impl std::error::Error for SmilesParseError {}

pub fn read_smiles_str(
    input: &str,
    _options: SmilesParseOptions,
) -> std::result::Result<SmallMolecule, SmilesParseError> {
    let chars = input.char_indices().collect::<Vec<_>>();
    let mut mol = Molecule::new();
    let mut current: Option<AtomId> = None;
    let mut stack = Vec::<AtomId>::new();
    let mut pending_bond = BondOrder::Single;
    let mut rings = BTreeMap::<char, (AtomId, BondOrder)>::new();
    let mut cursor = 0;
    while cursor < chars.len() {
        let (offset, ch) = chars[cursor];
        match ch {
            '(' => {
                let atom =
                    current.ok_or_else(|| SmilesParseError::new(offset, "branch without atom"))?;
                stack.push(atom);
                cursor += 1;
            }
            ')' => {
                current = Some(
                    stack
                        .pop()
                        .ok_or_else(|| SmilesParseError::new(offset, "unmatched branch close"))?,
                );
                cursor += 1;
            }
            '.' => {
                current = None;
                pending_bond = BondOrder::Single;
                cursor += 1;
            }
            '-' => {
                pending_bond = BondOrder::Single;
                cursor += 1;
            }
            '=' => {
                pending_bond = BondOrder::Double;
                cursor += 1;
            }
            '#' => {
                pending_bond = BondOrder::Triple;
                cursor += 1;
            }
            ':' => {
                pending_bond = BondOrder::Aromatic;
                cursor += 1;
            }
            '0'..='9' => {
                let atom = current
                    .ok_or_else(|| SmilesParseError::new(offset, "ring closure without atom"))?;
                if let Some((other, order)) = rings.remove(&ch) {
                    let order = if pending_bond == BondOrder::Single {
                        order
                    } else {
                        pending_bond
                    };
                    mol.add_bond(other, atom, order)
                        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
                    pending_bond = BondOrder::Single;
                } else {
                    rings.insert(ch, (atom, pending_bond));
                    pending_bond = BondOrder::Single;
                }
                cursor += 1;
            }
            '[' => {
                let (atom, next_cursor) = parse_bracket_atom(&chars, cursor)?;
                let atom_id = mol.add_atom(atom);
                if let Some(previous) = current {
                    mol.add_bond(previous, atom_id, pending_bond)
                        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
                }
                current = Some(atom_id);
                pending_bond = BondOrder::Single;
                cursor = next_cursor;
            }
            _ => {
                let (atom, next_cursor) = parse_organic_atom(&chars, cursor)?;
                let atom_id = mol.add_atom(atom);
                if let Some(previous) = current {
                    mol.add_bond(previous, atom_id, pending_bond)
                        .map_err(|error| SmilesParseError::new(offset, error.to_string()))?;
                }
                current = Some(atom_id);
                pending_bond = BondOrder::Single;
                cursor = next_cursor;
            }
        }
    }
    if !stack.is_empty() {
        return Err(SmilesParseError::new(input.len(), "unclosed branch"));
    }
    if !rings.is_empty() {
        return Err(SmilesParseError::new(input.len(), "unclosed ring closure"));
    }
    Ok(SmallMolecule { mol })
}

fn parse_organic_atom(
    chars: &[(usize, char)],
    cursor: usize,
) -> std::result::Result<(Atom, usize), SmilesParseError> {
    let (offset, ch) = chars[cursor];
    let mut symbol = ch.to_string();
    let mut aromatic = false;
    let mut next = cursor + 1;
    if matches!(ch, 'B' | 'C') && chars.get(cursor + 1).map(|(_, c)| *c) == Some('l')
        || ch == 'B' && chars.get(cursor + 1).map(|(_, c)| *c) == Some('r')
    {
        symbol.push(chars[cursor + 1].1);
        next += 1;
    } else if matches!(ch, 'b' | 'c' | 'n' | 'o' | 'p' | 's') {
        symbol = ch.to_ascii_uppercase().to_string();
        aromatic = true;
    }
    let element = Element::from_symbol(&symbol)
        .ok_or_else(|| SmilesParseError::new(offset, format!("unsupported atom `{ch}`")))?;
    let mut atom = Atom::new(element);
    atom.aromatic = aromatic;
    Ok((atom, next))
}

fn parse_bracket_atom(
    chars: &[(usize, char)],
    cursor: usize,
) -> std::result::Result<(Atom, usize), SmilesParseError> {
    let start = chars[cursor].0;
    let mut end = cursor + 1;
    while end < chars.len() && chars[end].1 != ']' {
        end += 1;
    }
    if end == chars.len() {
        return Err(SmilesParseError::new(start, "unclosed bracket atom"));
    }
    let text = chars[cursor + 1..end]
        .iter()
        .map(|(_, c)| *c)
        .collect::<String>();
    let mut index = 0;
    let isotope_digits = take_digits(&text, &mut index);
    let symbol = take_element_symbol(&text, &mut index)
        .ok_or_else(|| SmilesParseError::new(start, "bracket atom missing element"))?;
    let mut aromatic = false;
    let canonical_symbol = if symbol.chars().next().is_some_and(char::is_lowercase) {
        aromatic = true;
        let mut chars = symbol.chars();
        let first = chars.next().expect("symbol has first").to_ascii_uppercase();
        format!("{first}{}", chars.as_str())
    } else {
        symbol
    };
    let element = Element::from_symbol(&canonical_symbol)
        .ok_or_else(|| SmilesParseError::new(start, "unsupported bracket element"))?;
    let mut atom = Atom::new(element);
    atom.aromatic = aromatic;
    if !isotope_digits.is_empty() {
        atom.isotope = isotope_digits.parse::<u16>().ok();
    }
    while index < text.len() {
        let rest = &text[index..];
        if let Some(after_h) = rest.strip_prefix('H') {
            index += 1;
            let mut h_index = 0;
            let digits = take_digits(after_h, &mut h_index);
            atom.explicit_hydrogens = if digits.is_empty() {
                1
            } else {
                digits.parse().unwrap_or(0)
            };
            index += h_index;
        } else if rest.starts_with('+') || rest.starts_with('-') {
            let sign = if rest.starts_with('+') { 1 } else { -1 };
            index += 1;
            let mut repeats = 1i8;
            while text[index..].starts_with(if sign > 0 { '+' } else { '-' }) {
                repeats += 1;
                index += 1;
            }
            let mut charge_index = 0;
            let digits = take_digits(&text[index..], &mut charge_index);
            if !digits.is_empty() {
                repeats = digits.parse().unwrap_or(repeats);
                index += charge_index;
            }
            atom.formal_charge = sign * repeats;
        } else if let Some(after_colon) = rest.strip_prefix(':') {
            let mut map_index = 0;
            let digits = take_digits(after_colon, &mut map_index);
            atom.atom_map = digits.parse::<u32>().ok();
            index += 1 + map_index;
        } else {
            index += 1;
        }
    }
    Ok((atom, end + 1))
}

fn take_digits(text: &str, index: &mut usize) -> String {
    let start = *index;
    while *index < text.len() && text.as_bytes()[*index].is_ascii_digit() {
        *index += 1;
    }
    text[start..*index].to_owned()
}

fn take_element_symbol(text: &str, index: &mut usize) -> Option<String> {
    let bytes = text.as_bytes();
    if *index >= bytes.len() || !bytes[*index].is_ascii_alphabetic() {
        return None;
    }
    let mut symbol = String::new();
    symbol.push(bytes[*index] as char);
    *index += 1;
    if *index < bytes.len() && bytes[*index].is_ascii_lowercase() {
        symbol.push(bytes[*index] as char);
        *index += 1;
    }
    Some(symbol)
}

pub fn write_smiles(
    molecule: &SmallMolecule,
    _options: SmilesWriteOptions,
) -> std::result::Result<String, MolWriteError> {
    let mol = &molecule.mol;
    let mut visited = BTreeMap::<AtomId, ()>::new();
    let mut ring_numbers = BTreeMap::<(AtomId, AtomId), usize>::new();
    let mut next_ring = 1usize;
    let mut parts = Vec::new();
    for start in mol.atom_ids() {
        if visited.contains_key(&start) {
            continue;
        }
        parts.push(write_smiles_component(
            mol,
            start,
            None,
            &mut visited,
            &mut ring_numbers,
            &mut next_ring,
        )?);
    }
    Ok(parts.join("."))
}

fn write_smiles_component(
    mol: &Molecule,
    atom_id: AtomId,
    parent: Option<AtomId>,
    visited: &mut BTreeMap<AtomId, ()>,
    ring_numbers: &mut BTreeMap<(AtomId, AtomId), usize>,
    next_ring: &mut usize,
) -> std::result::Result<String, MolWriteError> {
    visited.insert(atom_id, ());
    let atom = mol
        .atom(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?;
    let mut out = smiles_atom(atom);
    let mut children = Vec::<(BondOrder, AtomId)>::new();
    for (_, bond) in mol
        .incident_bonds(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
    {
        let neighbor = bond.other_atom(atom_id);
        if Some(neighbor) == parent {
            continue;
        }
        if visited.contains_key(&neighbor) {
            let key = ordered_atom_pair(atom_id, neighbor);
            if let std::collections::btree_map::Entry::Vacant(entry) = ring_numbers.entry(key) {
                let number = *next_ring;
                *next_ring += 1;
                entry.insert(number);
                out.push_str(smiles_bond(bond.order));
                out.push_str(&number.to_string());
            }
        } else {
            children.push((bond.order, neighbor));
        }
    }
    children.sort_by_key(|(_, atom)| *atom);
    for (index, (order, child)) in children.into_iter().enumerate() {
        let child_text = format!(
            "{}{}",
            smiles_bond(order),
            write_smiles_component(mol, child, Some(atom_id), visited, ring_numbers, next_ring)?
        );
        if index == 0 {
            out.push_str(&child_text);
        } else {
            out.push('(');
            out.push_str(&child_text);
            out.push(')');
        }
    }
    Ok(out)
}

fn smiles_bond(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Aromatic => ":",
        BondOrder::Zero | BondOrder::Dative | BondOrder::Quadruple => "-",
    }
}

fn smiles_atom(atom: &Atom) -> String {
    let organic = atom.isotope.is_none()
        && atom.formal_charge == 0
        && atom.explicit_hydrogens == 0
        && atom.atom_map.is_none()
        && matches!(
            atom.element.symbol(),
            "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        );
    if organic {
        if atom.aromatic {
            atom.element.symbol().to_ascii_lowercase()
        } else {
            atom.element.symbol().to_owned()
        }
    } else {
        let mut out = String::from("[");
        if let Some(isotope) = atom.isotope {
            out.push_str(&isotope.to_string());
        }
        if atom.aromatic {
            out.push_str(&atom.element.symbol().to_ascii_lowercase());
        } else {
            out.push_str(atom.element.symbol());
        }
        if atom.explicit_hydrogens > 0 {
            out.push('H');
            if atom.explicit_hydrogens > 1 {
                out.push_str(&atom.explicit_hydrogens.to_string());
            }
        }
        if atom.formal_charge > 0 {
            out.push('+');
            if atom.formal_charge > 1 {
                out.push_str(&atom.formal_charge.to_string());
            }
        } else if atom.formal_charge < 0 {
            out.push('-');
            if atom.formal_charge < -1 {
                out.push_str(&(-atom.formal_charge).to_string());
            }
        }
        if let Some(map) = atom.atom_map {
            out.push(':');
            out.push_str(&map.to_string());
        }
        out.push(']');
        out
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MacroMolecule {
    pub mol: Molecule,
    pub hierarchy: BioHierarchy,
}

impl MacroMolecule {
    pub fn add_atom_site(
        &mut self,
        residue: ResidueId,
        atom: AtomId,
        metadata: AtomSiteMetadata,
    ) -> std::result::Result<AtomSiteId, BioHierarchyError> {
        self.mol
            .atom(atom)
            .map_err(|_| BioHierarchyError::InvalidAtomId(atom))?;
        self.hierarchy.add_atom_site(residue, atom, metadata)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ModelId(u32);

impl ModelId {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ChainId(u32);

impl ChainId {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ResidueId(u32);

impl ResidueId {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomSiteId(u32);

impl AtomSiteId {
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

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BioHierarchy {
    models: Vec<Model>,
    chains: Vec<Chain>,
    residues: Vec<Residue>,
    atom_sites: Vec<AtomSite>,
    atom_lookup: BTreeMap<AtomId, AtomSiteId>,
    pub props: PropMap,
}

impl BioHierarchy {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_model(&mut self, model_id: impl Into<String>) -> ModelId {
        let id = ModelId::new(self.models.len() as u32);
        self.models.push(Model {
            id,
            model_id: model_id.into(),
            chains: Vec::new(),
            props: PropMap::new(),
        });
        id
    }

    pub fn add_chain(
        &mut self,
        model: ModelId,
        label_id: impl Into<String>,
        author_id: Option<String>,
    ) -> std::result::Result<ChainId, BioHierarchyError> {
        self.model(model)?;
        let id = ChainId::new(self.chains.len() as u32);
        self.chains.push(Chain {
            id,
            model,
            label_id: label_id.into(),
            author_id,
            residues: Vec::new(),
            props: PropMap::new(),
        });
        self.models[model.index()].chains.push(id);
        Ok(id)
    }

    pub fn add_residue(
        &mut self,
        chain: ChainId,
        name: impl Into<String>,
        label_seq_id: Option<i32>,
        author_seq_id: Option<String>,
        insertion_code: Option<String>,
    ) -> std::result::Result<ResidueId, BioHierarchyError> {
        self.chain(chain)?;
        let id = ResidueId::new(self.residues.len() as u32);
        self.residues.push(Residue {
            id,
            chain,
            name: name.into(),
            label_seq_id,
            author_seq_id,
            insertion_code,
            atom_sites: Vec::new(),
            props: PropMap::new(),
        });
        self.chains[chain.index()].residues.push(id);
        Ok(id)
    }

    pub fn add_atom_site(
        &mut self,
        residue: ResidueId,
        atom: AtomId,
        metadata: AtomSiteMetadata,
    ) -> std::result::Result<AtomSiteId, BioHierarchyError> {
        self.residue(residue)?;
        if self.atom_lookup.contains_key(&atom) {
            return Err(BioHierarchyError::DuplicateAtomPlacement(atom));
        }
        let id = AtomSiteId::new(self.atom_sites.len() as u32);
        self.atom_sites.push(AtomSite {
            id,
            residue,
            atom,
            metadata,
            props: PropMap::new(),
        });
        self.residues[residue.index()].atom_sites.push(id);
        self.atom_lookup.insert(atom, id);
        Ok(id)
    }

    pub fn model(&self, id: ModelId) -> std::result::Result<&Model, BioHierarchyError> {
        self.models
            .get(id.index())
            .ok_or(BioHierarchyError::InvalidModelId(id))
    }

    pub fn chain(&self, id: ChainId) -> std::result::Result<&Chain, BioHierarchyError> {
        self.chains
            .get(id.index())
            .ok_or(BioHierarchyError::InvalidChainId(id))
    }

    pub fn residue(&self, id: ResidueId) -> std::result::Result<&Residue, BioHierarchyError> {
        self.residues
            .get(id.index())
            .ok_or(BioHierarchyError::InvalidResidueId(id))
    }

    pub fn atom_site(&self, id: AtomSiteId) -> std::result::Result<&AtomSite, BioHierarchyError> {
        self.atom_sites
            .get(id.index())
            .ok_or(BioHierarchyError::InvalidAtomSiteId(id))
    }

    pub fn atom_site_for_atom(&self, atom: AtomId) -> Option<&AtomSite> {
        self.atom_lookup
            .get(&atom)
            .and_then(|id| self.atom_sites.get(id.index()))
    }

    pub fn models(&self) -> impl Iterator<Item = (ModelId, &Model)> {
        self.models.iter().map(|model| (model.id, model))
    }

    pub fn chains(&self) -> impl Iterator<Item = (ChainId, &Chain)> {
        self.chains.iter().map(|chain| (chain.id, chain))
    }

    pub fn residues(&self) -> impl Iterator<Item = (ResidueId, &Residue)> {
        self.residues.iter().map(|residue| (residue.id, residue))
    }

    pub fn atom_sites(&self) -> impl Iterator<Item = (AtomSiteId, &AtomSite)> {
        self.atom_sites.iter().map(|site| (site.id, site))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    pub id: ModelId,
    pub model_id: String,
    pub chains: Vec<ChainId>,
    pub props: PropMap,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Chain {
    pub id: ChainId,
    pub model: ModelId,
    pub label_id: String,
    pub author_id: Option<String>,
    pub residues: Vec<ResidueId>,
    pub props: PropMap,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    pub id: ResidueId,
    pub chain: ChainId,
    pub name: String,
    pub label_seq_id: Option<i32>,
    pub author_seq_id: Option<String>,
    pub insertion_code: Option<String>,
    pub atom_sites: Vec<AtomSiteId>,
    pub props: PropMap,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomSite {
    pub id: AtomSiteId,
    pub residue: ResidueId,
    pub atom: AtomId,
    pub metadata: AtomSiteMetadata,
    pub props: PropMap,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct AtomSiteMetadata {
    pub label_atom_id: Option<String>,
    pub auth_atom_id: Option<String>,
    pub label_alt_id: Option<String>,
    pub occupancy: Option<f64>,
    pub b_factor: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BioHierarchyError {
    InvalidModelId(ModelId),
    InvalidChainId(ChainId),
    InvalidResidueId(ResidueId),
    InvalidAtomSiteId(AtomSiteId),
    InvalidAtomId(AtomId),
    DuplicateAtomPlacement(AtomId),
}

impl fmt::Display for BioHierarchyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidModelId(id) => write!(f, "invalid model id: {}", id.raw()),
            Self::InvalidChainId(id) => write!(f, "invalid chain id: {}", id.raw()),
            Self::InvalidResidueId(id) => write!(f, "invalid residue id: {}", id.raw()),
            Self::InvalidAtomSiteId(id) => write!(f, "invalid atom-site id: {}", id.raw()),
            Self::InvalidAtomId(id) => write!(f, "invalid hierarchy atom id: {id}"),
            Self::DuplicateAtomPlacement(id) => write!(f, "duplicate hierarchy placement for {id}"),
        }
    }
}

impl std::error::Error for BioHierarchyError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MmcifParseOptions {
    pub strict: bool,
}

impl Default for MmcifParseOptions {
    fn default() -> Self {
        Self { strict: true }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifParseError {
    pub line: usize,
    pub message: String,
}

impl MmcifParseError {
    fn new(line: usize, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }
}

impl fmt::Display for MmcifParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "mmCIF parse error at line {}: {}",
            self.line, self.message
        )
    }
}

impl std::error::Error for MmcifParseError {}

pub fn read_mmcif_str(
    input: &str,
    options: MmcifParseOptions,
) -> std::result::Result<MacroMolecule, MmcifParseError> {
    let tokens = tokenize_mmcif(input)?;
    let atom_site_loop = find_atom_site_loop(&tokens)
        .ok_or_else(|| MmcifParseError::new(1, "missing _atom_site loop"))?;
    build_macro_molecule_from_atom_site_loop(atom_site_loop, options)
}

#[derive(Debug, Clone)]
struct MmcifToken {
    text: String,
    line: usize,
}

#[derive(Debug, Clone)]
struct MmcifLoop<'a> {
    tags: Vec<&'a MmcifToken>,
    values: Vec<&'a MmcifToken>,
}

fn tokenize_mmcif(input: &str) -> std::result::Result<Vec<MmcifToken>, MmcifParseError> {
    let normalized = input.replace("\r\n", "\n").replace('\r', "\n");
    let lines = normalized.lines().collect::<Vec<_>>();
    let mut tokens = Vec::new();
    let mut line_index = 0usize;

    while line_index < lines.len() {
        let line_number = line_index + 1;
        let line = lines[line_index];
        if line.starts_with(';') {
            let mut text = String::new();
            line_index += 1;
            while line_index < lines.len() && !lines[line_index].starts_with(';') {
                if !text.is_empty() {
                    text.push('\n');
                }
                text.push_str(lines[line_index]);
                line_index += 1;
            }
            if line_index == lines.len() {
                return Err(MmcifParseError::new(
                    line_number,
                    "unterminated semicolon text",
                ));
            }
            tokens.push(MmcifToken {
                text,
                line: line_number,
            });
            line_index += 1;
            continue;
        }

        let bytes = line.as_bytes();
        let mut column = 0usize;
        while column < bytes.len() {
            while column < bytes.len() && bytes[column].is_ascii_whitespace() {
                column += 1;
            }
            if column == bytes.len() || bytes[column] == b'#' {
                break;
            }
            let start = column;
            let text = if bytes[column] == b'\'' || bytes[column] == b'"' {
                let quote = bytes[column];
                column += 1;
                let value_start = column;
                while column < bytes.len() && bytes[column] != quote {
                    column += 1;
                }
                if column == bytes.len() {
                    return Err(MmcifParseError::new(
                        line_number,
                        "unterminated quoted value",
                    ));
                }
                let value = &line[value_start..column];
                column += 1;
                value.to_owned()
            } else {
                while column < bytes.len()
                    && !bytes[column].is_ascii_whitespace()
                    && bytes[column] != b'#'
                {
                    column += 1;
                }
                line[start..column].to_owned()
            };
            tokens.push(MmcifToken {
                text,
                line: line_number,
            });
        }
        line_index += 1;
    }

    Ok(tokens)
}

fn find_atom_site_loop(tokens: &[MmcifToken]) -> Option<MmcifLoop<'_>> {
    let mut index = 0usize;
    while index < tokens.len() {
        if tokens[index].text != "loop_" {
            index += 1;
            continue;
        }
        index += 1;
        let tag_start = index;
        while index < tokens.len() && tokens[index].text.starts_with('_') {
            index += 1;
        }
        let tags = tokens[tag_start..index].iter().collect::<Vec<_>>();
        let value_start = index;
        while index < tokens.len()
            && tokens[index].text != "loop_"
            && !tokens[index].text.starts_with("data_")
            && !tokens[index].text.starts_with('_')
        {
            index += 1;
        }
        if tags.iter().any(|tag| tag.text.starts_with("_atom_site.")) {
            return Some(MmcifLoop {
                tags,
                values: tokens[value_start..index].iter().collect(),
            });
        }
    }
    None
}

fn build_macro_molecule_from_atom_site_loop(
    atom_loop: MmcifLoop<'_>,
    options: MmcifParseOptions,
) -> std::result::Result<MacroMolecule, MmcifParseError> {
    let width = atom_loop.tags.len();
    if width == 0 || atom_loop.values.len() % width != 0 {
        let line = atom_loop
            .values
            .first()
            .map(|token| token.line)
            .unwrap_or(1);
        return Err(MmcifParseError::new(line, "atom-site loop has ragged rows"));
    }

    let tag_index = atom_loop
        .tags
        .iter()
        .enumerate()
        .map(|(index, token)| (token.text.as_str(), index))
        .collect::<BTreeMap<_, _>>();
    let mut macro_mol = MacroMolecule::default();
    let mut models = BTreeMap::<String, ModelId>::new();
    let mut chains = BTreeMap::<(String, String), ChainId>::new();
    let mut residues = BTreeMap::<(ChainId, String, Option<i32>, Option<String>), ResidueId>::new();

    for row in atom_loop.values.chunks(width) {
        let line = row.first().map(|token| token.line).unwrap_or(1);
        let type_symbol = required_mmcif_value(row, &tag_index, "_atom_site.type_symbol", line)?;
        let element = Element::from_symbol(type_symbol).ok_or_else(|| {
            MmcifParseError::new(line, format!("unknown atom-site element `{type_symbol}`"))
        })?;
        let label_atom_id = optional_mmcif_value(row, &tag_index, "_atom_site.label_atom_id");
        let auth_atom_id = optional_mmcif_value(row, &tag_index, "_atom_site.auth_atom_id");
        let label_asym_id = optional_mmcif_value(row, &tag_index, "_atom_site.label_asym_id");
        let auth_asym_id = optional_mmcif_value(row, &tag_index, "_atom_site.auth_asym_id");
        let label_comp_id = optional_mmcif_value(row, &tag_index, "_atom_site.label_comp_id");
        let auth_comp_id = optional_mmcif_value(row, &tag_index, "_atom_site.auth_comp_id");
        let label_seq_id =
            optional_i32_mmcif_value(row, &tag_index, "_atom_site.label_seq_id", line)?;
        let auth_seq_id = optional_mmcif_value(row, &tag_index, "_atom_site.auth_seq_id");
        let insertion_code = optional_mmcif_value(row, &tag_index, "_atom_site.pdbx_PDB_ins_code");
        let model_key = optional_mmcif_value(row, &tag_index, "_atom_site.pdbx_PDB_model_num")
            .unwrap_or("1")
            .to_owned();
        let label_alt_id = optional_mmcif_value(row, &tag_index, "_atom_site.label_alt_id");
        let occupancy = optional_f64_mmcif_value(row, &tag_index, "_atom_site.occupancy", line)?;
        let b_factor =
            optional_f64_mmcif_value(row, &tag_index, "_atom_site.B_iso_or_equiv", line)?;

        let chain_label = label_asym_id
            .or(auth_asym_id)
            .ok_or_else(|| MmcifParseError::new(line, "missing atom-site chain identifier"))?;
        let residue_name = label_comp_id
            .or(auth_comp_id)
            .ok_or_else(|| MmcifParseError::new(line, "missing atom-site residue name"))?;
        if options.strict && label_atom_id.is_none() {
            return Err(MmcifParseError::new(
                line,
                "missing atom-site label atom id",
            ));
        }

        let model = *models
            .entry(model_key.clone())
            .or_insert_with(|| macro_mol.hierarchy.add_model(model_key.clone()));
        let chain_key = (model_key.clone(), chain_label.to_owned());
        let chain = if let Some(chain) = chains.get(&chain_key) {
            *chain
        } else {
            let chain = macro_mol
                .hierarchy
                .add_chain(
                    model,
                    chain_label.to_owned(),
                    auth_asym_id.map(str::to_owned),
                )
                .map_err(|error| MmcifParseError::new(line, error.to_string()))?;
            chains.insert(chain_key, chain);
            chain
        };
        let residue_key = (
            chain,
            residue_name.to_owned(),
            label_seq_id,
            insertion_code.map(str::to_owned),
        );
        let residue = if let Some(residue) = residues.get(&residue_key) {
            *residue
        } else {
            let residue = macro_mol
                .hierarchy
                .add_residue(
                    chain,
                    residue_name.to_owned(),
                    label_seq_id,
                    auth_seq_id.map(str::to_owned),
                    insertion_code.map(str::to_owned),
                )
                .map_err(|error| MmcifParseError::new(line, error.to_string()))?;
            residues.insert(residue_key, residue);
            residue
        };

        let atom = macro_mol.mol.add_atom(Atom::new(element));
        macro_mol
            .add_atom_site(
                residue,
                atom,
                AtomSiteMetadata {
                    label_atom_id: label_atom_id.map(str::to_owned),
                    auth_atom_id: auth_atom_id.map(str::to_owned),
                    label_alt_id: label_alt_id.map(str::to_owned),
                    occupancy,
                    b_factor,
                },
            )
            .map_err(|error| MmcifParseError::new(line, error.to_string()))?;
    }

    Ok(macro_mol)
}

fn required_mmcif_value<'a>(
    row: &'a [&'a MmcifToken],
    tag_index: &BTreeMap<&str, usize>,
    tag: &str,
    line: usize,
) -> std::result::Result<&'a str, MmcifParseError> {
    optional_mmcif_value(row, tag_index, tag)
        .ok_or_else(|| MmcifParseError::new(line, format!("missing required {tag}")))
}

fn optional_mmcif_value<'a>(
    row: &'a [&'a MmcifToken],
    tag_index: &BTreeMap<&str, usize>,
    tag: &str,
) -> Option<&'a str> {
    let value = row.get(*tag_index.get(tag)?)?.text.as_str();
    (!matches!(value, "." | "?")).then_some(value)
}

fn optional_i32_mmcif_value(
    row: &[&MmcifToken],
    tag_index: &BTreeMap<&str, usize>,
    tag: &str,
    line: usize,
) -> std::result::Result<Option<i32>, MmcifParseError> {
    optional_mmcif_value(row, tag_index, tag)
        .map(|value| {
            value
                .parse()
                .map_err(|_| MmcifParseError::new(line, format!("invalid integer {tag}")))
        })
        .transpose()
}

fn optional_f64_mmcif_value(
    row: &[&MmcifToken],
    tag_index: &BTreeMap<&str, usize>,
    tag: &str,
    line: usize,
) -> std::result::Result<Option<f64>, MmcifParseError> {
    optional_mmcif_value(row, tag_index, tag)
        .map(|value| {
            value
                .parse()
                .map_err(|_| MmcifParseError::new(line, format!("invalid float {tag}")))
        })
        .transpose()
}

pub type Result<T> = std::result::Result<T, MoleculeError>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoleculeError {
    InvalidAtomId(AtomId),
    InvalidBondId(BondId),
    InvalidConformerId(ConformerId),
    SelfBond(AtomId),
    DuplicateBond { a: AtomId, b: AtomId },
    UnsupportedFeature(&'static str),
}

impl fmt::Display for MoleculeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(id) => write!(f, "invalid atom id: {id}"),
            Self::InvalidBondId(id) => write!(f, "invalid bond id: {id}"),
            Self::InvalidConformerId(id) => write!(f, "invalid conformer id: {id}"),
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
    fn mol_v2000_preserves_coordinates_charges_isotopes_radicals_and_atom_maps() {
        let input = "\
charged radical
molecules validation
metadata fixture
  2  1  0  0  0  0            999 V2000
    0.1000    0.2000    0.3000 N   0  0  0  0  0  0  0  0  0  7  0  0
    1.4000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  CHG  1   1   1
M  ISO  1   2  13
M  RAD  1   1   2
M  END
";

        let small = read_mol_v2000_str(input).expect("mol should parse");
        let atom0 = small.mol.atom(AtomId::new(0)).expect("atom exists");
        let atom1 = small.mol.atom(AtomId::new(1)).expect("atom exists");
        assert_eq!(atom0.formal_charge, 1);
        assert_eq!(atom0.radical_electrons, 1);
        assert_eq!(atom0.atom_map, Some(7));
        assert_eq!(atom1.isotope, Some(13));
        let (_, conformer) = small.mol.first_conformer().expect("conformer exists");
        assert_eq!(
            conformer.position(AtomId::new(0)),
            Some(Point3::new(0.1, 0.2, 0.3))
        );
    }

    #[test]
    fn mol_and_sdf_v2000_writers_round_trip_metadata_and_fields() {
        let input = "\
ammonium_acetate_like
molecules validation
M CHG and M ISO fixture
  4  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6000    0.7000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6000   -0.7000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  CHG  2   1   1   4  -1
M  ISO  1   2  13
M  END
>  <fixture_id>
charged_isotope_records

$$$$
";

        let records =
            read_sdf_v2000_records(input, SdfParseOptions::default()).expect("sdf should parse");
        let molecules = records
            .iter()
            .map(|record| record.molecule.clone())
            .collect::<Vec<_>>();
        let sdf = write_sdf_v2000(&molecules).expect("sdf should write");
        let reparsed =
            read_sdf_v2000_records(&sdf, SdfParseOptions::default()).expect("written sdf parses");

        assert_eq!(reparsed.len(), 1);
        assert_eq!(
            reparsed[0]
                .molecule
                .mol
                .atom(AtomId::new(0))
                .expect("atom")
                .formal_charge,
            1
        );
        assert_eq!(
            reparsed[0].molecule.mol.props().get("sdf.field.fixture_id"),
            Some(&PropValue::String("charged_isotope_records".to_owned()))
        );
    }

    #[test]
    fn valence_and_sanitization_are_explicit() {
        let mut small = read_smiles_str("CCO", SmilesParseOptions).expect("smiles should parse");
        assert_eq!(small.mol.perception().valence, ComputedState::Absent);

        let report = sanitize_small_molecule(&mut small, SanitizeOptions::default())
            .expect("ethanol should sanitize");

        assert!(report.valence.expect("valence report").is_ok());
        assert_eq!(small.mol.perception().valence, ComputedState::Fresh);
        assert_eq!(small.mol.perception().rings, ComputedState::Fresh);
        assert_eq!(
            small
                .mol
                .atom(AtomId::new(2))
                .expect("oxygen")
                .implicit_hydrogens,
            Some(1)
        );
    }

    #[test]
    fn valence_reports_excess_common_valence() {
        let mut mol = Molecule::new();
        let c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("C")));
        for _ in 0..5 {
            let h = mol.add_atom(Atom::new(Element::from_symbol("H").expect("H")));
            mol.add_bond(c, h, BondOrder::Single).expect("bond");
        }

        let report = perceive_valence(&mut mol, ValenceModel::RdkitLike);

        assert_eq!(report.issues.len(), 1);
        assert!(!report.is_ok());
    }

    #[test]
    fn ring_set_reports_a_basis_for_fused_rings() {
        let (mut mol, _, _) = ring_molecule(
            &["C", "C", "C", "C", "C", "C"],
            &[
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
            ],
        );
        let a = mol.add_atom(carbon());
        let b = mol.add_atom(carbon());
        mol.add_bond(AtomId::new(0), a, BondOrder::Single)
            .expect("bond");
        mol.add_bond(a, b, BondOrder::Single).expect("bond");
        mol.add_bond(b, AtomId::new(3), BondOrder::Single)
            .expect("bond");

        let ring_set = perceive_ring_set(&mut mol);

        assert_eq!(ring_set.len(), 2);
        assert!(ring_set.rings().iter().all(|ring| ring.atoms.len() >= 4));
    }

    #[test]
    fn smiles_parses_branches_rings_brackets_and_fragments_without_sanitizing() {
        let small = read_smiles_str("C(C)O.C1=CC=CC=C1.[13NH4+:7]", SmilesParseOptions)
            .expect("smiles should parse");

        assert_eq!(small.mol.atom_count(), 10);
        assert_eq!(small.mol.bond_count(), 8);
        assert_eq!(small.mol.perception().valence, ComputedState::Absent);
        let bracket_atom = small.mol.atom(AtomId::new(9)).expect("bracket atom");
        assert_eq!(bracket_atom.isotope, Some(13));
        assert_eq!(bracket_atom.explicit_hydrogens, 4);
        assert_eq!(bracket_atom.formal_charge, 1);
        assert_eq!(bracket_atom.atom_map, Some(7));
    }

    #[test]
    fn smiles_writer_round_trips_graph_shape() {
        let small = read_smiles_str("CC(=O)O", SmilesParseOptions).expect("smiles should parse");
        let text = write_smiles(&small, SmilesWriteOptions).expect("smiles should write");
        let reparsed =
            read_smiles_str(&text, SmilesParseOptions).expect("written smiles should parse");

        assert_eq!(reparsed.mol.atom_count(), small.mol.atom_count());
        assert_eq!(reparsed.mol.bond_count(), small.mol.bond_count());
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

        perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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
        perceive_aromaticity(&mut cyclohexane, AromaticityModel::RdkitLike)
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
        perceive_aromaticity(&mut cyclobutadiene, AromaticityModel::RdkitLike)
            .expect("cyclobutadiene should be supported");
        assert!(atoms
            .iter()
            .all(|atom| !cyclobutadiene.atom(*atom).expect("atom exists").aromatic));
        assert!(bonds
            .iter()
            .all(|bond| !cyclobutadiene.bond(*bond).expect("bond exists").aromatic));
    }

    #[test]
    fn aromaticity_supports_heteroaromatic_ring() {
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

        perceive_aromaticity(&mut furan_like, AromaticityModel::RdkitLike)
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

        perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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

        perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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
        perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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
    fn bio_hierarchy_adds_models_chains_residues_and_atom_sites() {
        let mut hierarchy = BioHierarchy::new();
        let model = hierarchy.add_model("1");
        let chain = hierarchy
            .add_chain(model, "A", Some("authA".to_owned()))
            .expect("chain should be valid");
        let residue = hierarchy
            .add_residue(
                chain,
                "GLY",
                Some(10),
                Some("42".to_owned()),
                Some("A".to_owned()),
            )
            .expect("residue should be valid");
        let metadata = AtomSiteMetadata {
            label_atom_id: Some("CA".to_owned()),
            auth_atom_id: Some("CAY".to_owned()),
            label_alt_id: Some("B".to_owned()),
            occupancy: Some(0.5),
            b_factor: Some(12.25),
        };
        let site = hierarchy
            .add_atom_site(residue, AtomId::new(7), metadata.clone())
            .expect("atom site should be valid");

        assert_eq!(model.raw(), 0);
        assert_eq!(chain.raw(), 0);
        assert_eq!(residue.raw(), 0);
        assert_eq!(site.raw(), 0);
        assert_eq!(
            hierarchy.model(model).expect("model exists").chains,
            vec![chain]
        );
        assert_eq!(
            hierarchy.chain(chain).expect("chain exists").residues,
            vec![residue]
        );
        assert_eq!(
            hierarchy
                .residue(residue)
                .expect("residue exists")
                .atom_sites,
            vec![site]
        );
        assert_eq!(
            hierarchy
                .atom_site_for_atom(AtomId::new(7))
                .expect("site exists")
                .metadata,
            metadata
        );
    }

    #[test]
    fn bio_hierarchy_iteration_is_insertion_order() {
        let mut hierarchy = BioHierarchy::new();
        let first_model = hierarchy.add_model("1");
        let second_model = hierarchy.add_model("2");
        let first_chain = hierarchy.add_chain(first_model, "A", None).expect("chain");
        let second_chain = hierarchy.add_chain(second_model, "B", None).expect("chain");

        assert_eq!(
            hierarchy.models().map(|(id, _)| id).collect::<Vec<_>>(),
            vec![first_model, second_model]
        );
        assert_eq!(
            hierarchy.chains().map(|(id, _)| id).collect::<Vec<_>>(),
            vec![first_chain, second_chain]
        );
    }

    #[test]
    fn bio_hierarchy_rejects_missing_parents_and_duplicate_atom_placement() {
        let mut hierarchy = BioHierarchy::new();
        assert_eq!(
            hierarchy
                .add_chain(ModelId::new(99), "A", None)
                .expect_err("missing model should fail"),
            BioHierarchyError::InvalidModelId(ModelId::new(99))
        );

        let model = hierarchy.add_model("1");
        let chain = hierarchy.add_chain(model, "A", None).expect("chain");
        assert_eq!(
            hierarchy
                .add_residue(ChainId::new(99), "GLY", None, None, None)
                .expect_err("missing chain should fail"),
            BioHierarchyError::InvalidChainId(ChainId::new(99))
        );
        let residue = hierarchy
            .add_residue(chain, "GLY", None, None, None)
            .expect("residue");
        let atom = AtomId::new(2);
        hierarchy
            .add_atom_site(residue, atom, AtomSiteMetadata::default())
            .expect("first placement should work");
        assert_eq!(
            hierarchy
                .add_atom_site(residue, atom, AtomSiteMetadata::default())
                .expect_err("duplicate atom placement should fail"),
            BioHierarchyError::DuplicateAtomPlacement(atom)
        );
    }

    #[test]
    fn macro_molecule_validates_atom_site_atom_ids() {
        let mut macro_mol = MacroMolecule::default();
        let atom = macro_mol.mol.add_atom(carbon());
        let model = macro_mol.hierarchy.add_model("1");
        let chain = macro_mol
            .hierarchy
            .add_chain(model, "A", Some("authA".to_owned()))
            .expect("chain");
        let residue = macro_mol
            .hierarchy
            .add_residue(chain, "ALA", Some(1), Some("1".to_owned()), None)
            .expect("residue");

        macro_mol
            .add_atom_site(
                residue,
                atom,
                AtomSiteMetadata {
                    label_atom_id: Some("CA".to_owned()),
                    auth_atom_id: Some("CA".to_owned()),
                    label_alt_id: None,
                    occupancy: Some(1.0),
                    b_factor: Some(10.0),
                },
            )
            .expect("valid atom should attach");
        assert_eq!(
            macro_mol
                .add_atom_site(residue, AtomId::new(99), AtomSiteMetadata::default())
                .expect_err("missing atom should fail"),
            BioHierarchyError::InvalidAtomId(AtomId::new(99))
        );
    }

    #[test]
    fn core_atom_does_not_store_biomolecular_labels() {
        let atom = carbon();

        assert_eq!(atom.element.symbol(), "C");
        assert!(atom.props.is_empty());
    }

    #[test]
    fn mmcif_parse_builds_macro_molecule_hierarchy() {
        let input = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA CAY . GLY GLY A X 10 42 A 0.50 12.25 1
ATOM 2 O O O . GLY GLY A X 10 42 A 1.00 10.00 1
"#;

        let macro_mol =
            read_mmcif_str(input, MmcifParseOptions::default()).expect("mmCIF should parse");

        assert_eq!(macro_mol.mol.atom_count(), 2);
        assert_eq!(macro_mol.mol.bond_count(), 0);
        assert_eq!(macro_mol.hierarchy.models().count(), 1);
        assert_eq!(macro_mol.hierarchy.chains().count(), 1);
        assert_eq!(macro_mol.hierarchy.residues().count(), 1);
        assert_eq!(macro_mol.hierarchy.atom_sites().count(), 2);
        let (_, chain) = macro_mol.hierarchy.chains().next().expect("chain exists");
        assert_eq!(chain.label_id, "A");
        assert_eq!(chain.author_id, Some("X".to_owned()));
        let (_, residue) = macro_mol
            .hierarchy
            .residues()
            .next()
            .expect("residue exists");
        assert_eq!(residue.name, "GLY");
        assert_eq!(residue.label_seq_id, Some(10));
        assert_eq!(residue.author_seq_id, Some("42".to_owned()));
        assert_eq!(residue.insertion_code, Some("A".to_owned()));
        let site = macro_mol
            .hierarchy
            .atom_site_for_atom(AtomId::new(0))
            .expect("site exists");
        assert_eq!(site.metadata.label_atom_id, Some("CA".to_owned()));
        assert_eq!(site.metadata.auth_atom_id, Some("CAY".to_owned()));
        assert_eq!(site.metadata.occupancy, Some(0.5));
        assert_eq!(site.metadata.b_factor, Some(12.25));
    }

    #[test]
    fn mmcif_parse_handles_missing_values_and_quotes() {
        let input = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_alt_id
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
C "C A" ? "LIG" "AA" . ? ? ? ? . 2
"#;

        let macro_mol =
            read_mmcif_str(input, MmcifParseOptions::default()).expect("mmCIF should parse");
        let (_, model) = macro_mol.hierarchy.models().next().expect("model exists");
        let site = macro_mol
            .hierarchy
            .atom_site_for_atom(AtomId::new(0))
            .expect("site exists");

        assert_eq!(model.model_id, "2");
        assert_eq!(site.metadata.label_atom_id, Some("C A".to_owned()));
        assert_eq!(site.metadata.auth_atom_id, None);
        assert_eq!(site.metadata.label_alt_id, None);
        assert_eq!(site.metadata.occupancy, None);
        assert_eq!(site.metadata.b_factor, None);
    }

    #[test]
    fn mmcif_parse_rejects_missing_strict_atom_id_and_unknown_element() {
        let missing_atom_id = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_comp_id
_atom_site.label_asym_id
C GLY A
"#;
        let err = read_mmcif_str(missing_atom_id, MmcifParseOptions::default())
            .expect_err("strict mode should require label atom id");
        assert!(err.message.contains("label atom id"));

        let unknown_element = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
Xx CA GLY A
"#;
        let err = read_mmcif_str(unknown_element, MmcifParseOptions::default())
            .expect_err("unknown element should fail");
        assert!(err.message.contains("unknown atom-site element"));
    }

    #[test]
    fn mmcif_parse_does_not_infer_bonds_or_perception() {
        let input = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
C C1 BEN A
C C2 BEN A
"#;

        let macro_mol =
            read_mmcif_str(input, MmcifParseOptions::default()).expect("mmCIF should parse");

        assert_eq!(macro_mol.mol.atom_count(), 2);
        assert_eq!(macro_mol.mol.bond_count(), 0);
        assert_eq!(macro_mol.mol.perception().rings, ComputedState::Absent);
        assert_eq!(
            macro_mol.mol.perception().aromaticity,
            ComputedState::Absent
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
