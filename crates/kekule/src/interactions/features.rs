//! Interaction features perceived once per topology.

use std::fmt;
use std::ops::Range;
use std::sync::Arc;

use crate::chemistry::perceive_molecule;
use crate::geometry::{Point3, Vector3};
use crate::io::is_metal;
use crate::molecule::{AtomId, BondOrder, Element, Molecule};
use crate::structure::AsModelView;
use crate::topology::{InstanceAtomId, MoleculeDefinitionId, Topology, TopologyAtomIndex};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};

use super::InteractionError;

/// The sign of a charged group's charge. Closed: a charge is positive or
/// negative.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ChargeSign {
    /// A cation.
    Positive,
    /// An anion.
    Negative,
}

/// How a charged group was recognized.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum ChargedGroupKind {
    /// A carboxylic acid or carboxylate: a carbon bonded to two terminal
    /// oxygens and at most one other heavy atom, a carbon. Negative.
    Carboxylate,
    /// A phosphate, phosphonate, or phosphoric ester with at least two
    /// terminal oxygens on one phosphorus. Negative.
    Phosphate,
    /// A sulfate or sulfonate: at least three terminal oxygens on one
    /// sulfur. Negative.
    Sulfonate,
    /// An aliphatic amine or ammonium: a nitrogen outside aromatic rings
    /// whose bonds are all single and whose heavy neighbors are all carbons
    /// without multiple bonds or aromaticity. Positive.
    Amine,
    /// An amidine: a non-aromatic carbon bonded to exactly two non-aromatic
    /// nitrogens, one of them by a double bond, and otherwise only to carbon
    /// or hydrogen. Positive.
    Amidine,
    /// A guanidine: a non-aromatic carbon bonded to three non-aromatic
    /// nitrogens, one of them by a double bond. Positive.
    Guanidine,
    /// Any other atom with a formal charge, unless it is a metal or bonded
    /// to an atom of opposite charge (as in a nitro group).
    FormalCharge,
}

/// Indices into the shared atom pool of [`InteractionFeatures`].
type Span = Range<u32>;

#[derive(Debug, Clone)]
pub(super) struct DonorRecord {
    pub(super) atom: TopologyAtomIndex,
    /// Explicit hydrogens in the pool.
    pub(super) hydrogens: Span,
    pub(super) implicit: u8,
}

#[derive(Debug, Clone)]
pub(super) struct GroupRecord {
    pub(super) sign: ChargeSign,
    pub(super) kind: ChargedGroupKind,
    /// All of the group's atoms; the charged atoms come first.
    pub(super) atoms: Span,
    pub(super) charged: u32,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct HalogenRecord {
    pub(super) halogen: TopologyAtomIndex,
    pub(super) carbon: TopologyAtomIndex,
}

/// The interaction features of one topology: hydrogen-bond donors and
/// acceptors, charged groups, aromatic rings, halogen-bond donors, metals
/// and the atoms that can coordinate them, and hydrophobic atoms.
///
/// Features depend only on the topology, so a host perceives them once and
/// reuses them for every structure sharing its layout: a model, each sample
/// of an ensemble, and each frame of a trajectory. Each reusable definition
/// is classified once and mapped onto its instances.
///
/// Classification reads elements, formal charges, bonds, hydrogen counts,
/// and aromaticity. A definition with installed valence, ring-set, and
/// aromaticity perception uses it; any other is perceived on a private copy
/// with the default profile, and the topology is left unchanged. A
/// definition whose perception fails keeps its represented state: it has no
/// aromatic rings, and inferred hydrogen counts fall back to standard
/// valences; [`Self::unperceived_definitions`] lists it.
///
/// - **Donors**: nitrogen and oxygen atoms bearing hydrogens, with their
///   explicit hydrogen atoms and their implicit hydrogen count, except the
///   oxygens of carboxylates, phosphates, and sulfonates, which are
///   deprotonated at pH 7 however represented.
/// - **Acceptors**: oxygen atoms without positive charge; nitrogen atoms
///   without hydrogens or positive charge that are not amide, sulfonamide,
///   or aromatic nitrogens with three connections (pyrrole-type), and have
///   fewer than four connections.
/// - **Charged groups**: carboxylates, phosphates, sulfonates (charged at
///   pH 7 whether or not they are represented protonated), aliphatic amines,
///   amidines, guanidines, and other formal charges; see
///   [`ChargedGroupKind`]. Each atom belongs to at most one group. Templates
///   such as the PDB's represent histidine protonated, so its ring nitrogen
///   is a [`ChargedGroupKind::FormalCharge`] cation.
/// - **Aromatic rings**: the perceived ring set's rings whose atoms are all
///   aromatic; fused rings are separate rings.
/// - **Halogen-bond donors**: chlorine, bromine, and iodine bonded to one
///   carbon.
/// - **Metals**: alkali, alkaline-earth, transition, and post-transition
///   metals, lanthanides, and actinides.
/// - **Metal ligands**: nitrogen, oxygen, sulfur, and selenium atoms, except
///   nitrogens with four connections, and halide ions.
/// - **Hydrophobic atoms**: uncharged carbons whose heavy neighbors are all
///   carbons.
#[derive(Clone)]
pub struct InteractionFeatures {
    topology: Arc<Topology>,
    pub(super) pool: Vec<TopologyAtomIndex>,
    pub(super) donors: Vec<DonorRecord>,
    pub(super) acceptors: Vec<TopologyAtomIndex>,
    pub(super) groups: Vec<GroupRecord>,
    pub(super) rings: Vec<Span>,
    pub(super) halogens: Vec<HalogenRecord>,
    pub(super) metals: Vec<TopologyAtomIndex>,
    pub(super) ligands: Vec<TopologyAtomIndex>,
    pub(super) hydrophobic: Vec<TopologyAtomIndex>,
    unperceived: Vec<MoleculeDefinitionId>,
}

impl fmt::Debug for InteractionFeatures {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("InteractionFeatures")
            .field("atom_count", &self.topology.atom_count())
            .field("donors", &self.donors.len())
            .field("acceptors", &self.acceptors.len())
            .field("charged_groups", &self.groups.len())
            .field("aromatic_rings", &self.rings.len())
            .field("halogen_donors", &self.halogens.len())
            .field("metals", &self.metals.len())
            .field("metal_ligands", &self.ligands.len())
            .field("hydrophobic_atoms", &self.hydrophobic.len())
            .field("unperceived_definitions", &self.unperceived)
            .finish()
    }
}

impl InteractionFeatures {
    /// Perceives the features of `topology`.
    ///
    /// ```
    /// use kekule::interactions::InteractionFeatures;
    ///
    /// // Acetic acid and methylamine: a carboxylate and an amine.
    /// let topology = kekule::smiles::to_topology("CC(=O)O.CN")?;
    /// let features = InteractionFeatures::new(topology);
    /// assert_eq!(features.charged_groups().len(), 2);
    /// assert_eq!(features.donors().len(), 2);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        let topology = topology.into();
        let mut features = Self {
            topology: Arc::clone(&topology),
            pool: Vec::new(),
            donors: Vec::new(),
            acceptors: Vec::new(),
            groups: Vec::new(),
            rings: Vec::new(),
            halogens: Vec::new(),
            metals: Vec::new(),
            ligands: Vec::new(),
            hydrophobic: Vec::new(),
            unperceived: Vec::new(),
        };
        for definition in topology.definitions() {
            let molecule = definition.molecule();
            let perception = molecule.perception();
            let installed = perception.valence().is_some()
                && perception.aromaticity().is_some()
                && perception
                    .rings()
                    .is_some_and(|rings| rings.ring_set().is_some());
            let copy;
            let (molecule, perceived) = if installed {
                (molecule, true)
            } else {
                let mut perceived = molecule.clone();
                match perceive_molecule(&mut perceived) {
                    Ok(()) => {
                        copy = perceived;
                        (&copy, true)
                    }
                    Err(_) => {
                        features.unperceived.push(definition.id());
                        (molecule, false)
                    }
                }
            };
            let local = Local::classify(molecule, perceived);
            for instance in topology.instances_of(definition.id()) {
                let dense = |atom: AtomId| {
                    topology
                        .atom_index(InstanceAtomId::new(instance.id(), atom))
                        .expect("instance atoms are topology atoms")
                };
                features.extend(&local, dense);
            }
        }
        // Single-atom features sort by atom, so a search finds an atom's
        // features by bisection.
        let by_atom = |a: &TopologyAtomIndex, b: &TopologyAtomIndex| a.cmp(b);
        features.donors.sort_unstable_by_key(|donor| donor.atom);
        features
            .halogens
            .sort_unstable_by_key(|record| record.halogen);
        features.acceptors.sort_unstable_by(by_atom);
        features.metals.sort_unstable_by(by_atom);
        features.ligands.sort_unstable_by(by_atom);
        features.hydrophobic.sort_unstable_by(by_atom);
        features
    }

    fn extend(&mut self, local: &Local, dense: impl Fn(AtomId) -> TopologyAtomIndex) {
        let span = |pool: &mut Vec<TopologyAtomIndex>, atoms: &[AtomId]| {
            let start = pool.len() as u32;
            pool.extend(atoms.iter().map(|&atom| dense(atom)));
            start..pool.len() as u32
        };
        for (atom, hydrogens, implicit) in &local.donors {
            let hydrogens = span(&mut self.pool, hydrogens);
            self.donors.push(DonorRecord {
                atom: dense(*atom),
                hydrogens,
                implicit: *implicit,
            });
        }
        for (sign, kind, atoms, charged) in &local.groups {
            let atoms = span(&mut self.pool, atoms);
            self.groups.push(GroupRecord {
                sign: *sign,
                kind: *kind,
                atoms,
                charged: *charged,
            });
        }
        for ring in &local.rings {
            let ring = span(&mut self.pool, ring);
            self.rings.push(ring);
        }
        self.halogens.extend(
            local
                .halogens
                .iter()
                .map(|&(halogen, carbon)| HalogenRecord {
                    halogen: dense(halogen),
                    carbon: dense(carbon),
                }),
        );
        self.acceptors
            .extend(local.acceptors.iter().map(|&atom| dense(atom)));
        self.metals
            .extend(local.metals.iter().map(|&atom| dense(atom)));
        self.ligands
            .extend(local.ligands.iter().map(|&atom| dense(atom)));
        self.hydrophobic
            .extend(local.hydrophobic.iter().map(|&atom| dense(atom)));
    }

    /// The topology the features belong to.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// The topology the features belong to, shared.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub(super) fn span(&self, span: &Span) -> &[TopologyAtomIndex] {
        &self.pool[span.start as usize..span.end as usize]
    }

    /// Hydrogen-bond donors in dense order.
    pub fn donors(&self) -> impl ExactSizeIterator<Item = HydrogenBondDonor<'_>> + '_ {
        self.donors.iter().map(|donor| HydrogenBondDonor {
            atom: donor.atom,
            hydrogens: self.span(&donor.hydrogens),
            implicit: donor.implicit,
        })
    }

    /// Hydrogen-bond acceptors in dense order.
    pub fn acceptors(&self) -> &[TopologyAtomIndex] {
        &self.acceptors
    }

    /// Charged groups, in definition and then instance order.
    pub fn charged_groups(&self) -> impl ExactSizeIterator<Item = ChargedGroup<'_>> + '_ {
        self.groups.iter().map(|group| ChargedGroup {
            sign: group.sign,
            kind: group.kind,
            atoms: self.span(&group.atoms),
            charged: group.charged as usize,
        })
    }

    /// Aromatic rings, in definition and then instance order.
    pub fn aromatic_rings(&self) -> impl ExactSizeIterator<Item = AromaticRing<'_>> + '_ {
        self.rings.iter().map(|ring| AromaticRing {
            atoms: self.span(ring),
            topology: &self.topology,
        })
    }

    /// Halogen-bond donors in dense order of their halogens.
    pub fn halogen_donors(&self) -> impl ExactSizeIterator<Item = HalogenBondDonor> + '_ {
        self.halogens.iter().map(|record| HalogenBondDonor {
            halogen: record.halogen,
            carbon: record.carbon,
        })
    }

    /// Metal atoms in dense order.
    pub fn metals(&self) -> &[TopologyAtomIndex] {
        &self.metals
    }

    /// Atoms that can coordinate metals, in dense order.
    pub fn metal_ligands(&self) -> &[TopologyAtomIndex] {
        &self.ligands
    }

    /// Hydrophobic atoms in dense order.
    pub fn hydrophobic_atoms(&self) -> &[TopologyAtomIndex] {
        &self.hydrophobic
    }

    /// The definitions whose perception failed, classified from their
    /// represented state alone.
    pub fn unperceived_definitions(&self) -> &[MoleculeDefinitionId] {
        &self.unperceived
    }
}

/// A hydrogen-bond donor of [`InteractionFeatures`].
#[derive(Debug, Clone, Copy)]
pub struct HydrogenBondDonor<'a> {
    atom: TopologyAtomIndex,
    hydrogens: &'a [TopologyAtomIndex],
    implicit: u8,
}

impl<'a> HydrogenBondDonor<'a> {
    /// The donor atom.
    pub const fn atom(&self) -> TopologyAtomIndex {
        self.atom
    }

    /// Its explicit hydrogen atoms.
    pub const fn hydrogens(&self) -> &'a [TopologyAtomIndex] {
        self.hydrogens
    }

    /// Its implicit hydrogens, whose positions are unknown.
    pub const fn implicit_hydrogen_count(&self) -> usize {
        self.implicit as usize
    }
}

/// A charged group of [`InteractionFeatures`].
#[derive(Debug, Clone, Copy)]
pub struct ChargedGroup<'a> {
    sign: ChargeSign,
    kind: ChargedGroupKind,
    atoms: &'a [TopologyAtomIndex],
    charged: usize,
}

impl<'a> ChargedGroup<'a> {
    /// The sign of the group's charge.
    pub const fn sign(&self) -> ChargeSign {
        self.sign
    }

    /// How the group was recognized.
    pub const fn kind(&self) -> ChargedGroupKind {
        self.kind
    }

    /// The group's atoms: its charged atoms, then the others, such as a
    /// carboxylate's carbon.
    pub const fn atoms(&self) -> &'a [TopologyAtomIndex] {
        self.atoms
    }

    /// The atoms whose centroid is the group's charge center: a
    /// carboxylate's, phosphate's, or sulfonate's terminal oxygens, an
    /// amidine's or guanidine's nitrogens, or the charged atom.
    pub fn charged_atoms(&self) -> &'a [TopologyAtomIndex] {
        &self.atoms[..self.charged]
    }
}

/// An aromatic ring of [`InteractionFeatures`].
#[derive(Debug, Clone, Copy)]
pub struct AromaticRing<'a> {
    atoms: &'a [TopologyAtomIndex],
    topology: &'a Topology,
}

impl<'a> AromaticRing<'a> {
    /// The ring's atoms, in order around the ring.
    pub const fn atoms(&self) -> &'a [TopologyAtomIndex] {
        self.atoms
    }

    /// The ring's center, the centroid of its atoms, and its unit normal in
    /// `structure`.
    ///
    /// # Errors
    ///
    /// [`InteractionError::TopologyMismatch`] when `structure` does not
    /// share the features' topology layout, and
    /// [`InteractionError::DegenerateRing`] when the ring's atoms are
    /// collinear or coincident, so it has no plane.
    pub fn center_and_normal(
        &self,
        structure: &(impl AsModelView + ?Sized),
    ) -> Result<(Quantity<Point3>, Vector3), InteractionError> {
        let structure = structure.as_model_view();
        if !structure.topology().shares_layout(self.topology) {
            return Err(InteractionError::TopologyMismatch);
        }
        let positions = structure.positions().values();
        let (center, normal) =
            ring_geometry(self.atoms, positions.value()).ok_or(InteractionError::DegenerateRing)?;
        Ok((Quantity::new(center, CANONICAL_LENGTH_UNIT), normal))
    }
}

/// The centroid of a ring's atoms and its unit normal by Newell's method,
/// or `None` for a ring without a plane.
pub(super) fn ring_geometry(
    atoms: &[TopologyAtomIndex],
    positions: &[Point3],
) -> Option<(Point3, Vector3)> {
    let points = atoms.iter().map(|atom| positions[atom.index()]);
    let center = centroid(points.clone());
    let mut normal = Vector3::zero();
    let count = atoms.len();
    for (index, atom) in atoms.iter().enumerate() {
        let a = positions[atom.index()] - center;
        let b = positions[atoms[(index + 1) % count].index()] - center;
        normal += a.cross(b);
    }
    let length = normal.norm();
    (length.is_normal()).then(|| (center, normal / length))
}

/// The centroid of nonempty points.
pub(super) fn centroid(points: impl ExactSizeIterator<Item = Point3>) -> Point3 {
    let count = points.len() as f64;
    let sum = points.fold(Vector3::zero(), |sum, point| {
        sum + Vector3::new(point.x, point.y, point.z)
    });
    Point3::new(sum.x / count, sum.y / count, sum.z / count)
}

/// A halogen-bond donor of [`InteractionFeatures`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HalogenBondDonor {
    halogen: TopologyAtomIndex,
    carbon: TopologyAtomIndex,
}

impl HalogenBondDonor {
    /// The halogen atom.
    pub const fn halogen(&self) -> TopologyAtomIndex {
        self.halogen
    }

    /// The carbon it is bonded to.
    pub const fn carbon(&self) -> TopologyAtomIndex {
        self.carbon
    }
}

/// The features of one definition, by local atom.
#[derive(Debug, Default)]
struct Local {
    donors: Vec<(AtomId, Vec<AtomId>, u8)>,
    acceptors: Vec<AtomId>,
    groups: Vec<(ChargeSign, ChargedGroupKind, Vec<AtomId>, u32)>,
    rings: Vec<Vec<AtomId>>,
    halogens: Vec<(AtomId, AtomId)>,
    metals: Vec<AtomId>,
    ligands: Vec<AtomId>,
    hydrophobic: Vec<AtomId>,
}

/// One atom's chemistry as classification reads it.
#[derive(Debug, Clone)]
struct Site {
    element: Element,
    charge: i8,
    aromatic: bool,
    /// Heavy neighbors with their bond orders.
    heavy: Vec<(AtomId, BondOrder)>,
    explicit: Vec<AtomId>,
    hydrogens: usize,
}

impl Site {
    fn connections(&self) -> usize {
        self.heavy.len() + self.hydrogens
    }

    fn all_single(&self) -> bool {
        self.heavy
            .iter()
            .all(|(_, order)| *order == BondOrder::Single)
    }

    fn has_multiple(&self) -> bool {
        self.heavy
            .iter()
            .any(|(_, order)| matches!(order, BondOrder::Double | BondOrder::Triple))
    }
}

impl Local {
    fn classify(molecule: &Molecule, perceived: bool) -> Self {
        let aromaticity = molecule.perception().aromaticity();
        let sites = molecule
            .atom_ids()
            .map(|id| {
                let atom = molecule.atom(id).expect("listed atoms are live");
                let mut heavy = Vec::new();
                let mut explicit = Vec::new();
                for (_, bond) in molecule.incident_bonds(id).expect("listed atoms are live") {
                    let other = bond.other_atom(id);
                    let element = molecule
                        .atom(other)
                        .expect("bonds join live atoms")
                        .element();
                    // Coordination bonds to metals are interactions in their
                    // own right, not part of a group's covalent structure.
                    if element == Element::H {
                        explicit.push(other);
                    } else if !is_metal(element) || is_metal(atom.element()) {
                        heavy.push((other, bond.order));
                    }
                }
                let hydrogens = molecule.total_hydrogens(id).unwrap_or_else(|| {
                    explicit.len()
                        + standard_implicit(
                            atom.element(),
                            atom.formal_charge(),
                            &heavy,
                            explicit.len(),
                        )
                });
                Site {
                    element: atom.element(),
                    charge: atom.formal_charge(),
                    aromatic: aromaticity
                        .is_some_and(|aromaticity| aromaticity.atom_is_aromatic(id)),
                    heavy,
                    explicit,
                    hydrogens,
                }
            })
            .collect::<Vec<_>>();
        let site = |id: AtomId| &sites[id.index()];
        let mut local = Self {
            groups: charged_groups(&sites),
            ..Self::default()
        };
        // An acid's oxygens are deprotonated at pH 7, however represented.
        let mut anionic = vec![false; sites.len()];
        for (sign, _, atoms, charged) in &local.groups {
            if *sign == ChargeSign::Negative {
                for atom in &atoms[..*charged as usize] {
                    anionic[atom.index()] = true;
                }
            }
        }

        for (index, atom) in sites.iter().enumerate() {
            let id = AtomId::new(index as u32);
            let element = atom.element;
            if element == Element::H {
                continue;
            }
            // Donors and acceptors.
            if matches!(element, Element::N | Element::O) && atom.hydrogens > 0 && !anionic[index] {
                let implicit = atom.hydrogens.saturating_sub(atom.explicit.len());
                local.donors.push((
                    id,
                    atom.explicit.clone(),
                    implicit.min(u8::MAX as usize) as u8,
                ));
            }
            if element == Element::O && atom.charge <= 0 {
                local.acceptors.push(id);
            }
            if element == Element::N && is_nitrogen_acceptor(atom, &sites) {
                local.acceptors.push(id);
            }
            // Metals and their ligands.
            if is_metal(element) {
                local.metals.push(id);
            } else if matches!(element, Element::O | Element::S | Element::SE)
                || (element == Element::N && atom.connections() < 4)
                || (is_halogen(element) && atom.heavy.is_empty())
            {
                local.ligands.push(id);
            }
            // Halogen-bond donors.
            if matches!(element, Element::CL | Element::BR | Element::I)
                && atom.heavy.len() == 1
                && site(atom.heavy[0].0).element == Element::C
            {
                local.halogens.push((id, atom.heavy[0].0));
            }
            // Hydrophobic carbons.
            if element == Element::C
                && atom.charge == 0
                && atom
                    .heavy
                    .iter()
                    .all(|(other, _)| site(*other).element == Element::C)
            {
                local.hydrophobic.push(id);
            }
        }

        if perceived {
            if let Some(rings) = molecule
                .perception()
                .rings()
                .and_then(|rings| rings.ring_set())
            {
                for ring in rings.rings() {
                    if ring.atoms.iter().all(|&atom| site(atom).aromatic) {
                        local.rings.push(ring.atoms.clone());
                    }
                }
            }
        }
        local
    }
}

fn is_halogen(element: Element) -> bool {
    matches!(element, Element::F | Element::CL | Element::BR | Element::I)
}

/// The implicit hydrogens of an atom whose count is unknown, from the
/// standard valence of its element and charge.
fn standard_implicit(
    element: Element,
    charge: i8,
    heavy: &[(AtomId, BondOrder)],
    explicit: usize,
) -> usize {
    let valence: i32 = match element {
        Element::C => 4 - i32::from(charge).abs(),
        Element::N | Element::P => 3 + i32::from(charge),
        Element::O | Element::S | Element::SE => 2 + i32::from(charge),
        _ => return 0,
    };
    let used: i32 = heavy
        .iter()
        .map(|(_, order)| match order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            _ => 0,
        })
        .sum::<i32>()
        + explicit as i32;
    (valence - used).max(0) as usize
}

fn is_nitrogen_acceptor(atom: &Site, sites: &[Site]) -> bool {
    let site = |id: AtomId| &sites[id.index()];
    if atom.charge > 0 || atom.hydrogens > 0 || atom.connections() >= 4 {
        return false;
    }
    // A pyrrole-type nitrogen gives its lone pair to the ring.
    if atom.aromatic && atom.connections() == 3 {
        return false;
    }
    // An amide or sulfonamide nitrogen's lone pair is delocalized.
    let delocalized = atom.all_single()
        && atom.heavy.iter().any(|(other, _)| {
            let other = site(*other);
            matches!(other.element, Element::C | Element::S | Element::P)
                && !other.aromatic
                && other.heavy.iter().any(|(next, order)| {
                    *order == BondOrder::Double
                        && matches!(site(*next).element, Element::O | Element::S)
                })
        });
    !delocalized
}

/// Recognizes charged groups, each atom in at most one.
fn charged_groups(sites: &[Site]) -> Vec<(ChargeSign, ChargedGroupKind, Vec<AtomId>, u32)> {
    let site = |id: AtomId| &sites[id.index()];
    let mut claimed = vec![false; sites.len()];
    let mut groups = Vec::new();
    let terminal_oxygens = |atom: &Site| {
        atom.heavy
            .iter()
            .filter(|(other, _)| {
                let other = site(*other);
                other.element == Element::O && other.heavy.len() == 1
            })
            .map(|(other, _)| *other)
            .collect::<Vec<_>>()
    };
    let claim = |groups: &mut Vec<_>,
                 claimed: &mut Vec<bool>,
                 sign: ChargeSign,
                 kind: ChargedGroupKind,
                 charged: Vec<AtomId>,
                 others: &[AtomId]| {
        let atoms = charged.iter().chain(others).copied().collect::<Vec<_>>();
        if atoms.iter().any(|atom| claimed[atom.index()]) {
            return;
        }
        for atom in &atoms {
            claimed[atom.index()] = true;
        }
        groups.push((sign, kind, atoms, charged.len() as u32));
    };
    // Acids first, then bases, then the remaining formal charges.
    for (index, atom) in sites.iter().enumerate() {
        let id = AtomId::new(index as u32);
        let oxygens = terminal_oxygens(atom);
        let kind = match atom.element {
            Element::C
                if oxygens.len() == 2
                    && atom.heavy.len() <= 3
                    && atom.heavy.iter().all(|(other, _)| {
                        oxygens.contains(other) || site(*other).element == Element::C
                    }) =>
            {
                Some(ChargedGroupKind::Carboxylate)
            }
            Element::P if oxygens.len() >= 2 => Some(ChargedGroupKind::Phosphate),
            Element::S if oxygens.len() >= 3 => Some(ChargedGroupKind::Sulfonate),
            _ => None,
        };
        if let Some(kind) = kind {
            claim(
                &mut groups,
                &mut claimed,
                ChargeSign::Negative,
                kind,
                oxygens,
                &[id],
            );
        }
    }
    for (index, atom) in sites.iter().enumerate() {
        let id = AtomId::new(index as u32);
        if atom.element != Element::C || atom.aromatic {
            continue;
        }
        let nitrogens = atom
            .heavy
            .iter()
            .filter(|(other, _)| site(*other).element == Element::N && !site(*other).aromatic)
            .collect::<Vec<_>>();
        let double = nitrogens
            .iter()
            .any(|(_, order)| *order == BondOrder::Double);
        let others_plain = atom
            .heavy
            .iter()
            .filter(|(other, _)| site(*other).element != Element::N)
            .all(|(other, order)| {
                *order == BondOrder::Single && site(*other).element == Element::C
            });
        let kind = match nitrogens.len() {
            3 if double && atom.heavy.len() == 3 => Some(ChargedGroupKind::Guanidine),
            2 if double && others_plain => Some(ChargedGroupKind::Amidine),
            _ => None,
        };
        if let Some(kind) = kind {
            let charged = nitrogens.iter().map(|(other, _)| *other).collect();
            claim(
                &mut groups,
                &mut claimed,
                ChargeSign::Positive,
                kind,
                charged,
                &[id],
            );
        }
    }
    for (index, atom) in sites.iter().enumerate() {
        let id = AtomId::new(index as u32);
        let amine = atom.element == Element::N
            && !atom.aromatic
            && atom.charge >= 0
            && atom.all_single()
            && atom.heavy.iter().all(|(other, _)| {
                let other = site(*other);
                other.element == Element::C && !other.aromatic && !other.has_multiple()
            });
        if amine {
            claim(
                &mut groups,
                &mut claimed,
                ChargeSign::Positive,
                ChargedGroupKind::Amine,
                vec![id],
                &[],
            );
        }
    }
    for (index, atom) in sites.iter().enumerate() {
        let id = AtomId::new(index as u32);
        if atom.charge == 0 || atom.element == Element::H || is_metal(atom.element) {
            continue;
        }
        let paired = atom
            .heavy
            .iter()
            .any(|(other, _)| site(*other).charge.signum() == -atom.charge.signum());
        if paired {
            continue;
        }
        let sign = if atom.charge > 0 {
            ChargeSign::Positive
        } else {
            ChargeSign::Negative
        };
        claim(
            &mut groups,
            &mut claimed,
            sign,
            ChargedGroupKind::FormalCharge,
            vec![id],
            &[],
        );
    }
    groups
}
