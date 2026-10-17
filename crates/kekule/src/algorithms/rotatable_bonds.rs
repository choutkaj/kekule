use std::collections::{BTreeMap, BTreeSet};

use crate::chemistry::PerceptionError;
use crate::molecule::{
    AromaticityModel, AtomId, Bond, BondId, BondOrder, Element, ImplicitHydrogens, Molecule,
    RingBasisState, RingPerception, ValenceModel, ValenceValidation,
};

use super::rings::compute_ring_membership;

const HYDROGEN: u8 = 1;
const CARBON: u8 = 6;
const NITROGEN: u8 = 7;
const OXYGEN: u8 = 8;
const FLUORINE: u8 = 9;
const SULFUR: u8 = 16;
const CHLORINE: u8 = 17;
const BROMINE: u8 = 35;

/// Options controlling which represented single bonds are returned.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RotatableBondOptions {
    include_terminal_bonds: bool,
    include_resonance_restricted_bonds: bool,
    include_symmetric_bonds: bool,
    include_ring_bonds: bool,
}

impl RotatableBondOptions {
    /// Includes bonds whose endpoint has no other covalently connected
    /// non-hydrogen neighbor.
    pub const fn include_terminal_bonds(&self) -> bool {
        self.include_terminal_bonds
    }

    /// Sets [`Self::include_terminal_bonds`].
    #[must_use]
    pub const fn with_include_terminal_bonds(mut self, include_terminal_bonds: bool) -> Self {
        self.include_terminal_bonds = include_terminal_bonds;
        self
    }

    /// Includes resonance-restricted linkages such as amides and esters.
    pub const fn include_resonance_restricted_bonds(&self) -> bool {
        self.include_resonance_restricted_bonds
    }

    /// Sets [`Self::include_resonance_restricted_bonds`].
    #[must_use]
    pub const fn with_include_resonance_restricted_bonds(
        mut self,
        include_resonance_restricted_bonds: bool,
    ) -> Self {
        self.include_resonance_restricted_bonds = include_resonance_restricted_bonds;
        self
    }

    /// Includes axes involving symmetric groups such as trihalomethyl and
    /// tert-butyl-like centers.
    pub const fn include_symmetric_bonds(&self) -> bool {
        self.include_symmetric_bonds
    }

    /// Sets [`Self::include_symmetric_bonds`].
    #[must_use]
    pub const fn with_include_symmetric_bonds(mut self, include_symmetric_bonds: bool) -> Self {
        self.include_symmetric_bonds = include_symmetric_bonds;
        self
    }

    /// Includes represented single bonds that are members of graph cycles.
    pub const fn include_ring_bonds(&self) -> bool {
        self.include_ring_bonds
    }

    /// Sets [`Self::include_ring_bonds`].
    #[must_use]
    pub const fn with_include_ring_bonds(mut self, include_ring_bonds: bool) -> Self {
        self.include_ring_bonds = include_ring_bonds;
        self
    }
}

impl RotatableBondOptions {
    /// RDKit's strict rotatable-bond descriptor version 3.2.0, adapted to
    /// ignore graph hydrogens when evaluating terminal and symmetric groups.
    /// Resonance exclusions use the default RDKit-like aromaticity perception.
    pub const STRICT: Self = Self {
        include_terminal_bonds: false,
        include_resonance_restricted_bonds: false,
        include_symmetric_bonds: false,
        include_ring_bonds: false,
    };

    /// A permissive selection of all supported represented single-bond axes.
    ///
    /// This enables every optional category, including ring and terminal
    /// bonds. Unsupported bond orders, hydrogen axes, and axes adjacent to a
    /// triple bond remain excluded.
    pub const GENERAL: Self = Self {
        include_terminal_bonds: true,
        include_resonance_restricted_bonds: true,
        include_symmetric_bonds: true,
        include_ring_bonds: true,
    };
}

impl Default for RotatableBondOptions {
    fn default() -> Self {
        Self::STRICT
    }
}

/// A detached snapshot of the rotatable bonds detected in one molecule.
///
/// Bond identifiers remain in ascending stable-ID order. The snapshot is not
/// updated when its source molecule changes and should be recomputed after any
/// chemistry or topology mutation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RotatableBondSet {
    options: RotatableBondOptions,
    bonds: Vec<BondId>,
}

impl RotatableBondSet {
    /// Returns the options used to produce this snapshot.
    pub const fn options(&self) -> RotatableBondOptions {
        self.options
    }

    /// Returns detected bonds in ascending stable-ID order.
    pub fn bond_ids(&self) -> &[BondId] {
        &self.bonds
    }

    /// Returns whether the snapshot contains a bond identifier.
    pub fn contains(&self, bond: BondId) -> bool {
        self.bonds.binary_search(&bond).is_ok()
    }

    /// Returns the number of detected rotatable bonds.
    pub const fn len(&self) -> usize {
        self.bonds.len()
    }

    /// Returns whether no rotatable bonds were detected.
    pub const fn is_empty(&self) -> bool {
        self.bonds.is_empty()
    }
}

/// Detects rotatable bonds, the configurable represented single-bond axes of
/// canonical chemistry, without mutating the molecule or installing perception.
///
/// Ring membership is reused or computed on demand. Resonance exclusions reuse
/// installed RDKit-like valence and aromaticity with the default Figueras ring
/// model, or compute the default profile on a temporary copy. Strict valence
/// validation applies in either case: installed valence is reused only when
/// its provenance records strict validation. Incomplete or incompatible
/// perception is not reused. Valence and resource failures are
/// returned without changing the source graph or its perception. Disabling
/// resonance exclusions requires only graph ring membership and does not run
/// chemical perception.
///
/// # Errors
///
/// When resonance exclusions apply, the [`PerceptionError`] of strict valence
/// validation or of the default perception profile on the temporary copy.
pub fn rotatable_bonds(
    molecule: &Molecule,
    options: RotatableBondOptions,
) -> Result<RotatableBondSet, PerceptionError> {
    let mut perceived;
    let perception = molecule.perception();
    // Cached state is reused only when it was produced by the models and the
    // strict validation this descriptor would perceive with itself, so it
    // cannot hide a failure that an unperceived molecule would report.
    let valence_matches = perception.valence().is_some_and(|valence| {
        matches!(valence.model(), ValenceModel::RdkitLike)
            && valence.validation() == ValenceValidation::Strict
    });
    let aromaticity_matches = perception
        .aromaticity()
        .is_some_and(|aromaticity| aromaticity.model() == AromaticityModel::RdkitLike);
    let compatible = valence_matches
        && perception
            .rings()
            .and_then(RingPerception::basis)
            .is_some_and(RingBasisState::is_sssr_like)
        && aromaticity_matches;
    let molecule = if !options.include_resonance_restricted_bonds && !compatible {
        perceived = molecule.clone();
        perceived.perceive()?;
        &perceived
    } else {
        molecule
    };
    let computed_membership;
    let ring_membership = match molecule
        .perception()
        .rings()
        .map(|rings| rings.membership())
    {
        Some(membership) => membership,
        None => {
            computed_membership = compute_ring_membership(molecule);
            &computed_membership
        }
    };
    let profiles = atom_profiles(molecule, ring_membership, options);
    let bonds = molecule
        .bonds()
        .filter_map(|(bond_id, bond)| {
            let rotatable = rotatable_bond(bond_id, bond, ring_membership, &profiles, options);
            rotatable.then_some(bond_id)
        })
        .collect();

    Ok(RotatableBondSet { options, bonds })
}

#[derive(Debug, Clone, Copy, Default)]
struct AtomProfile {
    is_hydrogen: bool,
    heavy_degree: usize,
    incident_triple: bool,
    trihalomethyl_center: bool,
    tert_butyl_like_center: bool,
    resonance_restricted: bool,
}

fn rotatable_bond(
    bond_id: BondId,
    bond: &Bond,
    ring_membership: &crate::molecule::RingMembership,
    profiles: &BTreeMap<AtomId, AtomProfile>,
    options: RotatableBondOptions,
) -> bool {
    if bond.order != BondOrder::Single
        || (!options.include_ring_bonds && ring_membership.bond_in_ring(bond_id))
    {
        return false;
    }
    let (left, right) = bond.endpoints();
    endpoint_is_eligible(left, profiles, options)
        && endpoint_is_eligible(right, profiles, options)
        // RDKit 3.2.0 puts the recursive resonance exclusions on only one
        // side of its query. A match can orient either endpoint to that side,
        // so the bond disappears only when both endpoints are excluded.
        && (options.include_resonance_restricted_bonds
            || !(endpoint_is_resonance_restricted(left, profiles)
                && endpoint_is_resonance_restricted(right, profiles)))
}

fn endpoint_is_eligible(
    atom: AtomId,
    profiles: &BTreeMap<AtomId, AtomProfile>,
    options: RotatableBondOptions,
) -> bool {
    profiles.get(&atom).is_some_and(|profile| {
        !profile.is_hydrogen
            && (options.include_terminal_bonds || profile.heavy_degree > 1)
            && !profile.incident_triple
            && (options.include_symmetric_bonds
                || (!profile.trihalomethyl_center && !profile.tert_butyl_like_center))
    })
}

fn endpoint_is_resonance_restricted(
    atom: AtomId,
    profiles: &BTreeMap<AtomId, AtomProfile>,
) -> bool {
    profiles
        .get(&atom)
        .is_some_and(|profile| profile.resonance_restricted)
}

fn atom_profiles(
    molecule: &Molecule,
    ring_membership: &crate::molecule::RingMembership,
    options: RotatableBondOptions,
) -> BTreeMap<AtomId, AtomProfile> {
    let mut profiles = molecule
        .atom_ids()
        .map(|atom| {
            let mut profile = AtomProfile {
                is_hydrogen: molecule
                    .atom(atom)
                    .is_some_and(|atom| atom.element.atomic_number() == HYDROGEN),
                ..AtomProfile::default()
            };
            for (_, bond) in molecule
                .incident_bonds(atom)
                .expect("live atoms have valid adjacency")
            {
                if !is_ordinary_covalent(bond.order) {
                    continue;
                }
                let neighbor = bond.other_atom(atom);
                if molecule
                    .atom(neighbor)
                    .is_some_and(|neighbor| neighbor.element.atomic_number() != HYDROGEN)
                {
                    profile.heavy_degree += 1;
                }
                profile.incident_triple |= bond.order == BondOrder::Triple;
            }
            (atom, profile)
        })
        .collect::<BTreeMap<_, _>>();

    for atom in molecule.atom_ids() {
        let Some(value) = molecule.atom(atom) else {
            continue;
        };
        if value.element.atomic_number() != CARBON {
            continue;
        }
        let trihalomethyl_center = is_trihalomethyl_center(molecule, atom);
        let tert_butyl_like_center = molecule
            .incident_bonds(atom)
            .expect("live atoms have valid adjacency")
            .filter(|(_, bond)| bond.order == BondOrder::Single)
            .map(|(_, bond)| bond.other_atom(atom))
            .filter(|neighbor| is_methyl_like(molecule, *neighbor, &profiles))
            .take(3)
            .count()
            == 3;
        if let Some(profile) = profiles.get_mut(&atom) {
            profile.trihalomethyl_center = trihalomethyl_center;
            profile.tert_butyl_like_center = tert_butyl_like_center;
        }
    }
    if !options.include_resonance_restricted_bonds {
        for atom in resonance_restricted_atoms(
            molecule,
            ring_membership,
            &profiles,
            options.include_ring_bonds,
        ) {
            if let Some(profile) = profiles.get_mut(&atom) {
                profile.resonance_restricted = true;
            }
        }
    }
    profiles
}

fn is_trihalomethyl_center(molecule: &Molecule, atom: AtomId) -> bool {
    let mut fluorine = 0usize;
    let mut chlorine = 0usize;
    let mut bromine = 0usize;
    for (_, bond) in molecule
        .incident_bonds(atom)
        .expect("live atoms have valid adjacency")
    {
        if bond.order != BondOrder::Single {
            continue;
        }
        let Some(neighbor) = molecule.atom(bond.other_atom(atom)) else {
            continue;
        };
        match neighbor.element.atomic_number() {
            FLUORINE => fluorine += 1,
            CHLORINE => chlorine += 1,
            BROMINE => bromine += 1,
            _ => {}
        }
    }
    fluorine >= 3 || chlorine >= 3 || bromine >= 3
}

fn is_methyl_like(
    molecule: &Molecule,
    atom_id: AtomId,
    profiles: &BTreeMap<AtomId, AtomProfile>,
) -> bool {
    let Some(atom) = molecule.atom(atom_id) else {
        return false;
    };
    if atom.element.atomic_number() != CARBON
        || atom.formal_charge != 0
        || atom.radical.is_some()
        || profiles
            .get(&atom_id)
            .is_none_or(|profile| profile.heavy_degree != 1 || profile.incident_triple)
    {
        return false;
    }

    let mut heavy_single_bonds = 0usize;
    let mut graph_hydrogens = 0usize;
    for (_, bond) in molecule
        .incident_bonds(atom_id)
        .expect("live atoms have valid adjacency")
    {
        if !is_ordinary_covalent(bond.order) {
            continue;
        }
        let Some(neighbor) = molecule.atom(bond.other_atom(atom_id)) else {
            continue;
        };
        if neighbor.element.atomic_number() == HYDROGEN {
            if bond.order != BondOrder::Single {
                return false;
            }
            graph_hydrogens += 1;
        } else if bond.order == BondOrder::Single {
            heavy_single_bonds += 1;
        } else {
            return false;
        }
    }
    if heavy_single_bonds != 1 {
        return false;
    }

    match atom.hydrogens {
        ImplicitHydrogens::Inferred => graph_hydrogens <= 3,
        ImplicitHydrogens::Fixed(count) => graph_hydrogens.saturating_add(usize::from(count)) == 3,
    }
}

fn resonance_restricted_atoms(
    molecule: &Molecule,
    ring_membership: &crate::molecule::RingMembership,
    profiles: &BTreeMap<AtomId, AtomProfile>,
    include_ring_bonds: bool,
) -> BTreeSet<AtomId> {
    let mut restricted = BTreeSet::new();
    for center in molecule.atom_ids() {
        let Some(center_atom) = molecule.atom(center) else {
            continue;
        };
        if center_atom.element.atomic_number() != CARBON
            || molecule
                .perception()
                .aromaticity()
                .map(|aromaticity| aromaticity.atom_is_aromatic(center))
                == Some(true)
            || profiles
                .get(&center)
                .is_none_or(|profile| profile.heavy_degree != 3)
        {
            continue;
        }

        let mut has_heteroatom_double_bond = false;
        let mut attachments = Vec::new();
        for (bond_id, bond) in molecule
            .incident_bonds(center)
            .expect("live atoms have valid adjacency")
        {
            let neighbor_id = bond.other_atom(center);
            let Some(neighbor) = molecule.atom(neighbor_id) else {
                continue;
            };
            if bond.order == BondOrder::Double
                && is_nitrogen_oxygen_or_sulfur(neighbor.element)
                && molecule
                    .perception()
                    .aromaticity()
                    .map(|aromaticity| aromaticity.atom_is_aromatic(neighbor_id))
                    == Some(false)
            {
                // The strict query's [N,O,S] tests element and aliphaticity,
                // not charge. Charged imidates and amidines still have the
                // same resonance-restricted linkage.
                has_heteroatom_double_bond = true;
            }
            if bond.order != BondOrder::Single
                || (!include_ring_bonds && ring_membership.bond_in_ring(bond_id))
            {
                continue;
            }
            let atomic_number = neighbor.element.atomic_number();
            let heavy_degree = profiles
                .get(&neighbor_id)
                .map_or(0, |profile| profile.heavy_degree);
            if atomic_number == NITROGEN
                || (molecule
                    .perception()
                    .aromaticity()
                    .map(|aromaticity| aromaticity.atom_is_aromatic(neighbor_id))
                    == Some(false)
                    && (atomic_number == OXYGEN || (atomic_number == SULFUR && heavy_degree != 1)))
            {
                attachments.push(neighbor_id);
            }
        }

        if has_heteroatom_double_bond && !attachments.is_empty() {
            restricted.insert(center);
            restricted.extend(attachments);
        }
    }
    restricted
}

fn is_nitrogen_oxygen_or_sulfur(element: Element) -> bool {
    matches!(element.atomic_number(), NITROGEN | OXYGEN | SULFUR)
}

fn is_ordinary_covalent(order: BondOrder) -> bool {
    !matches!(order, BondOrder::Zero | BondOrder::Dative)
}
