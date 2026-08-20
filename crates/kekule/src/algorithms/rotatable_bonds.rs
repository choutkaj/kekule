use std::collections::{BTreeMap, BTreeSet};

use crate::core::{AtomId, Bond, BondId, BondOrder, Element, Molecule};

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
    /// Includes bonds whose endpoint has no other covalently connected
    /// non-hydrogen neighbor.
    pub include_terminal_bonds: bool,
    /// Includes resonance-restricted linkages such as amides and esters.
    pub include_resonance_restricted_bonds: bool,
    /// Includes axes involving symmetric groups such as trihalomethyl and
    /// tert-butyl-like centers.
    pub include_symmetric_bonds: bool,
    /// Includes represented single bonds that are members of graph cycles.
    pub include_ring_bonds: bool,
}

impl RotatableBondOptions {
    /// RDKit's strict rotatable-bond descriptor version 3.2.0, adapted to
    /// ignore graph hydrogens when evaluating terminal and symmetric groups.
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

/// Detects rotatable bonds using the supplied options.
pub fn detect_rotatable_bonds(
    molecule: &Molecule,
    options: RotatableBondOptions,
) -> RotatableBondSet {
    let computed_membership;
    let ring_membership = match molecule.ring_membership() {
        Some(membership) => membership,
        None => {
            computed_membership = compute_ring_membership(molecule);
            &computed_membership
        }
    };
    let profiles = atom_profiles(molecule, ring_membership, options.include_ring_bonds);
    let bonds = molecule
        .bonds()
        .filter_map(|(bond_id, bond)| {
            let rotatable = rotatable_bond(bond_id, bond, ring_membership, &profiles, options);
            rotatable.then_some(bond_id)
        })
        .collect();

    RotatableBondSet { options, bonds }
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
    ring_membership: &crate::core::RingMembership,
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
    ring_membership: &crate::core::RingMembership,
    include_ring_bonds: bool,
) -> BTreeMap<AtomId, AtomProfile> {
    let mut profiles = molecule
        .atom_ids()
        .map(|atom| {
            let mut profile = AtomProfile {
                is_hydrogen: molecule
                    .atom(atom)
                    .is_ok_and(|atom| atom.element.atomic_number() == HYDROGEN),
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
                    .is_ok_and(|neighbor| neighbor.element.atomic_number() != HYDROGEN)
                {
                    profile.heavy_degree += 1;
                }
                profile.incident_triple |= bond.order == BondOrder::Triple;
            }
            (atom, profile)
        })
        .collect::<BTreeMap<_, _>>();

    for atom in molecule.atom_ids() {
        let Some(value) = molecule.atom(atom).ok() else {
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
    for atom in resonance_restricted_atoms(molecule, ring_membership, &profiles, include_ring_bonds)
    {
        if let Some(profile) = profiles.get_mut(&atom) {
            profile.resonance_restricted = true;
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
        let Ok(neighbor) = molecule.atom(bond.other_atom(atom)) else {
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
    let Ok(atom) = molecule.atom(atom_id) else {
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
        let Ok(neighbor) = molecule.atom(bond.other_atom(atom_id)) else {
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

    let represented = usize::from(atom.hydrogens.explicit_count());
    if atom.hydrogens.allows_implicit() {
        graph_hydrogens.saturating_add(represented) <= 3
    } else {
        graph_hydrogens.saturating_add(represented) == 3
    }
}

fn resonance_restricted_atoms(
    molecule: &Molecule,
    ring_membership: &crate::core::RingMembership,
    profiles: &BTreeMap<AtomId, AtomProfile>,
    include_ring_bonds: bool,
) -> BTreeSet<AtomId> {
    let mut restricted = BTreeSet::new();
    for center in molecule.atom_ids() {
        let Ok(center_atom) = molecule.atom(center) else {
            continue;
        };
        if center_atom.element.atomic_number() != CARBON
            || profiles
                .get(&center)
                .is_none_or(|profile| profile.heavy_degree != 3)
        {
            continue;
        }

        let mut has_neutral_double_bond = false;
        let mut has_cationic_nitrogen_double_bond = false;
        let mut broad_attachments = Vec::new();
        let mut cationic_attachments = Vec::new();
        for (bond_id, bond) in molecule
            .incident_bonds(center)
            .expect("live atoms have valid adjacency")
        {
            let neighbor_id = bond.other_atom(center);
            let Ok(neighbor) = molecule.atom(neighbor_id) else {
                continue;
            };
            if bond.order == BondOrder::Double
                && is_nitrogen_oxygen_or_sulfur(neighbor.element)
                && !bond_is_in_localized_aromatic_cycle(molecule, bond_id, ring_membership)
            {
                has_neutral_double_bond |= neighbor.formal_charge == 0;
                has_cationic_nitrogen_double_bond |=
                    neighbor.element.atomic_number() == NITROGEN && neighbor.formal_charge > 0;
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
                || atomic_number == OXYGEN
                || (atomic_number == SULFUR && heavy_degree != 1)
            {
                broad_attachments.push(neighbor_id);
            }
            if atomic_number == NITROGEN && heavy_degree != 1 {
                cationic_attachments.push(neighbor_id);
            }
        }

        let mut attachments = Vec::new();
        if has_neutral_double_bond {
            attachments.extend(broad_attachments);
        }
        if has_cationic_nitrogen_double_bond {
            attachments.extend(cationic_attachments);
        }
        if !attachments.is_empty() {
            restricted.insert(center);
            restricted.extend(attachments);
        }
    }
    restricted
}

fn bond_is_in_localized_aromatic_cycle(
    molecule: &Molecule,
    bond_id: BondId,
    ring_membership: &crate::core::RingMembership,
) -> bool {
    // RDKit's SMARTS uses aliphatic `C`, while Kekule keeps canonical bond
    // orders localized and does not require aromaticity perception here.
    // Recognize the common localized five- and six-member conjugated cycles
    // directly so a ring C=N bond is not mistaken for an imide-like motif.
    ring_membership.bond_in_ring(bond_id)
        && bond_is_in_fully_conjugated_five_or_six_member_cycle(molecule, bond_id, ring_membership)
}

fn bond_is_in_fully_conjugated_five_or_six_member_cycle(
    molecule: &Molecule,
    bond_id: BondId,
    ring_membership: &crate::core::RingMembership,
) -> bool {
    let Ok(focus) = molecule.bond(bond_id) else {
        return false;
    };
    let (start, target) = focus.endpoints();
    if !atom_is_conjugation_capable(molecule, start)
        || !atom_is_conjugation_capable(molecule, target)
    {
        return false;
    }
    let mut visited = vec![false; molecule.atoms.len()];
    visited[start.index()] = true;
    conjugated_ring_path_exists(
        molecule,
        start,
        target,
        bond_id,
        ring_membership,
        0,
        &mut visited,
    )
}

fn conjugated_ring_path_exists(
    molecule: &Molecule,
    current: AtomId,
    target: AtomId,
    focus: BondId,
    ring_membership: &crate::core::RingMembership,
    depth: usize,
    visited: &mut [bool],
) -> bool {
    if depth == 5 {
        return false;
    }
    for (bond_id, bond) in molecule
        .incident_bonds(current)
        .expect("live atoms have valid adjacency")
    {
        if bond_id == focus || !ring_membership.bond_in_ring(bond_id) {
            continue;
        }
        let neighbor = bond.other_atom(current);
        let next_depth = depth + 1;
        if neighbor == target {
            if matches!(next_depth + 1, 5 | 6) {
                return true;
            }
            continue;
        }
        if visited[neighbor.index()] || !atom_is_conjugation_capable(molecule, neighbor) {
            continue;
        }
        visited[neighbor.index()] = true;
        if conjugated_ring_path_exists(
            molecule,
            neighbor,
            target,
            focus,
            ring_membership,
            next_depth,
            visited,
        ) {
            return true;
        }
        visited[neighbor.index()] = false;
    }
    false
}

fn atom_is_conjugation_capable(molecule: &Molecule, atom_id: AtomId) -> bool {
    let Ok(atom) = molecule.atom(atom_id) else {
        return false;
    };
    let explicit_pi_bonds = molecule
        .incident_bonds(atom_id)
        .expect("live atoms have valid adjacency")
        .filter(|(_, bond)| {
            matches!(
                bond.order,
                BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
            )
        })
        .count();
    if explicit_pi_bonds > 1 {
        return false;
    }
    if atom.formal_charge != 0 || atom.radical.is_some() {
        return true;
    }
    if matches!(atom.element.atomic_number(), NITROGEN | OXYGEN | SULFUR) {
        return true;
    }
    explicit_pi_bonds == 1
}

fn is_nitrogen_oxygen_or_sulfur(element: Element) -> bool {
    matches!(element.atomic_number(), NITROGEN | OXYGEN | SULFUR)
}

fn is_ordinary_covalent(order: BondOrder) -> bool {
    !matches!(order, BondOrder::Zero | BondOrder::Dative)
}
