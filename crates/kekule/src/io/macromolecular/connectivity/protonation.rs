//! The protonation state of a templated residue: which template or variant
//! its hydrogens, links, and metal coordination fit, and the formal charges
//! that follow. Readers and the mmCIF writer's test of what a reader
//! rebuilds share it.
//!
//! Hydrogens decide where the model has them on heteroatoms. A hydrogen
//! beyond the template's charges its atom by one, unless it fills the place
//! of a template heavy atom the model lacks, as `HO5'` fills a 5' end's
//! absent phosphorus. A hydrogen the model lacks deprotonates an acid
//! (carboxylic, phosphoric, or sulfonic acid, phenol, thiol); on a nitrogen
//! or an alcohol it is missing from the model, and the atom keeps its
//! template charge. A metal coordinating a thiol, an acid, or a nitrogen
//! takes the place of one hydrogen. Residues without hydrogens on their
//! heteroatoms take the template's state or, for amino acids and
//! nucleotides read with a dictionary template, the usual state at pH 7.
use std::collections::BTreeMap;

use crate::components::{Component, ComponentBondOrder as TemplateOrder};
use crate::molecule::Element;

/// What the model shows about one template heavy atom of a residue.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct AtomEvidence {
    /// Model hydrogens bonded to it.
    pub(crate) hydrogens: usize,
    /// The valence its bonds to other residues give it (see
    /// [`super::link_valence`]).
    pub(crate) links: usize,
    /// Metal coordination records between it and other residues.
    pub(crate) coordination: usize,
}

/// The evidence of one residue: its matched heavy atoms by template name,
/// and whether it has model hydrogens.
pub(crate) struct ResidueEvidence<'a> {
    pub(crate) heavy: BTreeMap<&'a str, AtomEvidence>,
    pub(crate) has_hydrogens: bool,
    /// Whether the residue is an amino acid or nucleotide read with a
    /// dictionary template, which takes the usual pH 7 state without
    /// hydrogens.
    pub(crate) usual_state: bool,
}

/// The state one heavy atom takes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct AtomState {
    /// The template atom name.
    pub(crate) name: String,
    /// The formal charge, in e.
    pub(crate) charge: i8,
    /// Whether the template states no charge, read as zero.
    pub(crate) charge_unknown: bool,
    /// The hydrogens the state places on the atom.
    pub(crate) hydrogens: usize,
    /// The part of the charge the model's hydrogens change.
    pub(crate) from_hydrogens: i8,
    /// Hydrogens the model lacks that no charge explains.
    pub(crate) missing: usize,
    /// Hydrogens metal coordination replaces.
    pub(crate) coordinated: usize,
}

/// The chosen template or variant with the state of each matched heavy
/// atom.
pub(crate) struct Protonation<'c> {
    pub(crate) candidate: &'c Component,
    pub(crate) atoms: Vec<AtomState>,
    /// Whether the residue had no hydrogens on its heteroatoms to decide by.
    pub(crate) without_evidence: bool,
}

/// How a heteroatom's hydrogens behave.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Role {
    /// The hydroxyl of a carboxylic, phosphoric, or sulfonic acid.
    StrongAcid,
    /// A phenol or a thiol.
    WeakAcid,
    /// A nitrogen.
    Base,
    /// An alcohol, water, or another heteroatom whose missing hydrogens
    /// are absent from the model.
    Neutral,
}

impl Role {
    const fn deprotonates(self) -> bool {
        matches!(self, Self::StrongAcid | Self::WeakAcid)
    }

    const fn coordinates_in_place_of_hydrogen(self) -> bool {
        matches!(self, Self::StrongAcid | Self::WeakAcid | Self::Base)
    }
}

/// The state of the first candidate that fits the evidence best, or `None`
/// when the model's hydrogens contradict every candidate. With hydrogens on
/// heteroatoms the best fit lacks the fewest hydrogens, then changes the
/// fewest charges; without, it changes the fewest charges for coordination,
/// then has the smallest net charge.
pub(crate) fn choose<'c>(
    candidates: impl IntoIterator<Item = &'c Component>,
    evidence: &ResidueEvidence<'_>,
) -> Option<Protonation<'c>> {
    let candidates = candidates.into_iter().collect::<Vec<_>>();
    let base = *candidates.first()?;
    let with_evidence = has_evidence(base, evidence);
    candidates
        .into_iter()
        .enumerate()
        .filter_map(|(order, candidate)| {
            let atoms = evaluate(candidate, evidence, with_evidence)?;
            let missing = atoms.iter().map(|atom| atom.missing).sum::<usize>();
            let changed = atoms
                .iter()
                .filter(|atom| atom.from_hydrogens != 0 || atom.coordinated != 0)
                .count();
            let coordination_charges = atoms
                .iter()
                .filter(|atom| atom.coordinated != 0 && atom.charge != 0)
                .count();
            let net = atoms
                .iter()
                .map(|atom| i32::from(atom.charge))
                .sum::<i32>()
                .unsigned_abs();
            let key = if with_evidence {
                (missing, changed, order)
            } else {
                (coordination_charges, net as usize, order)
            };
            Some((key, candidate, atoms))
        })
        .min_by_key(|(key, _, _)| *key)
        .map(|(_, candidate, atoms)| Protonation {
            candidate,
            atoms,
            without_evidence: !with_evidence,
        })
}

/// Whether the model has hydrogens to decide protonation by: hydrogens on
/// some heteroatom, or hydrogens where the template places none on
/// heteroatoms, as on an in-chain proline.
fn has_evidence(base: &Component, evidence: &ResidueEvidence<'_>) -> bool {
    if !evidence.has_hydrogens {
        return false;
    }
    let mut expects_any = false;
    for (&name, atom) in &evidence.heavy {
        let Some((index, template_atom)) = base.atom(name) else {
            continue;
        };
        if is_carbon(template_atom.element()) {
            continue;
        }
        if atom.hydrogens > 0 {
            return true;
        }
        expects_any |= expected_hydrogens(base, index, atom.links) > 0;
    }
    !expects_any
}

/// The state `candidate` gives each matched heavy atom, or `None` when the
/// model contradicts it.
fn evaluate(
    candidate: &Component,
    evidence: &ResidueEvidence<'_>,
    with_evidence: bool,
) -> Option<Vec<AtomState>> {
    // Variants only remove or move hydrogens: a model hydrogen a candidate
    // lacks counts on its heavy atom.
    let mut atoms = Vec::with_capacity(evidence.heavy.len());
    for (&name, atom) in &evidence.heavy {
        let (index, template_atom) = candidate.atom(name)?;
        let template_charge = template_atom.formal_charge();
        let base = template_charge.unwrap_or(0);
        let element = template_atom.element();
        let mut state = AtomState {
            name: template_atom.name().to_owned(),
            charge: base,
            charge_unknown: template_charge.is_none(),
            hydrogens: 0,
            from_hydrogens: 0,
            missing: 0,
            coordinated: 0,
        };
        if is_carbon(element) || is_hydrogen(element) {
            state.hydrogens = expected_hydrogens(candidate, index, atom.links);
            atoms.push(state);
            continue;
        }
        let role = role(candidate, index);
        let mut expected = expected_hydrogens(candidate, index, atom.links);
        if role.coordinates_in_place_of_hydrogen() {
            state.coordinated = atom.coordination.min(expected);
            expected -= state.coordinated;
        }
        let mut charge = i32::from(base) - count(state.coordinated);
        let mut hydrogens = expected;
        if with_evidence {
            let present = atom.hydrogens;
            if present > expected {
                // Hydrogens first fill the places of absent template atoms.
                let open = open_valence(candidate, index, evidence);
                let extra = present - expected - open.min(present - expected);
                if extra > 0 {
                    let ring_nitrogen =
                        element.atomic_number() == 7 && extra == 1 && in_ring(candidate, index);
                    if !only_single_bonds(candidate, index) && !ring_nitrogen {
                        return None;
                    }
                    state.from_hydrogens = i8::try_from(extra).ok()?;
                    charge += count(extra);
                }
                hydrogens = present;
            } else if present < expected {
                let lacking = expected - present;
                if role.deprotonates() {
                    state.from_hydrogens = -i8::try_from(lacking).ok()?;
                    charge -= count(lacking);
                    hydrogens = present;
                } else {
                    state.missing = lacking;
                }
            }
        } else if evidence.usual_state {
            match role {
                Role::StrongAcid if atom.links == 0 && expected > 0 => {
                    charge -= count(expected);
                    hydrogens = 0;
                }
                Role::Base if is_free_amine(candidate, index, atom.links) => {
                    charge += 1;
                    hydrogens += 1;
                }
                _ => {}
            }
        }
        state.charge = i8::try_from(charge).ok()?;
        state.hydrogens = hydrogens;
        atoms.push(state);
    }
    Some(atoms)
}

fn count(value: usize) -> i32 {
    i32::try_from(value).unwrap_or(i32::MAX)
}

/// Whether a template nitrogen is the free amine of an N-terminal amino
/// acid: the `N` of a peptide-linking component, with a hydrogen that
/// leaves when it links, single bonds only, no charge, and no link.
fn is_free_amine(component: &Component, index: usize, links: usize) -> bool {
    let atom = &component.atoms()[index];
    links == 0
        && atom.name() == "N"
        && atom.formal_charge() == Some(0)
        && component.linking() == crate::components::ComponentLinking::Peptide
        && only_single_bonds(component, index)
        && neighbours(component, index).any(|(other, _)| {
            let other = &component.atoms()[other];
            is_hydrogen(other.element()) && other.leaving()
        })
}

/// The valence of a template atom's bonds to template heavy atoms that the
/// residue lacks and that do not leave when it links: places a hydrogen
/// may fill.
fn open_valence(component: &Component, index: usize, evidence: &ResidueEvidence<'_>) -> usize {
    neighbours(component, index)
        .filter(|&(other, bond)| {
            let atom = &component.atoms()[other];
            !is_hydrogen(atom.element())
                && !atom.leaving()
                && !bond.metal_coordination()
                && !evidence.heavy.contains_key(atom.name())
        })
        .map(|(_, bond)| match bond.order() {
            TemplateOrder::Double => 2,
            TemplateOrder::Triple => 3,
            _ => 1,
        })
        .sum()
}

fn role(component: &Component, index: usize) -> Role {
    let atom = &component.atoms()[index];
    let number = atom.element().atomic_number();
    if number == 7 {
        return Role::Base;
    }
    // O, S, Se, and Te.
    if !matches!(number, 8 | 16 | 34 | 52) {
        return Role::Neutral;
    }
    let heavy = neighbours(component, index)
        .filter(|(other, bond)| {
            !is_hydrogen(component.atoms()[*other].element()) && !bond.metal_coordination()
        })
        .collect::<Vec<_>>();
    let [(other, bond)] = heavy[..] else {
        return Role::Neutral;
    };
    if bond.order() != TemplateOrder::Single || bond.aromatic() {
        return Role::Neutral;
    }
    let neighbour = &component.atoms()[other];
    match neighbour.element().atomic_number() {
        6 => {
            let carbonyl = neighbours(component, other).any(|(next, bond)| {
                next != index
                    && matches!(bond.order(), TemplateOrder::Double)
                    && !bond.aromatic()
                    && matches!(
                        component.atoms()[next].element().atomic_number(),
                        8 | 16 | 34
                    )
            });
            let aromatic = neighbour.aromatic()
                || neighbours(component, other).any(|(_, bond)| bond.aromatic());
            if carbonyl {
                Role::StrongAcid
            } else if aromatic || number != 8 {
                Role::WeakAcid
            } else {
                Role::Neutral
            }
        }
        // Phosphoric, arsenic, sulfuric, and selenic acids.
        15 | 33 | 16 | 34 => Role::StrongAcid,
        _ => Role::Neutral,
    }
}

fn neighbours(
    component: &Component,
    index: usize,
) -> impl Iterator<Item = (usize, &crate::components::ComponentBond)> {
    component.bonds().iter().filter_map(move |bond| {
        let [a, b] = bond.atoms();
        match (a == index, b == index) {
            (true, _) => Some((b, bond)),
            (_, true) => Some((a, bond)),
            _ => None,
        }
    })
}

fn only_single_bonds(component: &Component, index: usize) -> bool {
    neighbours(component, index).all(|(_, bond)| {
        bond.order() == TemplateOrder::Single && !bond.aromatic() && !bond.metal_coordination()
    })
}

/// The hydrogens a template atom carries when bonds to other residues give
/// it `links` valence. Each unit takes the place of one hydrogen, and an
/// atom that links loses all its leaving hydrogens: a peptide bond replaces
/// the amine's leaving `H2`, a disulfide replaces cysteine's `HG`, which the
/// dictionary does not mark leaving, and a Schiff base's double bond
/// replaces two of lysine's `NZ` hydrogens.
pub(crate) fn expected_hydrogens(component: &Component, index: usize, links: usize) -> usize {
    let (mut hydrogens, mut leaving) = (0, 0);
    for (other, _) in neighbours(component, index) {
        let atom = &component.atoms()[other];
        if is_hydrogen(atom.element()) {
            hydrogens += 1;
            leaving += usize::from(atom.leaving());
        }
    }
    if links == 0 {
        hydrogens
    } else {
        hydrogens.saturating_sub(leaving.max(links))
    }
}

/// Whether a template atom lies on a ring: some neighbour reaches it again
/// without the bond between them.
fn in_ring(component: &Component, index: usize) -> bool {
    let next = |atom: usize| neighbours(component, atom).map(|(other, _)| other);
    next(index).any(|start| {
        let mut seen = std::collections::BTreeSet::from([index, start]);
        let mut stack = vec![start];
        while let Some(atom) = stack.pop() {
            for other in next(atom) {
                if other == index && atom != start {
                    return true;
                }
                if seen.insert(other) {
                    stack.push(other);
                }
            }
        }
        false
    })
}

fn is_hydrogen(element: Element) -> bool {
    element.atomic_number() == 1
}

fn is_carbon(element: Element) -> bool {
    element.atomic_number() == 6
}
