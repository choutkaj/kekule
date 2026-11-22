//! The protonation state of a templated residue: which template or variant
//! its hydrogens, links, metal coordination, and stated charges fit, and the
//! formal charges that follow. Readers and the writers' tests of what a
//! reader rebuilds share it.
//!
//! Hydrogens decide where the model has them on heteroatoms. A hydrogen
//! beyond the template's charges its atom by one, unless it fills the place
//! of a template heavy atom the model lacks, as `HO5'` fills a 5' end's
//! absent phosphorus. A hydrogen the model lacks deprotonates a strong acid
//! (carboxylic, phosphoric, phosphonic, or sulfonic acid, sulfate). A
//! cationic nitrogen with at least its neutral count of hydrogens is the
//! neutral base. Elsewhere a lacking hydrogen is missing from the model, and
//! the atom keeps its template charge. A metal coordinating a thiol, an
//! acid, or a nitrogen takes the place of one hydrogen.
//!
//! Residues without hydrogens on their heteroatoms take the template's state
//! with their strong acids deprotonated, as at pH 7, when a dictionary
//! template applies. Polymer residues also take the ammonium of a true N
//! terminus and their least charged variant, such as neutral histidine.
//!
//! A stated charge is kept when the atom can hold it with some hydrogen
//! count, which it then sets. One that a resonance form of the atom's group
//! gives, such as `OD1` -1 of an aspartate, moves the double bond instead.
use std::collections::BTreeMap;

use crate::components::{Component, ComponentBond, ComponentBondOrder as TemplateOrder};
use crate::molecule::{BondOrder, Element};

/// What the model shows about one template heavy atom of a residue.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct AtomEvidence {
    /// Model hydrogen sites bonded to it: a hydrogen and a deuterium that
    /// share one exchangeable site count once.
    pub(crate) hydrogens: usize,
    /// The valence its bonds to other residues give it (see
    /// [`super::link_valence`]).
    pub(crate) links: usize,
    /// Metal coordination records between it and other residues.
    pub(crate) coordination: usize,
    /// Whether it bonds to a residue its polymer lacks, across a gap in the
    /// model.
    pub(crate) gap: bool,
    /// The charge the source states.
    pub(crate) stated: Option<i8>,
}

/// The evidence of one residue: its matched heavy atoms by template name,
/// and what is known about the residue as a whole.
#[derive(Debug, Clone)]
pub(crate) struct ResidueEvidence<'a> {
    pub(crate) heavy: BTreeMap<&'a str, AtomEvidence>,
    pub(crate) has_hydrogens: bool,
    /// Whether a dictionary template applies, not a force field's, whose
    /// residue names state their protonation.
    pub(crate) dictionary: bool,
    /// Whether the residue is a monomer of a polymer chain.
    pub(crate) polymer: bool,
}

/// What became of a charge the source states.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Stated {
    /// Kept, with the hydrogens it needs; the residue's state gives
    /// `derived`.
    Kept { derived: i8 },
    /// No hydrogen count holds it; the derived charge replaces it.
    Rejected { stated: i8 },
    /// A resonance form of the atom's group gives it.
    Resonance,
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
    /// Whether it bonds to another residue.
    pub(crate) linked: bool,
    /// Whether it bonds to a residue its polymer lacks, across a gap.
    pub(crate) gap: bool,
    /// Whether readers fix the atom's hydrogens at [`Self::hydrogens`]: it
    /// bonds across a gap, or anchors a leaving atom the residue lacks.
    pub(crate) fixed: bool,
    /// What became of the source's charge, when it differs from the
    /// residue's state.
    pub(crate) stated: Option<Stated>,
}

/// The chosen template or variant with the state of each matched heavy
/// atom.
pub(crate) struct Protonation<'c> {
    pub(crate) candidate: &'c Component,
    pub(crate) atoms: Vec<AtomState>,
    /// Whether the residue had no hydrogens on its heteroatoms to decide by.
    pub(crate) without_evidence: bool,
    /// Whether the residue took an assumed state: it has acids or bases,
    /// no hydrogens on its heteroatoms, and charges the source leaves
    /// unstated.
    pub(crate) assumed: bool,
    /// Candidate bonds, by template atom indices, whose order the state
    /// moves to another resonance form.
    reordered: Vec<([usize; 2], BondOrder)>,
}

impl Protonation<'_> {
    /// The order the state gives a bond of the candidate, or `None` for an
    /// aromatic bond that readers localize.
    pub(crate) fn bond_order(&self, bond: &ComponentBond) -> Option<BondOrder> {
        let [a, b] = bond.atoms();
        self.reordered
            .iter()
            .find(|(atoms, _)| *atoms == [a, b] || *atoms == [b, a])
            .map(|&(_, order)| order)
            .or_else(|| super::template_bond_order(bond))
    }

    /// Whether the state changes any bond order of the candidate.
    pub(crate) fn reorders(&self) -> bool {
        !self.reordered.is_empty()
    }
}

/// How a heteroatom's hydrogens behave.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Role {
    /// The hydroxyl of a carboxylic, phosphoric, phosphonic, or sulfonic
    /// acid, or of a sulfate.
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
    const fn coordinates_in_place_of_hydrogen(self) -> bool {
        matches!(self, Self::StrongAcid | Self::WeakAcid | Self::Base)
    }
}

/// The state of the first candidate that fits the evidence best, or `None`
/// when the model's hydrogens contradict every candidate. With hydrogens on
/// heteroatoms the best fit lacks the fewest hydrogens, then changes the
/// fewest charges; without, it charges the fewest atoms for coordination,
/// then, for a polymer residue, has the smallest net charge. Ties go to
/// the template. Stated charges then apply to the chosen state.
pub(crate) fn choose<'c>(
    candidates: impl IntoIterator<Item = &'c Component>,
    evidence: &ResidueEvidence<'_>,
) -> Option<Protonation<'c>> {
    let candidates = candidates.into_iter().collect::<Vec<_>>();
    let base = *candidates.first()?;
    let with_evidence = has_evidence(base, evidence);
    let (candidate, mut atoms, mut reordered, assumed) = candidates
        .into_iter()
        .enumerate()
        .filter_map(|(order, candidate)| {
            let (atoms, reordered, assumed) = evaluate(candidate, evidence, with_evidence)?;
            let missing = atoms.iter().map(|atom| atom.missing).sum::<usize>();
            let changed = atoms
                .iter()
                .filter(|atom| atom.from_hydrogens != 0 || atom.coordinated != 0)
                .count();
            let coordination_charges = atoms
                .iter()
                .filter(|atom| atom.coordinated != 0 && atom.charge != 0)
                .count();
            let net = if evidence.polymer {
                atoms
                    .iter()
                    .map(|atom| i32::from(atom.charge))
                    .sum::<i32>()
                    .unsigned_abs() as usize
            } else {
                0
            };
            let key = if with_evidence {
                (missing, changed, order)
            } else {
                (coordination_charges, net, order)
            };
            Some((key, (candidate, atoms, reordered, assumed)))
        })
        .min_by_key(|(key, _)| *key)
        .map(|(_, chosen)| chosen)?;
    apply_stated(
        candidate,
        evidence,
        with_evidence,
        &mut atoms,
        &mut reordered,
    );
    Some(Protonation {
        candidate,
        atoms,
        without_evidence: !with_evidence,
        assumed,
        reordered,
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
        expects_any |= expected_hydrogens(base, index, effective_links(atom)) > 0;
    }
    !expects_any
}

/// The links an atom has, counting a bond across a gap as one.
fn effective_links(atom: &AtomEvidence) -> usize {
    atom.links.max(usize::from(atom.gap))
}

type Evaluation = (Vec<AtomState>, Vec<([usize; 2], BondOrder)>, bool);

/// The state `candidate` gives each matched heavy atom, the bonds it
/// reorders, and whether it assumes a state; `None` when the model
/// contradicts the candidate.
fn evaluate(
    candidate: &Component,
    evidence: &ResidueEvidence<'_>,
    with_evidence: bool,
) -> Option<Evaluation> {
    // Variants only remove or move hydrogens: a model hydrogen a candidate
    // lacks counts on its heavy atom.
    let mut atoms = Vec::with_capacity(evidence.heavy.len());
    let mut assumed = false;
    for (&name, atom) in &evidence.heavy {
        let (index, template_atom) = candidate.atom(name)?;
        let template_charge = template_atom.formal_charge();
        let base = template_charge.unwrap_or(0);
        let element = template_atom.element();
        let links = effective_links(atom);
        let expected = expected_hydrogens(candidate, index, links);
        let fixed = atom.gap || (atom.links == 0 && anchors_missing(candidate, index, evidence));
        let mut state = AtomState {
            name: template_atom.name().to_owned(),
            charge: base,
            charge_unknown: template_charge.is_none(),
            hydrogens: expected,
            from_hydrogens: 0,
            missing: 0,
            coordinated: 0,
            linked: atom.links > 0,
            gap: atom.gap,
            fixed,
            stated: None,
        };
        if is_carbon(element) || is_hydrogen(element) {
            atoms.push(state);
            continue;
        }
        let role = role(candidate, index);
        let coordinates = role.coordinates_in_place_of_hydrogen();
        let mut charge = i32::from(base);
        let mut hydrogens = expected;
        if with_evidence {
            let mut expected = expected;
            if coordinates {
                state.coordinated = atom.coordination.min(expected);
                expected -= state.coordinated;
                charge -= count(state.coordinated);
            }
            hydrogens = expected;
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
                match role {
                    Role::StrongAcid => {
                        state.from_hydrogens = -i8::try_from(lacking).ok()?;
                        charge -= count(lacking);
                        hydrogens = present;
                    }
                    // A cationic nitrogen with at least its neutral count of
                    // hydrogens is the neutral base; one with none of them
                    // lacks them.
                    Role::Base if base > 0 && present > 0 => {
                        let released = lacking.min(usize::from(base.unsigned_abs()));
                        state.from_hydrogens = -i8::try_from(released).ok()?;
                        charge -= count(released);
                        state.missing = lacking - released;
                        hydrogens = expected - released;
                    }
                    _ => state.missing = lacking,
                }
            }
        } else {
            if evidence.dictionary {
                match role {
                    Role::StrongAcid if links == 0 && expected > 0 => {
                        charge -= count(expected);
                        hydrogens = 0;
                    }
                    Role::Base
                        if evidence.polymer && links == 0 && is_free_amine(candidate, index) =>
                    {
                        charge += 1;
                        hydrogens += 1;
                    }
                    _ => {}
                }
                assumed |= role != Role::Neutral && atom.stated.is_none();
            }
            // Coordination takes the place of a hydrogen the state leaves.
            if coordinates {
                state.coordinated = atom.coordination.min(hydrogens);
                hydrogens -= state.coordinated;
                charge -= count(state.coordinated);
            }
        }
        state.charge = i8::try_from(charge).ok()?;
        state.hydrogens = hydrogens;
        atoms.push(state);
    }
    let mut reordered = Vec::new();
    if with_evidence {
        release_resonance_partners(candidate, &mut atoms, &mut reordered);
    }
    Some((atoms, reordered, assumed))
}

/// A nitrogen that lacks a hydrogen next to a cationic nitrogen of the same
/// group, such as arginine's `NH1` beside the iminium `NH2`, takes the
/// double bond: the group is the neutral base.
fn release_resonance_partners(
    candidate: &Component,
    atoms: &mut [AtomState],
    reordered: &mut Vec<([usize; 2], BondOrder)>,
) {
    for lacking in 0..atoms.len() {
        let state = &atoms[lacking];
        let Some((index, atom)) = candidate.atom(&state.name) else {
            continue;
        };
        if atom.element().atomic_number() != 7
            || state.missing == 0
            || state.charge != 0
            || state.hydrogens == state.missing
        {
            continue;
        }
        let Some((partner, pair)) =
            resonance_partner(candidate, atoms, reordered, index, |other| {
                other.charge > 0 && other.missing == 0 && other.from_hydrogens == 0
            })
        else {
            continue;
        };
        set_swapped(candidate, reordered, pair);
        atoms[partner].charge -= 1;
        atoms[partner].from_hydrogens = -1;
        atoms[lacking].missing -= 1;
        atoms[lacking].hydrogens -= 1;
    }
}

/// The state of an atom of the same element as the template atom `index`,
/// bonded to the same centre, that `accepts`, with the two bonds to the
/// centre, one single and one double, that a resonance form swaps.
fn resonance_partner(
    candidate: &Component,
    atoms: &[AtomState],
    reordered: &[([usize; 2], BondOrder)],
    index: usize,
    accepts: impl Fn(&AtomState) -> bool,
) -> Option<(usize, [[usize; 2]; 2])> {
    let element = candidate.atoms()[index].element();
    let state_of = |template: usize| {
        let name = candidate.atoms()[template].name();
        atoms.iter().position(|state| state.name == name)
    };
    let own = state_of(index)?;
    for (centre, bond) in neighbours(candidate, index) {
        if bond.aromatic() || is_hydrogen(candidate.atoms()[centre].element()) {
            continue;
        }
        let order = current_order(reordered, bond)?;
        for (other, other_bond) in neighbours(candidate, centre) {
            if other == index
                || other_bond.aromatic()
                || candidate.atoms()[other].element() != element
            {
                continue;
            }
            let Some(partner) = state_of(other) else {
                continue;
            };
            let other_order = current_order(reordered, other_bond)?;
            let swaps = matches!(
                (order, other_order),
                (BondOrder::Single, BondOrder::Double) | (BondOrder::Double, BondOrder::Single)
            );
            if swaps && atoms[partner].hydrogens == atoms[own].hydrogens && accepts(&atoms[partner])
            {
                return Some((partner, [bond.atoms(), other_bond.atoms()]));
            }
        }
    }
    None
}

fn current_order(reordered: &[([usize; 2], BondOrder)], bond: &ComponentBond) -> Option<BondOrder> {
    let [a, b] = bond.atoms();
    reordered
        .iter()
        .find(|(atoms, _)| *atoms == [a, b] || *atoms == [b, a])
        .map(|&(_, order)| order)
        .or_else(|| super::template_bond_order(bond))
}

/// Applies stated charges that differ from the chosen state: a resonance
/// form that gives the charge moves the double bond; a charge some
/// hydrogen count holds is kept with that count; any other is rejected.
/// Resonance forms come first, so that a partner's stated zero does not
/// fix its hydrogens before the form moves.
fn apply_stated(
    candidate: &Component,
    evidence: &ResidueEvidence<'_>,
    with_evidence: bool,
    atoms: &mut [AtomState],
    reordered: &mut Vec<([usize; 2], BondOrder)>,
) {
    let stated_of = |state: &AtomState| {
        evidence
            .heavy
            .get(state.name.as_str())
            .and_then(|atom| atom.stated)
    };
    for position in 0..atoms.len() {
        let Some(stated) = stated_of(&atoms[position]) else {
            continue;
        };
        if stated == 0 || atoms[position].charge != 0 || atoms[position].stated.is_some() {
            continue;
        }
        let Some((index, _)) = candidate.atom(&atoms[position].name) else {
            continue;
        };
        let partner = resonance_partner(candidate, atoms, reordered, index, |other| {
            other.charge == stated
                && other.stated.is_none()
                && stated_of(other).is_none_or(|charge| charge == 0)
        });
        if let Some((partner, pair)) = partner {
            set_swapped(candidate, reordered, pair);
            atoms[partner].charge = 0;
            atoms[partner].stated = Some(Stated::Resonance);
            atoms[position].charge = stated;
            atoms[position].stated = Some(Stated::Resonance);
        }
    }
    for state in atoms.iter_mut() {
        let Some(atom) = evidence.heavy.get(state.name.as_str()) else {
            continue;
        };
        let Some(stated) = atom.stated else {
            continue;
        };
        if stated == state.charge || state.stated.is_some() {
            continue;
        }
        let Some((_, template_atom)) = candidate.atom(&state.name) else {
            continue;
        };
        let held = stated_hydrogens(
            template_atom.element(),
            state.charge,
            state.hydrogens,
            stated,
        )
        .filter(|&hydrogens| !with_evidence || hydrogens >= atom.hydrogens);
        match held {
            Some(hydrogens) => {
                state.stated = Some(Stated::Kept {
                    derived: state.charge,
                });
                state.charge = stated;
                state.hydrogens = hydrogens;
                state.from_hydrogens = 0;
                state.missing = if with_evidence {
                    hydrogens - atom.hydrogens
                } else {
                    0
                };
            }
            None => state.stated = Some(Stated::Rejected { stated }),
        }
    }
}

/// Records a resonance swap of two bonds to one centre, one single and one
/// double, by their template atoms.
fn set_swapped(
    candidate: &Component,
    reordered: &mut Vec<([usize; 2], BondOrder)>,
    pair: [[usize; 2]; 2],
) {
    let order_of = |reordered: &[([usize; 2], BondOrder)], atoms: [usize; 2]| {
        let bond = candidate.bonds().iter().find(|bond| {
            let [a, b] = bond.atoms();
            [a, b] == atoms || [b, a] == atoms
        })?;
        current_order(reordered, bond)
    };
    let (Some(first), Some(second)) = (order_of(reordered, pair[0]), order_of(reordered, pair[1]))
    else {
        return;
    };
    reordered.retain(|(atoms, _)| {
        let same = |pair: [usize; 2]| *atoms == pair || *atoms == [pair[1], pair[0]];
        !same(pair[0]) && !same(pair[1])
    });
    reordered.push((pair[0], second));
    reordered.push((pair[1], first));
}

/// The hydrogens a heteroatom or carbon holds with charge `stated` in place
/// of `charge` with `hydrogens`, or `None` when no hydrogen count holds it.
/// A proton adds one charge to N, O, S, Se, or Te; a carbon gains a charge
/// either way by losing one. An atom without hydrogens, such as a metal,
/// holds any charge.
fn stated_hydrogens(element: Element, charge: i8, hydrogens: usize, stated: i8) -> Option<usize> {
    let hydrogens = i32::try_from(hydrogens).ok()?;
    let change = i32::from(stated) - i32::from(charge);
    let held = match element.atomic_number() {
        7 | 8 | 16 | 34 | 52 => hydrogens + change,
        6 => hydrogens - change.abs(),
        _ if hydrogens == 0 => 0,
        _ => return None,
    };
    usize::try_from(held).ok()
}

fn count(value: usize) -> i32 {
    i32::try_from(value).unwrap_or(i32::MAX)
}

/// Whether a template nitrogen is the free amine of an N-terminal amino
/// acid: the `N` of a peptide-linking component, with a hydrogen that
/// leaves when it links, single bonds only, and no charge.
fn is_free_amine(component: &Component, index: usize) -> bool {
    let atom = &component.atoms()[index];
    atom.name() == "N"
        && atom.formal_charge() == Some(0)
        && component.linking() == crate::components::ComponentLinking::Peptide
        && only_single_bonds(component, index)
        && neighbours(component, index).any(|(other, _)| {
            let other = &component.atoms()[other];
            is_hydrogen(other.element()) && other.leaving()
        })
}

/// Whether a template atom anchors a leaving heavy atom the residue lacks,
/// as a C terminus without `OXT`.
fn anchors_missing(component: &Component, index: usize, evidence: &ResidueEvidence<'_>) -> bool {
    neighbours(component, index).any(|(other, _)| {
        let atom = &component.atoms()[other];
        atom.leaving() && !is_hydrogen(atom.element()) && !evidence.heavy.contains_key(atom.name())
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
        // Phosphoric, phosphonic, arsenic, sulfonic, sulfuric, and selenic
        // acids.
        15 | 33 | 16 | 34 => Role::StrongAcid,
        _ => Role::Neutral,
    }
}

fn neighbours(
    component: &Component,
    index: usize,
) -> impl Iterator<Item = (usize, &ComponentBond)> {
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
