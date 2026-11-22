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
//! acid, an aromatic ring nitrogen, or a cationic nitrogen takes the place
//! of one hydrogen.
//!
//! Residues without hydrogens on their heteroatoms take the template's state
//! with their strong acids deprotonated, as at pH 7, when a dictionary
//! template applies. Polymer residues also take the ammonium of a true N
//! terminus that no metal coordinates, and their least charged variant, such
//! as neutral histidine.
//!
//! A stated charge must fit an allowed valence of its atom. A resonance form
//! that gives it comes first: the charge moves along alternating single and
//! double bonds from an atom of the same sign, as `OD1` -1 of an aspartate
//! or `NE2` +1 of a histidine, or a double bond between a stated cation and
//! a stated anion becomes single, as in `[S+]([O-])` of DMSO. Otherwise the
//! atom keeps the charge with the hydrogens it has, in a residue with
//! hydrogens, or with the fewest that give an allowed valence, in one
//! without. A charge no such count holds is rejected. Residues without
//! hydrogens take a stated zero as no charge stated, since writers fill the
//! field with zero.
use std::collections::BTreeMap;

use crate::algorithms::allowed_valences_with_charge;
use crate::components::{Component, ComponentBond, ComponentBondOrder as TemplateOrder};
use crate::molecule::{BondOrder, Element};

/// What the model shows about one template heavy atom of a residue.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct AtomEvidence {
    /// Model hydrogen sites bonded to it (see
    /// [`super::hydrogen_sites`]).
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

    /// Charges each atom by the hydrogens it lacks, as where stated bond
    /// orders exclude the state the hydrogens fit: the model's hydrogens
    /// are then all there are.
    pub(crate) fn charge_lacking_hydrogens(&mut self) {
        for atom in &mut self.atoms {
            if atom.missing == 0 || atom.stated.is_some() {
                continue;
            }
            let lacking = i8::try_from(atom.missing).unwrap_or(i8::MAX);
            atom.charge = atom.charge.saturating_sub(lacking);
            atom.from_hydrogens = atom.from_hydrogens.saturating_sub(lacking);
            atom.hydrogens -= atom.missing;
            atom.missing = 0;
        }
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

/// Whether metal coordination takes the place of a hydrogen of a template
/// heteroatom: an acid, a thiol, a nitrogen of an aromatic ring, whose
/// variants move the hydrogen, or a cationic nitrogen. A neutral amine or
/// amide keeps its hydrogens.
fn coordinates_in_place_of_hydrogen(component: &Component, index: usize, role: Role) -> bool {
    match role {
        Role::StrongAcid | Role::WeakAcid => true,
        Role::Base => {
            let atom = &component.atoms()[index];
            atom.formal_charge().is_some_and(|charge| charge > 0)
                || (atom.aromatic() && in_ring(component, index))
        }
        Role::Neutral => false,
    }
}

/// The state of the first candidate that fits the evidence best, with the
/// stated charges applied, or `None` when the model's hydrogens contradict
/// every candidate. With hydrogens on heteroatoms the best fit lacks the
/// fewest hydrogens, then changes the fewest charges, then gives the most
/// stated charges by its own hydrogens or a resonance form. Without, it
/// gives the most stated charges so, then charges the fewest atoms for
/// coordination, then, for a polymer residue, charges the fewest atoms.
/// Ties go to the template.
pub(crate) fn choose<'c>(
    candidates: impl IntoIterator<Item = &'c Component>,
    evidence: &ResidueEvidence<'_>,
) -> Option<Protonation<'c>> {
    let candidates = candidates.into_iter().collect::<Vec<_>>();
    let base = *candidates.first()?;
    let with_evidence = has_evidence(base, evidence);
    let mut evidence = evidence.clone();
    if !with_evidence {
        for atom in evidence.heavy.values_mut() {
            atom.stated = atom.stated.filter(|&charge| charge != 0);
        }
    }
    let evidence = &evidence;
    let (candidate, atoms, reordered, assumed) = candidates
        .into_iter()
        .enumerate()
        .filter_map(|(order, candidate)| {
            let (mut atoms, mut reordered, assumed) = evaluate(candidate, evidence, with_evidence)?;
            let missing = atoms.iter().map(|atom| atom.missing).sum::<usize>();
            let changed = atoms
                .iter()
                .filter(|atom| atom.from_hydrogens != 0 || atom.coordinated != 0)
                .count();
            let coordination_charges = atoms
                .iter()
                .filter(|atom| atom.coordinated != 0 && atom.charge != 0)
                .count();
            // Fewest charged atoms, not the smallest net charge: a terminal
            // carboxylate must not pull a histidine into its cation.
            let charged = if evidence.polymer {
                atoms.iter().filter(|atom| atom.charge != 0).count()
            } else {
                0
            };
            apply_stated(
                candidate,
                evidence,
                with_evidence,
                &mut atoms,
                &mut reordered,
            );
            // Stated charges the candidate gives only with other hydrogens,
            // or not at all.
            let unfit = atoms
                .iter()
                .filter(|atom| {
                    matches!(
                        atom.stated,
                        Some(Stated::Kept { .. } | Stated::Rejected { .. })
                    )
                })
                .count();
            let key = if with_evidence {
                (missing, changed, unfit, order)
            } else {
                (unfit, coordination_charges, charged, order)
            };
            Some((key, (candidate, atoms, reordered, assumed)))
        })
        .min_by_key(|(key, _)| *key)
        .map(|(_, chosen)| chosen)?;
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
        let coordinates = coordinates_in_place_of_hydrogen(candidate, index, role);
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
                        let released = released_charge(
                            candidate, index, evidence, charge, lacking, present, atom.links,
                        );
                        state.from_hydrogens = i8::try_from(released - charge).ok()?;
                        charge = released;
                        hydrogens = present;
                    }
                    // A cationic nitrogen with at least its neutral count of
                    // hydrogens, none if that count is zero, is the neutral
                    // base; one with none of several lacks them.
                    Role::Base
                        if base > 0
                            && (present > 0 || lacking <= usize::from(base.unsigned_abs())) =>
                    {
                        let released = lacking.min(usize::from(base.unsigned_abs()));
                        state.from_hydrogens = -i8::try_from(released).ok()?;
                        charge -= count(released);
                        state.missing = lacking - released;
                        hydrogens = expected - released;
                    }
                    _ => state.missing = lacking,
                }
            } else if base != 0 {
                // A template charge the atom's bonds and hydrogens cannot
                // hold, as heme's `O2A` -1 with `H2A`, follows the hydrogens.
                let held = held_charge(candidate, index, evidence, charge, present + links);
                state.from_hydrogens = i8::try_from(held - charge).ok()?;
                charge = held;
            }
        } else {
            if evidence.dictionary {
                match role {
                    Role::StrongAcid if links == 0 && expected > 0 => {
                        charge =
                            released_charge(candidate, index, evidence, charge, expected, 0, 0);
                        hydrogens = 0;
                    }
                    Role::Base
                        if evidence.polymer
                            && links == 0
                            && atom.coordination == 0
                            && is_free_amine(candidate, index) =>
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
        let Some(order) = current_order(reordered, bond) else {
            continue;
        };
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
            let Some(other_order) = current_order(reordered, other_bond) else {
                continue;
            };
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

/// Applies stated charges that differ from the chosen state. Resonance
/// forms come first: a charge moves to the atom stated with it along
/// alternating bonds, then a double bond between a stated cation and a
/// stated anion becomes single. Each move must leave both atoms at an
/// allowed valence with the hydrogens they hold. A charge the atom then
/// holds at an allowed valence is kept, with the hydrogens it has where the
/// residue's hydrogens decide; any other is rejected.
fn apply_stated(
    candidate: &Component,
    evidence: &ResidueEvidence<'_>,
    with_evidence: bool,
    atoms: &mut [AtomState],
    reordered: &mut Vec<([usize; 2], BondOrder)>,
) {
    let context = Context {
        candidate,
        evidence,
        with_evidence,
    };
    while move_charge(&context, atoms, reordered) {}
    while separate_charges(&context, atoms, reordered) {}
    for state in atoms.iter_mut() {
        let Some(stated) = context.stated(state) else {
            continue;
        };
        if stated == state.charge || state.stated.is_some() {
            continue;
        }
        let held = candidate
            .atom(&state.name)
            .and_then(|(index, _)| context.held(index, stated, reordered));
        match held {
            Some(hydrogens) => {
                state.stated = Some(Stated::Kept {
                    derived: state.charge,
                });
                state.charge = stated;
                state.hydrogens = hydrogens;
                state.from_hydrogens = 0;
                state.missing = 0;
            }
            None => state.stated = Some(Stated::Rejected { stated }),
        }
    }
}

/// What the stated charges of one residue apply to.
struct Context<'a, 'c> {
    candidate: &'c Component,
    evidence: &'a ResidueEvidence<'a>,
    with_evidence: bool,
}

impl Context<'_, '_> {
    fn stated(&self, state: &AtomState) -> Option<i8> {
        self.evidence
            .heavy
            .get(state.name.as_str())
            .and_then(|atom| atom.stated)
    }

    /// The state of the template atom `index`, by position.
    fn position(&self, atoms: &[AtomState], index: usize) -> Option<usize> {
        let name = self.candidate.atoms()[index].name();
        atoms.iter().position(|state| state.name == name)
    }

    /// The hydrogens the template atom `index` holds with `charge` and the bond orders `reordered` gives: those
    /// the model has where the residue's hydrogens decide, else the fewest
    /// that give an allowed valence. `None` when no allowed valence holds
    /// it, or a bond of the atom has no localized order.
    fn held(
        &self,
        index: usize,
        charge: i8,
        reordered: &[([usize; 2], BondOrder)],
    ) -> Option<usize> {
        let atom = &self.candidate.atoms()[index];
        let evidence = self.evidence.heavy.get(atom.name())?;
        // A gap is no bond: a stated charge sets the hydrogens of an atom
        // without its neighbour.
        let mut valence = evidence.links;
        for (other, bond) in neighbours(self.candidate, index) {
            let other = &self.candidate.atoms()[other];
            if is_hydrogen(other.element())
                || bond.metal_coordination()
                || !self.evidence.heavy.contains_key(other.name())
            {
                continue;
            }
            valence += match current_order(reordered, bond)? {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Quadruple => 4,
                BondOrder::Zero | BondOrder::Dative => 0,
            };
        }
        let allowed = allowed_valences_with_charge(atom.element().atomic_number(), charge)?;
        // Metals hold any charge.
        if allowed.is_empty() {
            return Some(if self.with_evidence {
                evidence.hydrogens
            } else {
                0
            });
        }
        let allowed = allowed.iter().map(|&valence| usize::from(valence));
        if self.with_evidence {
            let open = open_valence(self.candidate, index, self.evidence);
            let total = valence + evidence.hydrogens;
            let mut allowed = allowed;
            return allowed
                .any(|allowed| total <= allowed && allowed <= total + open)
                .then_some(evidence.hydrogens);
        }
        allowed
            .filter(|&allowed| allowed >= valence)
            .min()
            .map(|allowed| allowed - valence)
    }
}

/// Moves one unit of charge to an atom stated with it from an atom of the
/// same sign along a path of alternating single and double bonds, the
/// shortest first. Returns whether a charge moved.
fn move_charge(
    context: &Context<'_, '_>,
    atoms: &mut [AtomState],
    reordered: &mut Vec<([usize; 2], BondOrder)>,
) -> bool {
    let candidate = context.candidate;
    for position in 0..atoms.len() {
        let Some(stated) = context.stated(&atoms[position]) else {
            continue;
        };
        let own = &atoms[position];
        if own.stated.is_some() || stated == own.charge {
            continue;
        }
        let delta: i8 = if stated > own.charge { 1 } else { -1 };
        let Some((index, _)) = candidate.atom(&own.name) else {
            continue;
        };
        for length in [2, 4, 6] {
            let mut path = vec![index];
            let found = find_path(context, atoms, reordered, &mut path, length, delta);
            if let Some(partner) = found {
                let bonds = path
                    .windows(2)
                    .map(|pair| [pair[0], pair[1]])
                    .collect::<Vec<_>>();
                let Some(swapped) = toggled(candidate, reordered, &bonds) else {
                    continue;
                };
                *reordered = swapped;
                let other = context
                    .position(atoms, partner)
                    .expect("paths end at matched atoms");
                let own_hydrogens = context.held(index, atoms[position].charge + delta, reordered);
                let partner_hydrogens =
                    context.held(partner, atoms[other].charge - delta, reordered);
                for (at, change, hydrogens) in [
                    (position, delta, own_hydrogens),
                    (other, -delta, partner_hydrogens),
                ] {
                    let state = &mut atoms[at];
                    state.charge += change;
                    state.hydrogens = hydrogens.unwrap_or(state.hydrogens);
                    state.missing = 0;
                    state.stated = Some(Stated::Resonance);
                }
                return true;
            }
        }
    }
    false
}

/// Extends `path`, which starts at a template atom, to `length` bonds that
/// alternate between single and double, ending at an atom that gives up
/// `delta` of its charge: one of the same sign whose stated charge, if any,
/// is its charge less `delta`. Both ends must then hold their charges at
/// allowed valences with the hydrogens they had, or, where the residue's
/// hydrogens decide, with the ones the model has. Returns that atom.
fn find_path(
    context: &Context<'_, '_>,
    atoms: &[AtomState],
    reordered: &[([usize; 2], BondOrder)],
    path: &mut Vec<usize>,
    length: usize,
    delta: i8,
) -> Option<usize> {
    let candidate = context.candidate;
    let last = *path.last()?;
    let previous = path.len().checked_sub(2).map(|at| path[at]);
    let previous_order = previous.and_then(|previous| {
        let bond = bond_between(candidate, previous, last)?;
        current_order(reordered, bond)
    });
    for (next, bond) in neighbours(candidate, last) {
        let atom = &candidate.atoms()[next];
        if path.contains(&next)
            || is_hydrogen(atom.element())
            || bond.metal_coordination()
            || !context.evidence.heavy.contains_key(atom.name())
        {
            continue;
        }
        let Some(order) = current_order(reordered, bond) else {
            continue;
        };
        let alternates = match previous_order {
            None => matches!(order, BondOrder::Single | BondOrder::Double),
            Some(BondOrder::Single) => order == BondOrder::Double,
            Some(BondOrder::Double) => order == BondOrder::Single,
            Some(_) => false,
        };
        if !alternates {
            continue;
        }
        path.push(next);
        if path.len() == length + 1 {
            if ends_move(context, atoms, reordered, path, delta) {
                return Some(next);
            }
        } else if let Some(found) = find_path(context, atoms, reordered, path, length, delta) {
            return Some(found);
        }
        path.pop();
    }
    None
}

/// Whether swapping the bonds along `path` moves `delta` of charge from its
/// last atom to its first, each then holding its charge.
fn ends_move(
    context: &Context<'_, '_>,
    atoms: &[AtomState],
    reordered: &[([usize; 2], BondOrder)],
    path: &[usize],
    delta: i8,
) -> bool {
    let (Some(&first), Some(&last)) = (path.first(), path.last()) else {
        return false;
    };
    let (Some(own), Some(partner)) = (
        context.position(atoms, first),
        context.position(atoms, last),
    ) else {
        return false;
    };
    let partner_state = &atoms[partner];
    let charge = partner_state.charge;
    if charge.signum() != delta.signum()
        || partner_state.stated.is_some()
        || context
            .stated(partner_state)
            .is_some_and(|stated| stated != charge - delta)
        || candidate_element(context, first) != candidate_element(context, last)
    {
        return false;
    }
    let bonds = path
        .windows(2)
        .map(|pair| [pair[0], pair[1]])
        .collect::<Vec<_>>();
    let Some(swapped) = toggled(context.candidate, reordered, &bonds) else {
        return false;
    };
    let keeps = |at: usize, index: usize, change: i8| {
        let state = &atoms[at];
        context
            .held(index, state.charge + change, &swapped)
            .is_some_and(|hydrogens| context.with_evidence || hydrogens == state.hydrogens)
    };
    keeps(own, first, delta) && keeps(partner, last, -delta)
}

fn candidate_element(context: &Context<'_, '_>, index: usize) -> Element {
    context.candidate.atoms()[index].element()
}

/// Makes a double bond between an atom stated more positive than its state
/// and one stated more negative single, each taking one unit of charge, as
/// the charge-separated sulfoxide of DMSO. Returns whether a bond changed.
fn separate_charges(
    context: &Context<'_, '_>,
    atoms: &mut [AtomState],
    reordered: &mut Vec<([usize; 2], BondOrder)>,
) -> bool {
    let candidate = context.candidate;
    let wants = |state: &AtomState| {
        if state
            .stated
            .is_some_and(|stated| stated != Stated::Resonance)
        {
            return 0;
        }
        context
            .stated(state)
            .map_or(0, |stated| (stated - state.charge).signum())
    };
    for bond in candidate.bonds() {
        let [a, b] = bond.atoms();
        let (Some(left), Some(right)) = (context.position(atoms, a), context.position(atoms, b))
        else {
            continue;
        };
        if current_order(reordered, bond) != Some(BondOrder::Double) {
            continue;
        }
        let (cation, anion, cation_index, anion_index) =
            match (wants(&atoms[left]), wants(&atoms[right])) {
                (1, -1) => (left, right, a, b),
                (-1, 1) => (right, left, b, a),
                _ => continue,
            };
        let mut swapped = reordered.clone();
        swapped.retain(|(atoms, _)| *atoms != [a, b] && *atoms != [b, a]);
        swapped.push(([a, b], BondOrder::Single));
        let held = |at: usize, index: usize, change: i8| {
            let state = &atoms[at];
            context
                .held(index, state.charge + change, &swapped)
                .filter(|&hydrogens| context.with_evidence || hydrogens == state.hydrogens)
        };
        let (Some(cation_hydrogens), Some(anion_hydrogens)) =
            (held(cation, cation_index, 1), held(anion, anion_index, -1))
        else {
            continue;
        };
        *reordered = swapped;
        for (at, change, hydrogens) in [(cation, 1, cation_hydrogens), (anion, -1, anion_hydrogens)]
        {
            let state = &mut atoms[at];
            state.charge += change;
            state.hydrogens = hydrogens;
            state.missing = 0;
            state.stated = Some(Stated::Resonance);
        }
        return true;
    }
    false
}

/// `reordered` with each of `bonds`, by template atoms, swapped between
/// single and double; `None` when one is neither.
fn toggled(
    candidate: &Component,
    reordered: &[([usize; 2], BondOrder)],
    bonds: &[[usize; 2]],
) -> Option<Vec<([usize; 2], BondOrder)>> {
    let mut swapped = reordered.to_vec();
    for &[a, b] in bonds {
        let bond = bond_between(candidate, a, b)?;
        let order = match current_order(reordered, bond)? {
            BondOrder::Single => BondOrder::Double,
            BondOrder::Double => BondOrder::Single,
            _ => return None,
        };
        swapped.retain(|(atoms, _)| *atoms != [a, b] && *atoms != [b, a]);
        swapped.push((bond.atoms(), order));
    }
    Some(swapped)
}

fn bond_between(candidate: &Component, a: usize, b: usize) -> Option<&ComponentBond> {
    candidate.bonds().iter().find(|bond| {
        let atoms = bond.atoms();
        atoms == [a, b] || atoms == [b, a]
    })
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

/// Whether a template atom bonds to a leaving heavy atom that `present`
/// names, as a C terminus's `C` to `OXT`.
pub(crate) fn bonds_leaving_atom(
    component: &Component,
    index: usize,
    present: &std::collections::BTreeSet<&str>,
) -> bool {
    neighbours(component, index).any(|(other, _)| {
        let atom = &component.atoms()[other];
        atom.leaving() && !is_hydrogen(atom.element()) && present.contains(atom.name())
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

/// The charge of the acid atom `index` with charge `charge` once it lacks
/// `lacking` of its template hydrogens, keeping `left` and the valence
/// `links` to other residues: one less per hydrogen, but not below the
/// charge at which its bonds give an allowed valence. A template that
/// already charges the atom while listing its hydrogen, as the carboxylates
/// of heme, keeps that charge.
fn released_charge(
    candidate: &Component,
    index: usize,
    evidence: &ResidueEvidence<'_>,
    charge: i32,
    lacking: usize,
    left: usize,
    links: usize,
) -> i32 {
    let Some(valence) = heavy_valence(candidate, index, evidence) else {
        return charge - count(lacking);
    };
    (0..=lacking)
        .rev()
        .map(|released| charge - count(released))
        .find(|&charge| holds(candidate, index, charge, valence + links + left))
        .unwrap_or(charge - count(lacking))
}

/// The charge of the template atom `index`, charged `charge`, whose
/// hydrogens and bonds to other residues give the valence `outside`: the
/// charge itself where its valence is allowed, else the one nearest it
/// towards zero that allows it. Unchanged where none does or a bond has no
/// localized order.
fn held_charge(
    candidate: &Component,
    index: usize,
    evidence: &ResidueEvidence<'_>,
    charge: i32,
    outside: usize,
) -> i32 {
    let Some(valence) = heavy_valence(candidate, index, evidence).map(|valence| valence + outside)
    else {
        return charge;
    };
    if holds(candidate, index, charge, valence) {
        return charge;
    }
    let step = -charge.signum();
    (1..=charge.abs())
        .map(|offset| charge + step * offset)
        .find(|&charge| holds(candidate, index, charge, valence))
        .unwrap_or(charge)
}

/// The valence the template bonds of atom `index` to heavy atoms give it:
/// all but metal coordination and bonds to leaving atoms the residue lacks.
/// `None` when a bond has no localized order.
fn heavy_valence(
    candidate: &Component,
    index: usize,
    evidence: &ResidueEvidence<'_>,
) -> Option<usize> {
    let mut valence = 0;
    for (other, bond) in neighbours(candidate, index) {
        let other = &candidate.atoms()[other];
        if is_hydrogen(other.element())
            || bond.metal_coordination()
            || (other.leaving() && !evidence.heavy.contains_key(other.name()))
        {
            continue;
        }
        valence += match current_order(&[], bond)? {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Zero | BondOrder::Dative => 0,
        };
    }
    Some(valence)
}

/// Whether the template atom `index` has `valence` among its allowed
/// valences with `charge`.
fn holds(candidate: &Component, index: usize, charge: i32, valence: usize) -> bool {
    let atomic_number = candidate.atoms()[index].element().atomic_number();
    i8::try_from(charge)
        .ok()
        .and_then(|charge| allowed_valences_with_charge(atomic_number, charge))
        .is_some_and(|allowed| allowed.iter().any(|&value| usize::from(value) == valence))
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
