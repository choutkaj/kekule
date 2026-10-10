//! Complete, bounded individualization/refinement for SMILES traversal ties.
//!
//! Color refinement alone is not a canonical labeling: inequivalent vertices
//! can retain the same color. Every unresolved choice is explored, except for
//! exact twins whose transposition is a proved automorphism. The least complete
//! colored adjacency and stereo certificate selects an order. Equal certificates describe
//! the same labeled graph, so choosing either cannot affect emitted SMILES.

use super::*;
use std::result::Result;

const MAX_SEARCH_STATES: usize = 100_000;
const MAX_REFINEMENT_WORK: usize = 50_000_000;
const MAX_PENDING_ATOMS: usize = 2_000_000;

pub(super) struct CanonicalOrder {
    // Writer priorities are ordered initial classes, not chemical symmetries.
    initial_classes: Vec<u32>,
    labels: Vec<usize>,
}

impl CanonicalOrder {
    pub(super) fn new(
        molecule: &Molecule,
        style: CanonicalAtomStyle,
    ) -> Result<Self, MolWriteError> {
        Self::with_limits(molecule, style, MAX_SEARCH_STATES, MAX_REFINEMENT_WORK)
    }

    fn with_limits(
        molecule: &Molecule,
        style: CanonicalAtomStyle,
        max_states: usize,
        max_work: usize,
    ) -> Result<Self, MolWriteError> {
        let atoms = molecule.atom_ids().collect::<Vec<_>>();
        let slots = molecule.graph.atom_slot_count();
        let mut dense = vec![usize::MAX; slots];
        for (index, atom) in atoms.iter().enumerate() {
            dense[atom.index()] = index;
        }
        let keys = atoms
            .iter()
            .map(|atom| {
                (
                    writer_atom_priority(molecule, *atom),
                    molecule
                        .atom(*atom)
                        .expect("canonical atom is live")
                        .radical
                        .map(|radical| (radical.electron_count(), radical.spin_multiplicity())),
                    canonical_smiles_atom_for_sort(molecule, *atom, style),
                )
            })
            .collect::<Vec<_>>();
        let mut sorted_keys = keys.clone();
        sorted_keys.sort_unstable();
        sorted_keys.dedup();
        let initial = keys
            .iter()
            .map(|key| {
                sorted_keys
                    .binary_search(key)
                    .expect("every atom key is indexed")
            })
            .collect::<Vec<_>>();
        let mut initial_classes = vec![u32::MAX; slots];
        for (atom, color) in atoms.iter().zip(&initial) {
            initial_classes[atom.index()] = *color as u32;
        }
        let adjacency = atoms
            .iter()
            .map(|atom| {
                let mut neighbors = smiles_incident_bonds_for_style(molecule, *atom, style)?
                    .into_iter()
                    .map(|(_, order, other)| (dense[other.index()], bond_order_code(order)))
                    .collect::<Vec<_>>();
                neighbors.sort_unstable();
                Ok(neighbors)
            })
            .collect::<Result<Vec<_>, MolWriteError>>()?;
        let element_indices: BTreeMap<_, _> = molecule
            .stereo_elements()
            .enumerate()
            .map(|(i, (id, _))| (id, i))
            .collect();
        let mut groups: Vec<_> = molecule
            .stereo_groups()
            .filter(|(_, group)| group.kind != StereoGroupKind::Absolute)
            .map(|(_, group)| {
                (
                    group.kind,
                    group.members.iter().map(|id| element_indices[id]).collect(),
                )
            })
            .collect();
        // Match the emitted projection: CX has one absolute set, and `r`
        // requires explicit shielding of otherwise ungrouped tetrahedra.
        let relative = groups
            .iter()
            .any(|(kind, _)| *kind == StereoGroupKind::Relative);
        let absolute: Vec<_> = molecule
            .stereo_elements()
            .filter_map(|(id, element)| {
                let is_absolute = element.group.is_some_and(|g| {
                    molecule.stereo_group(g).unwrap().kind == StereoGroupKind::Absolute
                });
                (is_absolute
                    || (relative
                        && element.group.is_none()
                        && matches!(&element.kind, StereoElementKind::Tetrahedral(_))))
                .then_some(element_indices[&id])
            })
            .collect();
        if !absolute.is_empty() {
            groups.push((StereoGroupKind::Absolute, absolute));
        }
        let stereo = molecule
            .stereo_elements()
            .map(|(_, element)| {
                let carrier = |c| match c {
                    StereoCarrier::Atom(a) => dense[a.index()],
                    StereoCarrier::ImplicitHydrogen => atoms.len(),
                    StereoCarrier::ImplicitLonePair => atoms.len() + 1,
                };
                match &element.kind {
                    StereoElementKind::Tetrahedral(value) => LabelStereo::Tetrahedral {
                        center: dense[value.center.index()],
                        carriers: value.carriers.iter().copied().map(carrier).collect(),
                        clockwise: value.orientation == Some(TetrahedralOrientation::Clockwise),
                    },
                    StereoElementKind::DoubleBond(value) => LabelStereo::DoubleBond {
                        left: dense[value.left.index()],
                        right: dense[value.right.index()],
                        left_carrier: carrier(value.left_carrier),
                        right_carrier: carrier(value.right_carrier),
                        together: value.orientation == Some(DoubleBondOrientation::Together),
                    },
                    StereoElementKind::Axis(_) => {
                        unreachable!("unsupported SMILES geometry is rejected before labeling")
                    }
                }
            })
            .collect();
        let stereo_atoms = atoms
            .iter()
            .map(|atom| {
                molecule
                    .stereo_elements()
                    .any(|(_, element)| element.references_atom(*atom))
            })
            .collect();
        let mut search = Search {
            adjacency,
            stereo,
            groups,
            stereo_atoms,
            initial,
            states: 0,
            work: 0,
            max_states,
            max_work,
            best: None,
            automorphisms: Vec::new(),
        };
        let order = search.run()?;
        let mut labels = vec![usize::MAX; slots];
        for (label, index) in order.into_iter().enumerate() {
            labels[atoms[index].index()] = label;
        }
        Ok(Self {
            initial_classes,
            labels,
        })
    }

    pub(super) fn rank(&self, atom: AtomId) -> (u32, usize) {
        (
            self.initial_classes[atom.index()],
            self.labels[atom.index()],
        )
    }
}

/// Writer-local priorities follow RDKit's map/degree/element/isotope/H/charge
/// ordering. Keep this separate from the public symmetry-class algorithm: the
/// writer orders classes, while canonical labeling still proves all remaining
/// ties using the complete adjacency/stereo certificate. This is not a claim of
/// exact RDKit stereo ranking or tie-breaking compatibility.
fn writer_atom_priority(molecule: &Molecule, id: AtomId) -> (u32, usize, u8, u16, usize, u32) {
    let atom = molecule.atom(id).expect("canonical atom is live");
    (
        atom.atom_map.unwrap_or(0),
        molecule
            .incident_bonds(id)
            .expect("canonical atom is live")
            .count(),
        atom.element.atomic_number(),
        atom.isotope.unwrap_or(0),
        molecule
            .implicit_hydrogens(id)
            .expect("canonical atom is live")
            .unwrap_or(0),
        // RDKit compares formal charge through an unsigned 32-bit value.
        // Sign-extend before conversion so negative charges follow neutral
        // and positive atoms while retaining their order among themselves.
        i32::from(atom.formal_charge) as u32,
    )
}

type Certificate = (
    Vec<(usize, Vec<(usize, u8)>)>,
    Vec<Vec<usize>>,
    Vec<(u8, Vec<Vec<usize>>)>,
);

enum LabelStereo {
    Tetrahedral {
        center: usize,
        carriers: Vec<usize>,
        clockwise: bool,
    },
    DoubleBond {
        left: usize,
        right: usize,
        left_carrier: usize,
        right_carrier: usize,
        together: bool,
    },
}

/// The order here is a serialization preference, not a chemical ordering of
/// stereo relationships. Correlated groups never refine by absolute handedness.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum StereoRefinementKind {
    None,
    Ungrouped,
    Absolute,
    Or,
    And,
    Relative,
    Racemic,
}

struct Search {
    stereo: Vec<LabelStereo>,
    groups: Vec<(StereoGroupKind, Vec<usize>)>,
    stereo_atoms: Vec<bool>,
    adjacency: Vec<Vec<(usize, u8)>>,
    initial: Vec<usize>,
    states: usize,
    work: usize,
    max_states: usize,
    max_work: usize,
    best: Option<(Certificate, Vec<usize>)>,
    automorphisms: Vec<Vec<usize>>,
}

struct Frame {
    colors: Vec<usize>,
    choices: Vec<usize>,
    next: usize,
}

impl Search {
    fn run(&mut self) -> Result<Vec<usize>, MolWriteError> {
        let mut pending = Vec::<Frame>::new();
        let mut colors = self.initial.clone();
        loop {
            self.states += 1;
            if self.states > self.max_states {
                return Err(limit("search states", self.max_states));
            }
            self.refine(&mut colors)?;
            let mut cells = BTreeMap::<usize, Vec<usize>>::new();
            for (atom, color) in colors.iter().copied().enumerate() {
                cells.entry(color).or_default().push(atom);
            }
            let cell = cells
                .values()
                .filter(|cell| cell.len() > 1)
                .min_by_key(|cell| (cell.len(), colors[cell[0]]));
            if let Some(cell) = cell {
                let mut choices = Vec::<usize>::new();
                for &atom in cell {
                    let mut redundant = false;
                    for &prior in &choices {
                        self.charge_work(
                            self.adjacency[prior]
                                .len()
                                .saturating_add(self.adjacency[atom].len())
                                .saturating_add(1),
                        )?;
                        if self.twins(prior, atom) {
                            redundant = true;
                            break;
                        }
                    }
                    if !redundant {
                        choices.push(atom);
                    }
                }
                if pending.len().saturating_add(1).saturating_mul(colors.len()) > MAX_PENDING_ATOMS
                {
                    return Err(limit("pending atom labels", MAX_PENDING_ATOMS));
                }
                pending.push(Frame {
                    colors,
                    choices,
                    next: 0,
                });
            } else {
                self.consider(&colors);
            }
            // An explicit stack keeps the search and its failure cleanup safe
            // even when a caller constructs a graph with a very deep partition.
            loop {
                let Some(frame) = pending.last_mut() else {
                    return Ok(self
                        .best
                        .take()
                        .expect("a completed search visits a leaf")
                        .1);
                };
                if let Some(&atom) = frame.choices.get(frame.next) {
                    let redundant = self.equivalent_to_prior_choice(
                        &frame.colors,
                        &frame.choices[..frame.next],
                        atom,
                    )?;
                    frame.next += 1;
                    if redundant {
                        continue;
                    }
                    colors = frame.colors.clone();
                    let chosen = colors[atom];
                    for (index, color) in colors.iter_mut().enumerate() {
                        *color = *color * 2 + usize::from(*color == chosen && index != atom);
                    }
                    break;
                }
                pending.pop();
            }
        }
    }

    fn equivalent_to_prior_choice(
        &mut self,
        colors: &[usize],
        prior: &[usize],
        atom: usize,
    ) -> Result<bool, MolWriteError> {
        if prior.is_empty() || self.automorphisms.is_empty() {
            return Ok(false);
        }
        self.charge_work(colors.len().saturating_mul(self.automorphisms.len()))?;
        // Only automorphisms preserving every current partition cell can map
        // two branches of this frame. Their compositions preserve it as well.
        let generators = self
            .automorphisms
            .iter()
            .filter(|permutation| {
                permutation
                    .iter()
                    .enumerate()
                    .all(|(from, &to)| colors[from] == colors[to])
            })
            .collect::<Vec<_>>();
        let mut seen = vec![false; colors.len()];
        let mut pending = prior.to_vec();
        for &value in prior {
            seen[value] = true;
        }
        let mut work = 0usize;
        while let Some(value) = pending.pop() {
            for permutation in &generators {
                work = work.saturating_add(1);
                let other = permutation[value];
                if other == atom {
                    self.charge_work(work)?;
                    return Ok(true);
                }
                if !seen[other] {
                    seen[other] = true;
                    pending.push(other);
                }
            }
        }
        self.charge_work(work)?;
        Ok(false)
    }

    fn refine(&mut self, colors: &mut Vec<usize>) -> Result<(), MolWriteError> {
        loop {
            let cost = colors
                .len()
                .saturating_add(self.adjacency.iter().map(Vec::len).sum::<usize>())
                .saturating_add(self.stereo.len())
                .saturating_add(
                    self.groups
                        .iter()
                        .map(|(_, members)| members.len())
                        .sum::<usize>(),
                );
            self.charge_work(cost)?;
            let stereo_keys = self.stereo_refinement_keys(colors);
            let signatures = self
                .adjacency
                .iter()
                .enumerate()
                .map(|(atom, neighbors)| {
                    let mut neighborhood = neighbors
                        .iter()
                        .map(|&(other, bond)| (bond, colors[other]))
                        .collect::<Vec<_>>();
                    // RDKit compares neighbor lists in descending bond/class
                    // order. The full certificate still resolves residual ties.
                    neighborhood.sort_unstable_by(|a, b| b.cmp(a));
                    (colors[atom], &stereo_keys[atom], neighborhood)
                })
                .collect::<Vec<_>>();
            let mut ordered = signatures.clone();
            ordered.sort_unstable();
            ordered.dedup();
            let refined = signatures
                .iter()
                .map(|signature| {
                    ordered
                        .binary_search(signature)
                        .expect("every refinement signature is indexed")
                })
                .collect::<Vec<_>>();
            if *colors == refined {
                return Ok(());
            }
            *colors = refined;
        }
    }

    fn stereo_refinement_keys(&self, colors: &[usize]) -> Vec<(StereoRefinementKind, usize, u8)> {
        // Build each group signature once. Copying its member list into every
        // atom key would make a large correlated group quadratic in memory.
        let group_keys = self
            .groups
            .iter()
            .map(|(kind, members)| {
                let kind = match kind {
                    StereoGroupKind::Absolute => StereoRefinementKind::Absolute,
                    StereoGroupKind::Or => StereoRefinementKind::Or,
                    StereoGroupKind::And => StereoRefinementKind::And,
                    StereoGroupKind::Relative => StereoRefinementKind::Relative,
                    StereoGroupKind::Racemic => StereoRefinementKind::Racemic,
                };
                let mut ranks = members
                    .iter()
                    .filter_map(|i| match &self.stereo[*i] {
                        LabelStereo::Tetrahedral { center, .. } => Some(colors[*center]),
                        _ => None,
                    })
                    .collect::<Vec<_>>();
                ranks.sort_unstable();
                ranks.dedup();
                (kind, ranks)
            })
            .collect::<Vec<_>>();
        let mut sorted = group_keys.iter().collect::<Vec<_>>();
        sorted.sort_unstable();
        sorted.dedup();
        let mut memberships = vec![(StereoRefinementKind::Ungrouped, 0); self.stereo.len()];
        for ((_, members), key) in self.groups.iter().zip(&group_keys) {
            let class = sorted.binary_search(&key).expect("group key is indexed");
            for &member in members {
                memberships[member] = (key.0, class);
            }
        }
        let mut keys = vec![(StereoRefinementKind::None, 0, 0); colors.len()];
        for (index, element) in self.stereo.iter().enumerate() {
            let LabelStereo::Tetrahedral {
                center,
                carriers,
                clockwise,
            } = element
            else {
                continue;
            };
            let (kind, members) = memberships[index];
            let mut chirality = 0;
            // A correlated group and its joint inverse are the same assertion.
            // Its absolute orientation must not split the refinement classes.
            if matches!(
                kind,
                StereoRefinementKind::Ungrouped | StereoRefinementKind::Absolute
            ) {
                let ranks = carriers
                    .iter()
                    .map(|atom| colors.get(*atom).copied().unwrap_or(*atom))
                    .collect::<Vec<_>>();
                let mut sorted = ranks.clone();
                sorted.sort_unstable();
                sorted.dedup();
                if sorted.len() == ranks.len() {
                    let odd = (0..ranks.len())
                        .flat_map(|i| (i + 1..ranks.len()).map(move |j| (i, j)))
                        .filter(|&(i, j)| ranks[i] > ranks[j])
                        .count()
                        % 2
                        != 0;
                    // Core Clockwise spells SMILES @, which RDKit names CCW.
                    chirality = 1 + u8::from(*clockwise == odd);
                }
            }
            keys[*center] = (kind, members, chirality);
        }
        keys
    }

    fn charge_work(&mut self, cost: usize) -> Result<(), MolWriteError> {
        self.work = self.work.saturating_add(cost);
        if self.work > self.max_work {
            Err(limit("refinement work", self.max_work))
        } else {
            Ok(())
        }
    }

    fn twins(&self, left: usize, right: usize) -> bool {
        !self.stereo_atoms[left]
            && !self.stereo_atoms[right]
            && self.adjacency[left]
                .iter()
                .filter(|(atom, _)| *atom != right)
                .eq(self.adjacency[right]
                    .iter()
                    .filter(|(atom, _)| *atom != left))
    }

    fn consider(&mut self, colors: &[usize]) {
        let mut order = (0..colors.len()).collect::<Vec<_>>();
        order.sort_unstable_by_key(|&atom| colors[atom]);
        let certificate = order
            .iter()
            .map(|&atom| {
                let mut neighbors = self.adjacency[atom]
                    .iter()
                    .map(|&(other, bond)| (colors[other], bond))
                    .collect::<Vec<_>>();
                neighbors.sort_unstable();
                (self.initial[atom], neighbors)
            })
            .collect::<Vec<_>>();
        let mut stereo = self
            .stereo
            .iter()
            .map(|element| {
                let label = |atom: usize| colors.get(atom).copied().unwrap_or(atom);
                match element {
                    LabelStereo::Tetrahedral {
                        center,
                        carriers,
                        clockwise,
                    } => {
                        let mut carriers = carriers.iter().copied().map(label).collect::<Vec<_>>();
                        let odd = (0..carriers.len())
                            .flat_map(|i| (i + 1..carriers.len()).map(move |j| (i, j)))
                            .filter(|&(i, j)| carriers[i] > carriers[j])
                            .count()
                            % 2
                            != 0;
                        carriers.sort_unstable();
                        let mut value = vec![0, label(*center)];
                        value.extend(carriers);
                        value.push(usize::from(*clockwise == odd));
                        value
                    }
                    LabelStereo::DoubleBond {
                        left,
                        right,
                        left_carrier,
                        right_carrier,
                        together,
                    } => {
                        let first_carrier =
                            |atom: usize, other: usize| {
                                self.adjacency[atom].iter()
                        .filter(|&&(neighbor, bond)| neighbor != other && bond == 1)
                        .map(|&(neighbor, _)| label(neighbor)).min()
                        .expect("SMILES stereo endpoints have an explicit single-bond carrier")
                            };
                        let mut endpoints = [label(*left), label(*right)];
                        endpoints.sort_unstable();
                        let inverted = (label(*left_carrier) != first_carrier(*left, *right))
                            != (label(*right_carrier) != first_carrier(*right, *left));
                        vec![
                            1,
                            endpoints[0],
                            endpoints[1],
                            usize::from(*together != inverted),
                        ]
                    }
                }
            })
            .collect::<Vec<_>>();
        let mut groups = Vec::new();
        for (kind, members) in &self.groups {
            // A non-absolute group represents a correlated choice of a whole
            // configuration and its inverse, not independent choices per site.
            if *kind != StereoGroupKind::Absolute {
                let first = *members
                    .iter()
                    .min_by_key(|&&i| &stereo[i][..stereo[i].len() - 1])
                    .expect("validated nonempty group");
                if stereo[first].last() == Some(&1) {
                    for &i in members {
                        *stereo[i].last_mut().unwrap() ^= 1;
                    }
                }
            }
            let mut values = members
                .iter()
                .map(|&i| stereo[i][..stereo[i].len() - 1].to_vec())
                .collect::<Vec<_>>();
            values.sort();
            groups.push((*kind as u8, values));
        }
        groups.sort();
        stereo.sort_unstable();
        let certificate = (certificate, stereo, groups);
        if let Some((best, representative)) = &self.best {
            if *best == certificate {
                // Equal complete certificates prove a graph-and-stereo
                // automorphism. Reuse it to avoid equivalent subtrees.
                let mut permutation = vec![0; order.len()];
                for (&from, &to) in representative.iter().zip(&order) {
                    permutation[from] = to;
                }
                if permutation.iter().enumerate().any(|(from, &to)| from != to)
                    && !self.automorphisms.contains(&permutation)
                    && self
                        .automorphisms
                        .len()
                        .saturating_add(1)
                        .saturating_mul(order.len())
                        <= MAX_PENDING_ATOMS
                {
                    self.automorphisms.push(permutation);
                }
            }
        }
        if self
            .best
            .as_ref()
            .is_none_or(|(best, _)| certificate < *best)
        {
            self.best = Some((certificate, order));
        }
    }
}

fn limit(resource: &str, bound: usize) -> MolWriteError {
    MolWriteError::resource_limit(format!(
        "canonical SMILES {resource} limit {bound} exceeded before canonical labeling completed"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn incomplete_search_never_publishes_a_candidate() {
        let molecule = crate::tests::read_smiles("C1CCCCC1").expect("cyclohexane");
        let error = CanonicalOrder::with_limits(
            &molecule,
            CanonicalAtomStyle::Aromatic,
            1,
            MAX_REFINEMENT_WORK,
        )
        .err()
        .expect("search is incomplete");
        assert_eq!(
            error.kind(),
            crate::smiles::MolWriteErrorKind::ResourceLimit
        );
        assert!(error.message().contains("search states"));
        let error = CanonicalOrder::with_limits(
            &molecule,
            CanonicalAtomStyle::Aromatic,
            MAX_SEARCH_STATES,
            1,
        )
        .err()
        .expect("refinement is incomplete");
        assert_eq!(
            error.kind(),
            crate::smiles::MolWriteErrorKind::ResourceLimit
        );
        assert!(error.message().contains("refinement work"));
    }
}
