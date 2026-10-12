//! Symmetry-equivalent atom permutations for symmetry-corrected RMSD.

use std::collections::VecDeque;
use std::fmt;

use crate::core::Element;
use crate::topology::{AtomSelection, Topology, TopologyAtomIndex, TopologyPerceptionError};

/// An [`AtomSelection`] with every permutation of its atoms that a symmetry
/// of the molecular graph allows, for symmetry-corrected fitting and RMSD.
///
/// Pass it wherever [`super::FitAtoms`] is accepted. Fits and RMSDs then
/// report the minimum over all permutations, so equivalent atoms (the
/// hydrogens of a methyl group, the ortho carbons of a flipped phenyl ring,
/// the oxygens of a carboxylate) are paired with whichever counterpart is
/// closest. Ties keep the earliest permutation; the identity is first.
///
/// # Equivalence
///
/// A permutation is admitted when it extends to an automorphism of the
/// covalent graph of the selected atoms' molecule instances that maps the
/// selection onto itself. The automorphism preserves element, isotope, total
/// hydrogen count, and connectivity, and never exchanges atoms between
/// molecule instances, so identical molecules (two ligand copies, waters,
/// homo-oligomer chains) are not swapped. Bond orders, formal charges, and
/// aromaticity are ignored, so localized resonance forms of carboxylates,
/// nitro groups, guanidinium, and Kekulé rings stay symmetric. Inferred
/// hydrogen counts of an unperceived topology are read from a perceived
/// snapshot; the topology itself is unchanged. Stereo is not considered: a
/// permutation that mirrors a stereocenter is admitted, and simply gives a
/// larger RMSD for a rigid structure.
///
/// # Cost
///
/// The search is exact: graph refinement only rules out mappings. Each
/// permutation is fitted or measured separately, so fused costs grow with
/// [`Self::permutation_count`]. The count multiplies across independent
/// symmetric groups, and selecting explicit hydrogens of many methyl groups
/// quickly reaches the limits; select heavy atoms, or omit symmetric side
/// chains of large biomolecules, to keep the count small. Exceeding a limit
/// is an error, never a truncated result.
///
/// The value is bound to its selection's topology layout and is reusable
/// across every model, member, and frame of that layout.
#[derive(Debug, Clone)]
pub struct SymmetricSelection {
    selection: AtomSelection,
    /// For each permutation and each selected atom in selection order, the
    /// selection position of the moving atom paired with it.
    permutations: Vec<Box<[usize]>>,
}

impl SymmetricSelection {
    /// Finds the permutations of `selection` within [`SymmetryLimits::default`].
    pub fn new(selection: &AtomSelection) -> Result<Self, SymmetryError> {
        Self::with_limits(selection, SymmetryLimits::default())
    }

    /// Finds the permutations of `selection` within explicit limits.
    pub fn with_limits(
        selection: &AtomSelection,
        limits: SymmetryLimits,
    ) -> Result<Self, SymmetryError> {
        let topology = selection.topology();
        let search = match Search::new(topology, selection.indices()) {
            Ok(search) => search,
            // Inferred hydrogen counts are known only after perception.
            Err(_) => {
                let perceived = topology.perceived().map_err(SymmetryError::Perception)?;
                Search::new(&perceived, selection.indices())
                    .map_err(|atom| SymmetryError::UnknownHydrogenCount { atom })?
            }
        };
        let permutations = search.run(limits)?;
        Ok(Self {
            selection: selection.clone(),
            permutations,
        })
    }

    /// The selection whose atoms are permuted.
    pub fn selection(&self) -> &AtomSelection {
        &self.selection
    }

    /// The number of selected atoms.
    pub fn len(&self) -> usize {
        self.selection.len()
    }

    pub fn is_empty(&self) -> bool {
        self.selection.is_empty()
    }

    /// The number of admitted permutations, including the identity.
    pub fn permutation_count(&self) -> usize {
        self.permutations.len()
    }

    /// The admitted permutations, identity first. Entry `p` of a permutation
    /// is the selection position of the moving atom paired with the
    /// reference atom at selection position `p`.
    pub fn permutations(&self) -> impl ExactSizeIterator<Item = &[usize]> {
        self.permutations.iter().map(AsRef::as_ref)
    }

    /// Moving/reference index pairs for every permutation, in permutation
    /// order and reference selection order.
    pub(crate) fn pair_sets(&self) -> Vec<Vec<(TopologyAtomIndex, TopologyAtomIndex)>> {
        let indices = self.selection.indices();
        self.permutations
            .iter()
            .map(|permutation| {
                permutation
                    .iter()
                    .zip(indices)
                    .map(|(moving, reference)| (indices[*moving], *reference))
                    .collect()
            })
            .collect()
    }
}

/// Bounds on the symmetry search of a [`SymmetricSelection`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SymmetryLimits {
    /// Most admitted permutations, including the identity.
    pub max_permutations: usize,
    /// Most candidate atom mappings examined by the search.
    pub max_search_steps: usize,
}

impl Default for SymmetryLimits {
    fn default() -> Self {
        Self {
            max_permutations: 4096,
            max_search_steps: 1_000_000,
        }
    }
}

/// Failure to find the symmetry permutations of a selection.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SymmetryError {
    /// The selection admits more permutations than allowed.
    PermutationLimit { limit: usize },
    /// The search examined more candidate mappings than allowed.
    SearchLimit { limit: usize },
    /// Perceiving inferred hydrogen counts failed.
    Perception(TopologyPerceptionError),
    /// An atom's hydrogen count is unknown even after perception.
    UnknownHydrogenCount { atom: TopologyAtomIndex },
}

impl fmt::Display for SymmetryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Perception(error) => {
                write!(formatter, "cannot perceive hydrogen counts: {error}")
            }
            Self::UnknownHydrogenCount { atom } => {
                write!(formatter, "hydrogen count of atom {atom} is unknown")
            }
            Self::PermutationLimit { limit } => write!(
                formatter,
                "selection admits more than {limit} symmetry permutations"
            ),
            Self::SearchLimit { limit } => write!(
                formatter,
                "symmetry search exceeded {limit} candidate mappings"
            ),
        }
    }
}

impl std::error::Error for SymmetryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Perception(error) => Some(error),
            _ => None,
        }
    }
}

const UNMAPPED: usize = usize::MAX;

/// Initial atom invariant; refinement only splits these classes.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct AtomKey {
    instance: usize,
    selected: bool,
    element: Element,
    isotope: Option<u16>,
    hydrogens: usize,
}

/// The covalent graph of the molecule instances holding selected atoms, in
/// local numbering.
struct Search {
    /// Local atom to topology atom.
    atoms: Vec<usize>,
    adjacency: Vec<Vec<usize>>,
    colors: Vec<usize>,
    /// Local atoms in search order: the selected atoms first, each later atom
    /// adjacent to an earlier one where the graph allows.
    order: Vec<usize>,
    /// Selection position of each local atom, if selected.
    positions: Vec<Option<usize>>,
    selected: usize,
}

impl Search {
    /// Fails with the first relevant atom whose hydrogen count is unknown.
    fn new(
        topology: &Topology,
        selection: &[TopologyAtomIndex],
    ) -> Result<Self, TopologyAtomIndex> {
        let mut instances = selection
            .iter()
            .map(|index| {
                topology
                    .atom_at(*index)
                    .expect("selection is bound to this topology")
                    .id()
                    .molecule()
            })
            .collect::<Vec<_>>();
        instances.sort_unstable();
        instances.dedup();

        let mut local = vec![UNMAPPED; topology.atom_count()];
        let mut atoms = Vec::new();
        let mut keys = Vec::new();
        for instance in &instances {
            let molecule = topology
                .molecule(*instance)
                .expect("selected atom has an instance");
            for atom in molecule.atoms() {
                local[atom.index().index()] = atoms.len();
                atoms.push(atom.index().index());
                keys.push(AtomKey {
                    instance: instance.index(),
                    selected: false,
                    element: atom.atom().element,
                    isotope: atom.atom().isotope,
                    hydrogens: atom.total_hydrogens().ok_or(atom.index())?,
                });
            }
        }
        let mut positions = vec![None; atoms.len()];
        for (position, index) in selection.iter().enumerate() {
            let atom = local[index.index()];
            positions[atom] = Some(position);
            keys[atom].selected = true;
        }
        let mut adjacency = vec![Vec::new(); atoms.len()];
        for instance in &instances {
            let molecule = topology.molecule(*instance).expect("collected instance");
            for bond in molecule.bonds() {
                let [a, b] = bond.atoms().map(|atom| local[atom.index().index()]);
                adjacency[a].push(b);
                adjacency[b].push(a);
            }
        }
        let colors = refined(&compressed(&keys), &adjacency);
        let order = search_order(selection, &local, &adjacency);
        Ok(Self {
            atoms,
            adjacency,
            colors,
            order,
            positions,
            selected: selection.len(),
        })
    }

    /// Enumerates each distinct restriction to the selection of the graph
    /// automorphisms that fix the selection setwise, identity first.
    fn run(self, limits: SymmetryLimits) -> Result<Vec<Box<[usize]>>, SymmetryError> {
        let identity = (0..self.selected).collect::<Box<[usize]>>();
        if self.selected == 0 {
            return Ok(vec![identity]);
        }
        let mut classes = vec![Vec::new(); self.atoms.len()];
        for (atom, color) in self.colors.iter().enumerate() {
            classes[*color].push(atom);
        }
        let count = self.order.len();
        // Candidates for the atom at `depth`, popped from the back; the atom
        // itself is last so the identity is found first.
        let candidates_at = |depth: usize| {
            let atom = self.order[depth];
            let mut class = classes[self.colors[atom]].clone();
            class.retain(|candidate| *candidate != atom);
            class.reverse();
            class.push(atom);
            class
        };
        let mut image = vec![UNMAPPED; count];
        let mut used = vec![false; count];
        // One candidate list per depth; every depth below the top is mapped,
        // and the top is the depth being mapped.
        let mut candidates = vec![candidates_at(0)];
        let mut permutations = Vec::new();
        let mut steps = 0_usize;
        loop {
            if self.advance(&mut candidates, &mut image, &mut used, &mut steps, limits)? {
                if candidates.len() < count {
                    candidates.push(candidates_at(candidates.len()));
                    continue;
                }
                if permutations.len() == limits.max_permutations {
                    return Err(SymmetryError::PermutationLimit {
                        limit: limits.max_permutations,
                    });
                }
                permutations.push(self.restriction(&image));
                // Other extensions of this restriction add nothing; resume
                // at the last selected atom.
                while candidates.len() > self.selected {
                    self.unassign(candidates.len() - 1, &mut image, &mut used);
                    candidates.pop();
                }
                self.unassign(candidates.len() - 1, &mut image, &mut used);
                continue;
            }
            // No candidate fits at the top depth; retry the one below.
            candidates.pop();
            match candidates.len() {
                0 => break,
                depth => self.unassign(depth - 1, &mut image, &mut used),
            }
        }
        debug_assert_eq!(permutations.first(), Some(&identity));
        Ok(permutations)
    }

    /// Assigns the next consistent candidate at the deepest level.
    fn advance(
        &self,
        candidates: &mut [Vec<usize>],
        image: &mut [usize],
        used: &mut [bool],
        steps: &mut usize,
        limits: SymmetryLimits,
    ) -> Result<bool, SymmetryError> {
        let depth = candidates.len() - 1;
        let atom = self.order[depth];
        while let Some(target) = candidates[depth].pop() {
            if *steps == limits.max_search_steps {
                return Err(SymmetryError::SearchLimit {
                    limit: limits.max_search_steps,
                });
            }
            *steps += 1;
            if self.consistent(atom, target, image, used) {
                image[atom] = target;
                used[target] = true;
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn unassign(&self, depth: usize, image: &mut [usize], used: &mut [bool]) {
        let atom = self.order[depth];
        if image[atom] != UNMAPPED {
            used[image[atom]] = false;
            image[atom] = UNMAPPED;
        }
    }

    /// Whether mapping `atom` to `target` keeps every edge and non-edge to
    /// already mapped atoms.
    fn consistent(&self, atom: usize, target: usize, image: &[usize], used: &[bool]) -> bool {
        if used[target] {
            return false;
        }
        let mut mapped = 0;
        for neighbor in &self.adjacency[atom] {
            let neighbor_image = image[*neighbor];
            if neighbor_image != UNMAPPED {
                if !self.adjacency[target].contains(&neighbor_image) {
                    return false;
                }
                mapped += 1;
            }
        }
        let target_mapped = self.adjacency[target]
            .iter()
            .filter(|neighbor| used[**neighbor])
            .count();
        mapped == target_mapped
    }

    fn restriction(&self, image: &[usize]) -> Box<[usize]> {
        let mut permutation = vec![0; self.selected];
        for (atom, position) in self.positions.iter().enumerate() {
            if let Some(position) = position {
                permutation[*position] =
                    self.positions[image[atom]].expect("automorphisms fix the selection setwise");
            }
        }
        permutation.into_boxed_slice()
    }
}

/// Dense ranks of `values`, preserving their order.
fn compressed<T: Ord + Clone>(values: &[T]) -> Vec<usize> {
    let mut ordered = values.to_vec();
    ordered.sort();
    ordered.dedup();
    values
        .iter()
        .map(|value| ordered.binary_search(value).expect("collected value"))
        .collect()
}

/// Color refinement by sorted neighbor colors until the partition is stable.
fn refined(colors: &[usize], adjacency: &[Vec<usize>]) -> Vec<usize> {
    let mut colors = colors.to_vec();
    let mut classes = colors.iter().max().map_or(0, |max| max + 1);
    loop {
        let signatures = adjacency
            .iter()
            .enumerate()
            .map(|(atom, neighbors)| {
                let mut environment = neighbors
                    .iter()
                    .map(|neighbor| colors[*neighbor])
                    .collect::<Vec<_>>();
                environment.sort_unstable();
                (colors[atom], environment)
            })
            .collect::<Vec<_>>();
        let next = compressed(&signatures);
        let next_classes = next.iter().max().map_or(0, |max| max + 1);
        if next_classes == classes {
            return colors;
        }
        colors = next;
        classes = next_classes;
    }
}

/// Selected atoms first, by breadth-first walks through selected atoms
/// seeded in selection order, then the remaining atoms breadth-first from
/// them, so each atom is adjacent to an earlier one where possible.
fn search_order(
    selection: &[TopologyAtomIndex],
    local: &[usize],
    adjacency: &[Vec<usize>],
) -> Vec<usize> {
    let count = adjacency.len();
    let mut selected = vec![false; count];
    for index in selection {
        selected[local[index.index()]] = true;
    }
    let mut visited = vec![false; count];
    let mut order = Vec::with_capacity(count);
    let mut queue = VecDeque::new();
    for index in selection {
        let start = local[index.index()];
        if visited[start] {
            continue;
        }
        visited[start] = true;
        queue.push_back(start);
        while let Some(atom) = queue.pop_front() {
            order.push(atom);
            for neighbor in &adjacency[atom] {
                if selected[*neighbor] && !visited[*neighbor] {
                    visited[*neighbor] = true;
                    queue.push_back(*neighbor);
                }
            }
        }
    }
    let mut queue = order.iter().copied().collect::<VecDeque<_>>();
    let mut next_start = 0;
    loop {
        while let Some(atom) = queue.pop_front() {
            for neighbor in &adjacency[atom] {
                if !visited[*neighbor] {
                    visited[*neighbor] = true;
                    order.push(*neighbor);
                    queue.push_back(*neighbor);
                }
            }
        }
        while next_start < count && visited[next_start] {
            next_start += 1;
        }
        if next_start == count {
            return order;
        }
        visited[next_start] = true;
        order.push(next_start);
        queue.push_back(next_start);
    }
}
