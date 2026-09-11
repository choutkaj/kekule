use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;

use crate::core::*;

pub(crate) fn compute_ring_membership(mol: &Molecule) -> RingMembership {
    let mut graph = vec![Vec::<(AtomId, BondId)>::new(); mol.graph.atom_slot_count()];
    let mut live_bonds = Vec::new();
    for (bond_id, bond) in mol.bonds() {
        if matches!(bond.order, BondOrder::Zero | BondOrder::Dative) {
            continue;
        }
        graph[bond.a.index()].push((bond.b, bond_id));
        graph[bond.b.index()].push((bond.a, bond_id));
        live_bonds.push(bond_id);
    }

    let mut discovery = vec![None; mol.graph.atom_slot_count()];
    let mut low = vec![0usize; mol.graph.atom_slot_count()];
    let mut bridge = vec![false; mol.graph.bond_slot_count()];
    let mut time = 0usize;
    for atom_id in mol.atom_ids() {
        if discovery[atom_id.index()].is_none() {
            ring_dfs_iterative(
                atom_id,
                &graph,
                &mut discovery,
                &mut low,
                &mut bridge,
                &mut time,
            );
        }
    }

    let mut membership = RingMembership {
        atom_flags: vec![false; mol.graph.atom_slot_count()],
        bond_flags: vec![false; mol.graph.bond_slot_count()],
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

/// Installs graph-theoretic cycle membership, excluding zero and dative bonds.
///
/// This linear-time traversal is independent of the model's selected rings.
/// Recomputing membership preserves valence and clears the ring set,
/// aromaticity, and dependent stereo perception.
pub fn perceive_ring_membership(mol: &mut Molecule) -> RingMembership {
    let membership = compute_ring_membership(mol);
    mol.install_ring_membership(membership.clone());
    membership
}

pub(crate) fn bond_in_ring_smaller_than(mol: &Molecule, bond_id: BondId, ring_size: usize) -> bool {
    let Ok(bond) = mol.bond(bond_id) else {
        return false;
    };
    if ring_size <= 1 || matches!(bond.order, BondOrder::Zero | BondOrder::Dative) {
        return false;
    }
    let max_path_edges = ring_size - 2;
    let mut seen = vec![false; mol.graph.atom_slot_count()];
    let mut queue = VecDeque::from([(bond.a(), 0usize)]);
    seen[bond.a().index()] = true;
    while let Some((atom, depth)) = queue.pop_front() {
        if atom == bond.b() {
            return true;
        }
        if depth == max_path_edges {
            continue;
        }
        let Ok(incident) = mol.incident_bonds(atom) else {
            continue;
        };
        for (next_bond, next) in incident
            .filter(|(_, edge)| !matches!(edge.order, BondOrder::Zero | BondOrder::Dative))
            .map(|(next_bond, edge)| (next_bond, edge.other_atom(atom)))
        {
            if next_bond == bond_id || seen.get(next.index()).copied().unwrap_or(true) {
                continue;
            }
            seen[next.index()] = true;
            queue.push_back((next, depth + 1));
        }
    }
    false
}

fn ring_dfs_iterative(
    start: AtomId,
    graph: &[Vec<(AtomId, BondId)>],
    discovery: &mut [Option<usize>],
    low: &mut [usize],
    bridge: &mut [bool],
    time: &mut usize,
) {
    struct Frame {
        atom: AtomId,
        parent_bond: Option<BondId>,
        next_edge: usize,
    }

    discovery[start.index()] = Some(*time);
    low[start.index()] = *time;
    *time += 1;
    let mut stack = vec![Frame {
        atom: start,
        parent_bond: None,
        next_edge: 0,
    }];
    while let Some(frame) = stack.last_mut() {
        if frame.next_edge >= graph[frame.atom.index()].len() {
            let finished = stack.pop().expect("bridge DFS frame should exist");
            if let (Some(parent), Some(parent_bond)) = (stack.last(), finished.parent_bond) {
                low[parent.atom.index()] = low[parent.atom.index()].min(low[finished.atom.index()]);
                if low[finished.atom.index()]
                    > discovery[parent.atom.index()].expect("parent atom is discovered")
                {
                    bridge[parent_bond.index()] = true;
                }
            }
            continue;
        }

        let atom = frame.atom;
        let parent_bond = frame.parent_bond;
        let (neighbor, bond_id) = graph[atom.index()][frame.next_edge];
        frame.next_edge += 1;
        if Some(bond_id) == parent_bond {
            continue;
        }
        if discovery[neighbor.index()].is_none() {
            discovery[neighbor.index()] = Some(*time);
            low[neighbor.index()] = *time;
            *time += 1;
            stack.push(Frame {
                atom: neighbor,
                parent_bond: Some(bond_id),
                next_edge: 0,
            });
        } else {
            low[atom.index()] =
                low[atom.index()].min(discovery[neighbor.index()].expect("neighbor discovered"));
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Resource bounds for selected-ring perception.
///
/// Total work includes search workspace initialization, search visits,
/// candidate generation, recovery graph and path copies, repeated recovery
/// bond scans, and ring comparisons during pruning and symmetrization.
pub struct RingPerceptionOptions {
    pub max_atoms: usize,
    pub max_bonds: usize,
    pub max_candidates: usize,
    pub max_path_expansions: usize,
    pub max_equivalent_shortest_paths: usize,
    pub max_cycle_size: usize,
    pub max_total_work: usize,
}

impl Default for RingPerceptionOptions {
    fn default() -> Self {
        Self {
            max_atoms: 1_000_000,
            max_bonds: 2_000_000,
            max_candidates: 100_000,
            max_path_expansions: 2_000_000,
            max_equivalent_shortest_paths: 100_000,
            max_cycle_size: 4_096,
            max_total_work: 5_000_000,
        }
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RingPerceptionError {
    /// A configured resource bound was exceeded.
    ResourceLimit {
        /// The bounded resource.
        resource: &'static str,
        /// The first observed value beyond the bound.
        observed: usize,
        /// The configured bound.
        limit: usize,
    },
    /// The candidate basis did not cover every bond known to be cyclic.
    IncompleteRingCoverage {
        /// Cyclic bonds not covered by the candidate basis.
        uncovered_bonds: Vec<BondId>,
    },
}

impl fmt::Display for RingPerceptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ResourceLimit {
                resource,
                observed,
                limit,
                ..
            } => write!(
                f,
                "ring perception {resource} limit exceeded: observed {observed}, limit {limit}"
            ),
            Self::IncompleteRingCoverage {
                uncovered_bonds, ..
            } => write!(
                f,
                "ring perception did not cover {} cyclic bond(s)",
                uncovered_bonds.len()
            ),
        }
    }
}

impl std::error::Error for RingPerceptionError {}

/// Installs RDKit-like Figueras rings and symmetric alternatives.
///
/// Incomplete candidate searches use deterministic depth-first cycles. This
/// selected-ring model does not promise a mathematical minimum or full-rank
/// cycle basis. Zero and dative bonds are excluded throughout.
///
/// Successful perception preserves valence and clears aromaticity and dependent
/// stereo. Failure preserves the previously installed perception state.
pub fn perceive_ring_set(mol: &mut Molecule) -> std::result::Result<RingSet, RingPerceptionError> {
    perceive_ring_set_with_options(mol, RingPerceptionOptions::default())
}

/// Like [`perceive_ring_set`], with explicit resource bounds for every search
/// phase and ring comparison. Bound failures never install partial state.
pub fn perceive_ring_set_with_options(
    mol: &mut Molecule,
    options: RingPerceptionOptions,
) -> std::result::Result<RingSet, RingPerceptionError> {
    let mut tracker = RingWorkTracker::new(options, mol.atom_count(), mol.bond_count())?;
    let membership = compute_ring_membership(mol);
    let (mut rings, mut extras) = figueras_sssr_candidates(mol, &membership, &mut tracker)?;
    if !uncovered_ring_bonds(mol, &membership, &rings).is_empty() {
        // RDKit can exceptionally omit a cyclic bond even after reaching its
        // candidate count. Keep true cycle membership and the installed ring
        // coverage invariant rather than copying that target defect.
        rings = depth_first_ring_basis(mol, &mut tracker)?;
        extras.clear();
        let uncovered_bonds = uncovered_ring_bonds(mol, &membership, &rings);
        if !uncovered_bonds.is_empty() {
            return Err(RingPerceptionError::IncompleteRingCoverage { uncovered_bonds });
        }
    }
    let rings = symmetrize_ring_set(rings, extras, mol.graph.bond_slot_count(), &mut tracker)?;
    let ring_set = RingSet::from_rings(rings);
    mol.install_ring_basis(
        membership,
        RingBasisModel::FiguerasSssrLike,
        ring_set.clone(),
    );
    Ok(ring_set)
}

#[derive(Clone)]
struct ActiveRingGraph {
    adjacency: Vec<Vec<(AtomId, BondId)>>,
    active_bonds: Vec<bool>,
    atom_degrees: Vec<usize>,
    copy_work: usize,
}

impl ActiveRingGraph {
    fn new(mol: &Molecule) -> Self {
        let mut adjacency = vec![Vec::new(); mol.graph.atom_slot_count()];
        let mut active_bonds = vec![false; mol.graph.bond_slot_count()];
        let mut atom_degrees = vec![0usize; mol.graph.atom_slot_count()];
        let mut adjacency_entries = 0usize;
        for (bond_id, bond) in mol.bonds() {
            if matches!(bond.order, BondOrder::Zero | BondOrder::Dative) {
                continue;
            }
            adjacency[bond.a.index()].push((bond.b, bond_id));
            adjacency[bond.b.index()].push((bond.a, bond_id));
            active_bonds[bond_id.index()] = true;
            atom_degrees[bond.a.index()] += 1;
            atom_degrees[bond.b.index()] += 1;
            adjacency_entries = adjacency_entries.saturating_add(2);
        }
        // Trimming only changes flags/degrees, so this allocation-size estimate
        // remains valid and constant-time to read for every recovery clone.
        let copy_work = adjacency
            .len()
            .saturating_mul(2)
            .saturating_add(active_bonds.len())
            .saturating_add(adjacency_entries);
        Self {
            adjacency,
            active_bonds,
            atom_degrees,
            copy_work,
        }
    }

    fn active_neighbors(&self, atom: AtomId) -> impl Iterator<Item = (AtomId, BondId)> + '_ {
        self.adjacency[atom.index()]
            .iter()
            .copied()
            .filter(|(_, bond)| self.active_bonds[bond.index()])
    }

    fn trim_atom(&mut self, atom: AtomId, changed: &mut VecDeque<AtomId>) {
        let incident = self.adjacency[atom.index()].clone();
        for (other, bond) in incident {
            if !self.active_bonds[bond.index()] {
                continue;
            }
            if self.atom_degrees[other.index()] <= 2 {
                changed.push_back(other);
            }
            self.active_bonds[bond.index()] = false;
            self.atom_degrees[other.index()] = self.atom_degrees[other.index()].saturating_sub(1);
            self.atom_degrees[atom.index()] = self.atom_degrees[atom.index()].saturating_sub(1);
        }
    }
}

fn figueras_sssr_candidates(
    mol: &Molecule,
    membership: &RingMembership,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<(Vec<Ring>, Vec<Ring>), RingPerceptionError> {
    let mut graph = ActiveRingGraph::new(mol);
    let fragments = active_fragments(mol, &graph);
    let mut seen_invariants = BTreeSet::<Vec<AtomId>>::new();
    let mut all_sssr = Vec::new();
    let mut all_extras = Vec::new();
    let mut discovered = RingMembership {
        atom_flags: vec![false; mol.graph.atom_slot_count()],
        bond_flags: vec![false; mol.graph.bond_slot_count()],
    };

    for fragment in fragments {
        if fragment.len() < 3 {
            continue;
        }
        let active_degree_sum = fragment
            .iter()
            .map(|atom| graph.atom_degrees[atom.index()])
            .sum::<usize>();
        let expected = (active_degree_sum / 2 + 1).saturating_sub(fragment.len());
        if expected == 0 {
            continue;
        }

        let mut changed = fragment
            .iter()
            .copied()
            .filter(|atom| graph.atom_degrees[atom.index()] < 2)
            .collect::<VecDeque<_>>();
        let mut done = vec![false; mol.graph.atom_slot_count()];
        let mut atoms_done = 0usize;
        let mut fragment_candidates = Vec::<Ring>::new();

        while atoms_done <= fragment.len().saturating_sub(3) {
            while let Some(atom) = changed.pop_front() {
                if done[atom.index()] {
                    continue;
                }
                done[atom.index()] = true;
                atoms_done += 1;
                graph.trim_atom(atom, &mut changed);
            }

            let d2_nodes = pick_degree_two_nodes(&fragment, &graph);
            if !d2_nodes.is_empty() {
                find_rings_from_degree_two_nodes(
                    &d2_nodes,
                    &mut graph,
                    &mut fragment_candidates,
                    &mut seen_invariants,
                    &mut discovered,
                    tracker,
                )?;
                for atom in d2_nodes {
                    if !done[atom.index()] {
                        done[atom.index()] = true;
                        atoms_done += 1;
                    }
                    graph.trim_atom(atom, &mut changed);
                }
            } else if atoms_done <= fragment.len().saturating_sub(3) {
                let Some(root) = fragment
                    .iter()
                    .copied()
                    .find(|atom| graph.atom_degrees[atom.index()] == 3)
                else {
                    break;
                };
                find_rings_from_degree_three_node(
                    root,
                    &graph,
                    &mut fragment_candidates,
                    &mut seen_invariants,
                    tracker,
                )?;
                if !done[root.index()] {
                    done[root.index()] = true;
                    atoms_done += 1;
                }
                graph.trim_atom(root, &mut changed);
            }
        }

        if fragment_candidates.len() < expected {
            recover_connecting_cycles(
                mol,
                &fragment,
                &mut fragment_candidates,
                &mut seen_invariants,
                &mut discovered,
                tracker,
            )?;
            if fragment_candidates.len() < expected {
                // RDKit falls back for the whole graph when Figueras cannot
                // reach the fragment's cyclomatic count. Preserve that model
                // choice instead of substituting a different shortest basis.
                return Ok((depth_first_ring_basis(mol, tracker)?, Vec::new()));
            }
        }
        let (kept, extras) = if fragment_candidates.len() > expected {
            remove_extra_rings(fragment_candidates, mol.graph.bond_slot_count(), tracker)?
        } else {
            (fragment_candidates, Vec::new())
        };
        all_sssr.extend(kept);
        all_extras.extend(extras);
    }

    // Candidate search uses the full active graph like RDKit. True membership
    // remains authoritative for the caller's final coverage validation.
    debug_assert!(all_sssr
        .iter()
        .flat_map(|ring| &ring.bonds)
        .all(|bond| membership.bond_in_ring(*bond)));
    Ok((all_sssr, all_extras))
}

fn active_fragments(mol: &Molecule, graph: &ActiveRingGraph) -> Vec<Vec<AtomId>> {
    let mut seen = vec![false; mol.graph.atom_slot_count()];
    let mut fragments = Vec::new();
    for start in mol.atom_ids() {
        if seen[start.index()] {
            continue;
        }
        seen[start.index()] = true;
        let mut stack = vec![start];
        let mut fragment = Vec::new();
        while let Some(atom) = stack.pop() {
            fragment.push(atom);
            for (neighbor, _) in graph.active_neighbors(atom) {
                if !seen[neighbor.index()] {
                    seen[neighbor.index()] = true;
                    stack.push(neighbor);
                }
            }
        }
        fragment.sort();
        fragments.push(fragment);
    }
    fragments
}

fn pick_degree_two_nodes(fragment: &[AtomId], graph: &ActiveRingGraph) -> Vec<AtomId> {
    let mut forbidden = vec![false; graph.atom_degrees.len()];
    let mut roots = Vec::new();
    for root in fragment.iter().copied() {
        if graph.atom_degrees[root.index()] != 2 || forbidden[root.index()] {
            continue;
        }
        roots.push(root);
        forbidden[root.index()] = true;
        let mut stack = vec![root];
        while let Some(atom) = stack.pop() {
            for (neighbor, _) in graph.active_neighbors(atom) {
                if !forbidden[neighbor.index()] && graph.atom_degrees[neighbor.index()] == 2 {
                    forbidden[neighbor.index()] = true;
                    stack.push(neighbor);
                }
            }
        }
    }
    roots
}

fn find_rings_from_degree_two_nodes(
    roots: &[AtomId],
    graph: &mut ActiveRingGraph,
    candidates: &mut Vec<Ring>,
    seen_invariants: &mut BTreeSet<Vec<AtomId>>,
    discovered: &mut RingMembership,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<(), RingPerceptionError> {
    let mut duplicate_roots = BTreeMap::<Vec<AtomId>, Vec<AtomId>>::new();
    let mut duplicate_map = BTreeMap::<AtomId, Vec<AtomId>>::new();

    for root in roots {
        let atom_rings = smallest_rings_bfs(*root, graph, &BTreeSet::new(), tracker)?;
        for atoms in &atom_rings {
            let invariant = ring_invariant(atoms);
            let prior_roots = duplicate_roots.entry(invariant.clone()).or_default();
            if seen_invariants.insert(invariant) {
                let ring = atom_ring_to_ring(atoms, graph, tracker)?;
                mark_discovered_ring(&ring, discovered);
                candidates.push(ring);
            } else {
                for other in prior_roots.iter().copied() {
                    duplicate_map.entry(*root).or_default().push(other);
                    duplicate_map.entry(other).or_default().push(*root);
                }
            }
            prior_roots.push(*root);
        }
        if atom_rings.is_empty() {
            let mut changed = VecDeque::from([*root]);
            while let Some(atom) = changed.pop_front() {
                graph.trim_atom(atom, &mut changed);
            }
        }
    }

    recover_duplicate_degree_two_candidates(
        graph,
        &duplicate_roots,
        &duplicate_map,
        candidates,
        seen_invariants,
        tracker,
    )
}

fn recover_duplicate_degree_two_candidates(
    graph: &ActiveRingGraph,
    duplicate_roots: &BTreeMap<Vec<AtomId>, Vec<AtomId>>,
    duplicate_map: &BTreeMap<AtomId, Vec<AtomId>>,
    candidates: &mut Vec<Ring>,
    seen_invariants: &mut BTreeSet<Vec<AtomId>>,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<(), RingPerceptionError> {
    for roots in duplicate_roots.values() {
        if roots.len() <= 1 {
            continue;
        }
        let mut recovered = Vec::<Vec<AtomId>>::new();
        let mut minimum_size = usize::MAX;
        for root in roots {
            tracker.add_work(graph.copy_work)?;
            let mut reduced = graph.clone();
            let mut changed = VecDeque::new();
            for duplicate in duplicate_map.get(root).into_iter().flatten().copied() {
                reduced.trim_atom(duplicate, &mut changed);
            }
            let atom_rings = smallest_rings_bfs(*root, &reduced, &BTreeSet::new(), tracker)?;
            for atoms in atom_rings {
                minimum_size = minimum_size.min(atoms.len());
                recovered.push(atoms);
            }
        }
        for atoms in recovered
            .into_iter()
            .filter(|atoms| atoms.len() == minimum_size)
        {
            if seen_invariants.insert(ring_invariant(&atoms)) {
                candidates.push(atom_ring_to_ring(&atoms, graph, tracker)?);
            }
        }
    }
    Ok(())
}

fn find_rings_from_degree_three_node(
    root: AtomId,
    graph: &ActiveRingGraph,
    candidates: &mut Vec<Ring>,
    seen_invariants: &mut BTreeSet<Vec<AtomId>>,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<(), RingPerceptionError> {
    let smallest = smallest_rings_bfs(root, graph, &BTreeSet::new(), tracker)?;
    store_unique_atom_rings(&smallest, graph, candidates, seen_invariants, tracker)?;
    if smallest.len() >= 3 {
        return Ok(());
    }
    let neighbors = graph
        .active_neighbors(root)
        .map(|(atom, _)| atom)
        .take(3)
        .collect::<Vec<_>>();
    if neighbors.len() < 3 {
        return Ok(());
    }

    if smallest.len() == 2 {
        if let Some(common) = neighbors
            .iter()
            .copied()
            .find(|neighbor| smallest[0].contains(neighbor) && smallest[1].contains(neighbor))
        {
            let forbidden = BTreeSet::from([common]);
            let rings = smallest_rings_bfs(root, graph, &forbidden, tracker)?;
            store_unique_atom_rings(&rings, graph, candidates, seen_invariants, tracker)?;
        }
    } else if smallest.len() == 1 {
        let absent = neighbors
            .iter()
            .copied()
            .filter(|neighbor| !smallest[0].contains(neighbor))
            .collect::<Vec<_>>();
        if absent.len() == 1 {
            let included = neighbors
                .iter()
                .copied()
                .filter(|neighbor| *neighbor != absent[0])
                .collect::<Vec<_>>();
            for forbidden_neighbor in included.into_iter().rev() {
                let forbidden = BTreeSet::from([forbidden_neighbor]);
                let rings = smallest_rings_bfs(root, graph, &forbidden, tracker)?;
                store_unique_atom_rings(&rings, graph, candidates, seen_invariants, tracker)?;
            }
        }
    }
    Ok(())
}

fn store_unique_atom_rings(
    rings: &[Vec<AtomId>],
    graph: &ActiveRingGraph,
    candidates: &mut Vec<Ring>,
    seen_invariants: &mut BTreeSet<Vec<AtomId>>,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<(), RingPerceptionError> {
    for atoms in rings {
        if seen_invariants.insert(ring_invariant(atoms)) {
            candidates.push(atom_ring_to_ring(atoms, graph, tracker)?);
        }
    }
    Ok(())
}

fn smallest_rings_bfs(
    root: AtomId,
    graph: &ActiveRingGraph,
    forbidden: &BTreeSet<AtomId>,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<Vec<Vec<AtomId>>, RingPerceptionError> {
    const WHITE: u8 = 0;
    const GRAY: u8 = 1;
    const BLACK: u8 = 2;
    tracker.add_work(graph.atom_degrees.len().saturating_mul(3))?;
    let mut colors = vec![WHITE; graph.atom_degrees.len()];
    for atom in forbidden {
        colors[atom.index()] = BLACK;
    }
    let mut parents = vec![None; graph.atom_degrees.len()];
    let mut depths = vec![0usize; graph.atom_degrees.len()];
    let mut queue = VecDeque::from([root]);
    let mut rings = Vec::new();
    let mut current_size = usize::MAX;
    while let Some(current) = queue.pop_front() {
        colors[current.index()] = BLACK;
        let depth = depths[current.index()].saturating_add(1);
        if depth > current_size {
            break;
        }
        for (neighbor, _) in graph.active_neighbors(current) {
            tracker.record_path_expansion()?;
            if colors[neighbor.index()] == BLACK || parents[current.index()] == Some(neighbor) {
                continue;
            }
            if colors[neighbor.index()] == WHITE {
                parents[neighbor.index()] = Some(current);
                colors[neighbor.index()] = GRAY;
                depths[neighbor.index()] = depth;
                queue.push_back(neighbor);
                continue;
            }

            let mut ring = vec![neighbor];
            let mut parent = parents[neighbor.index()];
            while let Some(atom) = parent {
                if atom == root {
                    break;
                }
                ring.push(atom);
                parent = parents[atom.index()];
            }
            ring.insert(0, current);
            parent = parents[current.index()];
            while let Some(atom) = parent {
                if ring.contains(&atom) {
                    ring.clear();
                    break;
                }
                ring.insert(0, atom);
                parent = parents[atom.index()];
            }
            if ring.len() > 1 {
                if ring.len() <= current_size {
                    tracker.check("cycle size", ring.len(), tracker.options.max_cycle_size)?;
                    tracker.record_shortest_path()?;
                    current_size = ring.len();
                    rings.push(ring);
                } else {
                    return Ok(rings);
                }
            }
        }
    }
    Ok(rings)
}

fn ring_invariant(atoms: &[AtomId]) -> Vec<AtomId> {
    let mut invariant = atoms.to_vec();
    invariant.sort();
    invariant.dedup();
    invariant
}

fn atom_ring_to_ring(
    atoms: &[AtomId],
    graph: &ActiveRingGraph,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<Ring, RingPerceptionError> {
    let mut bonds = Vec::with_capacity(atoms.len());
    for index in 0..atoms.len() {
        let left = atoms[index];
        let right = atoms[(index + 1) % atoms.len()];
        let bond = graph.adjacency[left.index()]
            .iter()
            .find_map(|(neighbor, bond)| (*neighbor == right).then_some(*bond))
            .expect("BFS ring edges must exist in the molecular graph");
        bonds.push(bond);
    }
    bonds.sort();
    tracker.record_candidate()?;
    Ok(Ring {
        atoms: atoms.to_vec(),
        bonds,
    })
}

fn remove_extra_rings(
    mut rings: Vec<Ring>,
    bond_slots: usize,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<(Vec<Ring>, Vec<Ring>), RingPerceptionError> {
    rings.sort_by_key(|ring| ring.bonds.len());
    let mut available = vec![true; rings.len()];
    let mut keep = vec![false; rings.len()];
    let mut union = vec![false; bond_slots];

    for index in 0..rings.len() {
        tracker.add_work(rings[index].bonds.len())?;
        if ring_is_subset_of(&rings[index], &union) {
            available[index] = false;
        }
        if !available[index] {
            continue;
        }
        add_ring_to_union(&rings[index], &mut union);
        keep[index] = true;
        let mut consider = BTreeSet::new();
        for other in index + 1..rings.len() {
            tracker.add_work(1)?;
            if available[other] && rings[other].bonds.len() == rings[index].bonds.len() {
                consider.insert(other);
            }
        }
        while !consider.is_empty() {
            let mut best = None;
            let mut best_overlap = None;
            for other in consider.iter().copied() {
                tracker.add_work(rings[other].bonds.len())?;
                let overlap = rings[other]
                    .bonds
                    .iter()
                    .filter(|bond| union[bond.index()])
                    .count();
                if best_overlap.is_none_or(|current| overlap > current) {
                    best = Some(other);
                    best_overlap = Some(overlap);
                }
            }
            let best = best.expect("nonempty candidate set has a best overlap");
            consider.remove(&best);
            if ring_is_subset_of(&rings[best], &union) {
                available[best] = false;
            } else {
                keep[best] = true;
                available[best] = false;
                add_ring_to_union(&rings[best], &mut union);
            }
        }
    }

    let mut kept = Vec::new();
    let mut extras = Vec::new();
    for (index, ring) in rings.into_iter().enumerate() {
        if keep[index] {
            kept.push(ring);
        } else {
            extras.push(ring);
        }
    }
    Ok((kept, extras))
}

fn ring_is_subset_of(ring: &Ring, union: &[bool]) -> bool {
    ring.bonds.iter().all(|bond| union[bond.index()])
}

fn add_ring_to_union(ring: &Ring, union: &mut [bool]) {
    for bond in &ring.bonds {
        union[bond.index()] = true;
    }
}

fn mark_discovered_ring(ring: &Ring, membership: &mut RingMembership) {
    for atom in &ring.atoms {
        membership.atom_flags[atom.index()] = true;
    }
    for bond in &ring.bonds {
        membership.bond_flags[bond.index()] = true;
    }
}

fn recover_connecting_cycles(
    mol: &Molecule,
    fragment: &[AtomId],
    rings: &mut Vec<Ring>,
    invariants: &mut BTreeSet<Vec<AtomId>>,
    discovered: &mut RingMembership,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<(), RingPerceptionError> {
    let graph = ActiveRingGraph::new(mol);
    let fragment_atoms = fragment.iter().copied().collect::<BTreeSet<_>>();
    let mut dead = vec![false; mol.graph.bond_slot_count()];
    loop {
        tracker.add_work(mol.graph.bond_slot_count())?;
        let candidate = mol.bonds().find(|(id, bond)| {
            graph.active_bonds[id.index()]
                && fragment_atoms.contains(&bond.a())
                && !dead[id.index()]
                && !discovered.bond_in_ring(*id)
                && discovered.atom_in_ring(bond.a())
                && discovered.atom_in_ring(bond.b())
        });
        let Some((bond_id, bond)) = candidate else {
            return Ok(());
        };
        let mut queue = VecDeque::from([vec![bond.a()]]);
        let mut recovered = None;
        'search: while let Some(path) = queue.pop_front() {
            let current = *path.last().expect("search path is nonempty");
            for (neighbor, _) in graph.active_neighbors(current) {
                tracker.record_path_expansion()?;
                if neighbor == bond.b() {
                    if current == bond.a() {
                        continue;
                    }
                    let mut atoms = path.clone();
                    atoms.push(neighbor);
                    if invariants.contains(&ring_invariant(&atoms)) {
                        continue;
                    }
                    tracker.check("cycle size", atoms.len(), tracker.options.max_cycle_size)?;
                    tracker.record_shortest_path()?;
                    recovered = Some(atom_ring_to_ring(&atoms, &graph, tracker)?);
                    break 'search;
                }
                if discovered.atom_in_ring(neighbor) && !path.contains(&neighbor) {
                    tracker.check("cycle size", path.len() + 1, tracker.options.max_cycle_size)?;
                    // Bound copied path storage as well as adjacency visits.
                    tracker.add_work(path.len() + 1)?;
                    let mut next = path.clone();
                    next.push(neighbor);
                    queue.push_back(next);
                }
            }
        }
        if let Some(ring) = recovered {
            invariants.insert(ring_invariant(&ring.atoms));
            mark_discovered_ring(&ring, discovered);
            rings.push(ring);
        } else {
            dead[bond_id.index()] = true;
        }
    }
}

fn depth_first_ring_basis(
    mol: &Molecule,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<Vec<Ring>, RingPerceptionError> {
    struct Frame {
        atom: AtomId,
        parent: Option<AtomId>,
        next_edge: usize,
    }
    let graph = ActiveRingGraph::new(mol);
    let mut colors = vec![0; mol.graph.atom_slot_count()];
    let mut positions = vec![0; mol.graph.atom_slot_count()];
    let mut rings = Vec::new();
    let mut stack = Vec::<Frame>::new();
    for root in mol.atom_ids() {
        if colors[root.index()] != 0 {
            continue;
        }
        if graph.atom_degrees[root.index()] < 2 {
            colors[root.index()] = 2;
            continue;
        }
        colors[root.index()] = 1;
        positions[root.index()] = 0;
        stack.push(Frame {
            atom: root,
            parent: None,
            next_edge: 0,
        });
        while let Some(frame) = stack.last_mut() {
            if frame.next_edge == graph.adjacency[frame.atom.index()].len() {
                colors[frame.atom.index()] = 2;
                stack.pop();
                continue;
            }
            let atom = frame.atom;
            let parent = frame.parent;
            let (neighbor, _) = graph.adjacency[atom.index()][frame.next_edge];
            frame.next_edge += 1;
            tracker.record_path_expansion()?;
            if colors[neighbor.index()] == 0 {
                if graph.atom_degrees[neighbor.index()] < 2 {
                    colors[neighbor.index()] = 2;
                    continue;
                }
                colors[neighbor.index()] = 1;
                positions[neighbor.index()] = stack.len();
                stack.push(Frame {
                    atom: neighbor,
                    parent: Some(atom),
                    next_edge: 0,
                });
            } else if colors[neighbor.index()] == 1 && parent != Some(neighbor) {
                let atoms = stack[positions[neighbor.index()]..]
                    .iter()
                    .rev()
                    .map(|frame| frame.atom)
                    .collect::<Vec<_>>();
                tracker.check("cycle size", atoms.len(), tracker.options.max_cycle_size)?;
                tracker.add_work(atoms.len())?;
                rings.push(atom_ring_to_ring(&atoms, &graph, tracker)?);
            }
        }
    }
    Ok(rings)
}

fn uncovered_ring_bonds(
    mol: &Molecule,
    membership: &RingMembership,
    rings: &[Ring],
) -> Vec<BondId> {
    let covered = rings
        .iter()
        .flat_map(|ring| ring.bonds.iter().copied())
        .collect::<BTreeSet<_>>();
    mol.bond_ids()
        .filter(|bond| membership.bond_in_ring(*bond) && !covered.contains(bond))
        .collect()
}

struct RingWorkTracker {
    options: RingPerceptionOptions,
    candidate_cycles: usize,
    equivalent_shortest_paths: usize,
    path_expansions: usize,
    total_work: usize,
}

impl RingWorkTracker {
    fn new(
        options: RingPerceptionOptions,
        atom_count: usize,
        bond_count: usize,
    ) -> std::result::Result<Self, RingPerceptionError> {
        let total_work = atom_count.saturating_add(bond_count);
        let tracker = Self {
            options,
            candidate_cycles: 0,
            equivalent_shortest_paths: 0,
            path_expansions: 0,
            total_work,
        };
        tracker.check("atoms", atom_count, options.max_atoms)?;
        tracker.check("bonds", bond_count, options.max_bonds)?;
        tracker.check("total work", total_work, options.max_total_work)?;
        Ok(tracker)
    }

    fn check(
        &self,
        resource: &'static str,
        observed: usize,
        limit: usize,
    ) -> std::result::Result<(), RingPerceptionError> {
        if observed > limit {
            Err(RingPerceptionError::ResourceLimit {
                resource,
                observed,
                limit,
            })
        } else {
            Ok(())
        }
    }

    fn add_work(&mut self, amount: usize) -> std::result::Result<(), RingPerceptionError> {
        self.total_work = self.total_work.saturating_add(amount);
        self.check("total work", self.total_work, self.options.max_total_work)
    }

    fn record_path_expansion(&mut self) -> std::result::Result<(), RingPerceptionError> {
        self.path_expansions = self.path_expansions.saturating_add(1);
        self.check(
            "path expansions",
            self.path_expansions,
            self.options.max_path_expansions,
        )?;
        self.add_work(1)
    }

    fn record_shortest_path(&mut self) -> std::result::Result<(), RingPerceptionError> {
        self.equivalent_shortest_paths = self.equivalent_shortest_paths.saturating_add(1);
        self.check(
            "equivalent shortest paths",
            self.equivalent_shortest_paths,
            self.options.max_equivalent_shortest_paths,
        )?;
        self.add_work(1)
    }

    fn record_candidate(&mut self) -> std::result::Result<(), RingPerceptionError> {
        self.candidate_cycles = self.candidate_cycles.saturating_add(1);
        self.check(
            "candidate cycles",
            self.candidate_cycles,
            self.options.max_candidates,
        )?;
        self.add_work(1)
    }
}

fn sssr_bond_counts(sssr: &[Ring], bond_slots: usize) -> Vec<usize> {
    let mut bond_counts = vec![0usize; bond_slots];
    for ring in sssr {
        for bond in &ring.bonds {
            bond_counts[bond.index()] = bond_counts[bond.index()].saturating_add(1);
        }
    }
    bond_counts
}

fn symmetrize_ring_set(
    mut rings: Vec<Ring>,
    extras: Vec<Ring>,
    bond_slots: usize,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<Vec<Ring>, RingPerceptionError> {
    let basis_len = rings.len();
    let bond_counts = sssr_bond_counts(&rings, bond_slots);
    let mut selected_bonds = rings
        .iter()
        .map(|ring| ring.bonds.clone())
        .collect::<BTreeSet<_>>();
    for extra in extras {
        // Only the original SSSR can witness a replacement. Previously added
        // symmetric rings must not recursively justify further additions.
        for ring in &rings[..basis_len] {
            tracker.add_work(1)?;
            if ring.bonds.len() != extra.bonds.len() {
                continue;
            }
            tracker.add_work(ring.bonds.len().saturating_mul(extra.bonds.len()))?;
            if can_replace_one_sssr_ring(&extra, ring, &bond_counts) {
                if selected_bonds.insert(extra.bonds.clone()) {
                    rings.push(extra);
                }
                break;
            }
        }
    }
    Ok(rings)
}

fn can_replace_one_sssr_ring(extra: &Ring, ring: &Ring, bond_counts: &[usize]) -> bool {
    if ring.bonds.len() != extra.bonds.len() {
        return false;
    }
    let mut shares_bond = false;
    for bond in &ring.bonds {
        let included = extra.bonds.contains(bond);
        shares_bond |= included;
        if bond_counts[bond.index()] == 1 && !included {
            return false;
        }
    }
    shares_bond
}

pub(crate) fn ordered_atom_pair(a: AtomId, b: AtomId) -> (AtomId, AtomId) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

#[cfg(test)]
mod tests;
