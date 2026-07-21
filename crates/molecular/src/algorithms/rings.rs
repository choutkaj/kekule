use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;

use crate::core::*;

pub(super) fn compute_ring_membership(mol: &Molecule) -> (RingMembership, usize) {
    let mut graph = vec![Vec::<(AtomId, BondId)>::new(); mol.atoms.len()];
    let mut live_bonds = Vec::new();
    for (bond_id, bond) in mol.bonds() {
        if matches!(bond.order, BondOrder::Zero | BondOrder::Dative) {
            continue;
        }
        graph[bond.a.index()].push((bond.b, bond_id));
        graph[bond.b.index()].push((bond.a, bond_id));
        live_bonds.push(bond_id);
    }

    let mut discovery = vec![None; mol.atoms.len()];
    let mut low = vec![0usize; mol.atoms.len()];
    let mut bridge = vec![false; mol.bonds.len()];
    let mut time = 0usize;
    let mut stack_peak = 0usize;

    for atom_id in mol.atom_ids().collect::<Vec<_>>() {
        if discovery[atom_id.index()].is_none() {
            stack_peak = stack_peak.max(ring_dfs_iterative(
                atom_id,
                &graph,
                &mut discovery,
                &mut low,
                &mut bridge,
                &mut time,
            ));
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
    (membership, stack_peak)
}

pub(super) fn bond_in_ring_smaller_than(mol: &Molecule, bond_id: BondId, ring_size: usize) -> bool {
    let Ok(bond) = mol.bond(bond_id) else {
        return false;
    };
    if ring_size <= 1 {
        return false;
    }
    let max_path_edges = ring_size - 2;
    let mut seen = vec![false; mol.atoms.len()];
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
        for (next_bond, next) in
            incident.map(|(next_bond, edge)| (next_bond, edge.other_atom(atom)))
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
) -> usize {
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
    let mut stack_peak = stack.len();

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
            stack_peak = stack_peak.max(stack.len());
        } else {
            low[atom.index()] =
                low[atom.index()].min(discovery[neighbor.index()].expect("neighbor discovered"));
        }
    }
    stack_peak
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
    ResourceLimit {
        resource: &'static str,
        observed: usize,
        limit: usize,
        work: RingWork,
    },
    IncompleteRingCoverage {
        uncovered_bonds: Vec<BondId>,
        work: RingWork,
    },
}

impl RingPerceptionError {
    pub fn work(&self) -> RingWork {
        match self {
            Self::ResourceLimit { work, .. } | Self::IncompleteRingCoverage { work, .. } => *work,
        }
    }
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

pub fn perceive_ring_set(mol: &mut Molecule) -> std::result::Result<RingSet, RingPerceptionError> {
    perceive_ring_set_with_options(mol, RingPerceptionOptions::default())
}

pub fn perceive_ring_set_with_options(
    mol: &mut Molecule,
    options: RingPerceptionOptions,
) -> std::result::Result<RingSet, RingPerceptionError> {
    let mut tracker = RingWorkTracker::new(options, mol.atom_count(), mol.bond_count())?;
    let (membership, bridge_stack_peak) = compute_ring_membership(mol);
    tracker.observe_stack(bridge_stack_peak);
    let (mut rings, extras) = figueras_sssr_candidates(mol, &membership, &mut tracker)?;
    if !uncovered_ring_bonds(mol, &membership, &rings).is_empty() {
        rings = complete_ring_coverage(mol, &membership, rings, &mut tracker)?;
    }
    if !rings.is_empty() {
        let bond_counts = sssr_bond_counts(&rings, mol.bonds.len());
        let mut selected_bonds = rings
            .iter()
            .map(|ring| ring.bonds.clone())
            .collect::<BTreeSet<_>>();
        for ring in extras {
            if can_replace_sssr_ring(&ring, &rings, &bond_counts)
                && selected_bonds.insert(ring.bonds.clone())
            {
                rings.push(ring);
            }
        }
    }
    let ring_set = RingSet {
        rings,
        work: tracker.work,
    };
    mol.install_rings(membership, ring_set.clone());
    Ok(ring_set)
}

#[derive(Clone)]
struct ActiveRingGraph {
    adjacency: Vec<Vec<(AtomId, BondId)>>,
    active_bonds: Vec<bool>,
    atom_degrees: Vec<usize>,
}

impl ActiveRingGraph {
    fn new(mol: &Molecule) -> Self {
        let mut adjacency = vec![Vec::new(); mol.atoms.len()];
        let mut active_bonds = vec![false; mol.bonds.len()];
        let mut atom_degrees = vec![0usize; mol.atoms.len()];
        for (bond_id, bond) in mol.bonds() {
            if matches!(bond.order, BondOrder::Zero | BondOrder::Dative) {
                continue;
            }
            adjacency[bond.a.index()].push((bond.b, bond_id));
            adjacency[bond.b.index()].push((bond.a, bond_id));
            active_bonds[bond_id.index()] = true;
            atom_degrees[bond.a.index()] += 1;
            atom_degrees[bond.b.index()] += 1;
        }
        Self {
            adjacency,
            active_bonds,
            atom_degrees,
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
        let mut done = vec![false; mol.atoms.len()];
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

        let (kept, extras) = remove_extra_rings(fragment_candidates, mol.bonds.len());
        all_sssr.extend(kept);
        all_extras.extend(extras);
    }

    // Candidate search uses the full active graph like RDKit. Ring membership is
    // still the authoritative filter for the bounded recovery path below.
    debug_assert!(all_sssr
        .iter()
        .flat_map(|ring| &ring.bonds)
        .all(|bond| membership.bond_in_ring(*bond)));
    Ok((all_sssr, all_extras))
}

fn active_fragments(mol: &Molecule, graph: &ActiveRingGraph) -> Vec<Vec<AtomId>> {
    let mut seen = vec![false; mol.atoms.len()];
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
    while let Some(root) = fragment
        .iter()
        .copied()
        .find(|atom| graph.atom_degrees[atom.index()] == 2 && !forbidden[atom.index()])
    {
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
                candidates.push(atom_ring_to_ring(atoms, graph, tracker)?);
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
            for forbidden_neighbor in included {
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
    let mut colors = vec![WHITE; graph.atom_degrees.len()];
    for atom in forbidden {
        colors[atom.index()] = BLACK;
    }
    let mut parents = vec![None; graph.atom_degrees.len()];
    let mut depths = vec![0usize; graph.atom_degrees.len()];
    let mut queue = VecDeque::from([root]);
    let mut rings = Vec::new();
    let mut current_size = usize::MAX;
    tracker.observe_queue(queue.len());

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
                tracker.observe_queue(queue.len());
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

fn remove_extra_rings(mut rings: Vec<Ring>, bond_slots: usize) -> (Vec<Ring>, Vec<Ring>) {
    rings.sort_by_key(|ring| ring.bonds.len());
    let mut available = vec![true; rings.len()];
    let mut keep = vec![false; rings.len()];
    let mut union = vec![false; bond_slots];

    for index in 0..rings.len() {
        if ring_is_subset_of(&rings[index], &union) {
            available[index] = false;
        }
        if !available[index] {
            continue;
        }
        add_ring_to_union(&rings[index], &mut union);
        keep[index] = true;
        let mut consider = ((index + 1)..rings.len())
            .filter(|other| {
                available[*other] && rings[*other].bonds.len() == rings[index].bonds.len()
            })
            .collect::<BTreeSet<_>>();
        while !consider.is_empty() {
            let mut best = None;
            let mut best_overlap = None;
            for other in consider.iter().copied() {
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
    (kept, extras)
}

fn ring_is_subset_of(ring: &Ring, union: &[bool]) -> bool {
    ring.bonds.iter().all(|bond| union[bond.index()])
}

fn add_ring_to_union(ring: &Ring, union: &mut [bool]) {
    for bond in &ring.bonds {
        union[bond.index()] = true;
    }
}

fn complete_ring_coverage(
    mol: &Molecule,
    membership: &RingMembership,
    mut rings: Vec<Ring>,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<Vec<Ring>, RingPerceptionError> {
    let mut graph = BTreeMap::<AtomId, Vec<(AtomId, BondId)>>::new();
    let mut ring_bonds = Vec::new();
    for (bond_id, bond) in mol.bonds() {
        if membership.bond_in_ring(bond_id) {
            graph.entry(bond.a()).or_default().push((bond.b(), bond_id));
            graph.entry(bond.b()).or_default().push((bond.a(), bond_id));
            ring_bonds.push((bond_id, bond.a(), bond.b()));
        }
    }
    for edges in graph.values_mut() {
        edges.sort_by_key(|(atom, bond)| (*atom, *bond));
    }

    let bit_len = mol.bonds.len().div_ceil(64);
    let mut basis_rows = BTreeMap::<usize, Vec<u64>>::new();
    let mut selected = BTreeSet::<Vec<BondId>>::new();
    rings.retain(|ring| {
        let independent =
            add_independent_cycle(cycle_bond_bits(&ring.bonds, bit_len), &mut basis_rows);
        if independent {
            selected.insert(ring.bonds.clone());
        }
        independent
    });
    if uncovered_ring_bonds(mol, membership, &rings).is_empty() {
        return Ok(rings);
    }
    for (closing_bond, a, b) in ring_bonds.iter().copied() {
        for mut ring in shortest_cycles_excluding(&graph, a, b, closing_bond, tracker)? {
            ring.bonds.push(closing_bond);
            ring.bonds.sort();
            ring.bonds.dedup();
            if selected.contains(&ring.bonds) {
                continue;
            }
            tracker.record_candidate()?;
            if add_independent_cycle(cycle_bond_bits(&ring.bonds, bit_len), &mut basis_rows) {
                selected.insert(ring.bonds.clone());
                rings.push(ring);
                if uncovered_ring_bonds(mol, membership, &rings).is_empty() {
                    return Ok(rings);
                }
            }
        }
    }

    // Edge-local shortest cycles do not necessarily span the complete cycle
    // space in highly bridged cage graphs. Deterministic BFS spanning trees
    // provide a guaranteed fundamental-cycle basis. Consider trees rooted at
    // every ring atom and take the shortest independent fallback, keeping the
    // usual small-ring candidates preferred while guaranteeing completeness.
    for ring in fundamental_cycle_candidates(&graph, &ring_bonds, tracker)? {
        if selected.contains(&ring.bonds) {
            continue;
        }
        if add_independent_cycle(cycle_bond_bits(&ring.bonds, bit_len), &mut basis_rows) {
            selected.insert(ring.bonds.clone());
            rings.push(ring);
            if uncovered_ring_bonds(mol, membership, &rings).is_empty() {
                return Ok(rings);
            }
        }
    }
    Err(RingPerceptionError::IncompleteRingCoverage {
        uncovered_bonds: uncovered_ring_bonds(mol, membership, &rings),
        work: tracker.work,
    })
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

fn fundamental_cycle_candidates(
    graph: &BTreeMap<AtomId, Vec<(AtomId, BondId)>>,
    ring_bonds: &[(BondId, AtomId, AtomId)],
    tracker: &mut RingWorkTracker,
) -> std::result::Result<Vec<Ring>, RingPerceptionError> {
    let atom_slots = graph
        .keys()
        .map(|atom| atom.index())
        .max()
        .map_or(0, |maximum| maximum + 1);
    let mut unique = BTreeMap::<Vec<BondId>, Vec<AtomId>>::new();

    for &root in graph.keys() {
        let mut parent_atom = vec![None; atom_slots];
        let mut parent_bond = vec![None; atom_slots];
        let mut depth = vec![None; atom_slots];
        let mut tree_bonds = BTreeSet::<BondId>::new();
        let mut queue = VecDeque::from([root]);
        depth[root.index()] = Some(0);
        tracker.observe_queue(queue.len());

        while let Some(atom) = queue.pop_front() {
            let atom_depth = depth[atom.index()].expect("queued atom has a depth");
            for (neighbor, bond) in graph.get(&atom).into_iter().flatten().copied() {
                tracker.record_path_expansion()?;
                if depth[neighbor.index()].is_some() {
                    continue;
                }
                depth[neighbor.index()] = Some(atom_depth + 1);
                parent_atom[neighbor.index()] = Some(atom);
                parent_bond[neighbor.index()] = Some(bond);
                tree_bonds.insert(bond);
                queue.push_back(neighbor);
                tracker.observe_queue(queue.len());
            }
        }

        for &(closing_bond, a, b) in ring_bonds {
            if tree_bonds.contains(&closing_bond)
                || depth[a.index()].is_none()
                || depth[b.index()].is_none()
            {
                continue;
            }
            let mut left = a;
            let mut right = b;
            let mut bonds = vec![closing_bond];
            let mut atoms = BTreeSet::from([a, b]);

            while left != right {
                let left_depth = depth[left.index()].expect("tree vertex depth");
                let right_depth = depth[right.index()].expect("tree vertex depth");
                if left_depth >= right_depth {
                    let bond = parent_bond[left.index()].expect("non-root tree vertex bond");
                    left = parent_atom[left.index()].expect("non-root tree vertex parent");
                    bonds.push(bond);
                    atoms.insert(left);
                }
                if right_depth >= left_depth && left != right {
                    let bond = parent_bond[right.index()].expect("non-root tree vertex bond");
                    right = parent_atom[right.index()].expect("non-root tree vertex parent");
                    bonds.push(bond);
                    atoms.insert(right);
                }
            }

            bonds.sort();
            bonds.dedup();
            tracker.check("cycle size", bonds.len(), tracker.options.max_cycle_size)?;
            if bonds.len() >= 3 && !unique.contains_key(&bonds) {
                tracker.record_candidate()?;
                unique.insert(bonds, atoms.into_iter().collect());
            }
        }
    }

    let mut candidates = unique
        .into_iter()
        .map(|(bonds, atoms)| Ring { atoms, bonds })
        .collect::<Vec<_>>();
    candidates.sort_by(|left, right| {
        left.bonds
            .len()
            .cmp(&right.bonds.len())
            .then_with(|| left.bonds.cmp(&right.bonds))
    });
    Ok(candidates)
}

struct RingWorkTracker {
    options: RingPerceptionOptions,
    work: RingWork,
}

impl RingWorkTracker {
    fn new(
        options: RingPerceptionOptions,
        atom_count: usize,
        bond_count: usize,
    ) -> std::result::Result<Self, RingPerceptionError> {
        let work = RingWork {
            atom_count,
            bond_count,
            total_work: atom_count.saturating_add(bond_count),
            ..RingWork::default()
        };
        let tracker = Self { options, work };
        tracker.check("atoms", atom_count, options.max_atoms)?;
        tracker.check("bonds", bond_count, options.max_bonds)?;
        tracker.check("total work", work.total_work, options.max_total_work)?;
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
                work: self.work,
            })
        } else {
            Ok(())
        }
    }

    fn add_work(&mut self, amount: usize) -> std::result::Result<(), RingPerceptionError> {
        self.work.total_work = self.work.total_work.saturating_add(amount);
        self.check(
            "total work",
            self.work.total_work,
            self.options.max_total_work,
        )
    }

    fn record_path_expansion(&mut self) -> std::result::Result<(), RingPerceptionError> {
        self.work.path_expansions = self.work.path_expansions.saturating_add(1);
        self.check(
            "path expansions",
            self.work.path_expansions,
            self.options.max_path_expansions,
        )?;
        self.add_work(1)
    }

    fn record_shortest_path(&mut self) -> std::result::Result<(), RingPerceptionError> {
        self.work.equivalent_shortest_paths = self.work.equivalent_shortest_paths.saturating_add(1);
        self.check(
            "equivalent shortest paths",
            self.work.equivalent_shortest_paths,
            self.options.max_equivalent_shortest_paths,
        )?;
        self.add_work(1)
    }

    fn record_candidate(&mut self) -> std::result::Result<(), RingPerceptionError> {
        self.work.candidate_cycles = self.work.candidate_cycles.saturating_add(1);
        self.check(
            "candidate cycles",
            self.work.candidate_cycles,
            self.options.max_candidates,
        )?;
        self.add_work(1)
    }

    fn observe_queue(&mut self, size: usize) {
        self.work.queue_peak = self.work.queue_peak.max(size);
    }

    fn observe_stack(&mut self, size: usize) {
        self.work.stack_peak = self.work.stack_peak.max(size);
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

fn can_replace_sssr_ring(extra: &Ring, sssr: &[Ring], bond_counts: &[usize]) -> bool {
    sssr.iter()
        .any(|ring| can_replace_one_sssr_ring(extra, ring, bond_counts))
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

pub(super) fn union_components(components: &mut [usize], left: usize, right: usize) {
    let left_root = find_component(components, left);
    let right_root = find_component(components, right);
    if left_root != right_root {
        components[right_root] = left_root;
    }
}

pub(super) fn find_component(components: &mut [usize], index: usize) -> usize {
    let mut root = index;
    while components[root] != root {
        root = components[root];
    }
    let mut current = index;
    while components[current] != current {
        let parent = components[current];
        components[current] = root;
        current = parent;
    }
    root
}

fn shortest_cycles_excluding(
    graph: &BTreeMap<AtomId, Vec<(AtomId, BondId)>>,
    start: AtomId,
    goal: AtomId,
    excluded_bond: BondId,
    tracker: &mut RingWorkTracker,
) -> std::result::Result<Vec<Ring>, RingPerceptionError> {
    let mut queue = VecDeque::new();
    let mut distances = BTreeMap::<AtomId, usize>::new();
    distances.insert(start, 0);
    queue.push_back(start);
    tracker.observe_queue(queue.len());

    while let Some(atom) = queue.pop_front() {
        let distance = distances[&atom];
        for (neighbor, bond_id) in graph.get(&atom).into_iter().flatten().copied() {
            tracker.record_path_expansion()?;
            if bond_id == excluded_bond || distances.contains_key(&neighbor) {
                continue;
            }
            distances.insert(neighbor, distance + 1);
            queue.push_back(neighbor);
            tracker.observe_queue(queue.len());
        }
    }

    let Some(goal_distance) = distances.get(&goal).copied() else {
        return Ok(Vec::new());
    };
    let cycle_size = goal_distance.saturating_add(1);
    tracker.check("cycle size", cycle_size, tracker.options.max_cycle_size)?;

    #[derive(Clone)]
    struct PathState {
        current: AtomId,
        atoms: Vec<AtomId>,
        bonds: Vec<BondId>,
    }

    let mut rings = Vec::new();
    let mut stack = vec![PathState {
        current: goal,
        atoms: vec![goal],
        bonds: Vec::new(),
    }];
    tracker.observe_stack(stack.len());
    while let Some(state) = stack.pop() {
        if state.current == start {
            tracker.record_shortest_path()?;
            rings.push(Ring {
                atoms: state.atoms,
                bonds: state.bonds,
            });
            continue;
        }

        let Some(distance) = distances.get(&state.current).copied() else {
            continue;
        };
        let mut predecessors = graph
            .get(&state.current)
            .into_iter()
            .flatten()
            .copied()
            .filter(|(neighbor, bond_id)| {
                *bond_id != excluded_bond
                    && distances.get(neighbor).copied() == distance.checked_sub(1)
            })
            .collect::<Vec<_>>();
        predecessors.sort_by_key(|(neighbor, bond_id)| (*neighbor, *bond_id));
        for (neighbor, bond_id) in predecessors.into_iter().rev() {
            tracker.record_path_expansion()?;
            if state.atoms.len() >= tracker.options.max_cycle_size {
                continue;
            }
            let mut next = state.clone();
            next.current = neighbor;
            next.atoms.push(neighbor);
            next.bonds.push(bond_id);
            stack.push(next);
            tracker.observe_stack(stack.len());
        }
    }
    Ok(rings)
}

fn cycle_bond_bits(bonds: &[BondId], bit_len: usize) -> Vec<u64> {
    let mut bits = vec![0u64; bit_len];
    for bond in bonds {
        let index = bond.index();
        bits[index / 64] |= 1u64 << (index % 64);
    }
    bits
}

fn add_independent_cycle(mut row: Vec<u64>, basis_rows: &mut BTreeMap<usize, Vec<u64>>) -> bool {
    while let Some(pivot) = first_set_bit(&row) {
        if let Some(existing) = basis_rows.get(&pivot) {
            xor_bits(&mut row, existing);
        } else {
            basis_rows.insert(pivot, row);
            return true;
        }
    }
    false
}

fn first_set_bit(bits: &[u64]) -> Option<usize> {
    bits.iter().enumerate().find_map(|(block, value)| {
        (*value != 0).then_some(block * 64 + value.trailing_zeros() as usize)
    })
}

fn xor_bits(left: &mut [u64], right: &[u64]) {
    for (left, right) in left.iter_mut().zip(right) {
        *left ^= *right;
    }
}

pub(crate) fn ordered_atom_pair(a: AtomId, b: AtomId) -> (AtomId, AtomId) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inconsistent_ring_membership_returns_structured_coverage_error() {
        let mut mol = Molecule::new();
        let atoms = (0..2)
            .map(|_| mol.add_atom(Atom::new(Element::from_symbol("C").expect("carbon"))))
            .collect::<Vec<_>>();
        let bond = mol
            .add_bond(atoms[0], atoms[1], BondOrder::Single)
            .expect("chain bond");
        let (mut membership, _) = compute_ring_membership(&mol);
        membership.atom_flags[atoms[0].index()] = true;
        membership.atom_flags[atoms[1].index()] = true;
        membership.bond_flags[bond.index()] = true;
        let mut tracker = RingWorkTracker::new(
            RingPerceptionOptions::default(),
            mol.atom_count(),
            mol.bond_count(),
        )
        .expect("tracker");

        let error = complete_ring_coverage(&mol, &membership, Vec::new(), &mut tracker)
            .expect_err("a bridge cannot be covered by a cycle");

        assert!(matches!(
            error,
            RingPerceptionError::IncompleteRingCoverage { uncovered_bonds, .. }
                if uncovered_bonds == vec![bond]
        ));
    }
}
