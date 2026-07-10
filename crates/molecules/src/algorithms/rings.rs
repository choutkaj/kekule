use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;

use super::*;
use crate::core::*;

pub(super) fn compute_ring_membership(mol: &Molecule) -> (RingMembership, usize) {
    let mut graph = vec![Vec::<(AtomId, BondId)>::new(); mol.atoms.len()];
    let mut live_bonds = Vec::new();
    for (bond_id, bond) in mol.bonds() {
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
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

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RingWork {
    pub atom_count: usize,
    pub bond_count: usize,
    pub candidate_cycles: usize,
    pub equivalent_shortest_paths: usize,
    pub path_expansions: usize,
    pub queue_peak: usize,
    pub stack_peak: usize,
    pub total_work: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingPerceptionError {
    pub resource: &'static str,
    pub observed: usize,
    pub limit: usize,
    pub work: RingWork,
}

impl fmt::Display for RingPerceptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "ring perception {} limit exceeded: observed {}, limit {}",
            self.resource, self.observed, self.limit
        )
    }
}

impl std::error::Error for RingPerceptionError {}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingSet {
    rings: Vec<Ring>,
    work: RingWork,
}

impl RingSet {
    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }

    pub fn len(&self) -> usize {
        self.rings.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rings.is_empty()
    }

    pub fn work(&self) -> RingWork {
        self.work
    }
}

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

    let mut candidates = Vec::<Ring>::new();
    let mut seen_candidate_bonds = BTreeSet::<Vec<BondId>>::new();
    for (closing_bond, a, b) in ring_bonds {
        for mut ring in shortest_cycles_excluding(&graph, a, b, closing_bond, &mut tracker)? {
            ring.bonds.push(closing_bond);
            let mut bond_key = ring.bonds.clone();
            bond_key.sort();
            bond_key.dedup();
            ring.bonds = bond_key.clone();
            if seen_candidate_bonds.insert(bond_key) {
                tracker.record_candidate()?;
                candidates.push(ring);
            }
        }
    }

    candidates.sort_by_key(|ring| (ring.bonds.len(), ring.atoms.clone(), ring.bonds.clone()));
    let symmetric_extra_allowed = symmetric_extra_candidates(&candidates);
    let cyclomatic = mol.bond_count().saturating_add(connected_components(mol)) - mol.atom_count();
    let bit_len = mol.bonds.len().div_ceil(64);
    let mut basis_rows = BTreeMap::<usize, Vec<u64>>::new();
    let mut rings = Vec::new();
    let mut selected_bonds = BTreeSet::<Vec<BondId>>::new();
    for ring in &candidates {
        let bits = cycle_bond_bits(&ring.bonds, bit_len);
        if add_independent_cycle(bits, &mut basis_rows) {
            selected_bonds.insert(ring.bonds.clone());
            rings.push(ring.clone());
            if rings.len() == cyclomatic {
                break;
            }
        }
    }
    if !rings.is_empty() {
        for (index, ring) in candidates.into_iter().enumerate() {
            if symmetric_extra_allowed[index] && selected_bonds.insert(ring.bonds.clone()) {
                rings.push(ring);
            }
        }
    }
    let ring_set = RingSet {
        rings,
        work: tracker.work,
    };
    mol.ring_membership = Some(membership);
    mol.ring_set = Some(ring_set.clone());
    mol.perception.rings = ComputedState::Fresh;
    Ok(ring_set)
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
            Err(RingPerceptionError {
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

fn symmetric_extra_candidates(candidates: &[Ring]) -> Vec<bool> {
    let mut components = (0..candidates.len()).collect::<Vec<_>>();
    for left in 0..candidates.len() {
        for right in (left + 1)..candidates.len() {
            if rings_share_atom(&candidates[left], &candidates[right]) {
                union_components(&mut components, left, right);
            }
        }
    }

    let mut component_candidates = BTreeMap::<usize, Vec<usize>>::new();
    for index in 0..candidates.len() {
        let root = find_component(&mut components, index);
        component_candidates.entry(root).or_default().push(index);
    }

    let mut allowed = vec![false; candidates.len()];
    for indexes in component_candidates.values() {
        let mut atoms = BTreeSet::new();
        let mut bonds = BTreeSet::new();
        for index in indexes {
            atoms.extend(candidates[*index].atoms.iter().copied());
            bonds.extend(candidates[*index].bonds.iter().copied());
        }
        let rank = bonds.len().saturating_add(1).saturating_sub(atoms.len());
        if indexes.len() == rank + 1 {
            for index in indexes {
                allowed[*index] = true;
            }
        }
    }
    allowed
}

fn rings_share_atom(left: &Ring, right: &Ring) -> bool {
    left.atoms.iter().any(|atom| right.atoms.contains(atom))
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

fn connected_components(mol: &Molecule) -> usize {
    let mut seen = BTreeMap::<AtomId, ()>::new();
    let mut count = 0;
    for start in mol.atom_ids() {
        if seen.contains_key(&start) {
            continue;
        }
        count += 1;
        let mut stack = vec![start];
        seen.insert(start, ());
        while let Some(atom) = stack.pop() {
            if let Ok(neighbors) = mol.neighbors(atom) {
                for neighbor in neighbors {
                    if seen.insert(neighbor, ()).is_none() {
                        stack.push(neighbor);
                    }
                }
            }
        }
    }
    count
}
