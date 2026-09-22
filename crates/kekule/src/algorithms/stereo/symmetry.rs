//! Exact, bounded orientation-reversing graph automorphisms.
//!
//! Refinement only rules out mappings; equal colors are never a symmetry proof.
//! Hydrogen expansion is detached and preserves the owner's IDs and declarations.
use super::eligibility::AtomIdentity;
use super::*;
use std::collections::{BTreeMap, BTreeSet};
use std::result::Result;

const UNMAPPED: usize = usize::MAX;
const LONE_PAIR: usize = usize::MAX - 1;

struct AnalysisGraph {
    adjacency: Vec<Vec<(usize, u8)>>,
    colors: Vec<usize>,
    atoms: Vec<usize>,
    hydrogens: Vec<Option<usize>>,
    anchors: Vec<usize>,
}

fn compressed<T: Ord + Clone>(values: &[T]) -> Vec<usize> {
    let mut ordered = values.to_vec();
    ordered.sort();
    ordered.dedup();
    values
        .iter()
        .map(|value| ordered.binary_search(value).expect("collected signature"))
        .collect()
}

impl AnalysisGraph {
    fn new(
        mol: &Molecule,
        options: StereoPerceptionOptions,
    ) -> Result<Self, StereoPerceptionError> {
        let size = mol.atom_count()
            + mol
                .atom_ids()
                .map(|id| usize::from(atom_hydrogen_count(mol, id)))
                .sum::<usize>();
        if size > options.max_atoms {
            return Err(StereoPerceptionError::ResourceLimit {
                resource: "analysis atoms",
                observed: size,
                limit: options.max_atoms,
            });
        }
        let mut atoms = vec![UNMAPPED; mol.graph.atom_slot_count()];
        let mut signatures = Vec::with_capacity(size);
        for (id, atom) in mol.atoms() {
            atoms[id.index()] = signatures.len();
            signatures.push((
                AtomIdentity::of(atom),
                mol.atom_is_aromatic(id).ok().flatten() == Some(true),
            ));
        }
        let mut adjacency = vec![Vec::new(); size];
        for (id, bond) in mol.bonds() {
            let (a, b) = (atoms[bond.a().index()], atoms[bond.b().index()]);
            let order = if matches!(bond.order, BondOrder::Single | BondOrder::Double)
                && mol.bond_is_aromatic(id).ok().flatten() == Some(true)
            {
                7
            } else {
                match bond.order {
                    BondOrder::Zero => 0,
                    BondOrder::Single => 1,
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    BondOrder::Quadruple => 4,
                    BondOrder::Dative => 5,
                }
            };
            adjacency[a].push((b, order));
            adjacency[b].push((a, if order == 5 { 6 } else { order }));
        }
        let mut hydrogens = vec![None; atoms.len()];
        for id in mol.atom_ids() {
            for _ in 0..atom_hydrogen_count(mol, id) {
                let hydrogen = signatures.len();
                hydrogens[id.index()].get_or_insert(hydrogen);
                signatures.push((AtomIdentity::hydrogen(), false));
                adjacency[atoms[id.index()]].push((hydrogen, 1));
                adjacency[hydrogen].push((atoms[id.index()], 1));
            }
        }
        // Unsupported axis/group dependencies must not supply a false proof.
        // Holding their explicit reference frame fixed is conservative.
        let mut anchors = Vec::new();
        for (_, element) in mol.stereo_elements() {
            if let StereoElementKind::Tetrahedral(stereo) = &element.kind {
                if unclassified_tetrahedral_geometry(mol, stereo.center) {
                    anchors.push(atoms[stereo.center.index()]);
                    anchors.extend(stereo.carriers.iter().filter_map(|carrier| match carrier {
                        StereoCarrier::Atom(id) => Some(atoms[id.index()]),
                        StereoCarrier::ImplicitHydrogen => hydrogens[stereo.center.index()],
                        StereoCarrier::ImplicitLonePair => None,
                    }));
                }
            }
            if let StereoElementKind::Axis(axis) = &element.kind {
                let bond = mol.bond(axis.axis).expect("validated axis");
                anchors.extend([atoms[bond.a().index()], atoms[bond.b().index()]]);
                anchors.extend(axis.carriers.iter().filter_map(|carrier| match carrier {
                    StereoCarrier::Atom(id) => Some(atoms[id.index()]),
                    _ => None,
                }));
            }
        }
        anchors.sort_unstable();
        anchors.dedup();
        Ok(Self {
            adjacency,
            colors: compressed(&signatures),
            atoms,
            hydrogens,
            anchors,
        })
    }

    fn carrier(&self, center: AtomId, carrier: StereoCarrier) -> Option<usize> {
        match carrier {
            StereoCarrier::Atom(id) => Some(self.atoms[id.index()]),
            StereoCarrier::ImplicitHydrogen => self.hydrogens[center.index()],
            StereoCarrier::ImplicitLonePair => Some(LONE_PAIR),
        }
    }

    fn refined(&self, sites: &[Site], active: &[bool], focus: usize) -> Vec<usize> {
        let mut marks = vec![Vec::new(); self.colors.len()];
        for (i, site) in sites.iter().enumerate().filter(|(i, _)| active[*i]) {
            let mark = if i == focus || site.fixed {
                i + 3
            } else {
                site.focus.len()
            };
            for &atom in &site.focus {
                marks[atom].push(mark);
            }
        }
        for (i, &atom) in self.anchors.iter().enumerate() {
            marks[atom].push(sites.len() + i + 3);
        }
        let mut colors = compressed(&self.colors.iter().copied().zip(marks).collect::<Vec<_>>());
        loop {
            let signatures: Vec<_> = self
                .adjacency
                .iter()
                .enumerate()
                .map(|(i, neighbors)| {
                    let mut environment: Vec<_> = neighbors
                        .iter()
                        .map(|&(j, bond)| (bond, colors[j]))
                        .collect();
                    environment.sort_unstable();
                    (colors[i], environment)
                })
                .collect();
            let next = compressed(&signatures);
            if next == colors {
                return colors;
            }
            colors = next;
        }
    }
}

struct Site {
    focus: Vec<usize>,
    carriers: Vec<Vec<usize>>,
    parity: Option<bool>,
    fixed: bool,
    group: Option<(StereoGroupId, StereoGroupKind)>,
}

fn odd_permutation(from: &[usize], to: &[usize]) -> Option<bool> {
    let permutation: Vec<_> = from
        .iter()
        .map(|atom| to.iter().position(|other| atom == other))
        .collect::<Option<_>>()?;
    Some(
        permutation
            .iter()
            .enumerate()
            .map(|(i, value)| {
                permutation[i + 1..]
                    .iter()
                    .filter(|other| *other < value)
                    .count()
            })
            .sum::<usize>()
            % 2
            != 0,
    )
}

impl Site {
    fn new(mol: &Molecule, graph: &AnalysisGraph, candidate: &StereoCandidate) -> Self {
        let (focus, carriers, element) = match candidate {
            StereoCandidate::Tetrahedral { center, carriers } => (
                vec![graph.atoms[center.index()]],
                vec![carriers.iter().map(|&c| graph.carrier(*center, c).expect("eligible carrier")).collect()],
                mol.stereo_elements().find(|(_, element)| matches!(&element.kind, StereoElementKind::Tetrahedral(s) if s.center == *center)),
            ),
            StereoCandidate::DoubleBond { bond, left, right, left_carriers, right_carriers } => {
                let carriers = [(*left, left_carriers), (*right, right_carriers)].into_iter().map(|(center, carriers)| {
                    let mut nodes: Vec<_> = carriers.iter().map(|&c| graph.carrier(center, c).expect("eligible carrier")).collect();
                    if nodes.len() == 1 { nodes.push(LONE_PAIR); }
                    nodes
                }).collect();
                (vec![graph.atoms[left.index()], graph.atoms[right.index()]], carriers,
                 mol.stereo_elements().find(|(_, element)| matches!(&element.kind, StereoElementKind::DoubleBond(s) if s.bond == *bond)))
            }
        };
        let mut site = Self {
            focus,
            carriers,
            parity: None,
            fixed: true,
            group: None,
        };
        if let Some((id, element)) = element {
            site.parity = site.represented_parity(graph, element);
            site.group = mol.stereo_groups().find_map(|(group_id, group)| {
                (group.kind != StereoGroupKind::Absolute && group.members.contains(&id))
                    .then_some((group_id, group.kind))
            });
            site.fixed = site.parity.is_none();
        }
        site
    }

    fn represented_parity(&self, graph: &AnalysisGraph, element: &StereoElement) -> Option<bool> {
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                let carriers = stereo
                    .carriers
                    .iter()
                    .map(|&c| graph.carrier(stereo.center, c))
                    .collect::<Option<Vec<_>>>()?;
                Some(
                    (stereo.orientation? == TetrahedralOrientation::Clockwise)
                        ^ odd_permutation(&carriers, &self.carriers[0])?,
                )
            }
            StereoElementKind::DoubleBond(stereo) => {
                let mut parity = stereo.orientation? == DoubleBondOrientation::Opposite;
                for (center, carrier) in [
                    (stereo.left, stereo.left_carrier),
                    (stereo.right, stereo.right_carrier),
                ] {
                    let side = self
                        .focus
                        .iter()
                        .position(|&i| i == graph.atoms[center.index()])?;
                    parity ^= self.carriers[side]
                        .iter()
                        .position(|&i| Some(i) == graph.carrier(center, carrier))?
                        != 0;
                }
                Some(parity)
            }
            StereoElementKind::Axis(_) => None,
        }
    }

    fn ready(&self, mapping: &[usize]) -> bool {
        self.focus
            .iter()
            .chain(self.carriers.iter().flatten())
            .all(|&atom| atom == LONE_PAIR || mapping[atom] != UNMAPPED)
    }

    fn mapped_focus_matches(&self, target: &Self, mapping: &[usize]) -> bool {
        self.focus.len() == target.focus.len()
            && (self
                .focus
                .iter()
                .map(|&i| mapping[i])
                .eq(target.focus.iter().copied())
                || self
                    .focus
                    .iter()
                    .map(|&i| mapping[i])
                    .eq(target.focus.iter().rev().copied()))
    }

    fn mapped_parity(&self, target: &Self, mapping: &[usize]) -> Option<bool> {
        let reverse = mapping[self.focus[0]] != target.focus[0];
        let mut odd = false;
        for (side, carriers) in self.carriers.iter().enumerate() {
            let mapped: Vec<_> = carriers
                .iter()
                .map(|&i| if i == LONE_PAIR { i } else { mapping[i] })
                .collect();
            let target_side = if reverse { 1 - side } else { side };
            odd ^= odd_permutation(&mapped, &target.carriers[target_side])?;
        }
        Some(odd)
    }
}

pub(super) fn filter(
    mol: &Molecule,
    candidates: Vec<StereoCandidate>,
    options: StereoPerceptionOptions,
) -> Result<Vec<StereoCandidate>, StereoPerceptionError> {
    let graph = AnalysisGraph::new(mol, options)?;
    let mut sites: Vec<_> = candidates
        .iter()
        .map(|candidate| Site::new(mol, &graph, candidate))
        .collect();
    // Keep a group fixed when some of its members are outside the classified
    // candidate families. A partial group must never prove a false symmetry.
    for (id, group) in mol.stereo_groups() {
        let classified = sites
            .iter()
            .filter(|s| s.group.is_some_and(|(g, _)| g == id) && s.parity.is_some())
            .count();
        if group.kind != StereoGroupKind::Absolute && classified != group.members.len() {
            for site in &mut sites {
                if site.group.is_some_and(|(g, _)| g == id) {
                    site.fixed = true;
                }
            }
        }
    }
    let mut active = vec![true; sites.len()];
    let mut work = 0;
    loop {
        let mut removed = Vec::new();
        for focus in (0..sites.len()).filter(|&i| active[i]) {
            let colors = graph.refined(&sites, &active, focus);
            let site = &sites[focus];
            let can_reverse = site.carriers.iter().any(|carriers| {
                let mut ranks: Vec<_> = carriers
                    .iter()
                    .map(|&i| if i == LONE_PAIR { LONE_PAIR } else { colors[i] })
                    .collect();
                ranks.sort_unstable();
                ranks.windows(2).any(|pair| pair[0] == pair[1])
            });
            if can_reverse
                && Search::new(&graph, &sites, &active, focus, colors)
                    .reverses(&mut work, options.max_search_states)?
            {
                removed.push(focus);
            }
        }
        if removed.is_empty() {
            break;
        }
        for i in removed {
            active[i] = false;
        }
    }
    Ok(candidates
        .into_iter()
        .zip(active)
        .filter_map(|(candidate, active)| active.then_some(candidate))
        .collect())
}

struct Search<'a> {
    graph: &'a AnalysisGraph,
    sites: &'a [Site],
    active: &'a [bool],
    focus: usize,
    colors: Vec<usize>,
    classes: BTreeMap<usize, Vec<usize>>,
    mapping: Vec<usize>,
    used: Vec<bool>,
    priority: Vec<u8>,
    twins: Vec<usize>,
}

struct Frame {
    atom: usize,
    choices: Vec<usize>,
    next: usize,
}

impl<'a> Search<'a> {
    fn new(
        graph: &'a AnalysisGraph,
        sites: &'a [Site],
        active: &'a [bool],
        focus: usize,
        colors: Vec<usize>,
    ) -> Self {
        let mut classes = BTreeMap::<_, Vec<_>>::new();
        for (atom, &color) in colors.iter().enumerate() {
            classes.entry(color).or_default().push(atom);
        }
        let mut mapping = vec![UNMAPPED; colors.len()];
        let mut used = vec![false; colors.len()];
        for class in classes.values().filter(|class| class.len() == 1) {
            mapping[class[0]] = class[0];
            used[class[0]] = true;
        }
        let mut priority = vec![2; colors.len()];
        for (i, site) in sites.iter().enumerate().filter(|(i, _)| active[*i]) {
            for &atom in site.focus.iter().chain(site.carriers.iter().flatten()) {
                if atom != LONE_PAIR {
                    priority[atom] = priority[atom].min(u8::from(i != focus));
                }
            }
        }
        for &atom in &graph.anchors {
            priority[atom] = priority[atom].min(1);
        }
        // Swapping same-colored vertices with identical typed neighborhoods is
        // an automorphism. Exclude all stereo reference vertices: only then can
        // these swaps be composed with a solution without changing its parity.
        let twins = compressed(
            &graph
                .adjacency
                .iter()
                .enumerate()
                .map(|(i, neighbors)| {
                    let mut neighbors = neighbors.clone();
                    neighbors.sort_unstable();
                    (colors[i], neighbors, (priority[i] < 2).then_some(i))
                })
                .collect::<Vec<_>>(),
        );
        Self {
            graph,
            sites,
            active,
            focus,
            colors,
            classes,
            mapping,
            used,
            priority,
            twins,
        }
    }

    fn stereo_matches(&self) -> bool {
        let mut groups = BTreeMap::new();
        let mut inverse_groups = BTreeMap::new();
        for (i, site) in self
            .sites
            .iter()
            .enumerate()
            .filter(|(i, _)| self.active[*i])
        {
            if !site.ready(&self.mapping) {
                continue;
            }
            if i == self.focus || site.fixed {
                if !site.mapped_focus_matches(site, &self.mapping)
                    || site.mapped_parity(site, &self.mapping) != Some(i == self.focus)
                {
                    return false;
                }
                if let Some((id, _)) = site.group {
                    if groups
                        .insert(id, (id, false))
                        .is_some_and(|v| v != (id, false))
                        || inverse_groups.insert(id, id).is_some_and(|v| v != id)
                    {
                        return false;
                    }
                }
            } else {
                let Some(target) = self.sites.iter().enumerate().find_map(|(j, target)| {
                    (self.active[j]
                        && !target.fixed
                        && site.mapped_focus_matches(target, &self.mapping))
                    .then_some(target)
                }) else {
                    return false;
                };
                let Some(parity) = site.mapped_parity(target, &self.mapping) else {
                    return false;
                };
                let delta = site.parity.unwrap() ^ parity ^ target.parity.unwrap();
                match (site.group, target.group) {
                    (None, None) if !delta => {}
                    (Some((source, kind)), Some((dest, target_kind))) if kind == target_kind => {
                        if groups
                            .insert(source, (dest, delta))
                            .is_some_and(|v| v != (dest, delta))
                            || inverse_groups
                                .insert(dest, source)
                                .is_some_and(|v| v != source)
                        {
                            return false;
                        }
                    }
                    _ => return false,
                }
            }
        }
        true
    }

    fn feasible(&self, atom: usize, target: usize) -> bool {
        let mut source_open = Vec::new();
        for &(neighbor, bond) in &self.graph.adjacency[atom] {
            if self.mapping[neighbor] == UNMAPPED {
                source_open.push((bond, self.colors[neighbor]));
            } else if !self.graph.adjacency[target].contains(&(self.mapping[neighbor], bond)) {
                return false;
            }
        }
        let mut target_open: Vec<_> = self.graph.adjacency[target]
            .iter()
            .filter(|(i, _)| !self.used[*i])
            .map(|&(i, bond)| (bond, self.colors[i]))
            .collect();
        source_open.sort_unstable();
        target_open.sort_unstable();
        source_open == target_open
    }

    fn frame(&self) -> Option<Frame> {
        let atom = (0..self.mapping.len())
            .filter(|&i| self.mapping[i] == UNMAPPED)
            .min_by_key(|&i| {
                (
                    self.priority[i],
                    self.classes[&self.colors[i]]
                        .iter()
                        .filter(|&&j| !self.used[j])
                        .count(),
                    std::cmp::Reverse(
                        self.graph.adjacency[i]
                            .iter()
                            .filter(|(j, _)| self.mapping[*j] != UNMAPPED)
                            .count(),
                    ),
                )
            })?;
        let mut seen = BTreeSet::new();
        Some(Frame {
            atom,
            // Unused target twins are interchangeable while fixing every
            // existing assignment, so retain one representative per class.
            choices: self.classes[&self.colors[atom]]
                .iter()
                .copied()
                .filter(|&target| !self.used[target] && seen.insert(self.twins[target]))
                .collect(),
            next: 0,
        })
    }

    fn reverses(mut self, work: &mut usize, limit: usize) -> Result<bool, StereoPerceptionError> {
        if !self.stereo_matches() {
            return Ok(false);
        }
        let Some(first) = self.frame() else {
            return Ok(true);
        };
        let mut stack = vec![first];
        while let Some(frame) = stack.last_mut() {
            let atom = frame.atom;
            if self.mapping[atom] != UNMAPPED {
                self.used[self.mapping[atom]] = false;
                self.mapping[atom] = UNMAPPED;
            }
            let mut assigned = false;
            while frame.next < frame.choices.len() {
                let target = frame.choices[frame.next];
                frame.next += 1;
                if self.used[target] {
                    continue;
                }
                *work += 1;
                if *work > limit {
                    return Err(StereoPerceptionError::ResourceLimit {
                        resource: "search states",
                        observed: *work,
                        limit,
                    });
                }
                if !self.feasible(atom, target) {
                    continue;
                }
                self.mapping[atom] = target;
                self.used[target] = true;
                if self.stereo_matches() {
                    assigned = true;
                    break;
                }
                self.mapping[atom] = UNMAPPED;
                self.used[target] = false;
            }
            if assigned {
                if let Some(next) = self.frame() {
                    stack.push(next);
                } else {
                    return Ok(true);
                }
            } else {
                stack.pop();
            }
        }
        Ok(false)
    }
}
