// Rust adaptation of RDKit 2026.03.3 Resonance.cpp, Copyright (C) 2015
// Paolo Tosco. BSD 3-Clause; see LICENSE-RDKit in the crate root.
use super::*;
use crate::algorithms::{aromaticity::rdkit_outer_electrons, valence::explicit_valence};
use std::collections::VecDeque;

pub(super) struct Work {
    left: usize,
    limit: usize,
}
impl Work {
    pub(super) fn charge(&mut self, n: usize) -> Result<(), ResonanceError> {
        self.left = self
            .left
            .checked_sub(n)
            .ok_or(ResonanceError::ResourceLimit { limit: self.limit })?;
        Ok(())
    }
}

struct AtomInfo {
    id: AtomId,
    outer: i32,
    z: usize,
    degree: i32,
    valence: i32,
    charge: i32,
    aromatic: bool,
    adjacent: Vec<(usize, usize)>,
}
struct BondInfo {
    id: BondId,
    ends: [usize; 2],
    order: u8,
}
struct Group {
    atoms: Vec<AtomInfo>,
    bonds: Vec<BondInfo>,
    charge: i32,
    electrons: i32,
    distances: Vec<Vec<i32>>,
}

#[derive(Clone)]
struct State {
    tv: Vec<i32>,
    nb: Vec<i32>,
    fc: Vec<i32>,
    orders: Vec<u8>,
    last: Vec<bool>,
    definitive: Vec<bool>,
    stacked: Vec<bool>,
    bond_done: Vec<bool>,
    todo: Vec<usize>,
    electrons: i32,
    allowed_charge: i32,
}
impl State {
    fn new(g: &Group) -> Self {
        let n = g.atoms.len();
        Self {
            tv: g.atoms.iter().map(|a| a.degree).collect(),
            nb: vec![0; n],
            fc: vec![0; n],
            orders: vec![1; g.bonds.len()],
            last: vec![false; n],
            definitive: vec![false; n],
            stacked: vec![false; n],
            bond_done: vec![false; g.bonds.len()],
            todo: Vec::new(),
            electrons: g.electrons,
            allowed_charge: g.charge,
        }
    }
    fn push(&mut self, a: usize) {
        if !self.stacked[a] {
            self.stacked[a] = true;
            self.todo.push(a);
        }
    }
    fn pop(&mut self) -> Option<usize> {
        while let Some(a) = self.todo.pop() {
            self.stacked[a] = false;
            if !self.definitive[a] {
                return Some(a);
            }
        }
        None
    }
    fn octet(&self, a: usize) -> bool {
        self.nb[a] + self.tv[a] * 2 == 8
    }
    fn needed(&self, a: usize) -> i32 {
        8 - self.tv[a] * 2 - self.nb[a]
    }
    fn neighbor_charged(&self, g: &Group, a: usize, order: i32, outer: i32) -> bool {
        g.atoms[a].adjacent.iter().any(|&(b, n)| {
            ((self.bond_done[b] && !self.octet(n))
                || (!self.bond_done[b] && self.definitive[n] && g.atoms[n].outer < 5 - order))
                && (outer == 0 || g.atoms[n].outer == outer)
        })
    }
    fn allowed(&mut self, g: &Group, a: usize, b: usize, flags: ResonanceFlags) -> u8 {
        self.last[a] |= g.atoms[a]
            .adjacent
            .iter()
            .all(|&(other, _)| other == b || self.bond_done[other]);
        let mut mask = 0;
        for order in 1..=3 {
            let mut allowed = !self.definitive[a] && self.tv[a] <= 5 - order;
            let mut needs_charge = false;
            if allowed && self.last[a] {
                let right = i32::from(g.atoms[a].outer > 4);
                let increment = if right != 0 {
                    i32::from(!self.neighbor_charged(g, a, order, 4))
                } else {
                    let neighbor = self.neighbor_charged(g, a, order, 0);
                    i32::from(
                        !(neighbor || self.allowed_charge == 0 && g.charge != 0)
                            || (neighbor && order == 3 && g.atoms[a].outer < 5),
                    )
                };
                let e = g.atoms[a].outer + self.tv[a] - 1 + order;
                allowed = e + increment + right >= 8;
                if allowed && e < 8 {
                    if g.charge != 0 || flags.bits() & 24 != 0 || right != 0 {
                        needs_charge = true;
                    } else {
                        allowed = false;
                    }
                }
            }
            if allowed {
                mask |= (1 | if needs_charge { 2 } else { 0 }) << ((order - 1) * 2);
            }
        }
        mask
    }
    fn set_order(&mut self, g: &Group, b: usize, order: u8) -> Result<(), ResonanceError> {
        self.electrons -= i32::from(order) * 2;
        if self.electrons < 0 {
            return Err(ResonanceError::InvalidElectronCount);
        }
        self.orders[b] = order;
        self.bond_done[b] = true;
        for a in g.bonds[b].ends {
            self.tv[a] += i32::from(order) - 1;
            if self.last[a] && g.atoms[a].outer < 5 {
                let needed = self.needed(a);
                if needed != 0 && self.allowed_charge != 0 {
                    self.allowed_charge -= (needed / 2) * self.allowed_charge.signum();
                }
            }
        }
        Ok(())
    }
    fn key(&self, atoms: bool, bonds: bool) -> Vec<i32> {
        let mut key = Vec::new();
        if atoms {
            key.extend(
                self.tv
                    .iter()
                    .zip(&self.nb)
                    .map(|(&tv, &nb)| tv | (nb << 4)),
            );
        }
        if bonds {
            key.extend(self.orders.iter().map(|&o| i32::from(o)));
        }
        key
    }
}

#[derive(Clone)]
struct Candidate {
    state: State,
    metrics: [i64; 7],
    cation_right: bool,
    separated: bool,
}
#[derive(Default)]
struct Stats {
    minimum: Option<i64>,
    no_right: bool,
    no_separation: bool,
}
impl Candidate {
    fn metrics_accept(&self, stats: &mut Stats, flags: ResonanceFlags) -> bool {
        let minimum = stats.minimum.get_or_insert(self.metrics[0]);
        *minimum = (*minimum).min(self.metrics[0]);
        if flags.contains(ResonanceFlags::ALLOW_INCOMPLETE_OCTETS) || self.metrics[0] <= *minimum {
            stats.no_right |= !self.cation_right;
            stats.no_separation |= !self.separated;
        }
        (flags.contains(ResonanceFlags::ALLOW_INCOMPLETE_OCTETS) || self.metrics[0] <= *minimum)
            && (flags.contains(ResonanceFlags::UNCONSTRAINED_CATIONS)
                || !self.cation_right
                || !stats.no_right)
            && (flags.contains(ResonanceFlags::ALLOW_CHARGE_SEPARATION)
                || !self.separated
                || !stats.no_separation)
    }
}

fn charges_valid(g: &Group, s: &State, flags: ResonanceFlags) -> (bool, bool, bool) {
    let (
        mut positive,
        mut negative,
        mut right,
        mut incomplete,
        mut nitrogen,
        mut pos_left,
        mut neg_left,
    ) = (false, false, false, false, false, false, false);
    let mut ok = true;
    for (a, info) in g.atoms.iter().enumerate() {
        if !ok {
            break;
        }
        let q = s.fc[a];
        ok = (-2..=1).contains(&q);
        if !ok {
            break;
        }
        positive |= q > 0;
        negative |= q < 0;
        if info.outer > 4 {
            nitrogen |= info.outer == 5;
            incomplete |= !s.octet(a);
            right |= q > 0 && info.outer > 5;
        } else {
            pos_left |= q > 0;
            neg_left |= q < 0;
        }
        ok = !(incomplete && neg_left);
    }
    if pos_left && !flags.contains(ResonanceFlags::UNCONSTRAINED_CATIONS) {
        ok &= g.charge > 0 && !negative && !nitrogen;
    }
    if neg_left && !flags.contains(ResonanceFlags::UNCONSTRAINED_ANIONS) {
        ok &= g.charge < 0 && !positive;
    }
    if !flags.contains(ResonanceFlags::ALLOW_INCOMPLETE_OCTETS) {
        ok &= !(incomplete && negative);
    }
    let mut multiple = vec![false; g.atoms.len()];
    for (b, info) in g.bonds.iter().enumerate() {
        if !ok {
            break;
        }
        for a in info.ends {
            if s.orders[b] > 1 && g.atoms[a].aromatic {
                if multiple[a] {
                    ok = false;
                } else {
                    multiple[a] = true;
                }
            }
            if g.atoms[a].outer < 5 {
                ok &= s.orders[b] == 1 || s.fc[a] < 1;
            }
        }
        let [a, b] = info.ends;
        ok &= !(g.atoms[a].outer < 5 && s.fc[a] != 0 && g.atoms[b].outer < 5 && s.fc[b] != 0);
    }
    (ok, right, positive && negative)
}

fn candidate(g: &Group, state: State, right: bool, separated: bool) -> Candidate {
    const EN: [i32; 89] = [
        1000, 2300, 4160, 912, 1576, 2051, 2544, 3066, 3610, 4193, 4789, 869, 1293, 1613, 1916,
        2253, 2589, 2869, 3242, 734, 1034, 1190, 1380, 1530, 1650, 1750, 1800, 1840, 1880, 1850,
        1590, 1756, 1994, 2211, 2434, 2685, 2966, 706, 963, 1120, 1320, 1410, 1470, 1510, 1540,
        1560, 1590, 1870, 1520, 1656, 1824, 1984, 2158, 2359, 2582, 659, 881, 1000, 1000, 1000,
        1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1090, 1160, 1340, 1470,
        1600, 1650, 1680, 1720, 1920, 1760, 1789, 1854, 2010, 2190, 2390, 2600, 670, 890,
    ];
    let mut m = [0i64; 7];
    for (a, info) in g.atoms.iter().enumerate() {
        let q = state.fc[a];
        m[0] += i64::from(state.needed(a));
        m[1] += i64::from(q.abs());
        m[2] += i64::from(q) * i64::from(EN.get(info.z).copied().unwrap_or(1000));
        if q != 0 {
            m[5] += i64::from(info.id.raw());
            for b in a + 1..g.atoms.len() {
                if state.fc[b] != 0 {
                    m[if q * state.fc[b] > 0 { 3 } else { 4 }] -= i64::from(g.distances[a][b]);
                }
            }
        }
    }
    for (b, info) in g.bonds.iter().enumerate() {
        if state.orders[b] > 1 {
            m[6] += i64::from(info.id.raw());
        }
    }
    Candidate {
        state,
        metrics: m,
        cation_right: right,
        separated,
    }
}

fn store(
    g: &Group,
    mut s: State,
    map: &mut BTreeMap<Vec<i32>, Candidate>,
    stats: &mut Stats,
    flags: ResonanceFlags,
    diversity: &mut BTreeSet<[i64; 5]>,
) {
    for (a, info) in g.atoms.iter().enumerate() {
        s.fc[a] = info.outer - s.tv[a] - s.nb[a];
    }
    let (valid, right, separated) = charges_valid(g, &s, flags);
    if !valid {
        return;
    }
    let key = s.key(true, flags.contains(ResonanceFlags::KEKULE_ALL));
    let c = candidate(g, s, right, separated);
    let previous_stats = (stats.minimum, stats.no_right, stats.no_separation);
    if c.metrics_accept(stats, flags) {
        diversity.insert(c.metrics[..5].try_into().expect("five ranking metrics"));
        map.entry(key).or_insert(c);
    }
    if previous_stats == (stats.minimum, stats.no_right, stats.no_separation) {
        return;
    }
    // Revisit existing candidates when a more complete / less charge-separated
    // arrangement changes the accepted family, as in RDKit's purgeMaps.
    loop {
        let before = (
            stats.minimum,
            stats.no_right,
            stats.no_separation,
            map.len(),
        );
        map.retain(|_, c| c.metrics_accept(stats, flags));
        if before
            == (
                stats.minimum,
                stats.no_right,
                stats.no_separation,
                map.len(),
            )
        {
            break;
        }
    }
    *diversity = map
        .values()
        .map(|c| c.metrics[..5].try_into().expect("five ranking metrics"))
        .collect();
}

fn nonbonded(
    g: &Group,
    s: &State,
    map: &mut BTreeMap<Vec<i32>, Candidate>,
    stats: &mut Stats,
    flags: ResonanceFlags,
    work: &mut Work,
    diversity: &mut BTreeSet<[i64; 5]>,
) -> Result<(), ResonanceError> {
    let slots: Vec<_> = (0..g.atoms.len()).filter(|&a| s.needed(a) > 0).collect();
    let total: i32 = slots.iter().map(|&a| s.needed(a)).sum();
    if total < s.electrons {
        return Ok(());
    }
    let missing = usize::try_from((total - s.electrons + 1) / 2)
        .map_err(|_| ResonanceError::InvalidElectronCount)?;
    if missing > slots.len() {
        return Ok(());
    }
    let mut selected: Vec<_> = (0..missing).collect();
    loop {
        work.charge(g.atoms.len() + g.bonds.len() + 1)?;
        let mut c = s.clone();
        for (i, &a) in slots.iter().enumerate() {
            c.nb[a] = s.needed(a) - if selected.contains(&i) { 2 } else { 0 };
            c.electrons -= c.nb[a];
            if c.nb[a] < 0 || c.electrons < 0 {
                return Err(ResonanceError::InvalidElectronCount);
            }
        }
        store(g, c, map, stats, flags, diversity);
        // Increasing bit-mask order without the upstream 32-bit mask limit.
        let mut pos = 0;
        while pos < missing
            && selected[pos] + 1
                == if pos + 1 < missing {
                    selected[pos + 1]
                } else {
                    slots.len()
                }
        {
            pos += 1;
        }
        if pos == missing {
            break;
        }
        selected[pos] += 1;
        for (i, entry) in selected.iter_mut().enumerate().take(pos) {
            *entry = i;
        }
    }
    Ok(())
}

fn enumerate_group(
    g: &Group,
    options: ResonanceOptions,
    work: &mut Work,
) -> Result<Vec<Candidate>, ResonanceError> {
    let flags = options.flags;
    let all = flags.contains(ResonanceFlags::KEKULE_ALL);
    let empty = State::new(g);
    let mut initial = empty.clone();
    for (a, info) in g.atoms.iter().enumerate() {
        initial.tv[a] = info.valence;
        initial.fc[a] = info.charge;
        initial.nb[a] = info.outer - info.valence - info.charge;
    }
    initial.orders = g.bonds.iter().map(|b| b.order).collect();
    initial.electrons = 0;
    let (_, right, separated) = charges_valid(g, &initial, flags);
    let mut map = BTreeMap::from([(
        initial.key(true, all),
        candidate(g, initial, right, separated),
    )]);
    let mut stats = Stats::default();
    let mut diversity: BTreeSet<[i64; 5]> = map
        .values()
        .map(|c| c.metrics[..5].try_into().expect("five ranking metrics"))
        .collect();
    let mut seen = BTreeSet::new();
    let mut stack = vec![empty];
    while let Some(mut s) = stack.pop() {
        work.charge(g.atoms.len() + g.bonds.len() + 1)?;
        if s.todo.is_empty() {
            s.push(0);
        }
        let mut dead = false;
        while let Some(a) = s.pop() {
            work.charge(1)?;
            let mut next = None;
            for &(b, n) in &g.atoms[a].adjacent {
                if s.definitive[n] || s.bond_done[b] {
                    continue;
                }
                if next.is_none() {
                    next = Some((b, n));
                } else {
                    s.push(n);
                }
            }
            let Some((b, n)) = next else {
                continue;
            };
            let mut mask = 0xff;
            for atom in [a, n] {
                mask &= s.allowed(g, atom, b, flags);
                if s.last[atom] {
                    s.definitive[atom] = true;
                } else {
                    s.push(atom);
                }
            }
            work.charge(g.atoms.len() + g.bonds.len())?;
            let base = s.clone();
            let mut any = false;
            for order in 1u8..=3 {
                let bit = 1 << ((order - 1) * 2);
                if s.electrons >= i32::from(order) * 2
                    && mask & bit != 0
                    && !(mask & (bit << 1) != 0 && (g.atoms[a].outer < 5 || g.atoms[n].outer < 5))
                {
                    if !any {
                        s.set_order(g, b, order)?;
                        any = true;
                    } else {
                        work.charge(g.atoms.len() + g.bonds.len())?;
                        let mut fork = base.clone();
                        fork.set_order(g, b, order)?;
                        stack.push(fork);
                    }
                }
            }
            if !any {
                dead = true;
                break;
            }
        }
        if !dead && seen.insert(s.key(!all, all)) {
            nonbonded(g, &s, &mut map, &mut stats, flags, work, &mut diversity)?;
            if diversity.len() >= options.max_structures {
                break;
            }
        }
    }
    let mut result: Vec<_> = map.into_values().collect();
    result.sort_by_key(|c| c.metrics);
    Ok(result)
}

fn prepare(
    mol: &Molecule,
    group: &ResonanceGroup,
    orders: &BTreeMap<BondId, BondOrder>,
    charge: i32,
    work: &mut Work,
) -> Result<Group, ResonanceError> {
    let index: BTreeMap<_, _> = group
        .atoms
        .iter()
        .enumerate()
        .map(|(i, &a)| (a, i))
        .collect();
    let mut atoms = Vec::new();
    for &id in &group.atoms {
        let atom = mol.atom(id).expect("live group atom");
        let h = mol
            .implicit_hydrogens(id)
            .expect("live atom")
            .ok_or(ResonanceError::UnknownHydrogens(id))?;
        let degree = mol
            .incident_bonds(id)
            .expect("live atom")
            .filter(|(_, b)| b.order != BondOrder::Dative)
            .count()
            + h;
        atoms.push(AtomInfo {
            id,
            outer: i32::from(rdkit_outer_electrons(atom)),
            z: usize::from(atom.element.atomic_number()),
            degree: i32::try_from(degree).map_err(|_| ResonanceError::InvalidElectronCount)?,
            valence: i32::try_from(explicit_valence(mol, id) + h)
                .map_err(|_| ResonanceError::InvalidElectronCount)?,
            charge: i32::from(atom.formal_charge),
            aromatic: mol.perception().atom_is_aromatic(id) == Some(true),
            adjacent: Vec::new(),
        });
    }
    let mut bonds = Vec::new();
    for &id in &group.bonds {
        let b = mol.bond(id).expect("live group bond");
        let order = match orders[&id] {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => return Err(ResonanceError::UnsupportedBond(id)),
        };
        let ends = [index[&b.a()], index[&b.b()]];
        atoms[ends[0]].adjacent.push((bonds.len(), ends[1]));
        atoms[ends[1]].adjacent.push((bonds.len(), ends[0]));
        bonds.push(BondInfo { id, ends, order });
    }
    // Preserve the graph's neighbor traversal order used by the reference DFS.
    for info in &mut atoms {
        let incident: Vec<_> = mol
            .incident_bonds(info.id)
            .expect("live atom")
            .map(|(b, _)| b)
            .collect();
        info.adjacent.sort_by_key(|&(b, _)| {
            incident
                .iter()
                .position(|&id| id == bonds[b].id)
                .expect("incident bond")
        });
    }
    let electrons = bonds.iter().map(|b| i32::from(b.order) * 2).sum::<i32>()
        + atoms
            .iter()
            .map(|a| a.outer - a.valence - a.charge)
            .sum::<i32>();
    if electrons < 0 {
        return Err(ResonanceError::InvalidElectronCount);
    }
    let mut distances = Vec::new();
    for &root in &group.atoms {
        work.charge(mol.graph.atom_slot_count() + mol.graph.bond_slot_count())?;
        let mut d = vec![-1; mol.graph.atom_slot_count()];
        d[root.index()] = 0;
        let mut todo = VecDeque::from([root]);
        while let Some(a) = todo.pop_front() {
            for n in mol.neighbors(a).expect("live atom") {
                if d[n.index()] < 0 {
                    d[n.index()] = d[a.index()] + 1;
                    todo.push_back(n);
                }
            }
        }
        distances.push(group.atoms.iter().map(|a| d[a.index()]).collect());
    }
    Ok(Group {
        atoms,
        bonds,
        charge,
        electrons,
        distances,
    })
}

pub(super) fn enumerate(
    mol: &Molecule,
    groups: &[ResonanceGroup],
    options: ResonanceOptions,
) -> Result<Vec<ResonanceContributor>, ResonanceError> {
    let mut work = Work {
        left: options.max_total_work,
        limit: options.max_total_work,
    };
    if options.max_structures == 0 {
        return Ok(Vec::new());
    }
    work.charge(mol.graph.atom_slot_count() + mol.graph.bond_slot_count())?;
    // The reference counts charge over all conjugated atoms when initializing
    // each group's allowed-charge counter (not over just the current group).
    let charge = groups
        .iter()
        .flat_map(|g| &g.atoms)
        .map(|&a| i32::from(mol.atom(a).expect("live atom").formal_charge))
        .sum();
    let mut families = Vec::new();
    let mut infos = Vec::new();
    let orders = super::kekule::orders(mol, &mut work)?;
    for group in groups {
        let info = prepare(mol, group, &orders, charge, &mut work)?;
        families.push(enumerate_group(&info, options, &mut work)?);
        infos.push(info);
    }
    let mut len = 1usize;
    for family in &families {
        len = len.saturating_mul(family.len()).min(options.max_structures);
    }
    let mut depths: Vec<Vec<usize>> = families
        .iter()
        .map(|f| {
            let mut sizes = Vec::new();
            let mut previous = None;
            for c in f {
                if previous != Some(&c.metrics[..5]) {
                    sizes.push(1);
                    previous = Some(&c.metrics[..5]);
                } else {
                    *sizes.last_mut().expect("depth") += 1;
                }
            }
            sizes
        })
        .collect();
    if len == options.max_structures && !families.is_empty() {
        let mut s = vec![0; families.len()];
        let mut t = s.clone();
        let mut size = 0;
        while size < len {
            work.charge(1)?;
            size = 1;
            for g in 0..families.len() {
                if size >= len {
                    break;
                }
                if s[g] < depths[g].len() {
                    t[g] += depths[g][s[g]];
                    s[g] += 1;
                }
                size = size.saturating_mul(t[g]);
            }
        }
        for g in 0..families.len() {
            depths[g].truncate(s[g]);
            families[g].truncate(t[g]);
        }
    }
    let family_sizes: Vec<_> = families.iter().map(Vec::len).collect();
    let permutations = build_permutations(&family_sizes, &depths, len, &mut work)?;
    let mut result = Vec::new();
    for (_, chosen) in permutations {
        work.charge(mol.atom_count() + mol.bond_count())?;
        let mut c = ResonanceContributor {
            charges: mol.atoms().map(|(id, a)| (id, a.formal_charge)).collect(),
            orders: mol.bonds().map(|(id, b)| (id, b.order)).collect(),
        };
        for (g, &i) in chosen.iter().enumerate() {
            let s = &families[g][i].state;
            for (a, info) in infos[g].atoms.iter().enumerate() {
                c.charges.insert(
                    info.id,
                    i8::try_from(s.fc[a]).map_err(|_| ResonanceError::InvalidElectronCount)?,
                );
            }
            for (b, info) in infos[g].bonds.iter().enumerate() {
                c.orders.insert(
                    info.id,
                    match s.orders[b] {
                        1 => BondOrder::Single,
                        2 => BondOrder::Double,
                        3 => BondOrder::Triple,
                        _ => return Err(ResonanceError::UnsupportedBond(info.id)),
                    },
                );
            }
        }
        result.push(c);
    }
    Ok(result)
}

type Permutation = ((usize, usize, Vec<usize>, Vec<usize>), Vec<usize>);

fn build_permutations(
    family_sizes: &[usize],
    depths: &[Vec<usize>],
    len: usize,
    work: &mut Work,
) -> Result<Vec<Permutation>, ResonanceError> {
    // Reserve work before allocating any permutations: three group-sized
    // vectors, two depth reductions, all depth scans, and an n log n sorting
    // allowance whose comparisons can scan both group-sized key vectors.
    // Checked arithmetic makes an overflowing estimate fail closed as well.
    let estimate = (|| {
        let groups = family_sizes.len();
        let scans = depths
            .iter()
            .try_fold(0usize, |n, d| n.checked_add(d.len()))?;
        let building = groups.checked_mul(5)?.checked_add(1)?.checked_add(scans)?;
        let levels = usize::BITS - len.saturating_sub(1).leading_zeros();
        let sorting = groups
            .checked_mul(2)?
            .checked_add(2)?
            .checked_mul(levels as usize)?;
        len.checked_mul(building.checked_add(sorting)?)
    })()
    .ok_or(ResonanceError::ResourceLimit { limit: work.limit })?;
    work.charge(estimate)?;

    let mut permutations = Vec::new();
    for index in 0..len {
        let mut rest = index;
        let mut chosen = Vec::new();
        let mut ds = Vec::new();
        let mut widths = Vec::new();
        for (g, &size) in family_sizes.iter().enumerate() {
            if size == 0 {
                return Err(ResonanceError::InvalidElectronCount);
            }
            let i = rest % size;
            rest /= size;
            chosen.push(i);
            let mut width = i;
            let mut depth = 0;
            while width >= depths[g][depth] {
                width -= depths[g][depth];
                depth += 1;
            }
            ds.push(depth);
            widths.push(width);
        }
        permutations.push((
            (
                ds.iter().sum::<usize>(),
                ds.iter().max().copied().unwrap_or(0),
                ds,
                widths,
            ),
            chosen,
        ));
    }
    permutations.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(permutations)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn permutation_materialization_and_sorting_require_budget() {
        // Twenty independently resonating groups would produce 2^20 forms.
        let mut work = Work {
            left: 10_000,
            limit: 10_000,
        };
        assert!(matches!(
            build_permutations(&[2; 20], &vec![vec![2]; 20], 1_000_000, &mut work),
            Err(ResonanceError::ResourceLimit { limit: 10_000 })
        ));
        // Construction alone fits, but construction plus sorting does not.
        let mut work = Work {
            left: 60,
            limit: 60,
        };
        assert!(matches!(
            build_permutations(&[2, 2], &[vec![1, 1], vec![1, 1]], 4, &mut work),
            Err(ResonanceError::ResourceLimit { limit: 60 })
        ));
    }

    #[test]
    fn budgeted_permutations_preserve_reference_priority_order() {
        let mut work = Work {
            left: 1_000,
            limit: 1_000,
        };
        let permutations =
            build_permutations(&[2, 2], &[vec![1, 1], vec![1, 1]], 4, &mut work).unwrap();
        assert_eq!(
            permutations
                .into_iter()
                .map(|(_, chosen)| chosen)
                .collect::<Vec<_>>(),
            vec![vec![0, 0], vec![0, 1], vec![1, 0], vec![1, 1]]
        );
        assert!(work.left < work.limit);
        assert!(matches!(
            build_permutations(&[2], &[vec![2]], usize::MAX, &mut work),
            Err(ResonanceError::ResourceLimit { .. })
        ));
    }
}
