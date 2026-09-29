// Canonical aromatic seed preparation for resonance search. The source graph
// and installed aromaticity are never changed by this private operation.
// Partition refinement follows RDKit new_canon (Copyright 2014 Greg Landrum,
// from Roger Sayle's pseudocode); BSD 3-Clause, see LICENSE-RDKit.
// Aromatic seed preparation follows RDKit Kekulize.cpp, Copyright (C)
// 2001-2021 Greg Landrum and other RDKit contributors, under the same license.
use super::*;
use crate::core::{
    AxisOrientation, DoubleBondOrientation, StereoCarrier, StereoElementKind,
    TetrahedralOrientation,
};
use std::collections::VecDeque;

fn ranks(mol: &Molecule, work: &mut super::search::Work) -> Result<Vec<usize>, ResonanceError> {
    let atoms: Vec<_> = mol.atom_ids().collect();
    let n = mol.graph.atom_slot_count();
    let mut rank = vec![0; n];
    let mut order = atoms.clone();
    let mut active = vec![atoms[0]];
    let mut queued = BTreeSet::from([atoms[0]]);
    loop {
        while let Some(partition) = active.pop() {
            work.charge(
                mol.atom_count()
                    .saturating_mul(1 + mol.stereo_elements().count())
                    + mol.bond_count(),
            )?;
            queued.remove(&partition);
            let offset = rank[partition.index()];
            let end = (offset..order.len())
                .take_while(|&j| rank[order[j].index()] == offset)
                .last()
                .expect("nonempty partition")
                + 1;
            let keys: BTreeMap<_, _> = atoms
                .iter()
                .map(|&id| {
                    let a = mol.atom(id).expect("live atom");
                    let mut neighbors: Vec<_> = mol
                        .incident_bonds(id)
                        .expect("live atom")
                        .map(|(bid, b)| {
                            (
                                if mol.perception().bond_is_aromatic(bid) == Some(true) {
                                    12
                                } else {
                                    match b.order {
                                        BondOrder::Zero => 0,
                                        BondOrder::Single => 1,
                                        BondOrder::Double => 2,
                                        BondOrder::Triple => 3,
                                        BondOrder::Quadruple => 4,
                                        BondOrder::Dative => 17,
                                    }
                                },
                                bond_stereo(mol, bid),
                                rank[b.other_atom(id).index()],
                            )
                        })
                        .collect();
                    neighbors.sort_by(|a, b| b.cmp(a));
                    (
                        id,
                        (
                            rank[id.index()],
                            a.atom_map.unwrap_or(0),
                            neighbors.len(),
                            a.element.atomic_number(),
                            a.isotope.unwrap_or(0),
                            mol.implicit_hydrogens(id).expect("live atom").unwrap_or(0),
                            a.formal_charge as u32,
                            atom_stereo(mol, id, &rank),
                            neighbors,
                        ),
                    )
                })
                .collect();
            order[offset..end].sort_by_key(|a| &keys[a]);
            let mut start = offset;
            let mut changed = Vec::new();
            for i in offset..end {
                if i > offset && keys[&order[i]] != keys[&order[i - 1]] {
                    start = i;
                }
                if start != offset {
                    changed.push(order[i]);
                }
                rank[order[i].index()] = start;
            }
            activate(mol, &changed, &order, &rank, &mut active, &mut queued);
        }
        let tie = order
            .windows(2)
            .position(|w| rank[w[0].index()] == rank[w[1].index()]);
        let Some(i) = tie else {
            return Ok(rank);
        };
        let end = (i + 1..order.len())
            .take_while(|&j| rank[order[j].index()] == rank[order[i].index()])
            .last()
            .expect("tied class");
        rank[order[end].index()] = end;
        activate(mol, &[order[end]], &order, &rank, &mut active, &mut queued);
    }
}

fn bond_stereo(mol: &Molecule, bond: BondId) -> u8 {
    mol.stereo_elements()
        .find_map(|(_, e)| match &e.kind {
            StereoElementKind::Axis(s) if s.axis == bond => Some(match s.orientation {
                None => 1,
                Some(AxisOrientation::Clockwise) => 6,
                Some(AxisOrientation::CounterClockwise) => 7,
            }),
            StereoElementKind::DoubleBond(s) if s.bond == bond => Some(match s.orientation {
                None => 1,
                Some(DoubleBondOrientation::Together) => 4,
                Some(DoubleBondOrientation::Opposite) => 5,
            }),
            _ => None,
        })
        .unwrap_or(0)
}

fn atom_stereo(mol: &Molecule, atom: AtomId, rank: &[usize]) -> (bool, u8) {
    for (_, e) in mol.stereo_elements() {
        if let StereoElementKind::Tetrahedral(s) = &e.kind {
            if s.center != atom {
                continue;
            }
            let Some(o) = s.orientation else {
                return (false, 0);
            };
            let r: Vec<_> = s
                .carriers
                .iter()
                .filter_map(|c| match c {
                    StereoCarrier::Atom(a) => Some(rank[a.index()]),
                    _ => None,
                })
                .collect();
            if r.iter().collect::<BTreeSet<_>>().len() != r.len() {
                return (true, 0);
            }
            let odd = r
                .iter()
                .enumerate()
                .map(|(i, a)| r[i + 1..].iter().filter(|b| a > *b).count())
                .sum::<usize>()
                % 2
                != 0;
            return (
                true,
                if (o == TetrahedralOrientation::Clockwise) ^ odd {
                    2
                } else {
                    1
                },
            );
        }
    }
    (false, 0)
}

fn activate(
    mol: &Molecule,
    changed: &[AtomId],
    order: &[AtomId],
    rank: &[usize],
    active: &mut Vec<AtomId>,
    queued: &mut BTreeSet<AtomId>,
) {
    let touched: BTreeSet<_> = changed
        .iter()
        .flat_map(|&a| mol.neighbors(a).expect("live atom"))
        .map(|a| rank[a.index()])
        .collect();
    for p in touched {
        if p + 1 < order.len() && rank[order[p + 1].index()] == p && queued.insert(order[p]) {
            active.push(order[p]);
        }
    }
}

#[derive(Clone)]
struct Walk {
    demand: BTreeSet<AtomId>,
    done: BTreeSet<AtomId>,
    queue: VecDeque<AtomId>,
    doubles: BTreeSet<BondId>,
}

pub(super) fn orders(
    mol: &Molecule,
    work: &mut super::search::Work,
) -> Result<BTreeMap<BondId, BondOrder>, ResonanceError> {
    let aromatic: BTreeSet<_> = mol
        .bonds()
        .filter_map(|(id, _)| (mol.perception().bond_is_aromatic(id) == Some(true)).then_some(id))
        .collect();
    let mut result: BTreeMap<_, _> = mol.bonds().map(|(id, b)| (id, b.order)).collect();
    if aromatic.is_empty() {
        return Ok(result);
    }
    work.charge(
        mol.atom_count()
            .saturating_mul(mol.atom_count() + mol.bond_count()),
    )?;
    let rank = ranks(mol, work)?;
    let mut all = BTreeSet::new();
    let mut demand = BTreeSet::new();
    for &id in &aromatic {
        let b = mol.bond(id).expect("live bond");
        all.extend([b.a(), b.b()]);
        if b.order == BondOrder::Double {
            demand.extend([b.a(), b.b()]);
        }
        // An aromatic flag can coexist with an explicit triple bond (arynes).
        // Canonical Kekulization redistributes only the single/double system.
        if matches!(b.order, BondOrder::Single | BondOrder::Double) {
            result.insert(id, BondOrder::Single);
        }
    }
    let mut all: Vec<_> = all.into_iter().collect();
    all.sort_by_key(|a| rank[a.index()]);
    let mut stack = vec![Walk {
        demand,
        done: BTreeSet::new(),
        queue: VecDeque::new(),
        doubles: BTreeSet::new(),
    }];
    while let Some(mut state) = stack.pop() {
        work.charge(mol.atom_count() + mol.bond_count())?;
        let mut failed = false;
        while state.done.len() < all.len() {
            let a = state.queue.pop_front().unwrap_or_else(|| {
                *all.iter()
                    .find(|a| !state.done.contains(a))
                    .expect("unvisited atom")
            });
            state.done.insert(a);
            let mut neighbors: Vec<_> = mol
                .incident_bonds(a)
                .expect("live atom")
                .filter(|(id, b)| aromatic.contains(id) && !state.done.contains(&b.other_atom(a)))
                .map(|(id, b)| (id, b.other_atom(a)))
                .collect();
            neighbors.sort_by_key(|(_, a)| rank[a.index()]);
            for &(_, n) in &neighbors {
                if !state.queue.contains(&n) {
                    state.queue.push_back(n);
                }
            }
            if !state.demand.remove(&a) {
                continue;
            }
            let choices: Vec<_> = neighbors
                .into_iter()
                .filter(|(_, n)| state.demand.contains(n))
                .collect();
            if choices.is_empty() {
                failed = true;
                break;
            }
            for &(b, n) in choices.iter().skip(1).rev() {
                work.charge(mol.atom_count() + mol.bond_count())?;
                let mut fork = state.clone();
                fork.demand.remove(&n);
                fork.doubles.insert(b);
                stack.push(fork);
            }
            state.demand.remove(&choices[0].1);
            state.doubles.insert(choices[0].0);
        }
        if !failed && state.demand.is_empty() {
            for b in state.doubles {
                result.insert(b, BondOrder::Double);
            }
            return Ok(result);
        }
    }
    Err(ResonanceError::InvalidElectronCount)
}
