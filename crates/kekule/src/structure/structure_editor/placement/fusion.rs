//! Fusing a ring fragment onto a bond of a structure draft.
use super::{
    fragment_point, perpendicular_part, signed_angle, StructureEditor, CLEARANCE_REACH,
    CLEARANCE_TIE, CLEAR_CONTACT,
};
use crate::math::{Point3, Vector3};
use crate::molecule::{AtomId, BondOrder, Molecule};
use crate::structure::placement::{
    fusion_dropped, fusion_orders, single_heavy_neighbor, unit, Fragment, FusedFragment, FusionEnd,
    FusionOrders, NameAllocator, PlacementError, Rotation,
};
use crate::structure::StructureEditError;
use crate::topology::edit::{EditAtomId, EditBondId, TopologyEditor};
use crate::topology::{FusionEndPlan, FusionPlan};
use std::collections::BTreeMap;

/// The largest ring, in atoms, that counts toward a bond already shared by
/// two rings.
const MAX_SHARED_RING: usize = 8;
/// Placements whose ring bonds lie within this angle, in radians (10°), of
/// the best one's free sites count as equally aligned; the least crowded
/// wins.
const ALIGNMENT_TOLERANCE: f64 = 0.174_532_925_199_432_95;
/// Contacts closer than this, in nm, crowd a placement.
const CROWDING_RANGE: f64 = 0.3;
/// Fragment atoms within this distance, in nm, of the plane a mirror image
/// is taken through do not count as moved by it.
const MIRROR_TOLERANCE: f64 = 1.0e-4;
/// Crowding scores closer than this, in nm², tie; the earlier wins.
const CROWDING_TIE: f64 = 1.0e-6;

/// Where a target atom's new ring bond can go: toward one of its explicit
/// hydrogens, or one of its free sites.
type Site = (Vector3, Option<EditAtomId>);

/// One candidate placement of the fused ring.
struct Candidate {
    /// Which fragment atom merges into each target atom, and the orders.
    orientation: usize,
    rotation: Rotation,
    /// The largest angle between a ring bond of a target atom and the free
    /// site nearest it.
    deviation: f64,
    /// The smallest distance from a new atom to the atoms around it.
    clearance: f64,
    /// How far, summed and squared, new atoms come within
    /// [`CROWDING_RANGE`] of the atoms around them.
    crowding: f64,
    /// The replaced and the further removed hydrogens of each target atom.
    hydrogens: [(Option<EditAtomId>, Vec<EditAtomId>); 2],
}

/// The fragment's side of one orientation: which of its atoms merges into
/// each target atom, their kept heavy neighbors, and the bond orders.
struct Orientation {
    fused: [AtomId; 2],
    neighbors: [AtomId; 2],
    orders: FusionOrders,
    /// The fragment's coordinates, as in the library or mirrored.
    points: BTreeMap<AtomId, Point3>,
}

impl Orientation {
    /// The rotation that turns the fusion bond, first fused atom to second,
    /// along `axis`.
    fn align(&self, axis: Vector3) -> Rotation {
        let [a, b] = self.fused.map(|atom| self.points[&atom]);
        Rotation::between(unit(b - a).expect("fusion atoms are apart"), axis)
    }
}

impl StructureEditor {
    /// Fuses the ring `fragment` onto `bond`: the fragment's
    /// [fusion bond](Fragment::fusion_bond) is laid onto it, its two atoms
    /// merge into the bond's atoms, and the rest of the ring is added, so
    /// benzene on a bond of benzene gives naphthalene and cyclohexane on a
    /// bond of cyclohexane decalin.
    ///
    /// Bond orders come first. The double bonds of the fragment's Kekulé
    /// structure pair its π atoms again over the fused ring, together with
    /// each bond atom whose fragment atom is π and that has no other multiple
    /// bond; the pairing that includes the most bond atoms wins, then the one
    /// that keeps the bond's order, then the one that keeps the most fragment
    /// double bonds. A fused benzene ring is therefore aromatic on an
    /// aromatic bond whatever its Kekulé order, and turns a cyclohexane bond
    /// into the aromatic bond of tetralin when both its atoms have two
    /// hydrogens to give up. No other existing bond changes.
    ///
    /// Each bond atom gives up one hydrogen for each unit of bond valence it
    /// gains: the explicit hydrogen nearest its new ring bond, when that is
    /// nearer than a free site, takes the bond's place in represented stereo
    /// and is deleted; otherwise an implicit hydrogen is taken. Further
    /// hydrogens come from the implicit count first, then from the explicit
    /// hydrogens nearest the new bond. Fixed counts follow as
    /// [`Self::add_bond`] makes them follow. The new atoms keep the
    /// fragment's hydrogens; [`Self::add_hydrogens`] completes nothing here.
    ///
    /// Geometry: the fragment keeps its ideal coordinates, with the midpoint
    /// of its fusion bond on the bond's midpoint and the two bonds aligned.
    /// Its turn about the bond lays each bond atom's new ring bond along one
    /// of that atom's hydrogens or free sites, or halfway between the two
    /// atoms' best ones, so a ring fused onto an aromatic bond lies in the
    /// plane, away from the existing ring. Of the turns that align within 10°
    /// of the best and keep the new atoms 2 Å from the atoms around them, the
    /// one that crowds them least within 3 Å wins; when none keeps clear,
    /// 15° steps are tried too and the clearest is taken. Both ways round
    /// the fusion bond are tried, and for a fragment without represented
    /// stereo also its mirror image through the plane of the fusion bond and
    /// the ring's center, which turns a chair the other way, so cyclohexane
    /// on cyclohexane gives trans-decalin. With the same bond orders, the
    /// library's coordinates and the fusion bond's first fragment atom
    /// merging into the bond's first atom win ties. Every new atom joins the
    /// residue of the bond's first atom, or else of its second, under the
    /// next free names.
    ///
    /// # Errors
    ///
    /// Nothing changes on error: [`PlacementError::Edit`] for an invalid
    /// handle, [`PlacementError::NotRingFragment`] for a group,
    /// [`PlacementError::FusionAtHydrogen`],
    /// [`PlacementError::FusionBondOrder`] for a bond neither single nor
    /// double, [`PlacementError::SharedFusionBond`] for a bond already in two
    /// rings of up to eight atoms,
    /// [`PlacementError::FusionAtomWithoutHydrogen`], and
    /// [`PlacementError::DegenerateGeometry`] when the bond's atoms coincide,
    /// with the geometry errors of [`Self::free_directions`].
    pub fn fuse_fragment(
        &mut self,
        bond: EditBondId,
        fragment: &Fragment,
    ) -> Result<FusedFragment, PlacementError> {
        let revision = self.topology.structural_revision;
        let target = self.topology.require_bond(bond)?;
        let fusion = fragment
            .fusion_bond()
            .ok_or(PlacementError::NotRingFragment {
                fragment: fragment.name(),
            })?;
        let ends = [target.a(), target.b()];
        if ends
            .iter()
            .any(|&atom| self.element(atom).atomic_number() == 1)
        {
            return Err(PlacementError::FusionAtHydrogen { bond });
        }
        if !matches!(target.order, BondOrder::Single | BondOrder::Double) {
            return Err(PlacementError::FusionBondOrder {
                bond,
                order: target.order,
            });
        }
        if shared_by_two_rings(&self.topology, ends)? {
            return Err(PlacementError::SharedFusionBond { bond });
        }

        // The hydrogens each bond atom can give up.
        let mut explicit = [Vec::new(), Vec::new()];
        let mut implicit = [0usize; 2];
        let mut states = [FusionEnd {
            hydrogens: 0,
            multiple_elsewhere: false,
        }; 2];
        for (end, &atom) in ends.iter().enumerate() {
            for neighbor in self.topology.neighbors(atom)?.collect::<Vec<_>>() {
                if neighbor != ends[1 - end] && self.is_replaceable_hydrogen(atom, neighbor)? {
                    explicit[end].push(neighbor);
                }
            }
            implicit[end] = self.topology.materializable_hydrogens(atom)?.0;
            let hydrogens = explicit[end].len() + implicit[end];
            if hydrogens == 0 {
                return Err(PlacementError::FusionAtomWithoutHydrogen { bond, atom });
            }
            let mut multiple_elsewhere = false;
            for other in self.topology.incident_bonds(atom)? {
                multiple_elsewhere |=
                    other != bond && self.topology.require_bond(other)?.order.valence() > 1;
            }
            states[end] = FusionEnd {
                hydrogens,
                multiple_elsewhere,
            };
        }

        // Both ways round the fusion bond, in the library's coordinates and,
        // for a fragment without stereo, mirrored through the plane of the
        // fusion bond and the ring's center, which turns a chair the other
        // way.
        let molecule = fragment.molecule();
        let fusion_bond = molecule.bond(fusion).expect("the fusion bond is live");
        let dropped = fusion_dropped(molecule, [fusion_bond.a(), fusion_bond.b()]);
        let kept = molecule
            .atom_ids()
            .filter(|atom| !dropped.contains(atom))
            .collect::<Vec<_>>();
        let library = molecule
            .atom_ids()
            .map(|atom| (atom, fragment_point(fragment, atom)))
            .collect::<BTreeMap<_, _>>();
        let fragment_middle = {
            let [a, b] = [fusion_bond.a(), fusion_bond.b()].map(|atom| library[&atom]);
            a + (b - a) * 0.5
        };
        let mirrored = mirror_plane(
            molecule,
            &library,
            &kept,
            [fusion_bond.a(), fusion_bond.b()],
        )
        .filter(|_| molecule.stereo_elements().next().is_none())
        .map(|normal| {
            library
                .iter()
                .map(|(atom, point)| {
                    (
                        *atom,
                        *point - normal * (2.0 * (*point - fragment_middle).dot(normal)),
                    )
                })
                .collect::<BTreeMap<_, _>>()
        });
        let mut orientations = Vec::with_capacity(4);
        for points in std::iter::once(&library).chain(&mirrored) {
            for fused in [
                [fusion_bond.a(), fusion_bond.b()],
                [fusion_bond.b(), fusion_bond.a()],
            ] {
                orientations.push(Orientation {
                    fused,
                    neighbors: [0, 1].map(|end| {
                        single_heavy_neighbor(molecule, fused[end], fused[1 - end])
                            .expect("fusion atoms have one other heavy neighbor")
                    }),
                    orders: fusion_orders(molecule, fused, target.order, states)
                        .expect("a fusion bond leaves the rest of its ring paired"),
                    points: points.clone(),
                });
            }
        }
        // Bond orders first: the orientations whose pairing ranks best.
        let best_rank = orientations
            .iter()
            .map(|orientation| orientation.orders.rank)
            .max()
            .expect("two orientations");
        orientations.retain(|orientation| orientation.orders.rank == best_rank);

        // Geometry.
        let points = ends.map(|atom| self.point(atom));
        let axis = unit(points[1] - points[0])
            .ok_or(PlacementError::DegenerateGeometry { atom: ends[0] })?;
        let middle = points[0] + (points[1] - points[0]) * 0.5;
        let mut sites: [Vec<Site>; 2] = [Vec::new(), Vec::new()];
        for end in 0..2 {
            for &hydrogen in &explicit[end] {
                if let Some(direction) = unit(self.point(hydrogen) - points[end]) {
                    sites[end].push((direction, Some(hydrogen)));
                }
            }
            if implicit[end] > 0 {
                let free = self.free_directions(ends[end])?;
                sites[end].extend(free.into_iter().map(|direction| (direction, None)));
            }
        }
        let reach = kept
            .iter()
            .map(|atom| (library[atom] - fragment_middle).norm())
            .fold(0.0, f64::max)
            + CLEARANCE_REACH;
        // The atoms bonded to each bond atom, which the atom its new ring
        // bond reaches sits beside by bond angle alone.
        let mut beside: [Vec<(EditAtomId, Point3)>; 2] = [Vec::new(), Vec::new()];
        for end in 0..2 {
            for neighbor in self.topology.neighbors(ends[end])? {
                if neighbor != ends[1 - end] {
                    beside[end].push((neighbor, self.point(neighbor)));
                }
            }
        }
        let excluded = ends
            .iter()
            .copied()
            .chain(beside.iter().flatten().map(|(atom, _)| *atom))
            .map(Some)
            .collect::<Vec<_>>();
        let environment = self.environment(middle, reach, &excluded);
        let hydrogen_points = explicit
            .iter()
            .flatten()
            .map(|hydrogen| (*hydrogen, self.point(*hydrogen)))
            .collect::<BTreeMap<_, _>>();

        let evaluate = |orientation: usize, angle: f64| -> Candidate {
            let side = &orientations[orientation];
            let rotation = Rotation::about(axis, angle).after(&side.align(axis));
            let place =
                |atom: &AtomId| middle + rotation.apply(side.points[atom] - fragment_middle);
            let mut deviation = 0.0f64;
            let mut hydrogens: [(Option<EditAtomId>, Vec<EditAtomId>); 2] = Default::default();
            for end in 0..2 {
                let bond_direction =
                    unit(place(&side.neighbors[end]) - points[end]).unwrap_or(axis);
                let nearest = sites[end]
                    .iter()
                    .max_by(|a, b| a.0.dot(bond_direction).total_cmp(&b.0.dot(bond_direction)));
                deviation = deviation.max(nearest.map_or(std::f64::consts::PI, |site| {
                    site.0.dot(bond_direction).clamp(-1.0, 1.0).acos()
                }));
                hydrogens[end] = given_up(
                    &explicit[end],
                    implicit[end],
                    side.orders.gains[end],
                    nearest.and_then(|site| site.1),
                    |hydrogen| (hydrogen_points[&hydrogen] - points[end]).dot(bond_direction),
                );
            }
            let removed = hydrogens
                .iter()
                .flat_map(|(replaced, removed)| replaced.iter().chain(removed))
                .collect::<Vec<_>>();
            let mut clearance = f64::INFINITY;
            let mut crowding = 0.0;
            for atom in &kept {
                let point = place(atom);
                // The atom a bond atom's new ring bond reaches is not
                // measured against that atom's other neighbors.
                let bonded = side.neighbors.iter().position(|neighbor| neighbor == atom);
                let near = beside
                    .iter()
                    .enumerate()
                    .filter(|(end, _)| Some(*end) != bonded)
                    .flat_map(|(_, atoms)| atoms)
                    .filter(|(atom, _)| !removed.contains(&atom))
                    .map(|(_, point)| point);
                for other in environment.iter().chain(near) {
                    let distance = (point - *other).norm();
                    clearance = clearance.min(distance);
                    crowding += (CROWDING_RANGE - distance).max(0.0).powi(2);
                }
            }
            Candidate {
                orientation,
                rotation,
                deviation,
                clearance,
                crowding,
                hydrogens,
            }
        };

        // Turns that lay a ring bond along a site, or halfway between two.
        let mut aligned = Vec::new();
        for (index, side) in orientations.iter().enumerate() {
            let align = side.align(axis);
            let exact = [0, 1].map(|end| {
                let bond =
                    align.apply(side.points[&side.neighbors[end]] - side.points[&side.fused[end]]);
                sites[end]
                    .iter()
                    .map(|(site, _)| {
                        signed_angle(
                            perpendicular_part(bond, axis),
                            perpendicular_part(*site, axis),
                            axis,
                        )
                    })
                    .collect::<Vec<_>>()
            });
            let mut angles = exact.iter().flatten().copied().collect::<Vec<_>>();
            for &first in &exact[0] {
                for &second in &exact[1] {
                    let (sin, cos) = (first.sin() + second.sin(), first.cos() + second.cos());
                    if sin.hypot(cos) > 1.0e-6 {
                        angles.push(sin.atan2(cos));
                    }
                }
            }
            aligned.extend(angles.into_iter().map(|angle| evaluate(index, angle)));
        }
        let chosen = match choose(&aligned) {
            Some(index) => aligned.swap_remove(index),
            None => {
                // No aligned turn keeps clear: try 15° steps as well.
                let mut all = aligned;
                for index in 0..orientations.len() {
                    all.extend(
                        (0..24).map(|step| evaluate(index, f64::from(step).to_radians() * 15.0)),
                    );
                }
                let index = choose(&all).unwrap_or_else(|| clearest(&all));
                all.swap_remove(index)
            }
        };
        let side = &orientations[chosen.orientation];

        // The edit.
        self.slots
            .try_reserve_atoms(kept.len())
            .map_err(StructureEditError::from)?;
        let index_of = kept
            .iter()
            .enumerate()
            .map(|(index, atom)| (*atom, index))
            .chain([(side.fused[0], kept.len()), (side.fused[1], kept.len() + 1)])
            .collect::<BTreeMap<_, _>>();
        let plan = FusionPlan {
            ends: [0, 1].map(|end| FusionEndPlan {
                atom: ends[end],
                replaced: chosen.hydrogens[end].0,
                removed: chosen.hydrogens[end].1.clone(),
            }),
            shared: bond,
            shared_order: side.orders.shared,
            atoms: kept
                .iter()
                .map(|atom| molecule.atom(*atom).expect("live").clone())
                .collect(),
            bonds: side
                .orders
                .bonds
                .iter()
                .map(|(id, order)| {
                    let value = molecule.bond(*id).expect("live");
                    (index_of[&value.a()], index_of[&value.b()], *order)
                })
                .collect(),
        };
        let fused = self.topology.fuse_in_instance(&plan)?;
        self.sync_slots(revision);
        self.slots.clear_atom_rows(&fused.dead_atom_slots);
        self.slots.clear_bond_rows(&fused.dead_bond_slots);
        for (atom, new) in kept.iter().zip(&fused.atoms) {
            self.place_row(
                *new,
                middle + chosen.rotation.apply(side.points[atom] - fragment_middle),
            );
        }
        let parent = match self.topology.atom_site_for_atom(ends[0])? {
            Some(_) => ends[0],
            None => ends[1],
        };
        let mut names = NameAllocator::default();
        for atom in &fused.atoms {
            self.join_residue_of(parent, *atom, &mut names)?;
        }
        Ok(FusedFragment {
            atoms: kept
                .iter()
                .copied()
                .zip(fused.atoms.iter().copied())
                .collect(),
            fused: [0, 1].map(|end| (side.fused[end], ends[end])),
            bond,
            removed: plan
                .ends
                .iter()
                .flat_map(|end| end.replaced.iter().chain(&end.removed).copied())
                .collect(),
        })
    }
}

/// The hydrogens a target atom gives up for `gain` units of new bond
/// valence: the one whose site its new ring bond takes (`nearest`, an
/// explicit hydrogen, or `None` for a free site, which takes an implicit
/// one), then implicit ones, then the explicit ones of `explicit` nearest
/// the new bond by `toward`.
fn given_up<H: Copy + PartialEq>(
    explicit: &[H],
    implicit: usize,
    gain: usize,
    nearest: Option<H>,
    toward: impl Fn(H) -> f64,
) -> (Option<H>, Vec<H>) {
    let replaced = nearest.or_else(|| {
        (implicit == 0)
            .then(|| {
                explicit
                    .iter()
                    .copied()
                    .max_by(|a, b| toward(*a).total_cmp(&toward(*b)))
            })
            .flatten()
    });
    let mut remaining = gain.saturating_sub(1);
    let implicit_left = implicit.saturating_sub(usize::from(replaced.is_none()));
    remaining -= remaining.min(implicit_left);
    let mut others = explicit
        .iter()
        .copied()
        .filter(|hydrogen| Some(*hydrogen) != replaced)
        .collect::<Vec<_>>();
    others.sort_by(|a, b| toward(*b).total_cmp(&toward(*a)));
    others.truncate(remaining);
    (replaced, others)
}

/// The index of the candidate to take among those that keep clear: of
/// those aligned within [`ALIGNMENT_TOLERANCE`] of the best aligned, the
/// least crowded, the earliest on ties; `None` when none keeps clear.
fn choose(candidates: &[Candidate]) -> Option<usize> {
    let clear = |candidate: &&Candidate| candidate.clearance >= CLEAR_CONTACT;
    let best = candidates
        .iter()
        .filter(clear)
        .map(|candidate| candidate.deviation)
        .fold(f64::INFINITY, f64::min);
    let mut chosen: Option<usize> = None;
    for (index, candidate) in candidates.iter().enumerate() {
        if candidate.clearance < CLEAR_CONTACT || candidate.deviation > best + ALIGNMENT_TOLERANCE {
            continue;
        }
        if chosen
            .is_none_or(|chosen| candidate.crowding < candidates[chosen].crowding - CROWDING_TIE)
        {
            chosen = Some(index);
        }
    }
    chosen
}

/// The index of the clearest candidate, the earliest on ties.
fn clearest(candidates: &[Candidate]) -> usize {
    let mut chosen = 0;
    for (index, candidate) in candidates.iter().enumerate() {
        if candidate.clearance > candidates[chosen].clearance + CLEARANCE_TIE {
            chosen = index;
        }
    }
    chosen
}

/// The normal of the plane through the fusion bond `fused` and the center
/// of the fragment's `kept` heavy atoms, or `None` when mirroring through it
/// moves no atom, as for a flat ring.
fn mirror_plane(
    molecule: &Molecule,
    points: &BTreeMap<AtomId, Point3>,
    kept: &[AtomId],
    fused: [AtomId; 2],
) -> Option<Vector3> {
    let [a, b] = fused.map(|atom| points[&atom]);
    let axis = unit(b - a)?;
    let middle = a + (b - a) * 0.5;
    let heavy = kept
        .iter()
        .filter(|atom| {
            molecule
                .atom(**atom)
                .is_some_and(|atom| atom.element.atomic_number() != 1)
        })
        .collect::<Vec<_>>();
    let mut center = Vector3::zero();
    for atom in &heavy {
        center += points[*atom] - middle;
    }
    let normal = unit(axis.cross(perpendicular_part(center, axis)))?;
    let moved = points
        .values()
        .any(|point| (*point - middle).dot(normal).abs() > MIRROR_TOLERANCE);
    moved.then_some(normal)
}

/// Whether the bond between `ends` already belongs to two rings of up to
/// [`MAX_SHARED_RING`] atoms: two paths between its atoms, besides the bond,
/// that share no atom.
fn shared_by_two_rings(
    topology: &TopologyEditor,
    ends: [EditAtomId; 2],
) -> Result<bool, PlacementError> {
    fn walk(
        topology: &TopologyEditor,
        path: &mut Vec<EditAtomId>,
        goal: EditAtomId,
        interiors: &mut Vec<Vec<EditAtomId>>,
    ) -> Result<(), PlacementError> {
        let last = *path.last().expect("paths start at an atom");
        for neighbor in topology.neighbors(last)?.collect::<Vec<_>>() {
            if neighbor == goal {
                if path.len() > 1 {
                    interiors.push(path[1..].to_vec());
                }
            } else if path.len() < MAX_SHARED_RING - 1 && !path.contains(&neighbor) {
                path.push(neighbor);
                walk(topology, path, goal, interiors)?;
                path.pop();
            }
        }
        Ok(())
    }
    let mut interiors = Vec::new();
    walk(topology, &mut vec![ends[0]], ends[1], &mut interiors)?;
    Ok(interiors.iter().enumerate().any(|(index, first)| {
        interiors[index + 1..]
            .iter()
            .any(|second| first.iter().all(|atom| !second.contains(atom)))
    }))
}

#[cfg(test)]
mod tests;
