//! Settling atoms and bonds of a structure draft after a chemistry edit:
//! builder hydrogen counts at ideal positions, the angle between an atom's
//! two neighbors, and ideal bond lengths.
use super::{bond_length, draft_error, DomainChange, StructureEditor};
use crate::algorithms::rdkit_hydrogens_at_valence;
use crate::math::Point3;
use crate::molecule::{BondOrder, Element};
use crate::structure::placement::{
    builder_charges, builder_hydrogen_count, free_sites, perpendicular, unit, AtomShape,
    NameAllocator, PlacementError, Rotation, SettledAtom, SettledAtomsReport, SettledBond,
    SettledBondsReport, SiteQuery,
};
use crate::topology::edit::{EditAtomId, EditBondId};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// The tetrahedral angle in degrees, `acos(-1/3)`.
const TETRAHEDRAL: f64 = 109.471_220_634_490_7;
/// Angles within this many degrees of their shape's stay as they are, so a
/// measured geometry keeps its small deviations.
const ANGLE_TOLERANCE: f64 = 5.0;
/// Bond lengths within this, in nm, of their ideal ones stay as they are.
const LENGTH_TOLERANCE: f64 = 1.0e-4;
/// A position within this, in nm, of where an atom is does not move it.
const MOVE_TOLERANCE: f64 = 1.0e-9;

/// What settling one atom needs to know, perceived before anything changes.
#[derive(Debug)]
struct AtomPlan {
    atom: EditAtomId,
    /// Neighbors other than its hydrogens, in neighbor order.
    others: Vec<EditAtomId>,
    /// Its hydrogens: explicit hydrogens singly bonded to it alone.
    hydrogens: Vec<EditAtomId>,
    /// The hydrogens it takes.
    target: usize,
    shape: Option<AtomShape>,
    lone_pairs: usize,
    /// The length of a bond from it to a hydrogen, in nm.
    length: f64,
}

/// Where one atom's hydrogens go.
#[derive(Debug, Default)]
struct HydrogenPlan {
    /// Kept hydrogens and their new positions.
    kept: Vec<(EditAtomId, Point3)>,
    removed: Vec<EditAtomId>,
    /// Positions of the new hydrogens.
    new: Vec<Point3>,
    /// Every hydrogen's position, in order of the atom's sites.
    points: Vec<Point3>,
}

/// Planned positions over the draft's.
type Moves = BTreeMap<EditAtomId, Point3>;

impl StructureEditor {
    /// Settles each distinct atom in `atoms` after an edit changed its
    /// element, charge, or bonds, as a builder does:
    ///
    /// - It takes its builder hydrogen count, from its element, charge, and
    ///   the bond orders to its other neighbors
    ///   ([`builder_hydrogen_count`](crate::structure::placement::builder_hydrogen_count)),
    ///   so the hydrogens it has do not count: methane charged +1 loses one,
    ///   and back at 0 gains it again. Hydrogens are added or deleted to
    ///   reach it, those that stay keep their identity, and none stays
    ///   implicit. Outside the builder's states the valence model gives the
    ///   count from the same bonds.
    /// - When it has two other neighbors whose angle is more than 5° off the
    ///   one of its shape (180° linear, 120° trigonal, 109.47° tetrahedral),
    ///   the smaller side of the two, by atoms, turns about it to that angle:
    ///   a carbon made sp is straight, and a straight one made sp² or sp³
    ///   bends in a deterministic plane. When both sides lie in a ring with
    ///   it, nothing turns, and the entry reports it.
    /// - All its hydrogens move to the ideal positions of its shape, at the
    ///   element pair's bond length, from its other neighbors and lone pairs
    ///   alone ([`placement`](crate::structure::placement)): methyl anion
    ///   pyramidal, methyl cation and borane planar, ammonium tetrahedral.
    ///
    /// The hydrogens of an atom are explicit hydrogens singly bonded to it
    /// alone; a hydrogen that is another atom's is left to it. New hydrogens
    /// join their atom's residue under the next free names. No other atom
    /// moves, and bond lengths stay; [`Self::settle_bonds`] sets them.
    ///
    /// # Errors
    ///
    /// Nothing changes on error: [`PlacementError::Edit`] for an invalid
    /// handle, [`PlacementError::BuilderValence`] for a charge or valence
    /// outside the builder's states, [`PlacementError::InvalidValence`] for
    /// an atom outside them over its valence,
    /// [`PlacementError::UnsupportedGeometry`] for more than six domains,
    /// [`PlacementError::DegenerateGeometry`] when a neighbor coincides with
    /// the atom, [`PlacementError::StereoHydrogen`] when its stereo names an
    /// implicit hydrogen that does not become exactly one new atom, and
    /// [`PlacementError::UnknownBondLength`].
    pub fn settle_atoms(
        &mut self,
        atoms: impl IntoIterator<Item = EditAtomId>,
    ) -> Result<SettledAtomsReport, PlacementError> {
        let mut seen = BTreeSet::new();
        let mut plans = Vec::new();
        for atom in atoms {
            if !seen.insert(atom) {
                continue;
            }
            self.topology.require_atom(atom)?;
            if self.is_hydrogen_of_neighbor(atom)? {
                continue;
            }
            plans.push(self.atom_plan(atom)?);
        }

        let mut moves = Moves::new();
        let mut angles = Vec::with_capacity(plans.len());
        for plan in &plans {
            angles.push(self.settle_angle(plan, &mut moves)?);
        }
        let centers = plans.iter().map(|plan| plan.atom).collect::<BTreeSet<_>>();
        let mut placed = BTreeMap::new();
        let mut hydrogens = Vec::with_capacity(plans.len());
        for plan in &plans {
            let hydrogen_plan = self.hydrogen_plan(plan, &moves, &centers, &placed)?;
            moves.extend(hydrogen_plan.kept.iter().copied());
            placed.insert(plan.atom, hydrogen_plan.points.clone());
            hydrogens.push(hydrogen_plan);
        }

        let total = hydrogens.iter().map(|plan| plan.new.len()).sum();
        self.slots
            .try_reserve_atoms(total)
            .map_err(super::super::StructureEditError::from)?;
        let mut names = NameAllocator::default();
        let mut report = Vec::with_capacity(plans.len());
        for ((plan, hydrogen_plan), (turned, angle_in_ring)) in
            plans.iter().zip(&hydrogens).zip(angles)
        {
            let atom = plan.atom;
            let declared = self.topology.materializable_hydrogens(atom)?.0 > 0
                || self.topology.owes_hydrogens(atom);
            let mut added = Vec::new();
            if declared || !hydrogen_plan.removed.is_empty() || !hydrogen_plan.new.is_empty() {
                added = self.structural(|topology| {
                    topology.settle_hydrogens(atom, &hydrogen_plan.removed, hydrogen_plan.new.len())
                })?;
                for (hydrogen, point) in added.iter().zip(&hydrogen_plan.new) {
                    self.place_row(*hydrogen, *point);
                    self.join_residue_of(atom, *hydrogen, &mut names)?;
                }
            }
            let moved = hydrogen_plan
                .kept
                .iter()
                .filter(|(hydrogen, point)| {
                    (self.point(*hydrogen) - *point).norm() > MOVE_TOLERANCE
                })
                .map(|(hydrogen, _)| *hydrogen)
                .collect();
            report.push(SettledAtom {
                atom,
                added,
                removed: hydrogen_plan.removed.clone(),
                moved,
                turned,
                angle_in_ring,
            });
        }
        self.apply_moves(moves);
        Ok(SettledAtomsReport { atoms: report })
    }

    /// Gives each distinct bond in `bonds` the ideal length for its elements
    /// and order
    /// ([`ideal_bond_length`](crate::structure::placement::ideal_bond_length))
    /// by moving the smaller side of the bond, by atoms, along it, so a
    /// single bond made double shortens from 1.54 to 1.34 Å without moving
    /// the larger part of the molecule. A bond in a ring has no side that
    /// moves alone; it keeps its length, and its entry reports it. A length
    /// within 0.001 Å of the ideal one stays. Bonds are settled in order,
    /// each from the positions the earlier ones left. Angles stay;
    /// [`Self::settle_atoms`] sets them.
    ///
    /// # Errors
    ///
    /// Nothing changes on error: [`PlacementError::Edit`] for an invalid
    /// handle, [`PlacementError::UnknownBondLength`], and
    /// [`PlacementError::DegenerateGeometry`] when the bond's atoms coincide.
    pub fn settle_bonds(
        &mut self,
        bonds: impl IntoIterator<Item = EditBondId>,
    ) -> Result<SettledBondsReport, PlacementError> {
        let mut seen = BTreeSet::new();
        let mut moves = Moves::new();
        let mut report = Vec::new();
        for bond in bonds {
            if !seen.insert(bond) {
                continue;
            }
            let value = self.topology.require_bond(bond)?;
            let (a, b) = (value.a(), value.b());
            let ideal = bond_length(self.element(a), self.element(b), value.order)?;
            let along = at(self, &moves, b) - at(self, &moves, a);
            let direction = unit(along).ok_or(PlacementError::DegenerateGeometry { atom: a })?;
            let shift = ideal - along.norm();
            let mut settled = SettledBond {
                bond,
                moved: Vec::new(),
                in_ring: false,
            };
            if shift.abs() > LENGTH_TOLERANCE {
                match self.smaller_side([(a, b), (b, a)])? {
                    None => settled.in_ring = true,
                    Some((side, atoms)) => {
                        let shift = if side == 0 { direction } else { -direction } * shift;
                        for atom in &atoms {
                            let point = at(self, &moves, *atom) + shift;
                            moves.insert(*atom, point);
                        }
                        settled.moved = atoms;
                    }
                }
            }
            report.push(settled);
        }
        self.apply_moves(moves);
        Ok(SettledBondsReport { bonds: report })
    }

    /// Whether `atom` is a hydrogen of its one neighbor.
    fn is_hydrogen_of_neighbor(&self, atom: EditAtomId) -> Result<bool, PlacementError> {
        let neighbors = self.topology.neighbors(atom)?.collect::<Vec<_>>();
        match neighbors.as_slice() {
            [neighbor] => self.is_replaceable_hydrogen(*neighbor, atom),
            _ => Ok(false),
        }
    }

    /// Perceives how `atom` settles: its hydrogens and other neighbors, its
    /// builder hydrogen count, and its shape with that count.
    fn atom_plan(&self, atom: EditAtomId) -> Result<AtomPlan, PlacementError> {
        let value = self.topology.require_atom(atom)?.clone();
        let mut others = Vec::new();
        let mut hydrogens = Vec::new();
        let mut valence = 0;
        for bond in self.topology.incident_bonds(atom)?.collect::<Vec<_>>() {
            let bond = self.topology.require_bond(bond)?;
            let other = if bond.a() == atom { bond.b() } else { bond.a() };
            if self.is_replaceable_hydrogen(atom, other)? {
                hydrogens.push(other);
            } else {
                others.push(other);
                valence += usize::from(bond.order.valence());
            }
        }
        let radical = value
            .radical
            .map_or(0, |radical| usize::from(radical.electron_count()));
        let target = if builder_charges(value.element).is_some() {
            builder_hydrogen_count(value.element, value.formal_charge, valence + radical)
                .map_err(|error| PlacementError::BuilderValence { atom, error })?
        } else {
            usize::from(
                rdkit_hydrogens_at_valence(&value, valence, !others.is_empty())
                    .ok_or(PlacementError::InvalidValence { atom })?,
            )
        };
        let (domains, lone_pairs) = self.domains_with(
            atom,
            DomainChange {
                without: &hydrogens,
                substituents: target,
                valence: target,
                coordinate: 0,
                multiple: false,
            },
        )?;
        let shape = match domains {
            0 => None,
            _ => Some(
                AtomShape::from_domains(domains)
                    .ok_or(PlacementError::UnsupportedGeometry { atom, domains })?,
            ),
        };
        let length = match target {
            0 => 0.0,
            _ => bond_length(value.element, Element::H, BondOrder::Single)?,
        };
        self.topology
            .check_settled_stereo(atom, target.saturating_sub(hydrogens.len()))
            .map_err(|error| draft_error(error, atom))?;
        Ok(AtomPlan {
            atom,
            others,
            hydrogens,
            target,
            shape,
            lone_pairs,
            length,
        })
    }

    /// Turns the smaller side of an atom with two other neighbors to the
    /// angle of its shape, in `moves`; returns the side's atoms and whether
    /// the angle is off but in a ring.
    fn settle_angle(
        &self,
        plan: &AtomPlan,
        moves: &mut Moves,
    ) -> Result<(Vec<EditAtomId>, bool), PlacementError> {
        let &[first, second] = plan.others.as_slice() else {
            return Ok((Vec::new(), false));
        };
        let ideal = match plan.shape {
            Some(AtomShape::Linear) => 180.0_f64,
            Some(AtomShape::TrigonalPlanar) => 120.0,
            Some(AtomShape::Tetrahedral) => TETRAHEDRAL,
            _ => return Ok((Vec::new(), false)),
        };
        let center = at(self, moves, plan.atom);
        let (Some(towards_first), Some(towards_second)) = (
            unit(at(self, moves, first) - center),
            unit(at(self, moves, second) - center),
        ) else {
            return Ok((Vec::new(), false));
        };
        let angle = towards_first
            .cross(towards_second)
            .norm()
            .atan2(towards_first.dot(towards_second))
            .to_degrees();
        if (angle - ideal).abs() <= ANGLE_TOLERANCE {
            return Ok((Vec::new(), false));
        }
        let Some((side, atoms)) = self.smaller_side([(plan.atom, first), (plan.atom, second)])?
        else {
            return Ok((Vec::new(), true));
        };
        let (turning, fixed) = if side == 0 {
            (towards_first, towards_second)
        } else {
            (towards_second, towards_first)
        };
        let plane =
            unit(turning - fixed * turning.dot(fixed)).unwrap_or_else(|| perpendicular(fixed));
        let (sin, cos) = ideal.to_radians().sin_cos();
        let rotation = Rotation::between(turning, fixed * cos + plane * sin);
        for atom in &atoms {
            let point = center + rotation.apply(at(self, moves, *atom) - center);
            moves.insert(*atom, point);
        }
        Ok((atoms, false))
    }

    /// Where the hydrogens of the atom of `plan` go: the free sites of its
    /// shape from its other neighbors and lone pairs, kept hydrogens on the
    /// sites nearest them.
    fn hydrogen_plan(
        &self,
        plan: &AtomPlan,
        moves: &Moves,
        centers: &BTreeSet<EditAtomId>,
        placed: &BTreeMap<EditAtomId, Vec<Point3>>,
    ) -> Result<HydrogenPlan, PlacementError> {
        let atom = plan.atom;
        let Some(shape) = plan.shape.filter(|_| plan.target > 0) else {
            return Ok(HydrogenPlan {
                removed: plan.hydrogens.clone(),
                ..HydrogenPlan::default()
            });
        };
        let center = at(self, moves, atom);
        let neighbors = plan
            .others
            .iter()
            .map(|other| at(self, moves, *other))
            .collect::<Vec<_>>();
        let reference = match plan.others.as_slice() {
            [neighbor] => self.settled_reference(*neighbor, atom, moves, centers, placed)?,
            _ => None,
        };
        let mut sites = free_sites(SiteQuery {
            center,
            neighbors: &neighbors,
            reference,
            shape,
            lone_pairs: plan.lone_pairs,
        })
        .map_err(|_| PlacementError::DegenerateGeometry { atom })?;
        if sites.len() < plan.target {
            return Err(PlacementError::NoFreePosition { atom });
        }
        sites.truncate(plan.target);

        // Each kept hydrogen takes the free site nearest it, best matches
        // first; the hydrogens matched to none are deleted.
        let mut pairs = Vec::with_capacity(plan.hydrogens.len() * sites.len());
        for (index, hydrogen) in plan.hydrogens.iter().enumerate() {
            let toward = unit(at(self, moves, *hydrogen) - center);
            for (site_index, site) in sites.iter().enumerate() {
                let alignment = toward.map_or(-2.0, |toward| toward.dot(*site));
                pairs.push((alignment, index, site_index));
            }
        }
        pairs.sort_by(|a, b| b.0.total_cmp(&a.0).then((a.1, a.2).cmp(&(b.1, b.2))));
        let mut hydrogen_site = vec![None; plan.hydrogens.len()];
        let mut site_taken = vec![false; sites.len()];
        for (_, index, site_index) in pairs {
            if hydrogen_site[index].is_none() && !site_taken[site_index] {
                hydrogen_site[index] = Some(site_index);
                site_taken[site_index] = true;
            }
        }
        let points = sites
            .iter()
            .map(|site| center + *site * plan.length)
            .collect::<Vec<_>>();
        let mut result = HydrogenPlan {
            points: points.clone(),
            ..HydrogenPlan::default()
        };
        for (hydrogen, site) in plan.hydrogens.iter().zip(hydrogen_site) {
            match site {
                Some(site) => result.kept.push((*hydrogen, points[site])),
                None => result.removed.push(*hydrogen),
            }
        }
        result.new = site_taken
            .iter()
            .zip(&points)
            .filter(|(taken, _)| !**taken)
            .map(|(_, point)| *point)
            .collect();
        Ok(result)
    }

    /// The substituent of `neighbor` other than `atom` that the hydrogens of
    /// `atom` stagger against: a heavy atom with the most heavy neighbors,
    /// else a hydrogen where it settles or already is; `None` when a settled
    /// neighbor's hydrogens are not placed yet.
    fn settled_reference(
        &self,
        neighbor: EditAtomId,
        atom: EditAtomId,
        moves: &Moves,
        centers: &BTreeSet<EditAtomId>,
        placed: &BTreeMap<EditAtomId, Vec<Point3>>,
    ) -> Result<Option<Point3>, PlacementError> {
        let mut best: Option<(usize, EditAtomId)> = None;
        let mut hydrogen = None;
        for candidate in self.topology.neighbors(neighbor)?.collect::<Vec<_>>() {
            if candidate == atom {
                continue;
            }
            if self.element(candidate) == Element::H {
                hydrogen = hydrogen.or(Some(candidate));
                continue;
            }
            let rank = self
                .topology
                .neighbors(candidate)?
                .filter(|other| self.element(*other) != Element::H)
                .count();
            if best.is_none_or(|(score, _)| rank > score) {
                best = Some((rank, candidate));
            }
        }
        if let Some((_, heavy)) = best {
            return Ok(Some(at(self, moves, heavy)));
        }
        if centers.contains(&neighbor) {
            return Ok(placed
                .get(&neighbor)
                .and_then(|points| points.first().copied()));
        }
        Ok(hydrogen.map(|hydrogen| at(self, moves, hydrogen)))
    }

    /// Of the sides beyond each `(near, far)` pair, the atoms reached from
    /// `far` without crossing its bond to `near`, the one with fewer atoms
    /// that does not reach `near`, with its index; the first on a tie.
    /// `None` when both reach `near`, both bonds lying in a ring. The search
    /// stops when the smaller side is complete.
    fn smaller_side(
        &self,
        pairs: [(EditAtomId, EditAtomId); 2],
    ) -> Result<Option<(usize, Vec<EditAtomId>)>, PlacementError> {
        struct Walk {
            near: EditAtomId,
            far: EditAtomId,
            seen: BTreeSet<EditAtomId>,
            atoms: Vec<EditAtomId>,
            queue: VecDeque<EditAtomId>,
            in_ring: bool,
        }
        let mut walks = pairs.map(|(near, far)| Walk {
            near,
            far,
            seen: BTreeSet::from([far]),
            atoms: vec![far],
            queue: VecDeque::from([far]),
            in_ring: false,
        });
        loop {
            for (index, walk) in walks.iter_mut().enumerate() {
                if walk.in_ring {
                    continue;
                }
                let Some(current) = walk.queue.pop_front() else {
                    return Ok(Some((index, std::mem::take(&mut walk.atoms))));
                };
                for next in self.topology.neighbors(current)? {
                    if current == walk.far && next == walk.near {
                        continue;
                    }
                    if next == walk.near {
                        walk.in_ring = true;
                        break;
                    }
                    if walk.seen.insert(next) {
                        walk.atoms.push(next);
                        walk.queue.push_back(next);
                    }
                }
            }
            if walks.iter().all(|walk| walk.in_ring) {
                return Ok(None);
            }
        }
    }

    /// Writes the planned positions of the live atoms that move.
    fn apply_moves(&mut self, moves: Moves) {
        for (atom, point) in moves {
            if self.topology.atom(atom).is_some()
                && (self.point(atom) - point).norm() > MOVE_TOLERANCE
            {
                self.place_row(atom, point);
            }
        }
    }
}

/// The planned position of `atom`, or its draft position.
fn at(editor: &StructureEditor, moves: &Moves, atom: EditAtomId) -> Point3 {
    moves
        .get(&atom)
        .copied()
        .unwrap_or_else(|| editor.point(atom))
}
