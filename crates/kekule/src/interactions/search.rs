//! The geometric search of one structure, or two.

use std::sync::Arc;

use crate::geometry::{Point3, PointGrid, Vector3};
use crate::io::is_metal;
use crate::molecule::{AtomId, Element};
use crate::structure::{AsModelView, ModelView};
use crate::topology::{AtomSelection, InstanceAtomId, Topology, TopologyAtomIndex};

use super::criteria::Canonical;
use super::features::{centroid, ring_geometry, InteractionFeatures};
use super::{
    Interaction, InteractionAngle, InteractionCriteria, InteractionEnd, InteractionError,
    InteractionKind, InteractionRole, Interactions,
};

/// Candidate pairs at most this many are measured pair by pair; more search
/// a grid of one list.
const DIRECT_PAIRS: usize = 256;

/// One side of a search across two structures: features of a topology, a
/// structure sharing its layout, and the atoms taking part.
#[derive(Debug, Clone, Copy)]
pub struct InteractionSide<'a> {
    features: &'a InteractionFeatures,
    structure: ModelView<'a>,
    atoms: &'a AtomSelection,
}

impl<'a> InteractionSide<'a> {
    /// A side of `structure`, a model, ensemble sample, or trajectory frame,
    /// whose `atoms` take part, with `features` of its topology.
    ///
    /// # Errors
    ///
    /// [`InteractionError::TopologyMismatch`] when `structure` or `atoms`
    /// does not share the features' topology layout.
    pub fn new(
        features: &'a InteractionFeatures,
        structure: &'a (impl AsModelView + ?Sized),
        atoms: &'a AtomSelection,
    ) -> Result<Self, InteractionError> {
        let structure = structure.as_model_view();
        if !structure.topology().shares_layout(features.topology())
            || !atoms.topology().shares_layout(features.topology())
        {
            return Err(InteractionError::TopologyMismatch);
        }
        Ok(Self {
            features,
            structure,
            atoms,
        })
    }

    /// The side's features.
    pub const fn features(&self) -> &'a InteractionFeatures {
        self.features
    }

    /// The side's structure.
    pub const fn structure(&self) -> ModelView<'a> {
        self.structure
    }

    /// The atoms taking part.
    pub const fn atoms(&self) -> &'a AtomSelection {
        self.atoms
    }
}

/// Finds the interactions between atoms of `first` and atoms of `second`,
/// two sets of atoms of one `structure` with `features` of its topology,
/// under `criteria`.
///
/// The sets may overlap: interactions among atoms of both sets are found
/// once. An end takes part through any of its atoms: a ring with one
/// selected atom, a charged group, or a donor whose hydrogen is selected. An
/// interaction's first end belongs to `first`; when either end could come
/// first, the donor, anion, cation (before a ring), halogen, metal, or end
/// with the lower atom index does. Two ends in one molecule fewer than
/// [`InteractionCriteria::min_bond_separation`] bonds apart, or sharing an
/// atom, do not interact, except a metal and its ligand. An atom pair that
/// forms hydrogen bonds in both directions is reported once: by its explicit
/// hydrogen when one qualifies, otherwise with the first end donating.
///
/// # Errors
///
/// [`InteractionError::TopologyMismatch`] when `structure`, `first`, or
/// `second` does not share the features' topology layout, and
/// [`InteractionError::InvalidCriterion`] for an invalid criterion.
pub fn find_interactions_within(
    features: &InteractionFeatures,
    structure: &(impl AsModelView + ?Sized),
    first: &AtomSelection,
    second: &AtomSelection,
    criteria: &InteractionCriteria,
) -> Result<Interactions, InteractionError> {
    let canonical = criteria.canonical()?;
    let first = InteractionSide::new(features, structure, first)?;
    let second = InteractionSide::new(features, structure, second)?;
    Ok(Search::new(&first, &second, true, canonical).run(criteria))
}

/// Finds the interactions between atoms of two structures, each side with
/// features of its own topology, under `criteria`.
///
/// Each side's features come from its own topology, so the structures may
/// be unrelated, such as two objects a host shows. An interaction's first
/// end belongs to `first`. Bonds never join the sides, so no pair is
/// excluded by bond separation; searching one structure split into two
/// gives the interactions [`find_interactions_within`] gives between the
/// parts, when no bond joins them.
///
/// # Errors
///
/// [`InteractionError::InvalidCriterion`] for an invalid criterion.
pub fn find_interactions_between(
    first: &InteractionSide<'_>,
    second: &InteractionSide<'_>,
    criteria: &InteractionCriteria,
) -> Result<Interactions, InteractionError> {
    let canonical = criteria.canonical()?;
    Ok(Search::new(first, second, false, canonical).run(criteria))
}

/// One feature near a point, by its index in its features' list.
#[derive(Debug, Clone, Copy)]
struct Candidate {
    point: Point3,
    feature: u32,
    /// A ring's normal; zero for other features.
    normal: Vector3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Role {
    Donor,
    Acceptor,
    Anion,
    Cation,
    Ring,
    Halogen,
    Metal,
    Ligand,
    Hydrophobic,
}

impl Role {
    const fn public(self) -> InteractionRole {
        match self {
            Self::Donor => InteractionRole::Donor,
            Self::Acceptor => InteractionRole::Acceptor,
            Self::Anion => InteractionRole::Anion,
            Self::Cation => InteractionRole::Cation,
            Self::Ring => InteractionRole::Ring,
            Self::Halogen => InteractionRole::HalogenDonor,
            Self::Metal => InteractionRole::Metal,
            Self::Ligand => InteractionRole::MetalLigand,
            Self::Hydrophobic => InteractionRole::Hydrophobic,
        }
    }
}

/// One side as a search reads it.
struct Side<'a> {
    features: &'a InteractionFeatures,
    topology: &'a Topology,
    positions: &'a [Point3],
    mask: Vec<u64>,
    /// The number of selected atoms.
    count: usize,
    /// Every candidate of each role, listed once when this is the side
    /// with fewer atoms.
    lists: [Option<Vec<Candidate>>; ROLES],
    /// The atoms within reach of the other side's candidates, in dense
    /// order, when this is the side with more atoms.
    nearby: Option<Vec<TopologyAtomIndex>>,
}

const ROLES: usize = 9;

impl<'a> Side<'a> {
    fn new(side: &InteractionSide<'a>) -> Self {
        Self {
            features: side.features,
            topology: side.features.topology(),
            positions: side.structure.positions().values().value(),
            mask: side.atoms.to_mask(),
            count: side.atoms.len(),
            lists: Default::default(),
            nearby: None,
        }
    }

    fn selected(&self, atom: TopologyAtomIndex) -> bool {
        let index = atom.index();
        self.mask[index / 64] & (1 << (index % 64)) != 0
    }

    fn position(&self, atom: TopologyAtomIndex) -> Point3 {
        self.positions[atom.index()]
    }

    /// The atoms of a feature that decide membership and exclusion.
    fn atoms(&self, role: Role, feature: u32) -> &'a [TopologyAtomIndex] {
        let features = self.features;
        let feature = feature as usize;
        match role {
            Role::Donor => std::slice::from_ref(&features.donors[feature].atom),
            Role::Acceptor => std::slice::from_ref(&features.acceptors[feature]),
            Role::Anion | Role::Cation => features.span(&features.groups[feature].atoms),
            Role::Ring => features.span(&features.rings[feature]),
            Role::Halogen => std::slice::from_ref(&features.halogens[feature].halogen),
            Role::Metal => std::slice::from_ref(&features.metals[feature]),
            Role::Ligand => std::slice::from_ref(&features.ligands[feature]),
            Role::Hydrophobic => std::slice::from_ref(&features.hydrophobic[feature]),
        }
    }

    /// Whether the side's atoms include the feature.
    fn includes(&self, role: Role, feature: u32) -> bool {
        self.atoms(role, feature)
            .iter()
            .any(|&atom| self.selected(atom))
            || (role == Role::Donor
                && self
                    .features
                    .span(&self.features.donors[feature as usize].hydrogens)
                    .iter()
                    .any(|&atom| self.selected(atom)))
    }

    /// The side's included features of `role` with their points, only
    /// those inside `region` when one is given.
    fn candidates(&self, role: Role, region: Option<(Point3, Point3)>) -> Vec<Candidate> {
        let features = self.features;
        let inside = |point: &Point3| {
            region.is_none_or(|(low, high)| {
                point.x >= low.x
                    && point.y >= low.y
                    && point.z >= low.z
                    && point.x <= high.x
                    && point.y <= high.y
                    && point.z <= high.z
            })
        };
        let plain = |feature: usize, point: Point3| Candidate {
            point,
            feature: feature as u32,
            normal: Vector3::zero(),
        };
        // The features of a single-atom role, by their atoms sorted in dense
        // order: those of the nearby atoms when they are known.
        let by_atom = |count: usize, atom: &dyn Fn(usize) -> TopologyAtomIndex| -> Vec<Candidate> {
            let keep = |index: usize| {
                let point = self.position(atom(index));
                (inside(&point) && self.includes(role, index as u32)).then(|| plain(index, point))
            };
            match &self.nearby {
                Some(nearby) => {
                    let mut found = Vec::new();
                    let mut from = 0;
                    for &near in nearby {
                        // Both lists ascend, so each search starts where the
                        // last one ended.
                        from = partition(from, count, |index| atom(index) < near);
                        if from < count && atom(from) == near {
                            found.extend(keep(from));
                        }
                    }
                    found
                }
                None => (0..count).filter_map(keep).collect(),
            }
        };
        let single = |list: &[TopologyAtomIndex]| by_atom(list.len(), &|index| list[index]);
        match role {
            Role::Donor => by_atom(features.donors.len(), &|index| features.donors[index].atom),
            Role::Acceptor => single(&features.acceptors),
            Role::Anion | Role::Cation => {
                let sign = if role == Role::Anion {
                    super::ChargeSign::Negative
                } else {
                    super::ChargeSign::Positive
                };
                features
                    .groups
                    .iter()
                    .enumerate()
                    .filter(|(index, group)| {
                        group.sign == sign && self.includes(role, *index as u32)
                    })
                    .map(|(index, group)| {
                        let atoms = features.span(&group.atoms);
                        let charged = &atoms[..group.charged as usize];
                        plain(
                            index,
                            centroid(charged.iter().map(|&atom| self.position(atom))),
                        )
                    })
                    .filter(|candidate| inside(&candidate.point))
                    .collect()
            }
            // A ring's center is its centroid, tested before its normal is
            // computed.
            Role::Ring => (0..features.rings.len() as u32)
                .filter(|&index| self.includes(role, index))
                .filter_map(|index| {
                    let atoms = features.span(&features.rings[index as usize]);
                    if !inside(&centroid(atoms.iter().map(|&atom| self.position(atom)))) {
                        return None;
                    }
                    ring_geometry(atoms, self.positions).map(|(point, normal)| Candidate {
                        point,
                        feature: index,
                        normal,
                    })
                })
                .collect(),
            Role::Halogen => by_atom(features.halogens.len(), &|index| {
                features.halogens[index].halogen
            }),
            Role::Metal => single(&features.metals),
            Role::Ligand => single(&features.ligands),
            Role::Hydrophobic => single(&features.hydrophobic),
        }
    }

    /// The heavy neighbors of `atom` other than metals, whose coordination
    /// bonds do not orient a donor or acceptor.
    fn heavy_neighbors(&self, atom: TopologyAtomIndex) -> Vec<TopologyAtomIndex> {
        let id = self
            .topology
            .atom_id(atom)
            .expect("feature atoms are topology atoms");
        let molecule = self.topology.definition_molecule(id.instance());
        molecule
            .neighbors(id.atom())
            .expect("feature atoms are live")
            .filter(|&neighbor| {
                let element = molecule
                    .atom(neighbor)
                    .expect("neighbors are live")
                    .element();
                element != Element::H && !is_metal(element)
            })
            .map(|neighbor| {
                self.topology
                    .atom_index(InstanceAtomId::new(id.instance(), neighbor))
                    .expect("neighbors are topology atoms")
            })
            .collect()
    }
}

/// The angle at `vertex` between `a` and `b`, or `None` when either arm has
/// no length.
fn angle(a: Point3, vertex: Point3, b: Point3) -> Option<f64> {
    let (u, v) = (a - vertex, b - vertex);
    let lengths = u.norm() * v.norm();
    (lengths > 0.0 && lengths.is_finite()).then(|| (u.dot(v) / lengths).clamp(-1.0, 1.0).acos())
}

/// The distance of `vector` from the axis along unit `normal`.
fn lateral(vector: Vector3, normal: Vector3) -> f64 {
    (vector - normal * vector.dot(normal)).norm()
}

/// The first index in `from..to` for which `before` is false, where it is
/// true for a prefix.
fn partition(mut from: usize, mut to: usize, before: impl Fn(usize) -> bool) -> usize {
    while from < to {
        let middle = from + (to - from) / 2;
        if before(middle) {
            from = middle + 1;
        } else {
            to = middle;
        }
    }
    from
}

/// The bounding box of nonempty candidates.
fn bounds(candidates: &[Candidate]) -> (Point3, Point3) {
    let first = candidates[0].point;
    candidates
        .iter()
        .fold((first, first), |(low, high), candidate| {
            let point = candidate.point;
            (
                Point3::new(low.x.min(point.x), low.y.min(point.y), low.z.min(point.z)),
                Point3::new(
                    high.x.max(point.x),
                    high.y.max(point.y),
                    high.z.max(point.z),
                ),
            )
        })
}

/// Passes every pair of a `left` and a `right` candidate within `cutoff` to
/// `each`.
fn pairs(
    left: &[Candidate],
    right: &[Candidate],
    cutoff: f64,
    grid: &mut PointGrid<u32>,
    mut each: impl FnMut(&Candidate, &Candidate),
) {
    if left.is_empty() || right.is_empty() {
        return;
    }
    let within = |a: &Point3, b: &Point3| (*a - *b).norm() <= cutoff;
    // The larger list keeps only the candidates inside the smaller one's
    // bounding box, widened by the cutoff.
    let left_is_small = left.len() <= right.len();
    let (small, large) = if left_is_small {
        (left, right)
    } else {
        (right, left)
    };
    let (low, high) = bounds(small);
    let inside = |point: &Point3| {
        point.x >= low.x - cutoff
            && point.y >= low.y - cutoff
            && point.z >= low.z - cutoff
            && point.x <= high.x + cutoff
            && point.y <= high.y + cutoff
            && point.z <= high.z + cutoff
    };
    let near = large
        .iter()
        .filter(|candidate| inside(&candidate.point))
        .collect::<Vec<_>>();
    if near.is_empty() {
        return;
    }
    let mut ordered = |a: &Candidate, b: &Candidate| {
        if left_is_small {
            each(a, b);
        } else {
            each(b, a);
        }
    };
    if small.len().saturating_mul(near.len()) <= DIRECT_PAIRS {
        for a in small {
            for b in &near {
                if within(&a.point, &b.point) {
                    ordered(a, b);
                }
            }
        }
        return;
    }
    let (low, high) = bounds(small);
    grid.build(
        small
            .iter()
            .enumerate()
            .map(|(index, candidate)| (candidate.point, index as u32)),
        (low, high),
        cutoff,
    );
    // Each near candidate meets the small candidates around it, in the
    // small list's order.
    let mut found = Vec::new();
    for b in near {
        found.clear();
        grid.for_each_near(&b.point, within, |_, index| found.push(index));
        found.sort_unstable();
        for &index in &found {
            ordered(&small[index as usize], b);
        }
    }
}

/// One search: two sides and the criteria in canonical units.
struct Search<'a> {
    first: Side<'a>,
    second: Side<'a>,
    first_topology: Arc<Topology>,
    second_topology: Arc<Topology>,
    within: bool,
    criteria: Canonical,
    grid: PointGrid<u32>,
    found: Vec<Interaction>,
}

impl<'a> Search<'a> {
    fn new(
        first: &InteractionSide<'a>,
        second: &InteractionSide<'a>,
        within: bool,
        criteria: Canonical,
    ) -> Self {
        Self {
            first: Side::new(first),
            second: Side::new(second),
            first_topology: first.features.shared_topology(),
            second_topology: second.features.shared_topology(),
            within,
            criteria,
            grid: PointGrid::default(),
            found: Vec::new(),
        }
    }

    fn run(mut self, criteria: &InteractionCriteria) -> Interactions {
        let c = self.criteria;
        let searched = [
            (
                c.hydrogen_bond.enabled,
                [Role::Donor, Role::Acceptor],
                c.hydrogen_bond.distance,
            ),
            (
                c.salt_bridge.enabled,
                [Role::Anion, Role::Cation],
                c.salt_bridge.distance,
            ),
            (
                c.pi_stacking.enabled,
                [Role::Ring, Role::Ring],
                c.pi_stacking.distance,
            ),
            (
                c.cation_pi.enabled,
                [Role::Cation, Role::Ring],
                c.cation_pi.distance,
            ),
            (
                c.halogen_bond.enabled,
                [Role::Halogen, Role::Acceptor],
                c.halogen_bond.distance,
            ),
            (
                c.metal_coordination.enabled,
                [Role::Metal, Role::Ligand],
                c.metal_coordination.distance,
            ),
            (
                c.hydrophobic.enabled,
                [Role::Hydrophobic, Role::Hydrophobic],
                c.hydrophobic.distance,
            ),
        ];
        // The side with more atoms reads only its atoms within reach of the
        // other side's candidates: those inside the union of their bounding
        // boxes, each widened by its search's cutoff.
        let small_is_first = self.first.count <= self.second.count;
        let (small, large) = if small_is_first {
            (&mut self.first, &mut self.second)
        } else {
            (&mut self.second, &mut self.first)
        };
        let mut regions = Vec::new();
        for (enabled, roles, cutoff) in searched {
            for role in roles.into_iter().filter(|_| enabled) {
                if small.lists[role as usize].is_none() {
                    small.lists[role as usize] = Some(small.candidates(role, None));
                }
                let listed = small.lists[role as usize].as_deref().unwrap_or_default();
                if !listed.is_empty() {
                    let (low, high) = bounds(listed);
                    let widen = Vector3::new(cutoff, cutoff, cutoff);
                    regions.push((low - widen, high + widen));
                }
            }
        }
        large.nearby = Some(
            (0..large.positions.len() as u32)
                .map(TopologyAtomIndex::new)
                .filter(|atom| {
                    let point = large.positions[atom.index()];
                    regions.iter().any(|(low, high)| {
                        point.x >= low.x
                            && point.y >= low.y
                            && point.z >= low.z
                            && point.x <= high.x
                            && point.y <= high.y
                            && point.z <= high.z
                    })
                })
                .collect(),
        );
        if c.hydrogen_bond.enabled {
            self.search(Role::Donor, Role::Acceptor, c.hydrogen_bond.distance);
        }
        if c.salt_bridge.enabled {
            self.search(Role::Anion, Role::Cation, c.salt_bridge.distance);
        }
        if c.pi_stacking.enabled {
            self.search(Role::Ring, Role::Ring, c.pi_stacking.distance);
        }
        if c.cation_pi.enabled {
            self.search(Role::Cation, Role::Ring, c.cation_pi.distance);
        }
        if c.halogen_bond.enabled {
            self.search(Role::Halogen, Role::Acceptor, c.halogen_bond.distance);
        }
        if c.metal_coordination.enabled {
            self.search(Role::Metal, Role::Ligand, c.metal_coordination.distance);
        }
        if c.hydrophobic.enabled {
            self.search(Role::Hydrophobic, Role::Hydrophobic, c.hydrophobic.distance);
        }
        let mut found = std::mem::take(&mut self.found);
        found.sort_by(|a, b| {
            (
                a.kind,
                &a.first.atoms,
                &a.second.atoms,
                a.first.hydrogen,
                a.second.hydrogen,
            )
                .cmp(&(
                    b.kind,
                    &b.first.atoms,
                    &b.second.atoms,
                    b.first.hydrogen,
                    b.second.hydrogen,
                ))
        });
        // An atom pair bonded in both directions keeps one hydrogen bond: by
        // an explicit hydrogen, otherwise with the first end donating, then
        // the lower donor.
        let pair = |interaction: &Interaction| {
            let (first, second) = (interaction.first.atoms[0], interaction.second.atoms[0]);
            if self.within && second < first {
                (second, first)
            } else {
                (first, second)
            }
        };
        let rank = |interaction: &Interaction| {
            let explicit = interaction.hydrogen_distance.is_some();
            let first_donates = interaction.first.role == InteractionRole::Donor;
            let donor = if first_donates {
                interaction.first.atoms[0]
            } else {
                interaction.second.atoms[0]
            };
            (!explicit, !first_donates, donor)
        };
        let mut best = std::collections::HashMap::new();
        for (index, interaction) in found.iter().enumerate() {
            if interaction.kind == InteractionKind::HydrogenBond {
                best.entry(pair(interaction))
                    .and_modify(|kept: &mut usize| {
                        if rank(interaction) < rank(&found[*kept]) {
                            *kept = index;
                        }
                    })
                    .or_insert(index);
            }
        }
        let kept = found
            .into_iter()
            .enumerate()
            .filter(|(index, interaction)| {
                interaction.kind != InteractionKind::HydrogenBond
                    || best.get(&pair(interaction)) == Some(index)
            })
            .map(|(_, interaction)| interaction)
            .collect::<Vec<_>>();
        Interactions {
            first: self.first_topology,
            second: self.second_topology,
            within: self.within,
            hydrogen_bond: criteria.hydrogen_bond(),
            interactions: kept,
        }
    }

    /// Searches pairs of a `p` feature and a `q` feature within `cutoff`, in
    /// both orientations unless the roles are the same.
    fn search(&mut self, p: Role, q: Role, cutoff: f64) {
        let symmetric = p == q;
        let orientations: &[bool] = if symmetric { &[false] } else { &[false, true] };
        for &flipped in orientations {
            // The `p` features come from the first side unless flipped.
            let (p_list, q_list) = self.lists(p, q, flipped, cutoff);
            let mut pairs_found = Vec::new();
            pairs(&p_list, &q_list, cutoff, &mut self.grid, |a, b| {
                pairs_found.push((*a, *b));
            });
            for (a, b) in pairs_found {
                if self.within && self.skip(p, q, flipped, symmetric, &a, &b) {
                    continue;
                }
                let (p_side, q_side) = if flipped {
                    (&self.second, &self.first)
                } else {
                    (&self.first, &self.second)
                };
                if let Some(interaction) = evaluate(&self.criteria, p, q, p_side, q_side, &a, &b) {
                    let interaction = if flipped {
                        Interaction {
                            first: interaction.second,
                            second: interaction.first,
                            ..interaction
                        }
                    } else {
                        interaction
                    };
                    self.found.push(interaction);
                }
            }
        }
    }

    /// The candidates of `p` on its side and `q` on the other. The side
    /// with fewer selected atoms lists all of its candidates, once per role;
    /// the other lists only those within `cutoff` of their bounding box.
    fn lists(
        &mut self,
        p: Role,
        q: Role,
        flipped: bool,
        cutoff: f64,
    ) -> (Vec<Candidate>, Vec<Candidate>) {
        let small_is_first = self.first.count <= self.second.count;
        // Whether the `p` features come from the side with fewer atoms: the
        // first side unless flipped.
        let p_small = flipped != small_is_first;
        let (small_role, large_role) = if p_small { (p, q) } else { (q, p) };
        let (small, large) = if small_is_first {
            (&mut self.first, &self.second)
        } else {
            (&mut self.second, &self.first)
        };
        if small.lists[small_role as usize].is_none() {
            small.lists[small_role as usize] = Some(small.candidates(small_role, None));
        }
        let listed = small.lists[small_role as usize].clone().unwrap_or_default();
        let found = if listed.is_empty() {
            Vec::new()
        } else {
            let (low, high) = bounds(&listed);
            let widen = Vector3::new(cutoff, cutoff, cutoff);
            large.candidates(large_role, Some((low - widen, high + widen)))
        };
        if p_small {
            (listed, found)
        } else {
            (found, listed)
        }
    }

    /// Whether a pair within one structure is excluded or found in the
    /// other orientation.
    fn skip(
        &self,
        p: Role,
        q: Role,
        flipped: bool,
        symmetric: bool,
        a: &Candidate,
        b: &Candidate,
    ) -> bool {
        let side = &self.first;
        let (p_atoms, q_atoms) = (side.atoms(p, a.feature), side.atoms(q, b.feature));
        if symmetric {
            if a.feature == b.feature {
                return true;
            }
            // Both ends in both sets: the lower atom index comes first.
            let reversible =
                self.first.includes(q, b.feature) && self.second.includes(p, a.feature);
            if reversible && q_atoms < p_atoms {
                return true;
            }
        } else if flipped {
            // Found unflipped when the `p` end is also in the first set and
            // the `q` end in the second.
            if self.first.includes(p, a.feature) && self.second.includes(q, b.feature) {
                return true;
            }
        }
        if p_atoms.iter().any(|atom| q_atoms.contains(atom)) {
            return true;
        }
        if p == Role::Metal {
            return false;
        }
        !separated(
            side.topology,
            p_atoms,
            q_atoms,
            self.criteria.min_bond_separation,
        )
    }
}

/// Whether every atom of `a` is at least `min` bonds from every atom of `b`
/// in one molecule, not counting paths through metals; atoms of different
/// molecules are.
fn separated(
    topology: &Topology,
    a: &[TopologyAtomIndex],
    b: &[TopologyAtomIndex],
    min: usize,
) -> bool {
    if min <= 1 {
        return true;
    }
    let id = |atom: &TopologyAtomIndex| {
        topology
            .atom_id(*atom)
            .expect("feature atoms are topology atoms")
    };
    let instance = id(&a[0]).instance();
    let targets = b
        .iter()
        .map(id)
        .filter(|target| target.instance() == instance)
        .map(InstanceAtomId::atom)
        .collect::<Vec<AtomId>>();
    if targets.is_empty() {
        return true;
    }
    let molecule = topology.definition_molecule(instance);
    let metal = |atom: AtomId| is_metal(molecule.atom(atom).expect("live atom").element());
    let mut seen = a.iter().map(|atom| id(atom).atom()).collect::<Vec<_>>();
    let mut frontier = seen.clone();
    // Atoms fewer than `min` bonds away are within `min - 1` steps.
    for _ in 1..min {
        let mut next = Vec::new();
        for &atom in &frontier {
            if metal(atom) {
                continue;
            }
            for neighbor in molecule.neighbors(atom).expect("live atom") {
                if seen.contains(&neighbor) || metal(neighbor) {
                    continue;
                }
                if targets.contains(&neighbor) {
                    return false;
                }
                seen.push(neighbor);
                next.push(neighbor);
            }
        }
        frontier = next;
    }
    true
}

fn end(side: &Side<'_>, role: Role, candidate: &Candidate) -> InteractionEnd {
    InteractionEnd {
        role: role.public(),
        atoms: side.atoms(role, candidate.feature).to_vec(),
        hydrogen: None,
        center: candidate.point,
        normal: (role == Role::Ring).then_some(candidate.normal),
    }
}

/// The interaction of a `p` candidate of `p_side` and a `q` candidate of
/// `q_side`, the `p` end first, when it meets the criteria.
fn evaluate(
    criteria: &Canonical,
    p: Role,
    q: Role,
    p_side: &Side<'_>,
    q_side: &Side<'_>,
    a: &Candidate,
    b: &Candidate,
) -> Option<Interaction> {
    let distance = (a.point - b.point).norm();
    let build = |kind, angles, offset| Interaction {
        kind,
        first: end(p_side, p, a),
        second: end(q_side, q, b),
        distance,
        hydrogen_distance: None,
        offset,
        angles,
    };
    match (p, q) {
        (Role::Donor, Role::Acceptor) => hydrogen_bond(criteria, p_side, q_side, a, b),
        (Role::Anion, Role::Cation) => (distance <= criteria.salt_bridge.distance)
            .then(|| build(InteractionKind::SaltBridge, Vec::new(), None)),
        (Role::Ring, Role::Ring) => {
            let c = criteria.pi_stacking;
            let planes = a.normal.dot(b.normal).abs().clamp(0.0, 1.0).acos();
            let delta = b.point - a.point;
            let offset = lateral(delta, a.normal).min(lateral(delta, b.normal));
            if distance > c.distance || offset > c.offset {
                return None;
            }
            let kind = if planes <= c.parallel_angle {
                InteractionKind::PiStackingParallel
            } else if planes >= c.perpendicular_angle {
                InteractionKind::PiStackingTShaped
            } else {
                return None;
            };
            Some(build(
                kind,
                vec![(InteractionAngle::RingPlanes, planes)],
                Some(offset),
            ))
        }
        (Role::Cation, Role::Ring) => {
            let c = criteria.cation_pi;
            let delta = a.point - b.point;
            let offset = lateral(delta, b.normal);
            if distance > c.distance || offset > c.offset {
                return None;
            }
            let normal = if distance > 0.0 {
                (delta.dot(b.normal).abs() / distance)
                    .clamp(0.0, 1.0)
                    .acos()
            } else {
                0.0
            };
            Some(build(
                InteractionKind::CationPi,
                vec![(InteractionAngle::RingNormalCation, normal)],
                Some(offset),
            ))
        }
        (Role::Halogen, Role::Acceptor) => {
            let c = criteria.halogen_bond;
            if distance > c.distance {
                return None;
            }
            let record = &p_side.features.halogens[a.feature as usize];
            let donor = angle(p_side.position(record.carbon), a.point, b.point)?;
            if donor < c.donor_angle {
                return None;
            }
            let mut angles = vec![(InteractionAngle::CarbonHalogenAcceptor, donor)];
            let acceptor = q_side.features.acceptors[b.feature as usize];
            for neighbor in q_side.heavy_neighbors(acceptor) {
                let value = angle(a.point, b.point, q_side.position(neighbor))?;
                if value < c.min_acceptor_angle || value > c.max_acceptor_angle {
                    return None;
                }
                angles.push((InteractionAngle::HalogenAcceptorNeighbor, value));
            }
            Some(build(InteractionKind::HalogenBond, angles, None))
        }
        (Role::Metal, Role::Ligand) => (distance <= criteria.metal_coordination.distance)
            .then(|| build(InteractionKind::MetalCoordination, Vec::new(), None)),
        (Role::Hydrophobic, Role::Hydrophobic) => (distance <= criteria.hydrophobic.distance)
            .then(|| build(InteractionKind::HydrophobicContact, Vec::new(), None)),
        _ => unreachable!("searched role pairs are listed above"),
    }
}

fn hydrogen_bond(
    criteria: &Canonical,
    donor_side: &Side<'_>,
    acceptor_side: &Side<'_>,
    a: &Candidate,
    b: &Candidate,
) -> Option<Interaction> {
    let c = criteria.hydrogen_bond;
    let distance = (a.point - b.point).norm();
    if distance > c.distance {
        return None;
    }
    let record = &donor_side.features.donors[a.feature as usize];
    let acceptor = acceptor_side.features.acceptors[b.feature as usize];
    let mut angles = Vec::new();
    let mut hydrogen = None;
    let mut hydrogen_distance = None;
    let hydrogens = donor_side.features.span(&record.hydrogens);
    if record.implicit == 0 && !hydrogens.is_empty() {
        // The qualifying hydrogen with the widest D–H…A angle.
        let mut best: Option<(f64, TopologyAtomIndex, f64)> = None;
        for &atom in hydrogens {
            let point = donor_side.position(atom);
            let reach = (point - b.point).norm();
            let Some(value) = angle(a.point, point, b.point) else {
                continue;
            };
            if reach <= c.hydrogen_distance
                && value >= c.hydrogen_angle
                && best.is_none_or(|(widest, _, _)| value > widest)
            {
                best = Some((value, atom, reach));
            }
        }
        let (value, atom, reach) = best?;
        angles.push((InteractionAngle::DonorHydrogenAcceptor, value));
        hydrogen = Some(atom);
        hydrogen_distance = Some(reach);
    } else {
        for neighbor in donor_side.heavy_neighbors(record.atom) {
            let value = angle(donor_side.position(neighbor), a.point, b.point)?;
            if value < c.donor_angle {
                return None;
            }
            angles.push((InteractionAngle::NeighborDonorAcceptor, value));
        }
    }
    for neighbor in acceptor_side.heavy_neighbors(acceptor) {
        let value = angle(a.point, b.point, acceptor_side.position(neighbor))?;
        if value < c.acceptor_angle {
            return None;
        }
        angles.push((InteractionAngle::DonorAcceptorNeighbor, value));
    }
    let mut first = end(donor_side, Role::Donor, a);
    first.hydrogen = hydrogen;
    Some(Interaction {
        kind: InteractionKind::HydrogenBond,
        first,
        second: end(acceptor_side, Role::Acceptor, b),
        distance,
        hydrogen_distance,
        offset: None,
        angles,
    })
}
