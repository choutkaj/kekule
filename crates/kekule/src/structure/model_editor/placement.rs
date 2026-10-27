//! Builder placement on a model draft; see [`crate::structure::placement`].
use super::{atom_row, checked_point, ModelEditor};
use crate::algorithms::{
    AddHydrogensOptions, AddHydrogensScope, AddedHydrogenOrigin, RetainedHydrogenReason,
};
use crate::diagnostics::{LimitExceeded, Resource};
use crate::geometry::{Point3, Vector3};
use crate::molecule::{Atom, AtomId, BondOrder, Element, ImplicitHydrogens};
use crate::structure::placement::{
    free_sites, ideal_bond_length, name_prefix, next_free_number, perceive_domains, unit,
    AtomGeometry, AtomShape, AttachedAtom, AttachedFragment, AttachmentSite, DomainInput, Fragment,
    NameAllocator, PlacedFragment, PlacedHydrogens, PlacedHydrogensReport, PlacementError,
    RemovedHydrogens, RemovedHydrogensReport, Rotation, SiteQuery,
};
use crate::structure::Positions;
use crate::topology::edit::{EditAtomId, EditResidueId};
use crate::topology::AtomSiteMetadata;
use crate::topology::DraftHydrogenError;
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};
use std::collections::{BTreeMap, BTreeSet};

/// The residue a new molecule's atoms join.
const NEW_RESIDUE: &str = "UNL";

impl ModelEditor {
    /// The ideal geometry of `atom` in the draft, from its element, formal
    /// charge, radical electrons, bond orders, neighbors, and implicit
    /// hydrogens; see [`placement`](crate::structure::placement).
    ///
    /// # Errors
    ///
    /// [`PlacementError::Edit`] when `atom` names no live atom, and
    /// [`PlacementError::UnsupportedGeometry`] for more than six domains.
    pub fn atom_geometry(&self, atom: EditAtomId) -> Result<AtomGeometry, PlacementError> {
        let (domains, lone_pairs, neighbors, implicit) = self.domains(atom)?;
        let shape = match domains {
            0 => None,
            _ => Some(
                AtomShape::from_domains(domains)
                    .ok_or(PlacementError::UnsupportedGeometry { atom, domains })?,
            ),
        };
        Ok(AtomGeometry {
            shape,
            neighbors,
            implicit_hydrogens: implicit,
            lone_pairs,
        })
    }

    /// Unit directions from `atom` to where its next substituents go, as for
    /// its hydrogens: the free sites of its [`AtomGeometry`], those of its
    /// implicit hydrogens included, in order of preference (anti to its
    /// neighbor's main substituent first).
    ///
    /// # Errors
    ///
    /// [`PlacementError::Edit`] for an invalid handle,
    /// [`PlacementError::UnsupportedGeometry`] for more than six domains, and
    /// [`PlacementError::DegenerateGeometry`] when its neighbors give no
    /// orientation.
    pub fn free_directions(&self, atom: EditAtomId) -> Result<Vec<Vector3>, PlacementError> {
        let geometry = self.atom_geometry(atom)?;
        let Some(shape) = geometry.shape else {
            return Ok(Vec::new());
        };
        let neighbors = self.topology.neighbors(atom)?.collect::<Vec<_>>();
        let points = neighbors
            .iter()
            .map(|neighbor| self.point(*neighbor))
            .collect::<Vec<_>>();
        let reference = match neighbors.as_slice() {
            [neighbor] => self.stagger_reference(*neighbor, atom)?,
            _ => None,
        };
        free_sites(SiteQuery {
            center: self.point(atom),
            neighbors: &points,
            reference,
            shape,
            lone_pairs: geometry.lone_pairs,
        })
        .map_err(|_| PlacementError::DegenerateGeometry { atom })
    }

    /// Makes the implicit hydrogens of each distinct atom in `atoms` explicit,
    /// bonded at the ideal length along its free directions, as
    /// [`Molecule::add_hydrogens`](crate::molecule::Molecule::add_hydrogens)
    /// does for a molecule. Inferred counts come from the valence model on
    /// the draft's current graph, so a bond order or charge changed earlier in
    /// the same draft counts; no perception is needed. A fixed count becomes
    /// zero. Represented stereo naming an implicit hydrogen names the new atom.
    /// New hydrogens join their atom's residue under the next free names.
    ///
    /// # Errors
    ///
    /// Nothing changes on error: [`PlacementError::Edit`] for an invalid
    /// handle, [`PlacementError::InvalidValence`] for an atom over its valence,
    /// [`PlacementError::NoFreePosition`] when its sites cannot hold its
    /// hydrogens, [`PlacementError::StereoHydrogen`],
    /// [`PlacementError::LimitExceeded`] beyond
    /// [`AddHydrogensOptions::max_added_hydrogens`], and the geometry errors
    /// of [`Self::free_directions`].
    pub fn add_hydrogens(
        &mut self,
        atoms: impl IntoIterator<Item = EditAtomId>,
        options: AddHydrogensOptions,
    ) -> Result<PlacedHydrogensReport, PlacementError> {
        let revision = self.topology.structural_revision;
        let hydrogen = Element::from_atomic_number(1).expect("hydrogen is an element");
        let mut seen = BTreeSet::new();
        let mut plans = Vec::new();
        let mut total = 0usize;
        for atom in atoms {
            if !seen.insert(atom) {
                continue;
            }
            if self.topology.exceeds_valence(atom)? {
                return Err(PlacementError::InvalidValence { atom });
            }
            let (mut count, origin) = self.topology.materializable_hydrogens(atom)?;
            if origin == AddedHydrogenOrigin::Inferred
                && options.scope() == AddHydrogensScope::FixedOnly
            {
                count = 0;
            }
            self.topology
                .check_hydrogen_stereo(atom, count)
                .map_err(|error| draft_error(error, atom))?;
            let mut points = Vec::with_capacity(count);
            if count > 0 {
                let directions = self.free_directions(atom)?;
                if directions.len() < count {
                    return Err(PlacementError::NoFreePosition { atom });
                }
                let length = bond_length(self.element(atom), hydrogen, BondOrder::Single)?;
                let center = self.point(atom);
                points.extend(directions[..count].iter().map(|d| center + *d * length));
            }
            total = total.saturating_add(count);
            plans.push((atom, origin, points));
        }
        if total > options.max_added_hydrogens() {
            return Err(PlacementError::LimitExceeded(LimitExceeded::new(
                Resource::Hydrogens,
                options.max_added_hydrogens(),
                total,
            )));
        }
        self.slots
            .try_reserve_atoms(total)
            .map_err(super::ModelEditError::from)?;
        let mut names = NameAllocator::default();
        let mut report = Vec::with_capacity(plans.len());
        for (atom, origin, points) in plans {
            let added = self
                .topology
                .materialize_hydrogens(atom, points.len(), origin)?;
            self.sync_slots(revision);
            let hydrogens = added
                .iter()
                .map(|(hydrogen, _)| *hydrogen)
                .collect::<Vec<_>>();
            for (hydrogen, point) in hydrogens.iter().zip(points) {
                self.place_row(*hydrogen, point);
                self.join_residue_of(atom, *hydrogen, &mut names)?;
            }
            report.push(PlacedHydrogens {
                atom,
                hydrogens,
                origin,
            });
        }
        Ok(PlacedHydrogensReport { atoms: report })
    }

    /// Collapses the ordinary explicit hydrogens bonded to each distinct atom
    /// in `atoms` into its implicit count, as
    /// [`Molecule::remove_hydrogens`](crate::molecule::Molecule::remove_hydrogens)
    /// does: hydrogens with information a count cannot hold (isotopes, maps,
    /// charges, topology or structure atom properties, stereo roles that
    /// cannot become implicit) are kept and reported. The atom keeps inference
    /// when it reproduces the complete count and is fixed otherwise; no
    /// perception is needed. Collapse changes representation, not
    /// protonation: an atom over its valence keeps its total as a fixed
    /// count, and changing that count is a separate edit (delete the
    /// hydrogens and set the count, for example to inferred).
    ///
    /// # Errors
    ///
    /// Nothing changes on error: [`PlacementError::Edit`] for an invalid
    /// handle and [`PlacementError::HydrogenCountOverflow`].
    pub fn remove_hydrogens(
        &mut self,
        atoms: impl IntoIterator<Item = EditAtomId>,
    ) -> Result<RemovedHydrogensReport, PlacementError> {
        let revision = self.topology.structural_revision;
        let atoms = atoms.into_iter().collect::<Vec<_>>();
        let slots = &self.slots;
        let collapsed = self
            .topology
            .collapse_hydrogens(&atoms, |atom, bond| {
                let has_data = |table: &crate::properties::RawPropertyTable, slot| {
                    table.row_has_data(slot).unwrap_or(false)
                };
                if has_data(slots.properties.atoms().raw(), atom) {
                    Some(RetainedHydrogenReason::AtomProperties)
                } else if has_data(slots.properties.bonds().raw(), bond) {
                    Some(RetainedHydrogenReason::BondProperties)
                } else {
                    None
                }
            })
            .map_err(|error| match error {
                DraftHydrogenError::Topology(error) => error.into(),
                DraftHydrogenError::Stereo(atom) => PlacementError::StereoHydrogen { atom },
                DraftHydrogenError::Overflow(atom, count) => {
                    PlacementError::HydrogenCountOverflow { atom, count }
                }
            })?;
        self.sync_slots(revision);
        self.slots.clear_atom_rows(&collapsed.dead_atom_slots);
        self.slots.clear_bond_rows(&collapsed.dead_bond_slots);
        Ok(RemovedHydrogensReport {
            atoms: collapsed
                .parents
                .into_iter()
                .map(|entry| RemovedHydrogens {
                    atom: entry.parent,
                    removed: entry.removed,
                    retained: entry.retained,
                    hydrogens: entry.hydrogens,
                    implicit_hydrogens: entry.implicit_hydrogens,
                })
                .collect(),
        })
    }

    /// Bonds a new atom to `parent` with `order`, at the ideal length for
    /// their elements and that order, at `site`: in place of one of the
    /// parent's explicit hydrogens, which is deleted, or at the parent's free
    /// position nearest a direction. The new atom takes the replaced
    /// hydrogen's place in the parent's represented stereo, as a single bond
    /// takes an implicit hydrogen's; a fixed implicit count on the parent
    /// drops by the valence the bond takes from it. The atom joins the
    /// parent's residue under the next free name. Its own hydrogens are not
    /// added; [`Self::add_hydrogens`] completes them. Closing a ring is a
    /// plain [`Self::add_bond`].
    ///
    /// # Errors
    ///
    /// Nothing changes on error: [`PlacementError::Edit`] for an invalid
    /// handle, [`PlacementError::NotReplaceableHydrogen`],
    /// [`PlacementError::InvalidDirection`],
    /// [`PlacementError::NoFreePosition`],
    /// [`PlacementError::UnknownBondLength`], and the geometry errors of
    /// [`Self::free_directions`].
    pub fn attach_atom(
        &mut self,
        parent: EditAtomId,
        atom: Atom,
        order: BondOrder,
        site: AttachmentSite,
    ) -> Result<AttachedAtom, PlacementError> {
        let revision = self.topology.structural_revision;
        let parent_atom = self.topology.require_atom(parent)?.clone();
        let length = bond_length(parent_atom.element, atom.element, order)?;
        let center = self.point(parent);
        let (direction, replaced) = self.attachment_direction(parent, site)?;
        self.slots
            .try_reserve_atoms(1)
            .map_err(super::ModelEditError::from)?;
        let attached = self.attach_first(parent, parent_atom, atom, order, replaced, revision)?;
        self.place_row(attached.atom, center + direction * length);
        self.join_residue_of(parent, attached.atom, &mut NameAllocator::default())?;
        Ok(AttachedAtom {
            atom: attached.atom,
            bond: attached.bond,
            replaced,
        })
    }

    /// Adds `atom` at `position` (any length unit) as a new molecule, with a
    /// chain and a residue `UNL` of its own; the atom is named by its element
    /// and `1` (`C1`). [`Self::add_hydrogens`] then completes it, so a carbon
    /// becomes methane. The chain's label is the first of `A`, `B`, …, `Z`,
    /// `AA`, … that no draft chain uses as a label or author identifier.
    ///
    /// # Errors
    ///
    /// Nothing changes on error: [`PlacementError::Edit`] when the position is
    /// not a finite length or capacity is exhausted.
    pub fn place_atom(
        &mut self,
        atom: Atom,
        position: Quantity<Point3>,
    ) -> Result<EditAtomId, PlacementError> {
        checked_point(position, None)?;
        let id = self.add_atom(atom, position)?;
        let residue = self.new_residue()?;
        self.name_in(residue, id, &mut NameAllocator::default())?;
        Ok(id)
    }

    /// The next free name for an atom of `element` in `residue`, by the rule
    /// placement names new atoms with: the element's uppercase symbol and one
    /// more than the largest number following it in a label or author atom
    /// name of the residue.
    ///
    /// # Errors
    ///
    /// [`PlacementError::Edit`] when `residue` names no live residue.
    pub fn next_atom_name(
        &self,
        residue: EditResidueId,
        element: Element,
    ) -> Result<String, PlacementError> {
        self.topology.require_residue(residue)?;
        let prefix = name_prefix(element);
        let number = next_free_number(&self.topology, residue, &prefix);
        Ok(format!("{prefix}{number}"))
    }

    /// Renames `atom` by the placement rule for its current element, as after
    /// changing its element: its label and author atom names become
    /// [`Self::next_atom_name`] in its residue, and its type symbol the
    /// element's. Returns the new name, or `None` when the atom has no site.
    ///
    /// # Errors
    ///
    /// [`PlacementError::Edit`] when `atom` names no live atom.
    pub fn rename_atom(&mut self, atom: EditAtomId) -> Result<Option<String>, PlacementError> {
        let element = self.element_checked(atom)?;
        let Some(site) = self.topology.atom_site_for_atom(atom)? else {
            return Ok(None);
        };
        let site_value = self.topology.require_atom_site(site)?;
        let residue = site_value.residue();
        let metadata = site_value.metadata().clone();
        let name = self.next_atom_name(residue, element)?;
        self.set_atom_site_metadata(
            site,
            metadata
                .with_type_symbol(Some(name_prefix(element)))
                .with_label_atom_id(Some(name.clone()))
                .with_auth_atom_id(Some(name.clone())),
        )?;
        Ok(Some(name))
    }

    /// Attaches `fragment` to `parent` at `site`, bonding its attachment atom
    /// at the ideal length for their elements and the order of the fragment's
    /// bond to its leaving atom, which is left out with any branch beyond it.
    /// The fragment keeps its own ideal coordinates; its torsion about the new
    /// bond is, of a few staggered torsions, the one that keeps it clearest
    /// of the atoms around it: three 120° steps between tetrahedral ends,
    /// four 90° steps between trigonal ends (coplanar and perpendicular), and
    /// six 60° steps between a trigonal and a tetrahedral end (a substituent
    /// in the plane, or two straddling it) and on a linear parent, none on a
    /// linear fragment, starting with the fragment's main
    /// substituent anti to the parent's. When none keeps every atom 2 Å from
    /// the atoms around it, the clearest of 15° steps is taken instead. Sites, stereo, fixed counts, and names follow
    /// [`Self::attach_atom`]; every atom of the fragment joins the parent's
    /// residue. Choosing the torsion reads the draft's position rows once and
    /// checks the liveness of nearby rows only.
    ///
    /// # Errors
    ///
    /// Nothing changes on error; the errors of [`Self::attach_atom`].
    pub fn attach_fragment(
        &mut self,
        parent: EditAtomId,
        fragment: &Fragment,
        site: AttachmentSite,
    ) -> Result<AttachedFragment, PlacementError> {
        let revision = self.topology.structural_revision;
        let parent_atom = self.topology.require_atom(parent)?.clone();
        let molecule = fragment.molecule();
        let attachment = fragment.attachment_atom();
        let attachment_atom = molecule.atom(attachment).expect("live").clone();
        let order = fragment.attachment_order();
        let length = bond_length(parent_atom.element, attachment_atom.element, order)?;
        let (direction, replaced) = self.attachment_direction(parent, site)?;
        let center = self.point(parent);
        let anchor = center + direction * length;
        let origin = fragment_point(fragment, attachment);
        let toward_parent = unit(fragment_point(fragment, fragment.leaving_atom()) - origin)
            .ok_or(PlacementError::DegenerateGeometry { atom: parent })?;
        let align = Rotation::between(toward_parent, -direction);
        let kept = fragment.kept_atoms();

        // Turn the fragment's main substituent anti to the parent's, then
        // try the staggered torsions from there.
        let parent_reference = self.main_substituent(parent, replaced)?;
        let fragment_reference = main_fragment_substituent(fragment, &kept);
        let base = match (parent_reference, fragment_reference) {
            (Some(parent_reference), Some(fragment_reference)) => {
                let ours = align.apply(fragment_point(fragment, fragment_reference) - origin);
                let theirs = parent_reference - center;
                signed_angle(
                    perpendicular_part(ours, direction),
                    -perpendicular_part(theirs, direction),
                    direction,
                )
            }
            _ => 0.0,
        };
        let parent_shape = self.atom_geometry(parent)?.shape();
        let fragment_shape = fragment_attachment_shape(molecule, attachment);
        let steps: &[f64] = match (parent_shape, fragment_shape) {
            // A linear fragment turns into itself; a linear parent has no
            // substituents to stagger against, so the fragment turns freely.
            (_, Some(AtomShape::Linear | AtomShape::Terminal)) => &[0.0],
            (Some(AtomShape::Linear | AtomShape::Terminal), _) => {
                &[0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
            }
            (Some(AtomShape::TrigonalPlanar), Some(AtomShape::TrigonalPlanar)) => {
                &[0.0, 90.0, 180.0, 270.0]
            }
            (Some(AtomShape::TrigonalPlanar), _) | (_, Some(AtomShape::TrigonalPlanar)) => {
                &[0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
            }
            _ => &[0.0, 120.0, 240.0],
        };
        let reach = kept
            .iter()
            .map(|(_, point)| (*point - origin).norm())
            .fold(0.0, f64::max)
            + CLEARANCE_REACH;
        let environment = self.environment(anchor, reach, &[Some(parent), replaced]);
        let clearance = |rotation: &Rotation| {
            kept.iter()
                .filter(|(atom, _)| *atom != attachment)
                .map(|(_, point)| anchor + rotation.apply(*point - origin))
                .flat_map(|point| environment.iter().map(move |other| (point - *other).norm()))
                .fold(f64::INFINITY, f64::min)
        };
        let clearest = |angles: &mut dyn Iterator<Item = f64>| {
            let mut best: Option<(f64, Rotation)> = None;
            for angle in angles {
                let rotation = Rotation::about(direction, base + angle.to_radians()).after(&align);
                let score = clearance(&rotation);
                if best
                    .as_ref()
                    .is_none_or(|(clearest, _)| score > *clearest + CLEARANCE_TIE)
                {
                    best = Some((score, rotation));
                }
            }
            best.expect("at least one torsion")
        };
        let mut best = clearest(&mut steps.iter().copied());
        // When no staggered torsion keeps clear, a finer turn may.
        if best.0 < CLEAR_CONTACT && steps.len() > 1 {
            let finer = clearest(&mut (0..24).map(|step| f64::from(step) * 15.0));
            if finer.0 > best.0 + CLEARANCE_TIE {
                best = finer;
            }
        }
        let rotation = best.1;

        self.slots
            .try_reserve_atoms(kept.len())
            .map_err(super::ModelEditError::from)?;
        let first = self.attach_first(
            parent,
            parent_atom,
            attachment_atom,
            order,
            replaced,
            revision,
        )?;
        let mut atoms = Vec::with_capacity(kept.len());
        for (local, point) in &kept {
            let id = if *local == attachment {
                first.atom
            } else {
                let atom = molecule.atom(*local).expect("live").clone();
                self.topology.add_atom_beside(parent, atom)?
            };
            atoms.push((*local, id));
            self.sync_slots(revision);
            self.place_row(id, anchor + rotation.apply(*point - origin));
        }
        let handles = atoms.iter().copied().collect::<BTreeMap<_, _>>();
        for (_, bond) in molecule.bonds() {
            if let (Some(a), Some(b)) = (handles.get(&bond.a()), handles.get(&bond.b())) {
                self.topology.add_bond(*a, *b, bond.order)?;
            }
        }
        self.sync_slots(revision);
        let mut names = NameAllocator::default();
        for (_, atom) in &atoms {
            self.join_residue_of(parent, *atom, &mut names)?;
        }
        Ok(AttachedFragment {
            atoms,
            attachment: first.atom,
            bond: first.bond,
            replaced,
        })
    }

    /// Adds `fragment`'s whole capped molecule as a new molecule, its atoms
    /// centered on `position` (any length unit) in the library's orientation,
    /// with a chain and a residue `UNL` of its own whose atoms are named by
    /// element and number in molecule order (`C1`, …, `H1`, …).
    ///
    /// # Errors
    ///
    /// Nothing changes on error: [`PlacementError::Edit`] when the position is
    /// not a finite length or capacity is exhausted.
    pub fn place_fragment(
        &mut self,
        fragment: &Fragment,
        position: Quantity<Point3>,
    ) -> Result<PlacedFragment, PlacementError> {
        let target = checked_point(position, None)?;
        let points = fragment.points();
        let count = points.len() as f64;
        let centroid = points.iter().fold(Vector3::zero(), |sum, point| {
            sum + (*point - Point3::origin())
        }) / count;
        let shift = target - (Point3::origin() + centroid);
        let positions = Positions::new(Quantity::new(
            points
                .iter()
                .map(|point| *point + shift)
                .collect::<Vec<_>>(),
            CANONICAL_LENGTH_UNIT,
        ))
        .map_err(super::ModelEditError::from)?;
        let added = self.add_molecule(fragment.molecule().clone(), &positions)?;
        let residue = self.new_residue()?;
        let mut names = NameAllocator::default();
        let atoms = added
            .atoms()
            .iter()
            .map(|(local, id)| (*local, *id))
            .collect::<Vec<_>>();
        for (_, atom) in &atoms {
            self.name_in(residue, *atom, &mut names)?;
        }
        Ok(PlacedFragment { atoms })
    }

    /// The positions of live draft atoms within `reach` of `center`, except
    /// `excluded`. Reads the position rows once, without allocating per row;
    /// liveness is checked only for nearby rows, against the draft's deleted
    /// slots.
    fn environment(
        &self,
        center: Point3,
        reach: f64,
        excluded: &[Option<EditAtomId>],
    ) -> Vec<Point3> {
        let excluded = excluded
            .iter()
            .flatten()
            .filter_map(|atom| self.topology.atom_slot(*atom).ok())
            .collect::<Vec<_>>();
        let rows = self.slots.positions.values();
        rows.value()
            .iter()
            .enumerate()
            .take(self.topology.atom_slot_count())
            .filter(|(slot, point)| {
                (**point - center).norm_squared() <= reach * reach
                    && !excluded.contains(slot)
                    && self.topology.is_live_atom_slot(*slot)
            })
            .map(|(_, point)| *point)
            .collect()
    }

    /// The direction from `parent` to a new substituent at `site`, and the
    /// hydrogen it replaces.
    fn attachment_direction(
        &self,
        parent: EditAtomId,
        site: AttachmentSite,
    ) -> Result<(Vector3, Option<EditAtomId>), PlacementError> {
        let center = self.point(parent);
        match site {
            AttachmentSite::Hydrogen(hydrogen) => {
                if !self.is_replaceable_hydrogen(parent, hydrogen)? {
                    return Err(PlacementError::NotReplaceableHydrogen {
                        atom: parent,
                        hydrogen,
                    });
                }
                let direction = unit(self.point(hydrogen) - center)
                    .ok_or(PlacementError::DegenerateGeometry { atom: parent })?;
                Ok((direction, Some(hydrogen)))
            }
            AttachmentSite::Toward(toward) => {
                let toward = unit(toward)
                    .filter(|toward| toward.is_finite())
                    .ok_or(PlacementError::InvalidDirection)?;
                let direction = self
                    .free_directions(parent)?
                    .into_iter()
                    .max_by(|a, b| a.dot(toward).total_cmp(&b.dot(toward)))
                    .ok_or(PlacementError::NoFreePosition { atom: parent })?;
                Ok((direction, None))
            }
        }
    }

    /// Bonds `atom` to `parent` in its occurrence, replacing `replaced`, and
    /// consumes a fixed implicit count on the parent by the valence the bond
    /// takes from it. Its position is the caller's.
    fn attach_first(
        &mut self,
        parent: EditAtomId,
        parent_atom: Atom,
        atom: Atom,
        order: BondOrder,
        replaced: Option<EditAtomId>,
        revision: u64,
    ) -> Result<crate::topology::AttachedInOccurrence, PlacementError> {
        let taken = bond_valence_of(order).saturating_sub(usize::from(replaced.is_some()));
        let attached = self
            .topology
            .attach_in_occurrence(parent, atom, order, replaced)?;
        self.sync_slots(revision);
        if let Some((atom_slot, bond_slots)) = &attached.dead {
            self.slots.clear_atom_rows(&[*atom_slot]);
            self.slots.clear_bond_rows(bond_slots);
        }
        if let ImplicitHydrogens::Fixed(count) = parent_atom.hydrogens {
            if taken > 0 && count > 0 {
                let mut updated = parent_atom;
                updated.hydrogens = ImplicitHydrogens::Fixed(
                    count.saturating_sub(u8::try_from(taken).unwrap_or(u8::MAX)),
                );
                self.topology.replace_atom(parent, updated)?;
            }
        }
        Ok(attached)
    }

    /// The position of `atom`'s main substituent other than `except`: a
    /// heavy atom with the most heavy neighbors first.
    fn main_substituent(
        &self,
        atom: EditAtomId,
        except: Option<EditAtomId>,
    ) -> Result<Option<Point3>, PlacementError> {
        let mut best: Option<(usize, EditAtomId)> = None;
        for candidate in self.topology.neighbors(atom)? {
            if Some(candidate) == except {
                continue;
            }
            let rank = if self.element(candidate).atomic_number() == 1 {
                0
            } else {
                1 + self
                    .topology
                    .neighbors(candidate)?
                    .filter(|other| self.element(*other).atomic_number() != 1)
                    .count()
            };
            if best.is_none_or(|(score, _)| rank > score) {
                best = Some((rank, candidate));
            }
        }
        Ok(best.map(|(_, candidate)| self.point(candidate)))
    }

    /// The domain count, site-taking lone pairs, neighbor count, and implicit
    /// hydrogens of `atom`.
    fn domains(&self, atom: EditAtomId) -> Result<(usize, usize, usize, usize), PlacementError> {
        let value = self.topology.require_atom(atom)?;
        let (implicit, _) = self.topology.materializable_hydrogens(atom)?;
        let mut neighbors = 0;
        let mut bond_valence = implicit;
        let mut coordinate_bonds = 0;
        let mut multiple = false;
        let mut conjugated = false;
        for bond in self.topology.incident_bonds(atom)? {
            let bond = self.topology.require_bond(bond)?;
            neighbors += 1;
            bond_valence += bond_valence_of(bond.order);
            if matches!(bond.order, BondOrder::Zero | BondOrder::Dative) {
                coordinate_bonds += 1;
            }
            multiple |= bond_valence_of(bond.order) > 1;
            let other = if bond.a() == atom { bond.b() } else { bond.a() };
            if !conjugated
                && value.element.atomic_number() == 7
                && (5..=8).contains(&self.element(other).atomic_number())
            {
                conjugated = self
                    .topology
                    .incident_bonds(other)?
                    .filter_map(|bond| self.topology.bond(bond))
                    .any(|bond| bond_valence_of(bond.order) > 1);
            }
        }
        let (domains, lone_pairs) = perceive_domains(DomainInput {
            element: value.element,
            formal_charge: value.formal_charge,
            radical_electrons: value
                .radical
                .map_or(0, |radical| usize::from(radical.electron_count())),
            substituents: neighbors + implicit,
            bond_valence,
            coordinate_bonds,
            conjugated_nitrogen: conjugated && !multiple && value.formal_charge == 0,
        });
        Ok((domains, lone_pairs, neighbors, implicit))
    }

    /// The substituent of `neighbor` (other than `atom`) that new sites of
    /// `atom` stagger against.
    fn stagger_reference(
        &self,
        neighbor: EditAtomId,
        atom: EditAtomId,
    ) -> Result<Option<Point3>, PlacementError> {
        self.main_substituent(neighbor, Some(atom))
    }

    fn is_replaceable_hydrogen(
        &self,
        parent: EditAtomId,
        hydrogen: EditAtomId,
    ) -> Result<bool, PlacementError> {
        let Some(atom) = self.topology.atom(hydrogen) else {
            return Ok(false);
        };
        if atom.element.atomic_number() != 1 {
            return Ok(false);
        }
        let bonds = self.topology.incident_bonds(hydrogen)?.collect::<Vec<_>>();
        Ok(match bonds.as_slice() {
            [bond] => self.topology.require_bond(*bond).is_ok_and(|bond| {
                bond.order == BondOrder::Single && (bond.a() == parent || bond.b() == parent)
            }),
            _ => false,
        })
    }

    /// Places `atom` in the residue of `parent`, if any, under the next free
    /// name for its element, with the parent's chain identifiers.
    fn join_residue_of(
        &mut self,
        parent: EditAtomId,
        atom: EditAtomId,
        names: &mut NameAllocator,
    ) -> Result<(), PlacementError> {
        let Some(site) = self.topology.atom_site_for_atom(parent)? else {
            return Ok(());
        };
        let parent_site = self.topology.require_atom_site(site)?;
        let residue = parent_site.residue();
        let element = self.element(atom);
        let name = names.allocate(&self.topology, residue, element);
        let metadata = parent_site
            .metadata()
            .clone()
            .with_type_symbol(Some(name_prefix(element)))
            .with_label_atom_id(Some(name.clone()))
            .with_auth_atom_id(Some(name));
        self.topology.add_atom_site(residue, atom, metadata)?;
        Ok(())
    }

    /// A new chain with a `UNL` residue, for a new molecule.
    fn new_residue(&mut self) -> Result<EditResidueId, PlacementError> {
        let label = self.free_chain_label();
        let chain = self.add_chain(label.clone(), Some(label))?;
        Ok(self.add_residue(chain, NEW_RESIDUE, Some(1), Some("1".to_owned()), None)?)
    }

    /// Places `atom` in `residue` under the next free name for its element,
    /// with the residue's chain identifiers.
    fn name_in(
        &mut self,
        residue: EditResidueId,
        atom: EditAtomId,
        names: &mut NameAllocator,
    ) -> Result<(), PlacementError> {
        let chain = self.topology.require_residue(residue)?.chain();
        let label = self.topology.require_chain(chain)?.label_id().to_owned();
        let element = self.element(atom);
        let name = names.allocate(&self.topology, residue, element);
        self.add_atom_site(
            residue,
            atom,
            AtomSiteMetadata::default()
                .with_type_symbol(Some(name_prefix(element)))
                .with_label_asym_id(Some(label.clone()))
                .with_auth_asym_id(Some(label))
                .with_label_atom_id(Some(name.clone()))
                .with_auth_atom_id(Some(name)),
        )?;
        Ok(())
    }

    /// The first chain label of `A`…`Z`, `AA`… unused as a label or author
    /// identifier.
    fn free_chain_label(&self) -> String {
        let used = self
            .topology
            .chains()
            .flat_map(|(_, chain)| [Some(chain.label_id()), chain.author_id()])
            .flatten()
            .map(str::to_owned)
            .collect::<BTreeSet<_>>();
        (0usize..)
            .map(|mut index| {
                let mut label = Vec::new();
                loop {
                    label.push(b'A' + u8::try_from(index % 26).expect("below 26"));
                    if index < 26 {
                        break;
                    }
                    index = index / 26 - 1;
                }
                label.reverse();
                String::from_utf8(label).expect("ASCII letters")
            })
            .find(|label| !used.contains(label))
            .expect("finitely many chains")
    }

    fn element_checked(&self, atom: EditAtomId) -> Result<Element, PlacementError> {
        Ok(self.topology.require_atom(atom)?.element)
    }

    fn element(&self, atom: EditAtomId) -> Element {
        self.topology.atom(atom).expect("live atom").element
    }

    /// The canonical draft position of a live atom.
    fn point(&self, atom: EditAtomId) -> Point3 {
        self.position(atom)
            .expect("live atoms have positions")
            .into_unit(CANONICAL_LENGTH_UNIT)
            .expect("positions are lengths")
            .into_value()
    }

    fn place_row(&mut self, atom: EditAtomId, point: Point3) {
        let slot = self.topology.atom_slot(atom).expect("live atom");
        self.slots
            .set_position(atom_row(slot), Quantity::new(point, CANONICAL_LENGTH_UNIT))
            .expect("placed positions are finite");
    }

    /// Grows structure rows to the topology's slots and, when the topology
    /// changed structurally since `revision`, clears structure owner
    /// annotations, as [`Self::structural`] does.
    fn sync_slots(&mut self, revision: u64) {
        self.slots.resize_slots(
            self.topology.atom_slot_count(),
            self.topology.bond_slot_count(),
        );
        if self.topology.structural_revision != revision {
            self.slots.properties.owner_mut().clear();
        }
    }
}

/// How far beyond a fragment's extent, in nm, atoms count toward choosing
/// its torsion.
const CLEARANCE_REACH: f64 = 0.4;
/// The clearance, in nm, below which a fragment's staggered torsions give
/// way to finer ones.
const CLEAR_CONTACT: f64 = 0.2;
/// Torsions whose clearances differ by less than this, in nm, tie; the
/// earlier wins.
const CLEARANCE_TIE: f64 = 1.0e-4;

fn fragment_point(fragment: &Fragment, atom: AtomId) -> Point3 {
    fragment
        .position(atom)
        .expect("a fragment atom")
        .into_value()
}

/// The fragment's main substituent on its attachment atom among `kept`: a
/// heavy atom with the most heavy neighbors first.
fn main_fragment_substituent(fragment: &Fragment, kept: &[(AtomId, Point3)]) -> Option<AtomId> {
    let molecule = fragment.molecule();
    let heavy = |atom: AtomId| {
        molecule
            .atom(atom)
            .is_some_and(|a| a.element.atomic_number() != 1)
    };
    let mut best: Option<(usize, AtomId)> = None;
    for candidate in molecule
        .neighbors(fragment.attachment_atom())
        .into_iter()
        .flatten()
    {
        if !kept.iter().any(|(atom, _)| *atom == candidate) {
            continue;
        }
        let rank = if heavy(candidate) {
            1 + molecule
                .neighbors(candidate)
                .into_iter()
                .flatten()
                .filter(|other| heavy(*other))
                .count()
        } else {
            0
        };
        if best.is_none_or(|(score, _)| rank > score) {
            best = Some((rank, candidate));
        }
    }
    best.map(|(_, atom)| atom)
}

/// The shape of a fragment's attachment atom in its capped molecule.
fn fragment_attachment_shape(
    molecule: &crate::molecule::Molecule,
    atom: AtomId,
) -> Option<AtomShape> {
    let value = molecule.atom(atom)?;
    let mut neighbors = 0;
    let mut bond_valence = 0;
    let mut coordinate_bonds = 0;
    let mut multiple = false;
    let mut conjugated = false;
    for (_, bond) in molecule.incident_bonds(atom).into_iter().flatten() {
        neighbors += 1;
        bond_valence += bond_valence_of(bond.order);
        coordinate_bonds += usize::from(matches!(bond.order, BondOrder::Zero | BondOrder::Dative));
        multiple |= bond_valence_of(bond.order) > 1;
        let other = bond.other_atom(atom);
        if value.element.atomic_number() == 7
            && molecule
                .atom(other)
                .is_some_and(|a| (5..=8).contains(&a.element.atomic_number()))
        {
            conjugated |= molecule
                .incident_bonds(other)
                .into_iter()
                .flatten()
                .any(|(_, bond)| bond_valence_of(bond.order) > 1);
        }
    }
    let (domains, _) = perceive_domains(DomainInput {
        element: value.element,
        formal_charge: value.formal_charge,
        radical_electrons: value
            .radical
            .map_or(0, |radical| usize::from(radical.electron_count())),
        substituents: neighbors,
        bond_valence,
        coordinate_bonds,
        conjugated_nitrogen: conjugated && !multiple && value.formal_charge == 0,
    });
    AtomShape::from_domains(domains)
}

/// The part of `v` perpendicular to the unit vector `axis`.
fn perpendicular_part(v: Vector3, axis: Vector3) -> Vector3 {
    v - axis * v.dot(axis)
}

/// The angle turning `from` onto `to` about the unit vector `axis`, both
/// perpendicular to it; zero when either vanishes.
fn signed_angle(from: Vector3, to: Vector3, axis: Vector3) -> f64 {
    if from.norm() < 1.0e-9 || to.norm() < 1.0e-9 {
        return 0.0;
    }
    from.cross(to).dot(axis).atan2(from.dot(to))
}

fn draft_error(error: DraftHydrogenError, atom: EditAtomId) -> PlacementError {
    match error {
        DraftHydrogenError::Topology(error) => error.into(),
        DraftHydrogenError::Stereo(_) => PlacementError::StereoHydrogen { atom },
        DraftHydrogenError::Overflow(atom, count) => {
            PlacementError::HydrogenCountOverflow { atom, count }
        }
    }
}

/// The ideal bond length in the canonical length unit.
fn bond_length(a: Element, b: Element, order: BondOrder) -> Result<f64, PlacementError> {
    ideal_bond_length(a, b, order)
        .and_then(|length| length.into_unit(CANONICAL_LENGTH_UNIT).ok())
        .map(Quantity::into_value)
        .ok_or(PlacementError::UnknownBondLength {
            first: a,
            second: b,
        })
}

fn bond_valence_of(order: BondOrder) -> usize {
    match order {
        BondOrder::Zero | BondOrder::Dative => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
    }
}
