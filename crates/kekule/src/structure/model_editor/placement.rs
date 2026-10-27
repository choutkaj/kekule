//! Builder placement on a model draft; see [`crate::structure::placement`].
use super::{atom_row, checked_point, ModelEditor};
use crate::algorithms::{
    AddHydrogensOptions, AddHydrogensScope, AddedHydrogenOrigin, RetainedHydrogenReason,
};
use crate::diagnostics::{LimitExceeded, Resource};
use crate::geometry::{Point3, Vector3};
use crate::molecule::{Atom, BondOrder, Element, ImplicitHydrogens};
use crate::structure::placement::{
    free_sites, ideal_bond_length, name_prefix, next_free_number, perceive_domains, unit,
    AtomGeometry, AtomShape, AttachedAtom, AttachmentSite, DomainInput, NameAllocator,
    PlacedHydrogens, PlacedHydrogensReport, PlacementError, RemovedHydrogens,
    RemovedHydrogensReport, SiteQuery,
};
use crate::topology::edit::{EditAtomId, EditResidueId};
use crate::topology::AtomSiteMetadata;
use crate::topology::DraftHydrogenError;
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};
use std::collections::BTreeSet;

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
        let (direction, replaced) = match site {
            AttachmentSite::Hydrogen(hydrogen) => {
                if !self.is_replaceable_hydrogen(parent, hydrogen)? {
                    return Err(PlacementError::NotReplaceableHydrogen {
                        atom: parent,
                        hydrogen,
                    });
                }
                let direction = unit(self.point(hydrogen) - center)
                    .ok_or(PlacementError::DegenerateGeometry { atom: parent })?;
                (direction, Some(hydrogen))
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
                (direction, None)
            }
        };
        let taken = bond_valence_of(order).saturating_sub(usize::from(replaced.is_some()));
        self.slots
            .try_reserve_atoms(1)
            .map_err(super::ModelEditError::from)?;
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
        let element = atom.element;
        let id = self.add_atom(atom, position)?;
        let label = self.free_chain_label();
        let chain = self.add_chain(label.clone(), Some(label.clone()))?;
        let residue = self.add_residue(chain, NEW_RESIDUE, Some(1), Some("1".to_owned()), None)?;
        let name = format!("{}1", name_prefix(element));
        self.add_atom_site(
            residue,
            id,
            AtomSiteMetadata::default()
                .with_type_symbol(Some(name_prefix(element)))
                .with_label_asym_id(Some(label.clone()))
                .with_auth_asym_id(Some(label))
                .with_label_atom_id(Some(name.clone()))
                .with_auth_atom_id(Some(name)),
        )?;
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
    /// `atom` stagger against: a heavy atom with the most heavy neighbors
    /// first, then any, in neighbor order.
    fn stagger_reference(
        &self,
        neighbor: EditAtomId,
        atom: EditAtomId,
    ) -> Result<Option<Point3>, PlacementError> {
        let mut best: Option<(usize, EditAtomId)> = None;
        for candidate in self.topology.neighbors(neighbor)? {
            if candidate == atom {
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
