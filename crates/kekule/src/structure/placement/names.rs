//! Names for new atoms: a template's hydrogen name, or element symbol and the
//! next free number in a residue; for a fragment placed alone, its ring
//! locants.
use super::Fragment;
use crate::components::{Component, ComponentLibrary};
use crate::molecule::Element;
use crate::topology::edit::{EditAtomId, EditResidueId, TopologyEditor};
use std::collections::{BTreeMap, BTreeSet};

/// The uppercase symbol new atoms of `element` are named by.
pub(crate) fn name_prefix(element: Element) -> String {
    element.symbol().to_ascii_uppercase()
}

/// The number following `prefix` in `name`, compared without case, when the
/// rest of the name is ASCII digits.
fn numbered(name: &str, prefix: &str) -> Option<u64> {
    let head = name.get(..prefix.len())?;
    let digits = &name[prefix.len()..];
    (head.eq_ignore_ascii_case(prefix)
        && !digits.is_empty()
        && digits.bytes().all(|byte| byte.is_ascii_digit()))
    .then(|| digits.parse().unwrap_or(u64::MAX))
}

/// The names of one residue: those its atoms take, and the component
/// templates by which readers place named atoms.
#[derive(Debug)]
pub(crate) struct ResidueNames<'c> {
    taken: BTreeSet<String>,
    /// The template of the residue's component, then its protonation
    /// variants.
    templates: Vec<&'c Component>,
}

impl<'c> ResidueNames<'c> {
    /// The names of a residue of the component `templates` define, whose
    /// atoms take `taken`.
    pub(crate) fn new(
        templates: Option<(&'c Component, Vec<&'c Component>)>,
        taken: impl IntoIterator<Item = String>,
    ) -> Self {
        Self {
            taken: taken.into_iter().collect(),
            templates: templates
                .map(|(component, variants)| std::iter::once(component).chain(variants).collect())
                .unwrap_or_default(),
        }
    }

    /// A free name for a new atom of `element` bonded to the atom named
    /// `parent`, taken for later calls. A hydrogen takes the first name the
    /// templates bond to `parent`, else a terminal hydrogen's name readers
    /// bond to it. Any other atom takes [`Self::numbered_name`].
    pub(crate) fn allocate(&mut self, element: Element, parent: Option<&str>) -> String {
        let name = parent
            .filter(|_| element.atomic_number() == 1)
            .and_then(|parent| self.hydrogen_on(parent))
            .map_or_else(|| self.numbered_name(element), str::to_owned);
        self.taken.insert(name.clone());
        name
    }

    /// Takes `name` for a new atom when no atom takes it and no template
    /// places it, and returns whether it did.
    pub(crate) fn claim(&mut self, name: &str) -> bool {
        let free = !self.is_taken(name) && !self.placed(name);
        if free {
            self.taken.insert(name.to_owned());
        }
        free
    }

    /// The element's uppercase symbol and one more than the largest number
    /// following it in a name of the residue, skipping names the templates
    /// place: readers then bond the atom by distance rather than by a name
    /// that means another atom.
    pub(crate) fn numbered_name(&self, element: Element) -> String {
        let prefix = name_prefix(element);
        let start = self
            .taken
            .iter()
            .filter_map(|name| numbered(name, &prefix))
            .max()
            .map_or(1, |largest| largest.saturating_add(1));
        (start..)
            .map(|number| format!("{prefix}{number}"))
            .find(|name| !self.is_taken(name) && !self.placed(name))
            .expect("numbers are unbounded")
    }

    /// Whether an atom takes `name`, or another name of the template atom it
    /// names, as `1HB` takes `HB2`: readers match both to one atom.
    fn is_taken(&self, name: &str) -> bool {
        self.taken.contains(name)
            || self.templates.iter().any(|template| {
                template.atom(name).is_some_and(|(index, _)| {
                    self.taken.iter().any(|taken| {
                        template
                            .atom(taken)
                            .is_some_and(|(other, _)| other == index)
                    })
                })
            })
    }

    /// Whether readers place an atom named `name` by the templates: an atom
    /// name or alternate name, or a terminal hydrogen name.
    fn placed(&self, name: &str) -> bool {
        self.templates.iter().any(|template| {
            template.atom(name).is_some()
                || template
                    .linking()
                    .terminal_hydrogens()
                    .iter()
                    .any(|(hydrogen, _)| *hydrogen == name)
        })
    }

    /// The first free hydrogen name that readers place on the atom named
    /// `parent`: one that every template naming it bonds to `parent`, else
    /// a terminal hydrogen's that no template names.
    fn hydrogen_on(&self, parent: &str) -> Option<&'c str> {
        let free = |name: &str| !self.is_taken(name);
        let unnamed = |name: &str| {
            self.templates
                .iter()
                .all(|template| template.atom(name).is_none())
        };
        let templated = self
            .templates
            .iter()
            .flat_map(|template| template.atoms())
            .filter(|atom| atom.element().atomic_number() == 1)
            .map(|atom| atom.name())
            .find(|&name| {
                free(name)
                    && self.templates.iter().all(|template| {
                        template.atom(name).is_none() || template.bonds_named(name, parent)
                    })
            });
        templated.or_else(|| {
            self.templates
                .first()?
                .linking()
                .terminal_hydrogens()
                .iter()
                .find(|(name, anchors)| anchors.contains(&parent) && free(name) && unnamed(name))
                .map(|(name, _)| *name)
        })
    }
}

/// Allocates names within residues for one placement operation, reading each
/// residue's names once.
#[derive(Debug, Default)]
pub(crate) struct NameAllocator {
    residues: BTreeMap<EditResidueId, ResidueNames<'static>>,
}

impl NameAllocator {
    /// A free name for `atom`, of `element`, in `residue` by
    /// [`ResidueNames::allocate`], reserved for this allocator's later
    /// calls. The atom's parent is its first heavy neighbor in the residue.
    pub(crate) fn allocate(
        &mut self,
        topology: &TopologyEditor,
        residue: EditResidueId,
        atom: EditAtomId,
        element: Element,
    ) -> String {
        let parent = parent_name(topology, residue, atom);
        self.residues
            .entry(residue)
            .or_insert_with(|| residue_names(topology, residue))
            .allocate(element, parent.as_deref())
    }
}

/// The names of `residue`: its label and author atom names, and the bundled
/// templates of its component.
pub(crate) fn residue_names(
    topology: &TopologyEditor,
    residue: EditResidueId,
) -> ResidueNames<'static> {
    let library = ComponentLibrary::bundled();
    let templates = topology.residue(residue).and_then(|residue| {
        let id = residue.label_comp_id().unwrap_or(residue.name());
        Some((library.get(id)?, library.variants(id).collect()))
    });
    let taken = topology
        .residue_atom_ids(residue)
        .into_iter()
        .flatten()
        .filter_map(|atom| {
            let site = topology.atom_site_for_atom(atom).ok()??;
            topology.atom_site(site)
        })
        .flat_map(|site| {
            let metadata = site.metadata();
            [metadata.label_atom_id(), metadata.auth_atom_id()].map(|name| name.map(str::to_owned))
        })
        .flatten()
        .collect::<Vec<_>>();
    ResidueNames::new(templates, taken)
}

/// The label or else author atom name of the first heavy atom bonded to
/// `atom` in `residue`.
fn parent_name(
    topology: &TopologyEditor,
    residue: EditResidueId,
    atom: EditAtomId,
) -> Option<String> {
    topology
        .neighbors(atom)
        .ok()?
        .filter(|&neighbor| {
            topology
                .atom(neighbor)
                .is_some_and(|neighbor| neighbor.element.atomic_number() != 1)
        })
        .find_map(|neighbor| {
            let site = topology.atom_site_for_atom(neighbor).ok()??;
            let site = topology.atom_site(site)?;
            let metadata = site.metadata();
            (site.residue() == residue)
                .then(|| metadata.label_atom_id().or(metadata.auth_atom_id()))
                .flatten()
                .map(str::to_owned)
        })
}

/// The names of a fragment's atoms when it is placed alone, in molecule atom
/// order: ring atoms by element and locant (`C4A`), other heavy atoms by
/// element and the next number after the ring's, and hydrogens after the
/// atom carrying them: `H` and its locant or number for one hydrogen, with a
/// running digit for several (`H2`; `H21`, `H22`). Where atoms of two
/// elements share a number, the hydrogens of the one that is not carbon
/// name its element too (`HN11`).
pub(crate) fn fragment_names(fragment: &Fragment) -> Vec<String> {
    let molecule = fragment.molecule();
    let atoms = molecule.atom_ids().collect::<Vec<_>>();
    let element = |atom| molecule.atom(atom).expect("fragment atom").element;
    let is_hydrogen = |atom| element(atom).atomic_number() == 1;
    let mut names = vec![String::new(); atoms.len()];
    // Heavy atoms: their tag is the locant or number after the symbol.
    let mut tags = vec![String::new(); atoms.len()];
    let mut next = BTreeMap::<String, u64>::new();
    for (index, &atom) in atoms.iter().enumerate() {
        if let Some(locant) = fragment.locant(atom) {
            let prefix = name_prefix(element(atom));
            tags[index] = locant.to_ascii_uppercase();
            names[index] = format!("{prefix}{}", tags[index]);
            if let Some(number) = numbered(&names[index], &prefix) {
                let entry = next.entry(prefix).or_insert(1);
                *entry = (*entry).max(number.saturating_add(1));
            }
        }
    }
    for (index, &atom) in atoms.iter().enumerate() {
        if names[index].is_empty() && !is_hydrogen(atom) {
            let prefix = name_prefix(element(atom));
            let number = next.entry(prefix.clone()).or_insert(1);
            tags[index] = number.to_string();
            names[index] = format!("{prefix}{number}");
            *number += 1;
        }
    }
    // Each hydrogen's carrier is its one heavy neighbor.
    let carrier = |atom| {
        molecule
            .neighbors(atom)
            .into_iter()
            .flatten()
            .find(|&neighbor| !is_hydrogen(neighbor))
            .and_then(|neighbor| atoms.iter().position(|&id| id == neighbor))
    };
    let carriers = atoms
        .iter()
        .map(|&atom| is_hydrogen(atom).then(|| carrier(atom)).flatten())
        .collect::<Vec<_>>();
    let mut counts = vec![0usize; atoms.len()];
    for &carrier in carriers.iter().flatten() {
        counts[carrier] += 1;
    }
    let shared = |index: usize| {
        (0..atoms.len())
            .any(|other| other != index && counts[other] > 0 && tags[other] == tags[index])
    };
    let mut seen = vec![0usize; atoms.len()];
    let mut loose = 0;
    for (index, carrier) in carriers.iter().enumerate() {
        if !is_hydrogen(atoms[index]) {
            continue;
        }
        let Some(carrier) = *carrier else {
            loose += 1;
            names[index] = format!("H{loose}");
            continue;
        };
        let element = element(atoms[carrier]);
        let tag = if shared(carrier) && element.atomic_number() != 6 {
            format!("{}{}", name_prefix(element), tags[carrier])
        } else {
            tags[carrier].clone()
        };
        seen[carrier] += 1;
        names[index] = if counts[carrier] == 1 {
            format!("H{tag}")
        } else {
            format!("H{tag}{}", seen[carrier])
        };
    }
    names
}

#[cfg(test)]
mod tests;
