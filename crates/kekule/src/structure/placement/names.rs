//! Names for new atoms: element symbol and the next free number in a residue,
//! or, for a fragment placed alone, its ring locants.
use super::Fragment;
use crate::molecule::Element;
use crate::topology::edit::{EditResidueId, TopologyEditor};
use std::collections::BTreeMap;

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

/// Allocates names within residues for one placement operation, reading each
/// residue's names once.
#[derive(Debug, Default)]
pub(crate) struct NameAllocator {
    next: BTreeMap<(EditResidueId, String), u64>,
}

impl NameAllocator {
    /// The next free name for `element` in `residue`, reserved for this
    /// allocator's later calls.
    pub(crate) fn allocate(
        &mut self,
        topology: &TopologyEditor,
        residue: EditResidueId,
        element: Element,
    ) -> String {
        let prefix = name_prefix(element);
        let next = self
            .next
            .entry((residue, prefix.clone()))
            .or_insert_with(|| next_free_number(topology, residue, &prefix));
        let name = format!("{prefix}{next}");
        *next = next.saturating_add(1);
        name
    }
}

/// One more than the largest number following `prefix` in a label or author
/// atom name of `residue`.
pub(crate) fn next_free_number(
    topology: &TopologyEditor,
    residue: EditResidueId,
    prefix: &str,
) -> u64 {
    let Some(atoms) = topology.residue_atom_ids(residue) else {
        return 1;
    };
    atoms
        .filter_map(|atom| {
            let site = topology.atom_site_for_atom(atom).ok()??;
            topology.atom_site(site)
        })
        .flat_map(|site| {
            let metadata = site.metadata();
            [metadata.label_atom_id(), metadata.auth_atom_id()]
                .map(|name| name.and_then(|name| numbered(name, prefix)))
        })
        .flatten()
        .max()
        .map_or(1, |largest| largest.saturating_add(1))
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
mod tests {
    use super::numbered;

    #[test]
    fn only_symbol_and_digits_are_numbered() {
        assert_eq!(numbered("C21", "C"), Some(21));
        assert_eq!(numbered("c7", "C"), Some(7));
        assert_eq!(numbered("CL1", "C"), None);
        assert_eq!(numbered("CA", "C"), None);
        assert_eq!(numbered("C", "C"), None);
        assert_eq!(numbered("CL12", "CL"), Some(12));
        assert_eq!(numbered("H", "HE"), None);
    }

    #[test]
    fn fragments_placed_alone_take_their_locants_and_distinct_names() {
        use super::super::{fragment, fragments, FragmentKind};
        use super::fragment_names;
        let names = |name: &str| fragment_names(fragment(name).unwrap()).join(" ");
        assert_eq!(names("benzene"), "C1 C2 C3 C4 C5 C6 H1 H2 H3 H4 H5 H6");
        assert_eq!(names("pyridine"), "N1 C2 C3 C4 C5 C6 H2 H3 H4 H5 H6");
        assert_eq!(
            names("naphthalene"),
            "C1 C2 C3 C4 C4A C5 C6 C7 C8 C8A H1 H2 H3 H4 H5 H6 H7 H8"
        );
        assert_eq!(
            names("indole"),
            "N1 C2 C3 C3A C4 C5 C6 C7 C7A H1 H2 H3 H4 H5 H6 H7"
        );
        assert_eq!(
            names("piperidine"),
            "N1 C2 C3 C4 C5 C6 H1 H21 H22 H31 H32 H41 H42 H51 H52 H61 H62"
        );
        assert_eq!(names("methylamino"), "N1 C1 HN11 HN12 H11 H12 H13");
        for fragment in fragments() {
            let molecule = fragment.molecule();
            let placed = fragment_names(fragment);
            let distinct = placed.iter().collect::<std::collections::BTreeSet<_>>();
            assert_eq!(distinct.len(), placed.len(), "{}", fragment.name());
            assert!(placed.iter().all(|name| name.len() <= 4), "{placed:?}");
            for (atom, name) in molecule.atom_ids().zip(&placed) {
                let element = molecule.atom(atom).unwrap().element();
                match fragment.locant(atom) {
                    Some(locant) => assert_eq!(
                        *name,
                        format!("{}{}", element.symbol(), locant).to_ascii_uppercase()
                    ),
                    // Every heavy atom of a ring fragment is a ring atom.
                    None => assert!(
                        fragment.kind() == FragmentKind::Group || element.atomic_number() == 1,
                        "{} {name}",
                        fragment.name()
                    ),
                }
            }
        }
    }
}
