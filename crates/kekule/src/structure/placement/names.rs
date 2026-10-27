//! Names for new atoms: element symbol and the next free number in a residue.
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
        })
        .flatten()
        .filter_map(|name| numbered(name, prefix))
        .max()
        .map_or(1, |largest| largest.saturating_add(1))
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
}
