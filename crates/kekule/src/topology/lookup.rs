//! Exact, scope-aware hierarchy lookup without identifier fallback.

use super::{AtomSiteView, ChainView, HierarchyError, HierarchyIdKind, ResidueView, Topology};

fn unique<T>(
    mut items: impl Iterator<Item = T>,
    kind: HierarchyIdKind,
    query: impl FnOnce() -> String,
) -> Result<T, HierarchyError> {
    let Some(first) = items.next() else {
        return Err(HierarchyError::NotFound {
            kind,
            query: query(),
        });
    };
    let remaining = items.count();
    if remaining > 0 {
        return Err(HierarchyError::Ambiguous {
            kind,
            query: query(),
            matches: remaining + 1,
        });
    }
    Ok(first)
}

impl Topology {
    /// Finds one residue by its full author address across matching chains.
    /// Author chain IDs may be shared by several label chains, so uniqueness
    /// is checked on the complete (chain, sequence, insertion) address.
    pub fn residue_by_author(
        &self,
        chain: &str,
        sequence: &str,
        insertion: Option<&str>,
    ) -> Result<ResidueView<'_>, HierarchyError> {
        unique(
            self.residues().filter(|residue| {
                residue.chain().author_id() == Some(chain)
                    && residue.author_seq_id() == Some(sequence)
                    && residue.insertion_code() == insertion
            }),
            HierarchyIdKind::Residue,
            || format!("author chain {chain:?} residue {sequence:?} insertion {insertion:?}"),
        )
    }

    /// Finds exactly one chain by its label identifier.
    pub fn chain_by_label(&self, label: &str) -> Result<ChainView<'_>, HierarchyError> {
        unique(
            self.chains().filter(|chain| chain.label_id() == label),
            HierarchyIdKind::Chain,
            || format!("label {label:?}"),
        )
    }

    /// Finds exactly one chain by its author identifier, without label fallback.
    pub fn chain_by_author(&self, author: &str) -> Result<ChainView<'_>, HierarchyError> {
        unique(
            self.chains()
                .filter(|chain| chain.author_id() == Some(author)),
            HierarchyIdKind::Chain,
            || format!("author {author:?}"),
        )
    }
}

impl<'a> ChainView<'a> {
    /// Finds exactly one residue with this label sequence number in this chain.
    pub fn residue_by_label(self, sequence: i32) -> Result<ResidueView<'a>, HierarchyError> {
        unique(
            self.residues()
                .filter(|residue| residue.label_seq_id() == Some(sequence)),
            HierarchyIdKind::Residue,
            || format!("{} label {sequence}", self.id()),
        )
    }

    /// Finds exactly one residue by author sequence identifier and insertion
    /// code in this chain. `None` requires no insertion code; it is not a wildcard.
    pub fn residue_by_author(
        self,
        sequence: &str,
        insertion: Option<&str>,
    ) -> Result<ResidueView<'a>, HierarchyError> {
        unique(
            self.residues().filter(|residue| {
                residue.author_seq_id() == Some(sequence) && residue.insertion_code() == insertion
            }),
            HierarchyIdKind::Residue,
            || format!("{} author {sequence:?} insertion {insertion:?}", self.id()),
        )
    }
}

impl<'a> ResidueView<'a> {
    /// Finds exactly one site by its label atom name in this residue.
    pub fn atom_site_by_label(self, name: &str) -> Result<AtomSiteView<'a>, HierarchyError> {
        unique(
            self.atom_sites()
                .filter(|site| site.metadata().label_atom_id.as_deref() == Some(name)),
            HierarchyIdKind::AtomSite,
            || format!("{} label atom {name:?}", self.id()),
        )
    }

    /// Finds exactly one site by its author atom name, without label fallback.
    pub fn atom_site_by_author(self, name: &str) -> Result<AtomSiteView<'a>, HierarchyError> {
        unique(
            self.atom_sites()
                .filter(|site| site.metadata().auth_atom_id.as_deref() == Some(name)),
            HierarchyIdKind::AtomSite,
            || format!("{} author atom {name:?}", self.id()),
        )
    }
}
