//! Exact, scope-aware hierarchy lookup without identifier fallback.

use std::fmt;

use super::{AtomSiteView, ChainView, HierarchyIdKind, ResidueView, Topology};

/// A request for one hierarchy match found zero or multiple matches.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum HierarchyLookupError {
    NotFound {
        kind: HierarchyIdKind,
        query: String,
    },
    Ambiguous {
        kind: HierarchyIdKind,
        query: String,
        matches: usize,
    },
}

impl fmt::Display for HierarchyLookupError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFound { kind, query } => write!(f, "no {kind} matches {query}"),
            Self::Ambiguous {
                kind,
                query,
                matches,
            } => write!(f, "{matches} {kind} matches for {query}"),
        }
    }
}

impl std::error::Error for HierarchyLookupError {}

fn unique<T>(
    mut items: impl Iterator<Item = T>,
    kind: HierarchyIdKind,
    query: impl FnOnce() -> String,
) -> Result<T, HierarchyLookupError> {
    let Some(first) = items.next() else {
        return Err(HierarchyLookupError::NotFound {
            kind,
            query: query(),
        });
    };
    let remaining = items.count();
    if remaining > 0 {
        return Err(HierarchyLookupError::Ambiguous {
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
    ) -> Result<ResidueView<'_>, HierarchyLookupError> {
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
    pub fn chain_by_label(&self, label: &str) -> Result<ChainView<'_>, HierarchyLookupError> {
        unique(
            self.chains().filter(|chain| chain.label_id() == label),
            HierarchyIdKind::Chain,
            || format!("label {label:?}"),
        )
    }

    /// Finds exactly one chain by its author identifier, without label fallback.
    pub fn chain_by_author(&self, author: &str) -> Result<ChainView<'_>, HierarchyLookupError> {
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
    pub fn residue_by_label(self, sequence: i32) -> Result<ResidueView<'a>, HierarchyLookupError> {
        unique(
            self.residues()
                .filter(|residue| residue.label_seq_id() == Some(sequence)),
            HierarchyIdKind::Residue,
            || format!("chain {} label {sequence}", self.id()),
        )
    }

    /// Finds exactly one residue by author sequence identifier and insertion
    /// code in this chain. `None` requires no insertion code; it is not a wildcard.
    pub fn residue_by_author(
        self,
        sequence: &str,
        insertion: Option<&str>,
    ) -> Result<ResidueView<'a>, HierarchyLookupError> {
        unique(
            self.residues().filter(|residue| {
                residue.author_seq_id() == Some(sequence) && residue.insertion_code() == insertion
            }),
            HierarchyIdKind::Residue,
            || {
                format!(
                    "chain {} author {sequence:?} insertion {insertion:?}",
                    self.id()
                )
            },
        )
    }
}

impl<'a> ResidueView<'a> {
    /// Finds exactly one site by its label atom name in this residue.
    pub fn atom_site_by_label(self, name: &str) -> Result<AtomSiteView<'a>, HierarchyLookupError> {
        unique(
            self.atom_sites()
                .filter(|site| site.metadata().label_atom_id.as_deref() == Some(name)),
            HierarchyIdKind::AtomSite,
            || format!("residue {} label atom {name:?}", self.id()),
        )
    }

    /// Finds exactly one site by its author atom name, without label fallback.
    pub fn atom_site_by_author(self, name: &str) -> Result<AtomSiteView<'a>, HierarchyLookupError> {
        unique(
            self.atom_sites()
                .filter(|site| site.metadata().auth_atom_id.as_deref() == Some(name)),
            HierarchyIdKind::AtomSite,
            || format!("residue {} author atom {name:?}", self.id()),
        )
    }
}
