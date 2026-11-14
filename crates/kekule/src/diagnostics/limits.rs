use std::fmt;

use super::Count;

/// A quantity a bounded algorithm or reader counts against a limit.
///
/// The error carrying a [`LimitExceeded`] names the operation, so a resource
/// says only what was counted: [`Resource::Atoms`] in a ring perception error
/// counts the atoms ring perception was asked to analyse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Resource {
    /// Bytes of input text.
    InputBytes,
    /// Bytes of one input record or line.
    RecordBytes,
    /// Bytes of one value or token.
    ValueBytes,
    /// Input records.
    Records,
    /// Data blocks of a CIF document.
    Blocks,
    /// Rows of one table, such as `_atom_site` rows.
    Rows,
    /// Values of one table or document.
    Values,
    /// Atoms, or atom storage slots.
    Atoms,
    /// Combined atom and bond storage slots.
    GraphSlots,
    /// Pending atom labels of canonical labeling.
    Labels,
    /// Bonds, or bond storage slots.
    Bonds,
    /// Hydrogen atoms an operation would add.
    Hydrogens,
    /// Molecule components or instances.
    Components,
    /// Nodes of a query expression or of all expressions of a query.
    ExpressionNodes,
    /// Nesting depth of a query expression.
    ExpressionDepth,
    /// Nesting depth of recursive queries.
    RecursionDepth,
    /// Nesting depth of branches.
    BranchDepth,
    /// Ring closures, open at once or in total.
    RingClosures,
    /// Units of work of a bounded calculation, such as graph visits.
    Work,
    /// States of a search, such as partial mappings or labelings.
    SearchStates,
    /// Candidate cycles of ring perception.
    CandidateCycles,
    /// Path expansions of ring perception.
    PathExpansions,
    /// Equivalent shortest paths of ring perception.
    ShortestPaths,
    /// Ring size, in atoms.
    RingSize,
    /// Matches reported by a search.
    Matches,
    /// Permutations of symmetry-equivalent atoms.
    Permutations,
    /// Nodes of a CIP hierarchical digraph.
    DigraphNodes,
    /// Depth of a CIP hierarchical digraph.
    DigraphDepth,
    /// Resonance contributors.
    Structures,
    /// Periodic lattice images examined.
    LatticeImages,
    /// Candidate positions, such as solvent placement sites.
    Positions,
    /// Pairwise entries of a distance matrix.
    MatrixEntries,
    /// Residues, such as protein residues of a secondary-structure assignment.
    Residues,
    /// Frames of a trajectory.
    Frames,
    /// Fields or items of one record.
    Fields,
    /// Combinations of carriers or other choices an encoder considers.
    Combinations,
    /// Candidate pairs, such as query-target atom pairs or residue pairs.
    CandidatePairs,
    /// Beta ladders of a secondary-structure assignment.
    Ladders,
}

impl fmt::Display for Resource {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::InputBytes => "input bytes",
            Self::RecordBytes => "record bytes",
            Self::ValueBytes => "value bytes",
            Self::Records => "records",
            Self::Blocks => "data blocks",
            Self::Rows => "rows",
            Self::Values => "values",
            Self::Atoms => "atoms",
            Self::GraphSlots => "graph slots",
            Self::Labels => "pending labels",
            Self::Bonds => "bonds",
            Self::Hydrogens => "hydrogens",
            Self::Components => "components",
            Self::ExpressionNodes => "expression nodes",
            Self::ExpressionDepth => "expression depth",
            Self::RecursionDepth => "recursion depth",
            Self::BranchDepth => "branch depth",
            Self::RingClosures => "ring closures",
            Self::Work => "work",
            Self::SearchStates => "search states",
            Self::CandidateCycles => "candidate cycles",
            Self::PathExpansions => "path expansions",
            Self::ShortestPaths => "equivalent shortest paths",
            Self::RingSize => "ring size",
            Self::Matches => "matches",
            Self::Permutations => "permutations",
            Self::DigraphNodes => "digraph nodes",
            Self::DigraphDepth => "digraph depth",
            Self::Structures => "resonance structures",
            Self::LatticeImages => "lattice images",
            Self::Positions => "candidate positions",
            Self::MatrixEntries => "matrix entries",
            Self::Residues => "residues",
            Self::Frames => "frames",
            Self::Fields => "fields",
            Self::Combinations => "combinations",
            Self::CandidatePairs => "candidate pairs",
            Self::Ladders => "ladders",
        })
    }
}

impl Resource {
    /// The name of the limit on this resource, such as `work` in "work
    /// limit".
    const fn limit_name(self) -> &'static str {
        match self {
            Self::InputBytes => "input size",
            Self::RecordBytes => "record size",
            Self::ValueBytes => "value size",
            Self::Records => "record",
            Self::Blocks => "data block",
            Self::Rows => "row",
            Self::Values => "value",
            Self::Atoms => "atom",
            Self::GraphSlots => "graph slot",
            Self::Labels => "pending label",
            Self::Bonds => "bond",
            Self::Hydrogens => "hydrogen",
            Self::Components => "component",
            Self::ExpressionNodes => "expression node",
            Self::ExpressionDepth => "expression depth",
            Self::RecursionDepth => "recursion depth",
            Self::BranchDepth => "branch depth",
            Self::RingClosures => "ring closure",
            Self::Work => "work",
            Self::SearchStates => "search state",
            Self::CandidateCycles => "candidate cycle",
            Self::PathExpansions => "path expansion",
            Self::ShortestPaths => "shortest path",
            Self::RingSize => "ring size",
            Self::Matches => "match",
            Self::Permutations => "permutation",
            Self::DigraphNodes => "digraph node",
            Self::DigraphDepth => "digraph depth",
            Self::Structures => "resonance structure",
            Self::LatticeImages => "lattice image",
            Self::Positions => "candidate position",
            Self::MatrixEntries => "matrix entry",
            Self::Residues => "residue",
            Self::Frames => "frame",
            Self::Fields => "field",
            Self::Combinations => "combination",
            Self::CandidatePairs => "candidate pair",
            Self::Ladders => "ladder",
        }
    }

    /// The unit a count of this resource is given in.
    const fn unit(self) -> &'static str {
        match self {
            Self::InputBytes | Self::RecordBytes | Self::ValueBytes => "bytes",
            Self::Records => "records",
            Self::Blocks => "data blocks",
            Self::Rows => "rows",
            Self::Values => "values",
            Self::Atoms | Self::RingSize => "atoms",
            Self::GraphSlots => "slots",
            Self::Labels => "labels",
            Self::Bonds => "bonds",
            Self::Hydrogens => "hydrogens",
            Self::Components => "components",
            Self::ExpressionNodes | Self::DigraphNodes => "nodes",
            Self::ExpressionDepth
            | Self::RecursionDepth
            | Self::BranchDepth
            | Self::DigraphDepth => "levels",
            Self::RingClosures => "ring closures",
            Self::Work => "steps",
            Self::SearchStates => "states",
            Self::CandidateCycles => "cycles",
            Self::PathExpansions => "expansions",
            Self::ShortestPaths => "paths",
            Self::Matches => "matches",
            Self::Permutations => "permutations",
            Self::Structures => "structures",
            Self::LatticeImages => "images",
            Self::Positions => "positions",
            Self::MatrixEntries => "entries",
            Self::Residues => "residues",
            Self::Frames => "frames",
            Self::Fields => "fields",
            Self::Combinations => "combinations",
            Self::CandidatePairs => "pairs",
            Self::Ladders => "ladders",
        }
    }
}

/// A bounded operation stopped because it would exceed a resource limit.
///
/// `observed` is the value that exceeded `limit`. An operation that stops at
/// the first excess reports that first value, so `observed` may be smaller
/// than the input would eventually have required; it is always greater than
/// `limit`. No partial result accompanies the error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LimitExceeded {
    resource: Resource,
    limit: usize,
    observed: usize,
}

impl LimitExceeded {
    /// Reports that `observed` units of `resource` exceed `limit`.
    pub const fn new(resource: Resource, limit: usize, observed: usize) -> Self {
        Self {
            resource,
            limit,
            observed,
        }
    }

    /// Reports the first unit of `resource` beyond `limit`, for an operation
    /// that counts one unit at a time.
    pub const fn first_beyond(resource: Resource, limit: usize) -> Self {
        Self::new(resource, limit, limit.saturating_add(1))
    }

    /// Returns `Err` when `observed` exceeds `limit`.
    pub(crate) const fn check(
        resource: Resource,
        limit: usize,
        observed: usize,
    ) -> Result<(), Self> {
        if observed > limit {
            Err(Self::new(resource, limit, observed))
        } else {
            Ok(())
        }
    }

    /// The bounded resource.
    pub const fn resource(&self) -> Resource {
        self.resource
    }

    /// The configured or built-in bound.
    pub const fn limit(&self) -> usize {
        self.limit
    }

    /// The value that exceeded the bound.
    pub const fn observed(&self) -> usize {
        self.observed
    }

    /// The limit with the observed value and the bound, such as `work limit:
    /// 69,400,898 steps of 69,248,000`, for the error of an operation to name
    /// in a sentence of its own.
    ///
    /// ```
    /// use kekule::diagnostics::{LimitExceeded, Resource};
    ///
    /// let limit = LimitExceeded::new(Resource::Work, 69_248_000, 69_400_898);
    /// assert_eq!(
    ///     format!("Ring perception exceeded its {}.", limit.describe()),
    ///     "Ring perception exceeded its work limit: 69,400,898 steps of 69,248,000."
    /// );
    /// ```
    pub fn describe(&self) -> impl fmt::Display + '_ {
        DescribedLimit(self)
    }
}

struct DescribedLimit<'a>(&'a LimitExceeded);

impl fmt::Display for DescribedLimit<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let limit = self.0;
        let unit = limit.resource.unit();
        let unit = match (limit.observed, unit) {
            (1, "matches") => "match",
            (1, "entries") => "entry",
            (1, unit) => unit.strip_suffix('s').unwrap_or(unit),
            (_, unit) => unit,
        };
        write!(
            formatter,
            "{} limit: {} {unit} of {}",
            limit.resource.limit_name(),
            Count(limit.observed),
            Count(limit.limit)
        )
    }
}

/// Writes `Exceeded the work limit: 69,400,898 steps of 69,248,000.`; an
/// error that carries the limit names its operation with
/// [`LimitExceeded::describe`] instead.
impl fmt::Display for LimitExceeded {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "Exceeded the {}.", self.describe())
    }
}

impl std::error::Error for LimitExceeded {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::{DiagnosticKind, Format, FormatError};

    #[test]
    fn limits_report_resource_bound_and_first_excess() {
        assert_eq!(LimitExceeded::check(Resource::Atoms, 3, 3), Ok(()));
        let exceeded = LimitExceeded::check(Resource::Atoms, 3, 5).unwrap_err();
        assert_eq!(
            (exceeded.resource(), exceeded.limit(), exceeded.observed()),
            (Resource::Atoms, 3, 5)
        );
        assert_eq!(
            exceeded.to_string(),
            "Exceeded the atom limit: 5 atoms of 3."
        );
        let work = LimitExceeded::new(Resource::Work, 69_248_000, 69_400_898);
        assert_eq!(
            format!("Ring perception exceeded its {}.", work.describe()),
            "Ring perception exceeded its work limit: 69,400,898 steps of 69,248,000."
        );
        assert_eq!(
            LimitExceeded::first_beyond(Resource::Records, usize::MAX).observed(),
            usize::MAX
        );
        let error = FormatError::parse(Format::Sdf, DiagnosticKind::Limit, "too many records")
            .exceeding(LimitExceeded::first_beyond(Resource::Records, 2));
        assert_eq!(
            error.limit_exceeded(),
            Some(LimitExceeded::new(Resource::Records, 2, 3))
        );
    }
}
