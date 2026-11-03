use super::*;

/// A checked structural operation or publication failed. No partial owner is published.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TopologyEditError {
    /// Publication was requested with no atoms left.
    EmptyTopology,
    /// An atom handle is deleted or belongs to another editor.
    InvalidAtom(EditAtomId),
    /// A bond handle is deleted or belongs to another editor.
    InvalidBond(EditBondId),
    /// A chain handle is deleted or belongs to another editor.
    InvalidChain(EditChainId),
    /// A residue handle is deleted or belongs to another editor.
    InvalidResidue(EditResidueId),
    /// An atom-site handle is deleted or belongs to another editor.
    InvalidAtomSite(EditAtomSiteId),
    /// A source-topology atom has no live editing handle.
    InvalidSourceAtom(InstanceAtomId),
    /// A source-topology bond has no live editing handle.
    InvalidSourceBond(InstanceBondId),
    /// A source molecule instance is absent or no longer matches its edited atoms.
    InvalidSourceInstance(MoleculeInstanceId),
    /// A source-topology chain has no live editing handle.
    InvalidSourceChain(ChainId),
    /// A source-topology residue has no live editing handle.
    InvalidSourceResidue(ResidueId),
    /// A source-topology atom site has no live editing handle.
    InvalidSourceAtomSite(AtomSiteId),
    /// The atom already has a hierarchy site.
    DuplicateAtomPlacement(EditAtomId),
    /// A molecule edit failed; see the source.
    Molecule(MoleculeError),
    /// An edited molecule could not be published; see the source.
    Publication(MoleculePublicationError),
    /// The edited topology could not be published; see the source.
    Topology(TopologyBuildError),
    /// A hierarchy edit failed; see the source.
    Hierarchy(HierarchyError),
    /// A property edit failed; see the source.
    Property(PropertyError),
    /// Property rows of an appended topology could not be copied into a table.
    AppendProperty {
        /// The table whose rows could not be appended.
        domain: crate::properties::PropertyDomain,
        /// The cause, also the source.
        error: Box<PropertyError>,
    },
}

impl fmt::Display for TopologyEditError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use crate::diagnostics::Capitalized;
        match self {
            Self::EmptyTopology => f.write_str("The edited topology has no atoms."),
            Self::InvalidAtom(id) => deleted_or_foreign(f, id),
            Self::InvalidBond(id) => deleted_or_foreign(f, id),
            Self::InvalidChain(id) => deleted_or_foreign(f, id),
            Self::InvalidResidue(id) => deleted_or_foreign(f, id),
            Self::InvalidAtomSite(id) => deleted_or_foreign(f, id),
            Self::InvalidSourceAtom(id) => no_live_handle(f, id),
            Self::InvalidSourceBond(id) => no_live_handle(f, id),
            Self::InvalidSourceInstance(id) => write!(
                f,
                "Molecule {id} is not in the editor or no longer matches its edited atoms."
            ),
            Self::InvalidSourceChain(id) => no_live_handle(f, id),
            Self::InvalidSourceResidue(id) => no_live_handle(f, id),
            Self::InvalidSourceAtomSite(id) => no_live_handle(f, id),
            Self::DuplicateAtomPlacement(id) => {
                write!(f, "{} already has a hierarchy site.", Capitalized(id))
            }
            Self::Molecule(_) => f.write_str("Cannot edit the molecule."),
            Self::Publication(_) => f.write_str("Cannot publish the edited molecule."),
            Self::Topology(_) => f.write_str("Cannot publish the edited topology."),
            Self::Hierarchy(_) => f.write_str("Cannot edit the hierarchy."),
            Self::Property(error) => error.fmt(f),
            Self::AppendProperty { domain, .. } => {
                write!(f, "Cannot append the {domain} properties.")
            }
        }
    }
}

fn deleted_or_foreign(f: &mut fmt::Formatter<'_>, id: &dyn fmt::Display) -> fmt::Result {
    write!(
        f,
        "{} is deleted or belongs to another editor.",
        crate::diagnostics::Capitalized(id)
    )
}

fn no_live_handle(f: &mut fmt::Formatter<'_>, id: &dyn fmt::Display) -> fmt::Result {
    write!(
        f,
        "{} of the source topology is deleted in the editor.",
        crate::diagnostics::Capitalized(id)
    )
}

impl std::error::Error for TopologyEditError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Molecule(e) => Some(e),
            Self::Publication(e) => Some(e),
            Self::Topology(e) => Some(e),
            Self::Hierarchy(e) => Some(e),
            Self::Property(e) => e.source(),
            Self::AppendProperty { error, .. } => Some(error.as_ref()),
            _ => None,
        }
    }
}
macro_rules! convert {
    ($ty:ty, $variant:ident) => {
        impl From<$ty> for TopologyEditError {
            fn from(e: $ty) -> Self {
                Self::$variant(e)
            }
        }
    };
}
convert!(MoleculeError, Molecule);
convert!(MoleculePublicationError, Publication);
convert!(TopologyBuildError, Topology);
convert!(HierarchyError, Hierarchy);
convert!(PropertyError, Property);
impl From<super::super::components::ComponentBuildError> for TopologyEditError {
    fn from(e: super::super::components::ComponentBuildError) -> Self {
        use super::super::components::ComponentBuildError as E;
        match e {
            E::Molecule(e) => e.into(),
            E::Publication(e) => e.into(),
            E::Property(e) => e.into(),
            E::Topology(e) => e.into(),
        }
    }
}
