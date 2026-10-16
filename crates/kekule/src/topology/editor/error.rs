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
        match self {
            Self::EmptyTopology => f.write_str("cannot publish an empty topology"),
            Self::InvalidAtom(id) => write!(f, "deleted or foreign atom handle: {id}"),
            Self::InvalidBond(id) => write!(f, "deleted or foreign bond handle: {id}"),
            Self::InvalidChain(id) => write!(f, "deleted or foreign chain handle: {id}"),
            Self::InvalidResidue(id) => write!(f, "deleted or foreign residue handle: {id}"),
            Self::InvalidAtomSite(id) => write!(f, "deleted or foreign atom-site handle: {id}"),
            Self::InvalidSourceAtom(id) => {
                write!(f, "source atom has no live editing handle: {id}")
            }
            Self::InvalidSourceBond(id) => {
                write!(f, "source bond has no live editing handle: {id}")
            }
            Self::InvalidSourceInstance(id) => write!(f, "invalid source molecule instance: {id}"),
            Self::InvalidSourceChain(id) => {
                write!(f, "source chain has no live editing handle: {id}")
            }
            Self::InvalidSourceResidue(id) => {
                write!(f, "source residue has no live editing handle: {id}")
            }
            Self::InvalidSourceAtomSite(id) => {
                write!(f, "source site has no live editing handle: {id}")
            }
            Self::DuplicateAtomPlacement(id) => {
                write!(f, "atom already has a hierarchy site: {id}")
            }
            Self::Molecule(_) => f.write_str("cannot edit molecule"),
            Self::Publication(_) => f.write_str("cannot publish edited molecule"),
            Self::Topology(_) => f.write_str("cannot publish edited topology"),
            Self::Hierarchy(_) => f.write_str("cannot edit hierarchy"),
            Self::Property(_) => f.write_str("cannot edit properties"),
            Self::AppendProperty { domain, .. } => write!(f, "cannot append {domain} properties"),
        }
    }
}
impl std::error::Error for TopologyEditError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Molecule(e) => Some(e),
            Self::Publication(e) => Some(e),
            Self::Topology(e) => Some(e),
            Self::Hierarchy(e) => Some(e),
            Self::Property(e) => Some(e),
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
