use super::*;

/// A checked structural operation or publication failed. No partial owner is published.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TopologyEditError {
    EmptyTopology,
    InvalidAtom(EditAtomId),
    InvalidBond(EditBondId),
    InvalidChain(EditChainId),
    InvalidResidue(EditResidueId),
    InvalidAtomSite(EditAtomSiteId),
    InvalidSourceAtom(InstanceAtomId),
    InvalidSourceBond(InstanceBondId),
    InvalidSourceInstance(MoleculeInstanceId),
    InvalidSourceChain(ChainId),
    InvalidSourceResidue(ResidueId),
    InvalidSourceAtomSite(AtomSiteId),
    DuplicateAtomPlacement(EditAtomId),
    Molecule(MoleculeError),
    Publication(MoleculePublicationError),
    Topology(TopologyBuildError),
    Hierarchy(HierarchyError),
    Property(PropertyError),
    AppendProperty {
        domain: &'static str,
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
            Self::Molecule(e) => write!(f, "cannot edit molecule: {e}"),
            Self::Publication(e) => write!(f, "cannot publish edited molecule: {e}"),
            Self::Topology(e) => write!(f, "cannot publish edited topology: {e}"),
            Self::Hierarchy(e) => write!(f, "cannot edit hierarchy: {e}"),
            Self::Property(e) => write!(f, "cannot edit properties: {e}"),
            Self::AppendProperty { domain, error } => {
                write!(f, "cannot append {domain} properties: {error}")
            }
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

/// Failed publication retaining the exact draft for repair, including its handles.
#[derive(Debug)]
pub struct TopologyFinishError {
    pub(super) error: Box<TopologyEditError>,
    pub(super) editor: Box<TopologyEditor>,
}
impl TopologyFinishError {
    pub fn error(&self) -> &TopologyEditError {
        &self.error
    }
    pub fn editor(&self) -> &TopologyEditor {
        &self.editor
    }
    pub fn into_editor(self) -> TopologyEditor {
        *self.editor
    }
}
impl fmt::Display for TopologyFinishError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.error.fmt(f)
    }
}
impl std::error::Error for TopologyFinishError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(self.error.as_ref())
    }
}
