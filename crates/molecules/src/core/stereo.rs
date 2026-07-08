use super::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoElement {
    pub kind: StereoElementKind,
    pub specifiedness: StereoSpecifiedness,
    pub source: StereoSource,
    pub group: Option<StereoGroupId>,
    pub descriptor: Option<StereoDescriptor>,
}

impl StereoElement {
    pub fn specified(kind: StereoElementKind, source: StereoSource) -> Self {
        Self {
            kind,
            specifiedness: StereoSpecifiedness::Specified,
            source,
            group: None,
            descriptor: None,
        }
    }

    pub fn references_atom(&self, atom: AtomId) -> bool {
        self.kind.references_atom(atom)
    }

    pub fn references_bond(&self, bond: BondId) -> bool {
        self.kind.references_bond(bond)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StereoElementKind {
    Tetrahedral(TetrahedralStereo),
    DoubleBond(DoubleBondStereo),
    Axis(AxisStereo),
}

impl StereoElementKind {
    fn references_atom(&self, atom: AtomId) -> bool {
        match self {
            Self::Tetrahedral(stereo) => {
                stereo.center == atom
                    || stereo
                        .carriers
                        .iter()
                        .any(|carrier| matches!(carrier, StereoCarrier::Atom(id) if *id == atom))
            }
            Self::DoubleBond(stereo) => {
                stereo.left == atom
                    || stereo.right == atom
                    || matches!(stereo.left_carrier, StereoCarrier::Atom(id) if id == atom)
                    || matches!(stereo.right_carrier, StereoCarrier::Atom(id) if id == atom)
            }
            Self::Axis(stereo) => stereo
                .carriers
                .iter()
                .any(|carrier| matches!(carrier, StereoCarrier::Atom(id) if *id == atom)),
        }
    }

    fn references_bond(&self, bond: BondId) -> bool {
        match self {
            Self::Tetrahedral(_) => false,
            Self::DoubleBond(stereo) => stereo.bond == bond,
            Self::Axis(stereo) => stereo.axis == bond,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    pub center: AtomId,
    pub carriers: Vec<StereoCarrier>,
    pub orientation: TetrahedralOrientation,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleBondStereo {
    pub bond: BondId,
    pub left: AtomId,
    pub right: AtomId,
    pub left_carrier: StereoCarrier,
    pub right_carrier: StereoCarrier,
    pub orientation: DoubleBondOrientation,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AxisStereo {
    pub axis: BondId,
    pub carriers: Vec<StereoCarrier>,
    pub orientation: AxisOrientation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StereoCarrier {
    Atom(AtomId),
    ImplicitHydrogen,
    ImplicitLonePair,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TetrahedralOrientation {
    Clockwise,
    CounterClockwise,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DoubleBondOrientation {
    Together,
    Opposite,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AxisOrientation {
    Clockwise,
    CounterClockwise,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StereoSpecifiedness {
    Specified,
    Unknown,
    Unspecified,
    InvalidCleared,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StereoSource {
    Smiles,
    MolfileV2000,
    MolfileV3000,
    Coordinates2D,
    Coordinates3D,
    Reaction,
    User,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StereoDescriptor {
    R,
    S,
    LowerR,
    LowerS,
    E,
    Z,
    M,
    P,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    pub kind: StereoGroupKind,
    pub members: Vec<StereoElementId>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StereoGroupKind {
    Absolute,
    Relative,
    Racemic,
    And,
    Or,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoBondMark {
    pub bond: BondId,
    pub kind: StereoBondMarkKind,
    pub source: StereoSource,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StereoBondMarkKind {
    DirectionalUp,
    DirectionalDown,
    WedgeUp,
    WedgeDown,
    WedgeEither,
    DoubleBondEither,
}
