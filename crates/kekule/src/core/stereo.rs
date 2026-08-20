use super::*;

/// A canonical local stereochemical assertion.
///
/// Absence of an element means stereo is not asserted. A present element whose
/// kind has no orientation is an explicit assertion of unknown configuration.
/// Parser provenance and source-format marks are deliberately not canonical
/// stereo payload.
///
/// Source provenance is not available on canonical stereo elements:
///
/// ```compile_fail
/// fn source_is_not_canonical(element: &kekule::core::StereoElement) {
///     let _ = element.source;
/// }
/// ```
///
/// Placeholder specifiedness is likewise not canonical state:
///
/// ```compile_fail
/// fn placeholder_state_is_not_canonical(element: &kekule::core::StereoElement) {
///     let _ = element.specifiedness;
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoElement {
    pub kind: StereoElementKind,
    /// Owning relation group, assigned only by [`Molecule::add_stereo_group`].
    ///
    /// [`Molecule::add_stereo_element`] rejects values with this field set, and
    /// [`Molecule::remove_stereo_element`] clears it on the returned value.
    pub group: Option<StereoGroupId>,
}

impl StereoElement {
    /// Creates an ungrouped stereo assertion for checked molecule insertion.
    ///
    /// [`Molecule::add_stereo_element`] canonicalizes carrier and endpoint
    /// conventions before storing it.
    pub fn new(kind: StereoElementKind) -> Self {
        Self { kind, group: None }
    }

    /// Whether this element asserts a concrete local configuration.
    pub fn is_specified(&self) -> bool {
        self.kind.is_specified()
    }

    /// Whether stereo is asserted while its local configuration is unknown.
    pub fn is_explicitly_unknown(&self) -> bool {
        !self.is_specified()
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
    pub fn is_specified(&self) -> bool {
        match self {
            Self::Tetrahedral(stereo) => stereo.orientation.is_some(),
            Self::DoubleBond(stereo) => stereo.orientation.is_some(),
            Self::Axis(stereo) => stereo.orientation.is_some(),
        }
    }

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
    /// `None` represents an explicit assertion of unknown configuration.
    pub orientation: Option<TetrahedralOrientation>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleBondStereo {
    pub bond: BondId,
    pub left: AtomId,
    pub right: AtomId,
    pub left_carrier: StereoCarrier,
    pub right_carrier: StereoCarrier,
    /// `None` represents an explicit assertion of unknown configuration.
    pub orientation: Option<DoubleBondOrientation>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AxisStereo {
    pub axis: BondId,
    pub carriers: Vec<StereoCarrier>,
    /// `None` represents an explicit assertion of unknown configuration.
    pub orientation: Option<AxisOrientation>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StereoCarrier {
    Atom(AtomId),
    ImplicitHydrogen,
    ImplicitLonePair,
}

impl StereoCarrier {
    pub(crate) const fn canonical_order_key(self) -> (u8, u32) {
        match self {
            Self::Atom(atom) => (0, atom.raw()),
            Self::ImplicitHydrogen => (1, 0),
            Self::ImplicitLonePair => (2, 0),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TetrahedralOrientation {
    Clockwise,
    CounterClockwise,
}

impl TetrahedralOrientation {
    pub(crate) const fn inverted(self) -> Self {
        match self {
            Self::Clockwise => Self::CounterClockwise,
            Self::CounterClockwise => Self::Clockwise,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DoubleBondOrientation {
    Together,
    Opposite,
}

impl DoubleBondOrientation {
    pub(crate) const fn inverted(self) -> Self {
        match self {
            Self::Together => Self::Opposite,
            Self::Opposite => Self::Together,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AxisOrientation {
    Clockwise,
    CounterClockwise,
}

impl AxisOrientation {
    pub(crate) const fn inverted(self) -> Self {
        match self {
            Self::Clockwise => Self::CounterClockwise,
            Self::CounterClockwise => Self::Clockwise,
        }
    }
}

pub(crate) fn has_tetrahedral_stereo(molecule: &Molecule, center: AtomId) -> bool {
    molecule.stereo_elements().any(|(_, element)| {
        matches!(
            &element.kind,
            StereoElementKind::Tetrahedral(stereo) if stereo.center == center
        )
    })
}

pub(crate) fn has_double_bond_stereo(molecule: &Molecule, bond: BondId) -> bool {
    molecule.stereo_elements().any(|(_, element)| {
        matches!(
            &element.kind,
            StereoElementKind::DoubleBond(stereo) if stereo.bond == bond
        )
    })
}

pub(crate) fn has_axis_stereo(molecule: &Molecule, axis: BondId) -> bool {
    molecule.stereo_elements().any(|(_, element)| {
        matches!(
            &element.kind,
            StereoElementKind::Axis(stereo) if stereo.axis == axis
        )
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StereoDescriptor {
    R,
    S,
    LowerR,
    LowerS,
    SeqTrans,
    SeqCis,
    E,
    Z,
    M,
    P,
    LowerM,
    LowerP,
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
