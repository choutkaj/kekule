use super::*;

/// A canonical local stereochemical assertion.
///
/// Absence of an element means stereo is not asserted. A present element whose
/// kind has no orientation is an explicit assertion of unknown configuration.
/// Parser provenance and source-format marks are deliberately not canonical
/// stereo payload.
/// Each atom or bond focus has at most one assertion; changing a configuration
/// uses [`MoleculeEditor::replace_stereo_element`].
/// Structural equality compares the stored assertion and group ID. It does not
/// identify equivalent enhanced-group representatives under joint inversion;
/// canonical SMILES normalizes that equivalence for supported encodings.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoElement {
    /// The focus, carriers, and orientation of the assertion.
    pub kind: StereoElementKind,
    pub(crate) group: Option<StereoGroupId>,
}

impl StereoElement {
    /// Creates an ungrouped stereo assertion for checked molecule insertion.
    ///
    /// [`MoleculeEditor::add_stereo_element`] canonicalizes carrier and endpoint
    /// conventions before storing it.
    pub fn new(kind: StereoElementKind) -> Self {
        Self { kind, group: None }
    }

    /// Owning relation group, assigned only by [`MoleculeEditor::add_stereo_group`].
    ///
    /// An element read from a molecule reports its group;
    /// [`MoleculeEditor::add_stereo_element`] rejects a grouped element, and
    /// [`MoleculeEditor::remove_stereo_element`] returns the element ungrouped.
    pub const fn group(&self) -> Option<StereoGroupId> {
        self.group
    }

    /// Whether this element asserts a concrete local configuration.
    pub fn is_specified(&self) -> bool {
        self.kind.is_specified()
    }

    /// Whether stereo is asserted while its local configuration is unknown.
    pub fn is_explicitly_unknown(&self) -> bool {
        !self.is_specified()
    }

    /// Whether `atom` is this element's focus, an axis or double-bond
    /// endpoint, or one of its carriers.
    pub fn references_atom(&self, atom: AtomId) -> bool {
        self.kind.references_atom(atom)
    }

    /// Whether `bond` is this element's double-bond or axis focus.
    pub fn references_bond(&self, bond: BondId) -> bool {
        self.kind.references_bond(bond)
    }
}

/// The geometry of one stereo assertion and its focus.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum StereoElementKind {
    /// Tetrahedral configuration around a center atom.
    Tetrahedral(TetrahedralStereo),
    /// Cis/trans configuration about a double bond.
    DoubleBond(DoubleBondStereo),
    /// Axial configuration about a bond, as in an atropisomer.
    Axis(AxisStereo),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum StereoFocus {
    Atom(AtomId),
    Bond(BondId),
}

impl StereoElementKind {
    pub(crate) fn focus(&self) -> StereoFocus {
        match self {
            Self::Tetrahedral(stereo) => StereoFocus::Atom(stereo.center),
            Self::DoubleBond(stereo) => StereoFocus::Bond(stereo.bond),
            Self::Axis(stereo) => StereoFocus::Bond(stereo.axis),
        }
    }

    /// Whether the assertion has a concrete orientation; `false` for an
    /// explicitly unknown configuration.
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

/// Tetrahedral configuration around one center atom.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    /// The stereocenter.
    pub center: AtomId,
    /// Four distinct carriers in the order used by `orientation`. The explicit
    /// atom carriers are exactly the center's bonded neighbors.
    pub carriers: [StereoCarrier; 4],
    /// `None` represents an explicit assertion of unknown configuration.
    pub orientation: Option<TetrahedralOrientation>,
}

/// Cis/trans configuration about one double bond, stated for one reference
/// carrier at each endpoint.
///
/// Checked insertion requires a [`BondOrder::Double`] focus bond and stores
/// the endpoints in ascending ID order, exchanging the carriers with them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleBondStereo {
    /// The double bond.
    pub bond: BondId,
    /// One endpoint of `bond`.
    pub left: AtomId,
    /// The other endpoint of `bond`.
    pub right: AtomId,
    /// The reference at `left`: a neighbor other than `right`, or its implicit
    /// hydrogen.
    pub left_carrier: StereoCarrier,
    /// The reference at `right`: a neighbor other than `left`, or its implicit
    /// hydrogen.
    pub right_carrier: StereoCarrier,
    /// `None` represents an explicit assertion of unknown configuration.
    pub orientation: Option<DoubleBondOrientation>,
}

/// Axial configuration about one bond, stated for one reference atom at each
/// axis endpoint.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AxisStereo {
    /// The bond whose endpoints define the axis.
    pub axis: BondId,
    /// Two explicit atom references, one bonded exclusively to each axis
    /// endpoint. Checked insertion canonicalizes endpoint and reference order.
    pub carriers: [StereoCarrier; 2],
    /// `None` represents an explicit assertion of unknown configuration.
    pub orientation: Option<AxisOrientation>,
}

/// One reference around a stereo focus.
///
/// Implicit carriers stand for a non-graph hydrogen or a lone pair of the focus
/// atom they are listed for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum StereoCarrier {
    /// A graph atom bonded to the focus.
    Atom(AtomId),
    /// An implicit (non-graph) hydrogen of the focus atom.
    ImplicitHydrogen,
    /// A lone pair of the focus atom, as at a tetrahedral nitrogen or sulfur.
    /// Not accepted for double-bond or axis stereo.
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

/// The handedness of four tetrahedral carriers in their stored order.
///
/// This enum is closed: a tetrahedral arrangement has exactly two
/// handednesses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TetrahedralOrientation {
    /// Positive signed volume `(p0 - p3) × (p1 - p3) · (p2 - p3)`
    /// for carrier coordinates in their stored order.
    Clockwise,
    /// Negative signed volume for the stored carrier order.
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

/// Whether two double-bond reference carriers lie on the same side.
///
/// This enum is closed: two references are either together or opposite.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DoubleBondOrientation {
    /// The two references are on the same side (cis).
    Together,
    /// The two references are on opposite sides (trans).
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

/// The handedness of an axis and its endpoint references.
///
/// This enum is closed: the signed triple product has two signs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AxisOrientation {
    /// Positive `(right - left) · ((left_reference - left) ×
    /// (right_reference - right))`. Exchanging both endpoints and their
    /// references preserves this handedness.
    Clockwise,
    /// Negative signed triple product for the endpoint references.
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

/// A CIP stereodescriptor assigned to one stereo element.
///
/// Descriptors are derived perception, installed by CIP assignment; they are
/// never part of represented stereo.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum StereoDescriptor {
    /// `R` at a stereogenic center.
    R,
    /// `S` at a stereogenic center.
    S,
    /// `r` at a pseudoasymmetric center.
    LowerR,
    /// `s` at a pseudoasymmetric center.
    LowerS,
    /// `seqTrans` at a pseudoasymmetric double bond.
    SeqTrans,
    /// `seqCis` at a pseudoasymmetric double bond.
    SeqCis,
    /// `E` at a stereogenic double bond.
    E,
    /// `Z` at a stereogenic double bond.
    Z,
    /// `M` at a stereogenic axis.
    M,
    /// `P` at a stereogenic axis.
    P,
    /// `m` at a pseudoasymmetric axis.
    LowerM,
    /// `p` at a pseudoasymmetric axis.
    LowerP,
}

/// An enhanced stereo group: stereo elements whose specified configurations
/// are related as `kind` states.
///
/// A valid group has at least one member and no repeated member, and each
/// element belongs to at most one group.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    /// How the members' configurations are related.
    pub kind: StereoGroupKind,
    /// The member stereo elements; order carries no chemical meaning.
    pub members: Vec<StereoElementId>,
}

/// Relationship between the specified configurations of group members.
///
/// Group identifiers and member order are not chemical information. For a
/// non-absolute group the reference configuration may be inverted jointly at
/// every member without changing the relationship; inverting only some members
/// generally changes the represented diastereomer relationship. Separate groups
/// make independent choices. The kinds remain distinct even when they admit the
/// same pair of configurations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum StereoGroupKind {
    /// Every member has the specified absolute configuration.
    Absolute,
    /// Relative configuration only, without asserting purity or composition.
    Relative,
    /// An equal mixture of the represented configuration and its enantiomer.
    /// Standard CX AND groups do not specify this quantitative composition.
    Racemic,
    /// A mixture containing the correlated configurations; proportions unspecified.
    And,
    /// One of the correlated configurations, without specifying which one.
    Or,
}
