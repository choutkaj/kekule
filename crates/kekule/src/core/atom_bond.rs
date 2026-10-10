use super::*;

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub element: Element,
    pub isotope: Option<u16>,
    pub formal_charge: i8,
    pub radical: Option<AtomRadical>,
    pub hydrogens: ImplicitHydrogens,
    pub atom_map: Option<u32>,
}

impl Atom {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            isotope: None,
            formal_charge: 0,
            radical: None,
            hydrogens: ImplicitHydrogens::Inferred,
            atom_map: None,
        }
    }
}

/// Specifies how an atom's implicit (non-graph) hydrogen count is determined.
///
/// Explicit hydrogens are separate graph atoms and are never counted here.
/// Use [`Molecule::implicit_hydrogens`] for the resolved implicit count and
/// [`Molecule::total_hydrogens`] to include explicit hydrogen neighbors.
///
/// An atom's implicit count is either wholly fixed or wholly inferred. Chemical
/// edits preserve a fixed count and invalidate an inferred one; reperception
/// never overwrites a fixed count. In particular, SMILES `[C]` fixes the count
/// at zero, `[CH3]` fixes it at three, and `C` infers it.
///
/// [`Molecule::implicit_hydrogens`]: super::Molecule::implicit_hydrogens
/// [`Molecule::total_hydrogens`]: super::Molecule::total_hydrogens
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ImplicitHydrogens {
    /// The valence model supplies the complete implicit count. It is unknown
    /// until valence perception is installed.
    #[default]
    Inferred,
    /// Exactly this many non-graph hydrogens, known without perception and
    /// never changed by valence inference.
    Fixed(u8),
}

impl ImplicitHydrogens {
    /// Returns the fixed count, or `None` when the count is inferred.
    pub const fn fixed_count(self) -> Option<u8> {
        match self {
            Self::Inferred => None,
            Self::Fixed(count) => Some(count),
        }
    }

    /// Returns the count represented without perception: the fixed count, or
    /// zero for an inferred count. Inferred hydrogens are not included.
    pub const fn represented_count(self) -> u8 {
        match self {
            Self::Inferred => 0,
            Self::Fixed(count) => count,
        }
    }

    pub const fn is_inferred(self) -> bool {
        matches!(self, Self::Inferred)
    }

    /// Resolves the implicit count from an installed inferred count.
    /// Fixed counts are known without perception; an inferred count is
    /// unknown until perception supplies it.
    pub fn resolve(self, inferred: Option<u8>) -> Option<usize> {
        match self {
            Self::Fixed(count) => Some(usize::from(count)),
            Self::Inferred => inferred.map(usize::from),
        }
    }
}

/// Nonbonding radical-electron occupancy with an optional local spin assertion.
///
/// Electron count and spin multiplicity are distinct. In particular, the
/// two-electron singlet and triplet states encoded by molfiles both reserve two
/// electrons during valence inference. A count inferred from bracket SMILES
/// does not by itself assert the spin state. This is atom-local information,
/// not a molecular spin multiplicity or a prediction of the electronic ground state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AtomRadical {
    electron_count: u8,
    spin_multiplicity: Option<u8>,
}

impl AtomRadical {
    /// Constructs a nonzero radical-electron count and optional multiplicity `2S+1`.
    ///
    /// Returns `None` for zero electrons or a multiplicity incompatible with
    /// coupling that many spin-one-half electrons. An atom without radical
    /// occupancy uses `Atom::radical = None`; an unspecified spin uses
    /// `AtomRadical::new(electrons, None)`.
    pub const fn new(electron_count: u8, spin_multiplicity: Option<u8>) -> Option<Self> {
        if electron_count == 0 {
            return None;
        }
        if let Some(multiplicity) = spin_multiplicity {
            if multiplicity == 0
                || multiplicity as u16 > electron_count as u16 + 1
                || multiplicity % 2 == electron_count % 2
            {
                return None;
            }
        }
        Some(Self {
            electron_count,
            spin_multiplicity,
        })
    }

    /// Electrons reserved from bonding, including a pair in a singlet center.
    /// This is not the number of physically unpaired electrons.
    pub const fn electron_count(self) -> u8 {
        self.electron_count
    }

    /// Explicit atom-local spin multiplicity, or `None` when not specified.
    pub const fn spin_multiplicity(self) -> Option<u8> {
        self.spin_multiplicity
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    pub(crate) a: AtomId,
    pub(crate) b: AtomId,
    pub order: BondOrder,
}

impl Bond {
    pub fn new(a: AtomId, b: AtomId, order: BondOrder) -> Self {
        Self { a, b, order }
    }

    pub const fn a(&self) -> AtomId {
        self.a
    }

    pub const fn b(&self) -> AtomId {
        self.b
    }

    pub const fn endpoints(&self) -> (AtomId, AtomId) {
        (self.a, self.b)
    }
}

/// A canonical localized bond order stored in a [`Molecule`](super::Molecule).
///
/// Aromaticity is perceived state, not a represented bond order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Zero,
    Single,
    Double,
    Triple,
    Quadruple,
    Dative,
}
