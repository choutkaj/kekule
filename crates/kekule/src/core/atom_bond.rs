use super::*;

#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub element: Element,
    pub isotope: Option<u16>,
    pub formal_charge: i8,
    pub radical: Option<AtomRadical>,
    pub hydrogens: HydrogenDeclaration,
    pub atom_map: Option<u32>,
}

impl Atom {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            isotope: None,
            formal_charge: 0,
            radical: None,
            hydrogens: HydrogenDeclaration::default(),
            atom_map: None,
        }
    }
}

/// The complete non-graph hydrogen statement represented on an atom.
///
/// Graph hydrogen atoms are separate atoms and are not counted here. Hydrogens
/// inferred by valence perception are stored in [`Perception`] rather
/// than this declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HydrogenDeclaration {
    /// The represented count is present and valence perception may infer
    /// additional implicit hydrogens.
    Infer { explicit: u8 },
    /// Exactly this many non-graph hydrogens are represented.
    Fixed(u8),
}

impl HydrogenDeclaration {
    pub const fn explicit_count(self) -> u8 {
        match self {
            Self::Infer { explicit } | Self::Fixed(explicit) => explicit,
        }
    }

    pub const fn allows_implicit(self) -> bool {
        matches!(self, Self::Infer { .. })
    }

    /// Returns the same inference policy with a different represented count.
    pub const fn with_explicit_count(self, explicit: u8) -> Self {
        match self {
            Self::Infer { .. } => Self::Infer { explicit },
            Self::Fixed(_) => Self::Fixed(explicit),
        }
    }
}

impl Default for HydrogenDeclaration {
    fn default() -> Self {
        Self::Infer { explicit: 0 }
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
