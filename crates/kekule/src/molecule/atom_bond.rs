use super::*;

/// The represented chemistry of one atom.
///
/// Construct an atom with [`Atom::new`] and the `with_*` builders; change an
/// atom in a molecule through [`AtomMut`](super::AtomMut).
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub(crate) element: Element,
    pub(crate) isotope: Option<NonZeroU16>,
    pub(crate) formal_charge: i8,
    pub(crate) radical: Option<AtomRadical>,
    pub(crate) hydrogens: ImplicitHydrogens,
    pub(crate) atom_map: Option<NonZeroU32>,
}

impl Atom {
    /// A neutral atom of `element` at natural isotopic abundance, without
    /// radical electrons or atom map, whose implicit hydrogens are inferred.
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

    /// Returns a copy with this element.
    #[must_use]
    pub fn with_element(mut self, element: Element) -> Self {
        self.element = element;
        self
    }

    /// Returns a copy with this mass number, or natural abundance for `None`.
    #[must_use]
    pub fn with_isotope(mut self, isotope: Option<NonZeroU16>) -> Self {
        self.isotope = isotope;
        self
    }

    /// Returns a copy with this formal charge.
    #[must_use]
    pub fn with_formal_charge(mut self, formal_charge: i8) -> Self {
        self.formal_charge = formal_charge;
        self
    }

    /// Returns a copy with these radical electrons, or none for `None`.
    #[must_use]
    pub fn with_radical(mut self, radical: Option<AtomRadical>) -> Self {
        self.radical = radical;
        self
    }

    /// Returns a copy with this implicit-hydrogen declaration.
    #[must_use]
    pub fn with_hydrogens(mut self, hydrogens: ImplicitHydrogens) -> Self {
        self.hydrogens = hydrogens;
        self
    }

    /// Returns a copy with this reaction atom map, or none for `None`.
    #[must_use]
    pub fn with_atom_map(mut self, atom_map: Option<NonZeroU32>) -> Self {
        self.atom_map = atom_map;
        self
    }

    /// The element.
    pub const fn element(&self) -> Element {
        self.element
    }

    /// The mass number, or `None` for natural isotopic abundance.
    pub const fn isotope(&self) -> Option<NonZeroU16> {
        self.isotope
    }

    /// The formal charge.
    pub const fn formal_charge(&self) -> i8 {
        self.formal_charge
    }

    /// The radical electrons, or `None` without radical occupancy.
    pub const fn radical(&self) -> Option<AtomRadical> {
        self.radical
    }

    /// How the implicit-hydrogen count is determined.
    pub const fn hydrogens(&self) -> ImplicitHydrogens {
        self.hydrogens
    }

    /// The reaction atom map, or `None` for an unmapped atom.
    pub const fn atom_map(&self) -> Option<NonZeroU32> {
        self.atom_map
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
#[non_exhaustive]
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
    /// occupancy has no `AtomRadical`; an unspecified spin uses
    /// `AtomRadical::new(electrons, None)`. A single electron is always a
    /// doublet, so `new(1, None)` records multiplicity 2.
    pub const fn new(electron_count: u8, spin_multiplicity: Option<u8>) -> Option<Self> {
        if electron_count == 0 {
            return None;
        }
        let spin_multiplicity = match (electron_count, spin_multiplicity) {
            (1, None) => Some(2),
            (_, multiplicity) => multiplicity,
        };
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
    /// A single electron always reports a doublet.
    pub const fn spin_multiplicity(self) -> Option<u8> {
        self.spin_multiplicity
    }

    /// Whether the multiplicity says more than the electron count does: a
    /// single electron's doublet is implied, any other assertion is not.
    pub(crate) const fn has_explicit_spin(self) -> bool {
        self.spin_multiplicity.is_some() && self.electron_count != 1
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
#[non_exhaustive]
pub enum BondOrder {
    Zero,
    Single,
    Double,
    Triple,
    Quadruple,
    Dative,
}
