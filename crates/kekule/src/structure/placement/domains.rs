//! Electron-domain perception: ideal bond lengths, atom shapes, and the domain
//! counts that placement reads.

use super::lengths;
use crate::molecule::{BondOrder, Element};
use crate::units::{Quantity, ANGSTROM};

/// The ideal length of a bond of `order` between atoms of elements `a` and
/// `b`: a typical length for common pairs (C–H 1.09 Å, C=O 1.22 Å, …), and
/// otherwise the sum of the elements' Cordero covalent radii, shortened for
/// multiple bonds by Pauling's relation `0.71 Å · log10(order)`. Zero-order
/// and dative bonds take the single-bond length.
///
/// Returns `None` when an element has no covalent radius.
///
/// ```
/// use kekule::{molecule::{BondOrder, Element}, structure::placement::ideal_bond_length};
/// use kekule::units::ANGSTROM;
/// let carbon = Element::from_symbol("C").unwrap();
/// let length = ideal_bond_length(carbon, carbon, BondOrder::Double).unwrap();
/// assert_eq!(length.into_unit(ANGSTROM)?.into_value(), 1.34);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn ideal_bond_length(a: Element, b: Element, order: BondOrder) -> Option<Quantity<f64>> {
    lengths::ideal_length_angstrom(a, b, order).map(|length| Quantity::new(length, ANGSTROM))
}

/// The arrangement of an atom's electron domains.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum AtomShape {
    /// One domain, as of a hydrogen.
    Terminal,
    /// Two domains at 180°.
    Linear,
    /// Three coplanar domains at 120°.
    TrigonalPlanar,
    /// Four domains at 109.47°.
    Tetrahedral,
    /// Five domains: two axial, three equatorial.
    TrigonalBipyramidal,
    /// Six domains at 90°.
    Octahedral,
}

impl AtomShape {
    /// The number of domains the shape arranges.
    pub const fn domains(self) -> usize {
        match self {
            Self::Terminal => 1,
            Self::Linear => 2,
            Self::TrigonalPlanar => 3,
            Self::Tetrahedral => 4,
            Self::TrigonalBipyramidal => 5,
            Self::Octahedral => 6,
        }
    }

    pub(crate) const fn from_domains(domains: usize) -> Option<Self> {
        Some(match domains {
            1 => Self::Terminal,
            2 => Self::Linear,
            3 => Self::TrigonalPlanar,
            4 => Self::Tetrahedral,
            5 => Self::TrigonalBipyramidal,
            6 => Self::Octahedral,
            _ => return None,
        })
    }
}

/// An atom's perceived ideal geometry; see the [module](crate::structure::placement) for the rules.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomDomains {
    pub(crate) shape: Option<AtomShape>,
    pub(crate) neighbors: usize,
    pub(crate) implicit_hydrogens: usize,
    pub(crate) lone_pairs: usize,
}

impl AtomDomains {
    /// The arrangement of the atom's domains, or `None` for an atom with none
    /// (an isolated atom without hydrogens or lone pairs).
    pub const fn shape(&self) -> Option<AtomShape> {
        self.shape
    }
    /// Graph neighbors, explicit hydrogens included.
    pub const fn neighbor_count(&self) -> usize {
        self.neighbors
    }
    /// Implicit hydrogens, which take free sites when made explicit.
    pub const fn implicit_hydrogens(&self) -> usize {
        self.implicit_hydrogens
    }
    /// Lone pairs that take sites of the shape; a conjugated nitrogen's lone
    /// pair takes none. Two radical electrons take one site, so a carbene,
    /// singlet (H–C–H about 102°) or triplet (about 134°), places bent at
    /// the trigonal 120°.
    pub const fn lone_pairs(&self) -> usize {
        self.lone_pairs
    }
    /// Sites left for new substituents: the shape's domains less neighbors
    /// and lone pairs. Implicit hydrogens count as free.
    pub fn free_site_count(&self) -> usize {
        self.shape.map_or(0, |shape| {
            shape
                .domains()
                .saturating_sub(self.neighbors + self.lone_pairs)
        })
    }
}

/// The electron-domain count of an atom from its chemistry.
#[derive(Debug, Clone, Copy)]
pub(crate) struct DomainInput {
    pub(crate) element: Element,
    pub(crate) formal_charge: i8,
    pub(crate) radical_electrons: usize,
    /// Neighbors (any bond order) plus implicit hydrogens.
    pub(crate) substituents: usize,
    /// Bond orders plus implicit hydrogens.
    pub(crate) bond_valence: usize,
    /// Zero-order and dative bonds, which may take a lone pair.
    pub(crate) coordinate_bonds: usize,
    /// A three-coordinate nitrogen bonded to an atom with a multiple bond.
    pub(crate) conjugated_nitrogen: bool,
}

/// The domain count, which may be zero or exceed six, and the lone pairs
/// that take sites among them.
pub(crate) fn perceive_domains(input: DomainInput) -> (usize, usize) {
    let lone_pairs = valence_electrons(input.element.atomic_number()).map_or(0, |valence| {
        let free = i64::from(valence)
            - i64::from(input.formal_charge)
            - i64::try_from(input.bond_valence).unwrap_or(i64::MAX / 4);
        let free = usize::try_from(free.max(0)).expect("non-negative");
        let radical = input.radical_electrons.min(free);
        // Each two radical electrons take a site, as a carbene's do: a
        // singlet's lone pair, or a triplet's in-plane electron (its other
        // one is in the p orbital). A single one, as of a methyl radical,
        // takes none.
        let pairs = (free - radical) / 2 + radical / 2;
        // A coordinate bond donates one of the atom's own pairs.
        pairs.saturating_sub(input.coordinate_bonds)
    });
    if input.conjugated_nitrogen && lone_pairs == 1 && input.substituents == 3 {
        return (3, 0);
    }
    (input.substituents + lone_pairs, lone_pairs)
}

/// Whether placement perceives the lone pairs of `element`, an s- or p-block
/// element, so its shape follows its electrons rather than its neighbors.
pub(crate) const fn models_lone_pairs(element: Element) -> bool {
    valence_electrons(element.atomic_number()).is_some()
}

/// Valence electrons of s- and p-block elements; `None` elsewhere.
const fn valence_electrons(atomic_number: u8) -> Option<u8> {
    Some(match atomic_number {
        1 | 2 => atomic_number,
        3..=10 => atomic_number - 2,
        11..=18 => atomic_number - 10,
        19 | 37 | 55 | 87 => 1,
        20 | 38 | 56 | 88 => 2,
        31..=36 => atomic_number - 28,
        49..=54 => atomic_number - 46,
        81..=86 => atomic_number - 78,
        113..=118 => atomic_number - 110,
        _ => return None,
    })
}

#[cfg(test)]
mod tests;
