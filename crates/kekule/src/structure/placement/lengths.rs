//! Ideal bond lengths by element pair and bond order.
use crate::molecule::{BondOrder, Element};

/// Typical lengths in Å of common bonds in organic and biological molecules,
/// rounded from the mean values of Allen et al. (1987) and standard
/// compilations; pairs are ordered by atomic number.
const TABLE: &[(u8, u8, BondOrder, f64)] = {
    use BondOrder::{Double, Single, Triple};
    &[
        (1, 1, Single, 0.74),
        (1, 5, Single, 1.19),
        (1, 6, Single, 1.09),
        (1, 7, Single, 1.01),
        (1, 8, Single, 0.96),
        (1, 9, Single, 0.92),
        (1, 14, Single, 1.48),
        (1, 15, Single, 1.42),
        (1, 16, Single, 1.34),
        (1, 17, Single, 1.27),
        (1, 34, Single, 1.46),
        (1, 35, Single, 1.41),
        (1, 53, Single, 1.61),
        (5, 6, Single, 1.56),
        (5, 8, Single, 1.36),
        (6, 6, Single, 1.54),
        (6, 6, Double, 1.34),
        (6, 6, Triple, 1.20),
        (6, 7, Single, 1.47),
        (6, 7, Double, 1.28),
        (6, 7, Triple, 1.16),
        (6, 8, Single, 1.43),
        (6, 8, Double, 1.22),
        (6, 8, Triple, 1.13),
        (6, 9, Single, 1.35),
        (6, 14, Single, 1.87),
        (6, 15, Single, 1.85),
        (6, 15, Double, 1.67),
        (6, 16, Single, 1.82),
        (6, 16, Double, 1.67),
        (6, 17, Single, 1.77),
        (6, 34, Single, 1.95),
        (6, 35, Single, 1.94),
        (6, 53, Single, 2.14),
        (7, 7, Single, 1.45),
        (7, 7, Double, 1.25),
        (7, 7, Triple, 1.10),
        (7, 8, Single, 1.40),
        (7, 8, Double, 1.21),
        (7, 15, Single, 1.70),
        (7, 16, Single, 1.63),
        (7, 16, Double, 1.54),
        (8, 8, Single, 1.48),
        (8, 8, Double, 1.21),
        (8, 14, Single, 1.63),
        (8, 15, Single, 1.61),
        (8, 15, Double, 1.48),
        (8, 16, Single, 1.58),
        (8, 16, Double, 1.43),
        (8, 17, Single, 1.42),
        (8, 17, Double, 1.43),
        (15, 16, Double, 1.95),
        (16, 16, Single, 2.05),
    ]
};

/// Pauling's empirical shortening per decade of bond order, in Å.
const PAULING_SHORTENING: f64 = 0.71;

/// The ideal length in Å of a bond of `order` between `a` and `b`: the
/// tabulated typical length, or the sum of Cordero covalent radii shortened
/// by Pauling's relation for multiple bonds. Zero-order and dative bonds take
/// the single-bond length. `None` when an element has no covalent radius.
pub(crate) fn ideal_length_angstrom(a: Element, b: Element, order: BondOrder) -> Option<f64> {
    let (low, high) = {
        let (x, y) = (a.atomic_number(), b.atomic_number());
        (x.min(y), x.max(y))
    };
    let order = match order {
        BondOrder::Zero | BondOrder::Dative => BondOrder::Single,
        other => other,
    };
    if let Some((_, _, _, length)) = TABLE
        .iter()
        .find(|(x, y, tabulated, _)| *x == low && *y == high && *tabulated == order)
    {
        return Some(*length);
    }
    let radius = |element: Element| {
        element
            .covalent_radius()?
            .into_unit(crate::units::ANGSTROM)
            .ok()
            .map(crate::units::Quantity::into_value)
    };
    let single = radius(a)? + radius(b)?;
    let multiplicity: f64 = match order {
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        _ => 1.0,
    };
    Some(single - PAULING_SHORTENING * multiplicity.log10())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn element(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    #[test]
    fn table_is_ordered_and_symmetric_lookup_works() {
        for (low, high, _, length) in TABLE {
            assert!(low <= high && *length > 0.5 && *length < 2.5);
        }
        assert_eq!(
            ideal_length_angstrom(element("O"), element("C"), BondOrder::Double),
            Some(1.22)
        );
        assert_eq!(
            ideal_length_angstrom(element("H"), element("C"), BondOrder::Dative),
            Some(1.09)
        );
    }

    #[test]
    fn fallback_shortens_multiple_bonds() {
        let single =
            ideal_length_angstrom(element("Si"), element("Si"), BondOrder::Single).unwrap();
        let double =
            ideal_length_angstrom(element("Si"), element("Si"), BondOrder::Double).unwrap();
        let triple =
            ideal_length_angstrom(element("Si"), element("Si"), BondOrder::Triple).unwrap();
        assert!((single - 2.22).abs() < 1e-9);
        assert!(single > double && double > triple);
    }
}
