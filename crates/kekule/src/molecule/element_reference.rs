use super::Element;
use crate::units::{Quantity, ANGSTROM, BOHR, CUBIC_BOHR};

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;

/// Selects one of the van der Waals radius reference columns.
///
/// The columns come from different definitions and are not interchangeable:
/// the free-atom radii of Charry and Tkatchenko are larger than the
/// crystallographic contact radii of Bondi (carbon 1.91 Å against 1.70 Å),
/// and Alvarez's contact-distribution radii are larger again for most
/// metals. Only the two Charry–Tkatchenko columns fall back to each other;
/// [`Bondi`](Self::Bondi) and [`Alvarez`](Self::Alvarez) never borrow values
/// from another column, so a caller that wants one scale completed by another
/// composes the lookups explicitly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[non_exhaustive]
pub enum VanDerWaalsRadiusSource {
    /// Free-atom radii reported by Charry and Tkatchenko, J. Chem. Theory
    /// Comput. 2024, DOI 10.1021/acs.jctc.4c00784, tabulated in bohr.
    ///
    /// Elements 1–118 except livermorium (116); an element without a
    /// free-atom value takes the [`Reference`](Self::Reference) value.
    #[default]
    FreeAtom,
    /// The compiled bonded/reference-radius column reported in the same work,
    /// tabulated in bohr.
    ///
    /// An element without a reference value, such as the lanthanides from
    /// cerium and every element after protactinium, takes the
    /// [`FreeAtom`](Self::FreeAtom) value.
    Reference,
    /// Bondi's crystallographic van der Waals radii, J. Phys. Chem. 1964, 68,
    /// 441, DOI 10.1021/j100785a001, with the usual completions, tabulated in
    /// ångström. These are the radii molecular graphics conventionally
    /// display.
    ///
    /// The 44 main-group elements, hydrogen to radium, are the consistent set
    /// of Mantina, Chamberlin, Valero, Cramer and Truhlar, J. Phys. Chem. A
    /// 2009, 113, 5806, DOI 10.1021/jp8111556, Table 12, transcribed
    /// verbatim: Bondi's values for 27 elements, hydrogen 1.10 Å from Rowland
    /// and Taylor, J. Phys. Chem. 1996, 100, 7384, DOI 10.1021/jp953141+,
    /// in place of Bondi's 1.20 Å, and Mantina's own values for the 16
    /// elements Bondi did not cover (Be, B, Al, Ca, Ge, Rb, Sr, Sb, Cs, Ba,
    /// Bi, Po, At, Rn, Fr, Ra). Strontium is 2.49 Å as in Table 12; the
    /// paper's abstract states 2.50 Å. The ten further elements Bondi gave
    /// (Ni, Cu, Zn, Pd, Ag, Cd, Pt, Au, Hg, U) are transcribed from the Bondi
    /// column of Alvarez 2013, Table 1, which draws on Bondi 1964 and on
    /// Bondi's 1966 values for metals in covalent compounds, J. Phys. Chem.
    /// 1966, 70, 3006, DOI 10.1021/j100881a503.
    ///
    /// Covers 54 elements; every other element, including the remaining
    /// transition metals and lanthanides, has no Bondi radius.
    Bondi,
    /// Alvarez's van der Waals radii from the distribution of intermolecular
    /// contact distances in crystal structures, Dalton Trans. 2013, 42,
    /// 8617, DOI 10.1039/C3DT50599E, Table 1 (the `rvdW` column), tabulated
    /// in ångström.
    ///
    /// Covers the 93 elements hydrogen to neodymium, samarium to bismuth,
    /// and actinium to einsteinium; promethium, polonium to radium, and
    /// elements after einsteinium have no Alvarez radius. Values the paper
    /// brackets as rough estimates from very small datasets (He, Ne, Ac,
    /// Pa, Cm, Bk, Cf, Es) are included as printed. Later revisions, such
    /// as Vogt and Alvarez's noble-gas radii, are not part of this column.
    Alvarez,
}

impl Element {
    /// Returns the static dipole polarizability, in cubic bohr as tabulated.
    ///
    /// Dipole polarizabilities are the values compiled by Charry and
    /// Tkatchenko, J. Chem. Theory Comput. 2024, DOI 10.1021/acs.jctc.4c00784,
    /// from Schwerdtfeger and Nagle, Mol. Phys. 2019,
    /// DOI 10.1080/00268976.2018.1535143.
    pub fn dipole_polarizability(self) -> Option<Quantity<f64>> {
        DIPOLE_POLARIZABILITY_BOHR3[usize::from(self.atomic_number())]
            .map(|value| Quantity::new(value, CUBIC_BOHR))
    }

    /// Returns a van der Waals radius from the selected column, in the unit
    /// the column is tabulated in: bohr for the Charry–Tkatchenko columns,
    /// ångström for Bondi and Alvarez.
    ///
    /// If the requested Charry–Tkatchenko column is unavailable for this
    /// element, the other Charry–Tkatchenko column is used, and `None` is
    /// returned only when Kekule has neither. The Bondi and Alvarez columns
    /// have no fallback and return `None` outside their coverage. Each
    /// [`VanDerWaalsRadiusSource`] variant documents its source and
    /// coverage.
    ///
    /// ```
    /// use kekule::molecule::{Element, VanDerWaalsRadiusSource};
    /// use kekule::units::ANGSTROM;
    ///
    /// let carbon = Element::C.van_der_waals_radius(VanDerWaalsRadiusSource::Bondi);
    /// assert_eq!(carbon.unwrap().value_in(ANGSTROM)?, 1.70);
    /// let iron = Element::from_symbol("Fe").unwrap();
    /// assert!(iron.van_der_waals_radius(VanDerWaalsRadiusSource::Bondi).is_none());
    /// let alvarez = iron.van_der_waals_radius(VanDerWaalsRadiusSource::Alvarez);
    /// assert_eq!(alvarez.unwrap().value_in(ANGSTROM)?, 2.44);
    /// # Ok::<(), kekule::units::UnitError>(())
    /// ```
    pub fn van_der_waals_radius(self, source: VanDerWaalsRadiusSource) -> Option<Quantity<f64>> {
        let index = usize::from(self.atomic_number());
        match source {
            VanDerWaalsRadiusSource::FreeAtom | VanDerWaalsRadiusSource::Reference => self
                .charry_tkatchenko_radius_bohr(source)
                .map(|radius| Quantity::new(radius, BOHR)),
            VanDerWaalsRadiusSource::Bondi => {
                BONDI_VDW_RADIUS_ANGSTROM[index].map(|radius| Quantity::new(radius, ANGSTROM))
            }
            VanDerWaalsRadiusSource::Alvarez => {
                ALVAREZ_VDW_RADIUS_ANGSTROM[index].map(|radius| Quantity::new(radius, ANGSTROM))
            }
        }
    }

    fn charry_tkatchenko_radius_bohr(self, source: VanDerWaalsRadiusSource) -> Option<f64> {
        let index = usize::from(self.atomic_number());
        match source {
            VanDerWaalsRadiusSource::Reference => {
                REFERENCE_VDW_RADIUS_BOHR[index].or(TKATCHENKO_VDW_RADIUS_BOHR[index])
            }
            _ => TKATCHENKO_VDW_RADIUS_BOHR[index].or(REFERENCE_VDW_RADIUS_BOHR[index]),
        }
    }

    /// The radius in ångström. Charry–Tkatchenko radii are converted from
    /// bohr as solvation has always converted them.
    pub(crate) fn van_der_waals_radius_angstrom(
        self,
        source: VanDerWaalsRadiusSource,
    ) -> Option<f64> {
        let index = usize::from(self.atomic_number());
        match source {
            VanDerWaalsRadiusSource::FreeAtom | VanDerWaalsRadiusSource::Reference => self
                .charry_tkatchenko_radius_bohr(source)
                .map(|radius| radius * BOHR_TO_ANGSTROM),
            VanDerWaalsRadiusSource::Bondi => BONDI_VDW_RADIUS_ANGSTROM[index],
            VanDerWaalsRadiusSource::Alvarez => ALVAREZ_VDW_RADIUS_ANGSTROM[index],
        }
    }
}

const DIPOLE_POLARIZABILITY_BOHR3: [Option<f64>; 119] = [
    None,        // 0: X
    Some(4.5),   // 1: H
    Some(1.38),  // 2: He
    Some(164.0), // 3: Li
    Some(37.7),  // 4: Be
    Some(20.5),  // 5: B
    Some(11.3),  // 6: C
    Some(7.4),   // 7: N
    Some(5.3),   // 8: O
    Some(3.74),  // 9: F
    Some(2.66),  // 10: Ne
    Some(163.0), // 11: Na
    Some(71.2),  // 12: Mg
    Some(57.8),  // 13: Al
    Some(37.3),  // 14: Si
    Some(25.0),  // 15: P
    Some(19.4),  // 16: S
    Some(14.6),  // 17: Cl
    Some(11.1),  // 18: Ar
    Some(290.0), // 19: K
    Some(161.0), // 20: Ca
    Some(97.0),  // 21: Sc
    Some(100.0), // 22: Ti
    Some(87.0),  // 23: V
    Some(83.0),  // 24: Cr
    Some(68.0),  // 25: Mn
    Some(62.0),  // 26: Fe
    Some(55.0),  // 27: Co
    Some(49.0),  // 28: Ni
    Some(47.0),  // 29: Cu
    Some(38.7),  // 30: Zn
    Some(50.0),  // 31: Ga
    Some(40.0),  // 32: Ge
    Some(30.0),  // 33: As
    Some(29.0),  // 34: Se
    Some(21.0),  // 35: Br
    Some(16.8),  // 36: Kr
    Some(320.0), // 37: Rb
    Some(197.0), // 38: Sr
    Some(162.0), // 39: Y
    Some(112.0), // 40: Zr
    Some(98.0),  // 41: Nb
    Some(87.0),  // 42: Mo
    Some(79.0),  // 43: Tc
    Some(72.0),  // 44: Ru
    Some(66.0),  // 45: Rh
    Some(26.1),  // 46: Pd
    Some(55.0),  // 47: Ag
    Some(46.0),  // 48: Cd
    Some(65.0),  // 49: In
    Some(53.0),  // 50: Sn
    Some(43.0),  // 51: Sb
    Some(38.0),  // 52: Te
    Some(32.9),  // 53: I
    Some(27.3),  // 54: Xe
    Some(401.0), // 55: Cs
    Some(272.0), // 56: Ba
    Some(215.0), // 57: La
    Some(205.0), // 58: Ce
    Some(216.0), // 59: Pr
    Some(208.0), // 60: Nd
    Some(200.0), // 61: Pm
    Some(192.0), // 62: Sm
    Some(184.0), // 63: Eu
    Some(158.0), // 64: Gd
    Some(170.0), // 65: Tb
    Some(165.0), // 66: Dy
    Some(156.0), // 67: Ho
    Some(150.0), // 68: Er
    Some(144.0), // 69: Tm
    Some(139.0), // 70: Yb
    Some(137.0), // 71: Lu
    Some(103.0), // 72: Hf
    Some(74.0),  // 73: Ta
    Some(68.0),  // 74: W
    Some(62.0),  // 75: Re
    Some(57.0),  // 76: Os
    Some(54.0),  // 77: Ir
    Some(48.0),  // 78: Pt
    Some(36.0),  // 79: Au
    Some(33.9),  // 80: Hg
    Some(50.0),  // 81: Tl
    Some(47.0),  // 82: Pb
    Some(48.0),  // 83: Bi
    Some(44.0),  // 84: Po
    Some(42.0),  // 85: At
    Some(35.0),  // 86: Rn
    Some(318.0), // 87: Fr
    Some(246.0), // 88: Ra
    Some(203.0), // 89: Ac
    Some(217.0), // 90: Th
    Some(154.0), // 91: Pa
    Some(129.0), // 92: U
    Some(151.0), // 93: Np
    Some(132.0), // 94: Pu
    Some(131.0), // 95: Am
    Some(144.0), // 96: Cm
    Some(125.0), // 97: Bk
    Some(122.0), // 98: Cf
    Some(118.0), // 99: Es
    Some(113.0), // 100: Fm
    Some(109.0), // 101: Md
    Some(110.0), // 102: No
    Some(320.0), // 103: Lr
    Some(112.0), // 104: Rf
    Some(42.0),  // 105: Db
    Some(40.0),  // 106: Sg
    Some(38.0),  // 107: Bh
    Some(36.0),  // 108: Hs
    Some(34.0),  // 109: Mt
    Some(32.0),  // 110: Ds
    Some(32.0),  // 111: Rg
    Some(28.0),  // 112: Cn
    Some(29.0),  // 113: Nh
    Some(31.0),  // 114: Fl
    Some(71.0),  // 115: Mc
    None,        // 116: Lv
    Some(76.0),  // 117: Ts
    Some(58.0),  // 118: Og
];

const TKATCHENKO_VDW_RADIUS_BOHR: [Option<f64>; 119] = [
    None, // 0: no element
    Some(3.164697),
    Some(2.672999),
    Some(5.289595),
    Some(4.2875),
    Some(3.9302),
    Some(3.6096),
    Some(3.398),
    Some(3.24),
    Some(3.0822),
    Some(2.935712),
    Some(5.285),
    Some(4.6952),
    Some(4.5574),
    Some(4.281),
    Some(4.043),
    Some(3.8993),
    Some(3.7441),
    Some(3.600377),
    Some(5.7384),
    Some(5.276),
    Some(4.907),
    Some(4.929),
    Some(4.832),
    Some(4.799),
    Some(4.664),
    Some(4.603),
    Some(4.525),
    Some(4.451),
    Some(4.425),
    Some(4.3036),
    Some(4.464),
    Some(4.324),
    Some(4.15),
    Some(4.13),
    Some(3.944),
    Some(3.819973),
    Some(5.8196),
    Some(5.43),
    Some(5.28),
    Some(5.009),
    Some(4.914),
    Some(4.832),
    Some(4.765),
    Some(4.703),
    Some(4.64),
    Some(4.0681),
    Some(4.525),
    Some(4.411),
    Some(4.635),
    Some(4.501),
    Some(4.369),
    Some(4.292),
    Some(4.2049),
    Some(4.0943),
    Some(6.0103),
    Some(5.686),
    Some(5.498),
    Some(5.461),
    Some(5.502),
    Some(5.472),
    Some(5.442),
    Some(5.41),
    Some(5.377),
    Some(5.262),
    Some(5.317),
    Some(5.294),
    Some(5.252),
    Some(5.223),
    Some(5.192),
    Some(5.166),
    Some(5.155),
    Some(4.95),
    Some(4.72),
    Some(4.66),
    Some(4.603),
    Some(4.548),
    Some(4.513),
    Some(4.438),
    Some(4.259),
    Some(4.223),
    Some(4.464),
    Some(4.425),
    Some(4.438),
    Some(4.383),
    Some(4.354),
    Some(4.242),
    Some(5.8144),
    Some(5.605),
    Some(5.453),
    Some(5.51),
    Some(5.242),
    Some(5.111),
    Some(5.228),
    Some(5.13),
    Some(5.12),
    Some(5.19),
    Some(5.09),
    Some(5.07),
    Some(5.05),
    Some(5.02),
    Some(4.99),
    Some(4.996),
    Some(5.82),
    Some(5.009),
    Some(4.354),
    Some(4.324),
    Some(4.292),
    Some(4.259),
    Some(4.225),
    Some(4.188),
    Some(4.19),
    Some(4.109),
    Some(4.13),
    Some(4.169),
    Some(4.69),
    None,
    Some(4.74),
    Some(4.56),
];

const REFERENCE_VDW_RADIUS_BOHR: [Option<f64>; 119] = [
    None, // 0: no element
    Some(3.16),
    Some(2.65),
    Some(4.97),
    Some(4.2141),
    Some(3.8739),
    Some(3.7039),
    Some(3.3826),
    Some(3.2314),
    Some(3.118),
    Some(2.91),
    Some(5.2345),
    Some(4.5731),
    Some(4.5353),
    Some(4.2708),
    Some(4.044),
    Some(3.8928),
    Some(3.8739),
    Some(3.55),
    Some(5.707),
    Some(5.2534),
    Some(4.9511),
    Some(4.6109),
    Some(4.2897),
    Some(4.2141),
    Some(4.2519),
    Some(4.2897),
    Some(4.2519),
    Some(4.2141),
    Some(4.2897),
    Some(4.233),
    Some(4.5542),
    Some(4.3842),
    Some(4.2519),
    Some(4.1196),
    Some(3.9684),
    Some(3.82),
    Some(5.9526),
    Some(5.5558),
    Some(5.1212),
    Some(4.8566),
    Some(4.6487),
    Some(4.5164),
    Some(4.4787),
    Some(4.4787),
    Some(4.3842),
    Some(4.4409),
    Some(4.4787),
    Some(4.4787),
    Some(4.781),
    Some(4.6487),
    Some(4.5542),
    Some(4.5542),
    Some(4.4598),
    Some(4.08),
    Some(6.2361),
    Some(5.7637),
    Some(5.3101),
    None, // 58: Ce
    None, // 59: Pr
    None, // 60: Nd
    None, // 61: Pm
    None, // 62: Sm
    None, // 63: Eu
    None, // 64: Gd
    None, // 65: Tb
    None, // 66: Dy
    None, // 67: Ho
    None, // 68: Er
    None, // 69: Tm
    None, // 70: Yb
    None, // 71: Lu
    Some(4.7621),
    Some(4.5731),
    Some(4.4598),
    Some(4.4409),
    Some(4.4031),
    Some(4.422),
    Some(4.4787),
    Some(4.5542),
    Some(4.2519),
    Some(4.781),
    Some(4.781),
    Some(4.7621),
    None,
    None,
    Some(4.23),
    None,
    None,
    None,
    Some(5.1967),
    Some(5.0078),
    None, // 92: U
    None, // 93: Np
    None, // 94: Pu
    None, // 95: Am
    None, // 96: Cm
    None, // 97: Bk
    None, // 98: Cf
    None, // 99: Es
    None, // 100: Fm
    None, // 101: Md
    None, // 102: No
    None, // 103: Lr
    None, // 104: Rf
    None, // 105: Db
    None, // 106: Sg
    None, // 107: Bh
    None, // 108: Hs
    None, // 109: Mt
    None, // 110: Ds
    None, // 111: Rg
    None, // 112: Cn
    None, // 113: Nh
    None, // 114: Fl
    None, // 115: Mc
    None, // 116: Lv
    None, // 117: Ts
    None, // 118: Og
];

const BONDI_VDW_RADIUS_ANGSTROM: [Option<f64>; 119] = [
    None,       // 0: no element
    Some(1.10), // 1: H (Rowland & Taylor)
    Some(1.40), // 2: He
    Some(1.81), // 3: Li
    Some(1.53), // 4: Be (Mantina)
    Some(1.92), // 5: B (Mantina)
    Some(1.70), // 6: C
    Some(1.55), // 7: N
    Some(1.52), // 8: O
    Some(1.47), // 9: F
    Some(1.54), // 10: Ne
    Some(2.27), // 11: Na
    Some(1.73), // 12: Mg
    Some(1.84), // 13: Al (Mantina)
    Some(2.10), // 14: Si
    Some(1.80), // 15: P
    Some(1.80), // 16: S
    Some(1.75), // 17: Cl
    Some(1.88), // 18: Ar
    Some(2.75), // 19: K
    Some(2.31), // 20: Ca (Mantina)
    None,       // 21: Sc
    None,       // 22: Ti
    None,       // 23: V
    None,       // 24: Cr
    None,       // 25: Mn
    None,       // 26: Fe
    None,       // 27: Co
    Some(1.63), // 28: Ni
    Some(1.40), // 29: Cu
    Some(1.39), // 30: Zn
    Some(1.87), // 31: Ga
    Some(2.11), // 32: Ge (Mantina)
    Some(1.85), // 33: As
    Some(1.90), // 34: Se
    Some(1.83), // 35: Br
    Some(2.02), // 36: Kr
    Some(3.03), // 37: Rb (Mantina)
    Some(2.49), // 38: Sr (Mantina)
    None,       // 39: Y
    None,       // 40: Zr
    None,       // 41: Nb
    None,       // 42: Mo
    None,       // 43: Tc
    None,       // 44: Ru
    None,       // 45: Rh
    Some(1.63), // 46: Pd
    Some(1.72), // 47: Ag
    Some(1.62), // 48: Cd
    Some(1.93), // 49: In
    Some(2.17), // 50: Sn
    Some(2.06), // 51: Sb (Mantina)
    Some(2.06), // 52: Te
    Some(1.98), // 53: I
    Some(2.16), // 54: Xe
    Some(3.43), // 55: Cs (Mantina)
    Some(2.68), // 56: Ba (Mantina)
    None,       // 57: La
    None,       // 58: Ce
    None,       // 59: Pr
    None,       // 60: Nd
    None,       // 61: Pm
    None,       // 62: Sm
    None,       // 63: Eu
    None,       // 64: Gd
    None,       // 65: Tb
    None,       // 66: Dy
    None,       // 67: Ho
    None,       // 68: Er
    None,       // 69: Tm
    None,       // 70: Yb
    None,       // 71: Lu
    None,       // 72: Hf
    None,       // 73: Ta
    None,       // 74: W
    None,       // 75: Re
    None,       // 76: Os
    None,       // 77: Ir
    Some(1.72), // 78: Pt
    Some(1.66), // 79: Au
    Some(1.70), // 80: Hg
    Some(1.96), // 81: Tl
    Some(2.02), // 82: Pb
    Some(2.07), // 83: Bi (Mantina)
    Some(1.97), // 84: Po (Mantina)
    Some(2.02), // 85: At (Mantina)
    Some(2.20), // 86: Rn (Mantina)
    Some(3.48), // 87: Fr (Mantina)
    Some(2.83), // 88: Ra (Mantina)
    None,       // 89: Ac
    None,       // 90: Th
    None,       // 91: Pa
    Some(1.86), // 92: U
    None,       // 93: Np
    None,       // 94: Pu
    None,       // 95: Am
    None,       // 96: Cm
    None,       // 97: Bk
    None,       // 98: Cf
    None,       // 99: Es
    None,       // 100: Fm
    None,       // 101: Md
    None,       // 102: No
    None,       // 103: Lr
    None,       // 104: Rf
    None,       // 105: Db
    None,       // 106: Sg
    None,       // 107: Bh
    None,       // 108: Hs
    None,       // 109: Mt
    None,       // 110: Ds
    None,       // 111: Rg
    None,       // 112: Cn
    None,       // 113: Nh
    None,       // 114: Fl
    None,       // 115: Mc
    None,       // 116: Lv
    None,       // 117: Ts
    None,       // 118: Og
];

const ALVAREZ_VDW_RADIUS_ANGSTROM: [Option<f64>; 119] = [
    None,       // 0: no element
    Some(1.20), // 1: H
    Some(1.43), // 2: He (rough estimate)
    Some(2.12), // 3: Li
    Some(1.98), // 4: Be
    Some(1.91), // 5: B
    Some(1.77), // 6: C
    Some(1.66), // 7: N
    Some(1.50), // 8: O
    Some(1.46), // 9: F
    Some(1.58), // 10: Ne (rough estimate)
    Some(2.50), // 11: Na
    Some(2.51), // 12: Mg
    Some(2.25), // 13: Al
    Some(2.19), // 14: Si
    Some(1.90), // 15: P
    Some(1.89), // 16: S
    Some(1.82), // 17: Cl
    Some(1.83), // 18: Ar
    Some(2.73), // 19: K
    Some(2.62), // 20: Ca
    Some(2.58), // 21: Sc
    Some(2.46), // 22: Ti
    Some(2.42), // 23: V
    Some(2.45), // 24: Cr
    Some(2.45), // 25: Mn
    Some(2.44), // 26: Fe
    Some(2.40), // 27: Co
    Some(2.40), // 28: Ni
    Some(2.38), // 29: Cu
    Some(2.39), // 30: Zn
    Some(2.32), // 31: Ga
    Some(2.29), // 32: Ge
    Some(1.88), // 33: As
    Some(1.82), // 34: Se
    Some(1.86), // 35: Br
    Some(2.25), // 36: Kr
    Some(3.21), // 37: Rb
    Some(2.84), // 38: Sr
    Some(2.75), // 39: Y
    Some(2.52), // 40: Zr
    Some(2.56), // 41: Nb
    Some(2.45), // 42: Mo
    Some(2.44), // 43: Tc
    Some(2.46), // 44: Ru
    Some(2.44), // 45: Rh
    Some(2.15), // 46: Pd
    Some(2.53), // 47: Ag
    Some(2.49), // 48: Cd
    Some(2.43), // 49: In
    Some(2.42), // 50: Sn
    Some(2.47), // 51: Sb
    Some(1.99), // 52: Te
    Some(2.04), // 53: I
    Some(2.06), // 54: Xe
    Some(3.48), // 55: Cs
    Some(3.03), // 56: Ba
    Some(2.98), // 57: La
    Some(2.88), // 58: Ce
    Some(2.92), // 59: Pr
    Some(2.95), // 60: Nd
    None,       // 61: Pm
    Some(2.90), // 62: Sm
    Some(2.87), // 63: Eu
    Some(2.83), // 64: Gd
    Some(2.79), // 65: Tb
    Some(2.87), // 66: Dy
    Some(2.81), // 67: Ho
    Some(2.83), // 68: Er
    Some(2.79), // 69: Tm
    Some(2.80), // 70: Yb
    Some(2.74), // 71: Lu
    Some(2.63), // 72: Hf
    Some(2.53), // 73: Ta
    Some(2.57), // 74: W
    Some(2.49), // 75: Re
    Some(2.48), // 76: Os
    Some(2.41), // 77: Ir
    Some(2.29), // 78: Pt
    Some(2.32), // 79: Au
    Some(2.45), // 80: Hg
    Some(2.47), // 81: Tl
    Some(2.60), // 82: Pb
    Some(2.54), // 83: Bi
    None,       // 84: Po
    None,       // 85: At
    None,       // 86: Rn
    None,       // 87: Fr
    None,       // 88: Ra
    Some(2.8),  // 89: Ac (rough estimate)
    Some(2.93), // 90: Th
    Some(2.88), // 91: Pa (rough estimate)
    Some(2.71), // 92: U
    Some(2.82), // 93: Np
    Some(2.81), // 94: Pu
    Some(2.83), // 95: Am
    Some(3.05), // 96: Cm (rough estimate)
    Some(3.4),  // 97: Bk (rough estimate)
    Some(3.05), // 98: Cf (rough estimate)
    Some(2.7),  // 99: Es (rough estimate)
    None,       // 100: Fm
    None,       // 101: Md
    None,       // 102: No
    None,       // 103: Lr
    None,       // 104: Rf
    None,       // 105: Db
    None,       // 106: Sg
    None,       // 107: Bh
    None,       // 108: Hs
    None,       // 109: Mt
    None,       // 110: Ds
    None,       // 111: Rg
    None,       // 112: Cn
    None,       // 113: Nh
    None,       // 114: Fl
    None,       // 115: Mc
    None,       // 116: Lv
    None,       // 117: Ts
    None,       // 118: Og
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::ANGSTROM;

    #[test]
    fn carbon_reference_values_match_archived_sources() {
        let carbon = Element::from_atomic_number(6).unwrap();

        assert_eq!(
            carbon.covalent_radius(),
            Some(Quantity::new(0.76, ANGSTROM))
        );
        assert_eq!(
            carbon.dipole_polarizability(),
            Some(Quantity::new(11.3, CUBIC_BOHR))
        );
        assert!(
            (carbon
                .van_der_waals_radius(VanDerWaalsRadiusSource::FreeAtom)
                .unwrap()
                .value_in(ANGSTROM)
                .unwrap()
                - 1.91)
                .abs()
                < 0.02
        );
    }

    #[test]
    fn requested_reference_radius_falls_back_to_free_atom_value() {
        let neodymium = Element::from_atomic_number(60).unwrap();
        assert!(REFERENCE_VDW_RADIUS_BOHR[60].is_none());
        assert!(
            (neodymium
                .van_der_waals_radius(VanDerWaalsRadiusSource::Reference)
                .unwrap()
                .value_in(ANGSTROM)
                .unwrap()
                - 2.90)
                .abs()
                < 0.03
        );
    }

    fn angstrom(element: u8, source: VanDerWaalsRadiusSource) -> Option<f64> {
        Element::from_atomic_number(element)
            .unwrap()
            .van_der_waals_radius(source)
            .map(|radius| radius.value_in(ANGSTROM).unwrap())
    }

    #[test]
    fn bondi_column_follows_mantina_and_alvarez_transcriptions() {
        use VanDerWaalsRadiusSource::Bondi;
        // Rowland and Taylor hydrogen, Bondi carbon and halogens as Mantina
        // Table 12 gives them, a Mantina completion, and a Bondi metal.
        for (element, radius) in [(1, 1.10), (6, 1.70), (35, 1.83), (38, 2.49), (80, 1.70)] {
            assert_eq!(angstrom(element, Bondi), Some(radius), "element {element}");
        }
        let covered = (1..=118).filter(|&z| angstrom(z, Bondi).is_some()).count();
        assert_eq!(covered, 54);
        // No fallback to another column: iron and the lanthanides have no
        // Bondi radius.
        assert_eq!(angstrom(26, Bondi), None);
        assert_eq!(angstrom(64, Bondi), None);
    }

    #[test]
    fn alvarez_column_follows_table_one() {
        use VanDerWaalsRadiusSource::Alvarez;
        for (element, radius) in [(1, 1.20), (6, 1.77), (26, 2.44), (54, 2.06), (99, 2.7)] {
            assert_eq!(
                angstrom(element, Alvarez),
                Some(radius),
                "element {element}"
            );
        }
        let missing = (1..=118)
            .filter(|&z| angstrom(z, Alvarez).is_none())
            .collect::<Vec<_>>();
        let expected = [61, 84, 85, 86, 87, 88]
            .into_iter()
            .chain(100..=118)
            .collect::<Vec<_>>();
        assert_eq!(missing, expected);
    }

    #[test]
    fn charry_tkatchenko_conversion_is_unchanged_for_solvation() {
        let oxygen = Element::from_atomic_number(8).unwrap();
        assert_eq!(
            oxygen.van_der_waals_radius_angstrom(VanDerWaalsRadiusSource::Reference),
            Some(3.2314 * BOHR_TO_ANGSTROM)
        );
        assert_eq!(
            oxygen.van_der_waals_radius_angstrom(VanDerWaalsRadiusSource::Bondi),
            Some(1.52)
        );
        assert_eq!(
            VanDerWaalsRadiusSource::default(),
            VanDerWaalsRadiusSource::FreeAtom
        );
    }

    #[test]
    fn missing_both_vdw_sources_is_explicit() {
        let livermorium = Element::from_atomic_number(116).unwrap();
        assert_eq!(
            livermorium.van_der_waals_radius(VanDerWaalsRadiusSource::FreeAtom),
            None
        );
    }
}
