use super::Element;

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;

/// Selects one of the neutral-atom van der Waals reference columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum VanDerWaalsRadiusSource {
    /// Free-atom radii reported by Charry and Tkatchenko.
    #[default]
    FreeAtom,
    /// The compiled bonded/reference-radius column reported in the same work.
    Reference,
}

/// Authoritative periodic-table reference values carried by Molecular.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElementReferenceData {
    pub covalent_radius_angstrom: Option<f64>,
    pub dipole_polarizability_bohr3: Option<f64>,
    pub free_atom_vdw_radius_bohr: Option<f64>,
    pub reference_vdw_radius_bohr: Option<f64>,
}

impl Element {
    /// Returns the archived element reference values.
    ///
    /// Covalent radii follow Cordero et al., Dalton Trans. 2008,
    /// DOI 10.1039/B801115J. The van der Waals radii follow Charry and
    /// Tkatchenko, J. Chem. Theory Comput. 2024,
    /// DOI 10.1021/acs.jctc.4c00784. Dipole polarizabilities are the values
    /// compiled there from Schwerdtfeger and Nagle, Mol. Phys. 2019,
    /// DOI 10.1080/00268976.2018.1535143.
    pub fn reference_data(self) -> ElementReferenceData {
        let index = usize::from(self.atomic_number());
        ElementReferenceData {
            covalent_radius_angstrom: self.covalent_radius_angstrom(),
            dipole_polarizability_bohr3: DIPOLE_POLARIZABILITY_BOHR3[index],
            free_atom_vdw_radius_bohr: TKATCHENKO_VDW_RADIUS_BOHR[index],
            reference_vdw_radius_bohr: REFERENCE_VDW_RADIUS_BOHR[index],
        }
    }

    pub fn dipole_polarizability_bohr3(self) -> Option<f64> {
        DIPOLE_POLARIZABILITY_BOHR3[usize::from(self.atomic_number())]
    }

    /// Returns a van der Waals radius in bohr.
    ///
    /// If the requested reference column is unavailable, the other archived
    /// column is used. `None` is returned only when Molecular has neither.
    pub fn van_der_waals_radius_bohr(self, source: VanDerWaalsRadiusSource) -> Option<f64> {
        let index = usize::from(self.atomic_number());
        match source {
            VanDerWaalsRadiusSource::FreeAtom => {
                TKATCHENKO_VDW_RADIUS_BOHR[index].or(REFERENCE_VDW_RADIUS_BOHR[index])
            }
            VanDerWaalsRadiusSource::Reference => {
                REFERENCE_VDW_RADIUS_BOHR[index].or(TKATCHENKO_VDW_RADIUS_BOHR[index])
            }
        }
    }

    pub fn van_der_waals_radius_angstrom(self, source: VanDerWaalsRadiusSource) -> Option<f64> {
        self.van_der_waals_radius_bohr(source)
            .map(|radius| radius * BOHR_TO_ANGSTROM)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn carbon_reference_values_match_archived_sources() {
        let carbon = Element::from_atomic_number(6).unwrap();
        let reference = carbon.reference_data();

        assert_eq!(reference.covalent_radius_angstrom, Some(0.76));
        assert_eq!(reference.dipole_polarizability_bohr3, Some(11.3));
        assert!(
            (carbon
                .van_der_waals_radius_angstrom(VanDerWaalsRadiusSource::FreeAtom)
                .unwrap()
                - 1.91)
                .abs()
                < 0.02
        );
    }

    #[test]
    fn requested_reference_radius_falls_back_to_free_atom_value() {
        let neodymium = Element::from_atomic_number(60).unwrap();
        assert!(neodymium
            .reference_data()
            .reference_vdw_radius_bohr
            .is_none());
        assert!(
            (neodymium
                .van_der_waals_radius_angstrom(VanDerWaalsRadiusSource::Reference)
                .unwrap()
                - 2.90)
                .abs()
                < 0.03
        );
    }

    #[test]
    fn missing_both_vdw_sources_is_explicit() {
        let livermorium = Element::from_atomic_number(116).unwrap();
        assert_eq!(
            livermorium.van_der_waals_radius_bohr(VanDerWaalsRadiusSource::FreeAtom),
            None
        );
    }
}
