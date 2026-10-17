use crate::units::{Quantity, ANGSTROM};
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Element {
    atomic_number: u8,
}

const ELEMENT_SYMBOLS: [&str; 119] = [
    "?", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

impl Element {
    /// H, atomic number 1.
    pub const H: Self = Self { atomic_number: 1 };
    /// He, atomic number 2.
    pub const HE: Self = Self { atomic_number: 2 };
    /// Li, atomic number 3.
    pub const LI: Self = Self { atomic_number: 3 };
    /// Be, atomic number 4.
    pub const BE: Self = Self { atomic_number: 4 };
    /// B, atomic number 5.
    pub const B: Self = Self { atomic_number: 5 };
    /// C, atomic number 6.
    pub const C: Self = Self { atomic_number: 6 };
    /// N, atomic number 7.
    pub const N: Self = Self { atomic_number: 7 };
    /// O, atomic number 8.
    pub const O: Self = Self { atomic_number: 8 };
    /// F, atomic number 9.
    pub const F: Self = Self { atomic_number: 9 };
    /// Ne, atomic number 10.
    pub const NE: Self = Self { atomic_number: 10 };
    /// Na, atomic number 11.
    pub const NA: Self = Self { atomic_number: 11 };
    /// Mg, atomic number 12.
    pub const MG: Self = Self { atomic_number: 12 };
    /// Al, atomic number 13.
    pub const AL: Self = Self { atomic_number: 13 };
    /// Si, atomic number 14.
    pub const SI: Self = Self { atomic_number: 14 };
    /// P, atomic number 15.
    pub const P: Self = Self { atomic_number: 15 };
    /// S, atomic number 16.
    pub const S: Self = Self { atomic_number: 16 };
    /// Cl, atomic number 17.
    pub const CL: Self = Self { atomic_number: 17 };
    /// Ar, atomic number 18.
    pub const AR: Self = Self { atomic_number: 18 };
    /// K, atomic number 19.
    pub const K: Self = Self { atomic_number: 19 };
    /// Ca, atomic number 20.
    pub const CA: Self = Self { atomic_number: 20 };
    /// Sc, atomic number 21.
    pub const SC: Self = Self { atomic_number: 21 };
    /// Ti, atomic number 22.
    pub const TI: Self = Self { atomic_number: 22 };
    /// V, atomic number 23.
    pub const V: Self = Self { atomic_number: 23 };
    /// Cr, atomic number 24.
    pub const CR: Self = Self { atomic_number: 24 };
    /// Mn, atomic number 25.
    pub const MN: Self = Self { atomic_number: 25 };
    /// Fe, atomic number 26.
    pub const FE: Self = Self { atomic_number: 26 };
    /// Co, atomic number 27.
    pub const CO: Self = Self { atomic_number: 27 };
    /// Ni, atomic number 28.
    pub const NI: Self = Self { atomic_number: 28 };
    /// Cu, atomic number 29.
    pub const CU: Self = Self { atomic_number: 29 };
    /// Zn, atomic number 30.
    pub const ZN: Self = Self { atomic_number: 30 };
    /// Ga, atomic number 31.
    pub const GA: Self = Self { atomic_number: 31 };
    /// Ge, atomic number 32.
    pub const GE: Self = Self { atomic_number: 32 };
    /// As, atomic number 33.
    pub const AS: Self = Self { atomic_number: 33 };
    /// Se, atomic number 34.
    pub const SE: Self = Self { atomic_number: 34 };
    /// Br, atomic number 35.
    pub const BR: Self = Self { atomic_number: 35 };
    /// Kr, atomic number 36.
    pub const KR: Self = Self { atomic_number: 36 };
    /// Rb, atomic number 37.
    pub const RB: Self = Self { atomic_number: 37 };
    /// Sr, atomic number 38.
    pub const SR: Self = Self { atomic_number: 38 };
    /// Y, atomic number 39.
    pub const Y: Self = Self { atomic_number: 39 };
    /// Zr, atomic number 40.
    pub const ZR: Self = Self { atomic_number: 40 };
    /// Nb, atomic number 41.
    pub const NB: Self = Self { atomic_number: 41 };
    /// Mo, atomic number 42.
    pub const MO: Self = Self { atomic_number: 42 };
    /// Tc, atomic number 43.
    pub const TC: Self = Self { atomic_number: 43 };
    /// Ru, atomic number 44.
    pub const RU: Self = Self { atomic_number: 44 };
    /// Rh, atomic number 45.
    pub const RH: Self = Self { atomic_number: 45 };
    /// Pd, atomic number 46.
    pub const PD: Self = Self { atomic_number: 46 };
    /// Ag, atomic number 47.
    pub const AG: Self = Self { atomic_number: 47 };
    /// Cd, atomic number 48.
    pub const CD: Self = Self { atomic_number: 48 };
    /// In, atomic number 49.
    pub const IN: Self = Self { atomic_number: 49 };
    /// Sn, atomic number 50.
    pub const SN: Self = Self { atomic_number: 50 };
    /// Sb, atomic number 51.
    pub const SB: Self = Self { atomic_number: 51 };
    /// Te, atomic number 52.
    pub const TE: Self = Self { atomic_number: 52 };
    /// I, atomic number 53.
    pub const I: Self = Self { atomic_number: 53 };
    /// Xe, atomic number 54.
    pub const XE: Self = Self { atomic_number: 54 };
    /// Cs, atomic number 55.
    pub const CS: Self = Self { atomic_number: 55 };
    /// Ba, atomic number 56.
    pub const BA: Self = Self { atomic_number: 56 };
    /// La, atomic number 57.
    pub const LA: Self = Self { atomic_number: 57 };
    /// Ce, atomic number 58.
    pub const CE: Self = Self { atomic_number: 58 };
    /// Pr, atomic number 59.
    pub const PR: Self = Self { atomic_number: 59 };
    /// Nd, atomic number 60.
    pub const ND: Self = Self { atomic_number: 60 };
    /// Pm, atomic number 61.
    pub const PM: Self = Self { atomic_number: 61 };
    /// Sm, atomic number 62.
    pub const SM: Self = Self { atomic_number: 62 };
    /// Eu, atomic number 63.
    pub const EU: Self = Self { atomic_number: 63 };
    /// Gd, atomic number 64.
    pub const GD: Self = Self { atomic_number: 64 };
    /// Tb, atomic number 65.
    pub const TB: Self = Self { atomic_number: 65 };
    /// Dy, atomic number 66.
    pub const DY: Self = Self { atomic_number: 66 };
    /// Ho, atomic number 67.
    pub const HO: Self = Self { atomic_number: 67 };
    /// Er, atomic number 68.
    pub const ER: Self = Self { atomic_number: 68 };
    /// Tm, atomic number 69.
    pub const TM: Self = Self { atomic_number: 69 };
    /// Yb, atomic number 70.
    pub const YB: Self = Self { atomic_number: 70 };
    /// Lu, atomic number 71.
    pub const LU: Self = Self { atomic_number: 71 };
    /// Hf, atomic number 72.
    pub const HF: Self = Self { atomic_number: 72 };
    /// Ta, atomic number 73.
    pub const TA: Self = Self { atomic_number: 73 };
    /// W, atomic number 74.
    pub const W: Self = Self { atomic_number: 74 };
    /// Re, atomic number 75.
    pub const RE: Self = Self { atomic_number: 75 };
    /// Os, atomic number 76.
    pub const OS: Self = Self { atomic_number: 76 };
    /// Ir, atomic number 77.
    pub const IR: Self = Self { atomic_number: 77 };
    /// Pt, atomic number 78.
    pub const PT: Self = Self { atomic_number: 78 };
    /// Au, atomic number 79.
    pub const AU: Self = Self { atomic_number: 79 };
    /// Hg, atomic number 80.
    pub const HG: Self = Self { atomic_number: 80 };
    /// Tl, atomic number 81.
    pub const TL: Self = Self { atomic_number: 81 };
    /// Pb, atomic number 82.
    pub const PB: Self = Self { atomic_number: 82 };
    /// Bi, atomic number 83.
    pub const BI: Self = Self { atomic_number: 83 };
    /// Po, atomic number 84.
    pub const PO: Self = Self { atomic_number: 84 };
    /// At, atomic number 85.
    pub const AT: Self = Self { atomic_number: 85 };
    /// Rn, atomic number 86.
    pub const RN: Self = Self { atomic_number: 86 };
    /// Fr, atomic number 87.
    pub const FR: Self = Self { atomic_number: 87 };
    /// Ra, atomic number 88.
    pub const RA: Self = Self { atomic_number: 88 };
    /// Ac, atomic number 89.
    pub const AC: Self = Self { atomic_number: 89 };
    /// Th, atomic number 90.
    pub const TH: Self = Self { atomic_number: 90 };
    /// Pa, atomic number 91.
    pub const PA: Self = Self { atomic_number: 91 };
    /// U, atomic number 92.
    pub const U: Self = Self { atomic_number: 92 };
    /// Np, atomic number 93.
    pub const NP: Self = Self { atomic_number: 93 };
    /// Pu, atomic number 94.
    pub const PU: Self = Self { atomic_number: 94 };
    /// Am, atomic number 95.
    pub const AM: Self = Self { atomic_number: 95 };
    /// Cm, atomic number 96.
    pub const CM: Self = Self { atomic_number: 96 };
    /// Bk, atomic number 97.
    pub const BK: Self = Self { atomic_number: 97 };
    /// Cf, atomic number 98.
    pub const CF: Self = Self { atomic_number: 98 };
    /// Es, atomic number 99.
    pub const ES: Self = Self { atomic_number: 99 };
    /// Fm, atomic number 100.
    pub const FM: Self = Self { atomic_number: 100 };
    /// Md, atomic number 101.
    pub const MD: Self = Self { atomic_number: 101 };
    /// No, atomic number 102.
    pub const NO: Self = Self { atomic_number: 102 };
    /// Lr, atomic number 103.
    pub const LR: Self = Self { atomic_number: 103 };
    /// Rf, atomic number 104.
    pub const RF: Self = Self { atomic_number: 104 };
    /// Db, atomic number 105.
    pub const DB: Self = Self { atomic_number: 105 };
    /// Sg, atomic number 106.
    pub const SG: Self = Self { atomic_number: 106 };
    /// Bh, atomic number 107.
    pub const BH: Self = Self { atomic_number: 107 };
    /// Hs, atomic number 108.
    pub const HS: Self = Self { atomic_number: 108 };
    /// Mt, atomic number 109.
    pub const MT: Self = Self { atomic_number: 109 };
    /// Ds, atomic number 110.
    pub const DS: Self = Self { atomic_number: 110 };
    /// Rg, atomic number 111.
    pub const RG: Self = Self { atomic_number: 111 };
    /// Cn, atomic number 112.
    pub const CN: Self = Self { atomic_number: 112 };
    /// Nh, atomic number 113.
    pub const NH: Self = Self { atomic_number: 113 };
    /// Fl, atomic number 114.
    pub const FL: Self = Self { atomic_number: 114 };
    /// Mc, atomic number 115.
    pub const MC: Self = Self { atomic_number: 115 };
    /// Lv, atomic number 116.
    pub const LV: Self = Self { atomic_number: 116 };
    /// Ts, atomic number 117.
    pub const TS: Self = Self { atomic_number: 117 };
    /// Og, atomic number 118.
    pub const OG: Self = Self { atomic_number: 118 };

    /// The element with this atomic number, from 1 to 118.
    pub fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        if (1..=118).contains(&atomic_number) {
            Some(Self { atomic_number })
        } else {
            None
        }
    }

    /// The element with this IUPAC symbol, matched case-sensitively, such as
    /// `"Cl"`.
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        let atomic_number = ELEMENT_SYMBOLS
            .iter()
            .position(|candidate| *candidate == symbol)?;
        if atomic_number == 0 {
            return None;
        }
        Some(Self {
            atomic_number: atomic_number as u8,
        })
    }

    pub const fn atomic_number(self) -> u8 {
        self.atomic_number
    }

    /// Returns the single-bond covalent radius when a reference value is
    /// available, in ångström as tabulated.
    ///
    /// Carbon uses the general-purpose sp3 value. Transition-metal values use
    /// the conservative low-spin entry where separate spin-state values exist.
    /// Covalent radii follow Cordero et al., Dalton Trans. 2008,
    /// DOI 10.1039/B801115J.
    pub const fn covalent_radius(self) -> Option<Quantity<f64>> {
        match self.covalent_radius_angstrom() {
            Some(radius) => Some(Quantity::new(radius, ANGSTROM)),
            None => None,
        }
    }

    /// The tabulated covalent radius in ångström, for internal geometry.
    pub(crate) const fn covalent_radius_angstrom(self) -> Option<f64> {
        match self.atomic_number {
            1 => Some(0.31),
            2 => Some(0.28),
            3 => Some(1.28),
            4 => Some(0.96),
            5 => Some(0.84),
            6 => Some(0.76),
            7 => Some(0.71),
            8 => Some(0.66),
            9 => Some(0.57),
            10 => Some(0.58),
            11 => Some(1.66),
            12 => Some(1.41),
            13 => Some(1.21),
            14 => Some(1.11),
            15 => Some(1.07),
            16 => Some(1.05),
            17 => Some(1.02),
            18 => Some(1.06),
            19 => Some(2.03),
            20 => Some(1.76),
            21 => Some(1.70),
            22 => Some(1.60),
            23 => Some(1.53),
            24 => Some(1.39),
            25 => Some(1.39),
            26 => Some(1.32),
            27 => Some(1.26),
            28 => Some(1.24),
            29 => Some(1.32),
            30 => Some(1.22),
            31 => Some(1.22),
            32 => Some(1.20),
            33 => Some(1.19),
            34 => Some(1.20),
            35 => Some(1.20),
            36 => Some(1.16),
            37 => Some(2.20),
            38 => Some(1.95),
            39 => Some(1.90),
            40 => Some(1.75),
            41 => Some(1.64),
            42 => Some(1.54),
            43 => Some(1.47),
            44 => Some(1.46),
            45 => Some(1.42),
            46 => Some(1.39),
            47 => Some(1.45),
            48 => Some(1.44),
            49 => Some(1.42),
            50 => Some(1.39),
            51 => Some(1.39),
            52 => Some(1.38),
            53 => Some(1.39),
            54 => Some(1.40),
            55 => Some(2.44),
            56 => Some(2.15),
            57 => Some(2.07),
            58 => Some(2.04),
            59 => Some(2.03),
            60 => Some(2.01),
            61 => Some(1.99),
            62 => Some(1.98),
            63 => Some(1.98),
            64 => Some(1.96),
            65 => Some(1.94),
            66 => Some(1.92),
            67 => Some(1.92),
            68 => Some(1.89),
            69 => Some(1.90),
            70 => Some(1.87),
            71 => Some(1.87),
            72 => Some(1.75),
            73 => Some(1.70),
            74 => Some(1.62),
            75 => Some(1.51),
            76 => Some(1.44),
            77 => Some(1.41),
            78 => Some(1.36),
            79 => Some(1.36),
            80 => Some(1.32),
            81 => Some(1.45),
            82 => Some(1.46),
            83 => Some(1.48),
            84 => Some(1.40),
            85 => Some(1.50),
            86 => Some(1.50),
            87 => Some(2.60),
            88 => Some(2.21),
            89 => Some(2.15),
            90 => Some(2.06),
            91 => Some(2.00),
            92 => Some(1.96),
            93 => Some(1.90),
            94 => Some(1.87),
            95 => Some(1.80),
            96 => Some(1.69),
            _ => None,
        }
    }

    pub fn symbol(self) -> &'static str {
        ELEMENT_SYMBOLS
            .get(self.atomic_number as usize)
            .copied()
            .unwrap_or("?")
    }
}

impl std::str::FromStr for Element {
    type Err = ElementError;

    /// Parses an IUPAC symbol; see [`Element::from_symbol`].
    fn from_str(symbol: &str) -> Result<Self, Self::Err> {
        Self::from_symbol(symbol).ok_or_else(|| ElementError::UnknownSymbol(symbol.to_owned()))
    }
}

impl TryFrom<u8> for Element {
    type Error = ElementError;

    /// Converts an atomic number; see [`Element::from_atomic_number`].
    fn try_from(atomic_number: u8) -> Result<Self, Self::Error> {
        Self::from_atomic_number(atomic_number)
            .ok_or(ElementError::InvalidAtomicNumber(atomic_number))
    }
}

/// A symbol or atomic number names no element.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ElementError {
    /// The text is not an IUPAC element symbol.
    UnknownSymbol(String),
    /// The atomic number is outside 1 to 118.
    InvalidAtomicNumber(u8),
}

impl fmt::Display for ElementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownSymbol(symbol) => write!(f, "unknown element symbol {symbol:?}"),
            Self::InvalidAtomicNumber(number) => {
                write!(f, "atomic number {number} is outside 1 to 118")
            }
        }
    }
}

impl std::error::Error for ElementError {}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}
