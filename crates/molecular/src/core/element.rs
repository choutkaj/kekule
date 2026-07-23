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
    pub fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        if (1..=118).contains(&atomic_number) {
            Some(Self { atomic_number })
        } else {
            None
        }
    }

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

    /// Returns the single-bond covalent radius in ångström when a reference
    /// value is available.
    ///
    /// Carbon uses the general-purpose sp3 value. Transition-metal values use
    /// the conservative low-spin entry where separate spin-state values exist.
    pub const fn covalent_radius_angstrom(self) -> Option<f64> {
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

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}
