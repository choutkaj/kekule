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
