use std::fmt;

/// Stable identity for a trajectory file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TrajectoryFormat {
    /// Multi-frame extended XYZ text.
    Xyz,
    /// CHARMM/NAMD/LAMMPS binary DCD.
    Dcd,
    /// GROMACS compressed XTC.
    Xtc,
    /// GROMACS full-precision TRR.
    Trr,
}

impl TrajectoryFormat {
    /// The conventional upper-case format name, such as `"XYZ"`; also its
    /// `Display` form.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Xyz => "XYZ",
            Self::Dcd => "DCD",
            Self::Xtc => "XTC",
            Self::Trr => "TRR",
        }
    }
}

impl fmt::Display for TrajectoryFormat {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}
