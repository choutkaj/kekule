use std::fmt;

use molecules::core::BondOrder;
use molecules::modeling::{InstanceAtomId, InstanceBondId, MoleculeInstanceId};

/// Failure while converting and parameterizing a molecular model with DREIDING.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum DreidingPrepareError {
    UnresolvedImplicitHydrogens {
        atom: InstanceAtomId,
    },
    CountedHydrogens {
        atom: InstanceAtomId,
        explicit: u8,
        implicit: u8,
    },
    RadicalAtom {
        atom: InstanceAtomId,
    },
    UnsupportedBondOrder {
        bond: InstanceBondId,
        order: BondOrder,
    },
    InconsistentAromaticBond {
        bond: InstanceBondId,
    },
    UnsupportedElement {
        atom: InstanceAtomId,
        symbol: String,
    },
    Parameterization {
        molecule: Option<MoleculeInstanceId>,
        message: String,
    },
    AtomTypeMismatch {
        molecule: MoleculeInstanceId,
        atom: InstanceAtomId,
        whole_model: String,
        component_model: String,
    },
    MissingVdwParameters {
        first: InstanceAtomId,
        second: InstanceAtomId,
    },
    InvalidPreparedData {
        interaction: &'static str,
        detail: String,
    },
}

impl fmt::Display for DreidingPrepareError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnresolvedImplicitHydrogens { atom } => write!(
                f,
                "atom {atom} has an unresolved implicit-hydrogen count; DREIDING preparation requires an explicit zero count or no-implicit-hydrogens assertion"
            ),
            Self::CountedHydrogens {
                atom,
                explicit,
                implicit,
            } => write!(
                f,
                "atom {atom} has {explicit} explicit-count and {implicit} implicit hydrogen(s); DREIDING preparation requires coordinate-bearing hydrogen atoms"
            ),
            Self::RadicalAtom { atom } => {
                write!(
                    f,
                    "atom {atom} is radical; DREIDING radical parameters are unsupported"
                )
            }
            Self::UnsupportedBondOrder { bond, order } => {
                write!(
                    f,
                    "bond {bond} has unsupported DREIDING bond order {order:?}"
                )
            }
            Self::InconsistentAromaticBond { bond } => write!(
                f,
                "bond {bond} has inconsistent aromatic flag and aromatic bond order"
            ),
            Self::UnsupportedElement { atom, symbol } => {
                write!(
                    f,
                    "atom {atom} element {symbol} cannot be converted for DREIDING"
                )
            }
            Self::Parameterization { molecule, message } => match molecule {
                Some(molecule) => write!(
                    f,
                    "DREIDING parameterization failed for {molecule}: {message}"
                ),
                None => write!(
                    f,
                    "DREIDING parameterization failed for the model: {message}"
                ),
            },
            Self::AtomTypeMismatch {
                molecule,
                atom,
                whole_model,
                component_model,
            } => write!(
                f,
                "DREIDING atom type for {atom} in {molecule} differs between whole-model ({whole_model}) and molecule-local ({component_model}) preparation"
            ),
            Self::MissingVdwParameters { first, second } => write!(
                f,
                "DREIDING produced no van der Waals parameters for pair {first}-{second}"
            ),
            Self::InvalidPreparedData {
                interaction,
                detail,
            } => write!(
                f,
                "DREIDING produced invalid {interaction} parameters: {detail}"
            ),
        }
    }
}

impl std::error::Error for DreidingPrepareError {}
