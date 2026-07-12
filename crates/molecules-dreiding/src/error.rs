use std::fmt;

use molecules::core::{AtomId, BondId, BondOrder};
use molecules::modeling::ComponentId;

/// Failure while converting and parameterizing a molecular model with DREIDING.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum DreidingPrepareError {
    UnresolvedImplicitHydrogens {
        atom: AtomId,
    },
    CountedHydrogens {
        atom: AtomId,
        explicit: u8,
        implicit: u8,
    },
    RadicalAtom {
        atom: AtomId,
    },
    UnsupportedBondOrder {
        bond: BondId,
        order: BondOrder,
    },
    InconsistentAromaticBond {
        bond: BondId,
    },
    CrossComponentBond {
        bond: BondId,
    },
    UnsupportedElement {
        atom: AtomId,
        symbol: String,
    },
    Parameterization {
        component: Option<ComponentId>,
        message: String,
    },
    AtomTypeMismatch {
        component: ComponentId,
        atom: AtomId,
        whole_model: String,
        component_model: String,
    },
    MissingVdwParameters {
        first: AtomId,
        second: AtomId,
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
            Self::CrossComponentBond { bond } => {
                write!(f, "bond {bond} connects different model components")
            }
            Self::UnsupportedElement { atom, symbol } => {
                write!(
                    f,
                    "atom {atom} element {symbol} cannot be converted for DREIDING"
                )
            }
            Self::Parameterization { component, message } => match component {
                Some(component) => write!(
                    f,
                    "DREIDING parameterization failed for {component}: {message}"
                ),
                None => write!(
                    f,
                    "DREIDING parameterization failed for the model: {message}"
                ),
            },
            Self::AtomTypeMismatch {
                component,
                atom,
                whole_model,
                component_model,
            } => write!(
                f,
                "DREIDING atom type for {atom} in {component} differs between whole-model ({whole_model}) and component ({component_model}) preparation"
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
