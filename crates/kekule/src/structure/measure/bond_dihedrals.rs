use std::{fmt, sync::Arc};

use crate::algorithms::{bond_reference_atoms, CipAssignmentOptions, CipRankingError};
use crate::structure::ModelView;
use crate::topology::{InstanceAtomId, InstanceBondId, Topology};
use crate::units::Quantity;

use super::{dihedral, MeasurementError};

/// Failure to prepare or measure a deterministic bond dihedral.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum BondDihedralError {
    InvalidBondId(InstanceBondId),
    /// A prepared definition can only measure views sharing its topology layout.
    TopologyMismatch,
    Ranking {
        bond: InstanceBondId,
        error: CipRankingError,
    },
    Measurement(MeasurementError),
}

impl fmt::Display for BondDihedralError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBondId(bond) => write!(f, "invalid dihedral bond: {bond}"),
            Self::TopologyMismatch => f.write_str("bond dihedral belongs to a different topology"),
            Self::Ranking { bond, error } => {
                write!(f, "dihedral reference ranking for {bond}: {error}")
            }
            Self::Measurement(error) => write!(f, "bond dihedral: {error}"),
        }
    }
}

impl std::error::Error for BondDihedralError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Ranking { error, .. } => Some(error),
            Self::Measurement(error) => Some(error),
            _ => None,
        }
    }
}

/// A deterministic reference quartet for a bond, bound to one topology layout.
///
/// For B-C, choose A among B's explicit neighbors other than C, and D among
/// C's explicit neighbors other than B. Highest CIP priority wins; complete
/// priority ties use the smallest local atom ID. B precedes C by atom ID.
/// Bond order and rotatability do not restrict this operation. References are
/// independent of coordinates and are not changed to avoid degenerate geometry.
/// Renumbering atoms can change references, so this is not a numbering-invariant
/// canonical angle. Reversing all four atoms preserves the dihedral sign.
///
/// Implicit hydrogens and phantom atoms have no stored positions and cannot be
/// references. Explicit hydrogens participate normally. When ranking is needed,
/// implicit hydrogen counts throughout the molecule must be known (from fixed
/// declarations or explicit valence perception). Ranking uses represented stereo
/// and path-dependent auxiliary descriptors without installing perception.
/// Unresolved auxiliary stereo and expansion limits return errors, never ID ties.
///
/// [`Self::atoms`] is `None` if an endpoint lacks an explicit reference or if
/// the selected outer atoms coincide (a three-atom closed walk). Coordinate
/// degeneracy is determined separately by [`Self::measure`].
///
/// Prepare once, then measure models, ensemble members, or trajectory frames
/// sharing this topology:
/// ```
/// use kekule::structure::{measure::BondDihedral, Model};
/// # fn example(model: &Model) -> Result<(), Box<dyn std::error::Error>> {
/// let topology = model.shared_topology();
/// let definitions = topology.bond_ids().iter()
///     .map(|&bond| BondDihedral::new(&topology, bond))
///     .collect::<Result<Vec<_>, _>>()?;
/// for definition in &definitions {
///     let angle = definition.measure(model.as_model_view())?;
///     // `angle` is None when this bond has no defined dihedral in this frame.
/// }
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct BondDihedral {
    topology: Arc<Topology>,
    bond: InstanceBondId,
    atoms: Option<[InstanceAtomId; 4]>,
}

impl BondDihedral {
    /// Prepares references with the default CIP expansion bounds.
    pub fn new(topology: &Arc<Topology>, bond: InstanceBondId) -> Result<Self, BondDihedralError> {
        Self::with_options(topology, bond, CipAssignmentOptions::default())
    }

    /// Prepares references with explicit CIP expansion bounds.
    pub fn with_options(
        topology: &Arc<Topology>,
        bond: InstanceBondId,
        options: CipAssignmentOptions,
    ) -> Result<Self, BondDihedralError> {
        let instance = topology
            .bond(bond)
            .ok_or(BondDihedralError::InvalidBondId(bond))?
            .molecule();
        let atoms = bond_reference_atoms(instance.molecule(), bond.bond(), options)
            .map_err(|error| BondDihedralError::Ranking { bond, error })?
            .map(|atoms| atoms.map(|atom| InstanceAtomId::new(bond.molecule(), atom)));
        Ok(Self {
            topology: Arc::clone(topology),
            bond,
            atoms,
        })
    }

    pub fn bond(&self) -> InstanceBondId {
        self.bond
    }

    /// Selected A-B-C-D atoms, independent of the geometry of any frame.
    pub fn atoms(&self) -> Option<[InstanceAtomId; 4]> {
        self.atoms
    }

    /// Measures the fixed quartet in radians, with ordinary quantity conversions.
    /// Returns `None` for absent references, zero-length consecutive coordinate
    /// differences, or collinear defining triples. Invalid topology and numerical
    /// overflow remain errors. Periodic cells are ignored; signed angles retain
    /// the ordinary -pi/pi branch cut. No angular unwrapping is applied.
    pub fn measure(&self, view: ModelView<'_>) -> Result<Option<Quantity<f64>>, BondDihedralError> {
        if !self.topology.shares_layout(view.topology()) {
            return Err(BondDihedralError::TopologyMismatch);
        }
        let Some([a, b, c, d]) = self.atoms else {
            return Ok(None);
        };
        match dihedral(view, a, b, c, d) {
            Ok(angle) => Ok(Some(angle)),
            Err(MeasurementError::DegenerateGeometry) => Ok(None),
            Err(error) => Err(BondDihedralError::Measurement(error)),
        }
    }
}

/// Selects and measures one deterministic bond dihedral using default CIP bounds.
/// See [`BondDihedral`] for reference selection and the meaning of `None`.
/// Reuse a prepared definition to avoid reranking on every trajectory frame.
pub fn bond_dihedral(
    view: ModelView<'_>,
    bond: InstanceBondId,
) -> Result<Option<Quantity<f64>>, BondDihedralError> {
    BondDihedral::new(&view.shared_topology(), bond)?.measure(view)
}

/// Lazily measures every bond in topology dense order, retaining bond identities
/// and individual errors. Bonds without a defined dihedral are not filtered out.
pub fn bond_dihedrals(
    view: ModelView<'_>,
) -> impl ExactSizeIterator<
    Item = (
        InstanceBondId,
        Result<Option<Quantity<f64>>, BondDihedralError>,
    ),
> + '_ {
    view.topology()
        .bond_ids()
        .iter()
        .copied()
        .map(move |bond| (bond, bond_dihedral(view, bond)))
}
