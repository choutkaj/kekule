//! Periodic box construction and force-field-independent TIP3P solvation.
//!
//! These operations create starting structures, not parameterized or equilibrated
//! simulation systems. Existing coordinates, chemistry, and entity annotations are
//! preserved. Solvation publishes one new topology; old selections and parameter
//! assignments remain bound to the original snapshot. Existing owner properties
//! are cleared under the ordinary append contract and listed in the optional report.
//!
mod boxes;
mod neighbors;
mod packing;
pub use boxes::*;

use super::{Model, ModelBuildError, PositionError};
use crate::{
    geometry::{PeriodicGeometryError, Point3},
    properties::PropertyKey,
    topology::{HierarchyError, InstanceAtomId, Topology},
    units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT, MOLAR, NANOMETER},
};
use std::{fmt, sync::Arc};

/// Built-in coordinate model. This does not install force-field parameters.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum WaterModel {
    #[default]
    Tip3p,
}
/// Monovalent cations supported by the solvent builder.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum PositiveIon {
    Lithium,
    #[default]
    Sodium,
    Potassium,
    Rubidium,
    Cesium,
}
impl PositiveIon {
    pub(crate) fn symbol(self) -> &'static str {
        match self {
            Self::Lithium => "Li",
            Self::Sodium => "Na",
            Self::Potassium => "K",
            Self::Rubidium => "Rb",
            Self::Cesium => "Cs",
        }
    }
}
/// Monovalent anions supported by the solvent builder.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum NegativeIon {
    Fluoride,
    #[default]
    Chloride,
    Bromide,
    Iodide,
}
impl NegativeIon {
    pub(crate) fn symbol(self) -> &'static str {
        match self {
            Self::Fluoride => "F",
            Self::Chloride => "Cl",
            Self::Bromide => "Br",
            Self::Iodide => "I",
        }
    }
}

/// Complete exclusion radii bound to one exact topology, in its dense atom order.
/// Zero radii are allowed, for example for atoms without a Lennard-Jones term.
#[derive(Debug, Clone)]
pub struct SolvationRadii {
    topology: Arc<Topology>,
    radii: Vec<f64>,
}
impl SolvationRadii {
    pub fn new(topology: Arc<Topology>, radii: Quantity<Vec<f64>>) -> Result<Self, SolvationError> {
        let radii = radii.into_unit(CANONICAL_LENGTH_UNIT)?.into_value();
        if radii.len() != topology.atom_count() {
            return Err(SolvationError::RadiusCountMismatch);
        }
        if radii.iter().any(|r| !r.is_finite() || *r < 0.0) {
            return Err(SolvationError::InvalidRadius);
        }
        Ok(Self { topology, radii })
    }
}

/// Solvation controls. All concentration and distance inputs accept compatible units.
///
/// Existing waters and ions are retained. Salt is additional, not a target inferred
/// from existing ions. No protonation or missing-atom reconstruction is performed.
/// For incomplete formal charges, supply `net_charge` explicitly. Radii default to
/// the element's bonded/reference van der Waals column (including its documented
/// fallback), not force-field radii; this is not OpenMM parameterization parity.
#[derive(Debug, Clone)]
pub struct SolventOptions {
    pub water_model: WaterModel,
    pub neutralize: bool,
    pub ionic_strength: Quantity<f64>,
    pub positive_ion: PositiveIon,
    pub negative_ion: NegativeIon,
    /// Stable deterministic ion-site shuffle. Identical inputs and seeds reproduce placement.
    pub seed: u64,
    /// Minimum periodic distance between added ions and other added/existing monatomic ions.
    pub ion_separation: Quantity<f64>,
    pub solute_radii: Option<SolvationRadii>,
    /// Integer charge in elementary-charge units; does not rewrite represented atom charges.
    pub net_charge: Option<i64>,
    /// Bound on tiled candidate waters examined, before cropping and clash removal.
    pub max_candidates: usize,
}
impl Default for SolventOptions {
    fn default() -> Self {
        Self {
            water_model: WaterModel::Tip3p,
            neutralize: true,
            ionic_strength: Quantity::new(0.0, MOLAR),
            positive_ion: PositiveIon::Sodium,
            negative_ion: NegativeIon::Chloride,
            seed: 0,
            ion_separation: Quantity::new(0.5, NANOMETER),
            solute_radii: None,
            net_charge: None,
            max_candidates: 10_000_000,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolvationChargeBasis {
    FormalCharges,
    Override,
}
/// Optional diagnostics from a successful operation. It is fine to discard this value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SolvationReport {
    pub waters_added: usize,
    pub positive_ions_added: usize,
    pub negative_ions_added: usize,
    /// Charge used for ion counting, in elementary-charge units.
    pub solute_charge: i64,
    pub charge_basis: SolvationChargeBasis,
    pub seed: u64,
    pub cleared_topology_properties: Vec<PropertyKey>,
    pub cleared_model_properties: Vec<PropertyKey>,
}

#[derive(Debug)]
#[non_exhaustive]
pub enum SolvationError {
    MissingCell,
    NotFullyPeriodic,
    RadiusTopologyMismatch,
    RadiusCountMismatch,
    InvalidRadius,
    MissingElementRadius(InstanceAtomId),
    InvalidIonicStrength,
    InvalidIonSeparation,
    CandidateLimit,
    CapacityOverflow,
    NoSolventSpace,
    InsufficientIonSites { requested: usize, placed: usize },
    Unit(UnitError),
    Geometry(PeriodicGeometryError),
    Position(PositionError),
    Build(ModelBuildError),
    Hierarchy(HierarchyError),
}
impl fmt::Display for SolvationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingCell => f.write_str("solvation requires a periodic cell"),
            Self::NotFullyPeriodic => f.write_str("solvation requires three periodic axes"),
            Self::RadiusTopologyMismatch => {
                f.write_str("solute radii belong to a different topology snapshot")
            }
            Self::RadiusCountMismatch => {
                f.write_str("provide exactly one radius per topology atom")
            }
            Self::InvalidRadius => f.write_str("radii must be finite and nonnegative"),
            Self::MissingElementRadius(atom) => {
                write!(f, "no elemental radius for {atom}; supply explicit radii")
            }
            Self::InvalidIonicStrength => {
                f.write_str("ionic strength must be finite and nonnegative")
            }
            Self::InvalidIonSeparation => f.write_str("ion separation must be finite and positive"),
            Self::CandidateLimit => {
                f.write_str("solvent tiling exceeds the configured candidate limit")
            }
            Self::CapacityOverflow => f.write_str("solvation count or allocation exceeds capacity"),
            Self::NoSolventSpace => f.write_str("no nonclashing water sites fit in this cell"),
            Self::InsufficientIonSites { requested, placed } => write!(
                f,
                "cannot place {requested} ions; found {placed} eligible sites"
            ),
            Self::Unit(e) => e.fmt(f),
            Self::Geometry(e) => e.fmt(f),
            Self::Position(e) => e.fmt(f),
            Self::Build(e) => e.fmt(f),
            Self::Hierarchy(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for SolvationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(e) => Some(e),
            Self::Geometry(e) => Some(e),
            Self::Position(e) => Some(e),
            Self::Build(e) => Some(e),
            Self::Hierarchy(e) => Some(e),
            _ => None,
        }
    }
}
macro_rules! from_error {
    ($source:ty, $variant:ident) => {
        impl From<$source> for SolvationError {
            fn from(error: $source) -> Self {
                Self::$variant(error)
            }
        }
    };
}
from_error!(UnitError, Unit);
from_error!(PeriodicGeometryError, Geometry);
from_error!(PositionError, Position);
from_error!(ModelBuildError, Build);
from_error!(HierarchyError, Hierarchy);

impl Model {
    /// Fills the existing cell with TIP3P water and optional monovalent ions.
    ///
    /// Existing positions and cell remain unchanged. Only newly generated waters
    /// may be replaced by ions. Water oxygen sites are tested against solute atoms
    /// using the sum of the solute radius and TIP3P oxygen exclusion radius.
    /// Periodic water-water contacts below that oxygen radius are removed without
    /// changing the bundled water geometry. Thus boundary cropping may reduce density.
    ///
    /// Neutralizing counterions are counted first. Additional salt pairs are
    /// `floor((water_sites - counterions) * ionic_strength / 55.4 M + 0.5)`.
    /// Ion sites use a deterministic shuffled greedy search; infeasible placement
    /// returns an error rather than publishing a partially solvated model.
    ///
    /// Publication is atomic. Existing entity IDs/order and annotations survive;
    /// topology/model owner annotations are cleared as for an ordinary append.
    /// The new topology invalidates reuse of old snapshot-bound selections and
    /// parameterizations. The returned report is optional to inspect.
    ///
    /// Extremely large coordinate offsets that cannot retain generated positions
    /// within 1e-8 nm are rejected. This preserves the water's internal geometry.
    ///
    /// ```
    /// use kekule::{smiles, structure::{Model, Positions, PeriodicBoxOptions,
    ///     BoxShape, SolventOptions}, units::{Quantity, NANOMETER, MOLAR}};
    /// let mut model = Model::from_molecule(&smiles::to_molecules("[Na+]")?.remove(0),
    ///     &Positions::zeros(1))?;
    /// model.add_periodic_box(&PeriodicBoxOptions::Padding {
    ///     padding: Quantity::new(1.0, NANOMETER), shape: BoxShape::Cube,
    /// })?;
    /// model.add_solvent(&SolventOptions {
    ///     ionic_strength: Quantity::new(0.15, MOLAR), ..Default::default()
    /// })?; // The successful report can simply be discarded.
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn add_solvent(
        &mut self,
        options: &SolventOptions,
    ) -> Result<SolvationReport, SolvationError> {
        packing::solvate(self, options)
    }
}

/// Numerically safer bounding-box midpoint than `(min + max) / 2`.
fn center(points: &[Point3]) -> Result<Point3, PeriodicGeometryError> {
    let mut low = points[0];
    let mut high = points[0];
    for p in &points[1..] {
        low = Point3::new(low.x.min(p.x), low.y.min(p.y), low.z.min(p.z));
        high = Point3::new(high.x.max(p.x), high.y.max(p.y), high.z.max(p.z));
    }
    let result = Point3::new(
        low.x * 0.5 + high.x * 0.5,
        low.y * 0.5 + high.y * 0.5,
        low.z * 0.5 + high.z * 0.5,
    );
    if !result.is_finite() {
        return Err(PeriodicGeometryError::NumericalFailure);
    }
    Ok(result)
}
