mod boxes;
mod neighbors;
mod packing;
pub use boxes::{add_periodic_box, BoxShape, PeriodicBoxError, PeriodicBoxOptions};

use crate::diagnostics::{LimitExceeded, Resource};
use crate::structure::{DenseArrayError, Model, ModelBuildError};
use crate::{
    geometry::{PeriodicGeometryError, Point3},
    properties::PropertyKey,
    topology::{HierarchyError, InstanceAtomId, Topology},
    units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT, MOLAR, NANOMETER},
};
use std::{fmt, sync::Arc};

/// Built-in coordinate model. This does not install force-field parameters.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum WaterModel {
    /// Three-site TIP3P water geometry, bundled from OpenMM's equilibrated box.
    #[default]
    Tip3p,
}
/// Monovalent cations supported by the solvent builder.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum PositiveIon {
    /// Li⁺.
    Lithium,
    /// Na⁺, the default.
    #[default]
    Sodium,
    /// K⁺.
    Potassium,
    /// Rb⁺.
    Rubidium,
    /// Cs⁺.
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
#[non_exhaustive]
pub enum NegativeIon {
    /// F⁻.
    Fluoride,
    /// Cl⁻, the default.
    #[default]
    Chloride,
    /// Br⁻.
    Bromide,
    /// I⁻.
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

/// Complete exclusion radii bound to one topology layout, in its dense atom order.
/// Zero radii are allowed, for example for atoms without a Lennard-Jones term.
#[derive(Debug, Clone)]
pub struct SolvationRadii {
    topology: Arc<Topology>,
    radii: Vec<f64>,
}
impl SolvationRadii {
    /// Binds one radius per atom, in `topology`'s dense atom order, accepting
    /// any length unit.
    ///
    /// # Errors
    ///
    /// Returns [`SolvationError::Unit`] when the unit is not a length,
    /// [`SolvationError::RadiusCountMismatch`] when the count differs from the
    /// atom count, and [`SolvationError::InvalidRadius`] for a negative or
    /// non-finite radius.
    pub fn new(
        topology: impl Into<Arc<Topology>>,
        radii: Quantity<Vec<f64>>,
    ) -> Result<Self, SolvationError> {
        let topology = topology.into();
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
    water_model: WaterModel,
    neutralize: bool,
    ionic_strength: Quantity<f64>,
    positive_ion: PositiveIon,
    negative_ion: NegativeIon,
    seed: u64,
    ion_separation: Quantity<f64>,
    solute_radii: Option<SolvationRadii>,
    net_charge: Option<i64>,
    max_candidates: usize,
}

impl SolventOptions {
    /// The water geometry to add; default [`WaterModel::Tip3p`].
    pub const fn water_model(&self) -> WaterModel {
        self.water_model
    }

    /// Sets [`Self::water_model`].
    #[must_use]
    pub fn with_water_model(mut self, water_model: WaterModel) -> Self {
        self.water_model = water_model;
        self
    }

    /// Whether counterions are added to cancel the solute charge; default `true`.
    pub const fn neutralize(&self) -> bool {
        self.neutralize
    }

    /// Sets [`Self::neutralize`].
    #[must_use]
    pub fn with_neutralize(mut self, neutralize: bool) -> Self {
        self.neutralize = neutralize;
        self
    }

    /// Additional salt concentration of monovalent ion pairs; default 0 M. Any
    /// concentration unit is accepted; [`solvate`] converts to mol/L.
    pub const fn ionic_strength(&self) -> Quantity<f64> {
        self.ionic_strength
    }

    /// Sets [`Self::ionic_strength`].
    #[must_use]
    pub fn with_ionic_strength(mut self, ionic_strength: Quantity<f64>) -> Self {
        self.ionic_strength = ionic_strength;
        self
    }

    /// The cation used for counterions and salt; default [`PositiveIon::Sodium`].
    pub const fn positive_ion(&self) -> PositiveIon {
        self.positive_ion
    }

    /// Sets [`Self::positive_ion`].
    #[must_use]
    pub fn with_positive_ion(mut self, positive_ion: PositiveIon) -> Self {
        self.positive_ion = positive_ion;
        self
    }

    /// The anion used for counterions and salt; default [`NegativeIon::Chloride`].
    pub const fn negative_ion(&self) -> NegativeIon {
        self.negative_ion
    }

    /// Sets [`Self::negative_ion`].
    #[must_use]
    pub fn with_negative_ion(mut self, negative_ion: NegativeIon) -> Self {
        self.negative_ion = negative_ion;
        self
    }

    /// Stable deterministic ion-site shuffle. Identical inputs and seeds reproduce placement.
    pub const fn seed(&self) -> u64 {
        self.seed
    }

    /// Sets [`Self::seed`].
    #[must_use]
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }

    /// Minimum periodic distance between added ions and other added/existing monatomic ions.
    pub const fn ion_separation(&self) -> Quantity<f64> {
        self.ion_separation
    }

    /// Sets [`Self::ion_separation`].
    #[must_use]
    pub fn with_ion_separation(mut self, ion_separation: Quantity<f64>) -> Self {
        self.ion_separation = ion_separation;
        self
    }

    /// Explicit solute exclusion radii, or `None` to use each element's reference
    /// van der Waals radius.
    pub fn solute_radii(&self) -> Option<&SolvationRadii> {
        self.solute_radii.as_ref()
    }

    /// Sets [`Self::solute_radii`].
    #[must_use]
    pub fn with_solute_radii(mut self, solute_radii: Option<SolvationRadii>) -> Self {
        self.solute_radii = solute_radii;
        self
    }

    /// Integer charge in elementary-charge units; does not rewrite represented atom charges.
    pub const fn net_charge(&self) -> Option<i64> {
        self.net_charge
    }

    /// Sets [`Self::net_charge`].
    #[must_use]
    pub fn with_net_charge(mut self, net_charge: Option<i64>) -> Self {
        self.net_charge = net_charge;
        self
    }

    /// Bound on tiled candidate waters examined, before cropping and clash removal.
    pub const fn max_candidates(&self) -> usize {
        self.max_candidates
    }

    /// Sets [`Self::max_candidates`].
    #[must_use]
    pub fn with_max_candidates(mut self, max_candidates: usize) -> Self {
        self.max_candidates = max_candidates;
        self
    }
}
impl SolventOptions {
    /// The [`Default`] bound of [`Self::max_candidates`].
    pub const DEFAULT_MAX_CANDIDATES: usize = 10_000_000;
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
            max_candidates: Self::DEFAULT_MAX_CANDIDATES,
        }
    }
}

/// Where the solute charge used for ion counting came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SolvationChargeBasis {
    /// The sum of the topology's represented formal charges.
    FormalCharges,
    /// [`SolventOptions::net_charge`].
    Override,
}
/// Optional diagnostics from a successful operation. It is fine to discard this value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SolvationReport {
    pub(crate) waters_added: usize,
    pub(crate) positive_ions_added: usize,
    pub(crate) negative_ions_added: usize,
    pub(crate) solute_charge: i64,
    pub(crate) charge_basis: SolvationChargeBasis,
    pub(crate) seed: u64,
    pub(crate) cleared_topology_properties: Vec<PropertyKey>,
    pub(crate) cleared_model_properties: Vec<PropertyKey>,
}

impl SolvationReport {
    /// Number of water molecules added, after ion replacement.
    pub const fn waters_added(&self) -> usize {
        self.waters_added
    }

    /// Number of cations added, counterions and salt together.
    pub const fn positive_ions_added(&self) -> usize {
        self.positive_ions_added
    }

    /// Number of anions added, counterions and salt together.
    pub const fn negative_ions_added(&self) -> usize {
        self.negative_ions_added
    }

    /// Charge used for ion counting, in elementary-charge units.
    pub const fn solute_charge(&self) -> i64 {
        self.solute_charge
    }

    /// Whether [`Self::solute_charge`] was summed or supplied.
    pub const fn charge_basis(&self) -> SolvationChargeBasis {
        self.charge_basis
    }

    /// The ion-placement seed that was used.
    pub const fn seed(&self) -> u64 {
        self.seed
    }

    /// Destination topology owner keys cleared because the topology changed.
    pub fn cleared_topology_properties(&self) -> &[PropertyKey] {
        &self.cleared_topology_properties
    }

    /// Destination structure owner keys cleared because the model changed.
    pub fn cleared_model_properties(&self) -> &[PropertyKey] {
        &self.cleared_model_properties
    }
}

/// Solvation options or the model prevent filling the cell with solvent.
#[derive(Debug)]
#[non_exhaustive]
pub enum SolvationError {
    /// The model has no periodic cell.
    MissingCell,
    /// The model's cell is not periodic along all three axes.
    NotFullyPeriodic,
    /// The supplied solute radii belong to a different topology layout.
    RadiusTopologyMismatch,
    /// The number of radii differs from the topology's atom count.
    RadiusCountMismatch,
    /// A radius is negative or not finite.
    InvalidRadius,
    /// No reference van der Waals radius is known for this atom's element.
    MissingElementRadius(InstanceAtomId),
    /// The ionic strength is negative or not finite.
    InvalidIonicStrength,
    /// The ion separation is not finite and positive.
    InvalidIonSeparation,
    /// Tiling would examine more candidate waters than
    /// [`SolventOptions::max_candidates`] allows.
    LimitExceeded(LimitExceeded),
    /// An ion count or allocation does not fit in `usize`.
    CapacityOverflow,
    /// The cell is too small for a water, or every water site clashes.
    NoSolventSpace,
    /// Too few water sites satisfy the ion separation for the requested ions.
    InsufficientIonSites {
        /// The number of ions to place.
        requested: usize,
        /// The number of eligible sites found.
        placed: usize,
    },
    /// A quantity could not be converted to the required unit; see the source.
    Unit(UnitError),
    /// A periodic-geometry calculation failed; see the source.
    Geometry(PeriodicGeometryError),
    /// Position data was rejected; see the source.
    Position(DenseArrayError),
    /// Building the solvated model failed; see the source.
    Build(ModelBuildError),
    /// Adding solvent molecules, chains or residues failed; see the source.
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
            Self::LimitExceeded(limit) => write!(f, "solvent tiling {limit}"),
            Self::CapacityOverflow => f.write_str("solvation count or allocation exceeds capacity"),
            Self::NoSolventSpace => f.write_str("no nonclashing water sites fit in this cell"),
            Self::InsufficientIonSites { requested, placed } => write!(
                f,
                "cannot place {requested} ions; found {placed} eligible sites"
            ),
            Self::Unit(_) => f.write_str("invalid quantity unit"),
            Self::Geometry(_) => f.write_str("periodic geometry failed"),
            Self::Position(_) => f.write_str("invalid positions"),
            Self::Build(_) => f.write_str("cannot build the solvated model"),
            Self::Hierarchy(_) => f.write_str("cannot build the solvent hierarchy"),
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
from_error!(DenseArrayError, Position);
from_error!(ModelBuildError, Build);
from_error!(HierarchyError, Hierarchy);

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
/// # Errors
///
/// The model is unchanged on error. The cell must exist
/// ([`SolvationError::MissingCell`]) and be periodic along all three axes
/// ([`SolvationError::NotFullyPeriodic`]). Invalid options return
/// [`SolvationError::InvalidIonicStrength`], [`SolvationError::InvalidIonSeparation`],
/// [`SolvationError::Unit`] or [`SolvationError::RadiusTopologyMismatch`];
/// a solute element without a reference radius returns
/// [`SolvationError::MissingElementRadius`]. A cell with no room for water returns
/// [`SolvationError::NoSolventSpace`], too few separated sites for the ions
/// [`SolvationError::InsufficientIonSites`], and exceeding
/// [`SolventOptions::max_candidates`] [`SolvationError::LimitExceeded`]. The
/// remaining variants report numerical, capacity or publication failures.
///
/// ```
/// use kekule::{smiles, solvation::{self, BoxShape, PeriodicBoxOptions, SolventOptions},
///     structure::{Model, Positions}, units::{Quantity, NANOMETER, MOLAR}};
/// let mut model = Model::from_molecule(smiles::to_molecules("[Na+]")?.remove(0),
///     &Positions::zeros(1))?;
/// solvation::add_periodic_box(&mut model, &PeriodicBoxOptions::Padding {
///     padding: Quantity::new(1.0, NANOMETER), shape: BoxShape::Cube,
/// })?;
/// solvation::solvate(
///     &mut model,
///     &SolventOptions::default().with_ionic_strength(Quantity::new(0.15, MOLAR)),
/// )?; // The successful report can simply be discarded.
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn solvate(
    model: &mut Model,
    options: &SolventOptions,
) -> Result<SolvationReport, SolvationError> {
    packing::solvate(model, options)
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
