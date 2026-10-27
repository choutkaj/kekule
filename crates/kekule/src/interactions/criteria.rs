//! The geometric criteria of each interaction type.

use std::fmt;

use crate::units::{Quantity, Unit, ANGSTROM, CANONICAL_ANGLE_UNIT, CANONICAL_LENGTH_UNIT, DEGREE};

use super::InteractionError;

/// Generates a getter and a `with_*` setter for a [`Quantity`] criterion.
macro_rules! quantity_criterion {
    ($(#[$doc:meta])* $field:ident, $setter:ident) => {
        $(#[$doc])*
        pub const fn $field(&self) -> Quantity<f64> {
            self.$field
        }

        #[doc = concat!("Sets [`Self::", stringify!($field), "`]; any unit of its dimension.")]
        #[must_use]
        pub fn $setter(mut self, $field: Quantity<f64>) -> Self {
            self.$field = $field;
            self
        }
    };
}

/// Generates the `enabled` getter and setter of one interaction type.
macro_rules! enabled_criterion {
    () => {
        /// Whether interactions of this type are searched; default `true`.
        pub const fn enabled(&self) -> bool {
            self.enabled
        }

        /// Sets [`Self::enabled`].
        #[must_use]
        pub fn with_enabled(mut self, enabled: bool) -> Self {
            self.enabled = enabled;
            self
        }
    };
}

/// One adjustable criterion, named by
/// [`InteractionError::InvalidCriterion`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum InteractionCriterion {
    /// [`HydrogenBondCriteria::distance`].
    HydrogenBondDistance,
    /// [`HydrogenBondCriteria::hydrogen_distance`].
    HydrogenBondHydrogenDistance,
    /// [`HydrogenBondCriteria::hydrogen_angle`].
    HydrogenBondHydrogenAngle,
    /// [`HydrogenBondCriteria::donor_angle`].
    HydrogenBondDonorAngle,
    /// [`HydrogenBondCriteria::acceptor_angle`].
    HydrogenBondAcceptorAngle,
    /// [`SaltBridgeCriteria::distance`].
    SaltBridgeDistance,
    /// [`PiStackingCriteria::distance`].
    PiStackingDistance,
    /// [`PiStackingCriteria::offset`].
    PiStackingOffset,
    /// [`PiStackingCriteria::parallel_angle`].
    PiStackingParallelAngle,
    /// [`PiStackingCriteria::perpendicular_angle`].
    PiStackingPerpendicularAngle,
    /// [`CationPiCriteria::distance`].
    CationPiDistance,
    /// [`CationPiCriteria::offset`].
    CationPiOffset,
    /// [`HalogenBondCriteria::distance`].
    HalogenBondDistance,
    /// [`HalogenBondCriteria::donor_angle`].
    HalogenBondDonorAngle,
    /// [`HalogenBondCriteria::min_acceptor_angle`].
    HalogenBondMinAcceptorAngle,
    /// [`HalogenBondCriteria::max_acceptor_angle`].
    HalogenBondMaxAcceptorAngle,
    /// [`MetalCoordinationCriteria::distance`].
    MetalCoordinationDistance,
    /// [`HydrophobicCriteria::distance`].
    HydrophobicDistance,
}

impl fmt::Display for InteractionCriterion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::HydrogenBondDistance => "hydrogen-bond donor-acceptor distance",
            Self::HydrogenBondHydrogenDistance => "hydrogen-bond hydrogen-acceptor distance",
            Self::HydrogenBondHydrogenAngle => "hydrogen-bond donor-hydrogen-acceptor angle",
            Self::HydrogenBondDonorAngle => "hydrogen-bond donor angle",
            Self::HydrogenBondAcceptorAngle => "hydrogen-bond acceptor angle",
            Self::SaltBridgeDistance => "salt-bridge distance",
            Self::PiStackingDistance => "pi-stacking distance",
            Self::PiStackingOffset => "pi-stacking offset",
            Self::PiStackingParallelAngle => "pi-stacking parallel angle",
            Self::PiStackingPerpendicularAngle => "pi-stacking perpendicular angle",
            Self::CationPiDistance => "cation-pi distance",
            Self::CationPiOffset => "cation-pi offset",
            Self::HalogenBondDistance => "halogen-bond distance",
            Self::HalogenBondDonorAngle => "halogen-bond donor angle",
            Self::HalogenBondMinAcceptorAngle => "halogen-bond minimum acceptor angle",
            Self::HalogenBondMaxAcceptorAngle => "halogen-bond maximum acceptor angle",
            Self::MetalCoordinationDistance => "metal-coordination distance",
            Self::HydrophobicDistance => "hydrophobic-contact distance",
        })
    }
}

/// A length criterion in the canonical unit (nm), finite and nonnegative.
fn length(value: Quantity<f64>, criterion: InteractionCriterion) -> Result<f64, InteractionError> {
    canonical(value, CANONICAL_LENGTH_UNIT, criterion, f64::INFINITY)
}

/// An angle criterion in the canonical unit (rad), between 0 and π.
fn angle(value: Quantity<f64>, criterion: InteractionCriterion) -> Result<f64, InteractionError> {
    canonical(value, CANONICAL_ANGLE_UNIT, criterion, std::f64::consts::PI)
}

fn canonical(
    value: Quantity<f64>,
    unit: Unit,
    criterion: InteractionCriterion,
    max: f64,
) -> Result<f64, InteractionError> {
    match value.value_in(unit) {
        Ok(value) if value.is_finite() && (0.0..=max).contains(&value) => Ok(value),
        _ => Err(InteractionError::InvalidCriterion(criterion)),
    }
}

fn angstrom(value: f64) -> Quantity<f64> {
    Quantity::new(value, ANGSTROM)
}

fn degrees(value: f64) -> Quantity<f64> {
    Quantity::new(value, DEGREE)
}

/// Hydrogen-bond criteria.
///
/// A donor is a nitrogen or oxygen atom bearing a hydrogen; an acceptor is
/// an oxygen without positive charge, or a nitrogen without hydrogens or
/// positive charge whose lone pair is not delocalized, as an amide's,
/// sulfonamide's, or pyrrole's is (see [`super::InteractionFeatures`]). A
/// pair qualifies when the donor–acceptor distance is at most
/// [`Self::distance`], and every heavy neighbor Y of the acceptor makes a
/// D…A–Y angle of at least [`Self::acceptor_angle`]. Bonds to metals do not
/// count as neighbors, here or at the donor.
///
/// When all of the donor's hydrogens are explicit graph atoms, one of them
/// must lie within [`Self::hydrogen_distance`] of the acceptor with a D–H…A
/// angle of at least [`Self::hydrogen_angle`]; the hydrogen with the widest
/// such angle is reported. Otherwise, as in most crystal structures, the
/// hydrogens' positions are unknown, and every heavy neighbor X of the donor
/// must make an X–D…A angle of at least [`Self::donor_angle`].
///
/// Defaults: a donor–acceptor distance of 3.5 Å, a hydrogen–acceptor
/// distance of 2.5 Å and a D–H…A angle of at least 120° (Baker & Hubbard,
/// *Prog. Biophys. Mol. Biol.* 44, 97–179, 1984), and X–D…A and D…A–Y angles
/// of at least 90° (HBPLUS: McDonald & Thornton, *J. Mol. Biol.* 238,
/// 777–793, 1994).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HydrogenBondCriteria {
    enabled: bool,
    distance: Quantity<f64>,
    hydrogen_distance: Quantity<f64>,
    hydrogen_angle: Quantity<f64>,
    donor_angle: Quantity<f64>,
    acceptor_angle: Quantity<f64>,
}

impl Default for HydrogenBondCriteria {
    fn default() -> Self {
        Self {
            enabled: true,
            distance: angstrom(3.5),
            hydrogen_distance: angstrom(2.5),
            hydrogen_angle: degrees(120.0),
            donor_angle: degrees(90.0),
            acceptor_angle: degrees(90.0),
        }
    }
}

impl HydrogenBondCriteria {
    enabled_criterion!();
    quantity_criterion!(
        /// The longest donor–acceptor distance; default 3.5 Å.
        distance,
        with_distance
    );
    quantity_criterion!(
        /// The longest hydrogen–acceptor distance, with explicit hydrogens;
        /// default 2.5 Å.
        hydrogen_distance,
        with_hydrogen_distance
    );
    quantity_criterion!(
        /// The smallest D–H…A angle, with explicit hydrogens; default 120°.
        hydrogen_angle,
        with_hydrogen_angle
    );
    quantity_criterion!(
        /// The smallest X–D…A angle at the donor's heavy neighbors, without
        /// explicit hydrogens; default 90°.
        donor_angle,
        with_donor_angle
    );
    quantity_criterion!(
        /// The smallest D…A–Y angle at the acceptor's heavy neighbors;
        /// default 90°.
        acceptor_angle,
        with_acceptor_angle
    );

    pub(super) fn canonical(&self) -> Result<HydrogenBond, InteractionError> {
        use InteractionCriterion as C;
        Ok(HydrogenBond {
            enabled: self.enabled,
            distance: length(self.distance, C::HydrogenBondDistance)?,
            hydrogen_distance: length(self.hydrogen_distance, C::HydrogenBondHydrogenDistance)?,
            hydrogen_angle: angle(self.hydrogen_angle, C::HydrogenBondHydrogenAngle)?,
            donor_angle: angle(self.donor_angle, C::HydrogenBondDonorAngle)?,
            acceptor_angle: angle(self.acceptor_angle, C::HydrogenBondAcceptorAngle)?,
        })
    }
}

/// Salt-bridge criteria: a negatively and a positively charged group whose
/// charge centers, the centroids of their charged atoms (such as a
/// carboxylate's oxygens or a guanidinium's nitrogens), lie within
/// [`Self::distance`].
///
/// Default 4.0 Å (Kumar & Nussinov, *ChemBioChem* 3, 604–617, 2002).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SaltBridgeCriteria {
    enabled: bool,
    distance: Quantity<f64>,
}

impl Default for SaltBridgeCriteria {
    fn default() -> Self {
        Self {
            enabled: true,
            distance: angstrom(4.0),
        }
    }
}

impl SaltBridgeCriteria {
    enabled_criterion!();
    quantity_criterion!(
        /// The longest distance between charge centers; default 4.0 Å.
        distance,
        with_distance
    );
}

/// π-stacking criteria for two aromatic rings whose centers lie within
/// [`Self::distance`] and whose offset, the distance from one ring's center
/// to the other center's projection onto that ring's plane (the smaller of
/// the two), is at most [`Self::offset`]. The rings stack in parallel when
/// the angle between their planes is at most [`Self::parallel_angle`] and
/// T-shaped when it is at least [`Self::perpendicular_angle`].
///
/// Defaults: 5.5 Å, 2.0 Å, 30°, and 60°, as PLIP's (Salentin et al.,
/// *Nucleic Acids Res.* 43, W443–W447, 2015).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PiStackingCriteria {
    enabled: bool,
    distance: Quantity<f64>,
    offset: Quantity<f64>,
    parallel_angle: Quantity<f64>,
    perpendicular_angle: Quantity<f64>,
}

impl Default for PiStackingCriteria {
    fn default() -> Self {
        Self {
            enabled: true,
            distance: angstrom(5.5),
            offset: angstrom(2.0),
            parallel_angle: degrees(30.0),
            perpendicular_angle: degrees(60.0),
        }
    }
}

impl PiStackingCriteria {
    enabled_criterion!();
    quantity_criterion!(
        /// The longest distance between ring centers; default 5.5 Å.
        distance,
        with_distance
    );
    quantity_criterion!(
        /// The largest offset; default 2.0 Å.
        offset,
        with_offset
    );
    quantity_criterion!(
        /// The largest angle between the planes of parallel rings; default
        /// 30°.
        parallel_angle,
        with_parallel_angle
    );
    quantity_criterion!(
        /// The smallest angle between the planes of T-shaped rings; default
        /// 60°.
        perpendicular_angle,
        with_perpendicular_angle
    );
}

/// Cation–π criteria: a positively charged group whose charge center lies
/// within [`Self::distance`] of an aromatic ring's center and over the ring,
/// within [`Self::offset`] of the ring's axis.
///
/// Defaults: 6.0 Å (Gallivan & Dougherty, *PNAS* 96, 9459–9464, 1999; PLIP)
/// and 2.0 Å (PLIP: Salentin et al. 2015).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CationPiCriteria {
    enabled: bool,
    distance: Quantity<f64>,
    offset: Quantity<f64>,
}

impl Default for CationPiCriteria {
    fn default() -> Self {
        Self {
            enabled: true,
            distance: angstrom(6.0),
            offset: angstrom(2.0),
        }
    }
}

impl CationPiCriteria {
    enabled_criterion!();
    quantity_criterion!(
        /// The longest distance from the ring center to the charge center;
        /// default 6.0 Å.
        distance,
        with_distance
    );
    quantity_criterion!(
        /// The largest distance of the charge center from the ring's axis;
        /// default 2.0 Å.
        offset,
        with_offset
    );
}

/// Halogen-bond criteria: a chlorine, bromine, or iodine bonded to carbon
/// whose distance to a hydrogen-bond acceptor is at most [`Self::distance`],
/// with a C–X…A angle of at least [`Self::donor_angle`] and, at every heavy
/// neighbor Y of the acceptor, an X…A–Y angle between
/// [`Self::min_acceptor_angle`] and [`Self::max_acceptor_angle`].
///
/// Defaults: 4.0 Å, 135°, and 90°–150°, PLIP's 165° ± 30° and 120° ± 30°
/// (Salentin et al. 2015, after Auffinger et al., *PNAS* 101, 16789–16794,
/// 2004).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HalogenBondCriteria {
    enabled: bool,
    distance: Quantity<f64>,
    donor_angle: Quantity<f64>,
    min_acceptor_angle: Quantity<f64>,
    max_acceptor_angle: Quantity<f64>,
}

impl Default for HalogenBondCriteria {
    fn default() -> Self {
        Self {
            enabled: true,
            distance: angstrom(4.0),
            donor_angle: degrees(135.0),
            min_acceptor_angle: degrees(90.0),
            max_acceptor_angle: degrees(150.0),
        }
    }
}

impl HalogenBondCriteria {
    enabled_criterion!();
    quantity_criterion!(
        /// The longest halogen–acceptor distance; default 4.0 Å.
        distance,
        with_distance
    );
    quantity_criterion!(
        /// The smallest C–X…A angle; default 135°.
        donor_angle,
        with_donor_angle
    );
    quantity_criterion!(
        /// The smallest X…A–Y angle; default 90°.
        min_acceptor_angle,
        with_min_acceptor_angle
    );
    quantity_criterion!(
        /// The largest X…A–Y angle; default 150°.
        max_acceptor_angle,
        with_max_acceptor_angle
    );
}

/// Metal-coordination criteria: a metal atom within [`Self::distance`] of a
/// nitrogen, oxygen, sulfur, or selenium atom, or of a halide ion. Bonds
/// between them, as mmCIF metal links give, do not exclude the pair.
///
/// Default 3.0 Å (PLIP: Salentin et al. 2015; Harding, *Acta Cryst.* D62,
/// 678–682, 2006, reports typical distances of 1.9–2.6 Å).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MetalCoordinationCriteria {
    enabled: bool,
    distance: Quantity<f64>,
}

impl Default for MetalCoordinationCriteria {
    fn default() -> Self {
        Self {
            enabled: true,
            distance: angstrom(3.0),
        }
    }
}

impl MetalCoordinationCriteria {
    enabled_criterion!();
    quantity_criterion!(
        /// The longest metal–ligand distance; default 3.0 Å.
        distance,
        with_distance
    );
}

/// Hydrophobic-contact criteria: two hydrophobic atoms, carbons bonded only
/// to carbon and hydrogen, within [`Self::distance`].
///
/// Default 4.0 Å, as PLIP's (Salentin et al. 2015).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HydrophobicCriteria {
    enabled: bool,
    distance: Quantity<f64>,
}

impl Default for HydrophobicCriteria {
    fn default() -> Self {
        Self {
            enabled: true,
            distance: angstrom(4.0),
        }
    }
}

impl HydrophobicCriteria {
    enabled_criterion!();
    quantity_criterion!(
        /// The longest distance between the atoms; default 4.0 Å.
        distance,
        with_distance
    );
}

/// The criteria of every interaction type, and which types are searched.
///
/// A plain value: every criterion is a [`Quantity`] in any unit of its
/// dimension, read by its getter and replaced by its `with_*` setter, so a
/// host can store and edit them. Values are checked when used: lengths must
/// be finite and nonnegative and angles between 0 and 180°, or the search
/// fails with [`InteractionError::InvalidCriterion`]; [`Self::validate`]
/// checks them beforehand.
///
/// Within one structure, two ends in one molecule fewer than
/// [`Self::min_bond_separation`] bonds apart, such as the atoms of one
/// residue's carboxylate and amine, never interact, except in metal
/// coordination; paths through metal atoms do not count. The default, 4,
/// keeps 1–5 pairs and beyond.
///
/// ```
/// use kekule::interactions::InteractionCriteria;
/// use kekule::units::{Quantity, ANGSTROM};
///
/// let criteria = InteractionCriteria::default();
/// let wider = criteria.hydrogen_bond().with_distance(Quantity::new(3.9, ANGSTROM));
/// let criteria = criteria
///     .with_hydrogen_bond(wider)
///     .with_hydrophobic(criteria.hydrophobic().with_enabled(false));
/// assert!(criteria.validate().is_ok());
/// assert!(!criteria.hydrophobic().enabled());
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InteractionCriteria {
    hydrogen_bond: HydrogenBondCriteria,
    salt_bridge: SaltBridgeCriteria,
    pi_stacking: PiStackingCriteria,
    cation_pi: CationPiCriteria,
    halogen_bond: HalogenBondCriteria,
    metal_coordination: MetalCoordinationCriteria,
    hydrophobic: HydrophobicCriteria,
    min_bond_separation: usize,
}

impl Default for InteractionCriteria {
    fn default() -> Self {
        Self {
            hydrogen_bond: HydrogenBondCriteria::default(),
            salt_bridge: SaltBridgeCriteria::default(),
            pi_stacking: PiStackingCriteria::default(),
            cation_pi: CationPiCriteria::default(),
            halogen_bond: HalogenBondCriteria::default(),
            metal_coordination: MetalCoordinationCriteria::default(),
            hydrophobic: HydrophobicCriteria::default(),
            min_bond_separation: Self::DEFAULT_MIN_BOND_SEPARATION,
        }
    }
}

/// Generates the getter and setter of one type's criteria.
macro_rules! type_criteria {
    ($field:ident, $setter:ident, $type:ty) => {
        #[doc = concat!("The [`", stringify!($type), "`].")]
        pub const fn $field(&self) -> $type {
            self.$field
        }

        #[doc = concat!("Sets [`Self::", stringify!($field), "`].")]
        #[must_use]
        pub fn $setter(mut self, $field: $type) -> Self {
            self.$field = $field;
            self
        }
    };
}

impl InteractionCriteria {
    /// The [`Default`] of [`Self::min_bond_separation`].
    pub const DEFAULT_MIN_BOND_SEPARATION: usize = 4;

    type_criteria!(hydrogen_bond, with_hydrogen_bond, HydrogenBondCriteria);
    type_criteria!(salt_bridge, with_salt_bridge, SaltBridgeCriteria);
    type_criteria!(pi_stacking, with_pi_stacking, PiStackingCriteria);
    type_criteria!(cation_pi, with_cation_pi, CationPiCriteria);
    type_criteria!(halogen_bond, with_halogen_bond, HalogenBondCriteria);
    type_criteria!(
        metal_coordination,
        with_metal_coordination,
        MetalCoordinationCriteria
    );
    type_criteria!(hydrophobic, with_hydrophobic, HydrophobicCriteria);

    /// The fewest bonds between two interacting ends in one molecule of one
    /// structure; default [`Self::DEFAULT_MIN_BOND_SEPARATION`]. Zero and
    /// one exclude only overlapping ends.
    pub const fn min_bond_separation(&self) -> usize {
        self.min_bond_separation
    }

    /// Sets [`Self::min_bond_separation`].
    #[must_use]
    pub fn with_min_bond_separation(mut self, min_bond_separation: usize) -> Self {
        self.min_bond_separation = min_bond_separation;
        self
    }

    /// Checks every criterion, enabled or not.
    ///
    /// # Errors
    ///
    /// [`InteractionError::InvalidCriterion`] for the first length that is
    /// not a finite, nonnegative length, or angle that is not an angle
    /// between 0 and 180°.
    pub fn validate(&self) -> Result<(), InteractionError> {
        self.canonical().map(|_| ())
    }

    pub(super) fn canonical(&self) -> Result<Canonical, InteractionError> {
        use InteractionCriterion as C;
        let (pi, cation, halogen) = (self.pi_stacking, self.cation_pi, self.halogen_bond);
        Ok(Canonical {
            hydrogen_bond: self.hydrogen_bond.canonical()?,
            salt_bridge: Distance {
                enabled: self.salt_bridge.enabled,
                distance: length(self.salt_bridge.distance, C::SaltBridgeDistance)?,
            },
            pi_stacking: PiStacking {
                enabled: pi.enabled,
                distance: length(pi.distance, C::PiStackingDistance)?,
                offset: length(pi.offset, C::PiStackingOffset)?,
                parallel_angle: angle(pi.parallel_angle, C::PiStackingParallelAngle)?,
                perpendicular_angle: angle(
                    pi.perpendicular_angle,
                    C::PiStackingPerpendicularAngle,
                )?,
            },
            cation_pi: CationPi {
                enabled: cation.enabled,
                distance: length(cation.distance, C::CationPiDistance)?,
                offset: length(cation.offset, C::CationPiOffset)?,
            },
            halogen_bond: HalogenBond {
                enabled: halogen.enabled,
                distance: length(halogen.distance, C::HalogenBondDistance)?,
                donor_angle: angle(halogen.donor_angle, C::HalogenBondDonorAngle)?,
                min_acceptor_angle: angle(
                    halogen.min_acceptor_angle,
                    C::HalogenBondMinAcceptorAngle,
                )?,
                max_acceptor_angle: angle(
                    halogen.max_acceptor_angle,
                    C::HalogenBondMaxAcceptorAngle,
                )?,
            },
            metal_coordination: Distance {
                enabled: self.metal_coordination.enabled,
                distance: length(
                    self.metal_coordination.distance,
                    C::MetalCoordinationDistance,
                )?,
            },
            hydrophobic: Distance {
                enabled: self.hydrophobic.enabled,
                distance: length(self.hydrophobic.distance, C::HydrophobicDistance)?,
            },
            min_bond_separation: self.min_bond_separation,
        })
    }
}

/// The criteria in canonical units (nm and rad).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Canonical {
    pub(super) hydrogen_bond: HydrogenBond,
    pub(super) salt_bridge: Distance,
    pub(super) pi_stacking: PiStacking,
    pub(super) cation_pi: CationPi,
    pub(super) halogen_bond: HalogenBond,
    pub(super) metal_coordination: Distance,
    pub(super) hydrophobic: Distance,
    pub(super) min_bond_separation: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct HydrogenBond {
    pub(super) enabled: bool,
    pub(super) distance: f64,
    pub(super) hydrogen_distance: f64,
    pub(super) hydrogen_angle: f64,
    pub(super) donor_angle: f64,
    pub(super) acceptor_angle: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Distance {
    pub(super) enabled: bool,
    pub(super) distance: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct PiStacking {
    pub(super) enabled: bool,
    pub(super) distance: f64,
    pub(super) offset: f64,
    pub(super) parallel_angle: f64,
    pub(super) perpendicular_angle: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct CationPi {
    pub(super) enabled: bool,
    pub(super) distance: f64,
    pub(super) offset: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct HalogenBond {
    pub(super) enabled: bool,
    pub(super) distance: f64,
    pub(super) donor_angle: f64,
    pub(super) min_acceptor_angle: f64,
    pub(super) max_acceptor_angle: f64,
}
