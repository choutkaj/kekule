use crate::structure::DenseArrayError;
use std::fmt;
use std::sync::Arc;

use crate::geometry::{PeriodicCell, RigidTransform};
use crate::molecule::Molecule;
use crate::properties::{OwnerProperties, StructureProperties, TopologyPropertiesMut};
use crate::topology::TopologySubsetError;
use crate::topology::{
    AtomSelection, Topology, TopologyBuildError, TopologyBuilder, TopologyPerceptionError,
};
use crate::units::{Quantity, UnitError, CANONICAL_ENERGY_UNIT, CANONICAL_TEMPERATURE_UNIT};

use super::scalars::AtomScalars;
use super::state::{
    is_bound, site_mutators, site_readers, structure_mutators, structure_readers, subset_rows,
    transform_cell, transform_positions, validate_binding,
};
use super::{AsModelView, Model, ModelView, Positions, StructureError, StructureKind};

/// A weighted, unordered sample of structures of one shared topology, such
/// as conformers, Monte Carlo samples, or alternate experimental models.
///
/// Every [`EnsembleMember`] carries a finite positive statistical weight and
/// may carry an energy. Weights are relative: only their ratios are
/// meaningful, so selecting, duplicating, or removing members needs no
/// renormalization. Member order is stable but carries no physical meaning;
/// use [`super::Trajectory`] when order is temporal. Statistics over an
/// ensemble must use its weights. Energies never change weights implicitly;
/// [`Self::assign_boltzmann_weights`] derives them explicitly.
///
/// Members are validated against the shared topology on insertion, so they
/// can be read as [`super::ModelView`]s without further checks. Superposition
/// and RMSD are inherent methods; see [`crate::alignment`].
///
/// ```
/// use kekule::structure::{Ensemble, EnsembleMember, Positions};
/// let topology = kekule::smiles::to_topology("O")?;
/// let atoms = topology.atom_count();
/// let mut ensemble = Ensemble::new(topology);
/// ensemble.push(EnsembleMember::new(Positions::zeros(atoms), 3.0)?)?;
/// ensemble.push(EnsembleMember::new(Positions::zeros(atoms), 1.0)?)?;
/// ensemble.normalize_weights()?;
/// assert_eq!(ensemble.get(0).unwrap().weight(), 0.75);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// Member editors preserve dimensions and cannot overwrite a whole member:
///
/// ```compile_fail,E0594
/// use kekule::structure::{Ensemble, EnsembleMember, Positions};
/// fn overwrite(ensemble: &mut Ensemble) {
///     *ensemble.get_mut(0).unwrap() = EnsembleMember::new(Positions::zeros(1), 1.0).unwrap();
/// }
/// ```
#[derive(Debug, Clone)]
pub struct Ensemble {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    members: Vec<EnsembleMember>,
}

impl Ensemble {
    /// An empty ensemble of one shared topology.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self {
            topology: topology.into(),
            properties: OwnerProperties::new(),
            members: Vec::new(),
        }
    }

    /// Validates and collects each member in order, binding its bond property
    /// rows. Fails on the first invalid member.
    pub fn from_items(
        topology: impl Into<Arc<Topology>>,
        members: impl IntoIterator<Item = EnsembleMember>,
    ) -> Result<Self, EnsembleError> {
        let mut ensemble = Self::new(topology);
        for member in members {
            ensemble.push(member)?;
        }
        Ok(ensemble)
    }

    /// Collects models that share one topology layout as equally weighted
    /// members (weight `1.0`), keeping the first model's snapshot. Positions,
    /// cells, occupancies, B-factors, and properties move without copying.
    pub fn from_models(models: impl IntoIterator<Item = Model>) -> Result<Self, EnsembleError> {
        let mut models = models.into_iter();
        let first = models.next().ok_or(EnsembleError::EmptySource)?;
        let mut ensemble = Self::new(first.shared_topology());
        ensemble.members.push(EnsembleMember::from_model(first));
        for model in models {
            if !model.topology().shares_layout(&ensemble.topology) {
                return Err(EnsembleError::TopologyMismatch);
            }
            ensemble.members.push(EnsembleMember::from_model(model));
        }
        Ok(ensemble)
    }

    /// Builds equally weighted conformers of one molecule from dense
    /// positions in molecule atom order.
    pub fn from_molecule_positions(
        molecule: Molecule,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<Self, EnsembleError> {
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(molecule)?;
        builder.add_instance(definition)?;
        Self::from_items(
            builder.build()?,
            positions.into_iter().map(EnsembleMember::unit),
        )
    }

    /// Members already bound to `topology`, for projections from another
    /// collection of the same topology.
    pub(crate) fn from_bound_members(
        topology: Arc<Topology>,
        properties: OwnerProperties,
        members: Vec<EnsembleMember>,
    ) -> Self {
        debug_assert!(members.iter().all(|member| is_bound(
            member.atom_count(),
            &member.properties,
            &topology
        )));
        Self {
            topology,
            properties,
            members,
        }
    }

    /// The topology shared by every member.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// A shared handle to the topology; cloning the [`Arc`] copies no topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Ensemble-level annotations.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    /// Mutable ensemble-level annotations. Member properties are edited through
    /// [`Self::get_mut`].
    pub fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }

    /// Number of members.
    pub fn len(&self) -> usize {
        self.members.len()
    }

    /// Whether the ensemble has no members.
    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    /// One topology-bound member by stable index.
    pub fn get(&self, index: usize) -> Option<EnsembleMemberView<'_>> {
        self.members.get(index).map(|member| EnsembleMemberView {
            topology: &self.topology,
            member,
        })
    }

    /// A dimension-preserving editor for one member. Whole-member replacement
    /// goes through [`Self::replace`].
    pub fn get_mut(&mut self, index: usize) -> Option<EnsembleMemberMut<'_>> {
        self.members
            .get_mut(index)
            .map(|member| EnsembleMemberMut { member })
    }

    /// Topology-bound members in stable order.
    pub fn iter(&self) -> EnsembleIter<'_> {
        EnsembleIter {
            topology: &self.topology,
            members: self.members.iter(),
        }
    }

    /// Validates and appends one member, binding its bond property rows.
    pub fn push(&mut self, mut member: EnsembleMember) -> Result<(), EnsembleError> {
        self.bind(&mut member)?;
        self.members.push(member);
        Ok(())
    }

    /// Validates and replaces one member, returning the previous one. The
    /// index is checked first; failure leaves the ensemble unchanged.
    pub fn replace(
        &mut self,
        index: usize,
        mut member: EnsembleMember,
    ) -> Result<EnsembleMember, EnsembleError> {
        self.check_index(index)?;
        self.bind(&mut member)?;
        Ok(std::mem::replace(&mut self.members[index], member))
    }

    /// Removes and returns one member, shifting later ones down.
    pub fn remove(&mut self, index: usize) -> Result<EnsembleMember, EnsembleError> {
        self.check_index(index)?;
        Ok(self.members.remove(index))
    }

    /// Replaces the positions of every member transactionally, in order.
    /// Other state is kept.
    pub fn replace_positions(
        &mut self,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<(), EnsembleError> {
        let positions = positions.into_iter().collect::<Vec<_>>();
        if positions.len() != self.members.len() {
            return Err(EnsembleError::MemberCountMismatch {
                expected: self.members.len(),
                actual: positions.len(),
            });
        }
        for values in &positions {
            if values.len() != self.topology.atom_count() {
                return Err(StructureError::Array(DenseArrayError::CountMismatch {
                    expected: self.topology.atom_count(),
                    actual: values.len(),
                })
                .into());
            }
        }
        for (member, positions) in self.members.iter_mut().zip(positions) {
            member.positions = positions;
        }
        Ok(())
    }

    /// Copies members in the requested order, sharing the topology. Ranges
    /// and strides can be written as `(start..end).step_by(stride)`;
    /// duplicates and reordering are allowed. Ensemble properties are kept.
    pub fn select(&self, indices: impl IntoIterator<Item = usize>) -> Result<Self, EnsembleError> {
        let members = indices
            .into_iter()
            .map(|index| {
                self.check_index(index)?;
                Ok(self.members[index].clone())
            })
            .collect::<Result<Vec<_>, EnsembleError>>()?;
        Ok(Self {
            topology: self.shared_topology(),
            properties: self.properties.clone(),
            members,
        })
    }

    /// Publishes one induced topology subset and projects every member onto
    /// it. Ensemble and per-member owner properties are dropped; weights and
    /// energies describe the sampled structure and are kept.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, EnsembleError> {
        let subset = crate::topology::Topology::subset(self.topology.clone(), selection)?;
        let (atoms, bonds) = subset_rows(&subset);
        let members = self
            .members
            .iter()
            .map(|member| member.project(&atoms, &bonds))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology: subset.shared_topology(),
            properties: OwnerProperties::new(),
            members,
        })
    }

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once and keeps every member and
    /// the ensemble properties without copying. Failure leaves the ensemble
    /// unchanged. The new snapshot shares the original layout, so selections,
    /// buffers, and prepared potentials stay usable.
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    /// Length-preserving mutable access to the shared topology's static
    /// annotations.
    ///
    /// A topology shared with other owners is first replaced by a snapshot of
    /// it ([`Arc::make_mut`]), so those owners are unaffected. Either way the
    /// layout is kept, so every member, selection, buffer, and prepared
    /// potential remains usable. See [`Topology::properties_mut`].
    pub fn topology_properties_mut(&mut self) -> TopologyPropertiesMut<'_> {
        Arc::make_mut(&mut self.topology).properties_mut()
    }

    /// Consumes the ensemble, transferring topology, ensemble properties, and
    /// every member without copying.
    pub fn into_parts(self) -> (Arc<Topology>, OwnerProperties, Vec<EnsembleMember>) {
        (self.topology, self.properties, self.members)
    }

    /// Consumes the ensemble, returning its members in order and dropping the
    /// topology and ensemble properties.
    pub fn into_items(self) -> Vec<EnsembleMember> {
        self.members
    }

    /// Scales every weight so the weights sum to one, keeping their ratios.
    ///
    /// Fails on an empty ensemble, which has no probabilities, or when a
    /// scaled weight would underflow to zero. Failure changes nothing.
    pub fn normalize_weights(&mut self) -> Result<(), EnsembleError> {
        let probabilities = self.probabilities()?;
        for (member, weight) in self.members.iter_mut().zip(probabilities) {
            member.weight = weight;
        }
        Ok(())
    }

    /// The probability of each member, in member order: its weight divided by
    /// the total weight. The ensemble is not changed.
    ///
    /// Fails on an empty ensemble or when a probability would underflow to
    /// zero.
    pub fn probabilities(&self) -> Result<Vec<f64>, EnsembleError> {
        if self.is_empty() {
            return Err(EnsembleError::EmptySource);
        }
        let scaled = self.scaled_weights();
        let total = scaled.iter().sum::<f64>();
        Ok(scaled
            .into_iter()
            .map(|weight| checked_weight(weight / total))
            .collect::<Result<Vec<_>, _>>()?)
    }

    /// Kish's effective sample size, `(Σw)² / Σw²`: the number of equally
    /// weighted members that would give the same statistical precision.
    ///
    /// It is the member count for equal weights and approaches one when a
    /// single member dominates. Fails on an empty ensemble.
    pub fn effective_sample_size(&self) -> Result<f64, EnsembleError> {
        if self.is_empty() {
            return Err(EnsembleError::EmptySource);
        }
        let scaled = self.scaled_weights();
        let total = scaled.iter().sum::<f64>();
        let squares = scaled.iter().map(|weight| weight * weight).sum::<f64>();
        Ok(total * total / squares)
    }

    /// Member indices by descending weight; equal weights keep member order.
    ///
    /// Select the ten most populated members with
    /// `ensemble.select(ensemble.ranked_by_weight().into_iter().take(10))`.
    pub fn ranked_by_weight(&self) -> Vec<usize> {
        let mut indices = (0..self.members.len()).collect::<Vec<_>>();
        indices.sort_by(|&a, &b| self.members[b].weight.total_cmp(&self.members[a].weight));
        indices
    }

    /// The fewest members, by descending weight, whose weights together
    /// reach `fraction` of the total, in [`Self::ranked_by_weight`] order.
    ///
    /// `fraction` must lie in `(0, 1]`; `1` selects every member, including
    /// those too light to change a floating-point total. Pass the result to
    /// [`Self::select`] to keep, for example, 95% of a Boltzmann population.
    /// Fails on an empty ensemble.
    pub fn covering_weight(&self, fraction: f64) -> Result<Vec<usize>, EnsembleError> {
        if !(fraction > 0.0 && fraction <= 1.0) {
            return Err(EnsembleError::InvalidFraction);
        }
        if self.is_empty() {
            return Err(EnsembleError::EmptySource);
        }
        let ranked = self.ranked_by_weight();
        if fraction == 1.0 {
            return Ok(ranked);
        }
        let scaled = self.scaled_weights();
        // Sum in ranked order so the running total ends exactly at the total.
        let total = ranked.iter().map(|&index| scaled[index]).sum::<f64>();
        let target = fraction * total;
        let mut covered = 0.0;
        let count = ranked
            .iter()
            .position(|&index| {
                covered += scaled[index];
                covered >= target
            })
            .map_or(ranked.len(), |position| position + 1);
        Ok(ranked[..count].to_vec())
    }

    /// Member indices by ascending energy; equal energies keep member order.
    ///
    /// Fails with [`EnsembleError::MissingEnergy`] if any member has none.
    pub fn ranked_by_energy(&self) -> Result<Vec<usize>, EnsembleError> {
        let energies = self.energies()?;
        let mut indices = (0..energies.len()).collect::<Vec<_>>();
        indices.sort_by(|&a, &b| {
            energies[a]
                .partial_cmp(&energies[b])
                .expect("member energies are finite")
        });
        Ok(indices)
    }

    /// Indices, in member order, of the members whose energy lies at most
    /// `window` above the lowest member energy.
    ///
    /// `window` is a finite, nonnegative energy, such as the 6 kcal/mol
    /// window common in conformer searches. Fails if any member has no
    /// energy.
    pub fn within_energy_window(&self, window: Quantity<f64>) -> Result<Vec<usize>, EnsembleError> {
        let window = window.value_in(CANONICAL_ENERGY_UNIT)?;
        if !(window.is_finite() && window >= 0.0) {
            return Err(EnsembleError::InvalidEnergyWindow);
        }
        let energies = self.energies()?;
        let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
        Ok(energies
            .iter()
            .enumerate()
            .filter(|(_, energy)| *energy - minimum <= window)
            .map(|(index, _)| index)
            .collect())
    }

    /// Replaces every weight by the Boltzmann factor of the member energy
    /// relative to the lowest one, `exp(-(E - E_min) / RT)`, so the lowest
    /// member gets weight `1`.
    ///
    /// This treats each member as one distinct state, such as a conformer
    /// from a conformer search. Members already drawn from the distribution,
    /// such as molecular dynamics or Monte Carlo samples, are equally
    /// weighted and must not be Boltzmann weighted again. Previous weights,
    /// including degeneracies of collapsed members, are discarded; weight
    /// equivalent members individually before collapsing them.
    ///
    /// `temperature` must be finite and positive. Fails if any member has no
    /// energy or a weight would underflow to zero (restrict the ensemble to
    /// an energy window first). Failure changes nothing.
    pub fn assign_boltzmann_weights(
        &mut self,
        temperature: Quantity<f64>,
    ) -> Result<(), EnsembleError> {
        let temperature = temperature.value_in(CANONICAL_TEMPERATURE_UNIT)?;
        let thermal_energy = MOLAR_GAS_CONSTANT * temperature;
        if !(temperature.is_finite() && thermal_energy.is_finite() && thermal_energy > 0.0) {
            return Err(EnsembleError::InvalidTemperature);
        }
        let energies = self.energies()?;
        let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
        let weights = energies
            .iter()
            .map(|energy| checked_weight((-(energy - minimum) / thermal_energy).exp()))
            .collect::<Result<Vec<_>, _>>()?;
        for (member, weight) in self.members.iter_mut().zip(weights) {
            member.weight = weight;
        }
        Ok(())
    }

    /// Weights divided by the largest one, so sums cannot overflow.
    pub(crate) fn scaled_weights(&self) -> Vec<f64> {
        let maximum = self
            .members
            .iter()
            .map(|member| member.weight)
            .fold(0.0_f64, f64::max);
        self.members
            .iter()
            .map(|member| member.weight / maximum)
            .collect()
    }

    /// Every member energy in [`CANONICAL_ENERGY_UNIT`].
    fn energies(&self) -> Result<Vec<f64>, EnsembleError> {
        self.members
            .iter()
            .enumerate()
            .map(|(index, member)| {
                member
                    .energy
                    .map(Quantity::into_value)
                    .ok_or(EnsembleError::MissingEnergy { index })
            })
            .collect()
    }

    /// Applies one rigid transform per member transactionally.
    pub(crate) fn transform_members(
        &mut self,
        transforms: &[RigidTransform],
    ) -> Result<(), StructureError> {
        debug_assert_eq!(transforms.len(), self.members.len());
        let members = self
            .members
            .iter()
            .zip(transforms)
            .map(|(member, transform)| member.transformed(*transform))
            .collect::<Result<Vec<_>, _>>()?;
        self.members = members;
        Ok(())
    }

    fn check_index(&self, index: usize) -> Result<(), EnsembleError> {
        if index >= self.members.len() {
            return Err(EnsembleError::IndexOutOfRange {
                index,
                len: self.members.len(),
            });
        }
        Ok(())
    }

    fn bind(&self, member: &mut EnsembleMember) -> Result<(), EnsembleError> {
        validate_binding(member.atom_count(), &member.properties, &self.topology)?;
        member.properties.bind_bond_rows(self.topology.bond_count());
        Ok(())
    }
}

impl<'a> IntoIterator for &'a Ensemble {
    type Item = EnsembleMemberView<'a>;
    type IntoIter = EnsembleIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// One ensemble member: positions in dense topology atom order, an optional
/// periodic cell, optional crystallographic occupancies and B-factors,
/// structure properties, a statistical weight, and an optional energy.
///
/// The member carries no topology; it is bound when inserted into an
/// [`Ensemble`].
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleMember {
    positions: Positions,
    cell: Option<PeriodicCell>,
    occupancies: AtomScalars,
    b_factors: AtomScalars,
    properties: StructureProperties,
    weight: f64,
    energy: Option<Quantity<f64>>,
}

impl EnsembleMember {
    /// Creates a non-periodic member without occupancies, B-factors,
    /// properties, or energy, with a finite positive relative weight.
    pub fn new(positions: Positions, weight: f64) -> Result<Self, StructureError> {
        let mut member = Self::unit(positions);
        member.weight = checked_weight(weight)?;
        Ok(member)
    }

    /// A member of an equally weighted collection.
    fn unit(positions: Positions) -> Self {
        let properties = StructureProperties::new(positions.len());
        Self {
            positions,
            cell: None,
            occupancies: AtomScalars::default(),
            b_factors: AtomScalars::default(),
            properties,
            weight: 1.0,
            energy: None,
        }
    }

    /// An equally weighted member holding a model's structure state; the
    /// model's topology is dropped.
    pub(crate) fn from_model(model: Model) -> Self {
        Self {
            positions: model.positions,
            cell: model.cell,
            occupancies: model.occupancies,
            b_factors: model.b_factors,
            properties: model.properties,
            weight: 1.0,
            energy: None,
        }
    }

    /// An equally weighted member without occupancies or B-factors, for
    /// frames reinterpreted as an ensemble.
    pub(crate) fn from_frame_state(
        positions: Positions,
        cell: Option<PeriodicCell>,
        properties: StructureProperties,
    ) -> Self {
        Self {
            positions,
            cell,
            occupancies: AtomScalars::default(),
            b_factors: AtomScalars::default(),
            properties,
            weight: 1.0,
            energy: None,
        }
    }

    /// Returns this member with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.cell = Some(cell);
        self
    }

    /// Number of atoms, the length of every dense per-atom array of the member.
    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// The relative statistical weight; always finite and positive.
    pub const fn weight(&self) -> f64 {
        self.weight
    }

    /// Sets a finite positive relative weight.
    pub fn set_weight(&mut self, weight: f64) -> Result<(), StructureError> {
        self.weight = checked_weight(weight)?;
        Ok(())
    }

    /// Returns this member with a finite energy.
    pub fn with_energy(mut self, energy: Quantity<f64>) -> Result<Self, StructureError> {
        self.set_energy(Some(energy))?;
        Ok(self)
    }

    /// The member's energy in [`CANONICAL_ENERGY_UNIT`], if known.
    ///
    /// Energies of one ensemble share one scale, such as the electronic or
    /// free energies of conformers at one level of theory; only differences
    /// between members are meaningful. Record the scale in the ensemble
    /// properties.
    pub const fn energy(&self) -> Option<Quantity<f64>> {
        self.energy
    }

    /// Sets or clears a finite energy, converted to the canonical energy
    /// unit. The weight is unchanged.
    pub fn set_energy(&mut self, energy: Option<Quantity<f64>>) -> Result<(), StructureError> {
        self.energy = energy
            .map(|energy| {
                let energy = energy.into_unit(CANONICAL_ENERGY_UNIT)?;
                if !energy.value().is_finite() {
                    return Err(StructureError::NonFiniteEnergy);
                }
                Ok(energy)
            })
            .transpose()?;
        Ok(())
    }

    /// Applies a rigid transform transactionally: positions move and the
    /// cell rotates.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), StructureError> {
        *self = self.transformed(*transform)?;
        Ok(())
    }

    fn model_view<'a>(&'a self, topology: &'a Arc<Topology>) -> ModelView<'a> {
        ModelView::bound(
            topology,
            &self.positions,
            self.cell.as_ref(),
            &self.occupancies,
            &self.b_factors,
            &self.properties,
            StructureKind::Member {
                weight: self.weight,
                energy: self.energy,
            },
        )
    }

    /// Projects every dense row through a subset correspondence.
    fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, StructureError> {
        Ok(Self {
            positions: self.positions.select_indices(atoms)?,
            cell: self.cell,
            occupancies: self.occupancies.select(atoms),
            b_factors: self.b_factors.select(atoms),
            properties: self.properties.project(atoms, bonds)?,
            weight: self.weight,
            energy: self.energy,
        })
    }

    /// Positions move and the cell rotates.
    fn transformed(&self, transform: RigidTransform) -> Result<Self, StructureError> {
        Ok(Self {
            cell: transform_cell(self.cell, transform)?,
            positions: transform_positions(&self.positions, transform)?,
            occupancies: self.occupancies.clone(),
            b_factors: self.b_factors.clone(),
            properties: self.properties.clone(),
            weight: self.weight,
            energy: self.energy,
        })
    }

    fn storage_ref(&self) -> &Self {
        self
    }

    fn storage(&mut self) -> &mut Self {
        self
    }

    structure_readers!();
    site_readers!();
    structure_mutators!();
    site_mutators!();
}

/// The molar gas constant `N_A k_B` in kJ/(mol K), exact since the 2019 SI.
const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24e-3;

fn checked_weight(weight: f64) -> Result<f64, StructureError> {
    if weight.is_finite() && weight > 0.0 {
        Ok(weight)
    } else {
        Err(StructureError::InvalidWeight)
    }
}

/// One topology-bound [`EnsembleMember`] borrowed from an [`Ensemble`].
///
/// Reads member state through `Deref` and converts to a [`ModelView`]
/// through [`AsModelView`] without copying.
#[derive(Debug, Clone, Copy)]
pub struct EnsembleMemberView<'a> {
    topology: &'a Arc<Topology>,
    member: &'a EnsembleMember,
}

impl std::ops::Deref for EnsembleMemberView<'_> {
    type Target = EnsembleMember;

    fn deref(&self) -> &Self::Target {
        self.member
    }
}

impl<'a> EnsembleMemberView<'a> {
    /// Binds a member that already has every row `topology` requires.
    pub fn new(
        topology: &'a Arc<Topology>,
        member: &'a EnsembleMember,
    ) -> Result<Self, StructureError> {
        validate_binding(member.atom_count(), &member.properties, topology)?;
        if !is_bound(member.atom_count(), &member.properties, topology) {
            return Err(StructureError::BondCountMismatch {
                expected: topology.bond_count(),
                actual: member.properties.bonds().len(),
            });
        }
        Ok(Self { topology, member })
    }

    /// The topology the member is bound to, with the view's full lifetime.
    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    /// A shared handle to the bound topology.
    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    /// The borrowed member with the view's full lifetime.
    pub fn payload(self) -> &'a EnsembleMember {
        self.member
    }

    /// Borrows the member as the common kernel input, with the view's full
    /// lifetime.
    pub fn as_model_view(self) -> ModelView<'a> {
        self.member.model_view(self.topology)
    }

    /// Materializes the member as an owned model, dropping its weight and
    /// energy.
    pub fn to_model(self) -> Model {
        self.as_model_view().to_model()
    }
}

impl super::model::sealed::Sealed for EnsembleMemberView<'_> {}

impl AsModelView for EnsembleMemberView<'_> {
    fn as_model_view(&self) -> ModelView<'_> {
        EnsembleMemberView::as_model_view(*self)
    }
}

/// Dimension-preserving mutable access to one member of an [`Ensemble`].
///
/// Reads use `Deref`. There is deliberately no `DerefMut`; replace a whole
/// member through [`Ensemble::replace`].
#[derive(Debug)]
pub struct EnsembleMemberMut<'a> {
    member: &'a mut EnsembleMember,
}

impl std::ops::Deref for EnsembleMemberMut<'_> {
    type Target = EnsembleMember;

    fn deref(&self) -> &Self::Target {
        self.member
    }
}

impl EnsembleMemberMut<'_> {
    /// Sets a finite positive relative weight.
    pub fn set_weight(&mut self, weight: f64) -> Result<(), StructureError> {
        self.member.set_weight(weight)
    }

    /// Sets or clears a finite energy; the weight is unchanged.
    pub fn set_energy(&mut self, energy: Option<Quantity<f64>>) -> Result<(), StructureError> {
        self.member.set_energy(energy)
    }

    /// Applies a rigid transform transactionally: positions move and the
    /// cell rotates.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), StructureError> {
        self.member.apply_rigid_transform(transform)
    }

    fn storage(&mut self) -> &mut EnsembleMember {
        self.member
    }

    structure_mutators!();
    site_mutators!();
}

/// Iterator over the topology-bound members of an [`Ensemble`].
#[derive(Debug, Clone)]
pub struct EnsembleIter<'a> {
    topology: &'a Arc<Topology>,
    members: std::slice::Iter<'a, EnsembleMember>,
}

impl<'a> Iterator for EnsembleIter<'a> {
    type Item = EnsembleMemberView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let topology = self.topology;
        self.members
            .next()
            .map(|member| EnsembleMemberView { topology, member })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.members.size_hint()
    }
}

impl ExactSizeIterator for EnsembleIter<'_> {}

/// Failure of an [`Ensemble`] operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum EnsembleError {
    /// A constructor received no models, or weights of an empty ensemble
    /// were normalized.
    EmptySource,
    /// A member index is outside the ensemble.
    IndexOutOfRange {
        /// The requested member index.
        index: usize,
        /// The number of members.
        len: usize,
    },
    /// A model belongs to a different topology layout.
    TopologyMismatch,
    /// An ensemble-wide replacement has the wrong number of members.
    MemberCountMismatch {
        /// The ensemble's member count.
        expected: usize,
        /// The number of replacements supplied.
        actual: usize,
    },
    /// An energy-based operation found a member without an energy.
    MissingEnergy {
        /// The index of the member without an energy.
        index: usize,
    },
    /// A temperature is not finite and positive.
    InvalidTemperature,
    /// An energy window is not finite and nonnegative.
    InvalidEnergyWindow,
    /// A weight fraction is not in `(0, 1]`.
    InvalidFraction,
    /// A member's structural state is invalid; see the source.
    Structure(StructureError),
    /// Projecting the topology onto a selection failed; see the source.
    Subset(TopologySubsetError),
    /// Building the topology for the ensemble failed; see the source.
    TopologyBuild(TopologyBuildError),
    /// A quantity could not be converted to the required unit; see the source.
    Unit(Box<UnitError>),
}

impl fmt::Display for EnsembleError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySource => formatter.write_str("ensemble source is empty"),
            Self::IndexOutOfRange { index, len } => {
                write!(formatter, "index {index} is out of range for {len} members")
            }
            Self::TopologyMismatch => {
                formatter.write_str("model belongs to a different topology layout")
            }
            Self::MemberCountMismatch { expected, actual } => write!(
                formatter,
                "ensemble requires {expected} replacement members, but received {actual}"
            ),
            Self::MissingEnergy { index } => write!(formatter, "member {index} has no energy"),
            Self::InvalidTemperature => {
                formatter.write_str("temperature must be finite and positive")
            }
            Self::InvalidEnergyWindow => {
                formatter.write_str("energy window must be finite and nonnegative")
            }
            Self::InvalidFraction => {
                formatter.write_str("weight fraction must be greater than 0 and at most 1")
            }
            Self::Structure(_) => formatter.write_str("invalid ensemble member"),
            Self::Subset(_) => formatter.write_str("cannot subset ensemble"),
            Self::TopologyBuild(_) => formatter.write_str("cannot build topology"),
            Self::Unit(_) => formatter.write_str("invalid quantity unit"),
        }
    }
}

impl std::error::Error for EnsembleError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structure(error) => Some(error),
            Self::Subset(error) => Some(error),
            Self::TopologyBuild(error) => Some(error),
            Self::Unit(error) => Some(error.as_ref()),
            _ => None,
        }
    }
}

impl From<StructureError> for EnsembleError {
    fn from(error: StructureError) -> Self {
        Self::Structure(error)
    }
}

impl From<TopologySubsetError> for EnsembleError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Subset(error)
    }
}

impl From<TopologyBuildError> for EnsembleError {
    fn from(error: TopologyBuildError) -> Self {
        Self::TopologyBuild(error)
    }
}

impl From<UnitError> for EnsembleError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
