//! What a run sampled: [`Sampling`], the sampling record of a structure set.

use crate::properties::{PropertyColumn, PropertyKey, PropertyTable, PropertyValue};
use crate::units::{Unit, CANONICAL_ENERGY_UNIT};

use super::{StructureSetError, ThermodynamicState};

/// The sampling record of a structure set: the state a run sampled, from
/// which MD step its frames are samples of that state, and which of the
/// set's columns hold the unbiased potential energy and the bias.
///
/// The record states facts the engine knows, never a weight. From
/// [`Self::from_step`] on, each geometry is a sample of [`Self::state`]
/// under the bias of its row: the sum of the [`Self::bias`] columns, a
/// static bias such as GaMD's boost or a per-frame one such as
/// metadynamics' `V(s, t) − c(t)`. An analysis derives weights of the
/// sampled state from the record
/// ([`Weights::from_bias`](crate::statistics::Weights::from_bias)). The
/// potential-energy column, when named, records each frame's unbiased
/// energy. The method is a label, such as `GaMD dual boost`; kekule gives
/// it no meaning.
///
/// A set has at most one record ([`StructureSet::sampling`](super::StructureSet::sampling)),
/// and every column it names is a real column of energies or absent. A
/// column is absent when all its cells are missing, as after selecting
/// geometries that have none of its values.
///
/// ```
/// use kekule::properties::PropertyKey;
/// use kekule::structure::{Sampling, ThermodynamicState};
/// use kekule::units::{Quantity, BAR, KELVIN};
///
/// let state = ThermodynamicState::isothermal_isobaric(
///     Quantity::new(300.0, KELVIN),
///     Quantity::new(1.0, BAR),
/// )?;
/// let sampling = Sampling::new(state)
///     .with_from_step(2_600_000)
///     .with_method("GaMD dual boost")
///     .with_potential_energy(PropertyKey::new("potential_energy")?)?
///     .with_bias([PropertyKey::new("bias.potential")?, PropertyKey::new("bias.dihedral")?])?;
/// assert_eq!(sampling.bias().len(), 2);
/// assert_eq!(sampling.columns().count(), 3);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Sampling {
    state: ThermodynamicState,
    from_step: u64,
    potential_energy: Option<PropertyKey>,
    method: Option<String>,
    bias: Vec<PropertyKey>,
}

impl Sampling {
    /// A run that sampled `state` from step 0, without a bias, a method, or
    /// a potential-energy column.
    pub const fn new(state: ThermodynamicState) -> Self {
        Self {
            state,
            from_step: 0,
            potential_energy: None,
            method: None,
            bias: Vec::new(),
        }
    }

    /// Returns this record whose frames are samples from MD step `step` on.
    #[must_use]
    pub const fn with_from_step(mut self, step: u64) -> Self {
        self.from_step = step;
        self
    }

    /// Returns this record with a method label.
    #[must_use]
    pub fn with_method(mut self, method: impl Into<String>) -> Self {
        self.method = Some(method.into());
        self
    }

    /// Returns this record whose `column` holds each frame's unbiased
    /// potential energy.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::RepeatedSamplingColumn`] when `column` is a bias
    /// column of the record.
    pub fn with_potential_energy(mut self, column: PropertyKey) -> Result<Self, StructureSetError> {
        if self.bias.contains(&column) {
            return Err(StructureSetError::RepeatedSamplingColumn { key: column });
        }
        self.potential_energy = Some(column);
        Ok(self)
    }

    /// Returns this record whose bias is the sum of `columns`, in this
    /// order; no columns means no bias.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::RepeatedSamplingColumn`] for a column given
    /// twice or that is the potential-energy column.
    pub fn with_bias(
        mut self,
        columns: impl IntoIterator<Item = PropertyKey>,
    ) -> Result<Self, StructureSetError> {
        let mut bias = Vec::new();
        for column in columns {
            if bias.contains(&column) || self.potential_energy.as_ref() == Some(&column) {
                return Err(StructureSetError::RepeatedSamplingColumn { key: column });
            }
            bias.push(column);
        }
        self.bias = bias;
        Ok(self)
    }

    /// The state the run sampled.
    pub const fn state(&self) -> &ThermodynamicState {
        &self.state
    }

    /// The first MD step whose frame is a sample of the state.
    pub const fn from_step(&self) -> u64 {
        self.from_step
    }

    /// The column of each frame's unbiased potential energy, when named.
    pub const fn potential_energy(&self) -> Option<&PropertyKey> {
        self.potential_energy.as_ref()
    }

    /// The method label, when given.
    pub fn method(&self) -> Option<&str> {
        self.method.as_deref()
    }

    /// The bias columns, whose sum is each frame's bias; empty without a
    /// bias.
    pub fn bias(&self) -> &[PropertyKey] {
        &self.bias
    }

    /// Every column the record names: the potential energy first, then the
    /// bias columns.
    pub fn columns(&self) -> impl Iterator<Item = &PropertyKey> {
        self.potential_energy.iter().chain(&self.bias)
    }

    /// Whether the record names `key`.
    pub(crate) fn names(&self, key: &PropertyKey) -> bool {
        self.columns().any(|column| column == key)
    }

    /// Checks that each column the record names is a real column of
    /// energies in `columns`, or absent.
    pub(crate) fn check_columns<R>(
        &self,
        columns: &PropertyTable<R>,
    ) -> Result<(), StructureSetError> {
        for key in self.columns() {
            match columns.get(key) {
                None => {}
                Some(PropertyColumn::Real { unit, .. }) if is_energy(*unit) => {}
                Some(_) => return Err(StructureSetError::SamplingColumn { key: key.clone() }),
            }
        }
        Ok(())
    }

    /// Checks that `value`, written to the column `key`, keeps the record's
    /// columns real energies.
    pub(crate) fn check_value(
        &self,
        key: &PropertyKey,
        value: Option<&PropertyValue>,
    ) -> Result<(), StructureSetError> {
        match value {
            Some(PropertyValue::Real { unit, .. }) if !is_energy(*unit) => {}
            Some(PropertyValue::Real { .. }) | None => return Ok(()),
            Some(_) => {}
        }
        if self.names(key) {
            Err(StructureSetError::SamplingColumn { key: key.clone() })
        } else {
            Ok(())
        }
    }
}

/// Whether `unit` is an energy unit, such as kJ/mol.
pub(crate) fn is_energy(unit: Unit) -> bool {
    unit.is_compatible(CANONICAL_ENERGY_UNIT)
}
