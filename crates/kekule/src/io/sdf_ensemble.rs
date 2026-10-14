//! Explicit assembly of interpreted SDF records into one ensemble.
use std::fmt;

use crate::structure::{Ensemble, EnsembleError, EnsembleMember};
use crate::units::{Quantity, KILOJOULE_PER_MOLE};

use super::sdf_document::SdfRecordInterpretation;
use super::sdf_write::{SDF_ENERGY_FIELD, SDF_WEIGHT_FIELD};

/// Failure to assemble SDF records into an ensemble.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum SdfEnsembleError {
    /// The records do not form an ensemble, for example because they are
    /// empty or do not share one topology layout.
    Ensemble(Box<EnsembleError>),
    /// A member metadata field of the record at `record` (counting from zero)
    /// is repeated or does not hold a valid value.
    InvalidField {
        record: usize,
        name: String,
        message: String,
    },
}

impl fmt::Display for SdfEnsembleError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Ensemble(error) => write!(formatter, "SDF records form no ensemble: {error}"),
            Self::InvalidField {
                record,
                name,
                message,
            } => write!(
                formatter,
                "SDF record {record} data field `{name}` {message}"
            ),
        }
    }
}

impl std::error::Error for SdfEnsembleError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Ensemble(error) => Some(error.as_ref()),
            Self::InvalidField { .. } => None,
        }
    }
}

/// Assembles interpreted SDF records, in order, into one ensemble on the
/// first record's topology. Every record must have the same topology layout
/// ([`Topology::same_layout`](crate::topology::Topology::same_layout)). Each
/// member takes its relative weight from the record's [`SDF_WEIGHT_FIELD`]
/// (1.0 when absent) and its energy from [`SDF_ENERGY_FIELD`] in kJ/mol (none
/// when absent), the fields the SDF writer emits for ensemble members.
pub fn sdf_ensemble<'a>(
    records: impl IntoIterator<Item = &'a SdfRecordInterpretation>,
) -> Result<Ensemble, SdfEnsembleError> {
    let records = records.into_iter().collect::<Vec<_>>();
    let topology = records
        .first()
        .ok_or_else(|| SdfEnsembleError::Ensemble(Box::new(EnsembleError::EmptySource)))?
        .model()
        .shared_topology();
    let mut members = Vec::with_capacity(records.len());
    for (index, record) in records.iter().enumerate() {
        if !record.model().topology().same_layout(&topology) {
            return Err(SdfEnsembleError::Ensemble(Box::new(
                EnsembleError::TopologyMismatch,
            )));
        }
        let invalid = |name: &str, message: String| SdfEnsembleError::InvalidField {
            record: index,
            name: name.to_owned(),
            message,
        };
        let number = |name: &str| -> Result<Option<f64>, SdfEnsembleError> {
            let mut values = record
                .data_fields()
                .iter()
                .filter(|field| field.name() == name);
            let Some(field) = values.next() else {
                return Ok(None);
            };
            if values.next().is_some() {
                return Err(invalid(name, "is repeated".to_owned()));
            }
            let value = field
                .value()
                .trim()
                .parse::<f64>()
                .map_err(|_| invalid(name, format!("is not a number: `{}`", field.value())))?;
            Ok(Some(value))
        };
        let weight = number(SDF_WEIGHT_FIELD)?.unwrap_or(1.0);
        let mut member = EnsembleMember::new(record.model().positions().clone(), weight)
            .map_err(|error| invalid(SDF_WEIGHT_FIELD, error.to_string()))?;
        if let Some(energy) = number(SDF_ENERGY_FIELD)? {
            member = member
                .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
                .map_err(|error| invalid(SDF_ENERGY_FIELD, error.to_string()))?;
        }
        members.push(member);
    }
    Ensemble::from_items(topology, members)
        .map_err(|error| SdfEnsembleError::Ensemble(Box::new(error)))
}
