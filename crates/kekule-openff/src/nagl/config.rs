//! The supported NAGL configuration, shared by bundle validation and execution.
use crate::{error, Result};
use kekule::core::Element;
use serde::Deserialize;
use std::collections::BTreeMap;

#[derive(Debug, Deserialize)]
#[serde(tag = "name", deny_unknown_fields)]
pub(super) enum Feature {
    #[serde(rename = "atomic_element")]
    Element { categories: Vec<String> },
    #[serde(rename = "atom_connectivity")]
    Connectivity { categories: Vec<usize> },
    #[serde(rename = "atom_average_formal_charge")]
    AverageFormalCharge,
    #[serde(rename = "atom_in_ring_of_size")]
    Ring { ring_size: usize },
}
impl Feature {
    pub fn width(&self) -> usize {
        match self {
            Self::Element { categories } => categories.len(),
            Self::Connectivity { categories } => categories.len(),
            _ => 1,
        }
    }
    fn validate(&self) -> Result<()> {
        fn unique<T: PartialEq>(values: &[T]) -> bool {
            !values.is_empty()
                && values
                    .iter()
                    .enumerate()
                    .all(|(i, v)| !values[..i].contains(v))
        }
        let valid = match self {
            Self::Element { categories } => {
                unique(categories) && categories.iter().all(|s| Element::from_symbol(s).is_some())
            }
            Self::Connectivity { categories } => {
                unique(categories) && categories.iter().all(|&n| n <= 6)
            }
            Self::Ring { ring_size } => (3..=6).contains(ring_size),
            Self::AverageFormalCharge => true,
        };
        if !valid {
            return Err(error("unsupported NAGL feature categories or ring size"));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "lowercase")]
pub(super) enum Activation {
    Relu,
    Sigmoid,
    Identity,
}
impl Activation {
    pub fn apply(self, x: f32) -> f32 {
        match self {
            Self::Relu => {
                if x < 0.0 {
                    0.0
                } else {
                    x
                }
            }
            Self::Sigmoid => 1.0 / (1.0 + (-x).exp()),
            Self::Identity => x,
        }
    }
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Layer {
    pub hidden_feature_size: usize,
    pub activation_function: Activation,
    pub dropout: f64,
    pub aggregator_type: Option<String>,
}
impl Layer {
    fn validate(&self, convolution: bool) -> Result<()> {
        if !(1..=2048).contains(&self.hidden_feature_size)
            || !(0.0..=1.0).contains(&self.dropout)
            || (convolution && self.aggregator_type.as_deref() != Some("mean"))
            || (!convolution && self.aggregator_type.is_some())
        {
            return Err(error(
                "unsupported NAGL layer: expected width 1..2048, dropout 0..1 and mean aggregation",
            ));
        }
        // Dropout is disabled in evaluation mode, including models trained with nonzero dropout.
        Ok(())
    }
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Convolution {
    pub architecture: String,
    pub layers: Vec<Layer>,
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Readout {
    pub pooling: String,
    pub layers: Vec<Layer>,
    pub postprocess: String,
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Config {
    pub version: String,
    pub atom_features: Vec<Feature>,
    pub bond_features: Vec<serde_json::Value>,
    pub convolution: Convolution,
    pub readouts: BTreeMap<String, Readout>,
}
impl Config {
    pub fn validate(&self) -> Result<()> {
        if self.version != "0.1"
            || !self.bond_features.is_empty()
            || self.convolution.architecture != "SAGEConv"
            || !(1..=16).contains(&self.convolution.layers.len())
            || self.readouts.len() != 1
        {
            return Err(error("unsupported NAGL configuration: expected version 0.1, SAGEConv, no bond features and one charge readout"));
        }
        if self.atom_features.is_empty()
            || self.atom_features.len() > 64
            || self.atom_features.iter().map(Feature::width).sum::<usize>() > 256
        {
            return Err(error("NAGL feature count exceeds supported limits"));
        }
        for feature in &self.atom_features {
            feature.validate()?;
        }
        for layer in &self.convolution.layers {
            layer.validate(true)?;
        }
        let (name, readout) = self.readouts.first_key_value().unwrap();
        if name.is_empty()
            || readout.pooling != "atoms"
            || readout.postprocess != "regularized_compute_partial_charges"
            || readout.layers.len() > 16
        {
            return Err(error("unsupported NAGL readout: expected atom pooling and regularized_compute_partial_charges"));
        }
        for layer in &readout.layers {
            layer.validate(false)?;
        }
        Ok(())
    }
}
