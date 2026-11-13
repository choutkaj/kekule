//! Evaluation-only GraphSAGE and charge readout. Chemistry stays in `features`.
use super::config::{Activation, Config, Feature};
use crate::{Error, ErrorKind, Result};

use ndarray::{Array1, Array2};
use serde::Deserialize;
use std::collections::BTreeMap;

fn error(detail: impl std::fmt::Display) -> Error {
    Error::new(ErrorKind::Structure, detail)
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Tensor {
    shape: Vec<usize>,
    offset: usize,
    length: usize,
}
#[derive(Debug)]
struct Linear {
    weight: Array2<f32>,
    bias: Array1<f32>,
}
impl Linear {
    fn apply(&self, x: &Array2<f32>) -> Array2<f32> {
        x.dot(&self.weight.t()) + &self.bias
    }
}
#[derive(Debug)]
struct Convolution {
    own: Linear,
    neighbor: Linear,
    activation: Activation,
}
#[derive(Debug)]
pub(super) struct Network {
    layers: Vec<Convolution>,
    readout: Vec<(Linear, Activation)>,
}
impl Network {
    pub(crate) fn load(
        config: &Config,
        mut tensors: BTreeMap<String, Tensor>,
        raw: &[u8],
    ) -> Result<Self> {
        let (words, remainder) = raw.as_chunks::<4>();
        if !remainder.is_empty() {
            return Err(error(format!(
                "The weights have {} bytes, not a whole number of float32 values.",
                raw.len()
            )));
        }
        let weights = words
            .iter()
            .map(|&b| f32::from_le_bytes(b))
            .collect::<Vec<_>>();
        if weights.iter().any(|v| !v.is_finite()) {
            return Err(error("The weights hold a nonfinite value."));
        }
        // Every byte belongs to one named tensor. No hidden, overlapping or unused tensors.
        let mut spans = tensors
            .values()
            .map(|t| (t.offset, t.length))
            .collect::<Vec<_>>();
        spans.sort_unstable();
        let mut end = 0usize;
        for (offset, length) in spans {
            if offset != end || length == 0 {
                return Err(error("The weight tensors overlap or leave gaps."));
            }
            end = end
                .checked_add(length)
                .ok_or_else(|| error("The weight tensors end beyond the addressable size."))?;
        }
        if end != weights.len() {
            return Err(error(format!(
                "The weight tensors cover {end} values of {}.",
                weights.len()
            )));
        }
        let mut tensor = |name: String, shape: &[usize]| -> Result<Vec<f32>> {
            let t = tensors
                .remove(&name)
                .ok_or_else(|| error(format!("The weight tensor {name:?} is missing.")))?;
            if t.shape != shape || t.length != shape.iter().product::<usize>() {
                return Err(error(format!(
                    "The weight tensor {name:?} has shape {:?}, not {shape:?}.",
                    t.shape
                )));
            }
            Ok(weights
                .get(t.offset..t.offset + t.length)
                .ok_or_else(|| {
                    error(format!(
                        "The weight tensor {name:?} ends beyond the weights."
                    ))
                })?
                .to_vec())
        };
        let mut linear =
            |prefix: &str, input: usize, output: usize, bias: bool| -> Result<Linear> {
                Ok(Linear {
                    weight: Array2::from_shape_vec(
                        (output, input),
                        tensor(format!("{prefix}.weight"), &[output, input])?,
                    )
                    .map_err(|e| {
                        Error::caused(
                            ErrorKind::Structure,
                            format!("Cannot shape the weight tensor \"{prefix}.weight\"."),
                            e,
                        )
                    })?,
                    bias: if bias {
                        Array1::from_vec(tensor(format!("{prefix}.bias"), &[output])?)
                    } else {
                        Array1::zeros(output)
                    },
                })
            };
        let mut input = config.atom_features.iter().map(Feature::width).sum();
        let mut layers = Vec::new();
        for (i, layer) in config.convolution.layers.iter().enumerate() {
            let output = layer.hidden_feature_size;
            layers.push(Convolution {
                own: linear(
                    &format!("convolution_module.gcn_layers.{i}.fc_self"),
                    input,
                    output,
                    true,
                )?,
                neighbor: linear(
                    &format!("convolution_module.gcn_layers.{i}.fc_neigh"),
                    input,
                    output,
                    false,
                )?,
                activation: layer.activation_function,
            });
            input = output;
        }
        let (name, specification) = config.readouts.first_key_value().unwrap();
        let mut readout = Vec::new();
        for (i, layer) in specification.layers.iter().enumerate() {
            readout.push((
                linear(
                    &format!("readout_modules.{name}.readout_layers.{}", 3 * i),
                    input,
                    layer.hidden_feature_size,
                    true,
                )?,
                layer.activation_function,
            ));
            input = layer.hidden_feature_size;
        }
        readout.push((
            linear(
                &format!(
                    "readout_modules.{name}.readout_layers.{}",
                    3 * specification.layers.len()
                ),
                input,
                3,
                true,
            )?,
            Activation::Identity,
        ));
        if !tensors.is_empty() {
            return Err(error(format!(
                "The weights hold unused tensors: {}.",
                tensors
                    .keys()
                    .map(|name| format!("{name:?}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            )));
        }
        Ok(Self { layers, readout })
    }

    pub(crate) fn charges(
        &self,
        mut x: Array2<f32>,
        neighbors: &[Vec<usize>],
        formal_charge: i64,
    ) -> Result<Vec<f64>> {
        for layer in &self.layers {
            let mut mean = Array2::zeros(x.raw_dim());
            for (i, adjacent) in neighbors.iter().enumerate() {
                for &j in adjacent {
                    for k in 0..x.ncols() {
                        mean[(i, k)] += x[(j, k)];
                    }
                }
                if !adjacent.is_empty() {
                    for k in 0..x.ncols() {
                        mean[(i, k)] /= adjacent.len() as f32;
                    }
                }
            }
            x = layer.own.apply(&x) + layer.neighbor.apply(&mean);
            x.mapv_inplace(|v| layer.activation.apply(v));
        }
        for (linear, activation) in &self.readout {
            x = linear.apply(&x);
            x.mapv_inplace(|v| activation.apply(v));
        }
        // Molecule-wide sums accumulate in f64 so their rounding error does not
        // grow with the number of atoms; per-atom terms stay in float32 as in
        // the network.
        let mut priors = 0f64;
        let mut inverse_sum = 0f64;
        let mut e_sum = 0f64;
        for row in x.rows() {
            if row.iter().any(|v| !v.is_finite()) || row[2] == 0.0 {
                return Err(Error::new(
                    ErrorKind::Charges,
                    "Charge equilibration produced a nonfinite value or a zero hardness.",
                ));
            }
            priors += f64::from(row[0]);
            inverse_sum += f64::from(1.0 / row[2]);
            e_sum += f64::from(row[1] / row[2]);
        }
        let fraction = ((priors - formal_charge as f64 - e_sum) / inverse_sum) as f32;
        if !fraction.is_finite() {
            return Err(Error::new(
                ErrorKind::Charges,
                "Charge equilibration is singular.",
            ));
        }
        Ok(x.rows()
            .into_iter()
            .map(|row| f64::from(row[0] - row[1] / row[2] - fraction / row[2]))
            .collect())
    }
}
