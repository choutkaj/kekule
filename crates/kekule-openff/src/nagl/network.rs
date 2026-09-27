//! Evaluation-only GraphSAGE and charge readout. Chemistry stays in `features`.
use super::config::{Activation, Config, Feature};
use crate::{error, Result};
use ndarray::{Array1, Array2};
use serde::Deserialize;
use std::collections::BTreeMap;

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
    pub fn load(
        config: &Config,
        mut tensors: BTreeMap<String, Tensor>,
        raw: &[u8],
    ) -> Result<Self> {
        if !raw.len().is_multiple_of(4) {
            return Err(error("NAGL weights must be little-endian float32"));
        }
        let weights = raw
            .chunks_exact(4)
            .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
            .collect::<Vec<_>>();
        if weights.iter().any(|v| !v.is_finite()) {
            return Err(error("nonfinite NAGL weights"));
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
                return Err(error("invalid NAGL tensor storage layout"));
            }
            end = end
                .checked_add(length)
                .ok_or_else(|| error("NAGL tensor overflow"))?;
        }
        if end != weights.len() {
            return Err(error("NAGL tensor storage size mismatch"));
        }
        let mut tensor = |name: String, shape: &[usize]| -> Result<Vec<f32>> {
            let t = tensors
                .remove(&name)
                .ok_or_else(|| error(format!("missing NAGL tensor {name}")))?;
            if t.shape != shape || t.length != shape.iter().product::<usize>() {
                return Err(error(format!("NAGL tensor shape mismatch: {name}")));
            }
            Ok(weights
                .get(t.offset..t.offset + t.length)
                .ok_or_else(|| error("NAGL tensor offset out of bounds"))?
                .to_vec())
        };
        let mut linear =
            |prefix: &str, input: usize, output: usize, bias: bool| -> Result<Linear> {
                Ok(Linear {
                    weight: Array2::from_shape_vec(
                        (output, input),
                        tensor(format!("{prefix}.weight"), &[output, input])?,
                    )
                    .map_err(error)?,
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
            return Err(error("unexpected NAGL tensors"));
        }
        Ok(Self { layers, readout })
    }

    pub fn charges(
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
        let mut priors = 0f32;
        let mut inverse_sum = 0f32;
        let mut e_sum = 0f32;
        for row in x.rows() {
            if row.iter().any(|v| !v.is_finite()) || row[2] == 0.0 {
                return Err(error("invalid NAGL charge-equilibration output"));
            }
            priors += row[0];
            inverse_sum += 1.0 / row[2];
            e_sum += row[1] / row[2];
        }
        let fraction = (priors - formal_charge as f32 - e_sum) / inverse_sum;
        if !fraction.is_finite() {
            return Err(error("singular NAGL charge equilibration"));
        }
        Ok(x.rows()
            .into_iter()
            .map(|row| f64::from(row[0] - row[1] / row[2] - fraction / row[2]))
            .collect())
    }
}
