mod features;
use crate::{error, explicit, identity, offxml::ASH_HASH, Result};
use kekule::{
    core::Molecule,
    query::{parse_smarts, QueryGraph},
    substructure::find_substructure_match,
    units::{Quantity, ELEMENTARY_CHARGE},
};
use ndarray::{Array1, Array2};
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::{collections::BTreeMap, io::Read, path::Path};

// Frozen bundle sizes. Check before allocating and cap the read as well, since
// a file can change between metadata inspection and reading.
fn read_exact_bundle_file(path: &Path, expected: usize) -> Result<Vec<u8>> {
    let file = std::fs::File::open(path).map_err(error)?;
    if !file.metadata().map_err(error)?.is_file()
        || file.metadata().map_err(error)?.len() != expected as u64
    {
        return Err(error(format!(
            "invalid Ash bundle file size: {} (expected {expected} bytes)",
            path.display()
        )));
    }
    let mut data = Vec::with_capacity(expected);
    file.take(expected as u64 + 1)
        .read_to_end(&mut data)
        .map_err(error)?;
    if data.len() != expected {
        return Err(error("Ash bundle file size changed while reading"));
    }
    Ok(data)
}

// OpenFF Toolkit's Molecule representation does not retain isotopic masses.
// Ash lookup keys and graph correspondence must use that same representation,
// while the caller's Kekule graph (and the general InChI adapter) retains them.
fn lookup_molecule(input: &Molecule) -> Result<Molecule> {
    if !input.atoms().any(|(_, a)| a.isotope.is_some()) {
        return explicit(input);
    }
    let mut editor = input.edit();
    for (id, atom) in input.atoms() {
        if atom.isotope.is_some() {
            editor.atom_mut(id).map_err(error)?.isotope = None;
        }
    }
    explicit(&editor.finish().map_err(error)?)
}

fn check_atom_limit(input: &Molecule) -> Result<()> {
    if input.atom_count() > 4096 {
        return Err(error("NAGL feature atom limit exceeded (4096)"));
    }
    Ok(())
}

/// Provenance of one molecule's assigned charges.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ChargeSource {
    Library { parameter_ids: Vec<String> },
    Lookup { inchi: String },
    Inference { model: String },
}
#[derive(Debug, Clone, PartialEq)]
pub struct ChargeAssignment {
    pub charges: Quantity<Vec<f64>>,
    pub source: ChargeSource,
}
#[derive(Debug, Deserialize)]
struct Tensor {
    shape: Vec<usize>,
    offset: usize,
    length: usize,
}
#[derive(Debug, Deserialize)]
struct Entry {
    inchi: String,
    mapped_smiles: String,
    charges: Vec<f32>,
}
#[derive(Debug, Deserialize)]
struct Domain {
    allowed_elements: Vec<u8>,
    forbidden_patterns: Vec<String>,
}
#[derive(Debug, Deserialize)]
struct Manifest {
    schema: u32,
    model: String,
    checkpoint_sha256: String,
    weights_sha256: String,
    tensors: BTreeMap<String, Tensor>,
    domain: Domain,
    lookup_tables: BTreeMap<String, Vec<Entry>>,
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
}
/// Native CPU inference for the frozen Ash 1.0.0 model shipped with Rosemary.
///
/// Load the data-only bundle produced by `benchmarks/openff/export_model.py`.
/// The checkpoint, manifest, and weight fingerprints are pinned. No Python,
/// PyTorch, RDKit, network, or arbitrary checkpoint deserialization is used here.
#[derive(Debug)]
pub struct NaglModel {
    name: String,
    layers: Vec<Convolution>,
    hidden: Linear,
    output: Linear,
    lookup: BTreeMap<String, Entry>,
    max_lookup_atoms: usize,
    elements: Vec<u8>,
    forbidden: Vec<QueryGraph>,
}
impl NaglModel {
    pub fn load(directory: impl AsRef<Path>) -> Result<Self> {
        let directory = directory.as_ref();
        let metadata = read_exact_bundle_file(&directory.join("model.json"), 2_050_045)?;
        if format!("{:x}", Sha256::digest(&metadata))
            != "f359ed50ada12a120195464ef0259837884e73826829e562652cb5e46ddb9425"
        {
            return Err(error("Ash model manifest checksum mismatch"));
        }
        let manifest: Manifest = serde_json::from_slice(&metadata).map_err(error)?;
        if manifest.schema != 1
            || manifest.checkpoint_sha256 != ASH_HASH
            || manifest.model != "openff-gnn-am1bcc-1.0.0.pt"
        {
            return Err(error("unsupported NAGL model"));
        }
        let raw = read_exact_bundle_file(&directory.join("weights.bin"), 10_852_364)?;
        if raw.len() != 10_852_364
            || format!("{:x}", Sha256::digest(&raw)) != manifest.weights_sha256
        {
            return Err(error("Ash weights checksum mismatch"));
        }
        let weights = raw
            .chunks_exact(4)
            .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
            .collect::<Vec<_>>();
        if weights.iter().any(|v| !v.is_finite()) {
            return Err(error("nonfinite NAGL weights"));
        }
        let linear = |prefix: &str, input: usize, output: usize, bias: bool| -> Result<Linear> {
            let tensor = |suffix: &str, shape: &[usize]| -> Result<Vec<f32>> {
                let t = manifest
                    .tensors
                    .get(&format!("{prefix}.{suffix}"))
                    .ok_or_else(|| error("missing NAGL tensor"))?;
                if t.shape != shape || t.length != shape.iter().product::<usize>() {
                    return Err(error("NAGL tensor shape mismatch"));
                }
                Ok(weights
                    .get(
                        t.offset
                            ..t.offset
                                .checked_add(t.length)
                                .ok_or_else(|| error("tensor overflow"))?,
                    )
                    .ok_or_else(|| error("tensor offset out of bounds"))?
                    .to_vec())
            };
            Ok(Linear {
                weight: Array2::from_shape_vec(
                    (output, input),
                    tensor("weight", &[output, input])?,
                )
                .map_err(error)?,
                bias: if bias {
                    Array1::from_vec(tensor("bias", &[output])?)
                } else {
                    Array1::zeros(output)
                },
            })
        };
        let mut layers = Vec::new();
        for i in 0..6 {
            let input = if i == 0 { 22 } else { 512 };
            layers.push(Convolution {
                own: linear(
                    &format!("convolution_module.gcn_layers.{i}.fc_self"),
                    input,
                    512,
                    true,
                )?,
                neighbor: linear(
                    &format!("convolution_module.gcn_layers.{i}.fc_neigh"),
                    input,
                    512,
                    false,
                )?,
            });
        }
        let hidden = linear(
            "readout_modules.am1bcc_charges.readout_layers.0",
            512,
            128,
            true,
        )?;
        let output = linear(
            "readout_modules.am1bcc_charges.readout_layers.3",
            128,
            3,
            true,
        )?;
        let forbidden = manifest
            .domain
            .forbidden_patterns
            .iter()
            .map(|p| parse_smarts(p).map_err(error))
            .collect::<Result<Vec<_>>>()?;
        let entries = manifest
            .lookup_tables
            .into_iter()
            .find(|(k, _)| k == "am1bcc_charges")
            .ok_or_else(|| error("missing charge lookup table"))?
            .1;
        let max_lookup_atoms = entries.iter().map(|e| e.charges.len()).max().unwrap_or(0);
        let lookup = entries.into_iter().map(|e| (e.inchi.clone(), e)).collect();
        Ok(Self {
            name: manifest.model,
            layers,
            hidden,
            output,
            lookup,
            max_lookup_atoms,
            elements: manifest.domain.allowed_elements,
            forbidden,
        })
    }
    pub fn lookup_entry_count(&self) -> usize {
        self.lookup.len()
    }
    /// Full fixed-H charge lookup identifier, accepting explicit hydrogens.
    /// Isotopes are ignored to match OpenFF Toolkit's charge representation.
    /// This diagnostic is limited to 1024 atoms by the InChI adapter. Charge
    /// assignment avoids InChI when the frozen table cannot contain the input.
    pub fn lookup_identifier(&self, molecule: &Molecule) -> Result<String> {
        identity::fixed_h_inchi(&lookup_molecule(molecule)?)
    }
    /// Ash features in molecule atom iteration order, with 22 columns.
    pub fn atom_features(&self, molecule: &Molecule) -> Result<Vec<Vec<f32>>> {
        check_atom_limit(molecule)?;
        let m = explicit(molecule)?;
        self.check_domain(&m)?;
        Ok(features::features(&m)?
            .rows()
            .into_iter()
            .map(|r| r.to_vec())
            .collect())
    }
    /// Lookup first; otherwise evaluate GraphSAGE and conserve molecular charge.
    pub fn assign_charges(&self, molecule: &Molecule) -> Result<ChargeAssignment> {
        check_atom_limit(molecule)?;
        let molecule = lookup_molecule(molecule)?;
        // A full-graph lookup hit must have one charge for every input atom.
        // The fingerprint-pinned table has at most 11 atoms per entry. Do not
        // let the identifier library's unrelated size limit block inference.
        if molecule.atom_count() > self.max_lookup_atoms {
            return self.infer_prepared(&molecule);
        }
        let inchi = identity::fixed_h_inchi(&molecule)?;
        if let Some(entry) = self.lookup.get(&inchi) {
            let mut components =
                kekule::smiles::to_molecules(&entry.mapped_smiles).map_err(error)?;
            if components.len() != 1 {
                return Err(error("lookup entry is disconnected"));
            }
            let reference = explicit(&components.remove(0))?;
            let mapping = identity::mapping(&molecule, &reference)?;
            let mut values = Vec::new();
            for atom in mapping {
                let map = reference
                    .atom(atom)
                    .map_err(error)?
                    .atom_map
                    .ok_or_else(|| error("unmapped lookup atom"))?;
                let q = entry
                    .charges
                    .get(map.checked_sub(1).ok_or_else(|| error("zero lookup map"))? as usize)
                    .ok_or_else(|| error("lookup charge array mismatch"))?;
                values.push(f64::from(*q));
            }
            return normalize_charges(
                values,
                molecule.formal_charge(),
                ChargeSource::Lookup { inchi },
            );
        }
        self.infer_prepared(&molecule)
    }
    /// Evaluate the network directly, bypassing lookup, for reference validation.
    pub fn infer_charges(&self, molecule: &Molecule) -> Result<ChargeAssignment> {
        check_atom_limit(molecule)?;
        self.infer_prepared(&explicit(molecule)?)
    }
    fn check_domain(&self, molecule: &Molecule) -> Result<()> {
        for (_, atom) in molecule.atoms() {
            if !self.elements.contains(&atom.element.atomic_number()) {
                return Err(error(format!(
                    "element {} is outside the Ash domain",
                    atom.element
                )));
            }
        }
        for pattern in &self.forbidden {
            if find_substructure_match(molecule, pattern)
                .map_err(error)?
                .is_some()
            {
                return Err(error("molecule matches a forbidden Ash domain pattern"));
            }
        }
        Ok(())
    }
    fn infer_prepared(&self, molecule: &Molecule) -> Result<ChargeAssignment> {
        self.check_domain(molecule)?;
        let mut x = features::features(molecule)?;
        let neighbors = features::adjacency(molecule)?;
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
            x.mapv_inplace(|v| v.max(0.0));
        }
        x = self.hidden.apply(&x);
        x.mapv_inplace(|v| 1.0 / (1.0 + (-v).exp()));
        let output = self.output.apply(&x);
        let mut priors = 0f32;
        let mut inverse_sum = 0f32;
        let mut e_sum = 0f32;
        for row in output.rows() {
            priors += row[0];
            inverse_sum += 1.0 / row[2];
            e_sum += row[1] / row[2];
        }
        let fraction = (priors - molecule.formal_charge() as f32 - e_sum) / inverse_sum;
        let values = output
            .rows()
            .into_iter()
            .map(|row| f64::from(row[0] - row[1] / row[2] - fraction / row[2]))
            .collect();
        normalize_charges(
            values,
            molecule.formal_charge(),
            ChargeSource::Inference {
                model: self.name.clone(),
            },
        )
    }
}
fn normalize_charges(
    mut values: Vec<f64>,
    formal: i64,
    source: ChargeSource,
) -> Result<ChargeAssignment> {
    if values.is_empty() || values.iter().any(|q| !q.is_finite()) {
        return Err(error("NAGL produced nonfinite or empty charges"));
    }
    let correction = (formal as f64 - values.iter().sum::<f64>()) / values.len() as f64;
    for q in &mut values {
        *q += correction;
    }
    Ok(ChargeAssignment {
        charges: Quantity::new(values, ELEMENTARY_CHARGE),
        source,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isotope_charge_key_matches_openff_without_mutating_input() {
        // PubChem CID 143783 from the locked pubchem-100k corpus. OpenFF
        // Toolkit 0.19.0 reports the isotope-free fixed-H identifier below.
        let mut input = kekule::smiles::to_molecules("[2H]C([2H])([2H])C([2H])([2H])SCC")
            .unwrap()
            .remove(0);
        input.perceive().unwrap();
        input.add_hydrogens().unwrap();
        let before = input.clone();
        let prepared = lookup_molecule(&input).unwrap();
        let expected = "InChI=1/C4H10S/c1-3-5-4-2/h3-4H2,1-2H3";
        assert_eq!(identity::fixed_h_inchi(&prepared).unwrap(), expected);
        assert_eq!(input, before);
        assert_eq!(
            prepared.atom_ids().collect::<Vec<_>>(),
            input.atom_ids().collect::<Vec<_>>()
        );
        assert!(identity::fixed_h_inchi(&explicit(&input).unwrap())
            .unwrap()
            .contains("/i1D3,3D2"));
    }

    #[test]
    fn bundle_reader_rejects_truncation_oversize_and_missing_files() {
        let path = std::env::temp_dir().join(format!("kekule-ash-size-{}.bin", std::process::id()));
        let file = std::fs::File::create(&path).unwrap();
        file.set_len(1_000_000_000).unwrap();
        drop(file);
        assert!(read_exact_bundle_file(&path, 8)
            .unwrap_err()
            .to_string()
            .contains("file size"));
        std::fs::write(&path, [1, 2, 3]).unwrap();
        assert!(read_exact_bundle_file(&path, 8).is_err());
        assert_eq!(read_exact_bundle_file(&path, 3).unwrap(), [1, 2, 3]);
        std::fs::remove_file(&path).unwrap();
        assert!(read_exact_bundle_file(&path, 3).is_err());
    }

    #[test]
    fn invalid_charge_outputs_fail_without_publishing() {
        let source = ChargeSource::Inference {
            model: "test".into(),
        };
        for values in [vec![], vec![f64::NAN], vec![f64::INFINITY]] {
            assert!(normalize_charges(values, 0, source.clone()).is_err());
        }
        let q = normalize_charges(vec![0.2, 0.3], 1, source).unwrap();
        assert!((q.charges.value().iter().sum::<f64>() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn inference_size_guard_does_not_require_perception() {
        let molecule = kekule::smiles::to_molecules(&"C".repeat(4097))
            .unwrap()
            .remove(0);
        assert!(check_atom_limit(&molecule)
            .unwrap_err()
            .to_string()
            .contains("4096"));
    }
}
