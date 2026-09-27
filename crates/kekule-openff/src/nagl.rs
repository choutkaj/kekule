mod bundle;
mod config;
mod features;
mod network;
use crate::{error, explicit, identity, Result};
use kekule::{
    core::Molecule,
    query::{parse_smarts, QueryGraph},
    substructure::find_substructure_match,
    units::{Quantity, ELEMENTARY_CHARGE},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

/// Checkpoint identity declared by an OFFXML handler or exported model bundle.
/// The checksum identifies the original checkpoint, not the converted weights file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelIdentity {
    model_file: String,
    checkpoint_sha256: String,
}
impl ModelIdentity {
    pub(crate) fn new(model_file: String, checkpoint_sha256: String) -> Result<Self> {
        if model_file.trim().is_empty()
            || checkpoint_sha256.len() != 64
            || !checkpoint_sha256.bytes().all(|b| b.is_ascii_hexdigit())
        {
            return Err(error(
                "NAGL requires a nonempty model_file and a 64-digit SHA-256 model_file_hash",
            ));
        }
        Ok(Self {
            model_file,
            checkpoint_sha256: checkpoint_sha256.to_ascii_lowercase(),
        })
    }
    pub fn model_file(&self) -> &str {
        &self.model_file
    }
    pub fn checkpoint_sha256(&self) -> &str {
        &self.checkpoint_sha256
    }
}

// OpenFF Toolkit's Molecule representation does not retain isotopic masses.
// NAGL lookup keys and graph correspondence must use that same representation,
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
    Lookup { inchi: String, model: ModelIdentity },
    Inference { model: ModelIdentity },
}
#[derive(Debug, Clone, PartialEq)]
pub struct ChargeAssignment {
    pub charges: Quantity<Vec<f64>>,
    pub source: ChargeSource,
}
/// Native CPU inference for a supported NAGL model configuration.
///
/// Load the data-only bundle produced by `benchmarks/openff/scripts/export_model.py`.
/// Chemistry preparation, features, network operations and tensor dimensions are
/// validated before use. No Python, PyTorch or network is used by this loader.
#[derive(Debug)]
pub struct NaglModel {
    identity: ModelIdentity,
    features: Vec<config::Feature>,
    network: network::Network,
    lookup: BTreeMap<String, bundle::Entry>,
    max_lookup_atoms: usize,
    elements: Vec<u8>,
    forbidden: Vec<QueryGraph>,
}
impl NaglModel {
    /// Load a schema-2 bundle, or the original checksum-pinned schema-1 Ash bundle.
    /// Unsupported configurations fail explicitly; no fallback model is selected.
    pub fn load(directory: impl AsRef<Path>) -> Result<Self> {
        let directory = directory.as_ref();
        let bundle =
            bundle::load(directory).map_err(|e| error(format!("{}: {e}", directory.display())))?;
        let forbidden = bundle
            .domain
            .forbidden_patterns
            .iter()
            .map(|p| parse_smarts(p).map_err(error))
            .collect::<Result<Vec<_>>>()?;
        let max_lookup_atoms = bundle
            .lookup
            .values()
            .map(|e| e.charges.len())
            .max()
            .unwrap_or(0);
        Ok(Self {
            identity: bundle.identity,
            features: bundle.config.atom_features,
            network: bundle.network,
            lookup: bundle.lookup,
            max_lookup_atoms,
            elements: bundle.domain.allowed_elements,
            forbidden,
        })
    }
    pub fn identity(&self) -> &ModelIdentity {
        &self.identity
    }
    pub fn lookup_entry_count(&self) -> usize {
        self.lookup.len()
    }
    /// Full fixed-H charge lookup identifier, accepting explicit hydrogens.
    /// Isotopes are ignored to match OpenFF Toolkit's charge representation.
    /// This diagnostic is limited to 1024 atoms by the InChI adapter. Charge
    /// assignment avoids InChI when the model's table cannot contain the input.
    pub fn lookup_identifier(&self, molecule: &Molecule) -> Result<String> {
        identity::fixed_h_inchi(&lookup_molecule(molecule)?)
    }
    /// Configured features in molecule atom iteration order and manifest column order.
    pub fn atom_features(&self, molecule: &Molecule) -> Result<Vec<Vec<f32>>> {
        check_atom_limit(molecule)?;
        let m = explicit(molecule)?;
        self.check_domain(&m)?;
        Ok(features::features(&m, &self.features)?
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
        // The model table bounds the size of any possible full-graph hit. Do not
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
            let tags = reference
                .atoms()
                .filter_map(|(_, a)| a.atom_map)
                .collect::<BTreeSet<_>>();
            if reference.atom_count() != entry.charges.len()
                || tags != (1..=entry.charges.len() as u32).collect()
            {
                return Err(error(
                    "lookup charge array requires consecutive unique atom maps",
                ));
            }
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
                ChargeSource::Lookup {
                    inchi,
                    model: self.identity.clone(),
                },
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
            if !self.elements.is_empty() && !self.elements.contains(&atom.element.atomic_number()) {
                return Err(error(format!(
                    "element {} is outside the model domain",
                    atom.element
                )));
            }
        }
        for pattern in &self.forbidden {
            if find_substructure_match(molecule, pattern)
                .map_err(error)?
                .is_some()
            {
                return Err(error("molecule matches a forbidden model domain pattern"));
            }
        }
        Ok(())
    }
    fn infer_prepared(&self, molecule: &Molecule) -> Result<ChargeAssignment> {
        self.check_domain(molecule)?;
        let x = features::features(molecule, &self.features)?;
        let neighbors = features::adjacency(molecule)?;
        let values = self
            .network
            .charges(x, &neighbors, molecule.formal_charge())?;
        normalize_charges(
            values,
            molecule.formal_charge(),
            ChargeSource::Inference {
                model: self.identity.clone(),
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
    fn invalid_charge_outputs_fail_without_publishing() {
        let source = ChargeSource::Inference {
            model: ModelIdentity::new("test".into(), "0".repeat(64)).unwrap(),
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
