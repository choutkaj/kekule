mod bundle;
mod config;
mod features;
mod network;
use crate::assignment::options;
use crate::error::atom_list;
use crate::{explicit, identity, Error, ErrorKind, Result};
use kekule::{
    molecule::{Element, Molecule},
    query::QueryGraph,
    smarts,
    substructure::PreparedTarget,
    units::{Quantity, ELEMENTARY_CHARGE},
};
use std::{collections::BTreeSet, path::Path};

/// Checkpoint identity declared by an OFFXML handler or exported model bundle.
/// The checksum identifies the original checkpoint, not the converted weights file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelIdentity {
    model_file: String,
    checkpoint_sha256: String,
}
impl ModelIdentity {
    /// A nonempty model file name and a 64-digit hexadecimal SHA-256, or
    /// `None`; callers report the failure in their own terms.
    pub(crate) fn new(model_file: String, checkpoint_sha256: String) -> Option<Self> {
        let valid = !model_file.trim().is_empty()
            && checkpoint_sha256.len() == 64
            && checkpoint_sha256.bytes().all(|b| b.is_ascii_hexdigit());
        valid.then(|| Self {
            model_file,
            checkpoint_sha256: checkpoint_sha256.to_ascii_lowercase(),
        })
    }

    /// The model checkpoint's file name, never blank.
    pub fn model_file(&self) -> &str {
        &self.model_file
    }

    /// The original checkpoint's SHA-256 as 64 lowercase hexadecimal digits.
    pub fn checkpoint_sha256(&self) -> &str {
        &self.checkpoint_sha256
    }
}

// OpenFF Toolkit's Molecule representation does not retain isotopic masses.
// NAGL lookup selection and graph correspondence must use that same
// representation, while the caller's Kekule graph retains them.
pub(crate) fn lookup_molecule(input: &Molecule) -> Result<Molecule> {
    if !input.atoms().any(|(_, a)| a.isotope().is_some()) {
        return explicit(input);
    }
    let mut editor = input.to_editor();
    for (id, atom) in input.atoms() {
        if atom.isotope().is_some() {
            editor
                .atom_mut(id)
                .expect("atom IDs come from this molecule")
                .set_isotope(None);
        }
    }
    explicit(&editor.finish().map_err(Error::chemistry)?)
}

fn model_error(detail: impl std::fmt::Display) -> Error {
    Error::new(ErrorKind::Model, detail)
}

/// Provenance of one molecule's assigned charges.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ChargeSource {
    /// Complete LibraryCharges coverage.
    Library {
        /// The distinct IDs of the LibraryCharge rules that supplied the
        /// charges, sorted; a rule without an `id` contributes an empty string.
        parameter_ids: Vec<String>,
    },
    /// An exact entry of the NAGL model's lookup table, corrected to
    /// conserve the formal charge.
    Lookup {
        /// The InChI string of the selected lookup entry, as stored in the bundle.
        inchi: String,
        /// The model whose lookup table supplied the charges.
        model: ModelIdentity,
    },
    /// Neural inference, corrected to conserve the formal charge.
    Inference {
        /// The model that inferred the charges.
        model: ModelIdentity,
    },
}
/// Partial charges for one molecule with their provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct ChargeAssignment {
    pub(crate) charges: Quantity<Vec<f64>>,
    pub(crate) source: ChargeSource,
}

impl ChargeAssignment {
    /// One finite charge per atom in e, in the molecule's atom order. They sum
    /// to the formal charge (library charges within 1e-6 e).
    pub fn charges(&self) -> &Quantity<Vec<f64>> {
        &self.charges
    }

    /// Where the charges came from.
    pub fn source(&self) -> &ChargeSource {
        &self.source
    }
}
/// Native CPU inference for a supported NAGL model configuration.
///
/// Load the data-only bundle produced by `tools/openff/export_model.py`.
/// Chemistry preparation, features, network operations and tensor dimensions are
/// validated before use. No Python, PyTorch or network is used by this loader.
#[derive(Debug)]
pub struct NaglModel {
    identity: ModelIdentity,
    features: Vec<config::Feature>,
    network: network::Network,
    lookup: Vec<bundle::Entry>,
    index: identity::LookupIndex,
    max_lookup_atoms: usize,
    elements: Vec<u8>,
    forbidden: Vec<QueryGraph>,
}
impl NaglModel {
    /// Load a schema-2 bundle, or the original checksum-pinned schema-1 Ash bundle.
    /// Unsupported configurations fail explicitly; no fallback model is selected.
    ///
    /// # Errors
    ///
    /// [`ErrorKind::Io`] if a bundle file cannot be read, and
    /// [`ErrorKind::Model`] if the bundle is malformed, corrupt or unsupported.
    pub fn load(directory: impl AsRef<Path>) -> Result<Self> {
        Self::from_bundle(bundle::load(directory.as_ref())?)
    }

    /// The bundled Ash model, `openff-gnn-am1bcc-1.0.0`, required by
    /// [`crate::ForceField::rosemary`].
    ///
    /// The model ships in the `kekule-openff-ash` crate (default `ash`
    /// feature) and is validated exactly like [`Self::load`]. Decoding takes a
    /// fraction of a second, so load it once and reuse it. The weights are
    /// CC BY 4.0; retain their attribution when redistributing them.
    ///
    /// # Errors
    ///
    /// [`ErrorKind::Model`] if the bundled data fail validation, which
    /// indicates a corrupted build.
    ///
    /// ```
    /// let model = kekule_openff::NaglModel::ash()?;
    /// let rosemary = kekule_openff::ForceField::rosemary()?;
    /// assert_eq!(rosemary.charge_model(), Some(model.identity()));
    /// # Ok::<(), kekule_openff::Error>(())
    /// ```
    #[cfg(feature = "ash")]
    #[cfg_attr(docsrs, doc(cfg(feature = "ash")))]
    pub fn ash() -> Result<Self> {
        use std::io::Read;
        let mut planes = Vec::new();
        flate2::read::ZlibDecoder::new(kekule_openff_ash::WEIGHT_PLANES_ZLIB)
            .take(bundle::WEIGHTS_LIMIT + 1)
            .read_to_end(&mut planes)
            .map_err(|e| Error::wrap(ErrorKind::Model, e))?;
        if planes.len() % 4 != 0 {
            return Err(model_error(
                "bundled Ash weights are not float32 byte planes",
            ));
        }
        let count = planes.len() / 4;
        let mut weights = vec![0; planes.len()];
        for (plane, bytes) in planes.chunks_exact(count.max(1)).enumerate() {
            for (value, byte) in bytes.iter().enumerate() {
                weights[value * 4 + plane] = *byte;
            }
        }
        let model = Self::from_bundle(bundle::parse(kekule_openff_ash::MANIFEST, &weights)?)?;
        if model.identity().model_file() != kekule_openff_ash::MODEL_FILE
            || model.identity().checkpoint_sha256() != kekule_openff_ash::CHECKPOINT_SHA256
        {
            return Err(model_error(
                "bundled Ash manifest declares another checkpoint",
            ));
        }
        Ok(model)
    }

    fn from_bundle(bundle: bundle::Bundle) -> Result<Self> {
        let forbidden = bundle
            .domain
            .forbidden_patterns
            .iter()
            .map(|p| {
                smarts::parse_str(p, Default::default())
                    .map_err(|e| Error::wrap(ErrorKind::Model, e))
            })
            .collect::<Result<Vec<_>>>()?;
        let lookup: Vec<_> = bundle.lookup.into_values().collect();
        let max_lookup_atoms = lookup.iter().map(|e| e.charges.len()).max().unwrap_or(0);
        // Entries input preparation rejects can never be selected. Ash's ten
        // pairs of entries sharing one molecule carry identical charges.
        let graphs = lookup.iter().enumerate().filter_map(|(position, entry)| {
            let mut components = kekule::smiles::to_molecules(&entry.mapped_smiles).ok()?;
            (components.len() == 1)
                .then(|| explicit(&components.remove(0)).ok())
                .flatten()
                .map(|graph| (graph, position))
        });
        let index = identity::LookupIndex::new(graphs)?;
        Ok(Self {
            identity: bundle.identity,
            features: bundle.config.atom_features,
            network: bundle.network,
            lookup,
            index,
            max_lookup_atoms,
            elements: bundle.domain.allowed_elements,
            forbidden,
        })
    }

    /// The checkpoint identity the bundle declares.
    pub fn identity(&self) -> &ModelIdentity {
        &self.identity
    }

    /// The elements the model charges, in ascending atomic number: those its
    /// domain allows that its element features also encode.
    ///
    /// A molecule with any other element fails [`Self::assign_charges`] with
    /// [`ErrorKind::UnsupportedMolecule`], so a host can refuse it before
    /// parameterizing. The model's forbidden patterns, valence, and
    /// connectivity limits can still reject molecules of these elements.
    ///
    /// ```
    /// use kekule::molecule::Element;
    /// let ash = kekule_openff::NaglModel::ash()?;
    /// assert!(ash.supports_element(Element::C));
    /// assert!(!ash.elements().contains(&Element::from_symbol("Mg").unwrap()));
    /// # Ok::<(), kekule_openff::Error>(())
    /// ```
    pub fn elements(&self) -> Vec<Element> {
        (1..=118)
            .filter_map(Element::from_atomic_number)
            .filter(|&element| self.supports_element(element))
            .collect()
    }

    /// Whether [`Self::elements`] contains `element`.
    pub fn supports_element(&self, element: Element) -> bool {
        (self.elements.is_empty() || self.elements.contains(&element.atomic_number()))
            && features::encodes_element(&self.features, element)
    }
    #[cfg(test)]
    pub(crate) fn lookup_entry_count(&self) -> usize {
        self.lookup.len()
    }
    #[cfg(test)]
    pub(crate) fn lookup_key(&self, molecule: &Molecule) -> Result<Option<&str>> {
        let molecule = lookup_molecule(molecule)?;
        if molecule.atom_count() > self.max_lookup_atoms {
            return Ok(None);
        }
        Ok(self
            .index
            .select(&molecule)?
            .map(|position| self.lookup[position].inchi.as_str()))
    }
    #[cfg(test)]
    pub(crate) fn atom_features(&self, molecule: &Molecule) -> Result<Vec<Vec<f32>>> {
        let m = explicit(molecule)?;
        self.check_domain(&m)?;
        Ok(features::features(&m, &self.features)?
            .rows()
            .into_iter()
            .map(|r| r.to_vec())
            .collect())
    }
    /// Charges for one explicit-hydrogen molecule: the model's lookup table
    /// first, otherwise neural inference, corrected to conserve the formal
    /// charge. Charges follow the molecule's atom order. This ignores
    /// LibraryCharges; use [`crate::ForceField::parameterize`] for force-field
    /// charge precedence. There is no molecule size limit.
    ///
    /// # Errors
    ///
    /// - [`ErrorKind::UnsupportedMolecule`] if the molecule has implicit
    ///   hydrogens or radicals, or is outside the model's element or pattern
    ///   domain.
    /// - [`ErrorKind::Chemistry`] if perceiving the private copy fails.
    /// - [`ErrorKind::Identity`] or [`ErrorKind::Model`] if a lookup hit
    ///   cannot be mapped onto its entry.
    /// - [`ErrorKind::ResourceLimit`] if a bounded search exceeds its limit.
    /// - [`ErrorKind::Charges`] if charge equilibration is singular or yields
    ///   nonfinite charges.
    pub fn assign_charges(&self, molecule: &Molecule) -> Result<ChargeAssignment> {
        let molecule = lookup_molecule(molecule)?;
        // A full-graph lookup hit must have one charge for every input atom.
        if molecule.atom_count() > self.max_lookup_atoms {
            return self.infer_prepared(&molecule);
        }
        if let Some(position) = self.index.select(&molecule)? {
            let entry = &self.lookup[position];
            let mut components = kekule::smiles::to_molecules(&entry.mapped_smiles)
                .map_err(|e| Error::wrap(ErrorKind::Model, e))?;
            if components.len() != 1 {
                return Err(model_error("lookup entry is disconnected"));
            }
            let reference = explicit(&components.remove(0))?;
            let tags = reference
                .atoms()
                .filter_map(|(_, a)| a.atom_map().map(std::num::NonZeroU32::get))
                .collect::<BTreeSet<_>>();
            if reference.atom_count() != entry.charges.len()
                || tags != (1..=entry.charges.len() as u32).collect()
            {
                return Err(model_error(
                    "lookup charge array requires consecutive unique atom maps",
                ));
            }
            let mapping = identity::mapping(&molecule, &reference)?;
            let mut values = Vec::new();
            for atom in mapping {
                let map = reference
                    .atom(atom)
                    .expect("atom IDs come from this molecule")
                    .atom_map()
                    .ok_or_else(|| model_error("unmapped lookup atom"))?;
                let q = entry
                    .charges
                    .get(map.get() as usize - 1)
                    .ok_or_else(|| model_error("lookup charge array mismatch"))?;
                values.push(f64::from(*q));
            }
            return normalize_charges(
                values,
                molecule.formal_charge(),
                ChargeSource::Lookup {
                    inchi: entry.inchi.clone(),
                    model: self.identity.clone(),
                },
            );
        }
        self.infer_prepared(&molecule)
    }
    #[cfg(test)]
    pub(crate) fn infer_charges(&self, molecule: &Molecule) -> Result<ChargeAssignment> {
        self.infer_prepared(&explicit(molecule)?)
    }
    fn check_domain(&self, molecule: &Molecule) -> Result<()> {
        let outside: Vec<_> = molecule
            .atoms()
            .filter(|(_, atom)| !self.supports_element(atom.element()))
            .collect();
        if !outside.is_empty() {
            let elements = outside
                .iter()
                .map(|(_, atom)| atom.element())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .map(|element| element.symbol())
                .collect::<Vec<_>>();
            let atoms: Vec<_> = outside.iter().map(|(id, _)| *id).collect();
            return Err(Error::new(
                ErrorKind::UnsupportedMolecule,
                format!(
                    "{} {} {} outside the model domain: {}",
                    if elements.len() == 1 {
                        "element"
                    } else {
                        "elements"
                    },
                    elements.join(", "),
                    if elements.len() == 1 { "is" } else { "are" },
                    atom_list(&atoms)
                ),
            )
            .with_atoms(atoms));
        }
        let target = PreparedTarget::new(molecule);
        for pattern in &self.forbidden {
            let Some(first) = target
                .first(pattern, Default::default())
                .map_err(Error::chemistry)?
            else {
                continue;
            };
            // Name every matched atom when the search stays within bounds.
            let atoms: BTreeSet<_> = match target.all(pattern, options(molecule.atom_count())) {
                Ok(matches) => matches
                    .iter()
                    .flat_map(|matched| matched.atoms().iter().copied())
                    .collect(),
                Err(_) => first.atoms().iter().copied().collect(),
            };
            let atoms: Vec<_> = atoms.into_iter().collect();
            return Err(Error::new(
                ErrorKind::UnsupportedMolecule,
                format!(
                    "molecule matches a forbidden model domain pattern at {}",
                    atom_list(&atoms)
                ),
            )
            .with_atoms(atoms));
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
        return Err(Error::new(
            ErrorKind::Charges,
            "NAGL produced nonfinite or empty charges",
        ));
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
    fn lookup_preparation_drops_isotopes_without_mutating_input() {
        // PubChem CID 143783 from the locked pubchem-100k corpus. The OpenFF
        // Toolkit's charge representation retains no isotopic masses.
        let mut input = kekule::smiles::to_molecules("[2H]C([2H])([2H])C([2H])([2H])SCC")
            .unwrap()
            .remove(0);
        input.perceive().unwrap();
        input.add_hydrogens(Default::default()).unwrap();
        let before = input.clone();
        let prepared = lookup_molecule(&input).unwrap();
        assert!(prepared.atoms().all(|(_, a)| a.isotope().is_none()));
        assert_eq!(input, before);
        assert_eq!(
            prepared.atom_ids().collect::<Vec<_>>(),
            input.atom_ids().collect::<Vec<_>>()
        );
    }

    #[cfg(feature = "ash")]
    #[test]
    fn lookup_selects_exact_entry_molecules_ignoring_isotopes() {
        // Heavy water selects the same Ash entry as water, as upstream does.
        let model = NaglModel::ash().unwrap();
        let select = |smiles: &str| {
            let m = kekule::smiles::to_molecules(smiles).unwrap().remove(0);
            model.lookup_key(&m).unwrap().map(str::to_owned)
        };
        assert_eq!(select("[H]O[H]").as_deref(), Some("InChI=1/H2O/h1H2"));
        assert_eq!(select("[2H]O[2H]"), select("[H]O[H]"));
        assert_eq!(
            select("[H]C([H])([H])[H]").as_deref(),
            Some("InChI=1/CH4/h1H4")
        );
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
}
