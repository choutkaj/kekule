//! SMIRNOFF parameter assignment and configurable NAGL charge inference.
//!
//! Parameterization consumes explicit-hydrogen molecular graphs. It perceives
//! private copies using MDL aromaticity and never rewrites the input topology.
//! Results retain the exact topology snapshot and use Kekule's canonical units:
//! nm, kJ/mol, radians, and elementary charges. Energies, gradients, and
//! minimization are evaluated separately by `kekule-potentials`.
//! Use [`ForceField::from_file`] or [`ForceField::from_offxml`] for custom rules
//! within the supported SMIRNOFF subset. [`NaglModel`] loads a compatible model
//! bundle; its checkpoint identity must match the OFFXML charge handler.
//! [`ForceField::append`] composes compatible compiled rule sets. Complete
//! library-charge systems can use [`ForceField::parameterize_without_nagl`]
//! without a model bundle. Explicit distance constraints may join nonbonded
//! atoms within one molecule; they do not change chemical connectivity.
//! Rosemary and Ash remain the independently validated presets.
//!
//! ```no_run
//! use kekule::{hydrogens, smiles};
//! use kekule_openff::{ForceField, NaglModel};
//!
//! let mut molecule = smiles::to_molecules("CCO")?.remove(0);
//! molecule.perceive()?;
//! hydrogens::add_hydrogens(&mut molecule)?;
//! let model = NaglModel::load("path/to/exported-ash")?;
//! let parameters = ForceField::rosemary()?.parameterize_molecule(&molecule, &model)?;
//! assert_eq!(parameters.charges().value().len(), molecule.atom_count());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
#![forbid(unsafe_code)]
#![warn(rustdoc::broken_intra_doc_links)]

mod assignment;
mod identity;
mod nagl;
mod offxml;
mod parameters;

pub use nagl::{ChargeAssignment, ChargeSource, ModelIdentity, NaglModel};
pub use offxml::ForceField;
pub use parameters::*;

/// An explicit input, format, coverage, or resource-limit failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error(String);
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for Error {}
pub(crate) type Result<T> = std::result::Result<T, Error>;
pub(crate) fn error(e: impl std::fmt::Display) -> Error {
    Error(e.to_string())
}

pub(crate) fn explicit(molecule: &kekule::core::Molecule) -> Result<kekule::core::Molecule> {
    let mut copy = preparation::normalize_valence(molecule)?;
    copy.perceive().map_err(error)?;
    for (id, atom) in copy.atoms() {
        if copy.implicit_hydrogens(id).map_err(error)? != Some(0) {
            return Err(error(format!(
                "atom {id} has implicit hydrogens; expand hydrogens before parameterization"
            )));
        }
        if atom.radical.is_some() {
            return Err(error(
                "radicals are outside the supported parameterization domain",
            ));
        }
    }
    Ok(copy)
}

mod preparation;

#[cfg(test)]
fn reference_records() -> Vec<serde_json::Value> {
    let input =
        flate2::read::GzDecoder::new(&include_bytes!("../tests/fixtures/audit.json.gz")[..]);
    let report: serde_json::Value = serde_json::from_reader(input).unwrap();
    report["records"].as_array().unwrap().clone()
}
