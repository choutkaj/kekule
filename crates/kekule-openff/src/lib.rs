//! SMIRNOFF parameter assignment and NAGL partial charges for Kekule.
//!
//! [`ForceField::rosemary`] compiles the bundled Rosemary force field, and
//! `NaglModel::ash` (default `ash` feature) loads the bundled Ash charge model it requires, so a
//! complete OpenFF parameterization needs no external files.
#![cfg_attr(
    feature = "ash",
    doc = r#"
```
use kekule::smiles;
use kekule_openff::{ForceField, NaglModel};

let mut molecule = smiles::to_molecules("CCO")?.remove(0);
molecule.perceive()?;
molecule.add_hydrogens()?;
let atoms = molecule.atom_count();
let model = NaglModel::ash()?;
let parameters = ForceField::rosemary()?.parameterize_molecule(molecule, &model)?;
assert_eq!(parameters.charges().value().len(), atoms);
# Ok::<(), Box<dyn std::error::Error>>(())
```
"#
)]
//!
//! For systems, [`ForceField::parameterize`] takes a shared topology and
//! parameterizes each reusable definition once. Charges come from complete
//! LibraryCharges, then the model's lookup table, then neural inference;
//! [`ChargeMethod::LibraryOnly`] parameterizes library-charged systems without
//! a model. There is no molecule size limit: work grows linearly with the
//! molecule, and bounded searches only stop pathological patterns.
//!
//! Parameterization consumes explicit-hydrogen molecules. It perceives private
//! copies with MDL aromaticity and never changes the input topology. Results
//! retain the exact topology snapshot and use Kekule's canonical units: nm,
//! kJ/mol, radians, and elementary charges. Energies, gradients, and
//! minimization live in `kekule-potentials`.
//!
//! [`ForceField::from_file`] and [`ForceField::from_offxml`] compile custom
//! rules within the supported SMIRNOFF subset, [`ForceField::append`] composes
//! compatible rule sets, and [`NaglModel::load`] loads another exported model
//! bundle, whose identity must match the force field's NAGLCharges handler.
//! Every failure is an [`Error`] with a stable [`ErrorKind`]. See `CONTRACT.md`
//! for the supported subset, assignment semantics, and validation record.
//!
//! The bundled Ash weights are CC BY 4.0 (see the `kekule-openff-ash` crate);
//! disable the default `ash` feature to load models only from files.
#![forbid(unsafe_code)]
#![warn(rustdoc::broken_intra_doc_links)]

mod assignment;
pub mod diagnostics;
mod error;
mod identity;
mod nagl;
mod offxml;
mod parameters;
mod preparation;

pub use assignment::ChargeMethod;
pub use error::{Error, ErrorKind};
pub use nagl::{ChargeAssignment, ChargeSource, ModelIdentity, NaglModel};
pub use offxml::ForceField;
pub use parameters::*;

pub(crate) use error::Result;
pub(crate) use preparation::explicit;

// Compile and run the contract's examples so they cannot drift from the API.
#[cfg(all(doctest, feature = "ash"))]
#[doc = include_str!("../CONTRACT.md")]
struct ContractExamples;

#[cfg(all(doctest, feature = "ash"))]
#[doc = include_str!("../README.md")]
struct ReadmeExamples;

#[cfg(test)]
fn reference_records() -> Vec<serde_json::Value> {
    let input =
        flate2::read::GzDecoder::new(&include_bytes!("../tests/fixtures/audit.json.gz")[..]);
    let report: serde_json::Value = serde_json::from_reader(input).unwrap();
    report["records"].as_array().unwrap().clone()
}
