//! Inspection hooks for validating a model against a reference
//! implementation. Parameterization does not need them.
use kekule::core::Molecule;

use crate::{identity, nagl::lookup_molecule, ChargeAssignment, NaglModel, Result};

/// The full fixed-H InChI that keys a model's charge lookup table, for an
/// explicit-hydrogen molecule. Isotopes are ignored, as in the OpenFF
/// Toolkit's charge representation.
///
/// The InChI library accepts at most 1,023 atoms. Charge assignment is not
/// limited by this: it computes an identifier only for molecules no larger
/// than the model's largest lookup entry, and model loading rejects lookup
/// entries above 1,023 atoms.
pub fn lookup_identifier(molecule: &Molecule) -> Result<String> {
    identity::fixed_h_inchi(&lookup_molecule(molecule)?)
}

/// The model's input features, one row per atom in molecule atom order and
/// one column per configured feature column.
pub fn atom_features(model: &NaglModel, molecule: &Molecule) -> Result<Vec<Vec<f32>>> {
    model.atom_features(molecule)
}

/// Neural-network charges, bypassing the lookup table.
pub fn infer_charges(model: &NaglModel, molecule: &Molecule) -> Result<ChargeAssignment> {
    model.infer_charges(molecule)
}

/// The number of entries in the model's charge lookup table.
pub fn lookup_entry_count(model: &NaglModel) -> usize {
    model.lookup_entry_count()
}
