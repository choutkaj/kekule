//! Inspection hooks for validating a model against a reference
//! implementation. Parameterization does not need them.
use kekule::core::Molecule;

use crate::{ChargeAssignment, NaglModel, Result};

/// The stored identifier of the lookup entry selected for an explicit-hydrogen
/// molecule, or `None` for a lookup miss. This is the entry upstream's fixed-H
/// InChI lookup selects; no InChI is computed. Isotopes are ignored, as in the
/// OpenFF Toolkit's charge representation.
pub fn lookup_key<'a>(model: &'a NaglModel, molecule: &Molecule) -> Result<Option<&'a str>> {
    model.lookup_key(molecule)
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
