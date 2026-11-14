//! Checks shared by the geometry owners: row indices, topology binding, and
//! property row counts.

use crate::properties::GeometryProperties;
use crate::topology::{Topology, TopologyAtomIndex};

use super::{DenseArrayError, GeometryError};

/// The row of one dense atom index in an owner with `len` atoms.
pub(crate) fn atom_row(len: usize, atom: TopologyAtomIndex) -> Result<usize, GeometryError> {
    let index = atom.index();
    if index >= len {
        return Err(GeometryError::Array(DenseArrayError::IndexOutOfRange {
            index,
            len,
        }));
    }
    Ok(index)
}

/// Checks that geometry rows can be bound to `topology` without changing
/// anything: one position per atom, and bond rows either matching or not yet
/// allocated.
pub(crate) fn validate_binding(
    atoms: usize,
    properties: &GeometryProperties,
    topology: &Topology,
) -> Result<(), GeometryError> {
    if atoms != topology.atom_count() {
        return Err(GeometryError::Array(DenseArrayError::CountMismatch {
            expected: topology.atom_count(),
            actual: atoms,
        }));
    }
    properties
        .validate_bond_rows(topology.bond_count())
        .map_err(|_| GeometryError::BondCountMismatch {
            expected: topology.bond_count(),
            actual: properties.bonds().len(),
        })
}

/// Whether geometry rows already cover every atom and bond of `topology`.
pub(crate) fn is_bound(atoms: usize, properties: &GeometryProperties, topology: &Topology) -> bool {
    atoms == topology.atom_count() && properties.bonds().len() == topology.bond_count()
}

/// Checks that a replacement keeps the current atom and bond row counts.
pub(crate) fn check_property_rows(
    atoms: usize,
    current: &GeometryProperties,
    replacement: &GeometryProperties,
) -> Result<(), GeometryError> {
    if replacement.atoms().len() != atoms {
        return Err(GeometryError::Array(DenseArrayError::CountMismatch {
            expected: atoms,
            actual: replacement.atoms().len(),
        }));
    }
    if replacement.bonds().len() != current.bonds().len() {
        return Err(GeometryError::BondCountMismatch {
            expected: current.bonds().len(),
            actual: replacement.bonds().len(),
        });
    }
    Ok(())
}
