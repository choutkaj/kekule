//! Checks shared by the geometry owners: row indices, topology binding, and
//! property row counts.

use crate::diagnostics::{Countable, IndexOutOfRange, LengthMismatch};
use crate::properties::GeometryProperties;
use crate::topology::{Topology, TopologyAtomIndex};

use super::{DenseArrayError, GeometryError};

/// The row of one dense atom index in an owner with `len` atoms.
pub(crate) fn atom_row(len: usize, atom: TopologyAtomIndex) -> Result<usize, GeometryError> {
    let index = atom.index();
    if index >= len {
        return Err(GeometryError::Array(Box::new(
            DenseArrayError::IndexOutOfRange(IndexOutOfRange::new(Countable::Row, index, len)),
        )));
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
        return Err(GeometryError::Array(Box::new(
            DenseArrayError::CountMismatch(LengthMismatch::new(
                Countable::Row,
                topology.atom_count(),
                atoms,
            )),
        )));
    }
    properties
        .validate_bond_rows(topology.bond_count())
        .map_err(|_| {
            GeometryError::BondCountMismatch(LengthMismatch::new(
                Countable::Row,
                topology.bond_count(),
                properties.bonds().len(),
            ))
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
        return Err(GeometryError::Array(Box::new(
            DenseArrayError::CountMismatch(LengthMismatch::new(
                Countable::Row,
                atoms,
                replacement.atoms().len(),
            )),
        )));
    }
    if replacement.bonds().len() != current.bonds().len() {
        return Err(GeometryError::BondCountMismatch(LengthMismatch::new(
            Countable::Row,
            current.bonds().len(),
            replacement.bonds().len(),
        )));
    }
    Ok(())
}
