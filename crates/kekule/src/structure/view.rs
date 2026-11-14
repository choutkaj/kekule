use std::sync::Arc;

use crate::math::{PeriodicCell, Point3};
use crate::properties::GeometryProperties;
use crate::topology::{InstanceAtomId, Topology, TopologyAtomIndex};
use crate::units::Quantity;

use super::scalars::{b_factor_quantities, b_factor_quantity, AtomScalars};
use super::state::is_bound;
use super::{Positions, SetContext, Structure};

pub(crate) mod sealed {
    /// Restricts [`super::AsStructureView`] to Kekule's structure owners and views.
    pub trait Sealed {}
}

/// Anything that can be borrowed as a [`StructureView`] without copying: structures,
/// views, geometries, and trajectory frames.
///
/// Each implementor also has an inherent `as_structure_view` method, so the trait
/// only needs importing for generic code. This trait is sealed; other
/// structure-like types, such as a [`crate::traj::FrameBuffer`], are passed
/// through their own `as_structure_view()`.
pub trait AsStructureView: sealed::Sealed {
    /// Borrows the structure as a [`StructureView`] without copying.
    fn as_structure_view(&self) -> StructureView<'_>;
}

impl sealed::Sealed for Structure {}
impl sealed::Sealed for StructureView<'_> {}

impl AsStructureView for Structure {
    fn as_structure_view(&self) -> StructureView<'_> {
        Structure::as_structure_view(self)
    }
}

impl AsStructureView for StructureView<'_> {
    fn as_structure_view(&self) -> StructureView<'_> {
        *self
    }
}

impl<'a> From<&'a Structure> for StructureView<'a> {
    fn from(structure: &'a Structure) -> Self {
        structure.as_structure_view()
    }
}

/// What a [`StructureView`] borrows: a structure, or one geometry of a structure
/// set or trajectory frame source with what its source records about it.
///
/// Kernels that only need coordinates ignore it. Writers and reductions read
/// it so that a geometry's place, velocities, and forces are
/// written, rejected, or deliberately omitted, never silently dropped.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum StructureKind<'a> {
    /// A [`Structure`], or a view materialized as one: a structure held on its
    /// own, with no set around it.
    Standalone,
    /// A geometry of a [`StructureSet`](super::StructureSet), a trajectory frame,
    /// or a frame buffer, with its [`SetContext`].
    InSet(SetContext<'a>),
}

/// A borrowed topology plus the geometry state bound to it: the common
/// input contract of coordinate-dependent kernels.
///
/// Structures, geometries, trajectory frames, and frame buffers all borrow as a
/// view, and [`Self::kind`] says whether it is a structure or a geometry.
/// Accessors return borrows with the view's full lifetime. A view of a
/// trajectory frame has no occupancies or B-factors. Views preserve exact
/// topology identity.
#[derive(Debug, Clone, Copy)]
pub struct StructureView<'a> {
    topology: &'a Arc<Topology>,
    positions: &'a Positions,
    cell: Option<&'a PeriodicCell>,
    occupancies: &'a AtomScalars,
    b_factors: &'a AtomScalars,
    properties: &'a GeometryProperties,
    kind: StructureKind<'a>,
}

impl<'a> StructureView<'a> {
    /// Only owners that already enforce binding may use this.
    pub(crate) fn bound(
        topology: &'a Arc<Topology>,
        positions: &'a Positions,
        cell: Option<&'a PeriodicCell>,
        occupancies: &'a AtomScalars,
        b_factors: &'a AtomScalars,
        properties: &'a GeometryProperties,
        kind: StructureKind<'a>,
    ) -> Self {
        debug_assert!(is_bound(positions.len(), properties, topology));
        Self {
            topology,
            positions,
            cell,
            occupancies,
            b_factors,
            properties,
            kind,
        }
    }

    /// Whether this view borrows a structure or a geometry, with the geometry's
    /// place, velocities, and forces.
    ///
    /// ```
    /// use kekule::{smiles, structure::{Structure, Positions, StructureKind}};
    /// let structure = Structure::new(smiles::to_topology("O")?, Positions::zeros(1))?;
    /// assert_eq!(structure.as_structure_view().kind(), StructureKind::Standalone);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub const fn kind(self) -> StructureKind<'a> {
        self.kind
    }

    /// The topology the viewed structure belongs to.
    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub(crate) const fn topology_arc(self) -> &'a Arc<Topology> {
        self.topology
    }

    /// A shared handle to the bound topology.
    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    /// Number of atoms, equal to the topology's atom count.
    pub fn atom_count(self) -> usize {
        self.positions.len()
    }

    /// Positions in dense topology atom order.
    pub fn positions(self) -> &'a Positions {
        self.positions
    }

    /// The periodic cell, or `None` for a non-periodic structure.
    pub fn cell(self) -> Option<&'a PeriodicCell> {
        self.cell
    }

    /// Structure annotations: owner values plus one row per dense atom and bond.
    pub fn properties(self) -> &'a GeometryProperties {
        self.properties
    }

    /// Occupancy of one atom; `None` when it was never recorded, which is
    /// always the case for a trajectory frame, or when `atom` is out of range.
    pub fn occupancy(self, atom: TopologyAtomIndex) -> Option<f64> {
        self.occupancies.get(atom)
    }

    /// Every occupancy in dense atom order, or `None` when none was recorded.
    pub fn occupancies(self) -> Option<&'a [Option<f64>]> {
        self.occupancies.values()
    }

    /// Isotropic B-factor of one atom; `None` when it was never recorded,
    /// which is always the case for a trajectory frame, or when `atom` is out
    /// of range.
    pub fn b_factor(self, atom: TopologyAtomIndex) -> Option<Quantity<f64>> {
        b_factor_quantity(self.b_factors.get(atom))
    }

    /// Every B-factor in dense atom order and square nanometres, or `None`
    /// when none was recorded.
    pub fn b_factors(self) -> Option<Quantity<&'a [Option<f64>]>> {
        b_factor_quantities(self.b_factors)
    }

    pub(crate) fn occupancy_column(self) -> &'a AtomScalars {
        self.occupancies
    }

    pub(crate) fn b_factor_column(self) -> &'a AtomScalars {
        self.b_factors
    }

    /// The position of one atom, or `None` when `atom` is not in the topology.
    pub fn position(self, atom: InstanceAtomId) -> Option<Quantity<Point3>> {
        self.positions.get(self.topology.atom_index(atom)?.index())
    }

    /// One atom with its topology identity and geometry state.
    pub fn atom(self, atom: InstanceAtomId) -> Option<StructureAtomView<'a>> {
        let atom = self.topology.atom(atom)?;
        Some(StructureAtomView {
            atom,
            structure: self,
        })
    }

    /// One atom by dense index.
    pub fn atom_at(self, index: TopologyAtomIndex) -> Option<StructureAtomView<'a>> {
        let atom = self.topology.atom_at(index)?;
        Some(StructureAtomView {
            atom,
            structure: self,
        })
    }

    /// Atoms with their geometry state, in dense order.
    pub fn atoms(self) -> impl ExactSizeIterator<Item = StructureAtomView<'a>> + 'a {
        self.topology.atoms().map(move |atom| StructureAtomView {
            atom,
            structure: self,
        })
    }

    /// Materializes this borrowed structure as an owned structure, dropping the
    /// set context reported by [`Self::kind`]. A view of a
    /// trajectory frame yields a structure without occupancies or B-factors.
    pub fn to_structure(self) -> Structure {
        Structure {
            topology: Arc::clone(self.topology),
            positions: self.positions.clone(),
            cell: self.cell.copied(),
            occupancies: self.occupancies.clone(),
            b_factors: self.b_factors.clone(),
            properties: self.properties.clone(),
        }
    }
}

/// One atom of a structure: its topology view (through `Deref`) plus the
/// coordinate-dependent state stored for it.
///
/// ```
/// use kekule::{smiles, structure::{Structure, Positions}};
/// let structure = Structure::new(smiles::to_topology("CO")?, Positions::zeros(2))?;
/// let oxygen = structure.atoms().find(|atom| atom.element().symbol() == "O").unwrap();
/// assert_eq!(oxygen.occupancy(), None);
/// assert_eq!(oxygen.position().into_value(), kekule::math::Point3::origin());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct StructureAtomView<'a> {
    atom: crate::topology::AtomView<'a>,
    structure: StructureView<'a>,
}

impl<'a> std::ops::Deref for StructureAtomView<'a> {
    type Target = crate::topology::AtomView<'a>;

    fn deref(&self) -> &Self::Target {
        &self.atom
    }
}

impl<'a> StructureAtomView<'a> {
    /// The topology view with the structure's full lifetime.
    pub const fn topology_atom(self) -> crate::topology::AtomView<'a> {
        self.atom
    }

    /// The atom's stored position in the canonical length unit (nm).
    pub fn position(self) -> Quantity<Point3> {
        self.structure
            .positions()
            .get(self.atom.index().index())
            .expect("bound structures cover every atom")
    }

    /// The atom's occupancy, or `None` when it was never recorded.
    pub fn occupancy(self) -> Option<f64> {
        self.structure.occupancy(self.atom.index())
    }

    /// The atom's isotropic B-factor in nm², or `None` when it was never
    /// recorded.
    pub fn b_factor(self) -> Option<Quantity<f64>> {
        self.structure.b_factor(self.atom.index())
    }

    /// One geometry annotation of this atom. Static annotations are
    /// available through [`crate::topology::AtomView::property`].
    pub fn geometry_property(
        self,
        key: &crate::properties::PropertyKey,
    ) -> Option<crate::properties::PropertyValue> {
        self.structure
            .properties()
            .atoms()
            .value(key, self.atom.index())
            .expect("bound structures cover every atom")
    }
}
