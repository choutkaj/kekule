use super::*;
use crate::structure::ModelView;
use crate::topology::{AppendMapping, Topology};
use std::collections::BTreeSet;

impl ModelEditor {
    /// Appends a complete borrowed model, including its existing coordinates.
    ///
    /// Accepts `&Model` or a [`ModelView`] from an ensemble member or trajectory
    /// frame without first creating an owned model. Positions are used as supplied:
    /// the caller is responsible for their placement in the destination coordinate
    /// system. No fitting, imaging, bond inference or geometry generation is done.
    ///
    /// Molecule definitions retain their reuse within this import, represented
    /// stereo, perception, classification and definition properties. Instances,
    /// atoms, bonds and hierarchy receive distinct editing identities. Chains are
    /// kept separate even when their labels coincide; labels and author identifiers
    /// are preserved. All entity property columns at definition, topology and model
    /// scope are transferred, including occupancy, B factors and missing values.
    /// Incompatible property types or units reject the entire append.
    ///
    /// Changed destination topology/model owner properties are cleared, and source
    /// topology/model owner properties are not imported: they do not automatically
    /// describe the combined system. [`ModelAppend::report`] lists these keys.
    /// Other owners and all source data remain unchanged. Transferring an arbitrary
    /// annotation does not guarantee its scientific validity after later edits.
    ///
    /// A source without a cell uses the destination cell. A cell-less empty editor
    /// adopts the source cell. Otherwise a periodic source requires equal periodic
    /// axes and canonical cell vectors within floating-point conversion roundoff
    /// (16 machine epsilons times the largest vector component). Callers can set
    /// the intended cell before retrying. Failure leaves the entire editor unchanged.
    ///
    /// ```
    /// use kekule::{smiles, structure::{Model, Positions}, geometry::Point3,
    ///     units::{Quantity, ANGSTROM}};
    /// let molecule = smiles::to_molecules("CO")?.pop().unwrap();
    /// let ligand = Model::from_molecule(&molecule, &Positions::new(Quantity::new(
    ///     vec![Point3::origin(), Point3::new(1.4, 0.0, 0.0)], ANGSTROM))?)?;
    /// let mut editor = ligand.edit();
    /// let appended = editor.append_model(&ligand)?;
    /// let source_atom = ligand.topology().atom_ids()[0];
    /// let draft_atom = appended.atom(source_atom)?;
    /// assert_eq!(editor.position(draft_atom)?, ligand.position(source_atom)?);
    /// let result = editor.finish_with_correspondence()?;
    /// let mapping = appended.published(&result)?;
    /// assert!(mapping.atom(source_atom).is_some());
    /// assert_eq!(result.model().atom_count(), 4);
    /// // Use editor.finish() instead when only the completed Model is needed.
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn append_model<'a>(
        &mut self,
        source: impl Into<ModelView<'a>>,
    ) -> Result<ModelAppend, ModelEditError> {
        let source = source.into();
        let cell = match (self.cell, source.cell().copied()) {
            (None, cell) if self.is_empty() => cell,
            (cell, None) => cell,
            (Some(left), Some(right)) if compatible_cells(left, right) => Some(left),
            _ => return Err(ModelEditError::IncompatibleAppendCell),
        };
        let report = ModelAppendReport {
            cleared_topology_properties: owner_keys(self.topology.properties()),
            cleared_model_properties: owner_keys(&self.properties),
            omitted_topology_properties: owner_keys(source.topology().properties()),
            omitted_model_properties: owner_keys(source.properties()),
        };
        let mut staged = self.clone();
        staged
            .positions
            .try_reserve(source.positions().len())
            .map_err(|_| ModelEditError::CapacityOverflow)?;
        let mapping = staged.topology.append_topology(source.shared_topology())?;
        staged
            .positions
            .extend_from_slice(source.positions().values().value());
        staged
            .properties
            .resize_atoms(staged.topology.atom_slot_count());
        staged
            .properties
            .resize_bonds(staged.topology.bond_slot_count());
        let atom_rows = source
            .topology()
            .atom_ids()
            .iter()
            .map(|id| staged.topology.atom_slot(mapping.atoms[id]))
            .collect::<Result<Vec<_>, _>>()?;
        let bond_rows = source
            .topology()
            .bond_ids()
            .iter()
            .map(|id| staged.topology.bond_slot(mapping.bonds[id]))
            .collect::<Result<Vec<_>, _>>()?;
        staged
            .properties
            .atoms_mut()
            .copy_rows_from(source.atom_properties(), &atom_rows)
            .map_err(|error| ModelEditError::AppendProperty {
                domain: "model atom",
                error: Box::new(error),
            })?;
        staged
            .properties
            .bonds_mut()
            .copy_rows_from(source.bond_properties(), &bond_rows)
            .map_err(|error| ModelEditError::AppendProperty {
                domain: "model bond",
                error: Box::new(error),
            })?;
        staged.properties.clear_owner();
        staged.cell = cell;
        *self = staged;
        Ok(ModelAppend { mapping, report })
    }
}

fn compatible_cells(left: PeriodicCell, right: PeriodicCell) -> bool {
    if left.periodic_axes() != right.periodic_axes() {
        return false;
    }
    let left = left.vectors().into_value();
    let right = right.vectors().into_value();
    let scale = left
        .iter()
        .chain(&right)
        .flat_map(|v| [v.x.abs(), v.y.abs(), v.z.abs()])
        .fold(0.0_f64, f64::max);
    let tolerance = 16.0 * f64::EPSILON * scale;
    left.iter().zip(&right).all(|(a, b)| {
        (a.x - b.x).abs() <= tolerance
            && (a.y - b.y).abs() <= tolerance
            && (a.z - b.z).abs() <= tolerance
    })
}

fn owner_keys(properties: &Properties) -> Vec<PropertyKey> {
    properties.iter().map(|(key, _)| key.clone()).collect()
}

/// Owner annotations excluded by one successful append, in property-key order.
/// Entity properties are transferred under [`ModelEditor::append_model`]'s contract.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ModelAppendReport {
    /// Destination topology owner keys cleared because the system changed.
    pub cleared_topology_properties: Vec<PropertyKey>,
    /// Destination realization owner keys cleared because the model changed.
    pub cleared_model_properties: Vec<PropertyKey>,
    /// Source topology owner keys that do not describe the combined system.
    pub omitted_topology_properties: Vec<PropertyKey>,
    /// Source realization owner keys that do not describe the combined model.
    pub omitted_model_properties: Vec<PropertyKey>,
}

/// One import's source-to-draft identities and retention report.
///
/// Each append has a separate context, even when the same source is appended
/// repeatedly. The source topology is retained; IDs supplied to this mapping are
/// interpreted only in that source. Handles remain stable through later edits;
/// a deleted handle will be rejected by the editor. Use [`Self::published`] to
/// resolve surviving entities after publication, including splits and merges.
#[derive(Debug, Clone)]
pub struct ModelAppend {
    mapping: AppendMapping,
    report: ModelAppendReport,
}

impl ModelAppend {
    pub fn source_topology(&self) -> &Topology {
        &self.mapping.source
    }
    pub fn report(&self) -> &ModelAppendReport {
        &self.report
    }
    pub fn atom(&self, source: InstanceAtomId) -> Result<EditAtomId, ModelEditError> {
        self.mapping
            .atoms
            .get(&source)
            .copied()
            .ok_or_else(|| TopologyEditError::InvalidSourceAtom(source).into())
    }
    pub fn bond(&self, source: InstanceBondId) -> Result<EditBondId, ModelEditError> {
        self.mapping
            .bonds
            .get(&source)
            .copied()
            .ok_or_else(|| TopologyEditError::InvalidSourceBond(source).into())
    }
    pub fn chain(&self, source: ChainId) -> Result<EditChainId, ModelEditError> {
        self.mapping
            .chains
            .get(&source)
            .copied()
            .ok_or_else(|| TopologyEditError::InvalidSourceChain(source).into())
    }
    pub fn residue(&self, source: ResidueId) -> Result<EditResidueId, ModelEditError> {
        self.mapping
            .residues
            .get(&source)
            .copied()
            .ok_or_else(|| TopologyEditError::InvalidSourceResidue(source).into())
    }
    pub fn atom_site(&self, source: AtomSiteId) -> Result<EditAtomSiteId, ModelEditError> {
        self.mapping
            .sites
            .get(&source)
            .copied()
            .ok_or_else(|| TopologyEditError::InvalidSourceAtomSite(source).into())
    }
    /// Borrows the mapping into a publication containing this exact import.
    /// Unrelated publications reject, even if they have equal layouts. Publications
    /// of a cloned draft containing this import are valid, including after all
    /// imported entities have been deleted. No topology bindings are transferred.
    pub fn published<'a>(
        &'a self,
        result: &'a ModelEdit,
    ) -> Result<ModelAppendCorrespondence<'a>, ModelEditError> {
        if !result.correspondence().contains_append(&self.mapping.token) {
            return Err(ModelEditError::ForeignAppend);
        }
        Ok(ModelAppendCorrespondence {
            mapping: &self.mapping,
            published: result.correspondence(),
        })
    }
}

/// Borrowed, transaction-specific correspondence for one published model import.
/// Missing/deleted entity IDs return `None`. An original molecular occurrence can
/// map to zero, one or several final occurrences after deletion, merging or splitting.
#[derive(Debug, Clone, Copy)]
pub struct ModelAppendCorrespondence<'a> {
    mapping: &'a AppendMapping,
    published: &'a TopologyEditCorrespondence,
}

impl ModelAppendCorrespondence<'_> {
    pub fn source_topology(&self) -> &Topology {
        &self.mapping.source
    }
    pub fn target_topology(&self) -> &Topology {
        self.published.target_topology()
    }
    pub fn atom(&self, source: InstanceAtomId) -> Option<InstanceAtomId> {
        self.mapping
            .atoms
            .get(&source)
            .and_then(|&id| self.published.atom(id))
    }
    pub fn bond(&self, source: InstanceBondId) -> Option<InstanceBondId> {
        self.mapping
            .bonds
            .get(&source)
            .and_then(|&id| self.published.bond(id))
    }
    pub fn chain(&self, source: ChainId) -> Option<ChainId> {
        self.mapping
            .chains
            .get(&source)
            .and_then(|&id| self.published.chain(id))
    }
    pub fn residue(&self, source: ResidueId) -> Option<ResidueId> {
        self.mapping
            .residues
            .get(&source)
            .and_then(|&id| self.published.residue(id))
    }
    pub fn atom_site(&self, source: AtomSiteId) -> Option<AtomSiteId> {
        self.mapping
            .sites
            .get(&source)
            .and_then(|&id| self.published.atom_site(id))
    }
    /// Returns distinct surviving occurrences in published instance order.
    /// An invalid source instance rejects; a fully deleted instance returns an empty list.
    pub fn instances(
        &self,
        source: MoleculeInstanceId,
    ) -> Result<Vec<MoleculeInstanceId>, ModelEditError> {
        if self.mapping.source.instance(source).is_err() {
            return Err(TopologyEditError::InvalidSourceInstance(source).into());
        }
        Ok(self
            .mapping
            .atoms
            .iter()
            .filter(|(id, _)| id.molecule() == source)
            .filter_map(|(_, &id)| self.published.atom(id).map(|atom| atom.molecule()))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect())
    }
}
