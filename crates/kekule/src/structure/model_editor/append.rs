use super::*;
use crate::properties::PropertyKey;
use crate::structure::AsModelView;
use crate::topology::{AppendMapping, Topology};
use crate::topology::{AtomSiteId, ChainId, InstanceAtomId, InstanceBondId, ResidueId};

impl ModelEditor {
    /// Appends a complete borrowed model, including its existing coordinates.
    ///
    /// Accepts any [`AsModelView`] source, such as a model, an ensemble member,
    /// or a trajectory frame, without first creating an owned model. Positions are used as supplied:
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
    /// let ligand = Model::from_molecule(molecule, &Positions::new(Quantity::new(
    ///     vec![Point3::origin(), Point3::new(1.4, 0.0, 0.0)], ANGSTROM))?)?;
    /// let mut editor = ligand.to_editor();
    /// let appended = editor.append_model(&ligand)?;
    /// let source_atom = ligand.topology().atom_ids()[0];
    /// let draft_atom = appended.atom(source_atom).expect("appended atom");
    /// assert_eq!(editor.position(draft_atom), ligand.position(source_atom));
    /// let model = editor.finish()?;
    /// assert_eq!(model.atom_count(), 4);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn append_model(
        &mut self,
        source: &(impl AsModelView + ?Sized),
    ) -> Result<ModelAppend, ModelEditError> {
        let source = source.as_model_view();
        let cell = match (self.cell().copied(), source.cell().copied()) {
            (None, cell) if self.topology.is_empty() => cell,
            (cell, None) => cell,
            (Some(left), Some(right)) if compatible_cells(left, right) => Some(left),
            _ => return Err(ModelEditError::IncompatibleAppendCell),
        };
        let report = ModelAppendReport {
            cleared_topology_properties: self
                .topology
                .properties()
                .owner()
                .keys()
                .cloned()
                .collect(),
            cleared_model_properties: self.slots.properties.owner().keys().cloned().collect(),
            omitted_topology_properties: source
                .topology()
                .properties()
                .owner()
                .keys()
                .cloned()
                .collect(),
            omitted_model_properties: source.properties().owner().keys().cloned().collect(),
        };
        let mut staged = self.clone();
        staged.slots.try_reserve_atoms(source.atom_count())?;
        let mapping = staged.topology.append_topology(source.shared_topology())?;
        staged.slots.resize_slots(
            staged.topology.atom_slot_count(),
            staged.topology.bond_slot_count(),
        );
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
        staged.slots.copy_atom_state_from(source, &atom_rows);
        let properties = &mut staged.slots.properties;
        properties
            .atoms_raw_mut()
            .copy_rows_from(source.properties().atoms().raw(), &atom_rows)
            .map_err(|error| ModelEditError::AppendProperty {
                domain: crate::properties::PropertyDomain::ModelAtoms,
                error: Box::new(error),
            })?;
        properties
            .bonds_raw_mut()
            .copy_rows_from(source.properties().bonds().raw(), &bond_rows)
            .map_err(|error| ModelEditError::AppendProperty {
                domain: crate::properties::PropertyDomain::ModelBonds,
                error: Box::new(error),
            })?;
        properties.owner_mut().clear();
        staged.slots.cell = cell;
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

/// Owner annotations excluded by one successful append, in property-key order.
/// Entity properties are transferred under [`ModelEditor::append_model`]'s contract.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelAppendReport {
    pub(crate) cleared_topology_properties: Vec<PropertyKey>,
    pub(crate) cleared_model_properties: Vec<PropertyKey>,
    pub(crate) omitted_topology_properties: Vec<PropertyKey>,
    pub(crate) omitted_model_properties: Vec<PropertyKey>,
}

impl ModelAppendReport {
    /// Destination topology owner keys cleared because the system changed.
    pub fn cleared_topology_properties(&self) -> &[PropertyKey] {
        &self.cleared_topology_properties
    }

    /// Destination structure owner keys cleared because the model changed.
    pub fn cleared_model_properties(&self) -> &[PropertyKey] {
        &self.cleared_model_properties
    }

    /// Source topology owner keys that do not describe the combined system.
    pub fn omitted_topology_properties(&self) -> &[PropertyKey] {
        &self.omitted_topology_properties
    }

    /// Source structure owner keys that do not describe the combined model.
    pub fn omitted_model_properties(&self) -> &[PropertyKey] {
        &self.omitted_model_properties
    }
}

/// One import's source-to-draft identities and retention report.
///
/// Each append has a separate context, even when the same source is appended
/// repeatedly. The source topology is retained; IDs supplied to this mapping are
/// interpreted only in that source. Handles remain stable through later edits;
/// a deleted handle will be rejected by the editor. These handles are only for
/// use within the draft, before [`ModelEditor::finish`].
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
    /// The draft handle of an appended source atom, or `None` when the source
    /// has no such atom.
    pub fn atom(&self, source: InstanceAtomId) -> Option<EditAtomId> {
        self.mapping.atoms.get(&source).copied()
    }
    /// The draft handle of an appended source bond, or `None` when the source
    /// has no such bond.
    pub fn bond(&self, source: InstanceBondId) -> Option<EditBondId> {
        self.mapping.bonds.get(&source).copied()
    }
    /// The draft handle of an appended source chain, or `None` when the source
    /// has no such chain.
    pub fn chain(&self, source: ChainId) -> Option<EditChainId> {
        self.mapping.chains.get(&source).copied()
    }
    /// The draft handle of an appended source residue, or `None` when the
    /// source has no such residue.
    pub fn residue(&self, source: ResidueId) -> Option<EditResidueId> {
        self.mapping.residues.get(&source).copied()
    }
    /// The draft handle of an appended source atom site, or `None` when the
    /// source has no such site.
    pub fn atom_site(&self, source: AtomSiteId) -> Option<EditAtomSiteId> {
        self.mapping.sites.get(&source).copied()
    }
}
