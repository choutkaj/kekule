use std::fmt;
use std::sync::Arc;

use crate::bio::MacroMolecule;
use crate::core::{Atom, AtomId, Bond, ConformerError, ConformerId, Molecule};
use crate::geometry::{PeriodicCell, Point3};
use crate::small::SmallMolecule;
use crate::topology::{
    InstanceAtomId, InstanceAtomSite, InstanceAtomSiteId, InstanceBondId, InstanceChain,
    InstanceChainId, InstanceResidue, InstanceResidueId, InstanceSmcraHierarchy,
    MoleculeDefinitionId, MoleculeInstance, MoleculeInstanceId, MoleculeInstanceMetadata, Topology,
    TopologyAtomIndex, TopologyBuildError, TopologyBuilder, TopologyError, TopologyMapping,
};
use crate::units::{Quantity, UnitError, MODEL_LENGTH_UNIT};

use super::{AtomData, AtomDataError, BondData, PositionError, Positions, TopologyRemapError};

/// One concrete realization of one immutable topology.
#[derive(Debug, Clone)]
pub struct Model {
    pub(super) topology: Arc<Topology>,
    pub(super) positions: Positions,
    pub(super) cell: Option<PeriodicCell>,
    pub(super) atom_data: AtomData,
    pub(super) bond_data: BondData,
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology)
            && self.positions == other.positions
            && self.cell == other.cell
            && self.atom_data == other.atom_data
            && self.bond_data == other.bond_data
    }
}

impl Model {
    /// Creates a non-periodic model with empty atom and bond data.
    pub fn new(topology: Arc<Topology>, positions: Positions) -> Result<Self, ModelError> {
        if !positions.is_compatible(&topology) {
            return Err(ModelError::TopologyMismatch);
        }
        let atom_data = AtomData::new(&topology);
        let bond_data = BondData::new(&topology);
        Ok(Self {
            topology,
            positions,
            cell: None,
            atom_data,
            bond_data,
        })
    }

    /// Creates a model from complete topology-bound state.
    pub fn with_atom_data(
        topology: Arc<Topology>,
        positions: Positions,
        cell: Option<PeriodicCell>,
        atom_data: AtomData,
    ) -> Result<Self, ModelError> {
        let bond_data = BondData::new(&topology);
        Self::with_data(topology, positions, cell, atom_data, bond_data)
    }

    /// Creates a model from complete topology-bound atom and bond state.
    pub fn with_data(
        topology: Arc<Topology>,
        positions: Positions,
        cell: Option<PeriodicCell>,
        atom_data: AtomData,
        bond_data: BondData,
    ) -> Result<Self, ModelError> {
        if !positions.is_compatible(&topology)
            || !atom_data.is_compatible(&topology)
            || !bond_data.is_compatible(&topology)
        {
            return Err(ModelError::TopologyMismatch);
        }
        Ok(Self {
            topology,
            positions,
            cell,
            atom_data,
            bond_data,
        })
    }

    pub fn builder() -> ModelBuilder {
        ModelBuilder::new()
    }

    pub fn from_small_molecule(
        molecule: &SmallMolecule,
        conformer: ConformerId,
    ) -> Result<Self, ModelBuildError> {
        let mut builder = ModelBuilder::new();
        builder.add_small_molecule(molecule, conformer)?;
        builder.build()
    }

    pub fn from_macro_molecule(
        molecule: &MacroMolecule,
        conformer: ConformerId,
    ) -> Result<Self, ModelBuildError> {
        let mut builder = ModelBuilder::new();
        builder.add_macro_molecule(molecule, conformer)?;
        builder.build()
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub fn atom(&self, atom: InstanceAtomId) -> Result<&Atom, TopologyError> {
        self.topology.atom(atom)
    }

    pub fn bond(&self, bond: InstanceBondId) -> Result<&Bond, TopologyError> {
        self.topology.bond(bond)
    }

    pub fn atoms(&self) -> impl ExactSizeIterator<Item = (InstanceAtomId, &Atom)> {
        self.topology.atoms()
    }

    pub fn bonds(&self) -> impl ExactSizeIterator<Item = (InstanceBondId, &Bond)> {
        self.topology.bonds()
    }

    pub fn instances(
        &self,
    ) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &MoleculeInstance)> {
        self.topology.instances()
    }

    pub fn hierarchy(
        &self,
        instance: MoleculeInstanceId,
    ) -> Result<Option<InstanceSmcraHierarchy<'_>>, TopologyError> {
        self.topology.hierarchy(instance)
    }

    pub fn chains(&self) -> impl Iterator<Item = InstanceChain<'_>> {
        self.topology.chains()
    }

    pub fn residues(&self) -> impl Iterator<Item = InstanceResidue<'_>> {
        self.topology.residues()
    }

    pub fn atom_sites(&self) -> impl Iterator<Item = InstanceAtomSite<'_>> {
        self.topology.atom_sites()
    }

    pub fn chain(&self, chain: InstanceChainId) -> Result<InstanceChain<'_>, TopologyError> {
        self.topology.chain(chain)
    }

    pub fn residue(
        &self,
        residue: InstanceResidueId,
    ) -> Result<InstanceResidue<'_>, TopologyError> {
        self.topology.residue(residue)
    }

    pub fn atom_site(
        &self,
        atom_site: InstanceAtomSiteId,
    ) -> Result<InstanceAtomSite<'_>, TopologyError> {
        self.topology.atom_site(atom_site)
    }

    pub fn atom_for_site(
        &self,
        atom_site: InstanceAtomSiteId,
    ) -> Result<InstanceAtomId, TopologyError> {
        self.topology.atom_for_site(atom_site)
    }

    pub fn atom_site_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceAtomSite<'_>>, TopologyError> {
        self.topology.atom_site_for_atom(atom)
    }

    pub fn residue_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceResidue<'_>>, TopologyError> {
        self.topology.residue_for_atom(atom)
    }

    pub fn chain_for_atom(
        &self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceChain<'_>>, TopologyError> {
        self.topology.chain_for_atom(atom)
    }

    pub fn residue_for_site(
        &self,
        atom_site: InstanceAtomSiteId,
    ) -> Result<InstanceResidue<'_>, TopologyError> {
        self.topology.residue_for_site(atom_site)
    }

    pub fn chain_for_residue(
        &self,
        residue: InstanceResidueId,
    ) -> Result<InstanceChain<'_>, TopologyError> {
        self.topology.chain_for_residue(residue)
    }

    pub const fn positions(&self) -> &Positions {
        &self.positions
    }

    pub fn position(&self, atom: InstanceAtomId) -> Result<Quantity<Point3>, PositionError> {
        self.positions.position(&self.topology, atom)
    }

    pub fn position_at(&self, index: TopologyAtomIndex) -> Result<Quantity<Point3>, PositionError> {
        self.positions.position_at(index)
    }

    pub fn set_position(
        &mut self,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), PositionError> {
        self.positions.set_position(&self.topology, atom, position)
    }

    pub fn set_positions<T>(&mut self, positions: Quantity<T>) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        self.positions.set_all(&self.topology, positions)
    }

    pub const fn cell(&self) -> Option<&PeriodicCell> {
        self.cell.as_ref()
    }

    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.cell = cell;
    }

    pub fn atom_data(&self) -> &AtomData {
        &self.atom_data
    }

    pub fn atom_data_mut(&mut self) -> &mut AtomData {
        &mut self.atom_data
    }

    pub fn set_atom_data(&mut self, atom_data: AtomData) -> Result<(), ModelError> {
        if !atom_data.is_compatible(&self.topology) {
            return Err(ModelError::TopologyMismatch);
        }
        self.atom_data = atom_data;
        Ok(())
    }

    pub fn bond_data(&self) -> &BondData {
        &self.bond_data
    }

    pub fn bond_data_mut(&mut self) -> &mut BondData {
        &mut self.bond_data
    }

    pub fn set_bond_data(&mut self, bond_data: BondData) -> Result<(), ModelError> {
        if !bond_data.is_compatible(&self.topology) {
            return Err(ModelError::TopologyMismatch);
        }
        self.bond_data = bond_data;
        Ok(())
    }

    pub fn occupancy(&self, atom: InstanceAtomId) -> Result<Option<f64>, AtomDataError> {
        self.atom_data.occupancy(&self.topology, atom)
    }

    pub fn b_factor(&self, atom: InstanceAtomId) -> Result<Option<Quantity<f64>>, AtomDataError> {
        self.atom_data.b_factor(&self.topology, atom)
    }

    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }

    pub fn view(&self) -> ModelView<'_> {
        ModelView {
            topology: &self.topology,
            positions: &self.positions,
            cell: self.cell.as_ref(),
            atom_data: &self.atom_data,
            bond_data: &self.bond_data,
        }
    }

    /// Remaps this model to an explicitly related target topology.
    ///
    /// The source model is unchanged. Positions, cell, [`AtomData`], and
    /// [`BondData`] are staged before the target model is returned.
    ///
    /// # Examples
    ///
    /// ```
    /// use kekule::core::{Atom, Element, Molecule};
    /// use kekule::geometry::Point3;
    /// use kekule::small::SmallMolecule;
    /// use kekule::structure::{Model, Positions};
    /// use kekule::topology::{
    ///     transform, MoleculeInstanceMetadata, MoleculeRole, TopologyBuilder,
    /// };
    /// use kekule::units::{Quantity, ANGSTROM};
    /// use std::sync::Arc;
    ///
    /// let mut ligand_builder = Molecule::builder();
    /// ligand_builder.add_atom(Atom::new(Element::from_symbol("C").unwrap()))?;
    /// let ligand = SmallMolecule::from_molecule(ligand_builder.build()?);
    /// let mut water_builder = Molecule::builder();
    /// water_builder.add_atom(Atom::new(Element::from_symbol("O").unwrap()))?;
    /// let water = SmallMolecule::from_molecule(water_builder.build()?);
    /// let mut builder = TopologyBuilder::new();
    /// let ligand_definition = builder.add_small_molecule_definition(&ligand)?;
    /// let water_definition = builder.add_small_molecule_definition(&water)?;
    /// builder.add_instance(ligand_definition, MoleculeInstanceMetadata::default())?;
    /// let mut solvent = MoleculeInstanceMetadata::default();
    /// solvent.insert_role(MoleculeRole::Solvent);
    /// let water_instance = builder.add_instance(water_definition, solvent)?;
    /// let topology = Arc::new(builder.build()?);
    /// let positions = Positions::new(
    ///     &topology,
    ///     Quantity::new(
    ///         vec![Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 0.0, 0.0)],
    ///         ANGSTROM,
    ///     ),
    /// )?;
    /// let model = Model::new(Arc::clone(&topology), positions)?;
    ///
    /// let edit = transform::remove_instances(&topology, [water_instance])?;
    /// let target = edit.shared_topology();
    /// let stripped = model.remap_to(&target, edit.mapping())?;
    /// assert_eq!(stripped.atom_count(), 1);
    /// assert_eq!(stripped.positions().values().value()[0].x, 0.0);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn remap_to(
        &self,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TopologyRemapError> {
        let positions = self.positions.remap_to(&self.topology, target, mapping)?;
        let atom_data = self.atom_data.remap_to(&self.topology, target, mapping)?;
        let bond_data = self.bond_data.remap_to(&self.topology, target, mapping)?;
        Ok(Self {
            topology: Arc::clone(target),
            positions,
            cell: self.cell,
            atom_data,
            bond_data,
        })
    }

    /// Copies one instance's current positions to a compatible local conformer.
    pub fn instance_to_conformer(
        &self,
        instance: MoleculeInstanceId,
        target: &mut Molecule,
        conformer: ConformerId,
    ) -> Result<(), InstanceToConformerError> {
        let source = self
            .topology
            .graph_for_instance(instance)
            .map_err(|_| InstanceToConformerError::InvalidMoleculeInstanceId(instance))?;
        let mut updated = target
            .conformer(conformer)
            .map_err(|_| InstanceToConformerError::InvalidConformerId(conformer))?
            .clone();

        for atom in source.atom_ids() {
            if target.atom(atom).is_err() {
                return Err(InstanceToConformerError::MissingTargetAtom(atom));
            }
        }
        for atom in target.atom_ids() {
            if source.atom(atom).is_err() {
                return Err(InstanceToConformerError::UnexpectedTargetAtom(atom));
            }
        }
        for atom in source.atom_ids() {
            let qualified = InstanceAtomId::new(instance, atom);
            updated.set_position(atom, self.position(qualified)?)?;
        }
        *target
            .conformer_mut(conformer)
            .expect("validated conformer remains live") = updated;
        Ok(())
    }
}

/// Borrowed topology, positions, cell, atom data, and bond data for structural
/// kernels.
#[derive(Debug, Clone, Copy)]
pub struct ModelView<'a> {
    topology: &'a Arc<Topology>,
    positions: &'a Positions,
    cell: Option<&'a PeriodicCell>,
    atom_data: &'a AtomData,
    bond_data: &'a BondData,
}

impl<'a> ModelView<'a> {
    pub fn new(
        topology: &'a Arc<Topology>,
        positions: &'a Positions,
        cell: Option<&'a PeriodicCell>,
        atom_data: &'a AtomData,
        bond_data: &'a BondData,
    ) -> Result<Self, ModelError> {
        if !positions.is_compatible(topology)
            || !atom_data.is_compatible(topology)
            || !bond_data.is_compatible(topology)
        {
            return Err(ModelError::TopologyMismatch);
        }
        Ok(Self {
            topology,
            positions,
            cell,
            atom_data,
            bond_data,
        })
    }

    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub(crate) const fn topology_arc(self) -> &'a Arc<Topology> {
        self.topology
    }

    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    pub fn atom(self, atom: InstanceAtomId) -> Result<&'a Atom, TopologyError> {
        self.topology.atom(atom)
    }

    pub fn bond(self, bond: InstanceBondId) -> Result<&'a Bond, TopologyError> {
        self.topology.bond(bond)
    }

    pub fn atoms(self) -> impl ExactSizeIterator<Item = (InstanceAtomId, &'a Atom)> + 'a {
        self.topology.atoms()
    }

    pub fn bonds(self) -> impl ExactSizeIterator<Item = (InstanceBondId, &'a Bond)> + 'a {
        self.topology.bonds()
    }

    pub fn instances(
        self,
    ) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &'a MoleculeInstance)> + 'a {
        self.topology.instances()
    }

    pub fn hierarchy(
        self,
        instance: MoleculeInstanceId,
    ) -> Result<Option<InstanceSmcraHierarchy<'a>>, TopologyError> {
        self.topology.hierarchy(instance)
    }

    pub fn chains(self) -> impl Iterator<Item = InstanceChain<'a>> + 'a {
        self.topology.chains()
    }

    pub fn residues(self) -> impl Iterator<Item = InstanceResidue<'a>> + 'a {
        self.topology.residues()
    }

    pub fn atom_sites(self) -> impl Iterator<Item = InstanceAtomSite<'a>> + 'a {
        self.topology.atom_sites()
    }

    pub fn chain(self, chain: InstanceChainId) -> Result<InstanceChain<'a>, TopologyError> {
        self.topology.chain(chain)
    }

    pub fn residue(self, residue: InstanceResidueId) -> Result<InstanceResidue<'a>, TopologyError> {
        self.topology.residue(residue)
    }

    pub fn atom_site(
        self,
        atom_site: InstanceAtomSiteId,
    ) -> Result<InstanceAtomSite<'a>, TopologyError> {
        self.topology.atom_site(atom_site)
    }

    pub fn atom_for_site(
        self,
        atom_site: InstanceAtomSiteId,
    ) -> Result<InstanceAtomId, TopologyError> {
        self.topology.atom_for_site(atom_site)
    }

    pub fn atom_site_for_atom(
        self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceAtomSite<'a>>, TopologyError> {
        self.topology.atom_site_for_atom(atom)
    }

    pub fn residue_for_atom(
        self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceResidue<'a>>, TopologyError> {
        self.topology.residue_for_atom(atom)
    }

    pub fn chain_for_atom(
        self,
        atom: InstanceAtomId,
    ) -> Result<Option<InstanceChain<'a>>, TopologyError> {
        self.topology.chain_for_atom(atom)
    }

    pub fn residue_for_site(
        self,
        atom_site: InstanceAtomSiteId,
    ) -> Result<InstanceResidue<'a>, TopologyError> {
        self.topology.residue_for_site(atom_site)
    }

    pub fn chain_for_residue(
        self,
        residue: InstanceResidueId,
    ) -> Result<InstanceChain<'a>, TopologyError> {
        self.topology.chain_for_residue(residue)
    }

    pub const fn positions(self) -> &'a Positions {
        self.positions
    }

    pub fn position(self, atom: InstanceAtomId) -> Result<Quantity<Point3>, PositionError> {
        self.positions.position(self.topology, atom)
    }

    pub fn position_at(self, index: TopologyAtomIndex) -> Result<Quantity<Point3>, PositionError> {
        self.positions.position_at(index)
    }

    pub const fn cell(self) -> Option<&'a PeriodicCell> {
        self.cell
    }

    pub const fn atom_data(self) -> &'a AtomData {
        self.atom_data
    }

    pub const fn bond_data(self) -> &'a BondData {
        self.bond_data
    }

    pub fn occupancy(self, atom: InstanceAtomId) -> Result<Option<f64>, AtomDataError> {
        self.atom_data.occupancy(self.topology, atom)
    }

    pub fn b_factor(self, atom: InstanceAtomId) -> Result<Option<Quantity<f64>>, AtomDataError> {
        self.atom_data.b_factor(self.topology, atom)
    }

    pub fn atom_count(self) -> usize {
        self.topology.atom_count()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ModelError {
    TopologyMismatch,
}

impl fmt::Display for ModelError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("model state belongs to a different topology")
            }
        }
    }
}

impl std::error::Error for ModelError {}

/// Convenience builder that assembles topology and one complete model.
///
/// Macro-molecule insertion validates only the explicitly selected conformer
/// while staging positions; topology insertion separately validates static
/// graph/hierarchy consistency and ignores every unselected conformer.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ModelBuilder {
    topology: TopologyBuilder,
    positions: Vec<Point3>,
}

impl ModelBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn topology_builder(&self) -> &TopologyBuilder {
        &self.topology
    }

    pub fn add_small_molecule_definition(
        &mut self,
        molecule: &SmallMolecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self.topology.add_small_molecule_definition(molecule)?)
    }

    pub fn add_macro_molecule_definition(
        &mut self,
        molecule: &MacroMolecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self.topology.add_macro_molecule_definition(molecule)?)
    }

    pub(crate) fn add_small_molecule_definition_unchecked_connectedness(
        &mut self,
        molecule: &SmallMolecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self
            .topology
            .add_small_molecule_definition_unchecked_connectedness(molecule)?)
    }

    pub(crate) fn add_macro_molecule_definition_unchecked_connectedness(
        &mut self,
        molecule: &MacroMolecule,
    ) -> Result<MoleculeDefinitionId, ModelBuildError> {
        Ok(self
            .topology
            .add_macro_molecule_definition_unchecked_connectedness(molecule)?)
    }

    pub fn add_instance<T>(
        &mut self,
        definition: MoleculeDefinitionId,
        positions: Quantity<T>,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, ModelBuildError>
    where
        T: AsRef<[Point3]>,
    {
        let graph = self.topology.definition(definition)?.graph();
        let staged = stage_instance_positions(graph, positions)?;
        self.positions
            .try_reserve(staged.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let instance = self.topology.add_instance(definition, metadata)?;
        self.positions.extend(staged);
        Ok(instance)
    }

    pub fn add_small_molecule(
        &mut self,
        molecule: &SmallMolecule,
        conformer: ConformerId,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        self.add_small_molecule_with_metadata(
            molecule,
            conformer,
            MoleculeInstanceMetadata::default(),
        )
    }

    pub fn add_small_molecule_with_metadata(
        &mut self,
        molecule: &SmallMolecule,
        conformer: ConformerId,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        if molecule.as_molecule().atom_count() == 0 {
            return Err(ModelBuildError::Topology(
                TopologyBuildError::EmptyMoleculeDefinition,
            ));
        }
        let staged = stage_conformer_positions(molecule.as_molecule(), conformer)?;
        self.positions
            .try_reserve(staged.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let (_, instance) = self
            .topology
            .add_small_molecule_instance(molecule, metadata)?;
        self.positions.extend(staged);
        Ok(instance)
    }

    pub(crate) fn add_small_molecule_with_metadata_unchecked_connectedness(
        &mut self,
        molecule: &SmallMolecule,
        conformer: ConformerId,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        if molecule.as_molecule().atom_count() == 0 {
            return Err(ModelBuildError::Topology(
                TopologyBuildError::EmptyMoleculeDefinition,
            ));
        }
        let staged = stage_conformer_positions(molecule.as_molecule(), conformer)?;
        self.positions
            .try_reserve(staged.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let (_, instance) = self
            .topology
            .add_small_molecule_instance_unchecked_connectedness(molecule, metadata)?;
        self.positions.extend(staged);
        Ok(instance)
    }

    pub fn add_macro_molecule(
        &mut self,
        molecule: &MacroMolecule,
        conformer: ConformerId,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        self.add_macro_molecule_with_metadata(
            molecule,
            conformer,
            MoleculeInstanceMetadata::default(),
        )
    }

    pub fn add_macro_molecule_with_metadata(
        &mut self,
        molecule: &MacroMolecule,
        conformer: ConformerId,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        if molecule.as_molecule().atom_count() == 0 {
            return Err(ModelBuildError::Topology(
                TopologyBuildError::EmptyMoleculeDefinition,
            ));
        }
        let staged = stage_conformer_positions(molecule.as_molecule(), conformer)?;
        self.positions
            .try_reserve(staged.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let (_, instance) = self
            .topology
            .add_macro_molecule_instance(molecule, metadata)?;
        self.positions.extend(staged);
        Ok(instance)
    }

    pub(crate) fn add_macro_molecule_with_metadata_unchecked_connectedness(
        &mut self,
        molecule: &MacroMolecule,
        conformer: ConformerId,
        metadata: MoleculeInstanceMetadata,
    ) -> Result<MoleculeInstanceId, ModelBuildError> {
        if molecule.as_molecule().atom_count() == 0 {
            return Err(ModelBuildError::Topology(
                TopologyBuildError::EmptyMoleculeDefinition,
            ));
        }
        let staged = stage_conformer_positions(molecule.as_molecule(), conformer)?;
        self.positions
            .try_reserve(staged.len())
            .map_err(|_| ModelBuildError::CapacityOverflow)?;
        let (_, instance) = self
            .topology
            .add_macro_molecule_instance_unchecked_connectedness(molecule, metadata)?;
        self.positions.extend(staged);
        Ok(instance)
    }

    pub fn build(self) -> Result<Model, ModelBuildError> {
        let topology = Arc::new(self.topology.build()?);
        let positions =
            Positions::new(&topology, Quantity::new(self.positions, MODEL_LENGTH_UNIT))?;
        Ok(Model::new(topology, positions).expect("builder creates exactly topology-bound state"))
    }
}

pub(super) fn stage_conformer_positions(
    graph: &Molecule,
    conformer: ConformerId,
) -> Result<Vec<Point3>, ModelBuildError> {
    let conformer = graph
        .conformer(conformer)
        .map_err(|_| ModelBuildError::InvalidConformerId(conformer))?;
    let mut positions = Vec::new();
    positions
        .try_reserve_exact(graph.atom_count())
        .map_err(|_| ModelBuildError::CapacityOverflow)?;
    for atom in graph.atom_ids() {
        let point = conformer
            .position(atom)
            .ok_or(ModelBuildError::MissingPosition { atom })?
            .to_unit(MODEL_LENGTH_UNIT)?
            .to_value();
        if !point.is_finite() {
            return Err(ModelBuildError::NonFinitePosition { atom });
        }
        positions.push(point);
    }
    Ok(positions)
}

fn stage_instance_positions<T>(
    graph: &Molecule,
    positions: Quantity<T>,
) -> Result<Vec<Point3>, ModelBuildError>
where
    T: AsRef<[Point3]>,
{
    let source = positions.value().as_ref();
    if source.len() != graph.atom_count() {
        return Err(ModelBuildError::InstancePositionCountMismatch {
            expected: graph.atom_count(),
            actual: source.len(),
        });
    }
    let factor = positions.unit().conversion_factor_to(MODEL_LENGTH_UNIT)?;
    source
        .iter()
        .copied()
        .zip(graph.atom_ids())
        .map(|(point, atom)| {
            let point = Point3::new(point.x * factor, point.y * factor, point.z * factor);
            if !point.is_finite() {
                return Err(ModelBuildError::NonFinitePosition { atom });
            }
            Ok(point)
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelBuildError {
    InvalidConformerId(ConformerId),
    MissingPosition { atom: AtomId },
    NonFinitePosition { atom: AtomId },
    InstancePositionCountMismatch { expected: usize, actual: usize },
    CapacityOverflow,
    Topology(TopologyBuildError),
    Position(PositionError),
    Unit(UnitError),
}

impl fmt::Display for ModelBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConformerId(id) => write!(formatter, "invalid source conformer: {id}"),
            Self::MissingPosition { atom } => {
                write!(
                    formatter,
                    "source conformer has no position for atom {atom}"
                )
            }
            Self::NonFinitePosition { atom } => {
                write!(formatter, "source position for atom {atom} is not finite")
            }
            Self::InstancePositionCountMismatch { expected, actual } => write!(
                formatter,
                "definition instance requires {expected} positions, but received {actual}"
            ),
            Self::CapacityOverflow => {
                formatter.write_str("model construction exceeds coordinate capacity")
            }
            Self::Topology(error) => write!(formatter, "cannot build topology: {error}"),
            Self::Position(error) => write!(formatter, "cannot build positions: {error}"),
            Self::Unit(error) => write!(formatter, "invalid source position unit: {error}"),
        }
    }
}

impl std::error::Error for ModelBuildError {}

impl From<TopologyBuildError> for ModelBuildError {
    fn from(error: TopologyBuildError) -> Self {
        Self::Topology(error)
    }
}

impl From<PositionError> for ModelBuildError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<UnitError> for ModelBuildError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum InstanceToConformerError {
    InvalidMoleculeInstanceId(MoleculeInstanceId),
    InvalidConformerId(ConformerId),
    MissingTargetAtom(AtomId),
    UnexpectedTargetAtom(AtomId),
    Position(PositionError),
    Conformer(ConformerError),
    Unit(UnitError),
}

impl fmt::Display for InstanceToConformerError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMoleculeInstanceId(id) => {
                write!(formatter, "invalid molecule instance: {id}")
            }
            Self::InvalidConformerId(id) => write!(formatter, "invalid target conformer: {id}"),
            Self::MissingTargetAtom(atom) => {
                write!(formatter, "target molecule is missing source atom {atom}")
            }
            Self::UnexpectedTargetAtom(atom) => {
                write!(formatter, "target molecule contains unexpected atom {atom}")
            }
            Self::Position(error) => write!(formatter, "cannot read model position: {error}"),
            Self::Conformer(error) => write!(formatter, "cannot update target conformer: {error}"),
            Self::Unit(error) => write!(formatter, "invalid target conformer unit: {error}"),
        }
    }
}

impl std::error::Error for InstanceToConformerError {}

impl From<PositionError> for InstanceToConformerError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<ConformerError> for InstanceToConformerError {
    fn from(error: ConformerError) -> Self {
        Self::Conformer(error)
    }
}

impl From<UnitError> for InstanceToConformerError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
