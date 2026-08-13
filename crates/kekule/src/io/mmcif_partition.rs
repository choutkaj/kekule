use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use crate::bio::{MacroMolecule, SmcraHierarchy};
use crate::core::{AtomId, Molecule};
use crate::small::SmallMolecule;
use crate::structure::{
    Configuration, Ensemble, EnsembleMember, Model, Positions, StructureObservation,
};
use crate::topology::{
    InstanceAtomId, MoleculeInstanceId, MoleculeRole, Topology, TopologyBuilder,
};
use crate::units::{Quantity, MODEL_LENGTH_UNIT};

use super::mmcif_connectivity as connectivity;
use super::mmcif_interpret as raw;
use super::MmcifDocument;

/// Final public mmCIF interpretation after authoritative connectivity has been
/// completed and every residual disconnected graph has been partitioned into
/// connected molecule instances.
#[derive(Debug, Clone)]
pub struct MmcifInterpretation {
    model: Model,
    report: raw::MmcifInterpretationReport,
}

impl MmcifInterpretation {
    pub fn model(&self) -> &Model {
        &self.model
    }

    pub fn report(&self) -> &raw::MmcifInterpretationReport {
        &self.report
    }

    pub fn topology(&self) -> &Topology {
        self.model.topology()
    }

    pub fn configuration(&self) -> &Configuration {
        self.model.configuration()
    }

    pub fn observation(&self) -> Option<&StructureObservation> {
        self.model.observation()
    }

    pub fn into_model(self) -> Model {
        self.model
    }

    pub fn into_parts(self) -> (Model, raw::MmcifInterpretationReport) {
        (self.model, self.report)
    }
}

/// Interprets one mmCIF coordinate model into connected molecule instances.
///
/// Authoritative component/polymer/branch/`_struct_conn` connectivity is
/// materialized first. Any remaining graph components are then represented as
/// separate molecule instances; no covalent bond is invented across an
/// unresolved experimental gap.
pub fn interpret_mmcif(
    document: &MmcifDocument,
    options: raw::MmcifInterpretOptions,
) -> Result<MmcifInterpretation, raw::MmcifInterpretError> {
    let interpretation = connectivity::interpret_mmcif(document, options)?;
    let (source, report) = interpretation.into_parts();
    let source_topology = source.shared_topology();
    let partition = partition_topology(&source_topology)?;
    let configuration = remap_configuration(
        source.configuration(),
        &source_topology,
        &partition.topology,
        &partition.source_atoms,
    )?;
    let observation = source
        .observation()
        .map(|observation| {
            remap_observation(
                observation,
                &source_topology,
                &partition.topology,
                &partition.source_atoms,
            )
        })
        .transpose()?;
    let model =
        Model::with_observation(Arc::clone(&partition.topology), configuration, observation)
            .map_err(interpret_error)?;
    let report = remap_report(report, &partition)?;
    Ok(MmcifInterpretation { model, report })
}

#[derive(Debug, Clone)]
pub struct MmcifEnsembleInterpretation {
    ensemble: Ensemble,
    reports: Vec<raw::MmcifInterpretationReport>,
}

impl MmcifEnsembleInterpretation {
    pub fn ensemble(&self) -> &Ensemble {
        &self.ensemble
    }

    pub fn reports(&self) -> &[raw::MmcifInterpretationReport] {
        &self.reports
    }

    pub fn into_parts(self) -> (Ensemble, Vec<raw::MmcifInterpretationReport>) {
        (self.ensemble, self.reports)
    }
}

/// Interprets an mmCIF ensemble using one shared connected-fragment topology.
pub fn interpret_mmcif_ensemble(
    document: &MmcifDocument,
    options: raw::MmcifEnsembleInterpretOptions,
) -> Result<MmcifEnsembleInterpretation, raw::MmcifEnsembleInterpretError> {
    let interpretation = connectivity::interpret_mmcif_ensemble(document, options)?;
    let (source, reports) = interpretation.into_parts();
    let source_topology = source.shared_topology();
    let partition = partition_topology(&source_topology).map_err(|error| {
        raw::MmcifEnsembleInterpretError::Model {
            model_id: reports
                .first()
                .and_then(raw::MmcifInterpretationReport::selected_model)
                .unwrap_or("<unknown>")
                .to_owned(),
            error,
        }
    })?;

    let topology = Arc::clone(&partition.topology);
    let mut ensemble = Ensemble::new(Arc::clone(&topology));
    for member in source.members() {
        let configuration = remap_configuration(
            member.configuration(),
            &source_topology,
            &topology,
            &partition.source_atoms,
        )
        .map_err(|error| raw::MmcifEnsembleInterpretError::Model {
            model_id: member
                .observation()
                .and_then(StructureObservation::source_model_id)
                .unwrap_or("<unknown>")
                .to_owned(),
            error,
        })?;
        let mut rebuilt = EnsembleMember::new(configuration);
        rebuilt
            .set_weight(member.weight())
            .map_err(|error| raw::MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
        if let Some(observation) = member.observation() {
            let observation = remap_observation(
                observation,
                &source_topology,
                &topology,
                &partition.source_atoms,
            )
            .map_err(|error| raw::MmcifEnsembleInterpretError::Model {
                model_id: observation
                    .source_model_id()
                    .unwrap_or("<unknown>")
                    .to_owned(),
                error,
            })?;
            rebuilt
                .set_observation(Some(observation))
                .map_err(|error| raw::MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
        }
        rebuilt.props_mut().clone_from(member.props());
        ensemble
            .push(rebuilt)
            .map_err(|error| raw::MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
    }

    let reports = reports
        .into_iter()
        .map(|report| {
            let model_id = report.selected_model().unwrap_or("<unknown>").to_owned();
            remap_report(report, &partition)
                .map_err(|error| raw::MmcifEnsembleInterpretError::Model { model_id, error })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(MmcifEnsembleInterpretation { ensemble, reports })
}

struct PartitionedTopology {
    topology: Arc<Topology>,
    /// Source atom corresponding to each dense atom of `topology`.
    source_atoms: Vec<InstanceAtomId>,
    /// Complete source-to-partitioned semantic atom mapping.
    atom_map: BTreeMap<InstanceAtomId, InstanceAtomId>,
}

fn partition_topology(source: &Topology) -> Result<PartitionedTopology, raw::MmcifInterpretError> {
    let mut builder = TopologyBuilder::new();
    let mut source_atoms = Vec::with_capacity(source.atom_count());
    let mut atom_map = BTreeMap::new();

    for (source_instance, source_metadata) in source.instances() {
        let definition = source
            .definition_for_instance(source_instance)
            .map_err(interpret_error)?;
        let components = definition.graph().connected_components();
        if components.is_empty() {
            return Err(interpret_error("mmCIF molecule instance has no atoms"));
        }

        if components.len() == 1 {
            let definition_id = if let Some(molecule) = definition.macro_molecule() {
                builder
                    .add_macro_molecule_definition(molecule)
                    .map_err(interpret_error)?
            } else if let Some(molecule) = definition.small_molecule() {
                builder
                    .add_small_molecule_definition(molecule)
                    .map_err(interpret_error)?
            } else {
                return Err(interpret_error("mmCIF definition has no molecule payload"));
            };
            let target_instance = builder
                .add_instance(definition_id, source_metadata.metadata().clone())
                .map_err(interpret_error)?;
            for atom in definition.graph().atom_ids() {
                let source_atom = InstanceAtomId::new(source_instance, atom);
                let target_atom = InstanceAtomId::new(target_instance, atom);
                source_atoms.push(source_atom);
                if atom_map.insert(source_atom, target_atom).is_some() {
                    return Err(interpret_error(
                        "duplicate mmCIF source atom during partition",
                    ));
                }
            }
            continue;
        }

        for component in components {
            let ExtractedConnectedGraph {
                graph,
                atom_map: local_map,
                ordered_source_atoms,
            } = extract_connected_graph(definition.graph(), &component)?;
            let definition_id = if let Some(molecule) = definition.macro_molecule() {
                let hierarchy = extract_hierarchy(molecule.hierarchy(), &local_map)?;
                let molecule =
                    MacroMolecule::try_from_parts(graph, hierarchy).map_err(interpret_error)?;
                builder
                    .add_macro_molecule_definition(&molecule)
                    .map_err(interpret_error)?
            } else if definition.small_molecule().is_some() {
                let molecule = SmallMolecule::from_graph(graph);
                builder
                    .add_small_molecule_definition(&molecule)
                    .map_err(interpret_error)?
            } else {
                return Err(interpret_error("mmCIF definition has no molecule payload"));
            };
            let target_instance = builder
                .add_instance(definition_id, source_metadata.metadata().clone())
                .map_err(interpret_error)?;
            for source_local in ordered_source_atoms {
                let target_local = local_map[&source_local];
                let source_atom = InstanceAtomId::new(source_instance, source_local);
                let target_atom = InstanceAtomId::new(target_instance, target_local);
                source_atoms.push(source_atom);
                if atom_map.insert(source_atom, target_atom).is_some() {
                    return Err(interpret_error(
                        "duplicate mmCIF source atom during partition",
                    ));
                }
            }
        }
    }

    let topology = Arc::new(builder.build().map_err(interpret_error)?);
    if topology.atom_count() != source_atoms.len() || atom_map.len() != source.atom_count() {
        return Err(interpret_error(
            "mmCIF fragment partition changed the represented atom count",
        ));
    }
    for (index, target_atom) in topology.atom_ids().iter().copied().enumerate() {
        let source_atom = source_atoms[index];
        if atom_map.get(&source_atom).copied() != Some(target_atom) {
            return Err(interpret_error(
                "mmCIF fragment partition produced inconsistent dense atom order",
            ));
        }
    }
    if topology.instances().any(|(instance, _)| {
        topology
            .connected_components(instance)
            .is_ok_and(|parts| parts.len() != 1)
    }) {
        return Err(interpret_error(
            "mmCIF fragment partition produced a disconnected molecule instance",
        ));
    }

    Ok(PartitionedTopology {
        topology,
        source_atoms,
        atom_map,
    })
}

struct ExtractedConnectedGraph {
    graph: Molecule,
    atom_map: BTreeMap<AtomId, AtomId>,
    ordered_source_atoms: Vec<AtomId>,
}

fn extract_connected_graph(
    source: &Molecule,
    component: &[AtomId],
) -> Result<ExtractedConnectedGraph, raw::MmcifInterpretError> {
    let selected = component.iter().copied().collect::<BTreeSet<_>>();
    let mut ordered = component.to_vec();
    ordered.sort_unstable();
    let mut builder = Molecule::builder();
    let mut atom_map = BTreeMap::new();
    let mut bond_props = Vec::new();

    for source_atom in ordered.iter().copied() {
        let atom = source.atom(source_atom).map_err(interpret_error)?.clone();
        let target_atom = builder.add_atom(atom).map_err(interpret_error)?;
        atom_map.insert(source_atom, target_atom);
    }
    for (_, bond) in source.bonds() {
        let (source_left, source_right) = bond.endpoints();
        if !selected.contains(&source_left) || !selected.contains(&source_right) {
            continue;
        }
        let target_left = atom_map[&source_left];
        let target_right = atom_map[&source_right];
        let target_bond = builder
            .add_bond(target_left, target_right, bond.order)
            .map_err(interpret_error)?;
        bond_props.push((target_bond, bond.props.clone()));
    }
    let mut graph = builder.build().map_err(interpret_error)?;
    for (target_bond, props) in bond_props {
        graph.bond_mut(target_bond).map_err(interpret_error)?.props = props;
    }
    graph.props_mut().clone_from(source.props());
    Ok(ExtractedConnectedGraph {
        graph,
        atom_map,
        ordered_source_atoms: ordered,
    })
}

fn extract_hierarchy(
    source: &SmcraHierarchy,
    atom_map: &BTreeMap<AtomId, AtomId>,
) -> Result<SmcraHierarchy, raw::MmcifInterpretError> {
    let mut hierarchy = SmcraHierarchy::new();
    hierarchy.props_mut().clone_from(source.props());

    for (_, source_chain) in source.chains() {
        let selected_residues = source_chain
            .residues()
            .iter()
            .filter_map(|residue_id| source.residue(*residue_id).ok())
            .filter(|residue| {
                residue.atom_sites().iter().any(|site_id| {
                    source
                        .atom_site(*site_id)
                        .ok()
                        .is_some_and(|site| atom_map.contains_key(&site.atom()))
                })
            })
            .collect::<Vec<_>>();
        if selected_residues.is_empty() {
            continue;
        }

        let chain = hierarchy
            .add_chain(
                source_chain.label_id().to_owned(),
                source_chain.author_id().map(str::to_owned),
            )
            .map_err(interpret_error)?;
        hierarchy
            .chain_props_mut(chain)
            .map_err(interpret_error)?
            .clone_from(source_chain.props());

        for source_residue in selected_residues {
            let residue = hierarchy
                .add_residue(
                    chain,
                    source_residue.name().to_owned(),
                    source_residue.label_seq_id(),
                    source_residue.author_seq_id().map(str::to_owned),
                    source_residue.insertion_code().map(str::to_owned),
                )
                .map_err(interpret_error)?;
            hierarchy
                .set_residue_component_ids(
                    residue,
                    source_residue.label_comp_id().map(str::to_owned),
                    source_residue.author_comp_id().map(str::to_owned),
                )
                .map_err(interpret_error)?;
            hierarchy
                .residue_props_mut(residue)
                .map_err(interpret_error)?
                .clone_from(source_residue.props());

            for source_site_id in source_residue.atom_sites() {
                let source_site = source.atom_site(*source_site_id).map_err(interpret_error)?;
                let Some(&target_atom) = atom_map.get(&source_site.atom()) else {
                    continue;
                };
                let site = hierarchy
                    .add_atom_site(residue, target_atom, source_site.metadata().clone())
                    .map_err(interpret_error)?;
                hierarchy
                    .atom_site_props_mut(site)
                    .map_err(interpret_error)?
                    .clone_from(source_site.props());
            }
        }
    }
    Ok(hierarchy)
}

fn remap_report(
    mut report: raw::MmcifInterpretationReport,
    partition: &PartitionedTopology,
) -> Result<raw::MmcifInterpretationReport, raw::MmcifInterpretError> {
    let source_instances = std::mem::take(&mut report.instances);
    let mut target_instances = Vec::new();

    for source in source_instances {
        let mut groups = BTreeMap::<MoleculeInstanceId, Vec<raw::MmcifAtomProvenance>>::new();
        for mut atom in source.atoms.clone() {
            let target = partition.atom_map.get(&atom.atom).copied().ok_or_else(|| {
                interpret_error("mmCIF provenance atom was lost during partition")
            })?;
            atom.atom = target;
            groups.entry(target.molecule()).or_default().push(atom);
        }
        for (molecule, atoms) in groups {
            target_instances.push(raw::MmcifInstanceProvenance {
                molecule,
                coordinate_model_id: source.coordinate_model_id.clone(),
                asym_ids: source.asym_ids.clone(),
                entity_ids: source.entity_ids.clone(),
                entity_kinds: source.entity_kinds.clone(),
                atoms,
            });
        }
    }
    target_instances.sort_by_key(raw::MmcifInstanceProvenance::molecule);
    report.instances = target_instances;
    report.template_bonds_pending = 0;
    report.macromolecules = partition
        .topology
        .instances()
        .filter(|(id, _)| {
            partition
                .topology
                .definition_for_instance(*id)
                .is_ok_and(|definition| definition.macro_molecule().is_some())
        })
        .count();
    report.small_molecules = partition
        .topology
        .instances()
        .filter(|(id, _)| {
            partition
                .topology
                .definition_for_instance(*id)
                .is_ok_and(|definition| definition.small_molecule().is_some())
        })
        .count();
    report.solvent_molecules = partition
        .topology
        .instances()
        .filter(|(_, molecule)| molecule.has_role(MoleculeRole::Solvent))
        .count();
    Ok(report)
}

fn remap_configuration(
    source: &Configuration,
    source_topology: &Arc<Topology>,
    target_topology: &Arc<Topology>,
    source_atoms: &[InstanceAtomId],
) -> Result<Configuration, raw::MmcifInterpretError> {
    let positions = source_atoms
        .iter()
        .copied()
        .map(|atom| {
            source
                .positions()
                .position(source_topology, atom)
                .map(|position| position.into_value())
                .map_err(interpret_error)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let positions = Positions::new(target_topology, Quantity::new(positions, MODEL_LENGTH_UNIT))
        .map_err(interpret_error)?;
    Ok(match source.cell().copied() {
        Some(cell) => Configuration::with_cell(positions, cell),
        None => Configuration::new(positions),
    })
}

fn remap_observation(
    source: &StructureObservation,
    source_topology: &Arc<Topology>,
    target_topology: &Arc<Topology>,
    source_atoms: &[InstanceAtomId],
) -> Result<StructureObservation, raw::MmcifInterpretError> {
    let atoms = source_atoms
        .iter()
        .copied()
        .map(|atom| {
            source
                .atom(source_topology, atom)
                .cloned()
                .map_err(interpret_error)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut observation =
        StructureObservation::new(target_topology, atoms).map_err(interpret_error)?;
    observation.set_source_model_id(source.source_model_id().map(str::to_owned));
    observation.props_mut().clone_from(source.props());
    Ok(observation)
}

fn interpret_error(error: impl std::fmt::Display) -> raw::MmcifInterpretError {
    raw::MmcifInterpretError {
        line: None,
        message: error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GAPPED_PEPTIDE: &str = r#"
data_gapped
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.nstd_linkage
1 'polypeptide(L)' no
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
A 1 1 GLY
A 1 2 GLY
A 1 3 GLY
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
GLY N CA sing
GLY CA C sing
GLY C O doub
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY A 1 1 0.0 0.0 0.0
ATOM 2 C CA GLY A 1 1 1.4 0.0 0.0
ATOM 3 C C GLY A 1 1 2.8 0.0 0.0
ATOM 4 O O GLY A 1 1 3.8 0.0 0.0
ATOM 5 N N GLY A 1 3 8.0 0.0 0.0
ATOM 6 C CA GLY A 1 3 9.4 0.0 0.0
ATOM 7 C C GLY A 1 3 10.8 0.0 0.0
ATOM 8 O O GLY A 1 3 11.8 0.0 0.0
"#;

    #[test]
    fn gapped_polymer_is_partitioned_into_connected_macro_instances() {
        let document = super::super::mmcif_document::parse_mmcif_str(
            GAPPED_PEPTIDE,
            super::super::mmcif_document::MmcifParseOptions::default(),
        )
        .expect("parse gapped peptide");
        let interpretation = interpret_mmcif(&document, raw::MmcifInterpretOptions::default())
            .expect("interpret gapped peptide");
        assert_eq!(interpretation.topology().instances().count(), 2);
        assert!(interpretation.topology().instances().all(|(instance, _)| {
            interpretation
                .topology()
                .connected_components(instance)
                .is_ok_and(|components| components.len() == 1)
        }));
        assert_eq!(interpretation.report().instances().len(), 2);
        assert!(interpretation
            .report()
            .instances()
            .iter()
            .all(|instance| { instance.asym_ids() == ["A"] && instance.entity_ids() == ["1"] }));
        assert_eq!(interpretation.report().template_bonds_pending(), 0);
    }
}
