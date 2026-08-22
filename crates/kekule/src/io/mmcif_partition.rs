use std::collections::BTreeMap;
use std::sync::Arc;

use crate::structure::{AtomData, Ensemble, EnsembleMember, Model, Positions};
use crate::topology::{InstanceAtomId, MoleculeInstanceId, Topology, TopologyBuilder};
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

    pub fn to_model(self) -> Model {
        self.model
    }

    pub fn to_parts(self) -> (Model, raw::MmcifInterpretationReport) {
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
    let (source, report) = interpretation.to_parts();
    let source_topology = source.shared_topology();
    let partition = partition_topology(&source_topology)?;
    let positions = remap_positions(
        source.positions(),
        &source_topology,
        &partition.topology,
        &partition.source_atoms,
    )?;
    let atom_data = remap_atom_data(
        source.atom_data(),
        &source_topology,
        &partition.topology,
        &partition.source_atoms,
    )?;
    let model = Model::with_atom_data(
        Arc::clone(&partition.topology),
        positions,
        source.cell().copied(),
        atom_data,
    )
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

    pub fn to_parts(self) -> (Ensemble, Vec<raw::MmcifInterpretationReport>) {
        (self.ensemble, self.reports)
    }
}

/// Interprets an mmCIF ensemble using one shared connected-fragment topology.
pub fn interpret_mmcif_ensemble(
    document: &MmcifDocument,
    options: raw::MmcifEnsembleInterpretOptions,
) -> Result<MmcifEnsembleInterpretation, raw::MmcifEnsembleInterpretError> {
    let interpretation = connectivity::interpret_mmcif_ensemble(document, options)?;
    let (source, reports) = interpretation.to_parts();
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
    for (member, report) in source.members().zip(&reports) {
        let model_id = report.selected_model().unwrap_or("<unknown>").to_owned();
        let positions = remap_positions(
            member.positions(),
            &source_topology,
            &topology,
            &partition.source_atoms,
        )
        .map_err(|error| raw::MmcifEnsembleInterpretError::Model {
            model_id: model_id.clone(),
            error,
        })?;
        let mut rebuilt = EnsembleMember::new(positions);
        rebuilt.set_cell(member.cell().copied());
        rebuilt
            .set_weight(member.weight())
            .map_err(|error| raw::MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
        let atom_data = remap_atom_data(
            member.atom_data(),
            &source_topology,
            &topology,
            &partition.source_atoms,
        )
        .map_err(|error| raw::MmcifEnsembleInterpretError::Model { model_id, error })?;
        rebuilt
            .set_atom_data(atom_data)
            .map_err(|error| raw::MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
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

    let definition_targets = source
        .definitions()
        .map(|(_, definition)| {
            builder
                .add_molecule_definition(definition.molecule())
                .map_err(interpret_error)
        })
        .collect::<Result<Vec<_>, _>>()?;
    for molecule in source.molecules() {
        let target_instance = builder
            .add_instance(definition_targets[molecule.definition_id().index()])
            .map_err(interpret_error)?;
        for (source_atom, _) in molecule.atoms() {
            let target_atom = InstanceAtomId::new(target_instance, source_atom.atom());
            source_atoms.push(source_atom);
            if atom_map.insert(source_atom, target_atom).is_some() {
                return Err(interpret_error(
                    "duplicate mmCIF source atom during topology reconstruction",
                ));
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
    Ok(PartitionedTopology {
        topology,
        source_atoms,
        atom_map,
    })
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
                .is_ok_and(|definition| definition.hierarchy().is_some())
        })
        .count();
    report.small_molecules = partition
        .topology
        .instances()
        .filter(|(id, _)| {
            partition
                .topology
                .definition_for_instance(*id)
                .is_ok_and(|definition| definition.hierarchy().is_none())
        })
        .count();
    report.solvent_molecules = report
        .instances
        .iter()
        .filter(|instance| {
            instance
                .entity_kinds()
                .contains(&raw::MmcifEntityKind::Water)
        })
        .count();
    Ok(report)
}

fn remap_positions(
    source: &Positions,
    source_topology: &Arc<Topology>,
    target_topology: &Arc<Topology>,
    source_atoms: &[InstanceAtomId],
) -> Result<Positions, raw::MmcifInterpretError> {
    let positions = source_atoms
        .iter()
        .copied()
        .map(|atom| {
            source
                .position(source_topology, atom)
                .map(|position| position.to_value())
                .map_err(interpret_error)
        })
        .collect::<Result<Vec<_>, _>>()?;
    Positions::new(target_topology, Quantity::new(positions, MODEL_LENGTH_UNIT))
        .map_err(interpret_error)
}

fn remap_atom_data(
    source: &AtomData,
    source_topology: &Arc<Topology>,
    target_topology: &Arc<Topology>,
    source_atoms: &[InstanceAtomId],
) -> Result<AtomData, raw::MmcifInterpretError> {
    let occupancies = source_atoms
        .iter()
        .copied()
        .map(|atom| {
            source
                .occupancy(source_topology, atom)
                .map_err(interpret_error)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let b_factors = source_atoms
        .iter()
        .copied()
        .map(|atom| {
            source
                .b_factor(source_topology, atom)
                .map_err(interpret_error)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut atom_data = AtomData::new(target_topology);
    atom_data
        .set_occupancies(occupancies)
        .map_err(interpret_error)?;
    atom_data
        .set_b_factors(Quantity::new(
            b_factors
                .into_iter()
                .map(|value| value.map(Quantity::to_value))
                .collect::<Vec<_>>(),
            crate::units::SQUARE_ANGSTROM,
        ))
        .map_err(interpret_error)?;
    Ok(atom_data)
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
        assert!(interpretation
            .topology()
            .molecules()
            .all(|molecule| molecule.molecule().validate_connected().is_ok()));
        assert_eq!(interpretation.report().instances().len(), 2);
        assert!(interpretation
            .report()
            .instances()
            .iter()
            .all(|instance| { instance.asym_ids() == ["A"] && instance.entity_ids() == ["1"] }));
        assert_eq!(interpretation.report().template_bonds_pending(), 0);
    }
}
