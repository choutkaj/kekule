//! Interpretation reports re-keyed to structures edited from their source.
use std::collections::BTreeMap;

use super::types::{
    AtomProvenance, BondProvenance, CoordinationProvenance, InstanceProvenance,
    StructureInterpretationReport,
};
use crate::structure::StructureAppend;
use crate::topology::edit::EditCorrespondence;
use crate::topology::{InstanceAtomId, MoleculeInstanceId, Topology};

impl StructureInterpretationReport {
    /// This report re-keyed to `edited`, which an editor published from
    /// `source` with `correspondence`, as
    /// [`StructureEditor::finish_with_correspondence`](crate::structure::StructureEditor::finish_with_correspondence)
    /// returns them. Writers given the result state the source's coordination
    /// and entities for the edited structure.
    ///
    /// Each surviving atom's provenance moves to its published atom and
    /// instance and takes the chain, residue, and atom identifiers of its
    /// edited atom site. Its source line and order, atom-site ID, residue key,
    /// occurrence, entity, and alternate location stay source facts. Deleted
    /// atoms lose their provenance, and added atoms have none. Bond provenance
    /// stays for bonds the edited structure keeps, and a coordination record
    /// while both partners survive. The instance counts follow the edited
    /// instances; the rest of the report, issues included, describes the
    /// source.
    #[must_use]
    pub fn follow_edit(
        &self,
        source: &Topology,
        edited: &Topology,
        correspondence: &EditCorrespondence,
    ) -> Self {
        self.follow(edited, |atom| {
            let index = correspondence.atoms().published(source.atom_index(atom)?)?;
            Some(edited.atom_at(index)?.id())
        })
    }

    /// This report re-keyed to `edited`, which an editor published with
    /// `correspondence` after appending the structure the report describes
    /// with `append`. Only the appended atoms have provenance: the PDB writer
    /// states their coordination, and the mmCIF writer refuses the other
    /// atoms, which lack provenance. Otherwise as [`Self::follow_edit`].
    #[must_use]
    pub fn follow_append(
        &self,
        append: &StructureAppend,
        edited: &Topology,
        correspondence: &EditCorrespondence,
    ) -> Self {
        self.follow(edited, |atom| {
            Some(correspondence.atom(append.atom(atom)?)?.0)
        })
    }

    /// This report re-keyed to `edited`, with `published` naming each source
    /// atom's published atom.
    fn follow(
        &self,
        edited: &Topology,
        published: impl Fn(InstanceAtomId) -> Option<InstanceAtomId>,
    ) -> Self {
        let mut instances = BTreeMap::<MoleculeInstanceId, InstanceProvenance>::new();
        for instance in &self.instances {
            for atom in &instance.atoms {
                let Some(id) = published(atom.atom) else {
                    continue;
                };
                let target = instances
                    .entry(id.instance())
                    .or_insert_with(|| InstanceProvenance {
                        molecule: id.instance(),
                        coordinate_model_id: instance.coordinate_model_id.clone(),
                        asym_ids: Vec::new(),
                        entity_ids: Vec::new(),
                        entity_kinds: Vec::new(),
                        atoms: Vec::new(),
                        bonds: Vec::new(),
                    });
                for kind in &instance.entity_kinds {
                    if !target.entity_kinds.contains(kind) {
                        target.entity_kinds.push(kind.clone());
                    }
                }
                target.atoms.push(followed_atom(atom, id, edited));
            }
            for bond in &instance.bonds {
                let [Some(a), Some(b)] = bond.atoms.map(&published) else {
                    continue;
                };
                let kept = a.instance() == b.instance()
                    && edited.instance(a.instance()).is_some_and(|instance| {
                        instance
                            .molecule()
                            .bond_between(a.atom(), b.atom())
                            .is_some()
                    });
                if let Some(target) = instances.get_mut(&a.instance()).filter(|_| kept) {
                    target.bonds.push(BondProvenance {
                        atoms: [a, b],
                        evidence: bond.evidence.clone(),
                    });
                }
            }
        }
        let instances = instances
            .into_values()
            .map(|mut instance| {
                instance
                    .atoms
                    .sort_by_key(|atom| edited.atom_index(atom.atom));
                instance.asym_ids = distinct(instance.atoms.iter().map(|atom| &atom.asym_id));
                instance.entity_ids = distinct(
                    instance
                        .atoms
                        .iter()
                        .filter_map(|atom| atom.entity_id.as_ref()),
                );
                instance
            })
            .collect();
        let coordination = self
            .coordination
            .iter()
            .filter_map(|record| {
                let [a, b] = record.atoms;
                Some(CoordinationProvenance {
                    atoms: [published(a)?, published(b)?],
                    connection_id: record.connection_id.clone(),
                    connection_type: record.connection_type.clone(),
                    source_line: record.source_line,
                })
            })
            .collect();
        let mut followed = Self {
            block_name: self.block_name.clone(),
            entity_definitions: self.entity_definitions,
            coordinate_models: self.coordinate_models,
            selected_model: self.selected_model.clone(),
            ignored_coordinate_models: self.ignored_coordinate_models.clone(),
            macromolecules: 0,
            small_molecules: 0,
            solvent_molecules: 0,
            applied_connections: self.applied_connections,
            connectivity_candidates: self.connectivity_candidates,
            instances,
            coordination,
            alternate_locations: self.alternate_locations.clone(),
            issues: self.issues.clone(),
            polymer_entities: self.polymer_entities.clone(),
            assemblies: self.assemblies.clone(),
        };
        followed.count_instances();
        followed
    }
}

/// `atom`'s provenance for published atom `id`, with the identifiers of its
/// edited atom site.
fn followed_atom(atom: &AtomProvenance, id: InstanceAtomId, edited: &Topology) -> AtomProvenance {
    let mut followed = atom.clone();
    followed.atom = id;
    let Some(site) = edited.atom(id).and_then(|atom| atom.atom_site()) else {
        return followed;
    };
    let owned = |value: Option<&str>| value.map(str::to_owned);
    let metadata = site.metadata();
    let residue = site.residue();
    followed.asym_id = residue.chain().label_id().to_owned();
    followed.label_asym_id = owned(metadata.label_asym_id());
    followed.auth_asym_id = owned(metadata.auth_asym_id());
    followed.label_sequence_id = residue.label_seq_id();
    followed.author_sequence_id = owned(residue.author_seq_id());
    followed.insertion_code = owned(residue.insertion_code());
    followed.component_id = residue.name().to_owned();
    followed.label_component_id = owned(residue.label_comp_id());
    followed.auth_component_id = owned(residue.author_comp_id());
    followed.label_atom_name = owned(metadata.label_atom_id());
    followed.auth_atom_name = owned(metadata.auth_atom_id());
    if let Some(name) = metadata.label_atom_id().or(metadata.auth_atom_id()) {
        followed.atom_name = name.to_owned();
    }
    if let Some(symbol) = metadata.type_symbol() {
        followed.type_symbol = symbol.to_owned();
    }
    followed
}

/// The distinct values in first-seen order.
fn distinct<'a>(values: impl Iterator<Item = &'a String>) -> Vec<String> {
    let mut seen = Vec::<String>::new();
    for value in values {
        if !seen.contains(value) {
            seen.push(value.clone());
        }
    }
    seen
}
