//! Shared connected-component publication for subsets and structural editors.
use super::TopologyBuilder;
use crate::core::{
    AtomId, BondId, MoleculeEditor, MoleculeError, MoleculePublicationError, StereoCarrier,
    StereoElement, StereoElementKind, StereoGroup,
};
use crate::properties::PropertyError;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

#[derive(Debug, Clone, PartialEq)]
pub(super) enum ComponentBuildError {
    Molecule(MoleculeError),
    Publication(MoleculePublicationError),
    Property(PropertyError),
    Topology(super::TopologyBuildError),
}
impl From<MoleculeError> for ComponentBuildError {
    fn from(e: MoleculeError) -> Self {
        Self::Molecule(e)
    }
}
impl From<MoleculePublicationError> for ComponentBuildError {
    fn from(e: MoleculePublicationError) -> Self {
        Self::Publication(e)
    }
}
impl From<PropertyError> for ComponentBuildError {
    fn from(e: PropertyError) -> Self {
        Self::Property(e)
    }
}
impl From<super::TopologyBuildError> for ComponentBuildError {
    fn from(e: super::TopologyBuildError) -> Self {
        Self::Topology(e)
    }
}

pub(super) struct ComponentDefinition {
    pub(super) id: super::MoleculeDefinitionId,
    pub(super) source_atoms: Vec<AtomId>,
    pub(super) source_bonds: Vec<BondId>,
}

pub(super) fn build_component_definitions(
    source_molecule: &crate::core::Molecule,
    selected_local: &BTreeSet<AtomId>,
    class: Option<(super::MoleculeClass, bool)>,
    builder: &mut TopologyBuilder,
) -> Result<Vec<ComponentDefinition>, ComponentBuildError> {
    let whole_instance_selected = selected_local.len() == source_molecule.atom_count()
        && source_molecule.connected_components().len() == 1;
    let mut definitions = Vec::new();
    let mut visited = BTreeSet::new();
    let mut components = Vec::new();
    let mut local_atoms = BTreeMap::new();
    let mut local_bonds = BTreeMap::new();
    for seed in source_molecule.atom_ids() {
        if !selected_local.contains(&seed) || !visited.insert(seed) {
            continue;
        }
        let mut queue = VecDeque::from([seed]);
        let mut component = Vec::new();
        while let Some(atom) = queue.pop_front() {
            component.push(atom);
            for neighbor in source_molecule.neighbors(atom)? {
                if selected_local.contains(&neighbor) && visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }
        component.sort_unstable();
        let mut editor = MoleculeEditor::new();
        for &atom in &component {
            let target_atom = editor.add_atom(source_molecule.atom(atom)?.clone())?;
            local_atoms.insert(atom, (components.len(), target_atom));
        }
        components.push(ComponentDraft {
            editor,
            source_atoms: component,
            source_bonds: Vec::new(),
        });
    }

    // Scan source bonds and stereo once for all components, rather than
    // scanning or cloning the whole source molecule for every fragment.
    for (source_bond, bond) in source_molecule.bonds() {
        let (Some(&(component, a)), Some(&(other_component, b))) =
            (local_atoms.get(&bond.a()), local_atoms.get(&bond.b()))
        else {
            continue;
        };
        debug_assert_eq!(component, other_component);
        let target_bond = components[component].editor.add_bond(a, b, bond.order)?;
        components[component].source_bonds.push(source_bond);
        local_bonds.insert(source_bond, (component, target_bond));
    }
    let mut local_stereo = BTreeMap::new();
    for (source_element, element) in source_molecule.stereo_elements() {
        if let Some((component, element)) = remap_subset_stereo(element, &local_atoms, &local_bonds)
        {
            let target_element = components[component].editor.add_stereo_element(element)?;
            local_stereo.insert(source_element, (component, target_element));
        }
    }
    for (_, group) in source_molecule.stereo_groups() {
        let mut members_by_component = BTreeMap::<_, Vec<_>>::new();
        for member in &group.members {
            if let Some(&(component, target_element)) = local_stereo.get(member) {
                members_by_component
                    .entry(component)
                    .or_default()
                    .push(target_element);
            }
        }
        for (component, members) in members_by_component {
            components[component].editor.add_stereo_group(StereoGroup {
                kind: group.kind,
                members,
            })?;
        }
    }

    for mut component in components {
        let properties = &mut component.editor.working_mut().properties;
        *properties.atoms_mut() = source_molecule.atom_properties().select_indices(
            &component
                .source_atoms
                .iter()
                .map(|atom| atom.index())
                .collect::<Vec<_>>(),
        )?;
        *properties.bonds_mut() = source_molecule.bond_properties().select_indices(
            &component
                .source_bonds
                .iter()
                .map(|bond| bond.index())
                .collect::<Vec<_>>(),
        )?;
        if whole_instance_selected {
            for (key, value) in source_molecule.properties().iter() {
                properties.insert(key.clone(), value.clone())?;
            }
        }
        let target_molecule = component.editor.finish()?;
        let definition = builder.add_molecule_definition_owned(target_molecule)?;
        if whole_instance_selected {
            if let Some((class, explicit)) = class {
                builder.preserve_molecule_class(definition, class, explicit)?;
            }
        }
        definitions.push(ComponentDefinition {
            id: definition,
            source_atoms: component.source_atoms,
            source_bonds: component.source_bonds,
        });
    }
    Ok(definitions)
}

#[derive(Debug, Default)]
struct ComponentDraft {
    editor: MoleculeEditor,
    source_atoms: Vec<AtomId>,
    source_bonds: Vec<BondId>,
}

// Retain an assertion only when all its represented references survive in the
// same component, matching atom/bond deletion's stereo-pruning semantics.
fn remap_subset_stereo(
    element: &StereoElement,
    atoms: &BTreeMap<AtomId, (usize, AtomId)>,
    bonds: &BTreeMap<BondId, (usize, BondId)>,
) -> Option<(usize, StereoElement)> {
    let component = match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => atoms.get(&stereo.center)?.0,
        StereoElementKind::DoubleBond(stereo) => bonds.get(&stereo.bond)?.0,
        StereoElementKind::Axis(stereo) => bonds.get(&stereo.axis)?.0,
    };
    let atom = |source| {
        atoms
            .get(&source)
            .and_then(|&(owner, target)| (owner == component).then_some(target))
    };
    let bond = |source| {
        bonds
            .get(&source)
            .and_then(|&(owner, target)| (owner == component).then_some(target))
    };
    let carrier = |source| match source {
        StereoCarrier::Atom(source) => atom(source).map(StereoCarrier::Atom),
        other => Some(other),
    };
    let mut kind = element.kind.clone();
    match &mut kind {
        StereoElementKind::Tetrahedral(stereo) => {
            stereo.center = atom(stereo.center)?;
            for value in &mut stereo.carriers {
                *value = carrier(*value)?;
            }
        }
        StereoElementKind::DoubleBond(stereo) => {
            stereo.bond = bond(stereo.bond)?;
            stereo.left = atom(stereo.left)?;
            stereo.right = atom(stereo.right)?;
            stereo.left_carrier = carrier(stereo.left_carrier)?;
            stereo.right_carrier = carrier(stereo.right_carrier)?;
        }
        StereoElementKind::Axis(stereo) => {
            stereo.axis = bond(stereo.axis)?;
            for value in &mut stereo.carriers {
                *value = carrier(*value)?;
            }
        }
    }
    Some((component, StereoElement::new(kind)))
}
