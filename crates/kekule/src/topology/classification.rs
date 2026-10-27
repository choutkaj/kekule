use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::molecule::{BondOrder, Molecule};

use super::{
    Hierarchy, InstanceAtomId, MoleculeClass, MoleculeDefinition, MoleculeDefinitionId,
    MoleculeInstance, ResidueClass, ResidueId,
};

pub(super) fn finalize(
    definitions: &mut [MoleculeDefinition],
    instances: &[MoleculeInstance],
    hierarchy: &mut Hierarchy,
    molecule_overrides: &BTreeMap<MoleculeDefinitionId, MoleculeClass>,
    preserved_molecule_classes: &BTreeMap<MoleculeDefinitionId, MoleculeClass>,
    appended_instances_start: Option<usize>,
    residue_overrides: &BTreeMap<ResidueId, ResidueClass>,
) {
    for index in 0..hierarchy.residues.len() {
        let residue_id = hierarchy.residues[index].id();
        let class = residue_overrides
            .get(&residue_id)
            .copied()
            .unwrap_or_else(|| {
                infer_residue_class(
                    &hierarchy.residues[index],
                    hierarchy,
                    definitions,
                    instances,
                )
            });
        hierarchy.residues[index].class = class;
    }

    // Only definitions that may be inferred need evidence: explicit classes
    // win, and a preserved class is reconsidered only when an appended
    // occurrence brings evidence. Unchanged definitions of a large system
    // cost nothing.
    let inferable =
        |definition: MoleculeDefinitionId| !molecule_overrides.contains_key(&definition);
    let mut has_new_evidence = vec![false; definitions.len()];
    if let Some(start) = appended_instances_start {
        for instance in instances.iter().skip(start) {
            let definition = instance.definition();
            if inferable(definition)
                && preserved_molecule_classes.contains_key(&definition)
                && !has_new_evidence[definition.index()]
            {
                let molecule = definitions[definition.index()].molecule();
                has_new_evidence[definition.index()] =
                    infer_instance_class(instance, molecule, hierarchy).is_some();
            }
        }
    }
    let mut evidence = vec![BTreeSet::new(); definitions.len()];
    for instance in instances {
        let definition = instance.definition();
        if !inferable(definition)
            || (preserved_molecule_classes.contains_key(&definition)
                && !has_new_evidence[definition.index()])
        {
            continue;
        }
        let molecule = definitions[definition.index()].molecule();
        if let Some(class) = infer_instance_class(instance, molecule, hierarchy) {
            evidence[definition.index()].insert(class);
        }
    }

    for definition in definitions {
        definition.class = molecule_overrides
            .get(&definition.id())
            .or_else(|| {
                // A resumed builder's new informative occurrence can change or
                // conflict with the evidence for a reused definition. Explicit
                // assignments still win; ordinary appends preserve cached classes.
                preserved_molecule_classes
                    .get(&definition.id())
                    .filter(|_| !has_new_evidence[definition.id().index()])
            })
            .copied()
            .unwrap_or_else(|| {
                let classes = &evidence[definition.id().index()];
                match classes.len() {
                    0 => MoleculeClass::SmallMolecule,
                    1 => *classes.iter().next().expect("one class is present"),
                    _ => MoleculeClass::Other,
                }
            });
    }
}

fn infer_residue_class(
    residue: &super::Residue,
    hierarchy: &Hierarchy,
    definitions: &[MoleculeDefinition],
    instances: &[MoleculeInstance],
) -> ResidueClass {
    let component = residue
        .label_comp_id()
        .or(residue.author_comp_id())
        .unwrap_or_else(|| residue.name())
        .trim()
        .to_ascii_uppercase();
    if AMINO_ACIDS.contains(&component.as_str()) {
        return ResidueClass::AminoAcid;
    }
    if DNA_NUCLEOTIDES.contains(&component.as_str()) {
        return ResidueClass::DnaNucleotide;
    }
    if RNA_NUCLEOTIDES.contains(&component.as_str()) {
        return ResidueClass::RnaNucleotide;
    }
    if CARBOHYDRATES.contains(&component.as_str()) {
        return ResidueClass::Carbohydrate;
    }
    if WATER_COMPONENTS.contains(&component.as_str()) {
        return ResidueClass::Water;
    }

    let atoms = residue
        .atom_sites()
        .iter()
        .filter_map(|site| hierarchy.atom_site(*site))
        .filter_map(|site| topology_atom(site.atom(), definitions, instances))
        .collect::<Vec<_>>();
    if atoms.len() == 1
        && (atoms[0].formal_charge != 0 || ION_COMPONENTS.contains(&component.as_str()))
    {
        return ResidueClass::Ion;
    }
    if is_explicit_water_atoms(&atoms) {
        return ResidueClass::Water;
    }
    ResidueClass::Other
}

fn topology_atom<'a>(
    atom: InstanceAtomId,
    definitions: &'a [MoleculeDefinition],
    instances: &[MoleculeInstance],
) -> Option<&'a crate::molecule::Atom> {
    let instance = instances.get(atom.instance().index())?;
    definitions
        .get(instance.definition().index())?
        .molecule()
        .atom(atom.atom())
}

fn infer_instance_class(
    instance: &MoleculeInstance,
    molecule: &Molecule,
    hierarchy: &Hierarchy,
) -> Option<MoleculeClass> {
    let residue_of = |atom: InstanceAtomId| {
        hierarchy
            .atom_site_for_atom(atom)
            .map(|site| site.residue())
    };
    let name_of = |atom: InstanceAtomId| {
        let site = hierarchy.atom_site_for_atom(atom)?;
        let name = site
            .metadata()
            .label_atom_id
            .as_deref()
            .or(site.metadata().auth_atom_id.as_deref())?;
        Some(normalize_atom_name(name))
    };
    let class_of = |residue: ResidueId| hierarchy.residues[residue.index()].class();
    let qualified_atoms = molecule
        .atom_ids()
        .map(|atom| InstanceAtomId::new(instance.id(), atom))
        .collect::<Vec<_>>();
    let touched_residues = qualified_atoms
        .iter()
        .filter_map(|&atom| residue_of(atom))
        .collect::<BTreeSet<_>>();

    if !touched_residues.is_empty()
        && qualified_atoms.iter().all(|&atom| {
            residue_of(atom).is_some_and(|residue| class_of(residue) == ResidueClass::Water)
        })
    {
        return Some(MoleculeClass::Water);
    }
    if is_explicit_water(molecule) {
        return Some(MoleculeClass::Water);
    }
    if molecule.atom_count() == 1
        && molecule
            .atoms()
            .next()
            .is_some_and(|(_, atom)| atom.formal_charge != 0)
    {
        return Some(MoleculeClass::Ion);
    }
    if !touched_residues.is_empty()
        && qualified_atoms.iter().all(|&atom| {
            residue_of(atom).is_some_and(|residue| class_of(residue) == ResidueClass::Ion)
        })
    {
        return Some(MoleculeClass::Ion);
    }

    let mut peptide_edges = Vec::new();
    let mut nucleotide_edges = Vec::new();
    for (_, bond) in molecule.bonds() {
        if bond.order != BondOrder::Single {
            continue;
        }
        let left = InstanceAtomId::new(instance.id(), bond.a());
        let right = InstanceAtomId::new(instance.id(), bond.b());
        let (Some(left_residue), Some(right_residue)) = (residue_of(left), residue_of(right))
        else {
            continue;
        };
        if left_residue == right_residue {
            continue;
        }
        let (left_name, right_name) = (name_of(left), name_of(right));
        let (left_name, right_name) = (left_name.as_deref(), right_name.as_deref());
        if matches!(
            (left_name, right_name),
            (Some("C"), Some("N")) | (Some("N"), Some("C"))
        ) {
            peptide_edges.push((left_residue, right_residue));
        }
        if matches!(
            (left_name, right_name),
            (Some("O3'"), Some("P")) | (Some("P"), Some("O3'"))
        ) {
            nucleotide_edges.push((left_residue, right_residue));
        }
    }

    let protein = component_has_at_least(&peptide_edges, &class_of, ResidueClass::AminoAcid, 2);
    let (dna, rna, nucleic_conflict) = nucleic_evidence(&nucleotide_edges, &class_of);
    if nucleic_conflict || (protein && (dna || rna)) {
        return Some(MoleculeClass::Other);
    }
    if protein {
        return Some(MoleculeClass::Protein);
    }
    if dna && rna {
        return Some(MoleculeClass::Other);
    }
    if dna {
        return Some(MoleculeClass::Dna);
    }
    if rna {
        return Some(MoleculeClass::Rna);
    }
    if touched_residues
        .iter()
        .any(|&residue| class_of(residue) == ResidueClass::Carbohydrate)
    {
        return Some(MoleculeClass::Carbohydrate);
    }
    None
}

fn is_explicit_water(molecule: &Molecule) -> bool {
    let atoms = molecule.atoms().map(|(_, atom)| atom).collect::<Vec<_>>();
    is_explicit_water_atoms(&atoms)
}

fn is_explicit_water_atoms(atoms: &[&crate::molecule::Atom]) -> bool {
    let oxygen_count = atoms
        .iter()
        .filter(|atom| atom.element.atomic_number() == 8)
        .count();
    let graph_hydrogens = atoms
        .iter()
        .filter(|atom| atom.element.atomic_number() == 1)
        .count();
    let declared_hydrogens = atoms
        .iter()
        .filter(|atom| atom.element.atomic_number() == 8)
        .map(|atom| usize::from(atom.hydrogens.represented_count()))
        .sum::<usize>();
    oxygen_count == 1
        && atoms
            .iter()
            .all(|atom| matches!(atom.element.atomic_number(), 1 | 8))
        && graph_hydrogens + declared_hydrogens == 2
}

fn component_has_at_least(
    edges: &[(ResidueId, ResidueId)],
    class_of: &impl Fn(ResidueId) -> ResidueClass,
    class: ResidueClass,
    minimum: usize,
) -> bool {
    residue_components(edges).into_iter().any(|component| {
        component
            .iter()
            .filter(|&&residue| class_of(residue) == class)
            .count()
            >= minimum
    })
}

fn nucleic_evidence(
    edges: &[(ResidueId, ResidueId)],
    class_of: &impl Fn(ResidueId) -> ResidueClass,
) -> (bool, bool, bool) {
    let mut dna = false;
    let mut rna = false;
    let mut conflict = false;
    for component in residue_components(edges) {
        let dna_count = component
            .iter()
            .filter(|&&residue| class_of(residue) == ResidueClass::DnaNucleotide)
            .count();
        let rna_count = component
            .iter()
            .filter(|&&residue| class_of(residue) == ResidueClass::RnaNucleotide)
            .count();
        if dna_count > 0 && rna_count > 0 && dna_count + rna_count >= 2 {
            conflict = true;
        }
        dna |= dna_count >= 2;
        rna |= rna_count >= 2;
    }
    (dna, rna, conflict)
}

fn residue_components(edges: &[(ResidueId, ResidueId)]) -> Vec<BTreeSet<ResidueId>> {
    let mut adjacency = BTreeMap::<ResidueId, BTreeSet<ResidueId>>::new();
    for &(left, right) in edges {
        adjacency.entry(left).or_default().insert(right);
        adjacency.entry(right).or_default().insert(left);
    }
    let mut visited = BTreeSet::new();
    let mut components = Vec::new();
    for &seed in adjacency.keys() {
        if !visited.insert(seed) {
            continue;
        }
        let mut component = BTreeSet::new();
        let mut queue = VecDeque::from([seed]);
        while let Some(residue) = queue.pop_front() {
            component.insert(residue);
            for &neighbor in &adjacency[&residue] {
                if visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }
        components.push(component);
    }
    components
}

fn normalize_atom_name(name: &str) -> String {
    name.trim().to_ascii_uppercase().replace(['*', '’'], "'")
}

const AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];
const DNA_NUCLEOTIDES: &[&str] = &["DA", "DC", "DG", "DT"];
const RNA_NUCLEOTIDES: &[&str] = &["A", "C", "G", "U"];
const WATER_COMPONENTS: &[&str] = &["HOH", "WAT", "H2O", "DOD"];
const CARBOHYDRATES: &[&str] = &[
    "GLC", "BGC", "GAL", "GLA", "MAN", "BMA", "FUC", "NAG", "NDG", "SIA",
];
const ION_COMPONENTS: &[&str] = &[
    "NA", "K", "CL", "CA", "MG", "ZN", "FE", "MN", "CU", "CO", "NI", "CD", "HG", "BR", "IOD",
];
