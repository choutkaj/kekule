use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::core::{BondOrder, HydrogenDeclaration, Molecule};

use super::{
    Hierarchy, InstanceAtomId, MoleculeClass, MoleculeDefinition, MoleculeDefinitionId,
    MoleculeInstance, ResidueClass, ResidueId,
};

pub(super) fn finalize(
    definitions: &mut [MoleculeDefinition],
    instances: &[MoleculeInstance],
    hierarchy: &mut Hierarchy,
    molecule_overrides: &BTreeMap<MoleculeDefinitionId, MoleculeClass>,
    residue_overrides: &BTreeMap<ResidueId, ResidueClass>,
) {
    let atom_residues = hierarchy
        .atom_sites()
        .map(|(_, site)| (site.atom(), site.residue()))
        .collect::<BTreeMap<_, _>>();
    let atom_names = hierarchy
        .atom_sites()
        .filter_map(|(_, site)| {
            let name = site
                .metadata()
                .label_atom_id
                .as_deref()
                .or(site.metadata().auth_atom_id.as_deref())?;
            Some((site.atom(), normalize_atom_name(name)))
        })
        .collect::<BTreeMap<_, _>>();

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

    let residue_classes = hierarchy
        .residues()
        .map(|(id, residue)| (id, residue.class()))
        .collect::<BTreeMap<_, _>>();
    let mut evidence = vec![BTreeSet::new(); definitions.len()];
    for instance in instances {
        let molecule = definitions[instance.definition().index()].molecule();
        if let Some(class) = infer_instance_class(
            instance,
            molecule,
            &atom_residues,
            &atom_names,
            &residue_classes,
        ) {
            evidence[instance.definition().index()].insert(class);
        }
    }

    for definition in definitions {
        definition.class = molecule_overrides
            .get(&definition.id())
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
        .filter_map(|site| hierarchy.atom_site(*site).ok())
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
) -> Option<&'a crate::core::Atom> {
    let instance = instances.get(atom.molecule().index())?;
    definitions
        .get(instance.definition().index())?
        .molecule()
        .atom(atom.atom())
        .ok()
}

fn infer_instance_class(
    instance: &MoleculeInstance,
    molecule: &Molecule,
    atom_residues: &BTreeMap<InstanceAtomId, ResidueId>,
    atom_names: &BTreeMap<InstanceAtomId, String>,
    residue_classes: &BTreeMap<ResidueId, ResidueClass>,
) -> Option<MoleculeClass> {
    let qualified_atoms = molecule
        .atom_ids()
        .map(|atom| InstanceAtomId::new(instance.id(), atom))
        .collect::<Vec<_>>();
    let touched_residues = qualified_atoms
        .iter()
        .filter_map(|atom| atom_residues.get(atom).copied())
        .collect::<BTreeSet<_>>();

    if !touched_residues.is_empty()
        && qualified_atoms.iter().all(|atom| {
            atom_residues
                .get(atom)
                .is_some_and(|residue| residue_classes.get(residue) == Some(&ResidueClass::Water))
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
        && qualified_atoms.iter().all(|atom| {
            atom_residues
                .get(atom)
                .is_some_and(|residue| residue_classes.get(residue) == Some(&ResidueClass::Ion))
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
        let (Some(&left_residue), Some(&right_residue)) =
            (atom_residues.get(&left), atom_residues.get(&right))
        else {
            continue;
        };
        if left_residue == right_residue {
            continue;
        }
        let left_name = atom_names.get(&left).map(String::as_str);
        let right_name = atom_names.get(&right).map(String::as_str);
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

    let protein =
        component_has_at_least(&peptide_edges, residue_classes, ResidueClass::AminoAcid, 2);
    let (dna, rna, nucleic_conflict) = nucleic_evidence(&nucleotide_edges, residue_classes);
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
        .any(|residue| residue_classes.get(residue) == Some(&ResidueClass::Carbohydrate))
    {
        return Some(MoleculeClass::Carbohydrate);
    }
    None
}

fn is_explicit_water(molecule: &Molecule) -> bool {
    let atoms = molecule.atoms().map(|(_, atom)| atom).collect::<Vec<_>>();
    is_explicit_water_atoms(&atoms)
}

fn is_explicit_water_atoms(atoms: &[&crate::core::Atom]) -> bool {
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
        .map(|atom| match atom.hydrogens {
            HydrogenDeclaration::Fixed(count) => usize::from(count),
            HydrogenDeclaration::Infer { explicit } => usize::from(explicit),
        })
        .sum::<usize>();
    oxygen_count == 1
        && atoms
            .iter()
            .all(|atom| matches!(atom.element.atomic_number(), 1 | 8))
        && graph_hydrogens + declared_hydrogens == 2
}

fn component_has_at_least(
    edges: &[(ResidueId, ResidueId)],
    classes: &BTreeMap<ResidueId, ResidueClass>,
    class: ResidueClass,
    minimum: usize,
) -> bool {
    residue_components(edges).into_iter().any(|component| {
        component
            .iter()
            .filter(|residue| classes.get(residue) == Some(&class))
            .count()
            >= minimum
    })
}

fn nucleic_evidence(
    edges: &[(ResidueId, ResidueId)],
    classes: &BTreeMap<ResidueId, ResidueClass>,
) -> (bool, bool, bool) {
    let mut dna = false;
    let mut rna = false;
    let mut conflict = false;
    for component in residue_components(edges) {
        let dna_count = component
            .iter()
            .filter(|residue| classes.get(residue) == Some(&ResidueClass::DnaNucleotide))
            .count();
        let rna_count = component
            .iter()
            .filter(|residue| classes.get(residue) == Some(&ResidueClass::RnaNucleotide))
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
