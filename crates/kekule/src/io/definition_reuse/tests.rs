use super::DefinitionReuse;
use crate::molecule::Molecule;
use crate::structure::StructureBuilder;

fn molecule(smiles: &str) -> Molecule {
    crate::smiles::to_molecules(smiles).unwrap().remove(0)
}

#[test]
fn only_identical_molecules_with_equal_labels_share_a_definition() {
    let mut builder = StructureBuilder::new();
    let mut reuse = DefinitionReuse::new();
    let mut definition = |smiles: &str, label: &str| {
        let molecule = molecule(smiles);
        let labels = vec![label.to_owned(); molecule.atom_count()];
        reuse.definition(&mut builder, molecule, labels).unwrap()
    };
    let ethanol = definition("CCO", "EOH");
    assert_eq!(definition("CCO", "EOH"), ethanol);
    // Each difference gives a definition of its own: atom order, charge,
    // isotope, stated hydrogens, atom map, bond order, and labels.
    let different = [
        ("OCC", "EOH"),
        ("CC[O-]", "EOH"),
        ("CC[18OH]", "EOH"),
        ("CC[OH]", "EOH"),
        ("[CH3:1]CO", "EOH"),
        ("C=CO", "EOH"),
        ("CCO", "ETH"),
    ]
    .map(|(smiles, label)| definition(smiles, label));
    let mut all = different.to_vec();
    all.push(ethanol);
    all.sort();
    all.dedup();
    assert_eq!(all.len(), different.len() + 1);

    // Stereo is chemistry: enantiomers do not share, equal stereo does.
    let left = definition("C[C@H](N)O", "ALO");
    let right = definition("C[C@@H](N)O", "ALO");
    assert_ne!(left, right);
    assert_eq!(definition("C[C@H](N)O", "ALO"), left);
}
