use super::*;

#[test]
fn residues_coordinated_or_linked_to_a_polymer_are_entities_of_their_own() {
    use AssignedEntityKind::{Branched, NonPolymer, Polymer, Water};
    let monomer = ResidueRole {
        monomer: true,
        chained: false,
        glycan: None,
    };
    let chained = ResidueRole {
        monomer: true,
        chained: true,
        glycan: None,
    };
    let attached = ResidueRole {
        monomer: false,
        chained: false,
        glycan: None,
    };
    let glycan = |residues| ResidueRole {
        monomer: false,
        chained: false,
        glycan: Some((ResidueId::new(0), residues)),
    };
    for (molecule, class, role, expected) in [
        (Polymer, ResidueClass::AminoAcid, monomer, Polymer),
        (Polymer, ResidueClass::RnaNucleotide, monomer, Polymer),
        // A modified residue is a monomer of the sequence; a covalently
        // bound ligand is not.
        (Polymer, ResidueClass::Other, monomer, Polymer),
        (Polymer, ResidueClass::Other, attached, NonPolymer),
        (Polymer, ResidueClass::Ion, attached, NonPolymer),
        (Polymer, ResidueClass::Carbohydrate, glycan(1), NonPolymer),
        (Polymer, ResidueClass::Carbohydrate, glycan(3), Branched),
        (Polymer, ResidueClass::Water, attached, Water),
        (Branched, ResidueClass::Carbohydrate, glycan(2), Branched),
        (Branched, ResidueClass::Ion, attached, NonPolymer),
        (Branched, ResidueClass::Water, attached, Water),
        // A peptide of D-amino acids is no protein, but its chained
        // residues are monomers; a lone one is not.
        (NonPolymer, ResidueClass::Other, chained, Polymer),
        (NonPolymer, ResidueClass::Other, monomer, NonPolymer),
        // A metal complex keeps its coordinated water.
        (NonPolymer, ResidueClass::Water, attached, NonPolymer),
        (NonPolymer, ResidueClass::Ion, attached, NonPolymer),
        (Water, ResidueClass::Water, attached, Water),
    ] {
        assert_eq!(
            residue_entity_kind(molecule, class, role),
            expected,
            "{molecule:?} {class:?} {role:?}"
        );
    }
}
