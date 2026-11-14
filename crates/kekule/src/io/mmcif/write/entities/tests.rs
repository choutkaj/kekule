use super::*;

#[test]
fn residues_coordinated_or_linked_to_a_polymer_are_entities_of_their_own() {
    use AssignedEntityKind::{Branched, NonPolymer, Polymer, Water};
    for (molecule, class, expected) in [
        (Polymer, ResidueClass::AminoAcid, Polymer),
        (Polymer, ResidueClass::RnaNucleotide, Polymer),
        // A modified residue cannot be told from a residue's class alone.
        (Polymer, ResidueClass::Other, Polymer),
        (Polymer, ResidueClass::Ion, NonPolymer),
        (Polymer, ResidueClass::Carbohydrate, NonPolymer),
        (Polymer, ResidueClass::Water, Water),
        (Branched, ResidueClass::Carbohydrate, Branched),
        (Branched, ResidueClass::Ion, NonPolymer),
        (Branched, ResidueClass::Water, Water),
        // A metal complex keeps its coordinated water.
        (NonPolymer, ResidueClass::Water, NonPolymer),
        (NonPolymer, ResidueClass::Ion, NonPolymer),
        (Water, ResidueClass::Water, Water),
    ] {
        assert_eq!(
            residue_entity_kind(molecule, class),
            expected,
            "{molecule:?} {class:?}"
        );
    }
}
