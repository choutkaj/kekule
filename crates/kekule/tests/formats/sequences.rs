use std::sync::Arc;

use kekule::components::ComponentLibrary;
use kekule::macromolecular::{PolymerType, SequenceSource, StructureInterpretation};
use kekule::sequence::{self, CodeOrigin, SequenceOptions};
use kekule::{mmcif, pdb};

const TRNA: &str = include_str!("../fixtures/mmcif/1EHZ.cif");
const GFP: &str = include_str!("../fixtures/mmcif/1EMA.cif");

fn interpret(text: &str) -> StructureInterpretation {
    mmcif::parse_str(text, Default::default())
        .expect("parses")
        .interpret(Default::default())
        .expect("interprets")
}

#[test]
fn modified_nucleotides_take_their_parents_codes() {
    let interpretation = interpret(TRNA);
    let entities = interpretation.report().polymer_entities();
    assert_eq!(entities.len(), 1);
    let trna = &entities[0];
    assert_eq!(trna.id(), "1");
    assert_eq!(trna.polymer_type(), Some(&PolymerType::Polyribonucleotide));
    assert_eq!(trna.monomers_source(), SequenceSource::Stated);
    assert_eq!(trna.sequence_source(), SequenceSource::Stated);
    assert_eq!(
        trna.sequence(),
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
    );
    let monomer = |number: i32| {
        trna.monomers()
            .iter()
            .find(|monomer| monomer.number() == number)
            .expect("monomer")
    };
    for (number, component, code) in [(10, "2MG", "G"), (37, "YYG", "G"), (39, "PSU", "U")] {
        let monomer = monomer(number);
        assert_eq!(monomer.component_id(), component);
        assert_eq!(monomer.code().as_str(), code);
        assert_eq!(monomer.code().origin(), CodeOrigin::Component);
    }

    // Every monomer is observed, so the chain's sequence is the entity's.
    let instance = &trna.instances()[0];
    assert_eq!(instance.asym_id(), "A");
    assert!(instance
        .positions()
        .iter()
        .all(|position| position.residue().is_some()));
    let topology = interpretation.topology();
    let chains = sequence::chain_sequences(topology, &SequenceOptions::default());
    assert_eq!(chains.len(), 1);
    assert_eq!(Some(chains[0].chain()), instance.chain());
    assert_eq!(chains[0].sequence(), trna.sequence());

    let residue = topology
        .residue(instance.positions()[9].residue().expect("observed"))
        .expect("residue");
    assert_eq!(residue.name(), "2MG");
    let code = sequence::residue_code(residue, &SequenceOptions::default()).expect("code");
    assert_eq!(code.as_str(), "G");
    // Magnesium ions and water have no code.
    assert!(topology
        .residues()
        .filter(|residue| matches!(residue.name(), "MG" | "HOH"))
        .all(|residue| sequence::residue_code(residue, &SequenceOptions::default()).is_none()));
}

#[test]
fn stated_sequences_place_the_observed_residues() {
    let interpretation = interpret(GFP);
    let gfp = &interpretation.report().polymer_entities()[0];
    assert_eq!(gfp.polymer_type(), Some(&PolymerType::PolypeptideL));
    assert_eq!(gfp.monomers().len(), 236);
    assert!(gfp.sequence().starts_with("MSKGEELFTG"));

    // The file states the chromophore as X; its dictionary code is TYG.
    let chromophore = &gfp.monomers()[64];
    assert_eq!(
        (chromophore.number(), chromophore.component_id()),
        (65, "CRO")
    );
    assert_eq!(chromophore.code().as_str(), "X");
    assert_eq!(chromophore.code().origin(), CodeOrigin::Sequence);
    let selenomethionine = &gfp.monomers()[0];
    assert_eq!(selenomethionine.component_id(), "MSE");
    assert_eq!(selenomethionine.code().as_str(), "M");

    // The first monomer and the last nine are not observed; the scheme still
    // numbers them.
    let positions = gfp.instances()[0].positions();
    let unobserved = positions
        .iter()
        .zip(gfp.monomers())
        .filter(|(position, _)| position.residue().is_none())
        .map(|(position, monomer)| (monomer.number(), position.author_seq_id()))
        .collect::<Vec<_>>();
    assert_eq!(unobserved.len(), 10);
    assert_eq!(unobserved[0], (1, Some("1")));
    assert_eq!(unobserved[1], (228, Some("230")));
    assert_eq!(unobserved[9], (236, Some("238")));

    let topology = interpretation.topology();
    let second = topology
        .residue(positions[1].residue().expect("observed"))
        .expect("residue");
    assert_eq!((second.name(), second.label_seq_id()), ("SER", Some(2)));

    // The chain lists observed residues, the chromophore with its own code.
    let chains = sequence::chain_sequences(topology, &SequenceOptions::default());
    let chain = &chains[0];
    assert!(chain.sequence().starts_with("SKGEELFTG"));
    assert!(chain.sequence().contains("TFTYGVQCF"));
    assert_eq!(chain.monomers().len(), 226);
}

const TOY: &str = r"
data_toy
loop_
_entity.id
_entity.type
1 polymer
2 polymer
3 water
loop_
_entity_poly.entity_id
_entity_poly.type
1 'polypeptide(L)'
2 'polypeptide(L)'
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 GLY n
1 2 ALA y
1 2 SER y
1 3 XSR n
loop_
_chem_comp.id
_chem_comp.type
_chem_comp.one_letter_code
_chem_comp.mon_nstd_parent_comp_id
GLY 'peptide linking' G ?
ALA 'L-peptide linking' ? ?
SER 'L-peptide linking' ? ?
XSR 'L-peptide linking' ? SER
ZZZ 'L-peptide linking' ? ?
HOH non-polymer ? ?
loop_
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
_atom_site.auth_seq_id
_atom_site.auth_asym_id
1 C CA GLY A 1 1 0.0 0.0 0.0 1 A
2 C CA ALA A 1 2 3.8 0.0 0.0 2 A
3 C CA XSR A 1 3 7.6 0.0 0.0 3 A
4 C CA ZZZ B 2 5 0.0 10.0 0.0 15 B
5 C CA GLY B 2 6 3.8 10.0 0.0 16 B
6 O O HOH C 3 . 0.0 20.0 0.0 1 A
";

#[test]
fn stated_components_and_derived_monomers() {
    let interpretation = interpret(TOY);
    let entities = interpretation.report().polymer_entities();
    assert_eq!(entities.len(), 2);

    let first = &entities[0];
    let codes = first
        .monomers()
        .iter()
        .map(|monomer| {
            (
                monomer.component_id(),
                monomer.code().as_str(),
                monomer.code().origin(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        codes,
        [
            ("GLY", "G", CodeOrigin::Component),
            ("ALA", "A", CodeOrigin::Component),
            ("XSR", "S", CodeOrigin::Parent),
        ]
    );
    assert_eq!(first.monomers()[1].alternatives(), ["SER"]);
    assert_eq!(first.sequence(), "GAS");
    assert_eq!(first.sequence_source(), SequenceSource::Derived);
    let instance = &first.instances()[0];
    assert_eq!(instance.author_chain_id(), Some("A"));
    assert_eq!(instance.positions()[2].author_seq_id(), Some("3"));

    // Without `_entity_poly_seq`, the atom sites give the monomers.
    let second = &entities[1];
    assert_eq!(second.monomers_source(), SequenceSource::Derived);
    assert_eq!(
        second
            .monomers()
            .iter()
            .map(|monomer| monomer.number())
            .collect::<Vec<_>>(),
        [5, 6]
    );
    assert_eq!(second.sequence(), "XG");
    assert!(second.monomers()[0].code().is_unknown());
    assert!(second.instances()[0]
        .positions()
        .iter()
        .all(|position| position.residue().is_some()));
}

#[test]
fn atom_sites_name_instances_a_partial_struct_asym_omits() {
    let text = TOY.replace(
        "1 3 XSR n
",
        "1 3 XSR n
2 5 ZZZ n
2 6 GLY n
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
",
    );
    let interpretation = interpret(&text);
    let entities = interpretation.report().polymer_entities();
    assert_eq!(entities[1].monomers_source(), SequenceSource::Stated);
    let instance = &entities[1].instances()[0];
    assert_eq!(instance.asym_id(), "B");
    assert!(instance.chain().is_some());
    assert!(instance
        .positions()
        .iter()
        .all(|position| position.residue().is_some()));
}

#[test]
fn chain_sequences_of_a_pdb_entry() {
    let text = "\
ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  MSE A   2       3.800   0.000   0.000  1.00  0.00           C
ATOM      3  CA  TPO A   3       7.600   0.000   0.000  1.00  0.00           C
TER
HETATM    4  O   HOH A 101       0.000  10.000   0.000  1.00  0.00           O
HETATM    5  PG  GTP A 102       0.000  20.000   0.000  1.00  0.00           P
END
";
    let interpretation = pdb::parse_str(text, Default::default())
        .expect("parses")
        .interpret(Default::default())
        .expect("interprets");
    assert!(interpretation.report().polymer_entities().is_empty());
    let chains = sequence::chain_sequences(interpretation.topology(), &SequenceOptions::default());
    assert_eq!(chains.len(), 1);
    assert_eq!(chains[0].sequence(), "GMT");
    assert_eq!(chains[0].monomers().len(), 3);
}

#[test]
fn caller_libraries_name_force_field_residues() {
    let amber =
        SequenceOptions::default().with_components(Some(Arc::new(ComponentLibrary::amber())));
    for (residue, code) in [
        ("HID", "H"),
        ("CYX", "C"),
        ("ASH", "D"),
        ("NALA", "A"),
        ("DA5", "A"),
    ] {
        let found = sequence::component_code(residue, &amber).expect(residue);
        assert_eq!(found.as_str(), code, "{residue}");
    }
    assert!(sequence::component_code("WAT", &amber).is_none());
}
