//! Biological assemblies that mmCIF and PDB entries state, and the models
//! built from them.

use kekule::diagnostics::Cancellation;
use kekule::geometry::Point3;
use kekule::macromolecular::{
    Assembly, AssemblyError, AssemblyOmissionReason, StructureInterpretIssue,
    StructureInterpretOptions, StructureInterpretation,
};
use kekule::units::{Quantity, ANGSTROM};
use kekule::{mmcif, pdb};

const INSULIN_RECEPTOR_KINASE_CIF: &str = include_str!("../fixtures/mmcif/1IRO.cif");
const INSULIN_RECEPTOR_KINASE_PDB: &str = include_str!("../fixtures/pdb/1IRO.pdb");
const ADENYLATE_KINASE: &str = include_str!("../fixtures/mmcif/1AKE.cif");

fn from_mmcif(text: &str) -> StructureInterpretation {
    mmcif::parse_str(text, Default::default())
        .unwrap()
        .interpret(StructureInterpretOptions::default())
        .unwrap()
}

fn from_pdb(text: &str) -> StructureInterpretation {
    pdb::parse_str(text, Default::default())
        .unwrap()
        .interpret(StructureInterpretOptions::default())
        .unwrap()
}

fn angstroms(point: Quantity<Point3>) -> Point3 {
    point.into_unit(ANGSTROM).unwrap().into_value()
}

fn chain_labels(model: &kekule::structure::Model) -> Vec<(String, Option<String>)> {
    model
        .topology()
        .chains()
        .map(|chain| {
            (
                chain.label_id().to_owned(),
                chain.author_id().map(str::to_owned),
            )
        })
        .collect()
}

fn assembly_issues(interpretation: &StructureInterpretation) -> Vec<&StructureInterpretIssue> {
    interpretation
        .report()
        .issues()
        .iter()
        .filter(|issue| matches!(issue, StructureInterpretIssue::AssemblyOmitted { .. }))
        .collect()
}

/// Copies sit where each operation puts the deposited atoms: rotated, then
/// translated, the first copy under the identity.
fn assert_copies_are_placed(
    interpretation: &StructureInterpretation,
    assembly: &Assembly,
    tolerance: f64,
) {
    let model = interpretation.model();
    let built = assembly.to_model(model, None).unwrap();
    let generator = &assembly.generators()[0];
    let operations = generator.operations();
    assert_eq!(built.atom_count(), model.atom_count() * operations.len());
    assert_eq!(
        assembly.atom_count(model.topology()),
        built.atom_count(),
        "the count matches the build"
    );
    let source = model.topology().atom_ids();
    let copies = built.topology().atom_ids();
    for (index, operation) in operations.iter().enumerate() {
        let rotation = operation.rotation();
        let shift = operation
            .translation()
            .into_unit(ANGSTROM)
            .unwrap()
            .into_value();
        for (offset, atom) in source.iter().enumerate().step_by(97) {
            let from = angstroms(model.position(*atom).unwrap());
            let turned =
                rotation.transform_vector(kekule::geometry::Vector3::new(from.x, from.y, from.z));
            let copy = angstroms(
                built
                    .position(copies[index * source.len() + offset])
                    .unwrap(),
            );
            for (got, want) in [
                (copy.x, turned.x + shift.x),
                (copy.y, turned.y + shift.y),
                (copy.z, turned.z + shift.z),
            ] {
                assert!(
                    (got - want).abs() < tolerance,
                    "operation {} atom {offset}: {got} against {want}",
                    operation.name()
                );
            }
        }
    }
}

#[test]
fn mmcif_assemblies_list_chains_and_crystal_operators() {
    let interpretation = from_mmcif(INSULIN_RECEPTOR_KINASE_CIF);
    assert!(assembly_issues(&interpretation).is_empty());
    let [assembly] = interpretation.report().assemblies() else {
        panic!("1IRO states one assembly");
    };
    assert_eq!(assembly.id(), "1");
    assert_eq!(assembly.details(), Some("author_defined_assembly"));
    assert_eq!(assembly.oligomeric_details(), Some("trimeric"));
    assert_eq!(assembly.oligomeric_count(), Some(3));
    let [generator] = assembly.generators() else {
        panic!("one generator");
    };
    assert_eq!(generator.chains(), ["A", "B", "C"]);
    let names = generator
        .operations()
        .iter()
        .map(|operation| operation.name())
        .collect::<Vec<_>>();
    assert_eq!(names, ["1", "2", "3"]);
    assert!(generator.operations()[0].is_identity());
    let shift = generator.operations()[1]
        .translation()
        .into_unit(ANGSTROM)
        .unwrap()
        .into_value();
    assert!((shift.x - 64.04).abs() < 1e-9 && shift.y == 0.0 && shift.z == 0.0);
    assert_copies_are_placed(&interpretation, assembly, 1e-9);
}

#[test]
fn copies_keep_identity_chains_and_name_the_others_by_operation() {
    let interpretation = from_mmcif(INSULIN_RECEPTOR_KINASE_CIF);
    let model = interpretation.model();
    let built = interpretation.report().assemblies()[0]
        .to_model(model, None)
        .unwrap();
    let source = chain_labels(model);
    let labels = chain_labels(&built);
    assert_eq!(labels.len(), source.len() * 3);
    assert_eq!(labels[..source.len()], source[..]);
    for (copy, (label, author)) in labels[source.len()..].iter().zip(source.iter().cycle()) {
        let suffix = if copy.0.ends_with("-2") { "-2" } else { "-3" };
        assert_eq!(copy.0, format!("{label}{suffix}"));
        assert_eq!(
            copy.1,
            author.as_ref().map(|author| format!("{author}{suffix}"))
        );
    }
    assert!(built.cell().is_none(), "copies are not a lattice");
    assert_eq!(
        built.topology().definition_count(),
        model.topology().definition_count(),
        "copies share their molecules' definitions"
    );
}

/// The chains of a model, sorted.
fn sorted_chains(model: &kekule::structure::Model) -> Vec<(String, Option<String>)> {
    let mut chains = chain_labels(model);
    chains.sort();
    chains
}

#[test]
fn written_assemblies_name_every_copy() {
    // mmCIF: the trimer's copies keep their names through writing, label and
    // author, so it reads back as built.
    let interpretation = from_mmcif(INSULIN_RECEPTOR_KINASE_CIF);
    let built = interpretation.report().assemblies()[0]
        .to_model(interpretation.model(), None)
        .unwrap();
    let (text, _) = mmcif::write(
        [mmcif::MmcifBlockSource::model(&built)],
        mmcif::MmcifWriteOptions::default(),
    )
    .unwrap();
    let read = from_mmcif(&text);
    assert_eq!(sorted_chains(read.model()), sorted_chains(&built));

    // PDB: the copies' chains, `A-2` and `A-3`, as segment identifiers,
    // which read back as chains.
    let interpretation = from_pdb(INSULIN_RECEPTOR_KINASE_PDB);
    let built = interpretation.report().assemblies()[0]
        .to_model(interpretation.model(), None)
        .unwrap();
    let options = pdb::PdbWriteOptions::default().with_chain_ids(pdb::PdbChainIdPolicy::SegmentId);
    let (text, _) = pdb::write(&built, options).unwrap();
    let read = from_pdb(&text);
    let mut labels = chain_labels(read.model())
        .into_iter()
        .map(|(label, _)| label)
        .collect::<Vec<_>>();
    labels.dedup();
    assert_eq!(labels, ["A", "A-2", "A-3"]);
    assert_eq!(read.model().atom_count(), built.atom_count());
}

#[test]
fn assemblies_pick_their_chains() {
    let interpretation = from_mmcif(ADENYLATE_KINASE);
    let model = interpretation.model();
    let assemblies = interpretation.report().assemblies();
    assert_eq!(
        assemblies.iter().map(Assembly::id).collect::<Vec<_>>(),
        ["1", "2"]
    );
    let built = assemblies
        .iter()
        .map(|assembly| assembly.to_model(model, None).unwrap())
        .collect::<Vec<_>>();
    assert_eq!(
        built[0].atom_count() + built[1].atom_count(),
        model.atom_count(),
        "the two assemblies partition the chains"
    );
    for (model, chains) in built.iter().zip([["A", "C", "E"], ["B", "D", "F"]]) {
        let labels = chain_labels(model)
            .into_iter()
            .map(|(label, _)| label)
            .collect::<Vec<_>>();
        assert_eq!(labels, chains);
    }
}

#[test]
fn pdb_biomolecules_read_as_assemblies() {
    let interpretation = from_pdb(INSULIN_RECEPTOR_KINASE_PDB);
    assert!(assembly_issues(&interpretation).is_empty());
    let [assembly] = interpretation.report().assemblies() else {
        panic!("1IRO states one biomolecule");
    };
    assert_eq!(assembly.id(), "1");
    assert_eq!(assembly.details(), Some("author_defined_assembly"));
    assert_eq!(assembly.oligomeric_details(), Some("TRIMERIC"));
    assert_eq!(assembly.oligomeric_count(), None);
    assert_eq!(assembly.generators()[0].chains(), ["A"]);
    assert_copies_are_placed(&interpretation, assembly, 1e-9);
    // Six decimals place the copies as the mmCIF entry's ten do, within
    // the coordinates' own precision.
    let reference = from_mmcif(INSULIN_RECEPTOR_KINASE_CIF);
    let stated = &reference.report().assemblies()[0].generators()[0].operations()[2];
    let rounded = &assembly.generators()[0].operations()[2];
    let [a, b, c] = stated.rotation().columns();
    let [x, y, z] = rounded.rotation().columns();
    for (stated, rounded) in [(a, x), (b, y), (c, z)] {
        assert!((stated.x - rounded.x).abs() < 1e-6);
        assert!((stated.y - rounded.y).abs() < 1e-6);
        assert!((stated.z - rounded.z).abs() < 1e-6);
    }
}

#[test]
fn pdb_biomolecules_stay_in_the_verbatim_records() {
    let document = pdb::parse_str(INSULIN_RECEPTOR_KINASE_PDB, Default::default()).unwrap();
    let remarks = document
        .other_records()
        .iter()
        .filter(|record| record.text().starts_with("REMARK 350"))
        .count();
    assert_eq!(
        remarks,
        INSULIN_RECEPTOR_KINASE_PDB
            .lines()
            .filter(|line| line.starts_with("REMARK 350"))
            .count()
    );
}

/// A two-chain block whose assemblies exercise the expression grammar and
/// every reason an assembly is left out.
fn block(assemblies: &str, operators: &str) -> String {
    format!(
        "data_toy
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
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
HETATM 1 O O HOH A 1 . 1.0 0.0 0.0 P 1
HETATM 2 O O HOH B 1 . 0.0 2.0 0.0 Q 1
#
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
{assemblies}
#
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 1 0 0 0 0 1 0 0 0 0 1 0
T 1 0 0 10 0 1 0 0 0 0 1 0
R 0 -1 0 0 1 0 0 0 0 0 1 0
{operators}
#
"
    )
}

#[test]
fn products_place_one_copy_per_combination() {
    let interpretation = from_mmcif(&block("1 '(1,R)(1,T)' A", ""));
    let assembly = &interpretation.report().assemblies()[0];
    let names = assembly.generators()[0]
        .operations()
        .iter()
        .map(|operation| operation.name())
        .collect::<Vec<_>>();
    assert_eq!(names, ["1-1", "1-T", "R-1", "R-T"]);
    let built = assembly.to_model(interpretation.model(), None).unwrap();
    let labels = chain_labels(&built);
    assert_eq!(
        labels,
        [
            ("A".to_owned(), Some("P".to_owned())),
            ("A-1-T".to_owned(), Some("P-1-T".to_owned())),
            ("A-R-1".to_owned(), Some("P-R-1".to_owned())),
            ("A-R-T".to_owned(), Some("P-R-T".to_owned())),
        ]
    );
    // R after T: the atom at x = 1 Å moves to x = 11 Å, then turns onto y.
    let atoms = built.topology().atom_ids();
    let moved = angstroms(built.position(atoms[3]).unwrap());
    assert!(moved.x.abs() < 1e-12 && (moved.y - 11.0).abs() < 1e-12);
}

#[test]
fn unreadable_assemblies_are_reported_not_listed() {
    let interpretation = from_mmcif(&block(
        "1 1 A\n2 '(1-3)' A\n3 M A\n4 S A\n5 '(1' A\n6 1 .\n7 U A",
        "M -1 0 0 0 0 1 0 0 0 0 1 0\nS 2 0 0 0 0 1 0 0 0 0 1 0\nU 1 0 0 ? 0 1 0 0 0 0 1 0",
    ));
    assert_eq!(
        interpretation
            .report()
            .assemblies()
            .iter()
            .map(Assembly::id)
            .collect::<Vec<_>>(),
        ["1"]
    );
    let reasons = assembly_issues(&interpretation)
        .into_iter()
        .map(|issue| match issue {
            StructureInterpretIssue::AssemblyOmitted {
                assembly, reason, ..
            } => (assembly.clone().unwrap(), reason.clone()),
            _ => unreachable!(),
        })
        .collect::<Vec<_>>();
    let operator = |id: &str| id.to_owned();
    assert_eq!(
        reasons,
        [
            (
                "2".to_owned(),
                AssemblyOmissionReason::UnknownOperator {
                    operator: operator("2")
                }
            ),
            (
                "3".to_owned(),
                AssemblyOmissionReason::NotRotation {
                    operator: operator("M")
                }
            ),
            (
                "4".to_owned(),
                AssemblyOmissionReason::NotRotation {
                    operator: operator("S")
                }
            ),
            (
                "5".to_owned(),
                AssemblyOmissionReason::Expression {
                    expression: "(1".to_owned()
                }
            ),
            ("6".to_owned(), AssemblyOmissionReason::NoGenerators),
            (
                "7".to_owned(),
                AssemblyOmissionReason::UnreadableOperator {
                    operator: operator("U")
                }
            ),
        ]
    );
    let message = assembly_issues(&interpretation)[0].to_string();
    assert!(
        message.starts_with("Assembly \"2\" at line ")
            && message.ends_with(" is not read. Operator \"2\" is not stated."),
        "{message}"
    );
}

#[test]
fn a_chain_copied_twice_under_one_name_is_refused() {
    let interpretation = from_mmcif(&block("1 1 A\n1 1,T A,B", ""));
    let assembly = &interpretation.report().assemblies()[0];
    assert_eq!(assembly.generators().len(), 2);
    let error = assembly.to_model(interpretation.model(), None).unwrap_err();
    assert!(
        matches!(&error, AssemblyError::RepeatedCopy { chain, .. } if chain == "A"),
        "{error:?}"
    );
    assert_eq!(
        error.to_string(),
        "Assembly \"1\" places two copies named chain \"A\"."
    );
}

#[test]
fn a_model_without_the_chains_builds_nothing() {
    let interpretation = from_mmcif(&block("1 1 Z", ""));
    let error = interpretation.report().assemblies()[0]
        .to_model(interpretation.model(), None)
        .unwrap_err();
    assert!(matches!(error, AssemblyError::NoAtoms { .. }), "{error:?}");
}

#[test]
fn building_stops_when_cancelled() {
    let interpretation = from_mmcif(INSULIN_RECEPTOR_KINASE_CIF);
    let cancellation = Cancellation::new(|| true);
    let error = interpretation.report().assemblies()[0]
        .to_model(interpretation.model(), Some(&cancellation))
        .unwrap_err();
    assert!(
        matches!(error, AssemblyError::Cancelled { .. }),
        "{error:?}"
    );
}

#[test]
fn reports_keep_assemblies_through_the_native_format() {
    let interpretation = from_mmcif(INSULIN_RECEPTOR_KINASE_CIF);
    let bytes = kekule::native::write(&interpretation, Default::default()).unwrap();
    let decoded = kekule::native::read(bytes.as_slice(), Default::default())
        .unwrap()
        .into_interpretation()
        .unwrap();
    assert_eq!(
        decoded.report().assemblies(),
        interpretation.report().assemblies()
    );
}
