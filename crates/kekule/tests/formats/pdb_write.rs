use kekule::diagnostics::Countable;
use kekule::diagnostics::LengthMismatch;
use std::collections::BTreeSet;
use std::sync::Arc;

use kekule::diagnostics::{
    AtomField, CellReason, DiagnosticKind, FormatErrorDetail, GeometryField, UnrepresentablePolicy,
    WriteOmission,
};
use kekule::macromolecular::{
    ModelSelection, StructureInterpretOptions, StructureInterpretation,
    StructureInterpretationReport,
};
use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{AtomId, BondOrder};
use kekule::pdb::{self, PdbChainIdPolicy, PdbSource, PdbWriteErrorKind, PdbWriteOptions};
use kekule::solvation::{self, BoxShape, PeriodicBoxOptions, SolventOptions};
use kekule::structure::{
    Geometry, Order, Positions, Structure, StructureBuilder, StructureSet, Velocities,
};
use kekule::topology::{InstanceAtomId, Topology};
use kekule::units::{Quantity, ANGSTROM, DEGREE, MOLAR, NANOMETER};
use kekule::{mmcif, smiles};

const CRAMBIN: &str = include_str!("../fixtures/pdb/1CRN.pdb");
const ADENYLATE_KINASE: &str = include_str!("../fixtures/pdb/1AKE.pdb");
const RUBREDOXIN: &str = include_str!("../fixtures/pdb/1IRO.pdb");
const KRAS: &str = include_str!("../fixtures/mmcif/9UOH.cif");

fn interpret(text: &str) -> Structure {
    pdb::parse_str(text, Default::default())
        .unwrap()
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
        .unwrap()
        .into_structure()
}

fn bonds(topology: &Topology) -> BTreeSet<(usize, usize, String)> {
    topology
        .bonds()
        .map(|bond| {
            let [left, right] = bond.atoms().map(|atom| atom.index().index());
            (
                left.min(right),
                left.max(right),
                format!("{:?}", bond.bond().order),
            )
        })
        .collect()
}

fn issue(error: &kekule::diagnostics::FormatError) -> &PdbWriteErrorKind {
    match error.detail() {
        Some(FormatErrorDetail::PdbWrite(issue)) => issue,
        other => panic!("not a PDB write issue: {other:?}"),
    }
}

#[test]
fn deposited_entries_round_trip_through_the_writer() {
    for (name, text) in [("1CRN", CRAMBIN), ("1AKE", ADENYLATE_KINASE)] {
        let structure = interpret(text);
        let written = pdb::write(&structure, Default::default())
            .map(|(text, _)| text)
            .unwrap();
        let reread = interpret(&written);
        let (left, right) = (reread.topology(), structure.topology());
        assert_eq!(left.atom_count(), right.atom_count(), "{name}");
        assert_eq!(bonds(left), bonds(right), "{name}");
        assert_eq!(left.instance_count(), right.instance_count(), "{name}");
        for (a, b) in left.residues().zip(right.residues()) {
            assert_eq!(
                (a.name(), a.author_seq_id(), a.insertion_code()),
                (b.name(), b.author_seq_id(), b.insertion_code()),
                "{name}"
            );
            assert_eq!(a.chain().author_id(), b.chain().author_id(), "{name}");
        }
        assert_eq!(reread.cell(), structure.cell(), "{name}");
        // PDB -> structure -> PDB is a fixed point.
        assert_eq!(
            pdb::write(&reread, Default::default())
                .map(|(text, _)| text)
                .unwrap(),
            written,
            "{name}"
        );
        // Disulfides and ligand bonds are CONECT records; polymer bonds are not.
        assert!(
            written.lines().any(|line| line.starts_with("CONECT")),
            "{name}"
        );
        assert!(
            written.lines().any(|line| line.starts_with("TER")),
            "{name}"
        );
    }
    // 1CRN's P 1 21 1 cell is not the structure's cell, so no CRYST1 is
    // written.
    let crambin = pdb::write(&interpret(CRAMBIN), Default::default())
        .map(|(text, _)| text)
        .unwrap();
    assert!(interpret(CRAMBIN).cell().is_none());
    assert!(!crambin.contains("CRYST1"));
    // Three disulfides, each listed from both atoms.
    assert_eq!(
        crambin
            .lines()
            .filter(|line| line.starts_with("CONECT"))
            .count(),
        6
    );
}

/// An entry's text with its chains named by segment, as CHARMM and NAMD
/// name them: chain `A` blank in the chain column and `PROA` in columns
/// 73-76.
fn named_by_segment(text: &str) -> String {
    let mut segmented = String::new();
    for line in text.lines() {
        let mut bytes = line.as_bytes().to_vec();
        if line.starts_with("ATOM") || line.starts_with("HETATM") {
            let chain = bytes[21];
            bytes[21] = b' ';
            bytes[72..76].copy_from_slice(&[b'P', b'R', b'O', chain]);
        } else if line.starts_with("TER") {
            bytes[21] = b' ';
        }
        segmented.push_str(std::str::from_utf8(&bytes).unwrap());
        segmented.push('\n');
    }
    segmented
}

#[test]
fn chains_named_by_segment_round_trip_under_the_segment_policy() {
    let structure = interpret(ADENYLATE_KINASE);
    let (named, _) = pdb::write(&structure, Default::default()).unwrap();
    let segmented = named_by_segment(&named);
    let reread = interpret(&segmented);
    let chains = reread
        .topology()
        .chains()
        .map(|chain| chain.author_id().unwrap_or(chain.label_id()).to_owned())
        .collect::<BTreeSet<_>>();
    assert_eq!(chains, ["PROA", "PROB"].map(str::to_owned).into());

    // The chains fail by default, and are written back as segments.
    let error = pdb::write(&reread, Default::default()).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::ChainIdTooLong { chain_id, .. } if chain_id == "PROA"
    ));
    assert!(error.message().contains("does not fit PDB's chain columns"));
    let options = PdbWriteOptions::default().with_chain_ids(PdbChainIdPolicy::SegmentId);
    let (written, report) = pdb::write(&reread, options).unwrap();
    assert_eq!(written, segmented);
    assert_eq!(report.omissions(), [WriteOmission::ChainIds { chains: 2 }]);
    // Each chain's polymer still ends in its own TER record.
    let ters = |text: &str| text.lines().filter(|line| line.starts_with("TER")).count();
    assert_eq!(ters(&written), 2);
    assert_eq!(ters(&written), ters(&named));
    // One-character chains are written as before under the policy.
    let (written, report) = pdb::write(&structure, options).unwrap();
    assert_eq!(written, named);
    assert!(report.is_empty());
}

/// Methane solvated in a 1.5 nm padded cube of water at 0.15 M salt.
fn solvated() -> Structure {
    let mut structure = Structure::from_molecule(
        smiles::to_molecules("C").unwrap().remove(0),
        &Positions::zeros(1),
    )
    .unwrap();
    solvation::add_periodic_box(
        &mut structure,
        &PeriodicBoxOptions::Padding {
            padding: Quantity::new(1.5, NANOMETER),
            shape: BoxShape::Cube,
        },
    )
    .unwrap();
    solvation::solvate(
        &mut structure,
        &SolventOptions::default().with_ionic_strength(Quantity::new(0.15, MOLAR)),
    )
    .unwrap();
    structure
}

#[test]
fn solvated_systems_export_their_solvent_chains_through_every_writer() {
    let structure = solvated();
    let solvent = structure
        .topology()
        .chains()
        .map(|chain| chain.label_id().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(solvent, ["SOL1", "ION1", "ION2"]);
    // Solvent counts no implicit hydrogen, as PDB and mmCIF readers infer.
    for instance in structure.topology().instances().skip(1) {
        let molecule = instance.definition().molecule();
        for (atom, _) in molecule.atoms() {
            assert_eq!(molecule.implicit_hydrogens(atom), Some(0));
        }
    }
    // By default PDB refuses the chain identifiers, and nothing else.
    let error = pdb::write(&structure, Default::default()).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::ChainIdTooLong { chain_id, .. } if chain_id == "SOL1"
    ));

    // PDB writes the solvent chains as segments, which read back as chains.
    let options = PdbWriteOptions::default().with_chain_ids(PdbChainIdPolicy::SegmentId);
    let (text, report) = pdb::write(&structure, options).unwrap();
    assert!(report.omissions().contains(&WriteOmission::ChainIds {
        chains: solvent.len()
    }));
    for line in text.lines().filter(|line| line.starts_with("HETATM")) {
        if &line[17..20] == "HOH" {
            assert_eq!((&line[21..22], &line[72..76]), (" ", "SOL1"), "{line}");
        }
    }
    let reread = interpret(&text);
    assert_eq!(
        reread.topology().atom_count(),
        structure.topology().atom_count()
    );
    let chains = reread
        .topology()
        .chains()
        .map(|chain| chain.author_id().unwrap_or(chain.label_id()).to_owned())
        .collect::<BTreeSet<_>>();
    assert!(
        solvent.iter().all(|chain| chains.contains(chain)),
        "{chains:?}"
    );
    assert_eq!(pdb::write(&reread, options).unwrap().0, text);

    // mmCIF names the chains as they are.
    let (cif, _) = mmcif::write(
        [mmcif::MmcifBlockSource::structure(&structure)],
        Default::default(),
    )
    .unwrap();
    let reread = mmcif::parse_str(&cif, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_structure();
    assert_eq!(
        reread.topology().atom_count(),
        structure.topology().atom_count()
    );
    let chains = reread
        .topology()
        .chains()
        .map(|chain| chain.author_id().unwrap_or(chain.label_id()).to_owned())
        .collect::<BTreeSet<_>>();
    assert!(
        solvent.iter().all(|chain| chains.contains(chain)),
        "{chains:?}"
    );
}

#[test]
fn chain_ids_longer_than_a_segment_fail_under_the_segment_policy() {
    let mut editor = solvated().to_editor();
    let (solvent, _) = editor.topology().chains().next().unwrap();
    editor
        .set_chain_identifiers(solvent, "SOLVENT", None)
        .unwrap();
    let renamed = editor.finish().unwrap();
    let options = PdbWriteOptions::default().with_chain_ids(PdbChainIdPolicy::SegmentId);
    let error = pdb::write(&renamed, options).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::ChainIdTooLong { chain_id, .. } if chain_id == "SOLVENT"
    ));
}

#[test]
fn bonds_without_a_template_are_conect_listings_by_order() {
    // Bracket atoms fix their implicit hydrogen counts, which PDB cannot
    // carry; here the counts are those a reader infers.
    let topology = smiles::to_topology("CC(=O)[O-].[Na+]").unwrap();
    let points = (0..5)
        .map(|index| Point3::new(f64::from(index) * 1.4, 0.0, 0.0))
        .collect::<Vec<_>>();
    let structure = Structure::new(
        topology,
        Positions::from_vec(Quantity::new(points, ANGSTROM)).unwrap(),
    )
    .unwrap();
    let omit = PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let (written, report) = pdb::write(&structure, omit).unwrap();
    // [O-] and [Na+] fix their counts.
    assert_eq!(
        report.omissions(),
        [WriteOmission::ImplicitHydrogens { atoms: 2 }]
    );
    let lines = written.lines().collect::<Vec<_>>();
    assert_eq!(
        lines[0],
        "HETATM    1  C1  MOL     1       0.000   0.000   0.000                       C"
    );
    assert_eq!(
        lines[3],
        "HETATM    4  O4  MOL     1       4.200   0.000   0.000                       O1-"
    );
    assert_eq!(
        lines[4],
        "HETATM    5 NA1  MOL     2       5.600   0.000   0.000                      NA1+"
    );
    // The double bond is listed twice from each atom.
    assert!(lines.contains(&"CONECT    2    1    3    3    4"));
    assert!(lines.contains(&"CONECT    3    2    2"));
    let reread = interpret(&written);
    assert_eq!(bonds(reread.topology()), bonds(structure.topology()));
}

#[test]
fn large_residue_numbers_use_hybrid_36() {
    let text = "\
HETATM    1  O   HOH A9999       0.000   0.000   0.000  1.00  0.00           O
HETATM    2  O   HOH AA000       3.000   0.000   0.000  1.00  0.00           O
END
";
    let structure = interpret(text);
    let written = pdb::write(&structure, Default::default())
        .map(|(text, _)| text)
        .unwrap();
    assert_eq!(written, text);
    let numbers = structure
        .topology()
        .residues()
        .map(|residue| residue.author_seq_id().unwrap().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(numbers, ["9999", "10000"]);
}

#[test]
fn identifiers_pdb_cannot_hold_fail_explicitly() {
    let cif = |chain: &str, component: &str| {
        format!(
            "data_x
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
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
HETATM 1 O O {component} {chain} 1 . 0.000 0.000 0.000 1 {chain} 1
"
        )
    };
    let structure = |text: &str| {
        mmcif::parse_str(text, Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap()
            .into_structure()
    };
    let error = pdb::write(&structure(&cif("AB", "HOH")), Default::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Unsupported);
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::ChainIdTooLong { chain_id, .. } if chain_id == "AB"
    ));
    let error = pdb::write(&structure(&cif("A", "LONG")), Default::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::ResidueNameTooLong { name, .. } if name == "LONG"
    ));
    // Omitting chemistry does not relax identifiers.
    let options = PdbWriteOptions::default().with_chemistry(UnrepresentablePolicy::Omit);
    assert!(pdb::write(&structure(&cif("AB", "HOH")), options)
        .map(|(text, _)| text)
        .is_err());
}

#[test]
fn chemistry_pdb_cannot_carry_fails_unless_omitted() {
    // A bracket atom fixes its implicit hydrogen count.
    let topology = smiles::to_topology("[CH4]").unwrap();
    let structure = Structure::new(topology, Positions::zeros(1)).unwrap();
    let error = pdb::write(&structure, Default::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::UnsupportedAtomField {
            field: AtomField::Hydrogens,
            ..
        }
    ));
    // The chemistry policy no longer covers implicit hydrogens.
    let chemistry = PdbWriteOptions::default().with_chemistry(UnrepresentablePolicy::Omit);
    assert!(pdb::write(&structure, chemistry).is_err());
    let hydrogens = PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let (_, report) = pdb::write(&structure, hydrogens).unwrap();
    assert_eq!(
        report.omissions(),
        [WriteOmission::ImplicitHydrogens { atoms: 1 }]
    );

    let topology = smiles::to_topology("C$C").unwrap();
    let structure = Structure::new(topology, Positions::zeros(2)).unwrap();
    let error = pdb::write(&structure, Default::default())
        .map(|(text, _)| text)
        .unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::UnsupportedBondOrder {
            order: BondOrder::Quadruple,
            ..
        }
    ));
    // Omitting implicit hydrogens does not omit other chemistry.
    assert!(pdb::write(&structure, hydrogens).is_err());
    let (written, report) = pdb::write(&structure, chemistry).unwrap();
    assert!(!written.contains("CONECT"));
    assert_eq!(report.omissions(), [WriteOmission::Bonds { bonds: 1 }]);
}

#[test]
fn residues_must_be_consecutive_in_dense_order() {
    let water = smiles::to_molecules("[H]O[H]").unwrap().remove(0);
    let mut builder = StructureBuilder::new();
    let positions = Positions::zeros(3);
    let first = builder.add_molecule(water.clone(), &positions).unwrap();
    let second = builder.add_molecule(water, &positions).unwrap();
    let atom = |instance, index| InstanceAtomId::new(instance, AtomId::new(index));
    builder
        .set_atom_order([
            atom(first, 0),
            atom(second, 0),
            atom(first, 1),
            atom(first, 2),
            atom(second, 1),
            atom(second, 2),
        ])
        .unwrap();
    let structure = builder.build().unwrap();
    let omit = PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let error = pdb::write(&structure, omit)
        .map(|(text, _)| text)
        .unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::NonContiguousResidue { .. }
    ));
}

/// Geometries of `positions` on `topology`, unordered.
fn geometries(
    topology: Arc<Topology>,
    positions: impl IntoIterator<Item = Positions>,
) -> StructureSet {
    StructureSet::from_geometries(topology, positions.into_iter().map(Geometry::new)).unwrap()
}

#[test]
fn structure_sets_write_one_model_per_geometry() {
    let topology = smiles::to_topology("C").unwrap();
    let positions = |x: f64| {
        Positions::from_vec(Quantity::new(vec![Point3::new(x, 0.0, 0.0)], ANGSTROM)).unwrap()
    };
    let set = geometries(Arc::clone(&topology), [positions(0.0), positions(1.0)]);
    let written = pdb::write(&set, Default::default())
        .map(|(text, _)| text)
        .unwrap();
    assert_eq!(written.matches("MODEL ").count(), 2);
    assert_eq!(written.matches("ENDMDL").count(), 2);
    assert_eq!(
        written,
        "MODEL        1\n\
         HETATM    1  C1  MOL     1       0.000   0.000   0.000                       C\n\
         ENDMDL\n\
         MODEL        2\n\
         HETATM    1  C1  MOL     1       1.000   0.000   0.000                       C\n\
         ENDMDL\n\
         END\n"
    );
    let reread = pdb::parse_str(&written, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap();
    assert_eq!(reread.structure_set().len(), 2);
    assert_eq!(reread.structure_set().order(), None);
    for (geometry, source) in reread.structure_set().iter().zip(&set) {
        assert_eq!(geometry.positions(), source.positions());
    }
    assert_eq!(
        pdb::write(reread.structure_set(), Default::default())
            .map(|(text, _)| text)
            .unwrap(),
        written
    );
}

#[test]
fn ordered_structure_sets_write_models_and_their_places_fail_unless_omitted() {
    let topology = smiles::to_topology("C").unwrap();
    let positions = |x: f64| {
        Positions::from_vec(Quantity::new(vec![Point3::new(x, 0.0, 0.0)], ANGSTROM)).unwrap()
    };
    let mut scan = StructureSet::from_geometries(
        Arc::clone(&topology),
        (0..3).map(|index| Geometry::new(positions(f64::from(index)))),
    )
    .unwrap();
    scan.set_order(Some(Order::Coordinate {
        unit: DEGREE,
        values: vec![0.0, 120.0, 240.0],
    }))
    .unwrap();
    let error = pdb::write(&scan, Default::default()).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::UnrepresentableGeometry {
            geometry: 1,
            field: GeometryField::OrderCoordinate
        }
    ));
    assert_eq!(
        error.to_string(),
        "Cannot write the PDB output. The order coordinate of geometry 1 has no PDB representation."
    );
    let options = PdbWriteOptions::default().with_geometry_metadata(UnrepresentablePolicy::Omit);
    let (written, report) = pdb::write(&scan, options).unwrap();
    assert_eq!(
        report.omissions(),
        [WriteOmission::GeometryState {
            field: GeometryField::OrderCoordinate,
            geometries: 3
        }]
    );
    assert_eq!(written.matches("MODEL ").count(), 3);
    assert_eq!(written.matches("ENDMDL").count(), 3);
    assert_eq!(
        written
            .lines()
            .filter(|line| line.starts_with("HETATM"))
            .collect::<Vec<_>>(),
        [0.0, 1.0, 2.0].map(|x| format!(
            "HETATM    1  C1  MOL     1    {x:>8.3}   0.000   0.000                       C"
        ))
    );
    // The models read back as an unordered set of the geometries' positions.
    let reread = pdb::parse_str(&written, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap();
    assert_eq!(reread.reports().len(), 3);
    let read = reread.structure_set();
    assert_eq!(read.len(), 3);
    assert_eq!(read.order(), None);
    for (geometry, source) in read.iter().zip(&scan) {
        assert_eq!(geometry.positions(), source.positions());
    }

    // A step order fails on its steps, and a lone geometry on its own.
    let mut history =
        StructureSet::from_geometries(topology, [Geometry::new(positions(0.0))]).unwrap();
    history.set_order(Some(Order::Step(vec![3]))).unwrap();
    let error = pdb::write(&history, Default::default()).unwrap_err();
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::UnrepresentableGeometry {
            geometry: 1,
            field: GeometryField::OrderStep
        }
    ));
    let lone = history.get(0).unwrap();
    let error = pdb::write(pdb::PdbSource::structure(&lone), Default::default()).unwrap_err();
    assert_eq!(
        error.detail(),
        Some(&FormatErrorDetail::UnrepresentableState(
            GeometryField::OrderStep
        ))
    );
    assert_eq!(
        error.to_string(),
        "Cannot write the PDB output. The order step of the geometry has no PDB representation."
    );
    let (_, report) = pdb::write(pdb::PdbSource::structure(&lone), options).unwrap();
    assert_eq!(
        report.omissions(),
        [WriteOmission::GeometryState {
            field: GeometryField::OrderStep,
            geometries: 1
        }]
    );

    // Plural state is named with a plural verb.
    let mut moving =
        StructureSet::from_geometries(history.shared_topology(), [Geometry::new(positions(0.0))])
            .unwrap();
    moving
        .get_mut(0)
        .unwrap()
        .set_velocities(Some(Velocities::zeros(1)))
        .unwrap();
    let error = pdb::write(&moving, Default::default()).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Cannot write the PDB output. The velocities of geometry 1 have no PDB representation."
    );
    let lone = moving.get(0).unwrap();
    let error = pdb::write(pdb::PdbSource::structure(&lone), Default::default()).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Cannot write the PDB output. The velocities of the geometry have no PDB representation."
    );

    // An empty structure set is no entry.
    let empty = StructureSet::new(smiles::to_topology("C").unwrap());
    let error = pdb::write(&empty, options).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert!(matches!(
        issue(&error),
        PdbWriteErrorKind::EmptyStructureSet
    ));
}

#[test]
fn cells_crystal_records_cannot_carry_fail_unless_omitted() {
    let topology = smiles::to_topology("C").unwrap();
    let cell = |lengths: Vector3, periodic: [bool; 3]| {
        PeriodicCell::orthorhombic(Quantity::new(lengths, ANGSTROM), periodic).unwrap()
    };
    let structure = |cell| {
        Structure::new(Arc::clone(&topology), Positions::zeros(1))
            .unwrap()
            .with_cell(cell)
    };
    let written = pdb::write(
        &structure(cell(Vector3::new(10.0, 20.0, 30.0), [true; 3])),
        Default::default(),
    )
    .map(|(text, _)| text)
    .unwrap();
    assert!(written
        .starts_with("CRYST1   10.000   20.000   30.000  90.00  90.00  90.00 P 1           1\n"));
    for (lengths, periodic, reason) in [
        (
            Vector3::new(10.0, 20.0, 30.0),
            [true, true, false],
            CellReason::PartiallyPeriodic,
        ),
        (
            Vector3::new(1.0, 1.0, 1.0),
            [true; 3],
            CellReason::Placeholder,
        ),
    ] {
        let source = structure(cell(lengths, periodic));
        let error = pdb::write(&source, Default::default())
            .map(|(text, _)| text)
            .unwrap_err();
        assert_eq!(
            issue(&error),
            &PdbWriteErrorKind::UnrepresentableCell { reason },
            "{reason:?}"
        );
        let options = PdbWriteOptions::default().with_cell(UnrepresentablePolicy::Omit);
        let (text, report) = pdb::write(&source, options).unwrap();
        assert!(!text.contains("CRYST1"));
        assert_eq!(report.omissions(), [WriteOmission::Cell { structures: 1 }]);
    }
}

#[test]
fn sink_failures_are_io_errors_after_validation() {
    struct Failing;
    impl std::io::Write for Failing {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("full"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let structure = interpret(CRAMBIN);
    let error = pdb::write_to(&mut Failing, &structure, Default::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Io);
    let mut sink = Vec::new();
    pdb::write_to(&mut sink, &structure, Default::default()).unwrap();
    assert_eq!(
        String::from_utf8(sink).unwrap(),
        pdb::write(&structure, Default::default())
            .map(|(text, _)| text)
            .unwrap()
    );
}

/// One partner of a coordination record: its author chain, residue number,
/// component, and atom name.
type Partner = (String, String, String, String);

/// The coordination a report states, each pair of partners once and in
/// sorted order.
fn coordination(report: &StructureInterpretationReport) -> BTreeSet<[Partner; 2]> {
    let partner = |atom| {
        let source = report
            .instances()
            .iter()
            .flat_map(|instance| instance.atoms())
            .find(|candidate| candidate.atom() == atom)
            .unwrap();
        (
            source.auth_asym_id().unwrap_or(source.asym_id()).to_owned(),
            source.author_sequence_id().unwrap_or_default().to_owned(),
            source
                .auth_component_id()
                .unwrap_or(source.component_id())
                .to_owned(),
            source
                .auth_atom_name()
                .unwrap_or(source.atom_name())
                .to_owned(),
        )
    };
    report
        .coordination()
        .iter()
        .map(|record| {
            let mut pair = record.atoms().map(partner);
            pair.sort();
            pair
        })
        .collect()
}

fn interpretation(text: &str) -> StructureInterpretation {
    pdb::parse_str(text, Default::default())
        .unwrap()
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
        .unwrap()
}

fn records<'a>(text: &'a str, name: &str) -> Vec<&'a str> {
    text.lines()
        .filter(|line| line.split_whitespace().next() == Some(name))
        .collect()
}

#[test]
fn metal_coordination_round_trips_as_link_records() {
    // Rubredoxin's iron is bound to four cysteines by LINK records.
    let source = interpretation(RUBREDOXIN);
    let iron = coordination(source.report());
    assert_eq!(iron.len(), 4);
    // The reader records each pair once, though LINK and CONECT both list it.
    assert_eq!(source.report().coordination().len(), 4);
    assert!(iron.iter().all(|pair| {
        pair.iter().any(|partner| partner.3 == "SG") && pair.iter().any(|partner| partner.2 == "FE")
    }));

    let (written, report) = pdb::write(&source, Default::default()).unwrap();
    assert!(report.omissions().is_empty());
    // The written records are the deposited ones, in the format's columns.
    let deposited = records(RUBREDOXIN, "LINK")
        .into_iter()
        .map(str::trim_end)
        .collect::<Vec<_>>();
    assert_eq!(records(&written, "LINK"), deposited);
    // They come before the coordinates.
    assert!(written.starts_with("LINK"));
    let reread = interpretation(&written);
    assert_eq!(coordination(reread.report()), iron);
    // PDB -> interpretation -> PDB is a fixed point.
    assert_eq!(pdb::write(&reread, Default::default()).unwrap().0, written);

    // Without reports the entry states no coordination, as before.
    let (bare, report) = pdb::write(source.structure(), Default::default()).unwrap();
    assert!(report.omissions().is_empty());
    assert!(records(&bare, "LINK").is_empty());
    assert!(coordination(interpretation(&bare).report()).is_empty());
    assert_eq!(
        written
            .lines()
            .filter(|line| !line.starts_with("LINK"))
            .collect::<Vec<_>>(),
        bare.lines().collect::<Vec<_>>()
    );
}

#[test]
fn mmcif_metal_coordination_round_trips_through_pdb() {
    // KRas G12C: Ser17 OG and GDP O2B coordinate Mg2+ (`metalc`).
    let source = mmcif::parse_str(KRAS, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let stated = coordination(source.report());
    assert_eq!(stated.len(), 2);
    // PDB's three residue-name columns cannot hold the inhibitor's code,
    // A1L9E, so it is written as the unknown ligand UNL.
    let inhibitor = source
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .find(|atom| atom.component_id() == "A1L9E")
        .unwrap()
        .atom();
    let mut editor = source.structure().to_editor();
    let handle = editor.topology().atom_handle(inhibitor).unwrap();
    let site = editor
        .topology()
        .atom_site_for_atom(handle)
        .unwrap()
        .unwrap();
    let residue = editor.topology().atom_site(site).unwrap().residue();
    let label = editor
        .topology()
        .residue(residue)
        .unwrap()
        .label_comp_id()
        .map(str::to_owned);
    editor
        .set_residue_component_ids(residue, label, Some("UNL".to_owned()))
        .unwrap();
    let structure = editor.finish().unwrap();
    let source =
        PdbSource::structure(&structure).with_reports(std::slice::from_ref(source.report()));
    let (written, report) = pdb::write(source, Default::default()).unwrap();
    assert!(report.omissions().is_empty(), "{:?}", report.omissions());
    assert_eq!(records(&written, "LINK").len(), 2);
    assert_eq!(coordination(interpretation(&written).report()), stated);
}

#[test]
fn coordination_the_entry_cannot_state_is_reported() {
    let source = interpretation(RUBREDOXIN);
    let reports = std::slice::from_ref(source.report());
    let iron = source
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .find(|atom| atom.atom_name() == "FE")
        .unwrap()
        .atom();
    let cysteine_6 = source
        .report()
        .instances()
        .iter()
        .flat_map(|instance| instance.atoms())
        .find(|atom| atom.atom_name() == "SG" && atom.author_sequence_id() == Some("6"))
        .unwrap()
        .atom();

    // A partner the structure no longer has.
    let mut editor = source.structure().to_editor();
    let handle = editor.topology().atom_handle(iron).unwrap();
    editor.delete_atom(handle).unwrap();
    let without_iron = editor.finish().unwrap();
    let (written, report) = pdb::write(
        PdbSource::structure(&without_iron).with_reports(reports),
        Default::default(),
    )
    .unwrap();
    assert!(records(&written, "LINK").is_empty());
    assert_eq!(
        report.omissions(),
        [WriteOmission::Coordination { records: 4 }]
    );

    // A partner whose atom site the report does not name so.
    let mut editor = source.structure().to_editor();
    let site = editor
        .topology()
        .atom_site_for_atom(editor.topology().atom_handle(cysteine_6).unwrap())
        .unwrap()
        .unwrap();
    let metadata = editor
        .topology()
        .atom_site(site)
        .unwrap()
        .metadata()
        .clone()
        .with_auth_atom_id(Some("SX".to_owned()));
    editor.set_atom_site_metadata(site, metadata).unwrap();
    let renamed = editor.finish().unwrap();
    let (written, report) = pdb::write(
        PdbSource::structure(&renamed).with_reports(reports),
        Default::default(),
    )
    .unwrap();
    assert_eq!(records(&written, "LINK").len(), 3);
    assert!(!written.contains("LINK         SG  CYS A   6"));
    assert_eq!(
        report.omissions(),
        [WriteOmission::Coordination { records: 1 }]
    );

    // A chain written as a segment identifier, which LINK cannot name.
    let mut editor = source.structure().to_editor();
    let sites = editor
        .topology()
        .atom_sites()
        .map(|(id, site)| (id, site.metadata().clone()))
        .collect::<Vec<_>>();
    for (id, metadata) in sites {
        editor
            .set_atom_site_metadata(id, metadata.with_auth_asym_id(Some("RUBR".to_owned())))
            .unwrap();
    }
    let segmented = editor.finish().unwrap();
    let options = PdbWriteOptions::default().with_chain_ids(PdbChainIdPolicy::SegmentId);
    let (written, report) = pdb::write(
        PdbSource::structure(&segmented).with_reports(reports),
        options,
    )
    .unwrap();
    assert!(records(&written, "LINK").is_empty());
    assert_eq!(
        report.omissions(),
        [
            WriteOmission::ChainIds { chains: 1 },
            WriteOmission::Coordination { records: 4 }
        ]
    );

    // Coordination between two metal ions, which a LINK would make a bond.
    let ions = mmcif::parse_str(
        "data_ions
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
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
HETATM 1 ZN ZN ZN A 1 . 0.000 0.000 0.000 1 A 1
HETATM 2 ZN ZN ZN B 1 . 2.600 0.000 0.000 2 A 1
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_atom_id
metalc1 metalc A ZN ZN B ZN ZN
",
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap();
    assert_eq!(ions.report().coordination().len(), 1);
    let (written, report) = pdb::write(&ions, Default::default()).unwrap();
    assert!(records(&written, "LINK").is_empty());
    assert_eq!(
        report.omissions(),
        [WriteOmission::Coordination { records: 1 }]
    );
    assert_eq!(interpretation(&written).topology().instance_count(), 2);

    // One report per structure.
    let error = pdb::write(
        PdbSource::structure(source.structure()).with_reports(&[]),
        Default::default(),
    )
    .unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        issue(&error),
        &PdbWriteErrorKind::ReportCountMismatch(LengthMismatch::new(Countable::Report, 1, 0))
    );
    assert_eq!(
        error.message(),
        "There are 0 interpretation reports instead of 1."
    );
}

#[test]
fn coordination_link_records_would_not_select_is_reported() {
    // CONECT names atoms by serial number, so it can state coordination to
    // water 10 beside water 10A. A LINK with a blank insertion code names
    // both, so readers would resolve neither.
    let record = |serial: u32, name: &str, residue: &str, sequence: &str, x: f64, element: &str| {
        format!(
            "HETATM{serial:>5} {name:<4} {residue:>3} A{sequence:>5}   {x:>8.3}{:>8.3}{:>8.3}  1.00  0.00          {element:>2}\n",
            0.0, 0.0
        )
    };
    let text = [
        record(1, "ZN", "ZN", "1 ", 0.0, "ZN"),
        record(2, " O", "HOH", "10 ", 2.1, "O"),
        record(3, " O", "HOH", "10A", 9.0, "O"),
        "CONECT    1    2\nCONECT    2    1\nEND\n".to_owned(),
    ]
    .concat();
    let source = interpretation(&text);
    assert_eq!(source.report().coordination().len(), 1);
    let (written, report) = pdb::write(&source, Default::default()).unwrap();
    assert!(records(&written, "LINK").is_empty());
    assert_eq!(
        report.omissions(),
        [WriteOmission::Coordination { records: 1 }]
    );

    // Without water 10A the record selects its water alone.
    let text = text
        .lines()
        .filter(|line| !line.contains("10A"))
        .map(|line| format!("{line}\n"))
        .collect::<String>();
    let source = interpretation(&text);
    let (written, report) = pdb::write(&source, Default::default()).unwrap();
    assert_eq!(records(&written, "LINK").len(), 1);
    assert!(report.omissions().is_empty());
    assert_eq!(
        coordination(interpretation(&written).report()),
        coordination(source.report())
    );
}

#[test]
fn structure_sets_state_the_coordination_every_geometry_states() {
    // The same entry without its first LINK record and its CONECT records,
    // which list the coordination again, states three of the four
    // iron-sulfur records.
    let full = interpretation(RUBREDOXIN);
    let first_link = RUBREDOXIN
        .lines()
        .position(|line| line.starts_with("LINK"))
        .unwrap();
    let partial = interpretation(
        &RUBREDOXIN
            .lines()
            .enumerate()
            .filter(|&(index, line)| index != first_link && !line.starts_with("CONECT"))
            .map(|(_, line)| format!("{line}\n"))
            .collect::<String>(),
    );
    assert_eq!(coordination(partial.report()).len(), 3);
    assert_eq!(
        partial.structure().atom_count(),
        full.structure().atom_count()
    );

    let set = geometries(
        full.structure().shared_topology(),
        [full.structure(), partial.structure()]
            .map(|structure| structure.as_structure_view().positions().clone()),
    );
    let options = PdbWriteOptions::default();

    // LINK records belong to the entry, so geometries that agree keep them all.
    let agreeing = [full.report().clone(), full.report().clone()];
    let (written, report) = pdb::write(
        PdbSource::structure_set(&set).with_reports(&agreeing),
        options,
    )
    .unwrap();
    assert_eq!(records(&written, "LINK").len(), 4);
    // They are the entry's own four records.
    let source_links = RUBREDOXIN
        .lines()
        .filter(|line| line.starts_with("LINK"))
        .map(str::trim_end)
        .collect::<Vec<_>>();
    assert_eq!(records(&written, "LINK"), source_links);
    assert_eq!(records(&written, "MODEL").len(), 2);
    assert!(written.starts_with("LINK"));
    assert!(!report
        .omissions()
        .iter()
        .any(|omission| matches!(omission, WriteOmission::Coordination { .. })));

    // A record one geometry does not state is left out of the entry.
    let differing = [full.report().clone(), partial.report().clone()];
    let (written, report) = pdb::write(
        PdbSource::structure_set(&set).with_reports(&differing),
        options,
    )
    .unwrap();
    assert_eq!(records(&written, "LINK").len(), 3);
    assert_eq!(records(&written, "LINK"), source_links[1..]);
    let differing_written = written.clone();
    assert!(report
        .omissions()
        .contains(&WriteOmission::Coordination { records: 1 }));
    let reread = pdb::parse_str(&written, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap();
    assert_eq!(reread.structure_set().len(), 2);
    for (geometry, source) in reread.structure_set().iter().zip(&set) {
        assert_eq!(geometry.positions(), source.positions());
    }
    for report in reread.reports() {
        assert_eq!(coordination(report), coordination(partial.report()));
    }
    // The order of the reports does not matter.
    let reversed = [partial.report().clone(), full.report().clone()];
    let written = pdb::write(
        PdbSource::structure_set(&set).with_reports(&reversed),
        options,
    )
    .unwrap()
    .0;
    assert_eq!(records(&written, "LINK").len(), 3);
    assert_eq!(written, differing_written);

    let error = pdb::write(
        PdbSource::structure_set(&set).with_reports(std::slice::from_ref(full.report())),
        options,
    )
    .unwrap_err();
    assert_eq!(
        issue(&error),
        &PdbWriteErrorKind::ReportCountMismatch(LengthMismatch::new(Countable::Report, 2, 1))
    );
}
