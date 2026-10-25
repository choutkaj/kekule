//! Macromolecular structures survive mmCIF output with their analyses intact.

use kekule::dssp::{self, DsspOptions};
use kekule::macromolecular::{self, ModelSelection, StructureInterpretOptions};
use kekule::mmcif::{self, MmcifBlockSource, MmcifWriteOptions};
use kekule::structure::Model;

use crate::support::fixture;

/// Every bond by the dense indices of its atoms, with its order, sorted.
fn bonds(model: &Model) -> Vec<(usize, usize, String)> {
    let mut bonds = model
        .topology()
        .bonds()
        .map(|bond| {
            let [a, b] = bond.atoms().map(|atom| atom.index().index());
            (a.min(b), a.max(b), format!("{:?}", bond.bond().order))
        })
        .collect::<Vec<_>>();
    bonds.sort();
    bonds
}

fn interpret(text: &str) -> (Model, macromolecular::StructureInterpretationReport) {
    (mmcif::parse_str(text, Default::default()).unwrap())
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
        .unwrap()
        .into_parts()
}

/// Atoms, bonds, hierarchy labels, coordinates at the written precision, and
/// DSSP secondary structure are unchanged by writing and re-reading an entry,
/// including its links between residues and to metals. Bonds the reader
/// takes from templates may come back in another order within a molecule,
/// which mmCIF does not carry.
#[test]
fn mmcif_round_trip_preserves_topology_coordinates_and_secondary_structure() {
    for path in ["1CRN", "1AKE", "1IRO", "6OIM"].map(|entry| fixture(&format!("mmcif/{entry}.cif")))
    {
        let entry = path.file_name().unwrap().to_string_lossy().into_owned();
        let (original, report) = interpret(&std::fs::read_to_string(&path).unwrap());
        let written = mmcif::write(
            [MmcifBlockSource::model(&original).with_reports(std::slice::from_ref(&report))],
            MmcifWriteOptions::default(),
        )
        .map(|(text, _)| text)
        .unwrap_or_else(|error| panic!("{entry}: {error}"));
        let (reread, reread_report) = interpret(&written);
        // The written block is a fixed point: the model read from it writes
        // the same text and reads back with the same layout.
        let rewritten = mmcif::write(
            [MmcifBlockSource::model(&reread).with_reports(std::slice::from_ref(&reread_report))],
            MmcifWriteOptions::default(),
        )
        .unwrap()
        .0;
        assert_eq!(rewritten, written, "{entry}: rewritten");
        assert!(
            interpret(&rewritten)
                .0
                .topology()
                .same_layout(reread.topology()),
            "{entry}: layout"
        );

        let atoms = |model: &Model| {
            model
                .topology()
                .atom_ids()
                .iter()
                .map(|&id| {
                    let atom = model.topology().atom(id).unwrap();
                    (atom.element().symbol(), atom.formal_charge())
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(atoms(&reread), atoms(&original), "{entry}: atoms");
        assert_eq!(bonds(&reread), bonds(&original), "{entry}: bonds");
        let residues = |model: &Model| {
            model
                .topology()
                .residues()
                .map(|residue| format!("{residue:?}"))
                .collect::<Vec<_>>()
        };
        assert_eq!(residues(&reread), residues(&original), "{entry}: residues");

        let digits = i32::try_from(MmcifWriteOptions::default().coordinate_precision()).unwrap();
        // Half a unit in the last written place, plus rounding slack.
        let precision = 0.5 * 10f64.powi(-digits) + 1e-12;
        for (left, right) in original
            .positions()
            .values()
            .value()
            .iter()
            .zip(reread.positions().values().value().iter())
        {
            let delta = (left.x - right.x)
                .abs()
                .max((left.y - right.y).abs())
                .max((left.z - right.z).abs());
            assert!(delta <= precision, "{entry}: {left:?} -> {right:?}");
        }

        let secondary_structure = |model: &Model| {
            dssp::assign(model, DsspOptions::default())
                .unwrap()
                .residues()
                .map(|residue| residue.secondary_structure().code())
                .collect::<String>()
        };
        assert_eq!(
            secondary_structure(&reread),
            secondary_structure(&original),
            "{entry}: DSSP"
        );
    }
}

/// The `_struct_conn` rows of written text, one value list per row.
fn struct_conn_rows(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .skip_while(|line| *line != "_struct_conn.pdbx_value_order")
        .skip(1)
        .take_while(|line| *line != "#")
        .map(|line| line.split_whitespace().map(str::to_owned).collect())
        .collect()
}

/// Like the PDB's own files, a written block states only the bonds that
/// templates and polymer linkage do not imply: links between residues other
/// than the standard peptide and phosphodiester links, metal coordination,
/// and the bonds of residues without a bundled template.
#[test]
fn mmcif_states_only_the_bonds_templates_and_polymer_linkage_do_not_imply() {
    let link = |kind, left: [&str; 4], right: [&str; 4]| {
        [kind]
            .into_iter()
            .chain(left)
            .chain(right)
            .collect::<Vec<_>>()
            .join(" ")
    };
    let sg = |seq| ["A", "CYS", seq, "SG"];
    let magnesium = ["B", "MG", ".", "MG"];
    let water = ["E", "HOH", ".", "O"];
    for (entry, links, untemplated) in [
        (
            "1CRN",
            vec![
                link("covale", sg("3"), sg("40")),
                link("covale", sg("4"), sg("32")),
                link("covale", sg("16"), sg("26")),
            ],
            vec![],
        ),
        ("1AKE", vec![], vec!["AP5"]),
        (
            "1IRO",
            ["6", "9", "39", "42"]
                .map(|seq| link("metalc", sg(seq), ["B", "FE", ".", "FE"]))
                .to_vec(),
            // The N-terminal ammonium has a hydrogen its template does not
            // name.
            vec!["MET"],
        ),
        (
            "6OIM",
            vec![
                link("covale", sg("26"), ["D", "MOV", ".", "C25"]),
                link("metalc", ["A", "SER", "31", "OG"], magnesium),
                link("metalc", magnesium, ["C", "GDP", ".", "O2B"]),
                link("metalc", magnesium, water),
                link("metalc", magnesium, water),
                link("metalc", magnesium, water),
                link("metalc", magnesium, water),
            ],
            vec!["MOV"],
        ),
    ] {
        let (model, report) =
            interpret(&std::fs::read_to_string(fixture(&format!("mmcif/{entry}.cif"))).unwrap());
        let written = mmcif::write(
            [MmcifBlockSource::model(&model).with_reports(std::slice::from_ref(&report))],
            MmcifWriteOptions::default(),
        )
        .unwrap()
        .0;
        let rows = struct_conn_rows(&written);
        let (within, between): (Vec<_>, Vec<_>) = rows
            .iter()
            .partition(|row| (&row[2], &row[3], &row[4]) == (&row[6], &row[7], &row[8]));
        let between = between
            .iter()
            .map(|row| row[1..10].join(" "))
            .collect::<Vec<_>>();
        assert_eq!(between, links, "{entry}");
        let components = within
            .iter()
            .map(|row| row[3].as_str())
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(
            components.into_iter().collect::<Vec<_>>(),
            untemplated,
            "{entry}"
        );
    }
}

/// Copies of KRas G12C side by side write and read back with every bond.
/// The block states only the ligand's bonds and the links of each copy.
#[test]
fn tiled_structure_round_trips_through_mmcif() {
    tiled_round_trip(8);
}

/// 66 copies, about 100,000 atoms as in a large export, round-trip in
/// seconds in release.
#[test]
#[ignore = "a 100,000-atom round trip; run in release with --ignored"]
fn large_tiled_structure_round_trips_through_mmcif() {
    tiled_round_trip(66);
}

fn tiled_round_trip(copies: u32) {
    use kekule::geometry::{Matrix3, RigidTransform, Vector3};
    use kekule::units::{Quantity, ANGSTROM};

    // The entry read through PDB, as an application that opened it there
    // does: its atom sites keep author identity alone.
    let (entry, _) = interpret(&std::fs::read_to_string(fixture("mmcif/6OIM.cif")).unwrap());
    let kras = kekule::pdb::parse_str(
        &kekule::pdb::write(&entry, Default::default()).unwrap().0,
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap()
    .into_model();
    let mut editor = kras.to_editor();
    for copy in 1..copies {
        let offset = Vector3::new(
            f64::from(copy % 6) * 80.0,
            f64::from(copy / 6 % 6) * 80.0,
            f64::from(copy / 36) * 80.0,
        );
        let mut shifted = kras.clone();
        shifted
            .apply_rigid_transform(
                &RigidTransform::new(Matrix3::identity(), Quantity::new(offset, ANGSTROM)).unwrap(),
            )
            .unwrap();
        let appended = editor.append_model(&shifted).unwrap();
        // mmCIF names each chain once, so every copy's chains are renamed.
        for chain in kras.topology().chains() {
            editor
                .set_chain_identifiers(
                    appended.chain(chain.id()).unwrap(),
                    format!("{}{copy}", chain.label_id()),
                    chain.author_id().map(str::to_owned),
                )
                .unwrap();
        }
    }
    let tiled = editor.finish().unwrap();
    assert_eq!(tiled.atom_count(), copies as usize * kras.atom_count());

    let written = mmcif::write([&tiled], MmcifWriteOptions::default())
        .unwrap()
        .0;
    let single = mmcif::write([&kras], MmcifWriteOptions::default())
        .unwrap()
        .0;
    assert_eq!(
        struct_conn_rows(&written).len(),
        copies as usize * struct_conn_rows(&single).len()
    );
    let reread = mmcif::parse_str(&written, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_model();
    assert_eq!(reread.atom_count(), tiled.atom_count());
    assert_eq!(
        reread.topology().instance_count(),
        tiled.topology().instance_count()
    );
    assert_eq!(bonds(&reread), bonds(&tiled));
}

/// A large structure without templates, as XYZ gives, carries every bond as
/// a `_struct_conn` row and reads back with them.
#[test]
fn large_untemplated_structure_round_trips_through_mmcif() {
    use kekule::diagnostics::UnrepresentablePolicy;

    // 8,000 methanes, 40,000 atoms, on a 4 Å grid.
    let hydrogens = [
        [0.63, 0.63, 0.63],
        [-0.63, -0.63, 0.63],
        [-0.63, 0.63, -0.63],
        [0.63, -0.63, -0.63],
    ];
    let mut text = String::from("40000\nmethanes\n");
    for index in 0..8000 {
        let origin = [index % 20, index / 20 % 20, index / 400].map(|step| f64::from(step) * 4.0);
        text.push_str(&format!("C {} {} {}\n", origin[0], origin[1], origin[2]));
        for offset in hydrogens {
            text.push_str(&format!(
                "H {} {} {}\n",
                origin[0] + offset[0],
                origin[1] + offset[1],
                origin[2] + offset[2]
            ));
        }
    }
    let model = kekule::xyz::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_model();
    assert_eq!(model.topology().bond_count(), 32_000);
    let options = MmcifWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit);
    let written = mmcif::write([&model], options).unwrap().0;
    assert_eq!(struct_conn_rows(&written).len(), 32_000);
    let reread = mmcif::parse_str(&written, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_model();
    assert_eq!(bonds(&reread), bonds(&model));
}
