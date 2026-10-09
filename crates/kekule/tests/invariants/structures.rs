//! Macromolecular structures survive mmCIF output with their analyses intact.

use kekule::dssp::{self, DsspOptions};
use kekule::mmcif::{
    self, MmcifBlockSource, MmcifInterpretOptions, MmcifModelSelection, MmcifWriteOptions,
};
use kekule::structure::Model;

use crate::support::{corpus, fixture};

fn interpret(text: &str) -> (Model, mmcif::MmcifInterpretationReport) {
    mmcif::interpret(
        &mmcif::parse_str(text).unwrap(),
        MmcifInterpretOptions {
            model_selection: MmcifModelSelection::First,
            ..MmcifInterpretOptions::default()
        },
    )
    .unwrap()
    .into_parts()
}

/// Atoms, bonds, hierarchy labels, coordinates at the written precision, and
/// DSSP secondary structure are unchanged by writing and re-reading an entry.
#[test]
fn mmcif_round_trip_preserves_topology_coordinates_and_secondary_structure() {
    for path in [
        corpus("smoke/data/rcsb/1CRN.cif"),
        fixture("mmcif/1AKE.cif"),
    ] {
        let entry = path.file_name().unwrap().to_string_lossy().into_owned();
        let (original, report) = interpret(&std::fs::read_to_string(&path).unwrap());
        let written = mmcif::write(
            [MmcifBlockSource::model(&original).with_reports(std::slice::from_ref(&report))],
            MmcifWriteOptions::default(),
        )
        .unwrap_or_else(|error| panic!("{entry}: {error}"));
        let (reread, _) = interpret(&written);

        let atoms = |model: &Model| {
            model
                .topology()
                .atom_ids()
                .iter()
                .map(|&id| {
                    let atom = model.topology().atom(id).unwrap();
                    (atom.element.symbol(), atom.formal_charge)
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(atoms(&reread), atoms(&original), "{entry}: atoms");
        let bonds = |model: &Model| {
            let mut bonds = model
                .topology()
                .bond_ids()
                .iter()
                .map(|&id| {
                    let bond = model.topology().bond(id).unwrap();
                    let index = |atom| {
                        let atom = kekule::topology::InstanceAtomId::new(id.molecule(), atom);
                        model.topology().atom_index(atom).unwrap().index()
                    };
                    let (a, b) = (index(bond.a()), index(bond.b()));
                    (a.min(b), a.max(b), format!("{:?}", bond.order))
                })
                .collect::<Vec<_>>();
            bonds.sort();
            bonds
        };
        assert_eq!(bonds(&reread), bonds(&original), "{entry}: bonds");
        let residues = |model: &Model| {
            model
                .topology()
                .residues()
                .map(|residue| format!("{residue:?}"))
                .collect::<Vec<_>>()
        };
        assert_eq!(residues(&reread), residues(&original), "{entry}: residues");

        let digits = i32::try_from(MmcifWriteOptions::default().coordinate_precision).unwrap();
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
            dssp::assign(model.as_model_view(), DsspOptions::default())
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
