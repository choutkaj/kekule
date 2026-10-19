//! PDB writer acceptance: every curated benchmark entry that Kekule
//! interprets and PDB can represent round-trips through the PDB writer.
//!
//! Each mmCIF entry is interpreted, written as PDB, and read back. The
//! reread model must keep every atom's element, charge, identity columns,
//! coordinates, occupancy, and B factor at PDB's printed precision, and
//! every bond with its order; writing it again must give the same text, so
//! PDB -> model -> PDB is a fixed point. NMR ensembles round-trip as `MODEL`
//! records. Entries PDB cannot represent must fail with a typed PDB write
//! issue, never write a file that reads back differently.

use std::collections::BTreeSet;
use std::io::Read;
use std::path::{Path, PathBuf};

use kekule::diagnostics::{DiagnosticKind, FormatErrorDetail};
use kekule::macromolecular::{ModelSelection, StructureInterpretOptions};
use kekule::mmcif;
use kekule::pdb::{self, PdbWriteIssue};
use kekule::structure::{Ensemble, Model, ModelView};
use kekule::topology::{Topology, TopologyAtomIndex};
use kekule::units::{ANGSTROM, SQUARE_ANGSTROM};

fn curated_bio() -> Vec<PathBuf> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../datasets/bio/curated/pdb");
    let mut files = std::fs::read_dir(&root)
        .expect("readable dataset directory")
        .map(|entry| entry.expect("dataset entry").path())
        .filter(|path| path.to_string_lossy().ends_with(".cif.gz"))
        .collect::<Vec<_>>();
    files.sort();
    files
}

fn read_gzip(path: &Path) -> String {
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(path).expect("dataset file"))
        .read_to_string(&mut text)
        .expect("gzip mmCIF");
    text
}

/// The identity columns a PDB record carries for one atom.
fn identity(topology: &Topology, index: usize) -> (String, Option<String>, String, String, String) {
    let atom = topology
        .atom_at(TopologyAtomIndex::new(index as u32))
        .expect("dense atom");
    let site = atom.atom_site().expect("interpreted atoms have sites");
    let residue = atom.residue().expect("interpreted atoms have residues");
    let chain = atom.chain().expect("interpreted atoms have chains");
    let name = site
        .metadata()
        .auth_atom_id()
        .or(site.metadata().label_atom_id())
        .unwrap_or_default()
        .to_owned();
    let residue_name = residue
        .author_comp_id()
        .unwrap_or(residue.name())
        .to_owned();
    let chain_id = site
        .metadata()
        .auth_asym_id()
        .or(chain.author_id())
        .unwrap_or(chain.label_id())
        .to_owned();
    let sequence = residue
        .author_seq_id()
        .map(str::to_owned)
        .or_else(|| residue.label_seq_id().map(|value| value.to_string()))
        .unwrap_or_default();
    (
        name,
        residue.insertion_code().map(str::to_owned),
        residue_name,
        chain_id,
        sequence
            .trim()
            .parse::<i64>()
            .expect("integer residue number")
            .to_string(),
    )
}

fn bonds(topology: &Topology) -> BTreeSet<(usize, usize, String)> {
    topology
        .bonds()
        .map(|bond| {
            let [left, right] = bond.atoms().map(|atom| atom.index().index());
            (
                left.min(right),
                left.max(right),
                format!("{:?}", bond.order),
            )
        })
        .collect()
}

fn close(left: f64, right: f64, tolerance: f64) -> bool {
    (left - right).abs() <= tolerance
}

/// The reread structure keeps everything PDB carries at its precision.
fn assert_same_structure(reread: ModelView<'_>, source: ModelView<'_>, label: &str) {
    let (left, right) = (reread.topology(), source.topology());
    assert_eq!(left.atom_count(), right.atom_count(), "{label}: atom count");
    for index in 0..left.atom_count() {
        let dense = TopologyAtomIndex::new(index as u32);
        let (a, b) = (left.atom_at(dense).unwrap(), right.atom_at(dense).unwrap());
        assert_eq!(
            a.atom().element(),
            b.atom().element(),
            "{label}: element {index}"
        );
        assert_eq!(
            a.atom().formal_charge(),
            b.atom().formal_charge(),
            "{label}: charge {index}"
        );
        assert_eq!(
            identity(left, index),
            identity(right, index),
            "{label}: identity {index}"
        );
        let p = reread
            .positions()
            .get(index)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap();
        let q = source
            .positions()
            .get(index)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap();
        assert!(
            close(p.x, q.x, 5.1e-4) && close(p.y, q.y, 5.1e-4) && close(p.z, q.z, 5.1e-4),
            "{label}: coordinates {index}"
        );
        match (reread.occupancy(dense), source.occupancy(dense)) {
            (Some(p), Some(q)) => assert!(close(p, q, 5.1e-3), "{label}: occupancy {index}"),
            (p, q) => assert_eq!(p.is_some(), q.is_some(), "{label}: occupancy {index}"),
        }
        let b_factor = |view: ModelView<'_>| {
            view.b_factor(dense)
                .map(|value| value.value_in(SQUARE_ANGSTROM).unwrap())
        };
        match (b_factor(reread), b_factor(source)) {
            (Some(p), Some(q)) => assert!(close(p, q, 5.1e-3), "{label}: B factor {index}"),
            (p, q) => assert_eq!(p.is_some(), q.is_some(), "{label}: B factor {index}"),
        }
    }
    let (reread_bonds, source_bonds) = (bonds(left), bonds(right));
    let gained = reread_bonds.difference(&source_bonds).collect::<Vec<_>>();
    let lost = source_bonds.difference(&reread_bonds).collect::<Vec<_>>();
    assert!(
        gained.is_empty() && lost.is_empty(),
        "{label}: bonds gained {gained:?}, lost {lost:?}"
    );
    assert_eq!(
        left.instance_count(),
        right.instance_count(),
        "{label}: molecules"
    );
}

/// The PDB write issues an entry may legitimately hit: identifiers or
/// values PDB has no columns for.
fn is_representability_limit(issue: &PdbWriteIssue) -> bool {
    matches!(
        issue,
        PdbWriteIssue::ChainIdTooLong { .. }
            | PdbWriteIssue::ResidueNameTooLong { .. }
            | PdbWriteIssue::InvalidAtomName { .. }
            | PdbWriteIssue::InsertionCodeTooLong { .. }
            | PdbWriteIssue::InvalidResidueNumber { .. }
            | PdbWriteIssue::TooManyAtoms { .. }
            | PdbWriteIssue::CoordinateOutOfRange { .. }
            | PdbWriteIssue::DictionaryBondConflict { .. }
    )
}

/// Writes, rereads, and rewrites one model; `None` when PDB cannot represent
/// it.
fn round_trip_model(model: &Model, label: &str) -> Result<Model, String> {
    let written = match pdb::write(model, Default::default()) {
        Ok(written) => written,
        Err(error) => {
            let Some(FormatErrorDetail::PdbWrite(issue)) = error.detail() else {
                panic!("{label}: untyped PDB write failure: {error}");
            };
            assert_eq!(
                error.kind(),
                DiagnosticKind::Unsupported,
                "{label}: {error}"
            );
            assert!(is_representability_limit(issue), "{label}: {error}");
            return Err(error.to_string());
        }
    };
    let reread = pdb::parse_str(&written, Default::default())
        .and_then(|document| document.interpret(Default::default()))
        .unwrap_or_else(|error| panic!("{label}: written PDB does not read back: {error}"))
        .into_model();
    assert_same_structure(reread.as_model_view(), model.as_model_view(), label);
    let again = pdb::write(&reread, Default::default())
        .unwrap_or_else(|error| panic!("{label}: reread model does not write: {error}"));
    assert!(
        again == written,
        "{label}: PDB -> model -> PDB changed the text"
    );
    Ok(reread)
}

fn round_trip_ensemble(ensemble: &Ensemble, label: &str) {
    let written = pdb::write(ensemble, Default::default())
        .unwrap_or_else(|error| panic!("{label}: ensemble does not write: {error}"));
    let reread = pdb::parse_str(&written, Default::default())
        .and_then(|document| document.interpret_ensemble(Default::default()))
        .unwrap_or_else(|error| panic!("{label}: written ensemble does not read back: {error}"))
        .into_ensemble();
    assert_eq!(reread.len(), ensemble.len(), "{label}: members");
    for (index, (left, right)) in reread.iter().zip(ensemble.iter()).enumerate() {
        assert_same_structure(
            left.as_model_view(),
            right.as_model_view(),
            &format!("{label} member {index}"),
        );
    }
    let again = pdb::write(&reread, Default::default())
        .unwrap_or_else(|error| panic!("{label}: reread ensemble does not write: {error}"));
    assert!(
        again == written,
        "{label}: ensemble text changed on rewrite"
    );
}

#[test]
fn curated_biomolecules_round_trip_through_pdb() {
    let mut written = 0;
    let mut unrepresentable = Vec::new();
    let mut ensembles = 0;
    for path in curated_bio() {
        let label = path.file_name().unwrap().to_string_lossy().into_owned();
        let Ok(document) = mmcif::parse_str(&read_gzip(&path), Default::default()) else {
            continue;
        };
        let first =
            StructureInterpretOptions::default().with_model_selection(ModelSelection::First);
        let Ok(interpretation) = document.interpret(first) else {
            continue;
        };
        match round_trip_model(interpretation.model(), &label) {
            Ok(_) => written += 1,
            Err(reason) => {
                unrepresentable.push(format!("{label}: {reason}"));
                continue;
            }
        }
        if let Ok(interpretation) = document.interpret_ensemble(Default::default()) {
            if interpretation.ensemble().len() > 1 {
                round_trip_ensemble(interpretation.ensemble(), &label);
                ensembles += 1;
            }
        }
    }
    eprintln!(
        "{written} entries round-tripped, {ensembles} as ensembles; unrepresentable: {unrepresentable:#?}"
    );
    assert!(written >= 15, "only {written} entries round-tripped");
    // mmCIF and PDB readings agree on metal coordination, so every curated
    // entry is representable.
    assert!(
        unrepresentable.is_empty(),
        "unrepresentable entries: {unrepresentable:#?}"
    );
}
