//! SMILES output is a faithful, numbering-independent serialization.

use kekule::smiles::{self, MolWriteErrorKind, SmilesWriteOptions};
use kekule::topology::Topology;

use crate::contract;
use crate::corpus::{assert_invariant, molecules, regression_mixtures};
use crate::support::{cip_labels, renumbered};

/// Canonical and isomeric output exist for every corpus molecule and re-read
/// to the same composition and canonical SMILES; ordinary output does too
/// unless the molecule carries stereo it cannot encode. Axial stereo is the one
/// documented gap: writers must reject it explicitly instead of dropping it.
#[test]
fn smiles_writers_round_trip_every_corpus_molecule() {
    let mut failures = Vec::new();
    for sample in molecules() {
        let molecule = sample.perceived();
        let mut fail = |message: String| failures.push((sample.label.clone(), message));
        if sample.has_axis_stereo() {
            for options in [
                SmilesWriteOptions::isomeric(),
                SmilesWriteOptions::canonical(),
            ] {
                match smiles::write(&molecule, options) {
                    Err(error) if error.kind() == MolWriteErrorKind::UnsupportedRepresentation => {}
                    other => fail(format!(
                        "axial stereo must be rejected explicitly, got {other:?}"
                    )),
                }
            }
            continue;
        }
        let canonical = match smiles::write(&molecule, SmilesWriteOptions::canonical()) {
            Ok(canonical) => canonical,
            Err(error) => {
                fail(format!("canonical writer failed: {error}"));
                continue;
            }
        };
        let has_stereo = molecule.stereo_elements().next().is_some();
        for (writer, options) in [
            ("canonical", SmilesWriteOptions::canonical()),
            ("isomeric", SmilesWriteOptions::isomeric()),
            ("ordinary", SmilesWriteOptions::ordinary()),
        ] {
            let written = match smiles::write(&molecule, options) {
                Ok(written) => written,
                Err(error)
                    if writer == "ordinary"
                        && has_stereo
                        && error.kind() == MolWriteErrorKind::UnsupportedRepresentation =>
                {
                    continue
                }
                Err(error) => {
                    fail(format!("{writer} writer failed: {error}"));
                    continue;
                }
            };
            let checked = std::panic::catch_unwind(|| {
                contract::assert_output(&molecule, &written, Some(&canonical))
            });
            if let Err(panic) = checked {
                fail(format!("{writer}: {}", panic_message(&panic)));
            }
        }
    }
    assert_invariant("smiles_round_trip", failures);
}

/// Canonical SMILES and CIP descriptors describe the molecule, not the order
/// in which its atoms, bonds and bond endpoints were stored.
#[test]
fn canonical_smiles_and_cip_labels_are_invariant_under_atom_renumbering() {
    let mut failures = Vec::new();
    let (mut renumberings, mut identities) = (0, 0);
    for sample in molecules() {
        let mut molecule = sample.perceived();
        let canonical = smiles::write(&molecule, SmilesWriteOptions::canonical()).ok();
        let labels = cip_labels(&mut molecule);
        for seed in [1, 2, 3] {
            let (mut copy, mapping) = renumbered(&sample.molecule, seed);
            if mapping.len() > 2 {
                renumberings += 1;
                if mapping
                    .iter()
                    .enumerate()
                    .all(|(old, new)| new.index() == old)
                {
                    identities += 1;
                }
            }
            copy.perceive()
                .unwrap_or_else(|error| panic!("{} seed {seed}: {error}", sample.label));
            let copy_canonical = smiles::write(&copy, SmilesWriteOptions::canonical()).ok();
            if copy_canonical != canonical {
                failures.push((
                    sample.label.clone(),
                    format!("seed {seed}: canonical {copy_canonical:?} != {canonical:?}"),
                ));
            }
            let expected = labels
                .iter()
                .map(|(focus, label)| (focus.mapped(&mapping), label.clone()))
                .collect();
            let actual = cip_labels(&mut copy);
            if actual != expected {
                failures.push((
                    sample.label.clone(),
                    format!("seed {seed}: CIP {actual:?} != {expected:?}"),
                ));
            }
        }
    }
    // A molecule may occasionally draw the identity permutation; most must not.
    assert!(
        identities * 10 < renumberings,
        "{identities} of {renumberings} kept numbering"
    );
    assert_invariant("renumbering", failures);
}

/// A topology's canonical SMILES sorts its components; the order in which a
/// mixture or salt was supplied must not matter.
#[test]
fn canonical_topology_smiles_is_invariant_under_component_order() {
    let mut failures = Vec::new();
    for (label, mut components) in regression_mixtures() {
        if components.len() < 2 {
            continue;
        }
        for component in &mut components {
            component.perceive().unwrap();
        }
        let canonical = |molecules: Vec<_>| {
            smiles::write(
                &Topology::from_molecules(molecules).unwrap(),
                SmilesWriteOptions::canonical(),
            )
            .unwrap()
        };
        let expected = canonical(components.clone());
        components.reverse();
        let reversed = canonical(components.clone());
        components.rotate_left(1);
        let rotated = canonical(components);
        if reversed != expected || rotated != expected {
            failures.push((
                label,
                format!("{expected} vs reversed {reversed} vs rotated {rotated}"),
            ));
        }
    }
    assert_invariant("component_order", failures);
}

fn panic_message(panic: &Box<dyn std::any::Any + Send>) -> String {
    panic
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| panic.downcast_ref::<&str>().map(|text| (*text).to_owned()))
        .unwrap_or_else(|| "non-string panic".to_owned())
}
