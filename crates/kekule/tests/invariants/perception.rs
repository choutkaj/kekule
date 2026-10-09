//! Perception derives state; it never rewrites what a molecule represents.

use kekule::core::Perception;
use kekule::smiles::{self, SmilesWriteOptions};

use crate::corpus::{assert_invariant, molecules};
use crate::support::{cip_labels, try_cip_labels};

/// Perception, reperception and CIP assignment leave represented chemistry
/// untouched, and repeating them reproduces the same derived state.
#[test]
fn perception_is_idempotent_and_never_rewrites_represented_chemistry() {
    let mut failures = Vec::new();
    for sample in molecules() {
        let mut fail = |message: &str| failures.push((sample.label.clone(), message.to_owned()));
        let mut molecule = sample.perceived();
        if molecule != sample.molecule {
            fail("perception changed the molecule");
        }
        let first = molecule.perception().clone();
        if first == Perception::default() {
            fail("perception installed nothing");
        }
        molecule.perceive().unwrap();
        if molecule.perception() != &first {
            fail("reperception differs");
        }
        let labels = cip_labels(&mut molecule);
        if molecule != sample.molecule {
            fail("CIP assignment changed the molecule");
        }
        if cip_labels(&mut molecule) != labels {
            fail("CIP reassignment differs");
        }
    }
    assert_invariant("perception", failures);
}

/// Materializing implicit hydrogens as atoms and collapsing them again keeps
/// the molecule's identity and stereo descriptors.
#[test]
fn explicit_hydrogen_round_trip_preserves_identity() {
    let mut failures = Vec::new();
    for sample in molecules() {
        let mut fail = |message: String| failures.push((sample.label.clone(), message));
        let mut original = sample.perceived();
        let canonical = smiles::write(&original, SmilesWriteOptions::canonical()).ok();
        let mut labels = cip_labels(&mut original).into_values().collect::<Vec<_>>();
        labels.sort();
        let mut explicit = original.clone();
        explicit
            .add_hydrogens()
            .unwrap_or_else(|error| panic!("{}: {error}", sample.label));
        explicit.perceive().unwrap();
        let mut collapsed = explicit.clone();
        collapsed
            .remove_hydrogens()
            .unwrap_or_else(|error| panic!("{}: {error}", sample.label));
        collapsed.perceive().unwrap();
        for (stage, molecule) in [("explicit", &mut explicit), ("collapsed", &mut collapsed)] {
            match try_cip_labels(molecule) {
                Ok(actual) => {
                    let mut actual = actual.into_values().collect::<Vec<_>>();
                    actual.sort();
                    if actual != labels {
                        fail(format!("{stage}: CIP {actual:?} != {labels:?}"));
                    }
                }
                Err(error) => fail(format!("{stage}: CIP failed: {error}")),
            }
        }
        let collapsed_canonical = smiles::write(&collapsed, SmilesWriteOptions::canonical()).ok();
        if collapsed_canonical != canonical {
            fail(format!(
                "collapsed canonical {collapsed_canonical:?} != {canonical:?}"
            ));
        }
    }
    assert_invariant("hydrogen_round_trip", failures);
}
