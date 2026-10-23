//! Cooperative cancellation of macromolecular interpretation.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use kekule::diagnostics::{Cancellation, DiagnosticKind, FormatError};
use kekule::macromolecular::{
    StructureEnsembleInterpretOptions, StructureEnsembleInterpretation, StructureInterpretOptions,
    StructureInterpretation,
};
use kekule::structure::Model;
use kekule::{mmcif, pdb};

const ADENYLATE_KINASE_CIF: &str = include_str!("../fixtures/mmcif/1AKE.cif");
const ADENYLATE_KINASE_PDB: &str = include_str!("../fixtures/pdb/1AKE.pdb");

/// A cancellation requested from the `after`th poll on, counting every poll.
fn requested_after(after: usize) -> (Cancellation, Arc<AtomicUsize>) {
    let polls = Arc::new(AtomicUsize::new(0));
    let counter = Arc::clone(&polls);
    let cancellation =
        Cancellation::new(move || counter.fetch_add(1, Ordering::Relaxed) + 1 >= after);
    (cancellation, polls)
}

/// Whether two independently interpreted models are equal: every model is
/// published with its own layout identity, so the layouts are compared by
/// content.
fn same_model(left: &Model, right: &Model) -> bool {
    let (left_view, right_view) = (left.as_model_view(), right.as_model_view());
    left.topology().same_layout(right.topology())
        && left_view.positions() == right_view.positions()
        && left_view.cell() == right_view.cell()
        && left_view.occupancies() == right_view.occupancies()
        && left_view.b_factors() == right_view.b_factors()
}

fn same_interpretation(left: &StructureInterpretation, right: &StructureInterpretation) -> bool {
    left.report() == right.report() && same_model(left.model(), right.model())
}

fn same_ensemble(
    left: &StructureEnsembleInterpretation,
    right: &StructureEnsembleInterpretation,
) -> bool {
    left.reports() == right.reports()
        && left.topology().same_layout(right.topology())
        && left.ensemble().len() == right.ensemble().len()
}

/// Runs `interpret` with a check that never turns true, which must not change
/// the result, then with checks that turn true at the first, a middle, and
/// the last poll, each of which must fail with a cancelled error.
fn assert_cancellable<T>(
    interpret: impl Fn(Option<Cancellation>) -> Result<T, FormatError>,
    same: impl Fn(&T, &T) -> bool,
) {
    let expected = interpret(None).unwrap();
    let (never, polls) = requested_after(usize::MAX);
    assert!(same(&interpret(Some(never)).unwrap(), &expected));
    let polls = polls.load(Ordering::Relaxed);
    assert!(polls > 100, "interpretation polled only {polls} times");
    for after in [1, polls / 3, polls / 2, polls] {
        let (cancellation, _) = requested_after(after);
        let Err(error) = interpret(Some(cancellation)) else {
            panic!("poll {after}: interpretation completed after cancellation");
        };
        assert_eq!(
            error.kind(),
            DiagnosticKind::Cancelled,
            "poll {after}: {error}"
        );
        assert!(
            std::error::Error::source(&error)
                .is_some_and(|source| source.is::<kekule::diagnostics::Cancelled>()),
            "poll {after}"
        );
    }
}

#[test]
fn mmcif_interpretation_stops_once_cancellation_is_requested() {
    let document = mmcif::parse_str(ADENYLATE_KINASE_CIF, Default::default()).unwrap();
    assert_cancellable(
        |cancellation| {
            document.interpret(StructureInterpretOptions::default().with_cancellation(cancellation))
        },
        same_interpretation,
    );
    assert_cancellable(
        |cancellation| {
            document.blocks()[0]
                .interpret(StructureInterpretOptions::default().with_cancellation(cancellation))
        },
        same_interpretation,
    );
}

#[test]
fn mmcif_ensemble_interpretation_stops_once_cancellation_is_requested() {
    let document = mmcif::parse_str(ADENYLATE_KINASE_CIF, Default::default()).unwrap();
    assert_cancellable(
        |cancellation| {
            document.interpret_ensemble(
                StructureEnsembleInterpretOptions::default().with_cancellation(cancellation),
            )
        },
        same_ensemble,
    );
    assert_cancellable(
        |cancellation| {
            document.blocks()[0].interpret_ensemble(
                StructureEnsembleInterpretOptions::default().with_cancellation(cancellation),
            )
        },
        same_ensemble,
    );
}

#[test]
fn pdb_interpretation_stops_once_cancellation_is_requested() {
    let document = pdb::parse_str(ADENYLATE_KINASE_PDB, Default::default()).unwrap();
    assert_cancellable(
        |cancellation| {
            document.interpret(StructureInterpretOptions::default().with_cancellation(cancellation))
        },
        same_interpretation,
    );
    assert_cancellable(
        |cancellation| {
            document.interpret_ensemble(
                StructureEnsembleInterpretOptions::default().with_cancellation(cancellation),
            )
        },
        same_ensemble,
    );
    let (cancellation, _) = requested_after(1);
    let error = document
        .interpret(StructureInterpretOptions::default().with_cancellation(Some(cancellation)))
        .unwrap_err();
    assert_eq!(error.format(), kekule::diagnostics::Format::Pdb);
}

#[test]
fn cancellation_is_an_identity_option() {
    let (cancellation, _) = requested_after(usize::MAX);
    let options =
        StructureInterpretOptions::default().with_cancellation(Some(cancellation.clone()));
    assert_eq!(options.cancellation(), Some(&cancellation));
    assert_eq!(
        options,
        StructureInterpretOptions::default().with_cancellation(Some(cancellation))
    );
    let (other, _) = requested_after(usize::MAX);
    assert_ne!(
        options,
        StructureInterpretOptions::default().with_cancellation(Some(other))
    );
    assert_eq!(StructureInterpretOptions::default().cancellation(), None);
    assert_eq!(
        StructureEnsembleInterpretOptions::default().cancellation(),
        None
    );
}
