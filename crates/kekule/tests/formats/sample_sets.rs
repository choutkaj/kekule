//! Readers that interpret sample sets give what the readers of ensembles,
//! series and trajectories give, converted.

use kekule::diagnostics::{DiagnosticKind, FormatErrorDetail, UnrepresentablePolicy};
use kekule::mmcif::{self, MmcifWriteOptions};
use kekule::native::{self, NativeWriteOptions};
use kekule::pdb::{self, PdbWriteOptions};
use kekule::sdf::{self, SdfRecordInterpretation, SdfWriteOptions};
use kekule::structure::{
    Ensemble, EnsemblePartition, Model, Order, Positions, SampleSet, Series, SeriesAxis, SeriesItem,
};
use kekule::topology::Topology;
use kekule::units::{Quantity, ANGSTROM, DEGREE, KILOJOULE_PER_MOLE};

const WATERS: &str = "3
first
O 0.0 0.0 0.0
H 0.957 0.0 0.0
H -0.24 0.927 0.0
3
second
O 0.0 0.0 0.1
H 0.957 0.0 0.1
H -0.24 0.927 0.1
";

/// The topology's exact native encoding, which compares every part of it.
fn topology_bytes(topology: &Topology) -> Vec<u8> {
    native::write(topology, NativeWriteOptions::default()).unwrap()
}

/// Asserts that two ensembles hold the same topology, samples, partition,
/// states, temperature, and properties.
fn assert_same_ensemble(actual: &Ensemble, expected: &Ensemble) {
    assert_eq!(
        topology_bytes(actual.topology()),
        topology_bytes(expected.topology())
    );
    assert_eq!(actual.properties(), expected.properties());
    assert_eq!(actual.len(), expected.len());
    for (left, right) in actual.iter().zip(expected.iter()) {
        assert_eq!(left.payload(), right.payload());
    }
    assert_eq!(actual.partition(), expected.partition());
    assert!(actual
        .states()
        .zip(expected.states())
        .all(|(left, right)| left.payload() == right.payload()));
    assert_eq!(actual.temperature(), expected.temperature());
}

/// Two waters as a two-model PDB entry.
fn pdb_waters() -> String {
    let ensemble = kekule::xyz::parse_str(WATERS, Default::default())
        .unwrap()
        .interpret_ensemble(Default::default())
        .unwrap()
        .into_ensemble();
    pdb::write(
        &ensemble,
        PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit),
    )
    .unwrap()
    .0
}

#[test]
fn pdb_and_mmcif_sample_sets_are_their_ensembles_of_one_unweighted_state() {
    let document = pdb::parse_str(&pdb_waters(), Default::default()).unwrap();
    let ensemble = document.interpret_ensemble(Default::default()).unwrap();
    let set = document.interpret_sample_set(Default::default()).unwrap();
    assert_eq!(set.reports(), ensemble.reports());
    assert_eq!(set.sample_set().len(), 2);
    assert_eq!(set.sample_set().order(), None);
    let distribution = set.sample_set().distribution().unwrap();
    assert_eq!(distribution.state_count(), 1);
    assert_eq!(distribution.weights(), None);
    assert_same_ensemble(
        &set.sample_set().clone().into_ensemble().unwrap(),
        ensemble.ensemble(),
    );
    // Both write the same PDB entry.
    assert_eq!(
        pdb::write(&set, PdbWriteOptions::default()).unwrap(),
        pdb::write(&ensemble, PdbWriteOptions::default()).unwrap()
    );

    let text = mmcif::write([&ensemble], MmcifWriteOptions::default())
        .unwrap()
        .0;
    assert_eq!(
        mmcif::write([&set], MmcifWriteOptions::default())
            .unwrap()
            .0,
        text
    );
    let document = mmcif::parse_str(&text, Default::default()).unwrap();
    let ensemble = document.interpret_ensemble(Default::default()).unwrap();
    let set = document.interpret_sample_set(Default::default()).unwrap();
    assert_eq!(set.reports(), ensemble.reports());
    assert_same_ensemble(
        &set.sample_set().clone().into_ensemble().unwrap(),
        ensemble.ensemble(),
    );
    let block = &document.blocks()[0];
    let from_block = block.interpret_sample_set(Default::default()).unwrap();
    assert_eq!(from_block.reports(), ensemble.reports());
    assert_same_ensemble(
        &from_block.into_sample_set().into_ensemble().unwrap(),
        ensemble.ensemble(),
    );
}

fn interpret_records(text: &str) -> Vec<SdfRecordInterpretation> {
    sdf::parse_str(text, Default::default())
        .unwrap()
        .records()
        .iter()
        .map(|record| record.interpret(Default::default()).unwrap())
        .collect()
}

/// Two methane conformers in two states, weighted 1:3, the second at
/// -2.5 kJ/mol.
fn weighted_methanes() -> Ensemble {
    let topology = kekule::smiles::to_topology("C").unwrap();
    let models = [0.0, 1.0].map(|x| {
        Model::new(
            topology.clone(),
            Positions::new(Quantity::new(
                [kekule::geometry::Point3::new(x, 0.0, 0.0)],
                ANGSTROM,
            ))
            .unwrap(),
        )
        .unwrap()
    });
    let mut ensemble =
        Ensemble::from_models(models, EnsemblePartition::state_per_sample(2)).unwrap();
    ensemble.get_mut(1).unwrap().set_weight(3.0).unwrap();
    ensemble
        .state_mut(1)
        .unwrap()
        .set_energy(Some(Quantity::new(-2.5, KILOJOULE_PER_MOLE)))
        .unwrap();
    ensemble
}

/// A three-item series of one methane carbon along `axis`.
fn methane_series(axis: SeriesAxis) -> Series {
    let topology = kekule::smiles::to_topology("C").unwrap();
    let mut series = Series::new(topology, axis);
    for index in 0..3u32 {
        let positions = Positions::zeros(1);
        let item = match axis {
            SeriesAxis::Step => SeriesItem::at_step(positions, u64::from(index) * 4),
            SeriesAxis::Coordinate(unit) => SeriesItem::at_coordinate(
                positions,
                Quantity::new(f64::from(index) * 30.0 - 60.0, unit),
            )
            .unwrap(),
            _ => unreachable!("a step or coordinate series"),
        };
        series.push(item).unwrap();
    }
    series
}

#[test]
fn sdf_sample_sets_are_the_ensembles_and_series_of_their_records() {
    // Weights, states and energies make the distribution.
    let records =
        interpret_records(&sdf::write(&weighted_methanes(), SdfWriteOptions::default()).unwrap());
    let set = sdf::sample_set(&records).unwrap();
    assert_eq!(set.order(), None);
    assert_same_ensemble(
        &set.clone().into_ensemble().unwrap(),
        &sdf::ensemble(&records).unwrap(),
    );
    assert_eq!(set.distribution().unwrap().state_count(), 2);

    // Steps or coordinates make the order.
    for axis in [SeriesAxis::Step, SeriesAxis::Coordinate(DEGREE)] {
        let series = methane_series(axis);
        let records = interpret_records(&sdf::write(&series, SdfWriteOptions::default()).unwrap());
        let set = sdf::sample_set(&records).unwrap();
        assert_eq!(
            set.order(),
            SampleSet::from_series(series.clone()).order(),
            "{axis:?}"
        );
        assert!(matches!(
            (axis, set.order()),
            (SeriesAxis::Step, Some(Order::Step(_)))
                | (SeriesAxis::Coordinate(_), Some(Order::Coordinate { .. }))
        ));
        // Without its order, the set is the ensemble of the same records.
        let ensemble = sdf::ensemble(&records).unwrap();
        let mut converted = set.clone();
        converted.set_order(None).unwrap();
        assert_same_ensemble(&converted.into_ensemble().unwrap(), &ensemble);
        // Without its distribution, it is the series.
        let mut set = set;
        set.set_distribution(None).unwrap();
        let read = set.into_series().unwrap();
        assert_eq!(read.axis(), series.axis());
        for (left, right) in read.iter().zip(series.iter()) {
            assert_eq!(left.payload(), right.payload());
        }
    }
}

#[test]
fn sdf_records_that_do_not_form_an_order_are_refused() {
    let series = methane_series(SeriesAxis::Step);
    let text = sdf::write(&series, SdfWriteOptions::default()).unwrap();
    let mut records = interpret_records(&text);
    // A record without the step its neighbours have.
    let stepless = interpret_records(
        &sdf::write(
            [sdf::SdfRecordSource::model(
                &series.get(0).unwrap().to_model(),
            )],
            SdfWriteOptions::default(),
        )
        .unwrap(),
    );
    records.push(stepless[0].clone());
    let error = sdf::sample_set(&records).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert!(
        error
            .to_string()
            .contains("missing. Other records have it."),
        "{error}"
    );
    assert!(sdf::ensemble(&records).is_ok());

    // Steps and coordinates together.
    let scan = methane_series(SeriesAxis::Coordinate(ANGSTROM));
    let mut mixed = interpret_records(&sdf::write(&scan, SdfWriteOptions::default()).unwrap());
    mixed.extend(interpret_records(&text));
    let error = sdf::sample_set(&mixed).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert!(matches!(
        error.detail(),
        Some(FormatErrorDetail::DataField { .. })
    ));
}

#[test]
fn xyz_sample_sets_have_no_distribution_where_ensembles_have_one_state() {
    let document = kekule::xyz::parse_str(WATERS, Default::default()).unwrap();
    let ensemble = document.interpret_ensemble(Default::default()).unwrap();
    let set = document.interpret_sample_set(Default::default()).unwrap();
    assert_eq!(set.reports(), ensemble.reports());
    assert_eq!(set.sample_set().order(), None);
    assert_eq!(set.sample_set().distribution(), None);
    assert_eq!(ensemble.ensemble().state_count(), 1);
    assert!(ensemble
        .ensemble()
        .iter()
        .all(|sample| sample.weight() == 1.0));
    assert_eq!(
        topology_bytes(set.topology()),
        topology_bytes(ensemble.topology())
    );
    for (left, right) in set.sample_set().iter().zip(ensemble.ensemble().iter()) {
        let (left, right) = (left.as_model_view(), right.as_model_view());
        assert_eq!(left.positions(), right.positions());
        assert_eq!(left.properties(), right.properties());
    }
    // With a distribution of one unweighted state it is the ensemble.
    let mut set = set.into_sample_set();
    set.set_distribution(Some(
        kekule::structure::Distribution::new(
            kekule::structure::Members::All,
            EnsemblePartition::one_state(2),
        )
        .unwrap(),
    ))
    .unwrap();
    assert_same_ensemble(&set.into_ensemble().unwrap(), ensemble.ensemble());
}
