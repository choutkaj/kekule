//! Energy surfaces: free energies and means of a distribution's samples on
//! one or two binned axes.

use kekule::properties::{PropertyColumn, PropertyKey};
use kekule::statistics::{
    Binning, CumulantOrder, Estimator, GridAxis, Landscape, LandscapeError, LandscapeOptions,
    SurfaceQuantity,
};
use kekule::structure::{
    Distribution, Ensemble, Geometry, Partition, Positions, Samples, Sampling, StructureSet,
    ThermodynamicState,
};
use kekule::units::{Quantity, DEGREE, DIMENSIONLESS, KELVIN, KILOJOULE_PER_MOLE, NANOMETER};

/// The molar gas constant in kJ/(mol K), exact since the 2019 SI.
const R: f64 = 8.314_462_618_153_24e-3;

const ESTIMATORS: [Estimator; 4] = [
    Estimator::Exponential,
    Estimator::Cumulant(CumulantOrder::First),
    Estimator::Cumulant(CumulantOrder::Second),
    Estimator::Cumulant(CumulantOrder::Third),
];

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} differs from {expected}"
    );
}

fn kelvin(value: f64) -> Quantity<f64> {
    Quantity::new(value, KELVIN)
}

fn at(value: f64) -> ThermodynamicState {
    ThermodynamicState::canonical(kelvin(value)).unwrap()
}

fn every_sample(count: usize) -> Distribution {
    Distribution::new(Samples::All, Partition::one_state(count)).unwrap()
}

fn options() -> LandscapeOptions {
    LandscapeOptions::new().with_fewest(1)
}

fn at_300() -> LandscapeOptions {
    options().with_temperature(Some(kelvin(300.0))).unwrap()
}

/// Unit bins from 0, `bins` of them.
fn unit_bins(bins: usize) -> Binning {
    Binning::new(0.0, 1.0, bins).unwrap()
}

/// A run sampled at 300 K whose reduced bias is `reduced`, one value per
/// geometry, read as one state with `estimator`.
fn reweighted(reduced: &[f64], estimator: Estimator) -> Ensemble {
    let topology = kekule::smiles::to_topology("O").unwrap();
    let atoms = topology.atom_count();
    let geometries = reduced
        .iter()
        .map(|_| Geometry::new(Positions::zeros(atoms)));
    let mut set = StructureSet::from_geometries(topology, geometries).unwrap();
    let boost = PropertyKey::new("boost").unwrap();
    let thermal_energy = *at(300.0).thermal_energy().value();
    set.columns_mut()
        .insert(
            boost.clone(),
            PropertyColumn::Real {
                unit: KILOJOULE_PER_MOLE,
                values: reduced
                    .iter()
                    .map(|value| Some(value * thermal_energy))
                    .collect(),
            },
        )
        .unwrap();
    set.set_sampling(Some(Sampling::new(at(300.0)).with_bias([boost]).unwrap()))
        .unwrap();
    Ensemble::reweighted(set, None, estimator).unwrap()
}

#[test]
fn equal_weights_count_the_samples() -> Result<(), LandscapeError> {
    // Counts 1, 2, 4, and 0 in four unit bins.
    let values = [0.5, 1.2, 1.7, 2.0, 2.1, 2.5, 2.9].map(Some);
    let axis = GridAxis::new(&values, unit_bins(4), NANOMETER);
    let surface = Landscape::free_energy(&every_sample(7), &[axis], &at_300())?;
    assert_eq!(surface.quantity(), SurfaceQuantity::FreeEnergy);
    assert_eq!(surface.shape(), [4]);
    assert_eq!(surface.value_unit(), KILOJOULE_PER_MOLE);
    assert_eq!(surface.units(), [NANOMETER]);
    assert_eq!(surface.temperature(), Some(kelvin(300.0)));
    let rt = R * 300.0;
    let cells = surface.cells();
    assert_eq!(
        cells.iter().map(|cell| cell.count()).collect::<Vec<_>>(),
        [1, 2, 4, 0]
    );
    assert_close(cells[0].value().unwrap(), rt * 4.0_f64.ln());
    assert_close(cells[1].value().unwrap(), rt * 2.0_f64.ln());
    assert_eq!(cells[2].value(), Some(0.0));
    assert_eq!(cells[2].effective_count(), 4.0);
    // An empty cell has nothing.
    assert_eq!(cells[3].value(), None);
    assert_eq!(cells[3].log_population(), None);
    assert_eq!(cells[3].representative(), None);
    assert_eq!(cells[3].effective_count(), 0.0);
    let probabilities = surface.probabilities();
    for (probability, count) in probabilities.iter().zip([1.0, 2.0, 4.0, 0.0]) {
        assert_close(*probability, count / 7.0);
    }
    Ok(())
}

#[test]
fn a_constant_reduced_bias_per_region_gives_the_exact_free_energy() -> Result<(), LandscapeError> {
    // Three samples with b = 0 in the first bin and two with b = 1.5 in the
    // second; every estimator sums a constant exactly.
    let reduced = [0.0, 0.0, 1.5, 0.0, 1.5];
    let values = [0.5, 0.25, 1.5, 0.75, 1.25].map(Some);
    let rt = R * 300.0;
    let low = 3.0_f64.ln();
    let high = 2.0_f64.ln() + 1.5;
    for estimator in ESTIMATORS {
        let ensemble = reweighted(&reduced, estimator);
        let axis = GridAxis::new(&values, unit_bins(2), NANOMETER);
        // The target, 300 K, is the free energy's temperature.
        let surface = Landscape::free_energy(ensemble.distribution(), &[axis], &options())?;
        assert_eq!(surface.temperature(), Some(kelvin(300.0)));
        let cells = surface.cells();
        assert_close(cells[0].value().unwrap(), rt * (high - low));
        assert_close(cells[1].value().unwrap(), 0.0);
        let total = low.exp() + high.exp();
        assert_close(surface.probabilities()[0], low.exp() / total);
        assert_close(cells[1].log_population().unwrap(), high - total.ln());
        assert_close(cells[0].effective_count(), 3.0);
        assert_close(cells[1].effective_count(), 2.0);
    }
    Ok(())
}

#[test]
fn a_bin_that_is_a_state_has_the_state_population() -> Result<(), Box<dyn std::error::Error>> {
    let reduced = [0.0, 0.7, 1.9, 0.3, 2.4, 1.1];
    let values = [0.2, 0.4, 1.6, 0.9, 1.1, 1.3].map(Some);
    for estimator in ESTIMATORS {
        let mut ensemble = reweighted(&reduced, estimator);
        ensemble
            .distribution_mut()
            .regroup(Partition::new(vec![0, 0, 1, 0, 1, 1], vec![0, 2])?)?;
        let distribution = ensemble.distribution();
        let axis = GridAxis::new(&values, unit_bins(2), NANOMETER);
        let surface = Landscape::free_energy(distribution, &[axis], &options())?;
        let states = distribution.state_probabilities()?;
        for (bin, state) in surface.probabilities().iter().zip(states) {
            assert_close(*bin, state);
        }
    }
    Ok(())
}

#[test]
fn periodic_axes_wrap_at_the_period() -> Result<(), LandscapeError> {
    let dihedral = Binning::periodic(-180.0, 360.0, 36)?;
    assert!(dihedral.is_periodic());
    assert_eq!(dihedral.bin(180.0), Some(0));
    assert_eq!(dihedral.bin(-180.0), Some(0));
    assert_eq!(dihedral.bin(540.0), Some(0));
    assert_eq!(dihedral.bin(-179.9), Some(0));
    assert_eq!(dihedral.bin(179.9), Some(35));
    assert_eq!(dihedral.bin(-180.1), Some(35));
    assert_eq!(dihedral.bin(f64::NAN), None);
    assert_eq!(dihedral.bin(f64::INFINITY), None);
    assert_eq!(dihedral.center(35), Some(175.0));
    assert_eq!(dihedral.center(36), None);

    // Bins centred at -135°, -45°, 45°, and 135°. The sample at 200° wraps to
    // -160°, 25° from the first centre, nearer than the one at -100°.
    let quarters = Binning::periodic(-180.0, 360.0, 4)?;
    let values = [Some(-100.0), Some(200.0), Some(-170.0)];
    let axis = GridAxis::new(&values, quarters, DEGREE);
    let surface = Landscape::free_energy(&every_sample(3), &[axis], &at_300())?;
    assert_eq!(surface.cells()[0].count(), 3);
    assert_eq!(surface.cells()[0].representative(), Some(1));
    assert_eq!(surface.cell_at(&[200.0]), Some(0));
    assert_eq!(surface.outside(), 0);
    Ok(())
}

#[test]
fn masking_hides_values_only() -> Result<(), LandscapeError> {
    // Counts 1, 2, and 4; the single sample weighs most.
    let values = [0.5, 1.2, 1.7, 2.0, 2.1, 2.5, 2.9].map(Some);
    let distribution = every_sample(7)
        .with_log_weights(vec![5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        .unwrap();
    let axis = GridAxis::new(&values, unit_bins(3), NANOMETER);
    let masked = Landscape::free_energy(&distribution, &[axis], &at_300().with_fewest(3))?;
    let unmasked = Landscape::free_energy(&distribution, &[axis], &at_300())?;
    assert_eq!(masked.fewest(), 3);
    assert_eq!(masked.probabilities(), unmasked.probabilities());
    let cells = masked.cells();
    assert_eq!(cells[0].value(), None);
    assert_eq!(cells[1].value(), None);
    // The lowest value is the most populated cell that is not masked.
    assert_eq!(cells[2].value(), Some(0.0));
    for (cell, other) in cells.iter().zip(unmasked.cells()) {
        assert_eq!(cell.count(), other.count());
        assert_eq!(cell.log_population(), other.log_population());
        assert_eq!(cell.effective_count(), other.effective_count());
        assert_eq!(cell.representative(), other.representative());
    }
    assert!(unmasked.cells()[0].value() == Some(0.0));
    assert!(unmasked.cells()[2].value().unwrap() > 0.0);
    // The default needs ten samples.
    assert_eq!(LandscapeOptions::default().fewest(), 10);
    let none = Landscape::free_energy(
        &distribution,
        &[axis],
        &LandscapeOptions::default().with_temperature(Some(kelvin(300.0)))?,
    )?;
    assert!(none.cells().iter().all(|cell| cell.value().is_none()));
    Ok(())
}

#[test]
fn representatives_are_nearest_the_centre() -> Result<(), LandscapeError> {
    // 0.25 and 0.75 lie a quarter bin from the centre, 0.4 nearer.
    let values = [0.25, 0.75, 0.4].map(Some);
    let axis = GridAxis::new(&values, unit_bins(1), NANOMETER);
    let nearest = Landscape::free_energy(&every_sample(3), &[axis], &at_300())?;
    assert_eq!(nearest.cells()[0].representative(), Some(2));

    // Equally near: the lower geometry index with equal weights.
    let values = [0.25, 0.75];
    let values = values.map(Some);
    let axis = GridAxis::new(&values, unit_bins(1), NANOMETER);
    let tie = Landscape::free_energy(&every_sample(2), &[axis], &at_300())?;
    assert_eq!(tie.cells()[0].representative(), Some(0));
    // The heavier sample otherwise, under every estimator.
    for estimator in ESTIMATORS {
        let ensemble = reweighted(&[0.0, 0.5], estimator);
        let tie = Landscape::free_energy(ensemble.distribution(), &[axis], &options())?;
        assert_eq!(tie.cells()[0].representative(), Some(1));
    }
    Ok(())
}

#[test]
fn the_temperature_is_the_target_s() {
    let values = [Some(0.5)];
    let axis = GridAxis::new(&values, unit_bins(1), NANOMETER);
    // Neither a target nor a temperature.
    let error = Landscape::free_energy(&every_sample(1), &[axis], &options()).unwrap_err();
    assert_eq!(error, LandscapeError::MissingTemperature);
    assert_eq!(
        error.to_string(),
        "The free energy has no temperature. Give one."
    );
    // A mean needs none.
    let mean = Landscape::mean(&every_sample(1), &[axis], &values, NANOMETER, &options()).unwrap();
    assert_eq!(mean.temperature(), None);
    // The given temperature without a target.
    let explicit = options().with_temperature(Some(kelvin(310.0))).unwrap();
    let surface = Landscape::free_energy(&every_sample(1), &[axis], &explicit).unwrap();
    assert_eq!(surface.temperature(), Some(kelvin(310.0)));
    // The target's own.
    let targeted = every_sample(1).with_target(at(300.0)).unwrap();
    let surface = Landscape::free_energy(&targeted, &[axis], &options()).unwrap();
    assert_eq!(surface.temperature(), Some(kelvin(300.0)));
    let same = options().with_temperature(Some(kelvin(300.0))).unwrap();
    let surface = Landscape::free_energy(&targeted, &[axis], &same).unwrap();
    assert_eq!(surface.temperature(), Some(kelvin(300.0)));
    // Another one conflicts, for a mean too.
    let error = Landscape::free_energy(&targeted, &[axis], &explicit).unwrap_err();
    assert_eq!(
        error,
        LandscapeError::ConflictingTemperature {
            given: 310.0,
            target: 300.0
        }
    );
    assert_eq!(
        error.to_string(),
        "The temperature 310 K differs from the target temperature 300 K."
    );
    assert!(Landscape::mean(&targeted, &[axis], &values, NANOMETER, &explicit).is_err());
    // Temperatures that are none.
    assert_eq!(
        options().with_temperature(Some(kelvin(-1.0))).unwrap_err(),
        LandscapeError::InvalidTemperature
    );
    let error = options()
        .with_temperature(Some(Quantity::new(1.0, KILOJOULE_PER_MOLE)))
        .unwrap_err();
    assert!(matches!(error, LandscapeError::Unit(_)));
    assert!(std::error::Error::source(&error).is_some());
    assert_eq!(
        error.to_string(),
        "The landscape temperature has the wrong unit."
    );
}

#[test]
fn covering_bins_hold_every_value() {
    assert_eq!(
        Binning::covering(&[Some(1.0)], 0),
        Err(LandscapeError::NoBins)
    );
    assert_eq!(Binning::covering(&[], 3), Err(LandscapeError::NoValues));
    assert_eq!(
        Binning::covering(&[None, Some(f64::NAN), Some(f64::INFINITY)], 3),
        Err(LandscapeError::NoValues)
    );
    assert_eq!(
        Binning::covering(&[Some(-1e308), Some(1e308)], 3),
        Err(LandscapeError::InvalidBinning)
    );
    // One value: unit bins around it.
    let single = Binning::covering(&[Some(2.0), Some(2.0)], 4).unwrap();
    assert_eq!((single.start(), single.width()), (0.0, 1.0));
    assert_eq!(single.bin(2.0), Some(2));
    // Non-finite values are ignored.
    let covering = Binning::covering(&[Some(f64::NAN), Some(0.0), None, Some(1.0)], 3).unwrap();
    assert_eq!(covering.bin(0.0), Some(0));
    assert_eq!(covering.bin(1.0), Some(2));
    assert_eq!(covering.bin(1.0 + 1e-9), None);
    assert_eq!(covering.bin(-1e-9), None);
    // The largest value lands in the last bin whatever the rounding.
    for bins in 1..60 {
        for (lowest, highest) in [
            (0.1, 0.7),
            (-3.3, 1e-3),
            (1e10, 1e10 + 7.77),
            (-180.0, 179.9),
            (0.0, 1.0),
            (1e-300, 3e-300),
        ] {
            let binning = Binning::covering(&[Some(highest), Some(lowest)], bins).unwrap();
            assert_eq!(binning.bin(lowest), Some(0), "{lowest} in {bins}");
            assert_eq!(binning.bin(highest), Some(bins - 1), "{highest} in {bins}");
            assert!(!binning.is_periodic());
        }
    }
}

#[test]
fn binnings_are_finite() {
    assert_eq!(Binning::new(0.0, 1.0, 0), Err(LandscapeError::NoBins));
    assert_eq!(Binning::periodic(0.0, 1.0, 0), Err(LandscapeError::NoBins));
    for (start, width) in [
        (0.0, 0.0),
        (0.0, -1.0),
        (f64::NAN, 1.0),
        (0.0, f64::INFINITY),
        (1e308, 1e308),
    ] {
        assert_eq!(
            Binning::new(start, width, 2),
            Err(LandscapeError::InvalidBinning)
        );
    }
    assert_eq!(
        Binning::periodic(0.0, 0.0, 4),
        Err(LandscapeError::InvalidBinning)
    );
    assert_eq!(
        LandscapeError::InvalidBinning.to_string(),
        "The axis bins are not finite or have no width."
    );
    let binning = Binning::new(-1.0, 0.5, 4).unwrap();
    assert_eq!(binning.end(), 1.0);
    assert_eq!(binning.bin(1.0), Some(3));
    assert_eq!(binning.bin(-1.0), Some(0));
    assert_eq!(binning.bin(-0.5), Some(1));
    assert_eq!(binning.bin(1.5), None);
}

#[test]
fn two_axes_run_the_first_axis_fastest() -> Result<(), LandscapeError> {
    // Three bins on x, two on y; samples at (x, y) bins (2, 0), (0, 1),
    // (1, 1), and (1, 1).
    let x = [2.5, 0.5, 1.5, 1.2].map(Some);
    let y = [0.5, 1.5, 1.5, 1.9].map(Some);
    let axes = [
        GridAxis::new(&x, unit_bins(3), NANOMETER),
        GridAxis::new(&y, unit_bins(2), DEGREE),
    ];
    let surface = Landscape::free_energy(&every_sample(4), &axes, &at_300())?;
    assert_eq!(surface.shape(), [3, 2]);
    assert_eq!(surface.binnings(), [unit_bins(3), unit_bins(2)]);
    assert_eq!(surface.units(), [NANOMETER, DEGREE]);
    let counts = surface
        .cells()
        .iter()
        .map(|cell| cell.count())
        .collect::<Vec<_>>();
    assert_eq!(counts, [0, 0, 1, 1, 2, 0]);
    assert_eq!(surface.cell_at(&[1.5, 1.5]), Some(4));
    assert_eq!(surface.cell_at(&[2.5, 0.5]), Some(2));
    assert_eq!(surface.cell_at(&[2.5]), None);
    assert_eq!(surface.cell_at(&[3.5, 0.5]), None);
    assert_eq!(surface.cells()[4].representative(), Some(2));
    assert_eq!(surface.cells()[4].value(), Some(0.0));
    Ok(())
}

#[test]
fn means_weigh_samples_by_their_exact_shares() -> Result<(), LandscapeError> {
    // One cell of weights 1 and 3, mapped values 10 and 20.
    let values = [Some(0.5), Some(0.6), Some(0.7), Some(1.5)];
    let of = [Some(10.0), Some(20.0), None, Some(f64::NAN)];
    let axis = GridAxis::new(&values, unit_bins(2), NANOMETER);
    for estimator in ESTIMATORS {
        let ensemble = reweighted(&[0.0, 3.0_f64.ln(), 0.0, 0.0], estimator);
        let surface = Landscape::mean(
            ensemble.distribution(),
            &[axis],
            &of,
            DIMENSIONLESS,
            &options(),
        )?;
        assert_eq!(surface.quantity(), SurfaceQuantity::Mean);
        assert_eq!(surface.value_unit(), DIMENSIONLESS);
        assert_eq!(surface.temperature(), Some(kelvin(300.0)));
        // Samples without a mapped value are excluded.
        assert_eq!(surface.excluded(), 2);
        assert_eq!(surface.cells()[0].count(), 2);
        assert_close(surface.cells()[0].value().unwrap(), 17.5);
        assert_close(surface.cells()[0].effective_count(), 1.6);
        assert_eq!(surface.cells()[1].value(), None);
    }
    Ok(())
}

#[test]
fn samples_without_a_cell_are_counted() -> Result<(), LandscapeError> {
    // Geometry 4 is no sample, so its value, out of range, is never read.
    let values = [Some(0.5), None, Some(f64::NAN), Some(7.0), Some(9.0)];
    let distribution =
        Distribution::new(Samples::Subset(vec![0, 1, 2, 3]), Partition::one_state(4)).unwrap();
    let axis = GridAxis::new(&values, unit_bins(2), NANOMETER);
    let surface = Landscape::free_energy(&distribution, &[axis], &at_300())?;
    assert_eq!(surface.excluded(), 2);
    assert_eq!(surface.outside(), 1);
    assert_eq!(surface.cells()[0].count(), 1);
    assert_eq!(surface.probabilities(), [1.0, 0.0]);
    // A missing value on either axis excludes a sample, and a value outside
    // the bins on either axis puts it outside.
    let other = [None, Some(0.5), Some(0.5), Some(0.5)];
    let axes = [axis, GridAxis::new(&other, unit_bins(1), NANOMETER)];
    let surface = Landscape::free_energy(&distribution, &axes, &at_300())?;
    assert_eq!(surface.excluded(), 3);
    assert_eq!(surface.outside(), 1);
    assert!(surface.cells().iter().all(|cell| cell.count() == 0));
    Ok(())
}

#[test]
fn inputs_must_cover_the_samples() {
    let short = [Some(0.5), Some(0.5)];
    let full = [Some(0.5); 3];
    let distribution = every_sample(3);
    let options = at_300();
    let error = Landscape::free_energy(
        &distribution,
        &[
            GridAxis::new(&full, unit_bins(1), NANOMETER),
            GridAxis::new(&short, unit_bins(1), NANOMETER),
        ],
        &options,
    )
    .unwrap_err();
    assert_eq!(
        error,
        LandscapeError::ShortAxis {
            axis: 1,
            geometry: 2
        }
    );
    assert_eq!(error.to_string(), "Axis 1 has no value for geometry 2.");
    let axis = GridAxis::new(&full, unit_bins(1), NANOMETER);
    let error = Landscape::mean(&distribution, &[axis], &short, NANOMETER, &options).unwrap_err();
    assert_eq!(error, LandscapeError::ShortValues { geometry: 2 });
    assert_eq!(
        error.to_string(),
        "The mapped values have no value for geometry 2."
    );
    for axes in [&[][..], &[axis; 3][..]] {
        let error = Landscape::free_energy(&distribution, axes, &options).unwrap_err();
        assert_eq!(error, LandscapeError::AxisCount(axes.len()));
    }
    assert_eq!(
        LandscapeError::AxisCount(3).to_string(),
        "A landscape has one or two axes. There are 3 axes."
    );
    let wide = GridAxis::new(&full, unit_bins(5_000), NANOMETER);
    let error = Landscape::free_energy(&distribution, &[wide, wide], &options).unwrap_err();
    assert!(matches!(error, LandscapeError::TooManyCells(_)));
    assert_eq!(
        error.to_string(),
        "The landscape exceeds its grid cell limit: 25,000,000 cells of 4,194,304."
    );
}
