//! Energy surfaces: the samples of weights binned on one or two axes, with a
//! free energy or a mean value per cell.

use std::fmt;

use super::{kish, log_sum_exp, Estimator, Weights};
use crate::diagnostics::{LimitExceeded, Resource};
use crate::structure::{StructureSetError, ThermodynamicState};
use crate::units::{Quantity, Unit, UnitError, CANONICAL_ENERGY_UNIT};

/// Equal bins along one axis: `bins` bins from `start` to `end`, in the
/// unit of the values they bin.
///
/// Every binning spans two edges and divides them into equal bins; the
/// width follows ([`Self::width`]). [`Self::between`] bins values between
/// the edges, [`Self::periodic`] wraps every finite value into the period
/// the edges span, and [`Self::covering`] takes the edges from the values.
///
/// Bin `i` holds the values from `start + i·width` up to the next edge; the
/// last bin also holds `end`. On `periodic(-180.0, 180.0, 36)` the value 180
/// wraps to -180 and lands in bin 0.
///
/// ```
/// use kekule::statistics::Binning;
///
/// let dihedral = Binning::periodic(-180.0, 180.0, 36)?;
/// assert_eq!(dihedral.width(), 10.0);
/// assert_eq!(dihedral.bin(-180.0), Some(0));
/// assert_eq!(dihedral.bin(180.0), Some(0));
/// assert_eq!(dihedral.bin(-175.0), Some(0));
/// assert_eq!(dihedral.bin(179.0), Some(35));
/// assert_eq!(dihedral.center(0), Some(-175.0));
///
/// let distance = Binning::between(0.0, 2.0, 20)?;
/// assert_eq!(distance.bin(2.0), Some(19));
/// assert_eq!(distance.bin(2.5), None);
///
/// let covering = Binning::covering(&[Some(0.1), None, Some(0.4)], 3)?;
/// assert_eq!(covering.bin(0.1), Some(0));
/// assert_eq!(covering.bin(0.4), Some(2));
/// assert_eq!(covering.bin(0.5), None);
/// # Ok::<(), kekule::statistics::LandscapeError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Binning {
    start: f64,
    end: f64,
    width: f64,
    bins: usize,
    periodic: bool,
}

impl Binning {
    /// `bins` bins from `start` to `end`. A value equal to `end` lands in
    /// the last bin, however the width rounds; a value outside the edges
    /// lands in none.
    ///
    /// # Errors
    ///
    /// [`LandscapeError::NoBins`] for no bins, and
    /// [`LandscapeError::InvalidBinning`] when `start` or `end` is not
    /// finite, or `end` is not above `start` by a positive width.
    pub fn between(start: f64, end: f64, bins: usize) -> Result<Self, LandscapeError> {
        Self::spanning(start, end, bins, false)
    }

    /// `bins` bins dividing the period from `start` to `end`, such as an
    /// angle from -180° to 180°. Values outside the period are wrapped into
    /// it, so `end` lands in bin 0.
    ///
    /// # Errors
    ///
    /// The errors of [`Self::between`].
    pub fn periodic(start: f64, end: f64, bins: usize) -> Result<Self, LandscapeError> {
        Self::spanning(start, end, bins, true)
    }

    /// `bins` bins from the smallest to the largest finite value, so that
    /// every finite value lands in a bin, the largest in the last one. When
    /// every value is the same, the bins are one unit wide, and the value
    /// lies at the centre of bin `bins / 2`.
    ///
    /// # Errors
    ///
    /// [`LandscapeError::NoBins`] for no bins,
    /// [`LandscapeError::NoValues`] without a finite value, and
    /// [`LandscapeError::InvalidBinning`] when the range of the values is
    /// not finite or too narrow to divide.
    pub fn covering(values: &[Option<f64>], bins: usize) -> Result<Self, LandscapeError> {
        if bins == 0 {
            return Err(LandscapeError::NoBins);
        }
        let (lowest, highest) = values
            .iter()
            .flatten()
            .filter(|value| value.is_finite())
            .fold(
                (f64::INFINITY, f64::NEG_INFINITY),
                |(lowest, highest), &value| (lowest.min(value), highest.max(value)),
            );
        if lowest > highest {
            return Err(LandscapeError::NoValues);
        }
        if lowest == highest {
            let start = lowest - ((bins / 2) as f64 + 0.5);
            return Self::checked(start, start + bins as f64, 1.0, bins, false);
        }
        Self::between(lowest, highest, bins)
    }

    fn spanning(start: f64, end: f64, bins: usize, periodic: bool) -> Result<Self, LandscapeError> {
        if bins == 0 {
            return Err(LandscapeError::NoBins);
        }
        Self::checked(start, end, (end - start) / bins as f64, bins, periodic)
    }

    fn checked(
        start: f64,
        end: f64,
        width: f64,
        bins: usize,
        periodic: bool,
    ) -> Result<Self, LandscapeError> {
        let finite = start.is_finite() && end.is_finite() && width.is_finite();
        if !(finite && width > 0.0 && end > start) {
            return Err(LandscapeError::InvalidBinning);
        }
        Ok(Self {
            start,
            end,
            width,
            bins,
            periodic,
        })
    }

    /// The lower edge of the first bin.
    pub const fn start(&self) -> f64 {
        self.start
    }

    /// The upper edge of the last bin.
    pub const fn end(&self) -> f64 {
        self.end
    }

    /// The width of every bin.
    pub const fn width(&self) -> f64 {
        self.width
    }

    /// The number of bins.
    pub const fn bins(&self) -> usize {
        self.bins
    }

    /// Whether values wrap into the period.
    pub const fn is_periodic(&self) -> bool {
        self.periodic
    }

    /// The period of a periodic binning, `end - start`.
    pub fn period(&self) -> Option<f64> {
        self.periodic.then_some(self.end - self.start)
    }

    /// The bin of `value`, or `None` for a value that is not finite or, on
    /// a binning that is not periodic, lies outside the bins.
    pub fn bin(&self, value: f64) -> Option<usize> {
        self.position(value)
            .map(|position| (position as usize).min(self.bins - 1))
    }

    /// The centre of bin `index`, or `None` past the last bin.
    pub fn center(&self, index: usize) -> Option<f64> {
        (index < self.bins).then_some(self.start + (index as f64 + 0.5) * self.width)
    }

    /// Where `value` lies, in bin widths from the start, wrapped into the
    /// period in value space on a periodic binning; `None` as [`Self::bin`]
    /// gives it.
    fn position(&self, value: f64) -> Option<f64> {
        if !value.is_finite() {
            return None;
        }
        if self.periodic {
            let period = self.end - self.start;
            let wrapped = (value - self.start).rem_euclid(period);
            wrapped
                .is_finite()
                .then_some(wrapped / period * self.bins as f64)
        } else {
            (value >= self.start && value <= self.end).then_some((value - self.start) / self.width)
        }
    }
}

/// One axis of a landscape: a value per geometry with its unit, and the
/// bins they fall into.
///
/// The values are indexed by geometry of the weights' set, so a host bins
/// set columns, the order, or its own measurements alike. Only the samples
/// are read; a sample whose value is missing or not finite lies in no cell.
#[derive(Debug, Clone, Copy)]
pub struct GridAxis<'a> {
    values: &'a [Option<f64>],
    unit: Unit,
    binning: Binning,
}

impl<'a> GridAxis<'a> {
    /// An axis of `values` by geometry, binned by `binning` in their unit.
    pub const fn new(values: Quantity<&'a [Option<f64>]>, binning: Binning) -> Self {
        Self {
            values: *values.value(),
            unit: values.unit(),
            binning,
        }
    }

    /// The value of each geometry, with its unit.
    pub const fn values(&self) -> Quantity<&'a [Option<f64>]> {
        Quantity::new(self.values, self.unit)
    }

    /// The bins, in the unit of the values.
    pub const fn binning(&self) -> Binning {
        self.binning
    }
}

/// Options of a [`Landscape`].
///
/// Unless [`Self::with_fewest`] sets it, a free-energy cell needs
/// [`Self::DEFAULT_FEWEST`] samples for a value, and a mean needs one: a
/// mean over a single structure, such as a point of a scan, is exact. Cell
/// populations are exact unless [`Self::with_estimator`] sets another
/// [`Estimator`].
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct LandscapeOptions {
    fewest: Option<usize>,
    temperature: Option<Quantity<f64>>,
    estimator: Estimator,
}

impl LandscapeOptions {
    /// The fewest samples a free-energy cell needs for a value unless set,
    /// PyReweighting's cutoff.
    pub const DEFAULT_FEWEST: usize = 10;

    /// The default options.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns these options with the fewest samples a cell needs for a
    /// value, for a free energy and a mean alike. A cell with fewer is
    /// masked: it keeps its count, population, and representative, but has
    /// no value.
    #[must_use]
    pub const fn with_fewest(mut self, fewest: usize) -> Self {
        self.fewest = Some(fewest);
        self
    }

    /// Returns these options with the estimator of each cell's population
    /// from the log weights of its samples; by default the exact
    /// [`Estimator::Exponential`]. Shares within a cell stay exact.
    #[must_use]
    pub const fn with_estimator(mut self, estimator: Estimator) -> Self {
        self.estimator = estimator;
        self
    }

    /// Returns these options with the temperature of the free energy, for
    /// weights without one. Weights with a temperature take only their own.
    ///
    /// # Errors
    ///
    /// [`LandscapeError::Unit`] for a quantity that is not a temperature, and
    /// [`LandscapeError::InvalidTemperature`] for one that is not finite and
    /// positive.
    pub fn with_temperature(
        mut self,
        temperature: Option<Quantity<f64>>,
    ) -> Result<Self, LandscapeError> {
        self.temperature = temperature
            .map(|temperature| {
                ThermodynamicState::canonical(temperature)
                    .map(|state| state.temperature())
                    .map_err(LandscapeError::from_temperature)
            })
            .transpose()?;
        Ok(self)
    }

    /// The fewest samples a cell needs for a value, if set.
    pub const fn fewest(&self) -> Option<usize> {
        self.fewest
    }

    /// The temperature given for weights without one, in
    /// [`CANONICAL_TEMPERATURE_UNIT`](crate::units::CANONICAL_TEMPERATURE_UNIT).
    pub const fn temperature(&self) -> Option<Quantity<f64>> {
        self.temperature
    }

    /// The estimator of each cell's population.
    pub const fn estimator(&self) -> Estimator {
        self.estimator
    }
}

/// What a [`Landscape`]'s cell values are.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SurfaceQuantity {
    /// The free energy of the cell's population relative to the most
    /// populated cell that is not masked ([`Landscape::free_energy`]).
    FreeEnergy,
    /// The probability-weighted mean of a value over the cell's samples
    /// ([`Landscape::mean`]).
    Mean,
}

/// One cell of a [`Landscape`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cell {
    count: usize,
    effective_count: f64,
    log_population: Option<f64>,
    value: Option<f64>,
    representative: Option<usize>,
}

impl Cell {
    /// The number of samples in the cell.
    pub const fn count(&self) -> usize {
        self.count
    }

    /// Kish's effective sample size of the cell over the exact weights of
    /// its samples; `0.0` for an empty cell.
    pub const fn effective_count(&self) -> f64 {
        self.effective_count
    }

    /// The natural logarithm of the cell's probability among the cells:
    /// the options' estimator over the log weights of its samples,
    /// normalized so that the cells' probabilities sum to one. `None` for an
    /// empty cell or one whose estimate is not finite.
    pub const fn log_population(&self) -> Option<f64> {
        self.log_population
    }

    /// The cell's value in [`Landscape::value_unit`], or `None` for an
    /// empty cell, one with fewer samples than [`Landscape::fewest`], or one
    /// without a finite value.
    pub const fn value(&self) -> Option<f64> {
        self.value
    }

    /// The geometry of the sample nearest the cell centre, or `None` for an
    /// empty cell.
    pub const fn representative(&self) -> Option<usize> {
        self.representative
    }
}

/// The samples of [`Weights`] binned on one or two axes, with a value per
/// cell: a free-energy surface or a map of a mean value.
///
/// **Cells.** A sample lies in the cell of its axis values. A sample whose
/// value on some axis is missing or not finite is excluded
/// ([`Self::excluded`]), and one whose value lies outside a binning that is
/// not periodic is outside ([`Self::outside`]); neither is in any cell.
/// Geometries that are not samples are not read.
///
/// **Order.** Cells are stored with the first axis fastest: on axes of `nx`
/// and `ny` bins, cell `(i, j)` is at index `j·nx + i`, from the lowest
/// bins of both axes, so the rows of a plotted grid run from the lowest
/// second-axis value up.
///
/// **Populations.** A cell's population is the log-sum-exp of its samples'
/// log weights, or the options' [`Estimator`] over them
/// ([`LandscapeOptions::with_estimator`]), and the cells are normalized
/// among themselves: [`Self::probabilities`] is the weighted histogram.
/// With the exact estimator, the default, a cell that holds exactly some
/// samples, with every sample in a cell, has their probability
/// ([`Weights::log_population`]). Within a cell, each sample's share is
/// exact, `exp(ℓ_i) / Σ exp(ℓ_j)`; the effective count is Kish's over the
/// same weights.
///
/// **Values.** A cell with fewer samples than [`Self::fewest`] is masked:
/// it has no value, but keeps its count, population, and representative,
/// and its probability stays in [`Self::probabilities`]. Unless the options
/// set it, a free energy needs [`LandscapeOptions::DEFAULT_FEWEST`] samples
/// per cell and a mean one. The free energy of a cell is
/// `F = -RT ln(P_cell / P_max)` in kJ/mol, with `P_max` the population of
/// the most populated cell that is not masked, so the lowest value is 0. A
/// mean is the average of a value over the cell's samples, each with its
/// exact share. [`Self::lowest_cell`] names the cell with the lowest value.
///
/// **Representative.** A cell's representative is the sample nearest its
/// centre, measured in bin widths on each axis, across the period on a
/// periodic axis; of equally near ones, the one with the larger weight,
/// then the lowest geometry index. It does not depend on the estimator.
///
/// **Temperature.** The free energy is at the weights' temperature
/// ([`Weights::temperature`]). Weights without one take
/// [`LandscapeOptions::with_temperature`]; weights with one refuse a
/// different temperature.
///
/// ```
/// use kekule::statistics::{Binning, GridAxis, Landscape, LandscapeOptions, Weights};
/// use kekule::units::{Quantity, DEGREE, KELVIN};
///
/// // Four equally weighted samples: three near 0°, one near 90°, on bins
/// // centred at 0°, 90°, 180°, and -90°.
/// let weights = Weights::equal(4);
/// let angles = [Some(-10.0), Some(5.0), Some(10.0), Some(95.0)];
/// let quarters = Binning::periodic(-45.0, 315.0, 4)?;
/// let axis = GridAxis::new(Quantity::new(&angles, DEGREE), quarters);
/// let options = LandscapeOptions::new()
///     .with_fewest(1)
///     .with_temperature(Some(Quantity::new(300.0, KELVIN)))?;
/// let surface = Landscape::free_energy(&weights, &[axis], &options)?;
/// assert_eq!(surface.shape(), [4]);
/// assert_eq!(surface.cell_at(&[0.0]), Some(0));
/// assert_eq!(surface.cell_at(&[90.0]), Some(1));
/// assert_eq!(surface.lowest_cell(), Some(0));
/// let [well, shoulder, ..] = surface.cells() else { unreachable!() };
/// assert_eq!(well.value(), Some(0.0));
/// // RT ln 3 above the well.
/// let rt = 8.314_462_618_153_24e-3 * 300.0;
/// assert!((shoulder.value().unwrap() - rt * 3.0_f64.ln()).abs() < 1e-12);
/// let probabilities = surface.probabilities();
/// assert!((probabilities[0] - 0.75).abs() < 1e-15);
/// assert_eq!(probabilities[2..], [0.0, 0.0]);
/// // The sample at 5° is nearest the well's centre, 0°.
/// assert_eq!(well.representative(), Some(1));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct Landscape {
    quantity: SurfaceQuantity,
    binnings: Vec<Binning>,
    units: Vec<Unit>,
    shape: Vec<usize>,
    value_unit: Unit,
    temperature: Option<Quantity<f64>>,
    fewest: usize,
    cells: Vec<Cell>,
    lowest_cell: Option<usize>,
    excluded: usize,
    outside: usize,
}

impl Landscape {
    /// The most cells a landscape holds, a grid of 2048 × 2048.
    pub const MAX_CELLS: usize = 1 << 22;

    /// The free-energy surface of the samples of `weights` on `axes`, in
    /// [`CANONICAL_ENERGY_UNIT`]. A cell with fewer samples than
    /// [`LandscapeOptions::with_fewest`] sets, or
    /// [`LandscapeOptions::DEFAULT_FEWEST`] unless set, is masked.
    ///
    /// # Errors
    ///
    /// [`LandscapeError::AxisCount`] for other than one or two axes,
    /// [`LandscapeError::TooManyCells`] beyond [`Self::MAX_CELLS`],
    /// [`LandscapeError::ShortAxis`] when an axis has no entry for a sample,
    /// [`LandscapeError::MissingTemperature`] when neither the weights nor
    /// the options give a temperature, and
    /// [`LandscapeError::ConflictingTemperature`] when the options give
    /// another one than the weights.
    pub fn free_energy(
        weights: &Weights,
        axes: &[GridAxis<'_>],
        options: &LandscapeOptions,
    ) -> Result<Self, LandscapeError> {
        Self::build(
            weights,
            axes,
            None,
            CANONICAL_ENERGY_UNIT,
            SurfaceQuantity::FreeEnergy,
            options,
        )
    }

    /// The map of the probability-weighted mean of `of`, a value per
    /// geometry, over the samples of each cell of `weights` on `axes`,
    /// in the unit of `of`. A sample whose value of `of` is missing or not
    /// finite is excluded, as one missing an axis value is. A cell is masked
    /// only below the fewest samples [`LandscapeOptions::with_fewest`] sets:
    /// the mean of a single structure, such as a point of a scan, is exact.
    ///
    /// # Errors
    ///
    /// The errors of [`Self::free_energy`] but
    /// [`LandscapeError::MissingTemperature`], and
    /// [`LandscapeError::ShortValues`] when `of` has no entry for a sample.
    pub fn mean(
        weights: &Weights,
        axes: &[GridAxis<'_>],
        of: Quantity<&[Option<f64>]>,
        options: &LandscapeOptions,
    ) -> Result<Self, LandscapeError> {
        Self::build(
            weights,
            axes,
            Some(*of.value()),
            of.unit(),
            SurfaceQuantity::Mean,
            options,
        )
    }

    fn build(
        weights: &Weights,
        axes: &[GridAxis<'_>],
        of: Option<&[Option<f64>]>,
        value_unit: Unit,
        quantity: SurfaceQuantity,
        options: &LandscapeOptions,
    ) -> Result<Self, LandscapeError> {
        if !(1..=2).contains(&axes.len()) {
            return Err(LandscapeError::AxisCount(axes.len()));
        }
        let temperature = resolve_temperature(weights, options)?;
        let thermal_energy = match (quantity, temperature) {
            (SurfaceQuantity::FreeEnergy, None) => {
                return Err(LandscapeError::MissingTemperature);
            }
            (SurfaceQuantity::FreeEnergy, Some(temperature)) => Some(
                *ThermodynamicState::canonical(temperature)
                    .map_err(LandscapeError::from_temperature)?
                    .thermal_energy()
                    .value(),
            ),
            (SurfaceQuantity::Mean, _) => None,
        };
        let shape = axes
            .iter()
            .map(|axis| axis.binning.bins)
            .collect::<Vec<_>>();
        let cell_count = shape
            .iter()
            .try_fold(1_usize, |count, &bins| count.checked_mul(bins))
            .unwrap_or(usize::MAX);
        LimitExceeded::check(Resource::Cells, Self::MAX_CELLS, cell_count)
            .map_err(LandscapeError::TooManyCells)?;

        let samples = weights.samples().collect::<Vec<_>>();
        let log_weights = samples
            .iter()
            .map(|&geometry| weights.log_weight(geometry).expect("a sample"))
            .collect::<Vec<_>>();
        let mut cell_of = Vec::with_capacity(samples.len());
        let (mut excluded, mut outside) = (0, 0);
        for &geometry in &samples {
            match sample_cell(axes, of, geometry)? {
                Placement::Cell(cell) => cell_of.push(Some(cell)),
                Placement::Excluded => {
                    excluded += 1;
                    cell_of.push(None);
                }
                Placement::Outside => {
                    outside += 1;
                    cell_of.push(None);
                }
            }
        }
        // Shifted by the largest log weight in a cell, cell populations keep
        // their digits under a large common offset, and a heavier sample in
        // no cell takes none.
        let shift = cell_of
            .iter()
            .zip(&log_weights)
            .filter(|(cell, _)| cell.is_some())
            .map(|(_, &log_weight)| log_weight)
            .reduce(f64::max)
            .unwrap_or(0.0);
        let grouping = Grouping::new(&cell_of, cell_count);
        let estimates = (0..cell_count)
            .map(|cell| {
                let values = grouping
                    .members(cell)
                    .iter()
                    .map(|&position| log_weights[position] - shift)
                    .collect::<Vec<_>>();
                options.estimator.log_sum(&values)
            })
            .collect::<Vec<_>>();
        let total = log_sum_exp(estimates.iter().flatten().copied());
        let log_populations = estimates
            .into_iter()
            .map(|estimate| estimate.map(|estimate| estimate - total))
            .collect::<Vec<_>>();

        let fewest = options.fewest.unwrap_or(match quantity {
            SurfaceQuantity::FreeEnergy => LandscapeOptions::DEFAULT_FEWEST,
            SurfaceQuantity::Mean => 1,
        });
        let reference = (0..cell_count)
            .filter(|&cell| grouping.members(cell).len() >= fewest)
            .filter_map(|cell| log_populations[cell])
            .fold(None, |highest: Option<f64>, value| {
                Some(highest.map_or(value, |highest| highest.max(value)))
            });
        let members = Members {
            samples: &samples,
            log_weights: &log_weights,
        };
        let cells = (0..cell_count)
            .map(|cell| {
                let positions = grouping.members(cell);
                let count = positions.len();
                let log_population = log_populations[cell];
                let value = if count == 0 || count < fewest {
                    None
                } else {
                    match (quantity, thermal_energy, of) {
                        (SurfaceQuantity::FreeEnergy, Some(thermal_energy), _) => log_population
                            .zip(reference)
                            .map(|(population, reference)| {
                                (reference - population) * thermal_energy
                            }),
                        (SurfaceQuantity::Mean, _, Some(of)) => {
                            members.weighted_mean(positions, of)
                        }
                        _ => None,
                    }
                };
                Cell {
                    count,
                    effective_count: if count == 0 {
                        0.0
                    } else {
                        kish(positions.iter().map(|&position| log_weights[position]))
                    },
                    log_population,
                    value: value.filter(|value| value.is_finite()),
                    representative: members.representative(axes, cell, positions),
                }
            })
            .collect::<Vec<_>>();
        let lowest_cell = find_lowest_cell(&cells);
        Ok(Self {
            quantity,
            binnings: axes.iter().map(|axis| axis.binning).collect(),
            units: axes.iter().map(|axis| axis.unit).collect(),
            shape,
            value_unit,
            temperature,
            fewest,
            cells,
            lowest_cell,
            excluded,
            outside,
        })
    }

    /// What the cell values are.
    pub const fn quantity(&self) -> SurfaceQuantity {
        self.quantity
    }

    /// The cells, the first axis fastest: cell `(i, j)` is at `j·nx + i`.
    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    /// The index of the cell with the lowest value, or `None` when no cell
    /// has a value. Of equal values, the first in cell order. On a free
    /// energy, it is the most populated cell with a value, the one at 0.
    pub const fn lowest_cell(&self) -> Option<usize> {
        self.lowest_cell
    }

    /// The number of bins on each axis.
    pub fn shape(&self) -> &[usize] {
        &self.shape
    }

    /// The bins of each axis.
    pub fn binnings(&self) -> &[Binning] {
        &self.binnings
    }

    /// The unit of each axis.
    pub fn units(&self) -> &[Unit] {
        &self.units
    }

    /// The unit of the cell values: [`CANONICAL_ENERGY_UNIT`] for a free
    /// energy, the mapped value's own for a mean.
    pub const fn value_unit(&self) -> Unit {
        self.value_unit
    }

    /// The weights' temperature, or the one the options give for weights
    /// without one, in
    /// [`CANONICAL_TEMPERATURE_UNIT`](crate::units::CANONICAL_TEMPERATURE_UNIT).
    /// `None` when neither gives one, which only a mean allows.
    pub const fn temperature(&self) -> Option<Quantity<f64>> {
        self.temperature
    }

    /// The fewest samples a cell needed for a value: the options' own, or
    /// the default for the quantity.
    pub const fn fewest(&self) -> usize {
        self.fewest
    }

    /// The number of samples in no cell for a missing or non-finite value.
    pub const fn excluded(&self) -> usize {
        self.excluded
    }

    /// The number of samples in no cell for a value outside the bins.
    pub const fn outside(&self) -> usize {
        self.outside
    }

    /// The probability of each cell, in cell order: the normalized cell
    /// populations, the weighted histogram of the samples. Empty cells have
    /// probability `0.0`; masked cells keep theirs.
    pub fn probabilities(&self) -> Vec<f64> {
        self.cells
            .iter()
            .map(|cell| cell.log_population.map_or(0.0, f64::exp))
            .collect()
    }

    /// The index of the cell holding `coordinates`, one value per axis in
    /// its unit, or `None` for another number of values or a value in no
    /// bin.
    pub fn cell_at(&self, coordinates: &[f64]) -> Option<usize> {
        if coordinates.len() != self.binnings.len() {
            return None;
        }
        let mut index = 0;
        let mut stride = 1;
        for (binning, &value) in self.binnings.iter().zip(coordinates) {
            index += binning.bin(value)? * stride;
            stride *= binning.bins;
        }
        Some(index)
    }

    /// The centre of cell `index`, one value per axis in its unit, or
    /// `None` past the last cell: the inverse of [`Self::cell_at`].
    pub fn center(&self, index: usize) -> Option<Vec<f64>> {
        if index >= self.cells.len() {
            return None;
        }
        let mut remainder = index;
        self.binnings
            .iter()
            .map(|binning| {
                let bin = remainder % binning.bins;
                remainder /= binning.bins;
                binning.center(bin)
            })
            .collect()
    }
}

/// Where a sample falls.
enum Placement {
    Cell(usize),
    Excluded,
    Outside,
}

/// The cell of the sample at `geometry`.
fn sample_cell(
    axes: &[GridAxis<'_>],
    of: Option<&[Option<f64>]>,
    geometry: usize,
) -> Result<Placement, LandscapeError> {
    let mut cell = Some(0);
    let mut stride = 1;
    let mut missing = false;
    for (axis_index, axis) in axes.iter().enumerate() {
        let value = axis.values.get(geometry).ok_or(LandscapeError::ShortAxis {
            axis: axis_index,
            geometry,
        })?;
        match value.filter(|value| value.is_finite()) {
            None => missing = true,
            Some(value) => {
                cell = cell
                    .zip(axis.binning.bin(value))
                    .map(|(cell, bin)| cell + bin * stride);
            }
        }
        stride *= axis.binning.bins;
    }
    if let Some(of) = of {
        let value = of
            .get(geometry)
            .ok_or(LandscapeError::ShortValues { geometry })?;
        missing |= !value.is_some_and(f64::is_finite);
    }
    Ok(match (missing, cell) {
        (true, _) => Placement::Excluded,
        (false, None) => Placement::Outside,
        (false, Some(cell)) => Placement::Cell(cell),
    })
}

/// The index of the cell with the lowest value, the first of equal ones.
fn find_lowest_cell(cells: &[Cell]) -> Option<usize> {
    let mut lowest: Option<(usize, f64)> = None;
    for (index, cell) in cells.iter().enumerate() {
        if let Some(value) = cell.value {
            if lowest.is_none_or(|(_, least)| value < least) {
                lowest = Some((index, value));
            }
        }
    }
    lowest.map(|(index, _)| index)
}

/// The temperature of a landscape over `weights`: theirs, or the options'
/// when they have none.
fn resolve_temperature(
    weights: &Weights,
    options: &LandscapeOptions,
) -> Result<Option<Quantity<f64>>, LandscapeError> {
    let Some(own) = weights.temperature() else {
        return Ok(options.temperature);
    };
    if let Some(given) = options.temperature {
        let given = *given.value();
        let own = *own.value();
        if (given - own).abs() > 1e-12 * own {
            return Err(LandscapeError::ConflictingTemperature {
                given,
                weights: own,
            });
        }
    }
    Ok(Some(own))
}

/// The sample positions of each cell, stored flat: every cell's members,
/// cell after cell, by increasing position.
struct Grouping {
    /// Where each cell's members start in `members`, and one past the last.
    starts: Vec<usize>,
    members: Vec<usize>,
}

impl Grouping {
    /// The cells of `cells` that `cell_of` puts each sample position in,
    /// leaving out the samples in no cell; every cell is below `cells`.
    fn new(cell_of: &[Option<usize>], cells: usize) -> Self {
        let mut starts = vec![0_usize; cells + 1];
        for cell in cell_of.iter().flatten() {
            starts[cell + 1] += 1;
        }
        for cell in 0..cells {
            starts[cell + 1] += starts[cell];
        }
        let mut members = vec![0_usize; starts[cells]];
        let mut next = starts[..cells].to_vec();
        for (position, cell) in cell_of.iter().enumerate() {
            if let Some(cell) = cell {
                members[next[*cell]] = position;
                next[*cell] += 1;
            }
        }
        Self { starts, members }
    }

    /// The sample positions of cell `cell`, by increasing position.
    fn members(&self, cell: usize) -> &[usize] {
        &self.members[self.starts[cell]..self.starts[cell + 1]]
    }
}

/// The samples of a landscape, by position: each one's geometry and log
/// weight.
struct Members<'a> {
    samples: &'a [usize],
    log_weights: &'a [f64],
}

impl Members<'_> {
    /// The mean of `of` over the samples at `positions`, each with its exact
    /// share among them.
    fn weighted_mean(&self, positions: &[usize], of: &[Option<f64>]) -> Option<f64> {
        let highest = positions
            .iter()
            .map(|&position| self.log_weights[position])
            .fold(f64::NEG_INFINITY, f64::max);
        let (mut total, mut mean) = (0.0, 0.0);
        for &position in positions {
            let value = of[self.samples[position]]?;
            let weight = (self.log_weights[position] - highest).exp();
            if weight == 0.0 {
                continue;
            }
            total += weight;
            mean += weight / total * (value - mean);
        }
        (total > 0.0).then_some(mean)
    }

    /// The sample of `positions` nearest the centre of `cell`, as
    /// [`Landscape`] chooses it.
    fn representative(
        &self,
        axes: &[GridAxis<'_>],
        cell: usize,
        positions: &[usize],
    ) -> Option<usize> {
        let mut best: Option<(f64, f64, usize)> = None;
        for &position in positions {
            let geometry = self.samples[position];
            let mut remainder = cell;
            let mut distance = 0.0;
            for axis in axes {
                let bin = remainder % axis.binning.bins;
                remainder /= axis.binning.bins;
                let offset = axis.values[geometry]
                    .and_then(|value| axis.binning.position(value))
                    .map_or(f64::INFINITY, |at| at - (bin as f64 + 0.5));
                distance += offset * offset;
            }
            let weight = self.log_weights[position];
            let nearer = best.is_none_or(|(best_distance, best_weight, _)| {
                distance < best_distance || (distance == best_distance && weight > best_weight)
            });
            if nearer {
                best = Some((distance, weight, geometry));
            }
        }
        best.map(|(_, _, geometry)| geometry)
    }
}

/// Why a [`Landscape`] or a [`Binning`] cannot be built.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum LandscapeError {
    /// A binning has no bins.
    NoBins,
    /// A binning's start or end is not finite, or the end is not above the
    /// start by a positive bin width.
    InvalidBinning,
    /// [`Binning::covering`] found no finite value.
    NoValues,
    /// A landscape has one or two axes; this many were given.
    AxisCount(usize),
    /// The values of an axis end before a sample's geometry. Give an axis a
    /// value, or `None`, for every geometry of the set.
    ShortAxis {
        /// The index of the axis.
        axis: usize,
        /// The sample's geometry.
        geometry: usize,
    },
    /// The values of [`Landscape::mean`] end before a sample's geometry.
    ShortValues {
        /// The sample's geometry.
        geometry: usize,
    },
    /// A free energy needs a temperature: the weights have none and the
    /// options give none ([`LandscapeOptions::with_temperature`]).
    MissingTemperature,
    /// The options give another temperature than the weights have.
    ConflictingTemperature {
        /// The given temperature, in kelvin.
        given: f64,
        /// The weights' temperature, in kelvin.
        weights: f64,
    },
    /// A temperature is not finite and positive.
    InvalidTemperature,
    /// A temperature is not a temperature; see the source.
    Unit(UnitError),
    /// The axes have more cells than [`Landscape::MAX_CELLS`].
    TooManyCells(LimitExceeded),
}

impl LandscapeError {
    /// The error of a temperature that [`ThermodynamicState`] refuses.
    fn from_temperature(error: StructureSetError) -> Self {
        match error {
            StructureSetError::Unit(error) => Self::Unit(error),
            _ => Self::InvalidTemperature,
        }
    }
}

impl fmt::Display for LandscapeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoBins => f.write_str("The axis has no bins."),
            Self::InvalidBinning => f.write_str("The axis bins are not finite or have no width."),
            Self::NoValues => f.write_str("There are no values to bin."),
            Self::AxisCount(axes) => write!(
                f,
                "A landscape has one or two axes. There are {}.",
                crate::diagnostics::counted(*axes, "axis", "axes")
            ),
            Self::ShortAxis { axis, geometry } => {
                write!(f, "Axis {axis} has no value for geometry {geometry}.")
            }
            Self::ShortValues { geometry } => {
                write!(
                    f,
                    "The mapped values have no value for geometry {geometry}."
                )
            }
            Self::MissingTemperature => {
                f.write_str("The free energy has no temperature. Give one.")
            }
            Self::ConflictingTemperature { given, weights } => write!(
                f,
                "The temperature {given} K differs from the temperature of the weights, {weights} K."
            ),
            Self::InvalidTemperature => {
                f.write_str("The landscape temperature is not positive or not finite.")
            }
            Self::Unit(_) => f.write_str("The landscape temperature has the wrong unit."),
            Self::TooManyCells(limit) => {
                write!(f, "The landscape exceeds its {}.", limit.describe())
            }
        }
    }
}

impl std::error::Error for LandscapeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}
