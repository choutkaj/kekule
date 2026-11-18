//! Energy surfaces: the samples of a distribution binned on one or two
//! axes, with a free energy or a mean value per cell.

use std::fmt;

use crate::diagnostics::{LimitExceeded, Resource};
use crate::structure::{Distribution, StructureSetError, ThermodynamicState};
use crate::units::{Quantity, Unit, UnitError, CANONICAL_ENERGY_UNIT};

/// Equal bins along one axis: `bins` bins of `width` from `start`, in the
/// unit of the values they bin.
///
/// Bin `i` holds the values from `start + i·width` up to the next edge; the
/// last bin also holds its upper edge. A periodic binning wraps every finite
/// value into its period first, so on `periodic(-180.0, 360.0, 36)` the
/// value 180 lands in bin 0 with -180.
///
/// ```
/// use kekule::statistics::Binning;
///
/// let dihedral = Binning::periodic(-180.0, 360.0, 36)?;
/// assert_eq!(dihedral.bin(-180.0), Some(0));
/// assert_eq!(dihedral.bin(180.0), Some(0));
/// assert_eq!(dihedral.bin(-175.0), Some(0));
/// assert_eq!(dihedral.bin(179.0), Some(35));
/// assert_eq!(dihedral.center(0), Some(-175.0));
///
/// let distance = Binning::covering(&[Some(0.1), None, Some(0.4)], 3)?;
/// assert_eq!(distance.bin(0.1), Some(0));
/// assert_eq!(distance.bin(0.4), Some(2));
/// assert_eq!(distance.bin(0.5), None);
/// # Ok::<(), kekule::statistics::LandscapeError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Binning {
    start: f64,
    width: f64,
    bins: usize,
    periodic: bool,
}

impl Binning {
    /// `bins` bins of `width` from `start`.
    ///
    /// # Errors
    ///
    /// [`LandscapeError::NoBins`] for no bins, and
    /// [`LandscapeError::InvalidBinning`] when `start` or `width` is not
    /// finite, `width` is not positive, or the last edge is not finite.
    pub fn new(start: f64, width: f64, bins: usize) -> Result<Self, LandscapeError> {
        Self::checked(start, width, bins, false)
    }

    /// `bins` bins dividing the period `period` from `start`, such as an
    /// angle from -180° over 360°. Values outside the period are wrapped
    /// into it.
    ///
    /// # Errors
    ///
    /// The errors of [`Self::new`] for the width `period / bins`.
    pub fn periodic(start: f64, period: f64, bins: usize) -> Result<Self, LandscapeError> {
        if bins == 0 {
            return Err(LandscapeError::NoBins);
        }
        Self::checked(start, period / bins as f64, bins, true)
    }

    /// `bins` bins from the smallest to the largest finite value, so that
    /// every finite value lands in a bin, the largest in the last one even
    /// when the width is not exact in binary. When every value is the same,
    /// the bins are one unit wide around it.
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
        let count = bins as f64;
        if lowest == highest {
            return Self::new(lowest - count / 2.0, 1.0, bins);
        }
        let span = highest - lowest;
        let mut width = span / count;
        if !(span.is_finite() && width > 0.0) {
            return Err(LandscapeError::InvalidBinning);
        }
        // `bin` divides by the width as here, so the largest value's
        // position must not pass the last edge.
        while span / width > count {
            width = width.next_up();
        }
        Self::new(lowest, width, bins)
    }

    fn checked(
        start: f64,
        width: f64,
        bins: usize,
        periodic: bool,
    ) -> Result<Self, LandscapeError> {
        if bins == 0 {
            return Err(LandscapeError::NoBins);
        }
        let end = start + width * bins as f64;
        if !(start.is_finite() && width.is_finite() && width > 0.0 && end.is_finite()) {
            return Err(LandscapeError::InvalidBinning);
        }
        Ok(Self {
            start,
            width,
            bins,
            periodic,
        })
    }

    /// The lower edge of the first bin.
    pub const fn start(&self) -> f64 {
        self.start
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

    /// The upper edge of the last bin, `start + bins·width`.
    pub fn end(&self) -> f64 {
        self.start + self.width * self.bins as f64
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

    /// Where `value` lies, in bin widths from the start, wrapped into
    /// `0..=bins` on a periodic binning; `None` as [`Self::bin`] gives it.
    fn position(&self, value: f64) -> Option<f64> {
        if !value.is_finite() {
            return None;
        }
        let position = (value - self.start) / self.width;
        let bins = self.bins as f64;
        if self.periodic {
            let wrapped = position.rem_euclid(bins);
            wrapped.is_finite().then_some(wrapped)
        } else {
            (value >= self.start && position <= bins).then_some(position)
        }
    }
}

/// One axis of a landscape: a value per geometry and the bins they fall
/// into, in `unit`.
///
/// The values are indexed by geometry of the distribution's set, so a host
/// bins set columns, the order, or its own measurements alike. Only the
/// distribution's samples are read; a sample whose value is missing or not
/// finite lies in no cell.
#[derive(Debug, Clone, Copy)]
pub struct GridAxis<'a> {
    values: &'a [Option<f64>],
    binning: Binning,
    unit: Unit,
}

impl<'a> GridAxis<'a> {
    /// An axis of `values` by geometry, binned by `binning`, in `unit`.
    pub const fn new(values: &'a [Option<f64>], binning: Binning, unit: Unit) -> Self {
        Self {
            values,
            binning,
            unit,
        }
    }

    /// The value of each geometry.
    pub const fn values(&self) -> &'a [Option<f64>] {
        self.values
    }

    /// The bins.
    pub const fn binning(&self) -> Binning {
        self.binning
    }

    /// The unit of the values and the bins.
    pub const fn unit(&self) -> Unit {
        self.unit
    }
}

/// Options of a [`Landscape`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LandscapeOptions {
    fewest: usize,
    temperature: Option<Quantity<f64>>,
}

impl Default for LandscapeOptions {
    /// At least [`Self::DEFAULT_FEWEST`] samples per cell, at the
    /// distribution's target temperature.
    fn default() -> Self {
        Self {
            fewest: Self::DEFAULT_FEWEST,
            temperature: None,
        }
    }
}

impl LandscapeOptions {
    /// The default fewest samples a cell needs for a value, PyReweighting's
    /// cutoff.
    pub const DEFAULT_FEWEST: usize = 10;

    /// The default options.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns these options with the fewest samples a cell needs for a
    /// value. A cell with fewer is masked: it keeps its count, population,
    /// and representative, but has no value.
    #[must_use]
    pub const fn with_fewest(mut self, fewest: usize) -> Self {
        self.fewest = fewest;
        self
    }

    /// Returns these options with the temperature of the free energy, for a
    /// distribution without a target. A distribution with a target takes
    /// only its own temperature.
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

    /// The fewest samples a cell needs for a value.
    pub const fn fewest(&self) -> usize {
        self.fewest
    }

    /// The temperature given for a distribution without a target, in
    /// [`CANONICAL_TEMPERATURE_UNIT`](crate::units::CANONICAL_TEMPERATURE_UNIT).
    pub const fn temperature(&self) -> Option<Quantity<f64>> {
        self.temperature
    }
}

/// What a [`Landscape`]'s cell values are.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SurfaceQuantity {
    /// The free energy of the cell's population relative to the most
    /// populated cell ([`Landscape::free_energy`]).
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
    /// the distribution's estimator over the exact log weights of its
    /// samples, normalized so that the cells' probabilities sum to one.
    /// `None` for an empty cell or one whose estimate is not finite.
    pub const fn log_population(&self) -> Option<f64> {
        self.log_population
    }

    /// The cell's value in [`Landscape::value_unit`], or `None` for an
    /// empty cell, one with fewer samples than
    /// [`LandscapeOptions::fewest`], or one without a finite value.
    pub const fn value(&self) -> Option<f64> {
        self.value
    }

    /// The geometry of the sample nearest the cell centre, or `None` for an
    /// empty cell.
    pub const fn representative(&self) -> Option<usize> {
        self.representative
    }
}

/// The samples of a [`Distribution`] binned on one or two axes, with a
/// value per cell: a free-energy surface or a map of a mean value.
///
/// **Cells.** A sample lies in the cell of its axis values. A sample whose
/// value on some axis is missing or not finite is excluded
/// ([`Self::excluded`]), and one whose value lies outside a binning that is
/// not periodic is outside ([`Self::outside`]); neither is in any cell.
/// Geometries that are not samples of the distribution are not read.
///
/// **Order.** Cells are stored with the first axis fastest: on axes of `nx`
/// and `ny` bins, cell `(i, j)` is at index `j·nx + i`, from the lowest
/// bins of both axes, so the rows of a plotted grid run from the lowest
/// second-axis value up.
///
/// **Populations.** A cell's population is the distribution's
/// [`estimator`](Distribution::estimator) over the exact log weights of its
/// samples, as [`Distribution::group_log_populations`] gives it, and the
/// cells are normalized among themselves: [`Self::probabilities`] is the
/// weighted histogram. A cell that holds exactly a state's samples, with
/// every sample in a cell, has the state's probability. Within a cell, each
/// sample's share is exact, `exp(ℓ_i) / Σ exp(ℓ_j)`; the effective count is
/// Kish's over the same weights.
///
/// **Values.** A cell with fewer samples than
/// [`LandscapeOptions::fewest`] is masked: it has no value, but keeps its
/// count, population, and representative, and its probability stays in
/// [`Self::probabilities`]. The free energy of a cell is
/// `F = -RT ln(P_cell / P_max)` in kJ/mol, with `P_max` the population of
/// the most populated cell that is not masked, so the lowest value is 0. A
/// mean is the average of a value over the cell's samples, each with its
/// exact share.
///
/// **Representative.** A cell's representative is the sample nearest its
/// centre, measured in bin widths on each axis, across the period on a
/// periodic axis; of equally near ones, the one with the larger exact
/// weight, then the lowest geometry index. It does not depend on the
/// estimator.
///
/// **Temperature.** The free energy is at the distribution's target
/// temperature. A distribution without a target takes
/// [`LandscapeOptions::with_temperature`]; one with a target refuses a
/// different temperature.
///
/// ```
/// use kekule::statistics::{Binning, GridAxis, Landscape, LandscapeOptions};
/// use kekule::structure::{Distribution, Partition, Samples};
/// use kekule::units::{Quantity, DEGREE, KELVIN};
///
/// // Four equally weighted samples: three near 0°, one near 90°, on bins
/// // centred at 0°, 90°, 180°, and -90°.
/// let distribution = Distribution::new(Samples::All, Partition::one_state(4))?;
/// let angles = [Some(-10.0), Some(5.0), Some(10.0), Some(95.0)];
/// let axis = GridAxis::new(&angles, Binning::periodic(-45.0, 360.0, 4)?, DEGREE);
/// let options = LandscapeOptions::new()
///     .with_fewest(1)
///     .with_temperature(Some(Quantity::new(300.0, KELVIN)))?;
/// let surface = Landscape::free_energy(&distribution, &[axis], &options)?;
/// assert_eq!(surface.shape(), [4]);
/// assert_eq!(surface.cell_at(&[0.0]), Some(0));
/// assert_eq!(surface.cell_at(&[90.0]), Some(1));
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
    excluded: usize,
    outside: usize,
}

impl Landscape {
    /// The most cells a landscape holds.
    pub const MAX_CELLS: usize = 1 << 24;

    /// The free-energy surface of `distribution` on `axes`, in
    /// [`CANONICAL_ENERGY_UNIT`].
    ///
    /// # Errors
    ///
    /// [`LandscapeError::AxisCount`] for other than one or two axes,
    /// [`LandscapeError::TooManyCells`] beyond [`Self::MAX_CELLS`],
    /// [`LandscapeError::ShortAxis`] when an axis has no entry for a sample,
    /// [`LandscapeError::MissingTemperature`] when neither the distribution
    /// nor the options give a temperature, and
    /// [`LandscapeError::ConflictingTemperature`] when the options give
    /// another one than the distribution's target.
    pub fn free_energy(
        distribution: &Distribution,
        axes: &[GridAxis<'_>],
        options: &LandscapeOptions,
    ) -> Result<Self, LandscapeError> {
        Self::build(
            distribution,
            axes,
            None,
            CANONICAL_ENERGY_UNIT,
            SurfaceQuantity::FreeEnergy,
            options,
        )
    }

    /// The map of the probability-weighted mean of `of`, a value per
    /// geometry in `unit`, over the samples of each cell of `distribution`
    /// on `axes`. A sample whose value of `of` is missing or not finite is
    /// excluded, as one missing an axis value is.
    ///
    /// # Errors
    ///
    /// The errors of [`Self::free_energy`] but
    /// [`LandscapeError::MissingTemperature`], and
    /// [`LandscapeError::ShortValues`] when `of` has no entry for a sample.
    pub fn mean(
        distribution: &Distribution,
        axes: &[GridAxis<'_>],
        of: &[Option<f64>],
        unit: Unit,
        options: &LandscapeOptions,
    ) -> Result<Self, LandscapeError> {
        Self::build(
            distribution,
            axes,
            Some(of),
            unit,
            SurfaceQuantity::Mean,
            options,
        )
    }

    fn build(
        distribution: &Distribution,
        axes: &[GridAxis<'_>],
        of: Option<&[Option<f64>]>,
        value_unit: Unit,
        quantity: SurfaceQuantity,
        options: &LandscapeOptions,
    ) -> Result<Self, LandscapeError> {
        if !(1..=2).contains(&axes.len()) {
            return Err(LandscapeError::AxisCount(axes.len()));
        }
        let temperature = resolve_temperature(distribution, options)?;
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

        let samples = distribution.sample_count();
        let mut group_of = Vec::with_capacity(samples);
        let (mut excluded, mut outside) = (0, 0);
        for position in 0..samples {
            let geometry = distribution.sample(position);
            let cell = sample_cell(axes, of, geometry)?;
            match cell {
                Placement::Cell(cell) => {
                    group_of.push(Some(
                        u32::try_from(cell).expect("cells are bounded below u32::MAX"),
                    ));
                }
                Placement::Excluded => {
                    excluded += 1;
                    group_of.push(None);
                }
                Placement::Outside => {
                    outside += 1;
                    group_of.push(None);
                }
            }
        }
        let log_populations = distribution
            .group_log_populations(&group_of, cell_count)
            .expect("every sample has a group below the cell count");

        // The samples of each cell, by sample position, cell after cell.
        let mut starts = vec![0_usize; cell_count + 1];
        for cell in group_of.iter().flatten() {
            starts[*cell as usize + 1] += 1;
        }
        for cell in 0..cell_count {
            starts[cell + 1] += starts[cell];
        }
        let mut members = vec![0_usize; starts[cell_count]];
        let mut next = starts.clone();
        for (position, cell) in group_of.iter().enumerate() {
            if let Some(cell) = cell {
                members[next[*cell as usize]] = position;
                next[*cell as usize] += 1;
            }
        }

        let log_weights = distribution.exact_log_weights();
        let fewest = options.fewest;
        let reference = (0..cell_count)
            .filter(|&cell| starts[cell + 1] - starts[cell] >= fewest)
            .filter_map(|cell| log_populations[cell])
            .fold(None, |highest: Option<f64>, value| {
                Some(highest.map_or(value, |highest| highest.max(value)))
            });
        let cells = (0..cell_count)
            .map(|cell| {
                let positions = &members[starts[cell]..starts[cell + 1]];
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
                            weighted_mean(distribution, &log_weights, positions, of)
                        }
                        _ => None,
                    }
                };
                Cell {
                    count,
                    effective_count: distribution
                        .group_effective_sample_size(positions)
                        .unwrap_or(0.0),
                    log_population,
                    value: value.filter(|value| value.is_finite()),
                    representative: representative(
                        distribution,
                        &log_weights,
                        axes,
                        cell,
                        positions,
                    ),
                }
            })
            .collect();
        Ok(Self {
            quantity,
            binnings: axes.iter().map(|axis| axis.binning).collect(),
            units: axes.iter().map(|axis| axis.unit).collect(),
            shape,
            value_unit,
            temperature,
            fewest,
            cells,
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

    /// The temperature of the free energy, or of the distribution's target,
    /// in [`CANONICAL_TEMPERATURE_UNIT`](crate::units::CANONICAL_TEMPERATURE_UNIT); `None` for a mean over a
    /// distribution without a target.
    pub const fn temperature(&self) -> Option<Quantity<f64>> {
        self.temperature
    }

    /// The fewest samples a cell needed for a value.
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

/// The temperature of a landscape over `distribution`: its target's, or
/// the options' when it has none.
fn resolve_temperature(
    distribution: &Distribution,
    options: &LandscapeOptions,
) -> Result<Option<Quantity<f64>>, LandscapeError> {
    let Some(target) = distribution.temperature() else {
        return Ok(options.temperature);
    };
    if let Some(given) = options.temperature {
        let given = *given.value();
        let target = *target.value();
        if (given - target).abs() > 1e-12 * target {
            return Err(LandscapeError::ConflictingTemperature { given, target });
        }
    }
    Ok(Some(target))
}

/// The mean of `of` over the samples at `positions`, each with its exact
/// share among them.
fn weighted_mean(
    distribution: &Distribution,
    log_weights: &[f64],
    positions: &[usize],
    of: &[Option<f64>],
) -> Option<f64> {
    let highest = positions
        .iter()
        .map(|&position| log_weights[position])
        .fold(f64::NEG_INFINITY, f64::max);
    let (mut total, mut mean) = (0.0, 0.0);
    for &position in positions {
        let value = of[distribution.sample(position)]?;
        let weight = (log_weights[position] - highest).exp();
        total += weight;
        mean += weight / total * (value - mean);
    }
    (total > 0.0).then_some(mean)
}

/// The sample of `positions` nearest the centre of `cell`, as
/// [`Landscape`] chooses it.
fn representative(
    distribution: &Distribution,
    log_weights: &[f64],
    axes: &[GridAxis<'_>],
    cell: usize,
    positions: &[usize],
) -> Option<usize> {
    let mut best: Option<(f64, f64, usize)> = None;
    for &position in positions {
        let geometry = distribution.sample(position);
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
        let weight = log_weights[position];
        let nearer = best.is_none_or(|(best_distance, best_weight, _)| {
            distance < best_distance || (distance == best_distance && weight > best_weight)
        });
        if nearer {
            best = Some((distance, weight, geometry));
        }
    }
    best.map(|(_, _, geometry)| geometry)
}

/// Why a [`Landscape`] or a [`Binning`] cannot be built.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum LandscapeError {
    /// A binning has no bins.
    NoBins,
    /// A binning's start or width is not finite, its width is not positive,
    /// or its last edge is not finite.
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
    /// A free energy needs a temperature: the distribution has no target and
    /// the options give none ([`LandscapeOptions::with_temperature`]).
    MissingTemperature,
    /// The options give another temperature than the distribution's target.
    ConflictingTemperature {
        /// The given temperature, in kelvin.
        given: f64,
        /// The target's temperature, in kelvin.
        target: f64,
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
            Self::ConflictingTemperature { given, target } => write!(
                f,
                "The temperature {given} K differs from the target temperature {target} K."
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
