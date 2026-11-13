//! Geometries processed as they are read: [`Processed`] and its
//! [`ProcessingSteps`].

use std::fmt;
use std::sync::{Arc, OnceLock};

use crate::alignment::{self, AlignmentOptions};
use crate::math::{PeriodicCell, Point3, Vector3};
use crate::properties::OwnerProperties;
use crate::structure::{
    subset_rows, DenseArrayError, Distribution, Geometry, GeometryError, GeometryView, Order,
    Positions, Structure, StructureError, StructureSetRead, StructureSetReadError,
};
use crate::topology::{AtomSelection, Topology, TopologySubsetError};

use super::periodic::{MoleculeImager, PeriodicError};

/// The steps a [`Processed`] reader runs on every geometry as it is read.
/// Each is optional; they run in this order:
///
/// 1. Make whole ([`Self::with_whole`]) or wrap ([`Self::with_wrap`]).
///    Making whole joins every molecule across the cell's boundaries, as
///    [`super::periodic::make_molecules_whole`] does. Wrapping makes every
///    molecule whole and images the molecules around the anchor atoms, as
///    [`super::periodic::image_molecules`] does: the molecules that hold an
///    anchor atom are imaged next to each other, the geometry is translated
///    so that the centroid of their atoms lies at the center of the cell
///    along every periodic axis, and every other molecule moves to its
///    periodic image nearest that centroid, so it lies within half a cell
///    of it. A geometry without a cell skips this step, where those two
///    functions fail with
///    [`PeriodicError::MissingCell`](super::periodic::PeriodicError::MissingCell),
///    so a trajectory whose frames have no box is read as stored;
///    [`Processed::skipped_periodic_steps`] tells such a geometry.
/// 2. Center ([`Self::with_center`]): the geometry translated so that the
///    centroid of the centered atoms lies at the center of its cell, half
///    the sum of all three cell vectors, periodic or not. Centering on the
///    anchor atoms after wrapping moves them to where wrapping put them, so
///    the wrapped molecules stay in the cell. A geometry without a cell is
///    translated so that the centroid lies where it lies in geometry 0 of the
///    reader below, as stored; in a reader where only some geometries have a
///    cell, the two kinds are centered on different targets.
/// 3. Fit ([`Self::with_fit`]): the geometry superposed onto a reference
///    through the fitted atoms, as
///    [`FrameSuperposer`](super::analysis::FrameSuperposer) does with
///    default options; positions move, and the cell, velocities, and forces
///    rotate.
/// 4. Smooth ([`Self::with_smoothing`]): positions averaged over a window of
///    geometries centered on the geometry, fewer at the ends, each processed by
///    steps 1–3. A geometry of the window with a cell has each of its
///    molecules moved to the periodic image nearest the same molecule of
///    the smoothed geometry first, as [`MoleculeImager::image_onto`] does, so
///    that a molecule crossing the cell between geometries is averaged whole.
///    Smoothing therefore expects whole molecules: make them whole or wrap
///    them. Everything but the positions is the geometry's own.
/// 5. Subset ([`Self::with_subset`]): the geometry projected onto the atoms of
///    a topology subset, which becomes the reader's topology.
///
/// Analyses read geometries without smoothing, which is for display: they read
/// a reader with [`Self::without_smoothing`].
///
/// ```
/// use kekule::topology::AtomSelection;
/// use kekule::traj::{FitReference, ProcessingSteps};
///
/// let topology = kekule::smiles::to_topology("CCO")?;
/// let all = AtomSelection::all(topology);
/// let steps = ProcessingSteps::default()
///     .with_whole(true)
///     .with_fit(all, FitReference::Geometry(0))
///     .with_smoothing(5)?;
/// assert_eq!(steps.smoothing(), Some(5));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Default)]
pub struct ProcessingSteps {
    whole: bool,
    wrap: Option<AtomSelection>,
    center: Option<AtomSelection>,
    fit: Option<(AtomSelection, FitReference)>,
    smoothing: Option<usize>,
    subset: Option<AtomSelection>,
}

/// What [`ProcessingSteps::with_fit`] superposes the geometries onto.
#[derive(Debug, Clone, PartialEq)]
pub enum FitReference {
    /// Geometry `n` of the same reader, after the steps before fitting.
    Geometry(usize),
    /// These positions, one per atom of the reader's topology.
    Positions(Positions),
}

/// One of the [`ProcessingSteps`], which [`ProcessingError`] names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ProcessingStep {
    /// Wrapping around anchor atoms.
    Wrap,
    /// Centering on atoms.
    Center,
    /// Fitting onto a reference.
    Fit,
    /// Projecting onto an atom subset.
    Subset,
}

impl fmt::Display for ProcessingStep {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Wrap => "wrapping",
            Self::Center => "centering",
            Self::Fit => "fitting",
            Self::Subset => "the subset",
        })
    }
}

impl ProcessingSteps {
    /// Makes every molecule whole across the cell's boundaries; wrapping
    /// does so too.
    #[must_use]
    pub fn with_whole(mut self, whole: bool) -> Self {
        self.whole = whole;
        self
    }

    /// Images every molecule whole around the molecules of `anchors`, whose
    /// combined center moves to the center of the cell.
    #[must_use]
    pub fn with_wrap(mut self, anchors: AtomSelection) -> Self {
        self.wrap = Some(anchors);
        self
    }

    /// Centers every geometry on the centroid of `atoms`.
    #[must_use]
    pub fn with_center(mut self, atoms: AtomSelection) -> Self {
        self.center = Some(atoms);
        self
    }

    /// Superposes every geometry onto `reference` through `atoms`.
    #[must_use]
    pub fn with_fit(mut self, atoms: AtomSelection, reference: FitReference) -> Self {
        self.fit = Some((atoms, reference));
        self
    }

    /// Averages positions over `window` geometries centered on each one.
    ///
    /// # Errors
    ///
    /// [`ProcessingError::InvalidWindow`] when `window` is not odd.
    pub fn with_smoothing(mut self, window: usize) -> Result<Self, ProcessingError> {
        if window.is_multiple_of(2) {
            return Err(ProcessingError::InvalidWindow { window });
        }
        self.smoothing = Some(window);
        Ok(self)
    }

    /// Projects every geometry onto the atoms of `atoms`, after the other
    /// steps.
    #[must_use]
    pub fn with_subset(mut self, atoms: AtomSelection) -> Self {
        self.subset = Some(atoms);
        self
    }

    /// The same steps without smoothing: what analyses read, through
    /// another [`Processed`] reader over the same reader below.
    ///
    /// ```
    /// use std::sync::Arc;
    ///
    /// use kekule::structure::{Positions, Geometry, StructureSetRead, StructureSet};
    /// use kekule::traj::{Processed, ProcessingSteps};
    ///
    /// let topology = kekule::smiles::to_topology("CCO")?;
    /// let mut set = StructureSet::new(topology.clone());
    /// for _ in 0..3 {
    ///     set.push(Geometry::new(Positions::zeros(topology.atom_count())), None)?;
    /// }
    /// let set: Arc<dyn StructureSetRead> = Arc::new(set);
    /// let steps = ProcessingSteps::default().with_whole(true).with_smoothing(3)?;
    /// let shown = Processed::new(Arc::clone(&set), steps.clone())?;
    /// // A background analysis owns its reader.
    /// let analysed: Arc<dyn StructureSetRead> =
    ///     Arc::new(Processed::new(set, steps.without_smoothing())?);
    /// assert_eq!(analysed.len(), shown.len());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[must_use]
    pub fn without_smoothing(mut self) -> Self {
        self.smoothing = None;
        self
    }

    /// Whether molecules are made whole without wrapping.
    pub const fn is_whole(&self) -> bool {
        self.whole
    }

    /// The anchor atoms of wrapping.
    pub const fn wrap(&self) -> Option<&AtomSelection> {
        self.wrap.as_ref()
    }

    /// The centered atoms.
    pub const fn center(&self) -> Option<&AtomSelection> {
        self.center.as_ref()
    }

    /// The fitted atoms and their reference.
    pub fn fit(&self) -> Option<(&AtomSelection, &FitReference)> {
        self.fit
            .as_ref()
            .map(|(atoms, reference)| (atoms, reference))
    }

    /// The smoothing window.
    pub const fn smoothing(&self) -> Option<usize> {
        self.smoothing
    }

    /// The atoms of the subset.
    pub const fn subset(&self) -> Option<&AtomSelection> {
        self.subset.as_ref()
    }

    /// Whether there is no step, so geometries are read as stored.
    pub const fn is_empty(&self) -> bool {
        !self.whole
            && self.wrap.is_none()
            && self.center.is_none()
            && self.fit.is_none()
            && self.smoothing.is_none()
            && self.subset.is_none()
    }
}

/// Failure to set up a [`Processed`] reader.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ProcessingError {
    /// A smoothing window is an odd number of geometries.
    InvalidWindow {
        /// The rejected window.
        window: usize,
    },
    /// The atoms of a step belong to another topology than the reader's.
    SelectionTopologyMismatch {
        /// The step.
        step: ProcessingStep,
    },
    /// The atoms of a step are none.
    EmptySelection {
        /// The step.
        step: ProcessingStep,
    },
    /// The reference geometry does not exist.
    ReferenceOutOfRange {
        /// The reference geometry index.
        index: usize,
        /// The number of geometries.
        len: usize,
    },
    /// The reference positions do not fit the reader's topology; see the
    /// source.
    Reference(Box<StructureError>),
    /// The subset cannot be made; see the source.
    Subset(TopologySubsetError),
    /// A step cannot process the first geometry or the reference; see the
    /// source.
    Geometry(StructureSetReadError),
}

impl fmt::Display for ProcessingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidWindow { window } => write!(
                formatter,
                "A smoothing window is an odd number of geometries, not {window}."
            ),
            Self::SelectionTopologyMismatch { step } => {
                write!(formatter, "The atoms of {step} belong to another topology.")
            }
            Self::EmptySelection { step } => {
                write!(formatter, "The atoms of {step} are none.")
            }
            Self::ReferenceOutOfRange { index, len } => write!(
                formatter,
                "The reference geometry {index} does not exist. There are {}.",
                crate::diagnostics::counted(*len, "geometry", "geometries")
            ),
            Self::Reference(_) => {
                formatter.write_str("The reference positions do not fit the topology.")
            }
            Self::Subset(error) => fmt::Display::fmt(error, formatter),
            Self::Geometry(error) => fmt::Display::fmt(error, formatter),
        }
    }
}

impl std::error::Error for ProcessingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Reference(error) => Some(error.as_ref()),
            Self::Subset(error) => std::error::Error::source(error),
            Self::Geometry(error) => std::error::Error::source(error),
            _ => None,
        }
    }
}

impl From<StructureSetReadError> for ProcessingError {
    fn from(error: StructureSetReadError) -> Self {
        Self::Geometry(error)
    }
}

/// The geometries of another reader with [`ProcessingSteps`] run on each as it
/// is read.
///
/// Reading geometry `i` reads geometry `i` of the reader below, and with
/// smoothing the geometries of its window, and runs the steps on them in the
/// order [`ProcessingSteps`] gives. A host caches geometries in front of this
/// reader when it reads them again. The order and the distribution are
/// those of the reader below.
///
/// Only geometry 0 is kept: [`Self::new`] processes it, [`StructureSetRead::peek`]
/// lends it, and [`StructureSetRead::read`] returns it, so a host's cache in front
/// of this reader takes it without processing it again. It stays as it was
/// processed then:
///
/// - Over a file later cut short, geometry 0 still reads while the geometries
///   after it fail.
/// - With smoothing over a reader whose length grows, geometry 0 keeps the
///   window the reader's length allowed then.
/// - It is returned shared, so a `Processed` reader on top of this one
///   copies it once more when its steps change it.
/// - Each reader keeps one whole geometry, with its velocities and forces:
///   about 72 MB at a million atoms, and as much again for a second reader
///   with [`ProcessingSteps::without_smoothing`].
///
/// The steps change the geometry they read in place. A geometry the reader
/// below owns alone, as a set in a file returns it, is not copied by
/// making whole, wrapping, centering, or fitting; one it shares, as a set
/// in memory, a host's cache, or geometry 0 of another `Processed` reader
/// returns it, is copied once, by the first of them that runs. Smoothing adds the geometries of its window into one sum
/// without copying them, and the subset copies the atoms it keeps.
///
/// Analyses read the geometries without smoothing, which is for display,
/// through a second reader with [`ProcessingSteps::without_smoothing`]. To
/// cache the geometries that smoothing averages, a host puts its cache in
/// front of a reader without smoothing and smooths through another
/// `Processed` reader, with smoothing alone, in front of the cache.
///
/// The reader below can be owned, borrowed (`&R`), or shared (`Arc<R>`),
/// such as an `Arc<dyn StructureSetRead>`. A host following a growing file takes
/// the reader below out with [`Self::into_parts`], refreshes it, and makes
/// the reader again with [`Self::new`].
///
/// ```
/// use kekule::math::Point3;
/// use kekule::structure::{Positions, Geometry, StructureSetRead, StructureSet};
/// use kekule::topology::AtomSelection;
/// use kekule::traj::{FitReference, Processed, ProcessingSteps};
/// use kekule::units::{Quantity, NANOMETER};
///
/// let topology = kekule::smiles::to_topology("CCO")?;
/// let mut set = StructureSet::new(topology.clone());
/// for x in [0.0, 1.0] {
///     let points = [Point3::new(x, 0.0, 0.0), Point3::new(x + 0.1, 0.0, 0.0), Point3::new(x, 0.1, 0.0)];
///     let positions = Positions::new(Quantity::new(points, NANOMETER))?;
///     set.push(Geometry::new(positions), None)?;
/// }
/// let all = AtomSelection::all(topology);
/// let steps = ProcessingSteps::default().with_fit(all, FitReference::Geometry(0));
/// let fitted = Processed::new(&set, steps)?;
/// let geometry = fitted.read(1)?;
/// let x = geometry.positions().values().value()[0].x;
/// assert!(x.abs() < 1.0e-12);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct Processed<R> {
    inner: R,
    steps: ProcessingSteps,
    /// The topology of the processed geometries: the subset's, or the inner
    /// reader's.
    topology: Arc<Topology>,
    /// Plans molecules for making whole, wrapping, and smoothing.
    imager: Option<MoleculeImager>,
    /// The molecules wrapping images around.
    anchors: Option<Vec<bool>>,
    /// Where centering puts the centered atoms in geometries without a cell:
    /// their centroid in geometry 0 of the reader below.
    rest: OnceLock<Vector3>,
    /// The structure fitting superposes onto.
    reference: Option<Structure>,
    /// The dense atom and bond rows the subset keeps.
    subset: Option<(Vec<usize>, Vec<usize>)>,
    /// Processed geometry 0, when there is one.
    first: Option<Geometry>,
}

impl<R: StructureSetRead> Processed<R> {
    /// The geometries of `inner` processed by `steps`. The steps are checked
    /// against the reader's topology, the reference made, and geometry 0
    /// processed, so steps that cannot process the geometries fail here.
    /// Geometry 0 is kept and lent.
    ///
    /// # Errors
    ///
    /// [`ProcessingError`] for atoms of another topology or no atoms, a
    /// reference geometry past the last or reference positions of another
    /// atom count, a subset that cannot be made, and a failure to read or
    /// process the reference or geometry 0.
    pub fn new(inner: R, steps: ProcessingSteps) -> Result<Self, ProcessingError> {
        let source = Arc::clone(inner.topology());
        for (step, atoms) in [
            (ProcessingStep::Wrap, steps.wrap.as_ref()),
            (ProcessingStep::Center, steps.center.as_ref()),
            (
                ProcessingStep::Fit,
                steps.fit.as_ref().map(|(atoms, _)| atoms),
            ),
            (ProcessingStep::Subset, steps.subset.as_ref()),
        ] {
            if let Some(atoms) = atoms {
                check_atoms(&source, atoms, step)?;
            }
        }
        let imager = (steps.whole || steps.wrap.is_some() || steps.smoothing.is_some())
            .then(|| MoleculeImager::new(Arc::clone(&source)));
        let anchors = match (&imager, &steps.wrap) {
            (Some(imager), Some(anchors)) => Some(
                imager
                    .anchor_groups(anchors)
                    .expect("the anchors were checked"),
            ),
            _ => None,
        };
        let (topology, subset) = match &steps.subset {
            Some(atoms) => {
                let subset = Topology::subset(Arc::clone(&source), atoms)
                    .map_err(ProcessingError::Subset)?;
                let rows = subset_rows(&subset);
                (subset.shared_topology(), Some(rows))
            }
            None => (Arc::clone(&source), None),
        };
        let mut processed = Self {
            inner,
            steps,
            topology,
            imager,
            anchors,
            rest: OnceLock::new(),
            reference: None,
            subset,
            first: None,
        };

        if let Some((_, reference)) = &processed.steps.fit {
            processed.reference = Some(processed.reference_structure(reference)?);
        }
        if !processed.inner.is_empty() {
            processed.first = Some(processed.processed(0)?);
        }
        Ok(processed)
    }

    /// The reader below.
    pub const fn inner(&self) -> &R {
        &self.inner
    }

    /// The steps.
    pub const fn steps(&self) -> &ProcessingSteps {
        &self.steps
    }

    /// The reader below and the steps.
    pub fn into_parts(self) -> (R, ProcessingSteps) {
        (self.inner, self.steps)
    }

    /// Whether this reader passed `geometry`, one it read, through making
    /// whole and wrapping unchanged because the geometry has no cell: the
    /// steps make whole or wrap, and the geometry has no cell. Every step
    /// keeps a geometry's cell or its lack of one, so the geometry read tells.
    /// The other steps ran on it as on any geometry.
    ///
    /// ```
    /// use kekule::structure::{Positions, Geometry, StructureSetRead, StructureSet};
    /// use kekule::traj::{Processed, ProcessingSteps};
    ///
    /// let topology = kekule::smiles::to_topology("O")?;
    /// let mut set = StructureSet::new(topology.clone());
    /// set.push(Geometry::new(Positions::zeros(topology.atom_count())), None)?;
    /// let whole = Processed::new(&set, ProcessingSteps::default().with_whole(true))?;
    /// assert!(whole.skipped_periodic_steps(&whole.read(0)?));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn skipped_periodic_steps(&self, geometry: &Geometry) -> bool {
        (self.steps.whole || self.steps.wrap.is_some()) && geometry.cell().is_none()
    }

    /// The structure fitting superposes onto.
    fn reference_structure(&self, reference: &FitReference) -> Result<Structure, ProcessingError> {
        match reference {
            FitReference::Geometry(index) => {
                let len = self.inner.len();
                if *index >= len {
                    return Err(ProcessingError::ReferenceOutOfRange { index: *index, len });
                }
                let mut geometry = self.read_inner(*index)?;
                self.prepare(*index, &mut geometry)?;
                Ok(self
                    .view(*index, &geometry)
                    .as_structure_view()
                    .to_structure())
            }
            FitReference::Positions(positions) => {
                Structure::new(Arc::clone(self.inner.topology()), positions.clone())
                    .map_err(|error| ProcessingError::Reference(Box::new(error)))
            }
        }
    }

    /// Geometry `index` of the reader below, checked to be bound to its
    /// topology, as every step expects.
    fn read_inner(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        let geometry = self.inner.read(index)?;
        GeometryView::new(&self.inner, index, &geometry)?;
        Ok(geometry)
    }

    /// A view of `geometry`, geometry `index` of the reader below.
    fn view<'a>(
        &'a self,
        index: usize,
        geometry: &'a Geometry,
    ) -> GeometryView<'a, dyn StructureSetRead> {
        GeometryView::of_reader(&self.inner, index, geometry)
    }

    /// Wraps or makes whole, then centers, geometry `index` of the reader
    /// below, in place.
    fn prepare(&self, index: usize, geometry: &mut Geometry) -> Result<(), StructureSetReadError> {
        if let (Some(imager), Some(cell)) = (&self.imager, geometry.cell().copied()) {
            if self.steps.whole || self.anchors.is_some() {
                imager
                    .whole_in_place(
                        index,
                        cell,
                        geometry.positions_in_place(),
                        self.anchors.as_deref(),
                    )
                    .map_err(periodic(index))?;
            }
        }
        if let Some(center) = &self.steps.center {
            let target = match geometry.cell() {
                Some(cell) => cell_center(cell),
                None => self.rest(center)?,
            };
            let shift = target - centroid(geometry.positions().values().value(), center);
            geometry
                .translate_in_place(shift)
                .map_err(geometry_error(index))?;
        }
        Ok(())
    }

    /// The centroid of the centered atoms in geometry 0 of the reader below.
    fn rest(&self, center: &AtomSelection) -> Result<Vector3, StructureSetReadError> {
        if let Some(rest) = self.rest.get() {
            return Ok(*rest);
        }
        let first = self.read_inner(0)?;
        let rest = centroid(first.positions().values().value(), center);
        Ok(*self.rest.get_or_init(|| rest))
    }

    /// Geometry `index` of the reader below after steps 1–3, before the
    /// subset.
    fn stepped(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        let mut geometry = self.read_inner(index)?;
        self.prepare(index, &mut geometry)?;
        if let (Some((atoms, _)), Some(reference)) = (&self.steps.fit, &self.reference) {
            let fitted = alignment::fit(
                self.view(index, &geometry).as_structure_view(),
                reference.as_structure_view(),
                atoms,
                AlignmentOptions::default(),
            )
            .map_err(|source| StructureSetReadError::Alignment {
                geometry: index,
                source: Box::new(source),
            })?;
            geometry
                .transform_in_place(&fitted.transform())
                .map_err(geometry_error(index))?;
        }
        Ok(geometry)
    }

    /// Geometry `index` after steps 1–3, averaged over the window of `half`
    /// geometries on each side.
    fn smoothed(&self, index: usize, half: usize) -> Result<Geometry, StructureSetReadError> {
        let mut center = self.stepped(index)?;
        let last = self.inner.len().saturating_sub(1);
        let window = index.saturating_sub(half)..=(index + half).min(last);
        if window.start() == window.end() {
            return Ok(center);
        }
        let mut sum: Vec<Vector3> = center
            .positions()
            .values()
            .value()
            .iter()
            .map(|point| *point - Point3::origin())
            .collect();
        let mut count = 1.0;
        let imager = self.imager.as_ref().expect("smoothing plans molecules");
        for other in window {
            if other == index {
                continue;
            }
            let geometry = self.stepped(other)?;
            let points = geometry.positions().values().into_value();
            match geometry.cell() {
                Some(cell) => imager
                    .add_imaged_onto(
                        other,
                        *cell,
                        points,
                        center.positions().values().value(),
                        &mut sum,
                    )
                    .map_err(periodic(other))?,
                None => add(&mut sum, points),
            }
            count += 1.0;
        }
        for (atom, (point, total)) in center.positions_in_place().iter_mut().zip(&sum).enumerate() {
            *point = Point3::origin() + *total / count;
            if !point.is_finite() {
                return Err(geometry_error(index)(GeometryError::Array(
                    DenseArrayError::NonFinite { index: atom },
                )));
            }
        }
        Ok(center)
    }

    /// Geometry `index` read from the reader below and processed.
    fn processed(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        let geometry = match self.steps.smoothing {
            Some(window) => self.smoothed(index, window / 2)?,
            None => self.stepped(index)?,
        };
        self.finish(index, geometry)
    }

    /// Projects geometry `index` onto the subset.
    fn finish(&self, index: usize, geometry: Geometry) -> Result<Geometry, StructureSetReadError> {
        match &self.subset {
            Some((atoms, bonds)) => {
                geometry
                    .projected(atoms, bonds)
                    .map_err(|source| StructureSetReadError::Geometry {
                        geometry: index,
                        source: Box::new(source),
                    })
            }
            None => Ok(geometry),
        }
    }
}

impl<R: StructureSetRead> StructureSetRead for Processed<R> {
    fn topology(&self) -> &Arc<Topology> {
        &self.topology
    }

    fn len(&self) -> usize {
        self.inner.len()
    }

    fn order(&self) -> Option<&Order> {
        self.inner.order()
    }

    fn distribution(&self) -> Option<&Distribution> {
        self.inner.distribution()
    }

    fn properties(&self) -> Option<&OwnerProperties> {
        self.inner.properties()
    }

    /// Geometry 0, which the reader keeps; `None` for every other geometry.
    fn peek(&self, index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        match (index, &self.first) {
            (0, Some(first)) => Some(GeometryView::of_reader(self, 0, first)),
            _ => None,
        }
    }

    /// Geometry `index`, processed. Geometry 0 is a clone of the kept geometry,
    /// which shares its storage.
    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        match (index, &self.first) {
            (0, Some(first)) => Ok(first.clone()),
            _ => self.processed(index),
        }
    }
}

/// Reports a periodic step's failure on geometry `geometry`.
fn periodic(geometry: usize) -> impl FnOnce(PeriodicError) -> StructureSetReadError {
    move |source| StructureSetReadError::Periodic {
        geometry,
        source: Box::new(source),
    }
}

/// Reports a step's failure to store its result in geometry `geometry`.
fn geometry_error(geometry: usize) -> impl FnOnce(GeometryError) -> StructureSetReadError {
    move |source| StructureSetReadError::Geometry {
        geometry,
        source: Box::new(source),
    }
}

/// Checks that the atoms of `step` are some atoms of `topology`.
fn check_atoms(
    topology: &Topology,
    atoms: &AtomSelection,
    step: ProcessingStep,
) -> Result<(), ProcessingError> {
    if !topology.shares_layout(atoms.topology()) {
        return Err(ProcessingError::SelectionTopologyMismatch { step });
    }
    if atoms.indices().is_empty() {
        return Err(ProcessingError::EmptySelection { step });
    }
    Ok(())
}

/// The centroid of the `atoms` of `points`, from the origin.
fn centroid(points: &[Point3], atoms: &AtomSelection) -> Vector3 {
    let indices = atoms.indices();
    let mut sum = Vector3::zero();
    for index in indices {
        sum += points[index.index()] - Point3::origin();
    }
    sum / indices.len() as f64
}

/// The center of a cell, from the origin.
fn cell_center(cell: &PeriodicCell) -> Vector3 {
    let [a, b, c] = cell.vectors().into_value();
    (a + b + c) * 0.5
}

/// Adds `points` to `sum`, atom by atom.
fn add(sum: &mut [Vector3], points: &[Point3]) {
    for (total, point) in sum.iter_mut().zip(points) {
        *total += *point - Point3::origin();
    }
}
