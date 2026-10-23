//! Interprets XYZ frames as structures whose connectivity is inferred from
//! interatomic distances.
use std::collections::{BTreeSet, HashMap};

use crate::diagnostics::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, LimitExceeded, Resource,
    SourceLocation, UnsupportedContent,
};
use crate::geometry::Point3;
use crate::molecule::{Atom, AtomId, BondOrder, Element, ImplicitHydrogens, MoleculeEditor};
use crate::structure::{Ensemble, EnsembleSample, EnsembleState, Model, ModelBuilder, Positions};
use crate::topology::{AtomSiteMetadata, InstanceAtomId, Topology, TopologyAtomIndex};
use crate::units::{Quantity, Unit, ANGSTROM};

use super::structure_connectivity::{
    covalent_bonding_radius, COVALENT_BOND_TOLERANCE_ANGSTROM, FALLBACK_COVALENT_RADIUS_ANGSTROM,
    MIN_COVALENT_BOND_DISTANCE_ANGSTROM,
};
use super::structure_interpret::is_metal;
use super::{XyzAtom, XyzDocument, XyzFrame};

/// Which frame [`XyzDocument::interpret`] interprets.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum XyzFrameSelection {
    /// The only frame; a document with several fails with
    /// [`DiagnosticKind::Inconsistent`].
    #[default]
    RequireSingle,
    /// The first frame in source order.
    First,
    /// The frame at this zero-based index.
    Index(usize),
}

/// Options shared by single-frame and ensemble interpretation.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Common {
    length_unit: Unit,
    extra_columns: UnsupportedContent,
    max_atom_contacts: usize,
}

impl Default for Common {
    fn default() -> Self {
        Self {
            length_unit: ANGSTROM,
            extra_columns: UnsupportedContent::Reject,
            max_atom_contacts: DEFAULT_XYZ_MAX_ATOM_CONTACTS,
        }
    }
}

impl Common {
    fn with_length_unit(mut self, length_unit: Unit) -> Result<Self, FormatError> {
        if length_unit.dimension() != ANGSTROM.dimension() {
            return Err(option_error(format!(
                "`{length_unit}` is not a length unit"
            )));
        }
        self.length_unit = length_unit;
        Ok(self)
    }

    fn with_max_atom_contacts(mut self, max_atom_contacts: usize) -> Result<Self, FormatError> {
        if max_atom_contacts == 0 {
            return Err(option_error("the contact limit must be at least one"));
        }
        self.max_atom_contacts = max_atom_contacts;
        Ok(self)
    }
}

fn option_error(message: impl Into<String>) -> FormatError {
    FormatError::interpret(Format::Xyz, DiagnosticKind::InvalidValue, message)
}

/// The largest coordinate magnitude, in ångström, that distance inference
/// indexes; larger ones fail rather than collapse into one grid cell.
const MAX_COORDINATE_ANGSTROM: f64 = 1.0e15;

/// The [`Default`] bound of [`XyzInterpretOptions::max_atom_contacts`].
///
/// No atom of a molecule or crystal has this many atoms within covalent
/// bonding distance, so reaching it means overlapping or corrupt coordinates.
pub const DEFAULT_XYZ_MAX_ATOM_CONTACTS: usize = 32;

macro_rules! common_options {
    () => {
        /// The length unit of the file's coordinates; ångström by default.
        pub const fn length_unit(&self) -> Unit {
            self.common.length_unit
        }

        /// Sets [`Self::length_unit`].
        ///
        /// # Errors
        ///
        /// A [`FormatError`] of kind [`DiagnosticKind::InvalidValue`] when
        /// `length_unit` is not a unit of length.
        pub fn with_length_unit(mut self, length_unit: Unit) -> Result<Self, FormatError> {
            self.common = self.common.with_length_unit(length_unit)?;
            Ok(self)
        }

        /// What happens to columns after the coordinates, such as extended
        /// XYZ per-atom properties, which Kekule does not interpret.
        /// [`UnsupportedContent::Reject`] by default;
        /// [`UnsupportedContent::Report`] omits them with one
        /// [`XyzInterpretIssue::ExtraColumnsOmitted`] per atom line.
        pub const fn extra_columns(&self) -> UnsupportedContent {
            self.common.extra_columns
        }

        /// Sets [`Self::extra_columns`].
        #[must_use]
        pub const fn with_extra_columns(mut self, extra_columns: UnsupportedContent) -> Self {
            self.common.extra_columns = extra_columns;
            self
        }

        /// The most atoms one atom may have within covalent bonding distance,
        /// bonded or too close to bond; [`DEFAULT_XYZ_MAX_ATOM_CONTACTS`] by
        /// default. It bounds the work of distance inference on overlapping
        /// coordinates.
        pub const fn max_atom_contacts(&self) -> usize {
            self.common.max_atom_contacts
        }

        /// Sets [`Self::max_atom_contacts`].
        ///
        /// # Errors
        ///
        /// A [`FormatError`] of kind [`DiagnosticKind::InvalidValue`] for zero.
        pub fn with_max_atom_contacts(
            mut self,
            max_atom_contacts: usize,
        ) -> Result<Self, FormatError> {
            self.common = self.common.with_max_atom_contacts(max_atom_contacts)?;
            Ok(self)
        }
    };
}

/// Options of [`XyzDocument::interpret`].
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct XyzInterpretOptions {
    frame: XyzFrameSelection,
    common: Common,
}

impl XyzInterpretOptions {
    /// Which frame to interpret; [`XyzFrameSelection::RequireSingle`] by
    /// default.
    pub const fn frame(&self) -> XyzFrameSelection {
        self.frame
    }

    /// Sets [`Self::frame`].
    #[must_use]
    pub const fn with_frame(mut self, frame: XyzFrameSelection) -> Self {
        self.frame = frame;
        self
    }

    common_options!();
}

/// Options of [`XyzDocument::interpret_ensemble`].
#[derive(Debug, Clone, Default, PartialEq)]
pub struct XyzEnsembleInterpretOptions {
    frames: Option<Vec<usize>>,
    common: Common,
}

impl XyzEnsembleInterpretOptions {
    /// Zero-based indices of the frames to interpret, in member order;
    /// `None`, the default, selects every frame in source order.
    pub fn frames(&self) -> Option<&[usize]> {
        self.frames.as_deref()
    }

    /// Sets [`Self::frames`].
    ///
    /// # Errors
    ///
    /// A [`FormatError`] with
    /// [`FormatErrorDetail::EmptyModelSelection`] for `Some` of an empty list.
    pub fn with_frames(mut self, frames: Option<Vec<usize>>) -> Result<Self, FormatError> {
        if frames.as_ref().is_some_and(Vec::is_empty) {
            return Err(interpret_error(
                DiagnosticKind::Inconsistent,
                "explicit frame selection is empty",
            )
            .with_detail(FormatErrorDetail::EmptyModelSelection));
        }
        self.frames = frames;
        Ok(self)
    }

    common_options!();
}

/// One bond inferred from an interatomic distance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XyzInferredBond {
    atoms: [TopologyAtomIndex; 2],
    order: BondOrder,
    distance: f64,
    cutoff: f64,
}

impl XyzInferredBond {
    /// The bonded atoms by dense index, which is their atom-line index in the
    /// frame, lower index first.
    pub const fn atoms(&self) -> [TopologyAtomIndex; 2] {
        self.atoms
    }

    /// [`BondOrder::Dative`] between a metal and a non-metal, otherwise
    /// [`BondOrder::Single`].
    pub const fn order(&self) -> BondOrder {
        self.order
    }

    /// The interatomic distance.
    pub fn distance(&self) -> Quantity<f64> {
        Quantity::new(self.distance, ANGSTROM)
    }

    /// The bonding cutoff for the pair: the sum of both covalent radii and
    /// the tolerance.
    pub fn cutoff(&self) -> Quantity<f64> {
        Quantity::new(self.cutoff, ANGSTROM)
    }
}

/// A non-fatal finding of XYZ interpretation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum XyzInterpretIssue {
    /// The element has no tabulated covalent radius, so distance inference
    /// used [`XyzInterpretationReport::fallback_covalent_radius`] for this
    /// atom.
    FallbackCovalentRadius {
        /// The atom.
        atom: TopologyAtomIndex,
        /// Its element.
        element: Element,
    },
    /// Two atoms are closer than
    /// [`XyzInterpretationReport::minimum_bond_distance`]; they are left
    /// unbonded.
    CloseContact {
        /// The atoms, lower index first.
        atoms: [TopologyAtomIndex; 2],
        /// Their distance.
        distance: Quantity<f64>,
    },
    /// The atom line's columns after the coordinates were omitted, as
    /// [`XyzInterpretOptions::extra_columns`] allows.
    ExtraColumnsOmitted {
        /// The atom.
        atom: TopologyAtomIndex,
        /// The atom line.
        location: SourceLocation,
    },
}

/// What interpretation inferred from one frame and why.
///
/// Connectivity comes from distances alone: two atoms are bonded when they
/// lie within the sum of their covalent radii (Cordero et al. 2008, as
/// [`Element::covalent_radius`] gives them) plus
/// [`Self::bond_tolerance`], and farther apart than
/// [`Self::minimum_bond_distance`]. This is the rule macromolecular
/// interpretation applies to residues without a template.
#[derive(Debug, Clone, PartialEq)]
pub struct XyzInterpretationReport {
    frame: usize,
    comment: String,
    bonds: Vec<XyzInferredBond>,
    issues: Vec<XyzInterpretIssue>,
}

impl XyzInterpretationReport {
    /// The zero-based index of the interpreted frame.
    pub const fn frame(&self) -> usize {
        self.frame
    }

    /// The frame's comment line verbatim. Kekule does not interpret it,
    /// including extended XYZ keys such as `Lattice`.
    pub fn comment(&self) -> &str {
        &self.comment
    }

    /// Every bond of the structure, in dense atom order of its endpoints;
    /// all of them are inferred.
    pub fn inferred_bonds(&self) -> &[XyzInferredBond] {
        &self.bonds
    }

    /// Fallback radii, close contacts, and omitted columns, in that order.
    pub fn issues(&self) -> &[XyzInterpretIssue] {
        &self.issues
    }

    /// The distance added to the sum of covalent radii: 0.45 Å.
    pub fn bond_tolerance(&self) -> Quantity<f64> {
        Quantity::new(COVALENT_BOND_TOLERANCE_ANGSTROM, ANGSTROM)
    }

    /// The distance at or below which atoms are too close to bond: 0.4 Å.
    pub fn minimum_bond_distance(&self) -> Quantity<f64> {
        Quantity::new(MIN_COVALENT_BOND_DISTANCE_ANGSTROM, ANGSTROM)
    }

    /// The covalent radius assumed for an element without a tabulated one:
    /// 0.77 Å.
    pub fn fallback_covalent_radius(&self) -> Quantity<f64> {
        Quantity::new(FALLBACK_COVALENT_RADIUS_ANGSTROM, ANGSTROM)
    }
}

/// One frame interpreted as a model, with its report.
#[derive(Debug, Clone)]
pub struct XyzInterpretation {
    model: Model,
    report: XyzInterpretationReport,
}

impl XyzInterpretation {
    /// The interpreted model.
    pub fn model(&self) -> &Model {
        &self.model
    }

    /// The model's topology.
    pub fn topology(&self) -> &Topology {
        self.model.topology()
    }

    /// What interpretation inferred and why.
    pub fn report(&self) -> &XyzInterpretationReport {
        &self.report
    }

    /// The interpreted model, discarding the report.
    pub fn into_model(self) -> Model {
        self.model
    }

    /// The model and its report.
    pub fn into_parts(self) -> (Model, XyzInterpretationReport) {
        (self.model, self.report)
    }
}

/// Several frames interpreted as one equally weighted ensemble, with one
/// report per member.
#[derive(Debug, Clone)]
pub struct XyzEnsembleInterpretation {
    ensemble: Ensemble,
    reports: Vec<XyzInterpretationReport>,
}

impl XyzEnsembleInterpretation {
    /// The interpreted ensemble, one member per selected frame.
    pub fn ensemble(&self) -> &Ensemble {
        &self.ensemble
    }

    /// One report per member, in member order.
    pub fn reports(&self) -> &[XyzInterpretationReport] {
        &self.reports
    }

    /// The topology all members share.
    pub fn topology(&self) -> &Topology {
        self.ensemble.topology()
    }

    /// The interpreted ensemble, discarding the reports.
    pub fn into_ensemble(self) -> Ensemble {
        self.ensemble
    }

    /// The ensemble and its per-member reports.
    pub fn into_parts(self) -> (Ensemble, Vec<XyzInterpretationReport>) {
        (self.ensemble, self.reports)
    }
}

fn interpret_error(kind: DiagnosticKind, message: impl Into<String>) -> FormatError {
    FormatError::interpret(Format::Xyz, kind, message)
}

/// Names the frame of a failure that carries no more specific detail, such
/// as an exceeded limit; its location names the line either way.
fn in_frame(error: FormatError, frame: usize) -> FormatError {
    if error.detail().is_some() {
        return error;
    }
    error.with_detail(FormatErrorDetail::Model {
        model_id: frame.to_string(),
    })
}

fn no_frames() -> FormatError {
    interpret_error(DiagnosticKind::Inconsistent, "input contains no frames")
        .with_detail(FormatErrorDetail::NoCoordinateModels)
}

pub(crate) fn interpret_xyz(
    document: &XyzDocument,
    options: &XyzInterpretOptions,
) -> Result<XyzInterpretation, FormatError> {
    let frames = document.frames();
    let index = match options.frame {
        _ if frames.is_empty() => return Err(no_frames()),
        XyzFrameSelection::RequireSingle if frames.len() > 1 => {
            return Err(interpret_error(
                DiagnosticKind::Inconsistent,
                format!(
                    "the input has {} frames; select one or interpret an ensemble",
                    frames.len()
                ),
            ))
        }
        XyzFrameSelection::RequireSingle | XyzFrameSelection::First => 0,
        XyzFrameSelection::Index(index) => {
            if index >= frames.len() {
                return Err(unknown_frame(index));
            }
            index
        }
    };
    let frame = Frame::analyse(&frames[index], index, &options.common)
        .map_err(|error| in_frame(error, index))?;
    let model = frame
        .build_model(options.common.length_unit)
        .map_err(|error| in_frame(error, index))?;
    Ok(XyzInterpretation {
        model,
        report: frame.report,
    })
}

fn unknown_frame(index: usize) -> FormatError {
    interpret_error(
        DiagnosticKind::Inconsistent,
        format!("frame {index} is not present"),
    )
    .with_detail(FormatErrorDetail::UnknownRequestedModel {
        model_id: index.to_string(),
    })
}

pub(crate) fn interpret_xyz_ensemble(
    document: &XyzDocument,
    options: &XyzEnsembleInterpretOptions,
) -> Result<XyzEnsembleInterpretation, FormatError> {
    let frames = document.frames();
    if frames.is_empty() {
        return Err(no_frames());
    }
    let selected = options
        .frames
        .clone()
        .unwrap_or_else(|| (0..frames.len()).collect());
    let mut seen = BTreeSet::new();
    for &index in &selected {
        if !seen.insert(index) {
            return Err(interpret_error(
                DiagnosticKind::Inconsistent,
                format!("frame {index} was requested more than once"),
            )
            .with_detail(FormatErrorDetail::DuplicateRequestedModel {
                model_id: index.to_string(),
            }));
        }
        if index >= frames.len() {
            return Err(unknown_frame(index));
        }
    }
    let common = &options.common;
    let first_index = selected[0];
    let first = Frame::analyse(&frames[first_index], first_index, common)
        .map_err(|error| in_frame(error, first_index))?;
    let model = first
        .build_model(common.length_unit)
        .map_err(|error| in_frame(error, first_index))?;
    let (first_model, first_report) = (model, first.report);
    let mut ensemble = Ensemble::new(first_model.shared_topology());
    let shared_bonds = bond_identity(&first_report.bonds);
    let mut reports = Vec::with_capacity(selected.len());
    push_sample(&mut ensemble, first.points, common.length_unit, first_index)?;
    reports.push(first_report);
    for &index in &selected[1..] {
        let frame = &frames[index];
        let mut same_atoms = frame.atoms().len() == first.elements.len();
        for (atom, &element) in frame.atoms().iter().zip(&first.elements) {
            let resolved = resolve_element(atom).map_err(|error| in_frame(error, index))?;
            same_atoms &= resolved == element;
        }
        if !same_atoms {
            return Err(interpret_error(
                DiagnosticKind::Inconsistent,
                format!("frame {index} has different atoms than frame {first_index}"),
            )
            .at(frame.location())
            .with_detail(FormatErrorDetail::InconsistentAtomSet {
                model_id: index.to_string(),
            }));
        }
        let analysed =
            Frame::analyse(frame, index, common).map_err(|error| in_frame(error, index))?;
        if bond_identity(&analysed.report.bonds) != shared_bonds {
            return Err(interpret_error(
                DiagnosticKind::Inconsistent,
                format!(
                    "frame {index} infers different bonds than frame {first_index} from its distances"
                ),
            )
            .at(frame.location())
            .with_detail(FormatErrorDetail::InconsistentTopology {
                model_id: index.to_string(),
            }));
        }
        push_sample(&mut ensemble, analysed.points, common.length_unit, index)?;
        reports.push(analysed.report);
    }
    drop(first_model);
    Ok(XyzEnsembleInterpretation { ensemble, reports })
}

fn bond_identity(bonds: &[XyzInferredBond]) -> Vec<([TopologyAtomIndex; 2], BondOrder)> {
    bonds.iter().map(|bond| (bond.atoms, bond.order)).collect()
}

/// Appends one frame as an equally weighted sample of the ensemble's one
/// state.
fn push_sample(
    ensemble: &mut Ensemble,
    points: Vec<Point3>,
    unit: Unit,
    frame: usize,
) -> Result<(), FormatError> {
    let chemistry = |error: &dyn std::fmt::Display| {
        in_frame(
            interpret_error(
                DiagnosticKind::Chemistry,
                format!("cannot add the frame to the ensemble: {error}"),
            ),
            frame,
        )
    };
    let positions =
        Positions::from_vec(Quantity::new(points, unit)).map_err(|error| chemistry(&error))?;
    let sample = EnsembleSample::new(positions, 1.0).map_err(|error| chemistry(&error))?;
    if ensemble.is_empty() {
        ensemble
            .push_state(sample, EnsembleState::new())
            .map(|_| ())
            .map_err(|error| chemistry(&error))
    } else {
        ensemble
            .push_sample(sample, 0)
            .map_err(|error| chemistry(&error))
    }
}

/// The element an atom line names: a symbol in any letter case, or an
/// atomic number.
fn resolve_element(atom: &XyzAtom) -> Result<Element, FormatError> {
    let token = atom.element();
    let element = if token.bytes().all(|byte| byte.is_ascii_digit()) {
        token
            .parse::<u8>()
            .ok()
            .and_then(Element::from_atomic_number)
    } else {
        let mut characters = token.chars();
        characters.next().and_then(|first| {
            let symbol = first
                .to_uppercase()
                .chain(characters.flat_map(char::to_lowercase))
                .collect::<String>();
            Element::from_symbol(&symbol)
        })
    };
    element.ok_or_else(|| {
        interpret_error(
            DiagnosticKind::InvalidValue,
            format!("`{token}` is not an element symbol or atomic number"),
        )
        .at(atom.location())
    })
}

/// One frame's elements, coordinates, and inferred connectivity.
struct Frame {
    elements: Vec<Element>,
    points: Vec<Point3>,
    report: XyzInterpretationReport,
}

impl Frame {
    fn analyse(frame: &XyzFrame, index: usize, common: &Common) -> Result<Self, FormatError> {
        let scale = Quantity::new(1.0, common.length_unit)
            .value_in(ANGSTROM)
            .expect("validated length unit converts to angstrom");
        let mut elements = Vec::with_capacity(frame.atoms().len());
        let mut points = Vec::with_capacity(frame.atoms().len());
        let mut fallback = Vec::new();
        let mut omitted = Vec::new();
        for (row, atom) in frame.atoms().iter().enumerate() {
            let element = resolve_element(atom)?;
            let atom_index = dense_index(row);
            if covalent_bonding_radius(element).1 {
                fallback.push(XyzInterpretIssue::FallbackCovalentRadius {
                    atom: atom_index,
                    element,
                });
            }
            if atom.extra().is_some() {
                match common.extra_columns {
                    UnsupportedContent::Report => {
                        omitted.push(XyzInterpretIssue::ExtraColumnsOmitted {
                            atom: atom_index,
                            location: atom.location(),
                        });
                    }
                    _ => {
                        return Err(interpret_error(
                            DiagnosticKind::Unsupported,
                            "the atom line has columns after its coordinates, which Kekule does not interpret",
                        )
                        .at(atom.location()))
                    }
                }
            }
            let position = atom.position();
            if [position.x, position.y, position.z]
                .into_iter()
                .any(|value| (value * scale).abs() > MAX_COORDINATE_ANGSTROM)
            {
                return Err(interpret_error(
                    DiagnosticKind::InvalidValue,
                    "a coordinate exceeds 1e15 Å, beyond what distance inference indexes",
                )
                .at(atom.location()));
            }
            elements.push(element);
            points.push(position);
        }
        let angstrom = points
            .iter()
            .map(|point| Point3::new(point.x * scale, point.y * scale, point.z * scale))
            .collect::<Vec<_>>();
        let (bonds, close) = infer_contacts(&elements, &angstrom, common.max_atom_contacts)
            .map_err(|limit| {
                interpret_error(
                    DiagnosticKind::Limit,
                    format!(
                        "atom {} has more atoms within bonding distance than allowed",
                        limit.1
                    ),
                )
                .at(frame.atoms()[limit.1].location())
                .exceeding(limit.0)
            })?;
        let mut issues = fallback;
        issues.extend(close);
        issues.extend(omitted);
        Ok(Self {
            elements,
            points,
            report: XyzInterpretationReport {
                frame: index,
                comment: frame.comment().to_owned(),
                bonds,
                issues,
            },
        })
    }

    /// Publishes one molecule per connected component, a synthetic hierarchy
    /// of one `UNL` residue per molecule in chain `A`, and atom-line order as
    /// dense atom order.
    fn build_model(&self, unit: Unit) -> Result<Model, FormatError> {
        let chemistry = |error: &dyn std::fmt::Display| {
            interpret_error(
                DiagnosticKind::Chemistry,
                format!("cannot publish the structure: {error}"),
            )
        };
        let components = components(self.elements.len(), &self.report.bonds);
        let mut builder = ModelBuilder::new();
        let chain = builder
            .hierarchy_mut()
            .add_chain("A", None)
            .map_err(|error| chemistry(&error))?;
        let mut order = vec![None; self.elements.len()];
        let mut bonds_by_component = vec![Vec::new(); components.len()];
        let mut component_of = vec![0; self.elements.len()];
        for (component, rows) in components.iter().enumerate() {
            for &row in rows {
                component_of[row] = component;
            }
        }
        for bond in &self.report.bonds {
            let row = bond.atoms[0].index();
            bonds_by_component[component_of[row]].push(bond);
        }
        for (component, rows) in components.iter().enumerate() {
            let mut editor = MoleculeEditor::new();
            let mut local = HashMap::with_capacity(rows.len());
            for &row in rows {
                // XYZ lists every atom, so no hydrogen is implicit.
                let atom =
                    Atom::new(self.elements[row]).with_hydrogens(ImplicitHydrogens::Fixed(0));
                let id = editor.add_atom(atom).map_err(|error| chemistry(&error))?;
                local.insert(row, id);
            }
            for bond in &bonds_by_component[component] {
                let [a, b] = bond.atoms.map(|atom| local[&atom.index()]);
                editor
                    .add_bond(a, b, bond.order)
                    .map_err(|error| chemistry(&error))?;
            }
            let molecule = editor.finish().map_err(|error| chemistry(&error))?;
            let positions = Positions::from_vec(Quantity::new(
                rows.iter().map(|&row| self.points[row]).collect(),
                unit,
            ))
            .map_err(|error| chemistry(&error))?;
            let atom_ids: Vec<AtomId> = rows.iter().map(|row| local[row]).collect();
            let instance = builder
                .add_molecule(molecule, &positions)
                .map_err(|error| chemistry(&error))?;
            let sequence = i32::try_from(component + 1).map_err(|error| chemistry(&error))?;
            let hierarchy = builder.hierarchy_mut();
            let residue = hierarchy
                .add_residue(
                    chain,
                    "UNL",
                    Some(sequence),
                    Some(sequence.to_string()),
                    None,
                )
                .map_err(|error| chemistry(&error))?;
            for (&row, atom) in rows.iter().zip(atom_ids) {
                let atom = InstanceAtomId::new(instance, atom);
                hierarchy
                    .add_atom_site(residue, atom, AtomSiteMetadata::default())
                    .map_err(|error| chemistry(&error))?;
                order[row] = Some(atom);
            }
        }
        builder
            .set_atom_order(
                order
                    .into_iter()
                    .map(|atom| atom.expect("every atom line belongs to one component")),
            )
            .map_err(|error| chemistry(&error))?;
        builder.build().map_err(|error| chemistry(&error))
    }
}

fn dense_index(row: usize) -> TopologyAtomIndex {
    TopologyAtomIndex::new(u32::try_from(row).expect("atom lines are bounded below u32::MAX"))
}

/// Atom lines grouped into connected components, each in line order, ordered
/// by their first line.
fn components(atom_count: usize, bonds: &[XyzInferredBond]) -> Vec<Vec<usize>> {
    let mut parent = (0..atom_count).collect::<Vec<_>>();
    fn root(parent: &mut [usize], mut atom: usize) -> usize {
        while parent[atom] != atom {
            parent[atom] = parent[parent[atom]];
            atom = parent[atom];
        }
        atom
    }
    for bond in bonds {
        let [a, b] = bond.atoms.map(TopologyAtomIndex::index);
        let (a, b) = (root(&mut parent, a), root(&mut parent, b));
        if a != b {
            parent[a.max(b)] = a.min(b);
        }
    }
    let mut slot = HashMap::new();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for atom in 0..atom_count {
        let root = root(&mut parent, atom);
        let index = *slot.entry(root).or_insert_with(|| {
            groups.push(Vec::new());
            groups.len() - 1
        });
        groups[index].push(atom);
    }
    groups
}

/// Bonds and too-close contacts between atoms within covalent bonding
/// distance, found through a uniform grid so the work stays linear in the
/// atom count. Fails with the limit and the offending atom line when an atom
/// has more than `max_contacts` contacts.
fn infer_contacts(
    elements: &[Element],
    points: &[Point3],
    max_contacts: usize,
) -> Result<(Vec<XyzInferredBond>, Vec<XyzInterpretIssue>), (LimitExceeded, usize)> {
    let largest_radius = elements
        .iter()
        .map(|&element| covalent_bonding_radius(element).0)
        .fold(0.0_f64, f64::max);
    let cell = 2.0 * largest_radius + COVALENT_BOND_TOLERANCE_ANGSTROM;
    let key = |point: &Point3| {
        // Coordinates are bounded by `MAX_COORDINATE_ANGSTROM`, so every key
        // is exact and distant atoms never share a cell; a cell then holds
        // few atoms that are not contacts, which keeps the work linear.
        (
            (point.x / cell).floor() as i64,
            (point.y / cell).floor() as i64,
            (point.z / cell).floor() as i64,
        )
    };
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (atom, point) in points.iter().enumerate() {
        grid.entry(key(point)).or_default().push(atom);
    }
    let mut contacts = vec![0usize; points.len()];
    let mut bonds = Vec::new();
    let mut close = Vec::new();
    let minimum_squared = MIN_COVALENT_BOND_DISTANCE_ANGSTROM.powi(2);
    for (left, point) in points.iter().enumerate() {
        let (x, y, z) = key(point);
        let mut neighbours = BTreeSet::new();
        for dx in -1..=1_i64 {
            for dy in -1..=1_i64 {
                for dz in -1..=1_i64 {
                    neighbours.insert((
                        x.saturating_add(dx),
                        y.saturating_add(dy),
                        z.saturating_add(dz),
                    ));
                }
            }
        }
        for neighbour in neighbours {
            let Some(atoms) = grid.get(&neighbour) else {
                continue;
            };
            for &right in atoms.iter().filter(|&&right| right > left) {
                let other = points[right];
                let distance_squared = (point.x - other.x).powi(2)
                    + (point.y - other.y).powi(2)
                    + (point.z - other.z).powi(2);
                let cutoff = covalent_bonding_radius(elements[left]).0
                    + covalent_bonding_radius(elements[right]).0
                    + COVALENT_BOND_TOLERANCE_ANGSTROM;
                if distance_squared > cutoff * cutoff {
                    continue;
                }
                for atom in [left, right] {
                    contacts[atom] += 1;
                    if contacts[atom] > max_contacts {
                        return Err((
                            LimitExceeded::new(
                                Resource::CandidatePairs,
                                max_contacts,
                                contacts[atom],
                            ),
                            atom,
                        ));
                    }
                }
                let atoms = [dense_index(left), dense_index(right)];
                let distance = distance_squared.sqrt();
                if distance_squared <= minimum_squared {
                    close.push((atoms, distance));
                } else {
                    let order = if is_metal(elements[left]) != is_metal(elements[right]) {
                        BondOrder::Dative
                    } else {
                        BondOrder::Single
                    };
                    bonds.push(XyzInferredBond {
                        atoms,
                        order,
                        distance,
                        cutoff,
                    });
                }
            }
        }
    }
    bonds.sort_by_key(|bond| bond.atoms);
    close.sort_by_key(|&(atoms, _)| atoms);
    let close = close
        .into_iter()
        .map(|(atoms, distance)| XyzInterpretIssue::CloseContact {
            atoms,
            distance: Quantity::new(distance, ANGSTROM),
        })
        .collect();
    Ok((bonds, close))
}
