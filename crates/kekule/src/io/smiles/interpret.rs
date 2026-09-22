use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::Range;

use crate::chemistry::{
    canonicalize_molecule_for_publication, localize_source_aromatic_bonds, NormalizationError,
    SourceStereoBondMark, SourceStereoBondMarkKind,
};
use crate::core::{
    Atom, AtomId, AtomRadical, BondId, BondOrder, Element, HydrogenDeclaration, Molecule,
    MoleculeEditor, StereoCarrier, StereoElement, StereoElementId, StereoElementKind,
    TetrahedralOrientation, TetrahedralStereo,
};
use crate::topology::{Topology, TopologyBuildError};

use super::cx::{CxFeatures, CxStereoGroup};
use super::parse::{
    PendingStereoCarrier, PendingTetrahedral, SmilesAtomSyntax, SmilesBondToken,
    SmilesChiralityToken, SmilesDirectionToken, SmilesDocument, SmilesProgram, SmilesStereoCarrier,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesInterpretError {
    pub(super) offset: usize,
    pub(super) message: String,
}

impl SmilesInterpretError {
    pub const fn offset(&self) -> usize {
        self.offset
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for SmilesInterpretError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "SMILES interpretation error at {}: {}",
            self.offset, self.message
        )
    }
}

impl std::error::Error for SmilesInterpretError {}

/// A single-molecule accessor was used for a component-aware interpretation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmilesComponentCountError {
    actual: usize,
}

impl SmilesComponentCountError {
    pub const fn actual(self) -> usize {
        self.actual
    }
}

impl fmt::Display for SmilesComponentCountError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "single-molecule SMILES access requires exactly one component, found {}",
            self.actual
        )
    }
}

impl std::error::Error for SmilesComponentCountError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesAtomMapping {
    atom: AtomId,
    source_index: usize,
    source_span: Range<usize>,
}

impl SmilesAtomMapping {
    /// Zero-based atom index in the complete source record, before partitioning.
    /// CXSMILES atom references use this index, not a component-local `AtomId`.
    pub const fn source_index(&self) -> usize {
        self.source_index
    }

    pub const fn atom(&self) -> AtomId {
        self.atom
    }

    pub fn source_span(&self) -> Range<usize> {
        self.source_span.clone()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmilesBondMapping {
    bond: BondId,
    source_offset: usize,
}

impl SmilesBondMapping {
    pub const fn bond(&self) -> BondId {
        self.bond
    }

    pub const fn source_offset(&self) -> usize {
        self.source_offset
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SmilesInterpretationReport {
    atom_mappings: Vec<SmilesAtomMapping>,
    bond_mappings: Vec<SmilesBondMapping>,
    created_stereo_elements: Vec<StereoElementId>,
}

impl SmilesInterpretationReport {
    pub fn atom_mappings(&self) -> &[SmilesAtomMapping] {
        &self.atom_mappings
    }

    pub fn bond_mappings(&self) -> &[SmilesBondMapping] {
        &self.bond_mappings
    }

    /// Canonical stereo elements decoded from directional source bond syntax.
    pub fn created_stereo_elements(&self) -> &[StereoElementId] {
        &self.created_stereo_elements
    }
}

/// Interpretation of one connected SMILES component.
#[derive(Debug, Clone, PartialEq)]
pub struct SmilesComponentInterpretation {
    source_span: Range<usize>,
    molecule: Molecule,
    report: SmilesInterpretationReport,
}

impl SmilesComponentInterpretation {
    /// Bounding span of the source fragments containing this component.
    /// Spans can overlap when branches interleave disconnected components.
    pub fn source_span(&self) -> Range<usize> {
        self.source_span.clone()
    }

    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }

    pub fn report(&self) -> &SmilesInterpretationReport {
        &self.report
    }

    pub fn into_molecule(self) -> Molecule {
        self.molecule
    }

    pub fn into_parts(self) -> (Molecule, SmilesInterpretationReport) {
        (self.molecule, self.report)
    }
}

/// Canonical interpretation of one SMILES document as connected molecules.
#[derive(Debug, Clone, PartialEq)]
pub struct SmilesInterpretation {
    components: Vec<SmilesComponentInterpretation>,
    name: Option<String>,
    cx_extension: Option<String>,
    omitted_cx_extension_span: Option<Range<usize>>,
}

impl SmilesInterpretation {
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// Original extension retained as format metadata, including vertical bars.
    pub fn cx_extension(&self) -> Option<&str> {
        self.cx_extension.as_deref()
    }

    /// The extension omitted by an explicit base-SMILES projection, if any.
    /// A normal, complete interpretation never silently omits an extension.
    pub fn omitted_cx_extension_span(&self) -> Option<Range<usize>> {
        self.omitted_cx_extension_span.clone()
    }

    pub fn components(&self) -> &[SmilesComponentInterpretation] {
        &self.components
    }

    pub fn molecules(&self) -> impl ExactSizeIterator<Item = &Molecule> + DoubleEndedIterator {
        self.components
            .iter()
            .map(SmilesComponentInterpretation::molecule)
    }

    pub fn into_molecules(self) -> Vec<Molecule> {
        self.components
            .into_iter()
            .map(SmilesComponentInterpretation::into_molecule)
            .collect()
    }

    /// Projects every source component into one topology occurrence.
    ///
    /// Component order becomes authoritative topology instance order. The
    /// projection fabricates no hierarchy and runs no perception.
    pub fn into_topology(self) -> Result<Topology, TopologyBuildError> {
        Topology::from_molecules(&self.into_molecules())
    }

    /// Convenience access for callers that require exactly one component.
    ///
    /// Prefer [`Self::components`] for general SMILES input. An input with
    /// several connected components returns a component-count error.
    pub fn molecule(&self) -> Result<&Molecule, SmilesComponentCountError> {
        Ok(self.single_component()?.molecule())
    }

    /// Convenience report access for an interpretation known to contain one component.
    pub fn report(&self) -> Result<&SmilesInterpretationReport, SmilesComponentCountError> {
        Ok(self.single_component()?.report())
    }

    /// Consumes an interpretation known to contain exactly one component.
    pub fn into_molecule(self) -> Result<Molecule, SmilesComponentCountError> {
        Ok(self.into_single_component()?.into_molecule())
    }

    /// Consumes an interpretation known to contain exactly one component and its report.
    pub fn into_parts(
        self,
    ) -> Result<(Molecule, SmilesInterpretationReport), SmilesComponentCountError> {
        Ok(self.into_single_component()?.into_parts())
    }

    fn single_component(
        &self,
    ) -> Result<&SmilesComponentInterpretation, SmilesComponentCountError> {
        match self.components.as_slice() {
            [component] => Ok(component),
            components => Err(SmilesComponentCountError {
                actual: components.len(),
            }),
        }
    }

    fn into_single_component(
        mut self,
    ) -> Result<SmilesComponentInterpretation, SmilesComponentCountError> {
        if self.components.len() != 1 {
            return Err(SmilesComponentCountError {
                actual: self.components.len(),
            });
        }
        Ok(self
            .components
            .pop()
            .expect("length was checked to contain one SMILES component"))
    }
}

/// Interprets each connected SMILES component independently.
///
/// Parsing remains record-level: [`SmilesDocument`] preserves the complete
/// source and fragment separators. Ring closures may join dot-separated
/// fragments. Interpretation turns each connected component into one [`Molecule`] with component-local atom and
/// bond identifiers while retaining mappings to the original source offsets.
pub fn interpret_smiles_document(
    document: &SmilesDocument,
) -> Result<SmilesInterpretation, SmilesInterpretError> {
    interpret_document(document, false)
}

pub(super) fn interpret_base_smiles_document(
    document: &SmilesDocument,
) -> Result<SmilesInterpretation, SmilesInterpretError> {
    interpret_document(document, true)
}

fn interpret_document(
    document: &SmilesDocument,
    base_only: bool,
) -> Result<SmilesInterpretation, SmilesInterpretError> {
    let cx = if base_only {
        CxFeatures::default()
    } else {
        CxFeatures::interpret(document)?
    };
    let component_count = document
        .program
        .atoms
        .iter()
        .map(|atom| atom.component)
        .max()
        .map_or(0, |last| last + 1);
    let fragment_spans = document
        .fragment_token_ranges()
        .iter()
        .map(|range| fragment_source_span(document, range.clone()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut source_spans = vec![document.source().len()..0; component_count];
    // Partition once. Dot-separated records can contain as many components as
    // atoms, so rescanning the complete program per component is quadratic.
    let mut programs = (0..component_count)
        .map(|_| ComponentProgram::default())
        .collect::<Vec<_>>();
    let mut fragment = 0;
    for (index, atom) in document.program.atoms.iter().enumerate() {
        programs[atom.component].atoms.push(index);
        while fragment_spans[fragment].end <= atom.span.start {
            fragment += 1;
        }
        let span = &mut source_spans[atom.component];
        span.start = span.start.min(fragment_spans[fragment].start);
        span.end = span.end.max(fragment_spans[fragment].end);
    }
    for (index, bond) in document.program.bonds.iter().enumerate() {
        programs[bond.component].bonds.push(index);
    }
    for pending in &document.program.tetrahedral {
        programs[document.program.atoms[pending.center].component]
            .tetrahedral
            .push(*pending);
    }
    for group in &cx.groups {
        let mut atoms_by_component = BTreeMap::<usize, Vec<usize>>::new();
        for &atom in &group.atoms {
            atoms_by_component
                .entry(document.program.atoms[atom].component)
                .or_default()
                .push(atom);
        }
        for (component, atoms) in atoms_by_component {
            programs[component].groups.push(CxStereoGroup {
                kind: group.kind,
                atoms,
                offset: group.offset,
            });
        }
    }
    let mut components = Vec::with_capacity(component_count);
    for (program, source_span) in programs.iter().zip(source_spans) {
        let (molecule, report) = interpret_smiles_program_component(
            &document.program,
            program,
            document.source(),
            &cx.radicals,
        )?;
        components.push(SmilesComponentInterpretation {
            source_span,
            molecule,
            report,
        });
    }
    Ok(SmilesInterpretation {
        components,
        name: document.name().map(str::to_owned),
        cx_extension: document.cx_extension().map(str::to_owned),
        omitted_cx_extension_span: if base_only {
            document.cx_extension_span()
        } else {
            None
        },
    })
}

fn fragment_source_span(
    document: &SmilesDocument,
    token_range: Range<usize>,
) -> Result<Range<usize>, SmilesInterpretError> {
    if token_range.is_empty() {
        return Err(SmilesInterpretError {
            offset: document.source().len(),
            message: "empty SMILES component".to_owned(),
        });
    }
    let first = document
        .tokens()
        .get(token_range.start)
        .ok_or_else(|| SmilesInterpretError {
            offset: document.source().len(),
            message: "component token range starts outside the SMILES document".to_owned(),
        })?;
    let last = document
        .tokens()
        .get(token_range.end.saturating_sub(1))
        .ok_or_else(|| SmilesInterpretError {
            offset: document.source().len(),
            message: "component token range ends outside the SMILES document".to_owned(),
        })?;
    Ok(first.span().start..last.span().end)
}

#[derive(Default)]
struct ComponentProgram {
    atoms: Vec<usize>,
    bonds: Vec<usize>,
    tetrahedral: Vec<PendingTetrahedral>,
    groups: Vec<CxStereoGroup>,
}

fn interpret_smiles_program_component(
    program: &SmilesProgram,
    component: &ComponentProgram,
    source: &str,
    radicals: &BTreeMap<usize, AtomRadical>,
) -> std::result::Result<(Molecule, SmilesInterpretationReport), SmilesInterpretError> {
    validate_smiles_source_aromaticity(program, component, source)?;
    let end_offset = source.len();

    let mut editor = crate::core::MoleculeEditor::new();
    let mut source_to_atom = BTreeMap::<usize, AtomId>::new();
    let mut atom_mappings = Vec::new();
    for &index in &component.atoms {
        let record = &program.atoms[index];
        let atom = interpret_smiles_atom(&record.syntax, record.span.start)?;
        let atom_id = editor
            .add_atom(atom)
            .map_err(|error| SmilesInterpretError {
                offset: record.span.start,
                message: format!("invalid represented atom: {error}"),
            })?;
        source_to_atom.insert(index, atom_id);
        atom_mappings.push(SmilesAtomMapping {
            atom: atom_id,
            source_index: index,
            source_span: record.span.clone(),
        });
    }
    let mut bond_mappings = Vec::new();
    let mut source_aromatic_bonds = BTreeSet::new();
    let mut source_stereo = Vec::new();
    let mut first_aromatic_offset = None;
    for &index in &component.bonds {
        let bond = &program.bonds[index];
        let left = source_to_atom
            .get(&bond.left)
            .copied()
            .ok_or_else(|| SmilesInterpretError {
                offset: bond.offset,
                message: "bond left endpoint is outside its SMILES component".to_owned(),
            })?;
        let right =
            source_to_atom
                .get(&bond.right)
                .copied()
                .ok_or_else(|| SmilesInterpretError {
                    offset: bond.offset,
                    message: "bond right endpoint is outside its SMILES component".to_owned(),
                })?;
        let (order, source_aromatic) = interpret_smiles_bond_token(bond.token);
        let bond_id = add_smiles_bond(&mut editor, left, right, order, bond.offset)?;
        if let Some(direction) = bond.direction {
            let source_from = bond.direction_from.ok_or_else(|| SmilesInterpretError {
                offset: bond.direction_offset.unwrap_or(bond.offset),
                message: "directional bond is missing its textual origin endpoint".to_owned(),
            })?;
            let from =
                source_to_atom
                    .get(&source_from)
                    .copied()
                    .ok_or_else(|| SmilesInterpretError {
                        offset: bond.direction_offset.unwrap_or(bond.offset),
                        message: "directional bond origin is outside its SMILES component"
                            .to_owned(),
                    })?;
            source_stereo.push(SourceStereoBondMark {
                bond: bond_id,
                from,
                kind: interpret_smiles_direction(direction),
            });
        }
        if source_aromatic {
            source_aromatic_bonds.insert(bond_id);
            first_aromatic_offset.get_or_insert(bond.offset);
        }
        bond_mappings.push(SmilesBondMapping {
            bond: bond_id,
            source_offset: bond.direction_offset.unwrap_or(bond.offset),
        });
    }
    localize_source_aromatic_bonds(editor.working_mut(), &source_aromatic_bonds).map_err(
        |error| SmilesInterpretError {
            offset: first_aromatic_offset.unwrap_or(end_offset),
            message: error.to_string(),
        },
    )?;

    for (&source_index, &atom_id) in &source_to_atom {
        if let Some(radical) = radicals.get(&source_index) {
            editor
                .atom_mut(atom_id)
                .expect("source atom is live")
                .radical = Some(*radical);
            continue;
        }
        let molecule = editor.working();
        let atom = molecule.atom(atom_id).expect("source atom is live");
        if let HydrogenDeclaration::Fixed(hydrogens) = atom.hydrogens {
            let electrons = crate::algorithms::rdkit_bracket_radical_electrons(
                atom,
                crate::algorithms::explicit_valence(molecule, atom_id),
                molecule
                    .incident_bonds(atom_id)
                    .expect("source atom is live")
                    .count(),
                hydrogens,
            );
            editor
                .atom_mut(atom_id)
                .expect("source atom is live")
                .radical = AtomRadical::new(electrons, None);
        }
    }

    add_smiles_tetrahedral_elements(
        &mut editor,
        &source_to_atom,
        &component.tetrahedral,
        &program.tetrahedral_carriers,
        end_offset,
    )?;
    let publication_report = canonicalize_molecule_for_publication(
        editor.working_mut(),
        None,
        &source_stereo,
    )
    .map_err(|error| SmilesInterpretError {
        offset: canonicalization_error_offset(&error, &atom_mappings, &bond_mappings, end_offset),
        message: format!("could not publish canonical molecule: {error}"),
    })?;
    debug_assert!(publication_report.warnings.is_empty());
    super::cx::install_stereo_groups(&mut editor, &source_to_atom, &component.groups)?;
    let molecule = editor.finish().map_err(|error| SmilesInterpretError {
        offset: atom_mappings
            .first()
            .map_or(end_offset, |mapping| mapping.source_span.start),
        message: error.to_string(),
    })?;
    Ok((
        molecule,
        SmilesInterpretationReport {
            atom_mappings,
            bond_mappings,
            created_stereo_elements: publication_report.created_stereo_elements,
        },
    ))
}

fn canonicalization_error_offset(
    error: &NormalizationError,
    atom_mappings: &[SmilesAtomMapping],
    bond_mappings: &[SmilesBondMapping],
    fallback: usize,
) -> usize {
    error
        .bond_location_hint()
        .and_then(|bond| {
            bond_mappings
                .iter()
                .find(|mapping| mapping.bond == bond)
                .map(|mapping| mapping.source_offset)
        })
        .or_else(|| {
            error.atom_location_hint().and_then(|atom| {
                atom_mappings
                    .iter()
                    .find(|mapping| mapping.atom == atom)
                    .map(|mapping| mapping.source_span.start)
            })
        })
        .or_else(|| {
            atom_mappings
                .first()
                .map(|mapping| mapping.source_span.start)
        })
        .unwrap_or(fallback)
}

/// Require atom- and bond-level source aromaticity assertions to cover the
/// same staged subgraph before any canonical molecule is constructed.
fn validate_smiles_source_aromaticity(
    program: &SmilesProgram,
    component: &ComponentProgram,
    source: &str,
) -> std::result::Result<(), SmilesInterpretError> {
    let mut source_aromatic_atoms = BTreeSet::new();
    for &index in &component.atoms {
        let record = &program.atoms[index];
        let imported = program.imported_aromatic_atoms.contains(&index);
        if record.syntax.aromatic != imported {
            return Err(SmilesInterpretError {
                offset: record.span.start,
                message: "inconsistent aromatic atom syntax state".to_owned(),
            });
        }
        if imported {
            source_aromatic_atoms.insert(index);
        }
    }

    let mut atoms_with_source_aromatic_bonds = BTreeSet::new();
    for bond in component
        .bonds
        .iter()
        .map(|&index| &program.bonds[index])
        .filter(|bond| bond.token == SmilesBondToken::Aromatic)
    {
        if !source_aromatic_atoms.contains(&bond.left)
            || !source_aromatic_atoms.contains(&bond.right)
        {
            return Err(SmilesInterpretError {
                offset: explicit_aromatic_bond_offset(source, bond.offset),
                message:
                    "source-aromatic bond requires source-aromatic atom syntax at both endpoints"
                        .to_owned(),
            });
        }
        atoms_with_source_aromatic_bonds.insert(bond.left);
        atoms_with_source_aromatic_bonds.insert(bond.right);
    }

    if let Some(index) = source_aromatic_atoms
        .difference(&atoms_with_source_aromatic_bonds)
        .next()
        .copied()
    {
        let record = &program.atoms[index];
        return Err(SmilesInterpretError {
            offset: record.span.start,
            message: "source-aromatic atom is not part of a source-aromatic bond".to_owned(),
        });
    }

    Ok(())
}

fn explicit_aromatic_bond_offset(source: &str, fallback: usize) -> usize {
    source
        .get(..fallback)
        .and_then(|prefix| prefix.char_indices().next_back())
        .filter(|(_, token)| *token == ':')
        .map_or(fallback, |(offset, _)| offset)
}

fn interpret_smiles_atom(
    syntax: &SmilesAtomSyntax,
    offset: usize,
) -> std::result::Result<Atom, SmilesInterpretError> {
    let element = Element::from_symbol(&syntax.symbol).ok_or_else(|| SmilesInterpretError {
        offset,
        message: format!("unsupported element symbol `{}`", syntax.symbol),
    })?;
    let mut atom = Atom::new(element);
    atom.isotope = syntax.isotope;
    atom.formal_charge = syntax.formal_charge;
    atom.hydrogens = if syntax.bracketed {
        HydrogenDeclaration::Fixed(syntax.specified_hydrogens)
    } else {
        HydrogenDeclaration::Infer {
            specified: syntax.specified_hydrogens,
        }
    };
    atom.atom_map = syntax.atom_map;
    Ok(atom)
}

const fn interpret_smiles_bond_token(token: SmilesBondToken) -> (BondOrder, bool) {
    match token {
        SmilesBondToken::Single => (BondOrder::Single, false),
        SmilesBondToken::Double => (BondOrder::Double, false),
        SmilesBondToken::Triple => (BondOrder::Triple, false),
        SmilesBondToken::Quadruple => (BondOrder::Quadruple, false),
        SmilesBondToken::Aromatic => (BondOrder::Single, true),
    }
}

const fn interpret_smiles_direction(direction: SmilesDirectionToken) -> SourceStereoBondMarkKind {
    match direction {
        SmilesDirectionToken::Up => SourceStereoBondMarkKind::DirectionalUp,
        SmilesDirectionToken::Down => SourceStereoBondMarkKind::DirectionalDown,
    }
}

fn add_smiles_bond(
    editor: &mut MoleculeEditor,
    left: AtomId,
    right: AtomId,
    order: BondOrder,
    offset: usize,
) -> std::result::Result<BondId, SmilesInterpretError> {
    editor
        .add_bond(left, right, order)
        .map_err(|error| SmilesInterpretError {
            offset,
            message: error.to_string(),
        })
}

fn add_smiles_tetrahedral_elements(
    editor: &mut MoleculeEditor,
    source_to_atom: &BTreeMap<usize, AtomId>,
    centers: &[PendingTetrahedral],
    carriers_by_center: &BTreeMap<usize, Vec<PendingStereoCarrier>>,
    offset: usize,
) -> std::result::Result<(), SmilesInterpretError> {
    for pending in centers {
        let Some(&center) = source_to_atom.get(&pending.center) else {
            continue;
        };
        let carriers = resolve_smiles_tetrahedral_carriers(
            editor.working(),
            center,
            source_to_atom,
            carriers_by_center
                .get(&pending.center)
                .cloned()
                .unwrap_or_default(),
            offset,
        )?;
        editor
            .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
                TetrahedralStereo {
                    center,
                    carriers,
                    orientation: Some(match pending.orientation {
                        SmilesChiralityToken::At => TetrahedralOrientation::Clockwise,
                        SmilesChiralityToken::AtAt => TetrahedralOrientation::CounterClockwise,
                    }),
                },
            )))
            .map_err(|error| SmilesInterpretError {
                offset,
                message: error.to_string(),
            })?;
    }
    Ok(())
}

fn resolve_smiles_tetrahedral_carriers(
    mol: &Molecule,
    center: AtomId,
    source_to_atom: &BTreeMap<usize, AtomId>,
    carriers: Vec<PendingStereoCarrier>,
    offset: usize,
) -> std::result::Result<Vec<StereoCarrier>, SmilesInterpretError> {
    let mut carriers = carriers
        .into_iter()
        .map(|carrier| match carrier {
            PendingStereoCarrier::Resolved(SmilesStereoCarrier::Atom(source)) => source_to_atom
                .get(&source)
                .copied()
                .map(StereoCarrier::Atom)
                .ok_or_else(|| SmilesInterpretError {
                    offset,
                    message: "tetrahedral carrier is outside its SMILES component".to_owned(),
                }),
            PendingStereoCarrier::Resolved(SmilesStereoCarrier::ImplicitHydrogen) => {
                Ok(StereoCarrier::ImplicitHydrogen)
            }
            PendingStereoCarrier::Ring { .. } => Err(SmilesInterpretError {
                offset,
                message: "unresolved tetrahedral ring carrier".to_owned(),
            }),
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    if carriers.len() == 3 && smiles_tetrahedral_center_can_have_lone_pair(mol, center) {
        // Trigonal-pyramidal SMILES with a bracket H use explicit neighbors,
        // then H, then the phantom lone-pair carrier (RDKit's convention).
        if let Some(index) = carriers
            .iter()
            .position(|carrier| *carrier == StereoCarrier::ImplicitHydrogen)
        {
            let hydrogen = carriers.remove(index);
            carriers.push(hydrogen);
        }
        carriers.push(StereoCarrier::ImplicitLonePair);
    }
    Ok(carriers)
}

fn smiles_tetrahedral_center_can_have_lone_pair(mol: &Molecule, center: AtomId) -> bool {
    mol.atom(center)
        .map(|atom| {
            matches!(
                atom.element.symbol(),
                "N" | "P" | "As" | "Sb" | "O" | "S" | "Se" | "Te"
            )
        })
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::super::parse::parse_smiles_document;
    use super::*;

    #[test]
    fn dot_smiles_interprets_as_connected_molecules() {
        let document = parse_smiles_document("CC(=O)[O-].[Na+]").expect("valid salt");
        let interpretation = interpret_smiles_document(&document).expect("interpret salt");
        assert_eq!(interpretation.components().len(), 2);
        assert!(interpretation
            .molecules()
            .all(|molecule| molecule.is_connected()));
        assert_eq!(interpretation.components()[0].source_span(), 0..10);
        assert_eq!(interpretation.components()[1].source_span(), 11..16);
    }

    #[test]
    fn component_mappings_retain_document_offsets() {
        let document = parse_smiles_document("C.[Na+]").expect("valid components");
        let interpretation = interpret_smiles_document(&document).expect("interpret components");
        assert_eq!(
            interpretation.components()[1].report().atom_mappings()[0].source_span(),
            2..7
        );
    }

    #[test]
    fn single_component_convenience_rejects_dot_smiles_without_panicking() {
        let document = parse_smiles_document("C.O").expect("valid components");
        let error = interpret_smiles_document(&document)
            .expect("interpret components")
            .into_molecule()
            .expect_err("single-component access must reject two components");
        assert_eq!(error.actual(), 2);
    }
}
