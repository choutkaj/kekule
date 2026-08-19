use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::Range;

use crate::chemistry::{
    canonicalize_molecule_for_publication, localize_source_aromatic_bonds, NormalizationError,
};
use crate::core::{
    Atom, AtomId, BondId, BondOrder, Element, Molecule, StereoBondMark, StereoBondMarkKind,
    StereoCarrier, StereoElement, StereoElementId, StereoElementKind, StereoSource,
    TetrahedralOrientation, TetrahedralStereo,
};
use crate::small::model::SmallMolecule;

use super::parse::{
    PendingStereoCarrier, PendingTetrahedral, SmilesAtomSyntax, SmilesBondToken,
    SmilesChiralityToken, SmilesDirectionToken, SmilesDocument, SmilesProgram, SmilesStereoCarrier,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesInterpretError {
    offset: usize,
    message: String,
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
    source_span: Range<usize>,
}

impl SmilesAtomMapping {
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

/// Interpretation of one dot-delimited connected SMILES component.
#[derive(Debug, Clone, PartialEq)]
pub struct SmilesComponentInterpretation {
    source_span: Range<usize>,
    molecule: SmallMolecule,
    report: SmilesInterpretationReport,
}

impl SmilesComponentInterpretation {
    pub fn source_span(&self) -> Range<usize> {
        self.source_span.clone()
    }

    pub fn molecule(&self) -> &SmallMolecule {
        &self.molecule
    }

    pub fn report(&self) -> &SmilesInterpretationReport {
        &self.report
    }

    pub fn into_molecule(self) -> SmallMolecule {
        self.molecule
    }

    pub fn into_parts(self) -> (SmallMolecule, SmilesInterpretationReport) {
        (self.molecule, self.report)
    }
}

/// Canonical interpretation of one SMILES document as connected molecules.
#[derive(Debug, Clone, PartialEq)]
pub struct SmilesInterpretation {
    components: Vec<SmilesComponentInterpretation>,
}

impl SmilesInterpretation {
    pub fn components(&self) -> &[SmilesComponentInterpretation] {
        &self.components
    }

    pub fn molecules(&self) -> impl ExactSizeIterator<Item = &SmallMolecule> + DoubleEndedIterator {
        self.components
            .iter()
            .map(SmilesComponentInterpretation::molecule)
    }

    pub fn into_molecules(self) -> Vec<SmallMolecule> {
        self.components
            .into_iter()
            .map(SmilesComponentInterpretation::into_molecule)
            .collect()
    }

    /// Convenience access for callers that require exactly one component.
    ///
    /// Prefer [`Self::components`] for general SMILES input. This method keeps
    /// the pre-component API convenient for callers whose input contract is
    /// already single-molecule and fails loudly rather than discarding data.
    pub fn molecule(&self) -> Result<&SmallMolecule, SmilesComponentCountError> {
        Ok(self.single_component()?.molecule())
    }

    /// Convenience report access for an interpretation known to contain one component.
    pub fn report(&self) -> Result<&SmilesInterpretationReport, SmilesComponentCountError> {
        Ok(self.single_component()?.report())
    }

    /// Consumes an interpretation known to contain exactly one component.
    pub fn into_molecule(self) -> Result<SmallMolecule, SmilesComponentCountError> {
        Ok(self.into_single_component()?.into_molecule())
    }

    /// Consumes an interpretation known to contain exactly one component and its report.
    pub fn into_parts(
        self,
    ) -> Result<(SmallMolecule, SmilesInterpretationReport), SmilesComponentCountError> {
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

/// Interprets each dot-delimited SMILES component independently.
///
/// Parsing remains record-level: [`SmilesDocument`] preserves the complete
/// source and component separators. Interpretation turns each syntactic
/// component into one connected [`SmallMolecule`] with component-local atom and
/// bond identifiers while retaining mappings to the original source offsets.
pub fn interpret_smiles_document(
    document: &SmilesDocument,
) -> Result<SmilesInterpretation, SmilesInterpretError> {
    let mut components = Vec::with_capacity(document.component_token_ranges().len());
    for (component_index, token_range) in document.component_token_ranges().iter().enumerate() {
        let source_span = component_source_span(document, token_range.clone())?;
        let local = interpret_smiles_component(document, component_index).map_err(|error| {
            SmilesInterpretError {
                offset: error.offset(),
                message: error.message().to_owned(),
            }
        })?;
        let (molecule, report) = local.into_parts();
        molecule
            .graph()
            .validate_connected()
            .map_err(|error| SmilesInterpretError {
                offset: source_span.start,
                message: error.to_string(),
            })?;

        let atom_mappings = report
            .atom_mappings()
            .iter()
            .map(|mapping| SmilesAtomMapping {
                atom: mapping.atom(),
                source_span: mapping.source_span(),
            })
            .collect();
        let bond_mappings = report
            .bond_mappings()
            .iter()
            .map(|mapping| SmilesBondMapping {
                bond: mapping.bond(),
                source_offset: mapping.source_offset(),
            })
            .collect();
        let created_stereo_elements = report.created_stereo_elements().to_vec();

        components.push(SmilesComponentInterpretation {
            source_span,
            molecule,
            report: SmilesInterpretationReport {
                atom_mappings,
                bond_mappings,
                created_stereo_elements,
            },
        });
    }
    Ok(SmilesInterpretation { components })
}

fn component_source_span(
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

#[derive(Debug, Clone, PartialEq)]
struct SmilesProgramInterpretation {
    molecule: SmallMolecule,
    report: SmilesInterpretationReport,
}

impl SmilesProgramInterpretation {
    fn into_parts(self) -> (SmallMolecule, SmilesInterpretationReport) {
        (self.molecule, self.report)
    }
}

fn interpret_smiles_component(
    document: &SmilesDocument,
    component: usize,
) -> std::result::Result<SmilesProgramInterpretation, SmilesInterpretError> {
    interpret_smiles_program_component(&document.program, component, document.source())
}

fn interpret_smiles_program_component(
    program: &SmilesProgram,
    component: usize,
    source: &str,
) -> std::result::Result<SmilesProgramInterpretation, SmilesInterpretError> {
    validate_smiles_source_aromaticity(program, component, source)?;
    let end_offset = source.len();

    let mut mol = Molecule::new();
    let mut source_to_atom = BTreeMap::<usize, AtomId>::new();
    let mut atom_mappings = Vec::new();
    for (index, record) in program.atoms.iter().enumerate() {
        if record.component != component {
            continue;
        }
        let atom = interpret_smiles_atom(&record.syntax, record.span.start)?;
        let atom_id = mol.add_atom(atom).map_err(|error| SmilesInterpretError {
            offset: record.span.start,
            message: format!("invalid represented atom: {error}"),
        })?;
        source_to_atom.insert(index, atom_id);
        atom_mappings.push(SmilesAtomMapping {
            atom: atom_id,
            source_span: record.span.clone(),
        });
    }
    let mut bond_mappings = Vec::new();
    let mut source_aromatic_bonds = BTreeSet::new();
    let mut first_aromatic_offset = None;
    for bond in &program.bonds {
        if bond.component != component {
            continue;
        }
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
        let bond_id = add_smiles_bond(
            &mut mol,
            left,
            right,
            order,
            bond.direction.map(interpret_smiles_direction),
            bond.offset,
        )?;
        if source_aromatic {
            source_aromatic_bonds.insert(bond_id);
            first_aromatic_offset.get_or_insert(bond.offset);
        }
        bond_mappings.push(SmilesBondMapping {
            bond: bond_id,
            source_offset: bond.offset,
        });
    }
    localize_source_aromatic_bonds(&mut mol, &source_aromatic_bonds).map_err(|error| {
        SmilesInterpretError {
            offset: first_aromatic_offset.unwrap_or(end_offset),
            message: error.to_string(),
        }
    })?;

    add_smiles_tetrahedral_elements(
        &mut mol,
        &source_to_atom,
        &program.tetrahedral,
        &program.tetrahedral_carriers,
        end_offset,
    )?;
    let publication_report =
        canonicalize_molecule_for_publication(&mut mol).map_err(|error| SmilesInterpretError {
            offset: canonicalization_error_offset(
                &error,
                &atom_mappings,
                &bond_mappings,
                end_offset,
            ),
            message: format!("could not publish canonical molecule: {error}"),
        })?;
    debug_assert!(publication_report.warnings.is_empty());
    Ok(SmilesProgramInterpretation {
        molecule: SmallMolecule::from_graph(mol),
        report: SmilesInterpretationReport {
            atom_mappings,
            bond_mappings,
            created_stereo_elements: publication_report.created_stereo_elements,
        },
    })
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
    component: usize,
    source: &str,
) -> std::result::Result<(), SmilesInterpretError> {
    let mut source_aromatic_atoms = BTreeSet::new();
    for (index, record) in program.atoms.iter().enumerate() {
        if record.component != component {
            continue;
        }
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
    for bond in program
        .bonds
        .iter()
        .filter(|bond| bond.component == component && bond.token == SmilesBondToken::Aromatic)
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
    atom.explicit_hydrogens = syntax.explicit_hydrogens;
    atom.no_implicit_hydrogens = syntax.bracketed;
    atom.atom_map = syntax.atom_map;
    Ok(atom)
}

const fn interpret_smiles_bond_token(token: SmilesBondToken) -> (BondOrder, bool) {
    match token {
        SmilesBondToken::Single => (BondOrder::Single, false),
        SmilesBondToken::Double => (BondOrder::Double, false),
        SmilesBondToken::Triple => (BondOrder::Triple, false),
        SmilesBondToken::Aromatic => (BondOrder::Single, true),
    }
}

const fn interpret_smiles_direction(direction: SmilesDirectionToken) -> StereoBondMarkKind {
    match direction {
        SmilesDirectionToken::Up => StereoBondMarkKind::DirectionalUp,
        SmilesDirectionToken::Down => StereoBondMarkKind::DirectionalDown,
    }
}

fn add_smiles_bond(
    mol: &mut Molecule,
    left: AtomId,
    right: AtomId,
    order: BondOrder,
    stereo: Option<StereoBondMarkKind>,
    offset: usize,
) -> std::result::Result<BondId, SmilesInterpretError> {
    let bond_id = mol
        .add_bond(left, right, order)
        .map_err(|error| SmilesInterpretError {
            offset,
            message: error.to_string(),
        })?;
    if let Some(kind) = stereo {
        mol.set_stereo_bond_mark(StereoBondMark {
            bond: bond_id,
            kind,
            source: StereoSource::Smiles,
        })
        .map_err(|error| SmilesInterpretError {
            offset,
            message: error.to_string(),
        })?;
    }
    Ok(bond_id)
}

fn add_smiles_tetrahedral_elements(
    mol: &mut Molecule,
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
            mol,
            center,
            source_to_atom,
            carriers_by_center
                .get(&pending.center)
                .cloned()
                .unwrap_or_default(),
            offset,
        )?;
        mol.add_stereo_element(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers,
                orientation: match pending.orientation {
                    SmilesChiralityToken::At => TetrahedralOrientation::Clockwise,
                    SmilesChiralityToken::AtAt => TetrahedralOrientation::CounterClockwise,
                },
            }),
            StereoSource::Smiles,
        ))
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
            ) && atom.explicit_hydrogens == 0
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
            .all(|molecule| molecule.graph().is_connected()));
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
