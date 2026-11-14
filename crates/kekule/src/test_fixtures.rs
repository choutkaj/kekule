//! Fixtures shared by unit tests in several modules: readers of test input,
//! small molecules, perception helpers, and source-backed molblocks.
//!
//! Tests import what they use by name; nothing here is glob-imported.

use crate::diagnostics::DiagnosticKind;
use crate::io::{
    canonicalize_molecule_for_publication, NormalizationError, NormalizationReport,
    SourceStereoBondMark,
};
use crate::math::Point3;
use crate::molecule::{
    AromaticityModel, Atom, AtomId, AtomRadical, BondId, BondOrder, Element, ImplicitHydrogens,
    Molecule, MoleculeEditor, Perception, StereoElement, StereoElementId, StereoElementKind,
    StereoGroup, StereoGroupId, ValenceModel,
};
use crate::perception::{rings as rings_api, rings::RingPerceptionOptions, valence as valence_api};
use crate::sdf::{SdfParseOptions, SdfRecordInterpretation};
use crate::smiles::SmilesInterpretationReport;
use crate::structure::{Positions, Structure, StructureSet};
use crate::{molfile, sdf, smiles as smiles_api};

/// The typed detail of a failed format operation.
pub(crate) fn error_detail<T>(
    result: Result<T, crate::diagnostics::FormatError>,
) -> Option<crate::diagnostics::FormatErrorDetail> {
    result.err()?.detail().cloned()
}

/// Each geometry of `set` as its atoms in topology order: the element symbol
/// and the position in Å, rounded to 1e-9 Å so that coordinates read from
/// text compare equal to their decimal literals.
pub(crate) fn geometry_atoms(set: &StructureSet) -> Vec<Vec<(&'static str, [f64; 3])>> {
    let angstroms = |nanometers: f64| (nanometers * 1.0e10).round() / 1.0e9;
    set.iter()
        .map(|geometry| {
            geometry
                .as_structure_view()
                .atoms()
                .map(|atom| {
                    let position = atom.position().into_value();
                    (
                        atom.element.symbol(),
                        [position.x, position.y, position.z].map(angstroms),
                    )
                })
                .collect()
        })
        .collect()
}

/// The messages of `error` and its sources, joined by `: `.
pub(crate) fn chain_text(error: &(dyn std::error::Error + 'static)) -> String {
    let mut text = error.to_string();
    let mut source = error.source();
    while let Some(cause) = source {
        text.push_str(": ");
        text.push_str(&cause.to_string());
        source = cause.source();
    }
    text
}

/// The mmCIF writer issue of a failed format operation.
pub(crate) fn mmcif_issue<T>(
    result: Result<T, crate::diagnostics::FormatError>,
) -> Option<crate::mmcif::MmcifWriteErrorKind> {
    match error_detail(result)? {
        crate::diagnostics::FormatErrorDetail::MmcifWrite(issue) => Some(issue),
        _ => None,
    }
}

/// The 0-based byte offset of a SMILES or SMARTS failure within its record.
pub(crate) fn byte_offset(error: &crate::diagnostics::FormatError) -> usize {
    error
        .location()
        .and_then(|location| location.column())
        .expect("SMILES and SMARTS failures carry a byte column")
        - 1
}

pub(crate) fn carbon() -> Atom {
    Atom::new(Element::from_symbol("C").expect("carbon should be available"))
}

pub(crate) fn oxygen() -> Atom {
    Atom::new(Element::from_symbol("O").expect("oxygen should be available"))
}

pub(crate) fn read_smiles(input: &str) -> Result<Molecule, Box<dyn std::error::Error>> {
    let document = smiles_api::parse_str(input, Default::default())?;
    Ok(document.interpret(Default::default())?.into_molecule()?)
}

pub(crate) fn read_smiles_components(
    input: &str,
) -> Result<Vec<Molecule>, Box<dyn std::error::Error>> {
    let document = smiles_api::parse_str(input, Default::default())?;
    Ok(document.interpret(Default::default())?.into_molecules())
}

pub(crate) fn read_smiles_component(
    input: &str,
    component: usize,
) -> Result<Molecule, Box<dyn std::error::Error>> {
    read_smiles_components(input)?
        .into_iter()
        .nth(component)
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("SMILES has no component at index {component}"),
            )
            .into()
        })
}

pub(crate) fn perceive(molecule: &mut Molecule) -> Result<(), Box<dyn std::error::Error>> {
    molecule.perceive()?;
    Ok(())
}

pub(crate) fn canonical_smiles_round_trip(molecule: &Molecule) -> (String, Molecule) {
    let written = smiles_api::write(molecule, smiles_api::SmilesWriteOptions::canonical())
        .unwrap_or_else(|error| panic!("canonical SMILES should write: {error}"));
    let mut reparsed = read_smiles(&written)
        .unwrap_or_else(|error| panic!("canonical output should parse: {written}: {error}"));
    perceive(&mut reparsed)
        .unwrap_or_else(|error| panic!("canonical output should perceive: {written}: {error}"));
    (written, reparsed)
}

pub(crate) fn read_smiles_with_report(
    input: &str,
) -> Result<(Molecule, SmilesInterpretationReport), Box<dyn std::error::Error>> {
    let document = smiles_api::parse_str(input, Default::default())?;
    Ok(document.interpret(Default::default())?.into_parts()?)
}

pub(crate) trait CanonicalizeFixture {
    fn canonicalize_fixture(&mut self) -> Result<NormalizationReport, NormalizationError>;

    fn canonicalize_fixture_with_source_stereo(
        &mut self,
        source_stereo: &[SourceStereoBondMark],
    ) -> Result<NormalizationReport, NormalizationError>;
}

/// Applies a structural edit to a published test molecule in place, without
/// publication's checks, as crate code once could on the molecule directly.
pub(crate) fn edit<R>(molecule: &mut Molecule, edit: impl FnOnce(&mut MoleculeEditor) -> R) -> R {
    let mut editor = molecule.clone().into_editor();
    let result = edit(&mut editor);
    *molecule = editor.into_unchecked();
    result
}

/// Runs `f` on a draft that has deleted nothing, through an unchecked
/// published copy whose IDs equal the draft's, keeping what `f` installs.
pub(crate) fn on_draft<R>(editor: &mut MoleculeEditor, f: impl FnOnce(&mut Molecule) -> R) -> R {
    assert!(!editor.has_deleted(), "draft IDs must survive publication");
    let mut molecule = std::mem::take(editor).into_unchecked();
    let result = f(&mut molecule);
    *editor = molecule.into_editor();
    result
}

impl CanonicalizeFixture for Molecule {
    fn canonicalize_fixture(&mut self) -> Result<NormalizationReport, NormalizationError> {
        edit(self, |editor| editor.canonicalize_fixture())
    }

    fn canonicalize_fixture_with_source_stereo(
        &mut self,
        source_stereo: &[SourceStereoBondMark],
    ) -> Result<NormalizationReport, NormalizationError> {
        edit(self, |editor| {
            editor.canonicalize_fixture_with_source_stereo(source_stereo)
        })
    }
}

impl CanonicalizeFixture for MoleculeEditor {
    fn canonicalize_fixture(&mut self) -> Result<NormalizationReport, NormalizationError> {
        canonicalize_molecule_for_publication(self, None, &[])
    }

    fn canonicalize_fixture_with_source_stereo(
        &mut self,
        source_stereo: &[SourceStereoBondMark],
    ) -> Result<NormalizationReport, NormalizationError> {
        canonicalize_molecule_for_publication(self, None, source_stereo)
    }
}

pub(crate) fn read_molfile(input: &str) -> Result<Molecule, Box<dyn std::error::Error>> {
    let document = molfile::parse_str(input, Default::default())?;
    let molecules = document.interpret(Default::default())?.into_molecules();
    exactly_one(molecules, "molfile")
}

/// The stage at which malformed input fails.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FailureStage {
    Parse,
    Interpret,
}

/// Where a malformed Molfile fails: its stage, diagnostic kind, and line.
pub(crate) fn molfile_failure(input: &str) -> (FailureStage, DiagnosticKind, Option<usize>) {
    match molfile::parse_str(input, Default::default()) {
        Err(error) => (
            FailureStage::Parse,
            error.kind(),
            Some(error.location().unwrap().line()),
        ),
        Ok(document) => {
            let error = document
                .interpret(Default::default())
                .expect_err("malformed Molfile must fail");
            (
                FailureStage::Interpret,
                error.kind(),
                error.location().map(|location| location.line()),
            )
        }
    }
}

pub(crate) fn read_molfile_with_report(
    input: &str,
) -> Result<(Molecule, molfile::MolfileInterpretationReport), Box<dyn std::error::Error>> {
    let document = molfile::parse_str(input, Default::default())?;
    let (structure, mut reports, _) = document.interpret(Default::default())?.into_parts();
    if reports.len() != 1 {
        return Err(format!("expected one molfile component, found {}", reports.len()).into());
    }
    let molecule = structure
        .topology()
        .instances()
        .next()
        .expect("one component")
        .molecule()
        .clone();
    Ok((molecule, reports.pop().expect("length checked")))
}

pub(crate) fn read_sdf_records(
    input: &str,
) -> Result<Vec<SdfRecordInterpretation>, Box<dyn std::error::Error>> {
    read_sdf_records_with_options(input, SdfParseOptions::default())
}

pub(crate) fn read_sdf_records_with_options(
    input: &str,
    options: SdfParseOptions,
) -> Result<Vec<SdfRecordInterpretation>, Box<dyn std::error::Error>> {
    let document = sdf::parse_str(input, options)?;
    Ok(document
        .records()
        .iter()
        .map(|record| record.interpret(Default::default()))
        .collect::<Result<_, _>>()?)
}

pub(crate) fn read_sdf_molecules(input: &str) -> Result<Vec<Molecule>, Box<dyn std::error::Error>> {
    read_sdf_records(input)?
        .into_iter()
        .map(|record| exactly_one(record.into_molecules(), "SDF record"))
        .collect::<Result<Vec<_>, _>>()
}

pub(crate) fn read_sdf_molecules_with_options(
    input: &str,
    options: SdfParseOptions,
) -> Result<Vec<Molecule>, Box<dyn std::error::Error>> {
    read_sdf_records_with_options(input, options)?
        .into_iter()
        .map(|record| exactly_one(record.into_molecules(), "SDF record"))
        .collect::<Result<Vec<_>, _>>()
}

fn exactly_one(
    mut molecules: Vec<Molecule>,
    source: &str,
) -> Result<Molecule, Box<dyn std::error::Error>> {
    if molecules.len() != 1 {
        return Err(format!(
            "expected one molecule from {source}, found {}",
            molecules.len()
        )
        .into());
    }
    Ok(molecules.pop().expect("length checked"))
}

pub(crate) trait MolfileInterpretationTestExt {
    fn molecule(&self) -> &Molecule;
    fn report(&self) -> &molfile::MolfileInterpretationReport;
    fn into_molecule(self) -> Molecule;
}

impl MolfileInterpretationTestExt for molfile::MolfileInterpretation {
    fn molecule(&self) -> &Molecule {
        assert_eq!(
            self.reports().len(),
            1,
            "test fixture must have one component"
        );
        self.molecules().next().expect("one component")
    }

    fn report(&self) -> &molfile::MolfileInterpretationReport {
        assert_eq!(
            self.reports().len(),
            1,
            "test fixture must have one component"
        );
        &self.reports()[0]
    }

    fn into_molecule(self) -> Molecule {
        let mut molecules = self.into_molecules();
        assert_eq!(molecules.len(), 1, "test fixture must have one component");
        molecules.pop().expect("length checked")
    }
}

pub(crate) trait SdfRecordTestExt {
    fn molecule(&self) -> &Molecule;
}

impl SdfRecordTestExt for SdfRecordInterpretation {
    fn molecule(&self) -> &Molecule {
        let mut molecules = self.molecules();
        let molecule = molecules.next().expect("test record must have a component");
        assert!(
            molecules.next().is_none(),
            "test record must have one component"
        );
        molecule
    }
}

pub(crate) fn test_structure(molecule: &Molecule) -> Structure {
    let positions = test_positions(vec![Point3::default(); molecule.atom_count()]);
    Structure::from_molecule(molecule.clone(), &positions).expect("test structure builds")
}

pub(crate) fn element_atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).expect("test element should be available"))
}

pub(crate) fn aromatic_carbon_no_hydrogens() -> Atom {
    let mut atom = carbon();
    atom.hydrogens = ImplicitHydrogens::Fixed(0);
    atom
}

pub(crate) fn charged_atom(symbol: &str, formal_charge: i8) -> Atom {
    let mut atom = element_atom(symbol);
    atom.formal_charge = formal_charge;
    atom
}

pub(crate) fn coordinate_axis_graph(three_dimensional: bool) -> (Molecule, Positions, BondId) {
    // Explicit sp2 endpoints: Br-C(=O)-C(=S)-Cl. This tests the opted-in
    // coordinate axis geometry without inventing aromatic perception flags.
    let mut mol = crate::molecule::MoleculeEditor::new();
    let left = mol
        .add_atom(aromatic_carbon_no_hydrogens())
        .expect("atom identifier capacity");
    let right = mol
        .add_atom(aromatic_carbon_no_hydrogens())
        .expect("atom identifier capacity");
    let left_reference = mol
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    let left_other = mol
        .add_atom(element_atom("O"))
        .expect("atom identifier capacity");
    let right_reference = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let right_other = mol
        .add_atom(element_atom("S"))
        .expect("atom identifier capacity");
    let axis = mol.add_bond(left, right, BondOrder::Single).expect("axis");
    mol.add_bond(left, left_reference, BondOrder::Single)
        .expect("left reference");
    mol.add_bond(left, left_other, BondOrder::Double)
        .expect("left other");
    mol.add_bond(right, right_reference, BondOrder::Single)
        .expect("right reference");
    mol.add_bond(right, right_other, BondOrder::Double)
        .expect("right other");
    let right_reference_point = if three_dimensional {
        Point3::new(1.0, 0.0, 1.0)
    } else {
        Point3::new(1.0, 1.0, 0.0)
    };
    let right_other_point = if three_dimensional {
        Point3::new(1.0, 0.0, -1.0)
    } else {
        Point3::new(1.0, -1.0, 0.0)
    };
    let positions = Positions::new(crate::units::Quantity::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            right_reference_point,
            right_other_point,
        ],
        crate::units::ANGSTROM,
    ))
    .unwrap();
    let mol = mol.finish().expect("connected axis molecule");
    (mol, positions, axis)
}

pub(crate) fn test_positions(points: Vec<Point3>) -> Positions {
    Positions::new(crate::units::Quantity::new(points, crate::units::ANGSTROM))
        .expect("finite test positions")
}

pub(crate) fn ring_molecule(
    symbols: &[&str],
    orders: &[BondOrder],
) -> (Molecule, Vec<AtomId>, Vec<BondId>) {
    assert_eq!(symbols.len(), orders.len());
    let mut mol = crate::molecule::MoleculeEditor::new();
    let atoms = symbols
        .iter()
        .map(|symbol| {
            mol.add_atom(Atom::new(
                Element::from_symbol(symbol).expect("test element should be available"),
            ))
            .expect("atom identifier capacity")
        })
        .collect::<Vec<_>>();
    let mut bonds = Vec::new();
    for index in 0..atoms.len() {
        let next = (index + 1) % atoms.len();
        bonds.push(
            mol.add_bond(atoms[index], atoms[next], orders[index])
                .expect("ring bond should be valid"),
        );
    }
    (mol.finish().expect("connected ring molecule"), atoms, bonds)
}

pub(crate) fn sorted_atom_ids(ids: impl IntoIterator<Item = AtomId>) -> Vec<AtomId> {
    let mut ids = ids.into_iter().collect::<Vec<_>>();
    ids.sort();
    ids
}

pub(crate) fn sorted_bond_ids(ids: impl IntoIterator<Item = BondId>) -> Vec<BondId> {
    let mut ids = ids.into_iter().collect::<Vec<_>>();
    ids.sort();
    ids
}

#[derive(Debug, Clone, PartialEq)]
struct RepresentedAtomSnapshot {
    element: Element,
    isotope: Option<u16>,
    formal_charge: i8,
    radical: Option<AtomRadical>,
    hydrogens: ImplicitHydrogens,
    atom_map: Option<u32>,
}

#[derive(Debug, Clone, PartialEq)]
struct RepresentedBondSnapshot {
    a: AtomId,
    b: AtomId,
    order: BondOrder,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RepresentedStereoElementSnapshot {
    kind: StereoElementKind,
    group: Option<StereoGroupId>,
}

/// Complete primary molecule state, excluding only `Perception`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RepresentedMoleculeSnapshot {
    atoms: Vec<Option<RepresentedAtomSnapshot>>,
    bonds: Vec<Option<RepresentedBondSnapshot>>,
    adjacency: Vec<Vec<BondId>>,
    stereo_elements: Vec<Option<RepresentedStereoElementSnapshot>>,
    stereo_groups: Vec<Option<StereoGroup>>,
}

pub(crate) fn represented_molecule_snapshot(molecule: &Molecule) -> RepresentedMoleculeSnapshot {
    RepresentedMoleculeSnapshot {
        atoms: molecule
            .atoms()
            .map(|(_, atom)| {
                Some(RepresentedAtomSnapshot {
                    element: atom.element,
                    isotope: atom.isotope.map(std::num::NonZero::get),
                    formal_charge: atom.formal_charge,
                    radical: atom.radical,
                    hydrogens: atom.hydrogens,
                    atom_map: atom.atom_map.map(std::num::NonZero::get),
                })
            })
            .collect(),
        bonds: molecule
            .bonds()
            .map(|(_, bond)| {
                Some(RepresentedBondSnapshot {
                    a: bond.a,
                    b: bond.b,
                    order: bond.order,
                })
            })
            .collect(),
        adjacency: molecule
            .atom_ids()
            .map(|atom| {
                molecule
                    .incident_bonds(atom)
                    .expect("live atom")
                    .map(|(id, _)| id)
                    .collect()
            })
            .collect(),
        stereo_elements: molecule
            .stereo_elements()
            .map(|(_, element)| {
                Some(RepresentedStereoElementSnapshot {
                    kind: element.kind.clone(),
                    group: element.group,
                })
            })
            .collect(),
        stereo_groups: molecule
            .stereo_groups()
            .map(|(_, group)| Some(group.clone()))
            .collect(),
    }
}

pub(crate) fn deterministic_text_mutations(seed: &str) -> Vec<String> {
    let mut mutations = vec![String::new(), seed.to_owned()];
    for index in 0..=seed.len().min(128) {
        for inserted in ["\0", "\n", "%", "[", "]", "é"] {
            let mut value = seed.to_owned();
            value.insert_str(index, inserted);
            mutations.push(value);
        }
        if index < seed.len() {
            let mut removed = seed.to_owned();
            removed.remove(index);
            mutations.push(removed);

            let mut replaced = seed.to_owned();
            replaced.replace_range(index..index + 1, "\u{7f}");
            mutations.push(replaced);
        }
    }
    mutations
}

pub(crate) fn mark_all_fresh(mol: &mut Molecule) {
    let _ = valence_api::perceive_valence(mol, ValenceModel::RdkitLike, Default::default());
    rings_api::perceive_rings(mol, RingPerceptionOptions::membership_only())
        .expect("membership perception cannot fail");
    let _ = mol
        .perception()
        .rings()
        .map(|rings| rings.membership())
        .expect("rings installed")
        .clone();
    mol.begin_aromaticity(AromaticityModel::RdkitLike);
}

pub(crate) fn assert_all_stale(mol: &Molecule) {
    assert_perception_stale(mol.perception());
}

pub(crate) fn assert_draft_stale(editor: &MoleculeEditor) {
    assert_perception_stale(editor.perception());
}

fn assert_perception_stale(perception: &Perception) {
    assert!(perception.valence().is_none());
    assert!(perception.rings().is_none());
    assert!(perception.aromaticity().is_none());
    assert!(perception.cip().is_none());
}

pub(crate) fn tetrahedral_marked_graph() -> (MoleculeEditor, AtomId, Vec<AtomId>, BondId) {
    let mut mol = crate::molecule::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let carriers = ["F", "Cl", "Br", "I"]
        .into_iter()
        .map(element_atom)
        .map(|atom| mol.add_atom(atom).expect("atom identifier capacity"))
        .collect::<Vec<_>>();
    let mut bonds = Vec::new();
    for carrier in &carriers {
        bonds.push(
            mol.add_bond(center, *carrier, BondOrder::Single)
                .expect("tetrahedral carrier bond"),
        );
    }
    (mol, center, carriers, bonds[0])
}

/// Bypasses checked insertion so validation can observe malformed stereo.
pub(crate) fn insert_unchecked_stereo(
    editor: &mut MoleculeEditor,
    element: StereoElement,
) -> StereoElementId {
    editor.graph_mut().push_stereo_element(element)
}

pub(crate) fn implicit_h_wedge_geometry_molblock() -> &'static str {
    "\
implicit H geometry wedge
kekule

  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  4  0  0  0  0  0  0
    1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000   -1.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  6  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
"
}

pub(crate) fn rdkit_rp6306_atrop_molblock() -> &'static str {
    include_str!("../tests/fixtures/corpus/smoke/data/rdkit_atropisomers/RP-6306_atrop1.mol")
}

pub(crate) fn rdkit_rp6306_atrop3_molblock() -> &'static str {
    include_str!("../tests/fixtures/corpus/smoke/data/rdkit_atropisomers/RP-6306_atrop3.mol")
}

pub(crate) fn rdkit_rp6306_atrop4_molblock() -> &'static str {
    include_str!("../tests/fixtures/corpus/smoke/data/rdkit_atropisomers/RP-6306_atrop4.mol")
}

pub(crate) fn rdkit_bms986142_atrop4_molblock() -> &'static str {
    include_str!("../tests/fixtures/corpus/smoke/data/rdkit_atropisomers/BMS-986142_atrop4.mol")
}

pub(crate) fn rdkit_bms986142_atrop5_molblock() -> &'static str {
    include_str!("../tests/fixtures/corpus/smoke/data/rdkit_atropisomers/BMS-986142_atrop5.mol")
}

pub(crate) fn rdkit_jdq443_atrop1_molblock() -> &'static str {
    include_str!("../tests/fixtures/corpus/smoke/data/rdkit_atropisomers/JDQ443_atrop1.mol")
}

pub(crate) fn rdkit_zm374979_atrop1_molblock() -> &'static str {
    include_str!("../tests/fixtures/corpus/smoke/data/rdkit_atropisomers/ZM374979_atrop1.mol")
}

pub(crate) fn rdkit_zm374979_atrop2_molblock() -> &'static str {
    include_str!("../tests/fixtures/corpus/smoke/data/rdkit_atropisomers/ZM374979_atrop2.mol")
}

pub(crate) fn rdkit_macrocycle8_ortho_wedge_molblock() -> &'static str {
    include_str!(
        "../tests/fixtures/corpus/smoke/data/rdkit_atropisomers/macrocycle-8-ortho-wedge.mol"
    )
}

pub(crate) fn rdkit_macrocycle8_ortho_hash_molblock() -> &'static str {
    include_str!(
        "../tests/fixtures/corpus/smoke/data/rdkit_atropisomers/macrocycle-8-ortho-hash.mol"
    )
}

// Test-only re-encoding of focused source fixtures. No chemistry is generated.
pub(crate) fn token(value: crate::mmcif::MmcifValue<'_>) -> String {
    let text = value.text();
    if value.is_missing() {
        text.to_owned()
    } else if text.is_empty()
        || text.chars().any(char::is_whitespace)
        || text.starts_with(['_', '#', ';', '\'', '"', '$', '[', ']'])
        || matches!(text, "." | "?")
    {
        assert!(
            !text.contains('\'') && !text.contains('\n'),
            "focused fixture token"
        );
        format!("'{text}'")
    } else {
        text.to_owned()
    }
}

pub(crate) fn singleton_loops_as_scalars(input: &str) -> String {
    let doc = crate::mmcif::parse_str(input, Default::default()).unwrap();
    let mut output = String::new();
    for block in doc.blocks() {
        output.push_str(&format!("data_{}\n", block.name()));
        for entry in block.entries() {
            match entry {
                crate::mmcif::MmcifEntry::Item(item) => {
                    output.push_str(&format!("{} {}\n", item.tag(), token(item.value())))
                }
                crate::mmcif::MmcifEntry::Loop(table) if table.row_count() == 1 => {
                    for tag in table.tags() {
                        output
                            .push_str(&format!("{tag} {}\n", token(table.value(0, tag).unwrap())));
                    }
                }
                crate::mmcif::MmcifEntry::Loop(table) => {
                    output.push_str(&format!("loop_\n{}\n", table.tags().join("\n")));
                    for row in 0..table.row_count() {
                        for value in table.row(row).unwrap() {
                            output.push_str(&format!("{}\n", token(value)));
                        }
                    }
                }
            }
        }
    }
    output
}
