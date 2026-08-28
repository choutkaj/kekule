use crate::chemistry::*;
use crate::core::*;
use crate::geometry::Point3;
use crate::perception::{
    aromaticity as aromaticity_api, aromaticity::*, rings as rings_api, rings::*,
    valence as valence_api, valence::*,
};
use crate::sdf::*;
use crate::smiles::*;
use crate::structure::{Model, Positions};
use crate::{
    canon, molfile, perception as perception_api, sdf, smiles as smiles_api, stereo as stereo_api,
    stereo::*,
};

pub(super) fn carbon() -> Atom {
    Atom::new(Element::from_symbol("C").expect("carbon should be available"))
}

pub(super) fn oxygen() -> Atom {
    Atom::new(Element::from_symbol("O").expect("oxygen should be available"))
}

pub(super) fn read_smiles(
    input: &str,
) -> std::result::Result<Molecule, Box<dyn std::error::Error>> {
    let document = smiles_api::parse_str(input)?;
    Ok(smiles_api::interpret(&document)?.to_molecule()?)
}

pub(super) fn read_smiles_components(
    input: &str,
) -> std::result::Result<Vec<Molecule>, Box<dyn std::error::Error>> {
    let document = smiles_api::parse_str(input)?;
    Ok(smiles_api::interpret(&document)?.to_molecules())
}

pub(super) fn read_smiles_component(
    input: &str,
    component: usize,
) -> std::result::Result<Molecule, Box<dyn std::error::Error>> {
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

pub(super) fn perceive(
    molecule: &mut Molecule,
) -> std::result::Result<(), Box<dyn std::error::Error>> {
    molecule.perceive()?;
    Ok(())
}

pub(super) fn canonical_smiles_round_trip(molecule: &Molecule) -> (String, Molecule) {
    let written = smiles_api::write_canonical(molecule)
        .unwrap_or_else(|error| panic!("canonical SMILES should write: {error}"));
    let mut reparsed = read_smiles(&written)
        .unwrap_or_else(|error| panic!("canonical output should parse: {written}: {error}"));
    perceive(&mut reparsed)
        .unwrap_or_else(|error| panic!("canonical output should perceive: {written}: {error}"));
    (written, reparsed)
}

pub(super) fn read_smiles_with_report(
    input: &str,
) -> std::result::Result<(Molecule, SmilesInterpretationReport), Box<dyn std::error::Error>> {
    let document = smiles_api::parse_str(input)?;
    Ok(smiles_api::interpret(&document)?.to_parts()?)
}

pub(super) trait CanonicalizeFixture {
    fn canonicalize_fixture(
        &mut self,
    ) -> std::result::Result<NormalizationReport, NormalizationError>;

    fn canonicalize_fixture_with_source_stereo(
        &mut self,
        source_stereo: &[SourceStereoBondMark],
    ) -> std::result::Result<NormalizationReport, NormalizationError>;
}

impl CanonicalizeFixture for Molecule {
    fn canonicalize_fixture(
        &mut self,
    ) -> std::result::Result<NormalizationReport, NormalizationError> {
        canonicalize_molecule_for_publication(self, None, &[])
    }

    fn canonicalize_fixture_with_source_stereo(
        &mut self,
        source_stereo: &[SourceStereoBondMark],
    ) -> std::result::Result<NormalizationReport, NormalizationError> {
        canonicalize_molecule_for_publication(self, None, source_stereo)
    }
}

pub(super) fn read_molfile(
    input: &str,
) -> std::result::Result<Molecule, Box<dyn std::error::Error>> {
    let document = molfile::parse_str(input)?;
    let molecules = molfile::interpret(&document)?.to_molecules();
    exactly_one(molecules, "molfile")
}

pub(super) fn read_molfile_with_report(
    input: &str,
) -> std::result::Result<(Molecule, molfile::MolfileInterpretationReport), Box<dyn std::error::Error>>
{
    let document = molfile::parse_str(input)?;
    let mut components = molfile::interpret(&document)?.to_components();
    if components.len() != 1 {
        return Err(format!("expected one molfile component, found {}", components.len()).into());
    }
    let (molecule, _positions, report) = components.pop().expect("length checked").to_parts();
    Ok((molecule, report))
}

pub(super) fn read_sdf_records(
    input: &str,
) -> std::result::Result<Vec<SdfRecordInterpretation>, Box<dyn std::error::Error>> {
    read_sdf_records_with_options(input, SdfParseOptions::default())
}

pub(super) fn read_sdf_records_with_options(
    input: &str,
    options: SdfParseOptions,
) -> std::result::Result<Vec<SdfRecordInterpretation>, Box<dyn std::error::Error>> {
    let document = sdf::parse_str_with_options(input, options)?;
    Ok(sdf::interpret(&document)?.to_records())
}

pub(super) fn read_sdf_molecules(
    input: &str,
) -> std::result::Result<Vec<Molecule>, Box<dyn std::error::Error>> {
    read_sdf_records(input)?
        .into_iter()
        .map(|record| exactly_one(record.to_molecules(), "SDF record"))
        .collect::<std::result::Result<Vec<_>, _>>()
}

pub(super) fn read_sdf_molecules_with_options(
    input: &str,
    options: SdfParseOptions,
) -> std::result::Result<Vec<Molecule>, Box<dyn std::error::Error>> {
    read_sdf_records_with_options(input, options)?
        .into_iter()
        .map(|record| exactly_one(record.to_molecules(), "SDF record"))
        .collect::<std::result::Result<Vec<_>, _>>()
}

fn exactly_one(
    mut molecules: Vec<Molecule>,
    source: &str,
) -> std::result::Result<Molecule, Box<dyn std::error::Error>> {
    if molecules.len() != 1 {
        return Err(format!(
            "expected one molecule from {source}, found {}",
            molecules.len()
        )
        .into());
    }
    Ok(molecules.pop().expect("length checked"))
}

pub(super) trait MolfileInterpretationTestExt {
    fn molecule(&self) -> &Molecule;
    fn report(&self) -> &molfile::MolfileInterpretationReport;
    fn to_molecule(self) -> Molecule;
}

impl MolfileInterpretationTestExt for molfile::MolfileInterpretation {
    fn molecule(&self) -> &Molecule {
        assert_eq!(
            self.components().len(),
            1,
            "test fixture must have one component"
        );
        self.components()[0].molecule()
    }

    fn report(&self) -> &molfile::MolfileInterpretationReport {
        assert_eq!(
            self.components().len(),
            1,
            "test fixture must have one component"
        );
        self.components()[0].report()
    }

    fn to_molecule(self) -> Molecule {
        let mut components = self.to_components();
        assert_eq!(components.len(), 1, "test fixture must have one component");
        components.pop().expect("length checked").to_molecule()
    }
}

pub(super) trait SdfRecordTestExt {
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

pub(super) fn test_model(molecule: &Molecule) -> Model {
    let positions = test_positions(vec![Point3::default(); molecule.atom_count()]);
    Model::from_molecule(molecule, &positions).expect("test model builds")
}

pub(super) fn element_atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).expect("test element should be available"))
}

pub(super) fn aromatic_carbon_no_hydrogens() -> Atom {
    let mut atom = carbon();
    atom.hydrogens = HydrogenDeclaration::Fixed(0);
    atom
}

pub(super) fn charged_atom(symbol: &str, formal_charge: i8) -> Atom {
    let mut atom = element_atom(symbol);
    atom.formal_charge = formal_charge;
    atom
}

pub(super) fn coordinate_axis_graph(three_dimensional: bool) -> (Molecule, Positions, BondId) {
    let mut mol = crate::core::MoleculeEditor::new();
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
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let right_reference = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let right_other = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let axis = mol.add_bond(left, right, BondOrder::Single).expect("axis");
    mol.add_bond(left, left_reference, BondOrder::Single)
        .expect("left reference");
    mol.add_bond(left, left_other, BondOrder::Single)
        .expect("left other");
    mol.add_bond(right, right_reference, BondOrder::Single)
        .expect("right reference");
    mol.add_bond(right, right_other, BondOrder::Single)
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
    let mut mol = mol.finish().expect("connected axis molecule");
    mol.begin_aromaticity(AromaticityModel::RdkitLike);
    mol.set_atom_aromatic(left, true);
    mol.set_atom_aromatic(right, true);
    (mol, positions, axis)
}

pub(super) fn test_positions(points: Vec<Point3>) -> Positions {
    Positions::new(crate::units::Quantity::new(points, crate::units::ANGSTROM))
        .expect("finite test positions")
}

pub(super) fn ring_molecule(
    symbols: &[&str],
    orders: &[BondOrder],
) -> (Molecule, Vec<AtomId>, Vec<BondId>) {
    assert_eq!(symbols.len(), orders.len());
    let mut mol = crate::core::MoleculeEditor::new();
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

pub(super) fn sorted_atom_ids(ids: impl IntoIterator<Item = AtomId>) -> Vec<AtomId> {
    let mut ids = ids.into_iter().collect::<Vec<_>>();
    ids.sort();
    ids
}

pub(super) fn sorted_bond_ids(ids: impl IntoIterator<Item = BondId>) -> Vec<BondId> {
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
    hydrogens: HydrogenDeclaration,
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
pub(super) struct RepresentedMoleculeSnapshot {
    atoms: Vec<Option<RepresentedAtomSnapshot>>,
    bonds: Vec<Option<RepresentedBondSnapshot>>,
    adjacency: Vec<Vec<BondId>>,
    stereo_elements: Vec<Option<RepresentedStereoElementSnapshot>>,
    stereo_groups: Vec<Option<StereoGroup>>,
}

pub(super) fn represented_molecule_snapshot(molecule: &Molecule) -> RepresentedMoleculeSnapshot {
    RepresentedMoleculeSnapshot {
        atoms: molecule
            .graph
            .atoms
            .iter()
            .map(|atom| {
                atom.as_ref().map(|atom| RepresentedAtomSnapshot {
                    element: atom.element,
                    isotope: atom.isotope,
                    formal_charge: atom.formal_charge,
                    radical: atom.radical,
                    hydrogens: atom.hydrogens,
                    atom_map: atom.atom_map,
                })
            })
            .collect(),
        bonds: molecule
            .graph
            .bonds
            .iter()
            .map(|bond| {
                bond.as_ref().map(|bond| RepresentedBondSnapshot {
                    a: bond.a,
                    b: bond.b,
                    order: bond.order,
                })
            })
            .collect(),
        adjacency: molecule.graph.adjacency.clone(),
        stereo_elements: molecule
            .graph
            .stereo_elements
            .iter()
            .map(|element| {
                element
                    .as_ref()
                    .map(|element| RepresentedStereoElementSnapshot {
                        kind: element.kind.clone(),
                        group: element.group,
                    })
            })
            .collect(),
        stereo_groups: molecule.graph.stereo_groups.clone(),
    }
}

pub(super) fn deterministic_text_mutations(seed: &str) -> Vec<String> {
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

pub(super) fn mark_all_fresh(mol: &mut Molecule) {
    let _ = valence_api::perceive_valence(mol, ValenceModel::RdkitLike);
    let _ = rings_api::perceive_ring_membership(mol);
    mol.begin_aromaticity(AromaticityModel::RdkitLike);
}

pub(super) fn assert_all_stale(mol: &Molecule) {
    assert!(!mol.perception().has_valence());
    assert!(!mol.perception().has_rings());
    assert!(!mol.perception().has_aromaticity());
    assert!(!mol.perception().has_stereo());
}

pub(super) fn implicit_h_wedge_geometry_molblock() -> &'static str {
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

pub(super) fn rdkit_rp6306_atrop_molblock() -> &'static str {
    include_str!("../../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/RP-6306_atrop1.mol")
}

pub(super) fn rdkit_rp6306_atrop3_molblock() -> &'static str {
    include_str!("../../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/RP-6306_atrop3.mol")
}

pub(super) fn rdkit_rp6306_atrop4_molblock() -> &'static str {
    include_str!("../../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/RP-6306_atrop4.mol")
}

pub(super) fn rdkit_bms986142_atrop4_molblock() -> &'static str {
    include_str!(
        "../../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/BMS-986142_atrop4.mol"
    )
}

pub(super) fn rdkit_bms986142_atrop5_molblock() -> &'static str {
    include_str!(
        "../../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/BMS-986142_atrop5.mol"
    )
}

pub(super) fn rdkit_jdq443_atrop1_molblock() -> &'static str {
    include_str!("../../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/JDQ443_atrop1.mol")
}

pub(super) fn rdkit_zm374979_atrop1_molblock() -> &'static str {
    include_str!("../../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/ZM374979_atrop1.mol")
}

pub(super) fn rdkit_zm374979_atrop2_molblock() -> &'static str {
    include_str!("../../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/ZM374979_atrop2.mol")
}

pub(super) fn rdkit_macrocycle8_ortho_wedge_molblock() -> &'static str {
    include_str!(
        "../../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/macrocycle-8-ortho-wedge.mol"
    )
}

pub(super) fn rdkit_macrocycle8_ortho_hash_molblock() -> &'static str {
    include_str!(
        "../../../../benchmarks/corpora/smoke/data/rdkit_atropisomers/macrocycle-8-ortho-hash.mol"
    )
}

mod canonical;
mod chemistry;
mod cip;
mod core_payload;
mod graph;
mod hierarchy;
mod hydrogens;
mod mmcif_contents;
mod normalization;
mod perception;
mod public_api;
mod query;
mod ring_limits;
mod rotatable_bonds;
mod smiles;
mod v2000;
mod v3000;
mod valence;
