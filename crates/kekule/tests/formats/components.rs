use kekule::components::{Component, ComponentBondOrder, ComponentLibrary, ComponentLinking};
use kekule::diagnostics::{DiagnosticKind, SourceLocation};
use kekule::mmcif;

/// Template facts of one component: type, linking, atom and bond counts,
/// leaving atoms, atoms whose charge is not zero (`None` is unknown), and the
/// number of metal-coordination bonds.
#[derive(Debug, PartialEq)]
struct Facts {
    kind: String,
    linking: ComponentLinking,
    atoms: usize,
    bonds: usize,
    leaving: Vec<String>,
    charges: Vec<(String, Option<i8>)>,
    metal_bonds: usize,
}

impl Facts {
    fn of(component: &Component) -> Self {
        Self {
            kind: component.kind().to_owned(),
            linking: component.linking(),
            atoms: component.atoms().len(),
            bonds: component.bonds().len(),
            leaving: component
                .atoms()
                .iter()
                .filter(|atom| atom.leaving())
                .map(|atom| atom.name().to_owned())
                .collect(),
            charges: component
                .atoms()
                .iter()
                .filter(|atom| atom.formal_charge() != Some(0))
                .map(|atom| (atom.name().to_owned(), atom.formal_charge()))
                .collect(),
            metal_bonds: component
                .bonds()
                .iter()
                .filter(|bond| bond.metal_coordination())
                .count(),
        }
    }

    fn new(
        kind: &str,
        linking: ComponentLinking,
        (atoms, bonds): (usize, usize),
        leaving: &[&str],
        charges: &[(&str, Option<i8>)],
        metal_bonds: usize,
    ) -> Self {
        Self {
            kind: kind.to_owned(),
            linking,
            atoms,
            bonds,
            leaving: leaving.iter().map(|name| (*name).to_owned()).collect(),
            charges: charges
                .iter()
                .map(|(name, charge)| ((*name).to_owned(), *charge))
                .collect(),
            metal_bonds,
        }
    }
}

#[test]
fn bundled_templates_match_their_dictionary_entries() {
    use ComponentLinking::{Nucleotide, Peptide};
    let none = ComponentLinking::None;
    let peptide = "L-PEPTIDE LINKING";
    let library = ComponentLibrary::bundled();
    assert_eq!(library.len(), 82);

    // Counted directly from the CCD entry files by an independent CIF reader;
    // tools/components records the files and their checksums.
    let leaving = ["OXT", "H2", "HXT"];
    let cases = [
        (
            "ALA",
            Facts::new(peptide, Peptide, (13, 12), &leaving, &[], 0),
        ),
        (
            "HIS",
            Facts::new(peptide, Peptide, (21, 21), &leaving, &[("ND1", Some(1))], 0),
        ),
        (
            "LYS",
            Facts::new(peptide, Peptide, (25, 24), &leaving, &[("NZ", Some(1))], 0),
        ),
        (
            "A",
            Facts::new(
                "RNA LINKING",
                Nucleotide,
                (37, 39),
                &["OP3", "HO3'"],
                &[],
                0,
            ),
        ),
        (
            "DA",
            Facts::new(
                "DNA LINKING",
                Nucleotide,
                (36, 38),
                &["OP3", "HO3'"],
                &[],
                0,
            ),
        ),
        ("HOH", Facts::new("NON-POLYMER", none, (3, 2), &[], &[], 0)),
        ("DOD", Facts::new("NON-POLYMER", none, (3, 2), &[], &[], 0)),
        (
            "ZN",
            Facts::new("NON-POLYMER", none, (1, 0), &[], &[("ZN", Some(2))], 0),
        ),
        (
            "HEM",
            Facts::new(
                "NON-POLYMER",
                none,
                (75, 82),
                &[],
                &[
                    ("O2A", Some(-1)),
                    ("O2D", Some(-1)),
                    ("NA", Some(-1)),
                    ("NC", Some(-1)),
                    ("FE", None),
                ],
                4,
            ),
        ),
        (
            "SF4",
            Facts::new(
                "NON-POLYMER",
                none,
                (8, 12),
                &[],
                &[
                    ("FE1", None),
                    ("FE2", None),
                    ("FE3", None),
                    ("FE4", None),
                    ("S1", Some(-2)),
                    ("S2", Some(-2)),
                    ("S3", Some(-2)),
                    ("S4", Some(-2)),
                ],
                12,
            ),
        ),
        (
            "CLA",
            Facts::new(
                "NON-POLYMER",
                none,
                (137, 145),
                &[],
                &[("MG", Some(2)), ("NB", Some(-1)), ("ND", Some(-1))],
                4,
            ),
        ),
        (
            "SO4",
            Facts::new(
                "NON-POLYMER",
                none,
                (5, 4),
                &[],
                &[("O3", Some(-1)), ("O4", Some(-1))],
                0,
            ),
        ),
    ];
    for (id, expected) in cases {
        let component = library.get(id).unwrap_or_else(|| panic!("{id} is bundled"));
        assert_eq!(component.id(), id);
        assert_eq!(component.variant(), None);
        assert_eq!(Facts::of(component), expected, "{id}");
    }

    // Deuterium is hydrogen of mass 2, and lookups ignore case.
    let (_, deuterium) = library.get("dod").unwrap().atom("D1").unwrap();
    assert_eq!(
        (deuterium.element().symbol(), deuterium.isotope()),
        ("H", Some(2))
    );
    // Atoms are found by their dictionary name or its older alternate.
    let histidine = library.get("HIS").unwrap();
    assert_eq!(histidine.atom("1HB").unwrap().1.name(), "HB2");
    assert_eq!(
        histidine.atom("1HB").unwrap().0,
        histidine.atom("HB2").unwrap().0
    );
    assert_eq!(
        histidine.atom("HB2").unwrap().1.alternate_name(),
        Some("1HB")
    );
}

#[test]
fn bundled_histidine_tautomers_are_neutral_and_localized() {
    use ComponentBondOrder::{Double, Single};
    let ring = |component: &Component| {
        let index = |name| component.atom(name).unwrap().0;
        [
            ("CG", "ND1"),
            ("ND1", "CE1"),
            ("CE1", "NE2"),
            ("NE2", "CD2"),
            ("CD2", "CG"),
        ]
        .map(|(a, b)| {
            let (a, b) = (index(a), index(b));
            component
                .bonds()
                .iter()
                .find(|bond| bond.atoms() == [a, b] || bond.atoms() == [b, a])
                .unwrap()
                .order()
        })
    };
    let library = ComponentLibrary::bundled();
    // The dictionary form: imidazolium with both ring hydrogens.
    assert_eq!(
        ring(library.get("HIS").unwrap()),
        [Single, Double, Single, Single, Double]
    );

    let variants = library
        .variants("HIS")
        .map(|variant| {
            let hydrogens = ["HD1", "HE2"].map(|name| variant.atom(name).is_some());
            let charged = variant
                .atoms()
                .iter()
                .filter(|atom| atom.formal_charge() != Some(0))
                .count();
            (
                variant.variant().unwrap(),
                variant.atoms().len(),
                variant.bonds().len(),
                hydrogens,
                charged,
                ring(variant),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        variants,
        [
            (
                "NE2-protonated",
                20,
                20,
                [false, true],
                0,
                [Single, Double, Single, Single, Double]
            ),
            (
                "ND1-protonated",
                20,
                20,
                [true, false],
                0,
                [Single, Single, Double, Single, Double]
            ),
        ]
    );
    assert_eq!(library.variants("ALA").len(), 0);
}

const TOY_COMPONENT: &str = "data_XYZ
_chem_comp.id XYZ
_chem_comp.name 'toy ligand'
_chem_comp.type NON-POLYMER
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
XYZ C1 C 0
XYZ O1 O 0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
XYZ C1 O1 DOUB
";

#[test]
fn component_libraries_read_dictionary_blocks_and_replace_by_identifier() {
    let document = mmcif::parse_str(TOY_COMPONENT).unwrap();
    let library = ComponentLibrary::from_mmcif(&document).unwrap();
    let toy = library.get("xyz").unwrap();
    assert_eq!(
        (toy.id(), toy.name(), toy.kind(), toy.linking()),
        ("XYZ", "toy ligand", "NON-POLYMER", ComponentLinking::None)
    );
    assert_eq!(toy.bonds()[0].atoms(), [0, 1]);
    assert_eq!(toy.bonds()[0].order(), ComponentBondOrder::Double);

    // A caller library extends the bundled one. Replacing a component drops
    // the variants derived from its previous definition.
    let mut combined = ComponentLibrary::bundled().clone();
    combined.extend(&library);
    assert_eq!(combined.len(), 83);
    assert_eq!(combined.variants("HIS").len(), 2);
    let histidine = combined.get("HIS").unwrap().clone();
    assert!(combined.insert(histidine).is_some());
    assert_eq!(combined.variants("HIS").len(), 0);
}

#[test]
fn component_definitions_fail_with_structured_kinds_and_lines() {
    for (case, source, kind, line) in [
        (
            "bond to an unknown atom",
            TOY_COMPONENT.replace("XYZ C1 O1 DOUB", "XYZ C1 N9 DOUB"),
            DiagnosticKind::InvalidReference,
            Some(17),
        ),
        (
            "duplicate atom",
            TOY_COMPONENT.replace("XYZ O1 O 0", "XYZ C1 O 0"),
            DiagnosticKind::Inconsistent,
            Some(11),
        ),
        (
            "unknown element",
            TOY_COMPONENT.replace("XYZ O1 O 0", "XYZ O1 Qq 0"),
            DiagnosticKind::InvalidValue,
            Some(11),
        ),
        (
            "unknown bond order",
            TOY_COMPONENT.replace("O1 DOUB", "O1 WEIRD"),
            DiagnosticKind::InvalidValue,
            Some(17),
        ),
        (
            "row of another component",
            TOY_COMPONENT.replace("XYZ O1 O 0", "ABC O1 O 0"),
            DiagnosticKind::Inconsistent,
            Some(11),
        ),
    ] {
        let document = mmcif::parse_str(&source).unwrap();
        let error = ComponentLibrary::from_mmcif(&document).unwrap_err();
        assert_eq!(
            (error.kind(), error.location()),
            (kind, line.map(SourceLocation::at_line)),
            "{case}: {error}"
        );
    }
}
