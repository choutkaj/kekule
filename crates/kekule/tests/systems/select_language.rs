//! The selection language: parsing, canonical printing, errors with spans,
//! host extensions, and a conformance suite whose expectations come from
//! per-atom reference predicates independent of the evaluator.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::sync::Arc;

use kekule::diagnostics::{DiagnosticKind, Format};
use kekule::geometry::Point3;
use kekule::molecule::Element;
use kekule::properties::{PropertyColumn, PropertyKey};
use kekule::select::{
    self, HostRequest, HostRole, HostTerm, HostTermError, HostTermKind, HostValue, ReferencePoints,
    SelectEvaluateError, SelectEvaluateOptions, SelectParseOptions, SelectionExpr,
    SelectionResolver,
};
use kekule::structure::{Model, Positions};
use kekule::topology::{AtomSelection, AtomView, MoleculeClass, ResidueClass, Topology};
use kekule::units::{Quantity, ANGSTROM, CANONICAL_LENGTH_UNIT, ELEMENTARY_CHARGE, NANOMETER};
use kekule::{mmcif, smiles};

use super::support::fixture;

fn parse(text: &str) -> SelectionExpr {
    select::parse(text, SelectParseOptions::default())
        .unwrap_or_else(|error| panic!("{text:?}: {error}"))
}

fn host_options() -> SelectParseOptions {
    SelectParseOptions::default()
        .with_host_keywords(["displayed", "selected"])
        .unwrap()
}

#[test]
fn canonical_printing_normalizes_and_round_trips() {
    let cases = [
        ("elem C", "element C"),
        ("  protein   and  name CA ", "protein and name CA"),
        (
            "(protein and name CA) and chain A",
            "protein and name CA and chain A",
        ),
        ("protein or (water or ion)", "protein or water or ion"),
        ("protein and (water or ion)", "protein and (water or ion)"),
        ("not (protein or water)", "not (protein or water)"),
        ("not (protein)", "not protein"),
        ("same residue as name CA", "byres name CA"),
        ("same chain as resid 5", "bychain resid 5"),
        ("same molecule as ligand", "bymolecule ligand"),
        ("within 5 A of resname LIG", "within 5 of resname LIG"),
        ("within 5.0 Å of resname LIG", "within 5 of resname LIG"),
        ("within 0.5 nm of ligand", "within 0.5 nm of ligand"),
        (
            "within 5 of resname LIG and not water",
            "within 5 of resname LIG and not water",
        ),
        (
            "within 5 of (resname LIG and not water)",
            "within 5 of (resname LIG and not water)",
        ),
        ("resid 10 to 50", "resid 10-50"),
        ("resid 52A-60 7", "resid 52A-60 7"),
        ("charge -2 to -1", "charge -2 to -1"),
        ("charge -2--1", "charge -2 to -1"),
        ("charge = -1", "charge == -1"),
        ("bfactor>50", "bfactor > 50"),
        ("bfactor 20-50.5", "bfactor 20-50.5"),
        ("name \"C1'\" C1'", "name C1' C1'"),
        ("resname \"and\" 'x y'", "resname \"and\" \"x y\""),
        ("smarts [OH]", "smarts \"[OH]\""),
        ("property q < -0.5 e", "property q < -0.5 e"),
        ("property q '12' 12", "property q \"12\" 12"),
        ("property flag", "property flag"),
        ("property tag a-b", "property tag a-b"),
        ("property q 1-2 e", "property q 1 e to 2 e"),
        ("property q -1 to 2 e", "property q -1 e to 2 e"),
        ("resclass amino_acid dna", "resclass amino_acid dna"),
        ("molclass small_molecule", "molclass small_molecule"),
        (
            "molclass water ion and water",
            "molclass water ion and water",
        ),
        ("bonded to element O", "bonded to element O"),
        ("all", "all"),
        ("none or all", "none or all"),
        ("#2-5 or $site", "#2-5 or $site"),
    ];
    for (input, canonical) in cases {
        let expression = parse(input);
        assert_eq!(expression.to_string(), canonical, "{input:?}");
        let reparsed = parse(canonical);
        assert_eq!(reparsed, expression, "{input:?}");
        assert_eq!(reparsed.to_string(), canonical, "{input:?}");
    }
    let hosted = select::parse("displayed and not selected", host_options()).unwrap();
    assert_eq!(hosted.to_string(), "displayed and not selected");
    let terms = hosted
        .host_terms()
        .iter()
        .map(|term| (term.kind(), term.text().to_owned()))
        .collect::<Vec<_>>();
    assert_eq!(
        terms,
        [
            (HostTermKind::Keyword, "displayed".to_owned()),
            (HostTermKind::Keyword, "selected".to_owned())
        ]
    );
}

#[test]
fn parse_errors_are_typed_and_located() {
    let cases: [(&str, DiagnosticKind, std::ops::Range<usize>); 12] = [
        ("protein water", DiagnosticKind::Syntax, 8..13),
        ("resname", DiagnosticKind::Syntax, 7..7),
        ("(protein", DiagnosticKind::Syntax, 0..1),
        ("protein)", DiagnosticKind::Syntax, 7..8),
        ("element Xx", DiagnosticKind::InvalidValue, 8..10),
        ("within -1 of protein", DiagnosticKind::InvalidValue, 7..9),
        ("bfactor > abc", DiagnosticKind::InvalidValue, 10..13),
        ("name < CA", DiagnosticKind::Syntax, 5..6),
        ("smarts \"C(\"", DiagnosticKind::InvalidValue, 7..11),
        ("frobnicate", DiagnosticKind::Syntax, 0..10),
        ("resname \"ALA", DiagnosticKind::Syntax, 8..12),
        ("", DiagnosticKind::Syntax, 0..0),
    ];
    for (input, kind, span) in cases {
        let error = select::parse(input, SelectParseOptions::default()).unwrap_err();
        assert_eq!(error.format(), Format::Selection, "{input:?}");
        assert_eq!(error.kind(), kind, "{input:?}: {error}");
        assert_eq!(error.span(), Some(span.clone()), "{input:?}: {error}");
        let location = error.location().expect("selection errors are located");
        assert_eq!(location.column(), Some(span.start + 1), "{input:?}");
    }
    let deep = format!("{}all", "not ".repeat(300));
    let error = select::parse(&deep, SelectParseOptions::default()).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Limit);
    assert!(SelectParseOptions::default()
        .with_host_keywords(["and"])
        .is_err());
    assert!(SelectParseOptions::default()
        .with_host_keywords(["two words"])
        .is_err());
}

fn model(name: &str) -> Model {
    let text = std::fs::read_to_string(fixture(name)).unwrap();
    let mut model = mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_model();
    model.perceive().unwrap();
    model
}

/// What a reference predicate may read about one atom.
struct Atom<'a> {
    view: AtomView<'a>,
    model: &'a Model,
}

impl Atom<'_> {
    fn element(&self) -> Element {
        self.view.atom().element()
    }

    fn name(&self) -> Option<&str> {
        let site = self.view.atom_site()?;
        let metadata = site.metadata();
        metadata.auth_atom_id().or(metadata.label_atom_id())
    }

    fn resname(&self) -> Option<&str> {
        self.view.residue().map(|residue| residue.name())
    }

    fn resid(&self) -> Option<(i64, Option<String>)> {
        let residue = self.view.residue()?;
        match residue.author_seq_id() {
            Some(text) => text
                .trim()
                .parse()
                .ok()
                .map(|number| (number, residue.insertion_code().map(str::to_owned))),
            None => residue
                .label_seq_id()
                .map(|number| (i64::from(number), None)),
        }
    }

    fn chain(&self) -> Option<&str> {
        let chain = self.view.chain()?;
        Some(chain.author_id().unwrap_or(chain.label_id()))
    }

    fn class(&self) -> MoleculeClass {
        self.view.instance().class()
    }

    fn residue_class(&self) -> Option<ResidueClass> {
        self.view.residue().map(|residue| residue.class())
    }

    fn b_factor(&self) -> Option<f64> {
        self.model
            .b_factor(self.view.index())
            .and_then(|value| value.value_in(kekule::units::SQUARE_ANGSTROM).ok())
    }

    fn position(&self) -> Point3 {
        self.model.positions().values().value()[self.view.index().index()]
    }

    fn is_protein_backbone(&self) -> bool {
        self.class() == MoleculeClass::Protein
            && self.residue_class() == Some(ResidueClass::AminoAcid)
            && matches!(self.name(), Some("N" | "CA" | "C" | "O"))
    }
}

/// A reference predicate: the expected membership of every atom.
type Predicate = Box<dyn Fn(&[Atom<'_>]) -> Vec<bool>>;

/// Atoms within `cutoff` ångströms of any atom satisfying `of`, by brute
/// force over every pair.
fn near(all: &[Atom<'_>], cutoff: f64, of: impl Fn(&Atom<'_>) -> bool) -> Vec<bool> {
    let cutoff = cutoff / 10.0;
    let references = all
        .iter()
        .filter(|atom| of(atom))
        .map(Atom::position)
        .collect::<Vec<_>>();
    all.iter()
        .map(|atom| {
            let here = atom.position();
            references.iter().any(|there| {
                let (x, y, z) = (here.x - there.x, here.y - there.y, here.z - there.z);
                (x * x + y * y + z * z).sqrt() <= cutoff
            })
        })
        .collect()
}

/// The conformance cases: an expression and its reference predicate.
fn cases() -> Vec<(&'static str, &'static str, Predicate)> {
    fn each(
        structure: &'static str,
        expression: &'static str,
        predicate: impl Fn(&Atom<'_>) -> bool + 'static,
    ) -> (&'static str, &'static str, Predicate) {
        (
            structure,
            expression,
            Box::new(move |all| all.iter().map(&predicate).collect()),
        )
    }
    fn whole(
        structure: &'static str,
        expression: &'static str,
        predicate: impl Fn(&[Atom<'_>]) -> Vec<bool> + 'static,
    ) -> (&'static str, &'static str, Predicate) {
        (structure, expression, Box::new(predicate))
    }
    let ligand = |atom: &Atom<'_>| atom.resname() == Some("AP5");
    vec![
        each("1AKE", "all", |_| true),
        each("1AKE", "none", |_| false),
        each("1AKE", "protein", |atom| {
            atom.class() == MoleculeClass::Protein
        }),
        each("1AKE", "water", |atom| atom.class() == MoleculeClass::Water),
        each("1AKE", "ligand", |atom| {
            atom.class() == MoleculeClass::SmallMolecule
        }),
        each("1AKE", "not protein and not water", |atom| {
            !matches!(atom.class(), MoleculeClass::Protein | MoleculeClass::Water)
        }),
        each("1AKE", "chain B", |atom| atom.chain() == Some("B")),
        each("1AKE", "resname AP5", ligand),
        each("1AKE", "resname GLY ALA", |atom| {
            matches!(atom.resname(), Some("GLY" | "ALA"))
        }),
        each("1AKE", "resid 10-50", |atom| {
            atom.resid().is_some_and(|(number, code)| {
                (10..=50).contains(&number) && (number != 50 || code.is_none())
            })
        }),
        each("1AKE", "resid > 200 and chain A", |atom| {
            atom.resid()
                .is_some_and(|(number, code)| number > 200 || (number == 200 && code.is_some()))
                && atom.chain() == Some("A")
        }),
        each("1AKE", "name CA", |atom| atom.name() == Some("CA")),
        each("1AKE", "name CA and chain B", |atom| {
            atom.name() == Some("CA") && atom.chain() == Some("B")
        }),
        each("1AKE", "backbone", |atom| atom.is_protein_backbone()),
        each("1AKE", "sidechain", |atom| {
            atom.class() == MoleculeClass::Protein
                && atom.residue_class() == Some(ResidueClass::AminoAcid)
                && !matches!(atom.name(), Some("N" | "CA" | "C" | "O" | "OXT"))
                && atom.element() != Element::H
        }),
        each("1AKE", "element C N", |atom| {
            matches!(atom.element(), Element::C | Element::N)
        }),
        each("1AKE", "heavy", |atom| atom.element() != Element::H),
        each("1AKE", "bfactor > 50", |atom| {
            atom.b_factor().is_some_and(|value| value > 50.0)
        }),
        whole("1AKE", "within 5 of resname AP5", move |all| {
            near(all, 5.0, ligand)
        }),
        whole(
            "1AKE",
            "within 0.5 nm of resname AP5 and protein",
            move |all| {
                near(all, 5.0, ligand)
                    .into_iter()
                    .zip(all)
                    .map(|(near, atom)| near && atom.class() == MoleculeClass::Protein)
                    .collect()
            },
        ),
        whole(
            "1AKE",
            "byres (within 4 of resname AP5 and protein)",
            move |all| {
                let near = near(all, 4.0, ligand);
                let residues = all
                    .iter()
                    .zip(&near)
                    .filter(|(atom, near)| **near && atom.class() == MoleculeClass::Protein)
                    .filter_map(|(atom, _)| atom.view.residue().map(|residue| residue.id()))
                    .collect::<std::collections::BTreeSet<_>>();
                all.iter()
                    .map(|atom| {
                        atom.view
                            .residue()
                            .is_some_and(|residue| residues.contains(&residue.id()))
                    })
                    .collect()
            },
        ),
        whole("1AKE", "bonded to name CA", |all| {
            all.iter()
                .map(|atom| {
                    atom.view
                        .neighbors()
                        .any(|neighbor| all[neighbor.index().index()].name() == Some("CA"))
                })
                .collect()
        }),
        each("1AKE", "aromatic and resname PHE", |atom| {
            atom.view.is_aromatic() == Some(true) && atom.resname() == Some("PHE")
        }),
        each("1AKE", "resclass amino_acid and not protein", |atom| {
            atom.residue_class() == Some(ResidueClass::AminoAcid)
                && atom.class() != MoleculeClass::Protein
        }),
        each("1AKE", "molclass water and index < 3500", |atom| {
            atom.class() == MoleculeClass::Water && atom.view.index().index() < 3500
        }),
        each("1CRN", "resname CYS and name SG", |atom| {
            atom.resname() == Some("CYS") && atom.name() == Some("SG")
        }),
        each("1CRN", "sidechain and resname THR", |atom| {
            atom.resname() == Some("THR")
                && !matches!(atom.name(), Some("N" | "CA" | "C" | "O" | "OXT"))
        }),
        each("1CRN", "bychain resid 1", |_| true),
        each("1CRN", "polar_hydrogen or hydrogen", |atom| {
            atom.element() == Element::H
        }),
    ]
}

fn ranges(indices: &[usize]) -> String {
    let mut text = String::new();
    let mut iter = indices.iter().copied().peekable();
    while let Some(start) = iter.next() {
        let mut end = start;
        while iter.peek() == Some(&(end + 1)) {
            end += 1;
            iter.next();
        }
        if !text.is_empty() {
            text.push(',');
        }
        if start == end {
            write!(text, "{start}").unwrap();
        } else {
            write!(text, "{start}-{end}").unwrap();
        }
    }
    if text.is_empty() {
        text.push('-');
    }
    text
}

#[test]
fn conformance_suite_matches_reference_predicates_and_the_committed_file() {
    let models = [
        ("1AKE", model("mmcif/1AKE.cif")),
        ("1CRN", model("mmcif/1CRN.cif")),
    ]
    .into_iter()
    .collect::<BTreeMap<_, _>>();
    let mut file = String::from(
        "# Selection-language conformance: structure, expression, and the expected\n\
         # dense atom indices as inclusive ranges (`-` for none). Structures are the\n\
         # perceived first models of tests/fixtures/mmcif/<structure>.cif. Expected\n\
         # sets come from per-atom reference predicates, not from the evaluator.\n",
    );
    for (structure, expression, predicate) in cases() {
        let model = &models[structure];
        let atoms = model
            .topology()
            .atoms()
            .map(|view| Atom { view, model })
            .collect::<Vec<_>>();
        let expected = predicate(&atoms)
            .into_iter()
            .enumerate()
            .filter_map(|(index, member)| member.then_some(index))
            .collect::<Vec<_>>();
        let selected = parse(expression)
            .evaluate_structure(model, SelectEvaluateOptions::default())
            .unwrap_or_else(|error| panic!("{expression}: {error}"));
        let actual = selected
            .indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        assert_eq!(actual, expected, "{structure}: {expression}");
        writeln!(file, "{structure}\t{expression}\t{}", ranges(&expected)).unwrap();
    }
    let path = fixture("select/conformance.txt");
    if std::env::var_os("KEKULE_BLESS").is_some() {
        std::fs::write(&path, &file).unwrap();
    }
    let committed = std::fs::read_to_string(&path)
        .unwrap()
        .replace("\r\n", "\n");
    assert_eq!(
        committed, file,
        "rerun with KEKULE_BLESS=1 after a reviewed change"
    );
}

#[test]
fn topology_evaluation_reports_terms_that_need_more() {
    let model = model("mmcif/1CRN.cif");
    let topology = model.shared_topology();
    let error = parse("protein and within 5 of name SG")
        .evaluate(topology.clone(), SelectEvaluateOptions::default())
        .unwrap_err();
    assert!(matches!(
        error,
        SelectEvaluateError::RequiresStructure { .. }
    ));
    assert_eq!(error.span(), 12..31);
    assert!(matches!(
        parse("bfactor > 10").evaluate(topology.clone(), SelectEvaluateOptions::default()),
        Err(SelectEvaluateError::RequiresStructure { .. })
    ));
    let unperceived = smiles::to_topology("c1ccccc1O").unwrap();
    let error = parse("aromatic")
        .evaluate(unperceived.clone(), SelectEvaluateOptions::default())
        .unwrap_err();
    assert!(matches!(
        error,
        SelectEvaluateError::MissingPerception { .. }
    ));
    assert!(matches!(
        parse("smarts \"cO\"").evaluate(unperceived, SelectEvaluateOptions::default()),
        Err(SelectEvaluateError::Smarts { .. })
    ));
    // Topology-only terms need no structure.
    let names = parse("name SG")
        .evaluate(topology, SelectEvaluateOptions::default())
        .unwrap();
    assert_eq!(names.len(), 6);
}

#[test]
fn chemistry_terms_read_installed_perception() {
    let topology = Arc::new(
        smiles::to_topology("c1ccccc1CC(=O)[OH].[Na+]")
            .unwrap()
            .perceived()
            .unwrap(),
    );
    let count = |text: &str| {
        parse(text)
            .evaluate(topology.clone(), SelectEvaluateOptions::default())
            .unwrap()
            .len()
    };
    assert_eq!(count("aromatic"), 6);
    assert_eq!(count("ring"), 6);
    assert_eq!(count("smarts \"[OX2H]\""), 1);
    assert_eq!(count("smarts \"C(=O)O\""), 3);
    assert_eq!(count("ion"), 1);
    assert_eq!(count("charge > 0"), 1);
    assert_eq!(count("element na"), 1);
    assert_eq!(count("not aromatic and heavy"), 5);
}

#[test]
fn hydrogen_terms_distinguish_polar_hydrogens() {
    let molecule = smiles::to_molecules("[H]OC([H])([H])N([H])[H]")
        .unwrap()
        .remove(0);
    let topology = Arc::new(Topology::from_molecule(molecule).unwrap());
    let count = |text: &str| {
        parse(text)
            .evaluate(topology.clone(), SelectEvaluateOptions::default())
            .unwrap()
            .len()
    };
    assert_eq!(count("hydrogen"), 5);
    assert_eq!(count("polar_hydrogen"), 3);
    assert_eq!(count("heavy"), 3);
    assert_eq!(count("bonded to element O"), 2);
}

#[test]
fn properties_compare_in_their_units_and_scopes() {
    let mut topology = smiles::to_topology("CCO").unwrap();
    let key = |name: &str| PropertyKey::new(name).unwrap();
    {
        let topology = Arc::make_mut(&mut topology);
        let mut properties = topology.properties_mut();
        let mut atoms = properties.atoms_mut();
        atoms
            .insert(
                key("q"),
                PropertyColumn::Real {
                    unit: ELEMENTARY_CHARGE,
                    values: vec![Some(-0.2), Some(0.1), Some(-0.6)],
                },
            )
            .unwrap();
        atoms
            .insert(
                key("radius"),
                PropertyColumn::Real {
                    unit: NANOMETER,
                    values: vec![Some(0.17), None, Some(0.15)],
                },
            )
            .unwrap();
        atoms
            .insert(
                key("tag"),
                PropertyColumn::String(vec![Some("x".into()), None, Some("y".into())]),
            )
            .unwrap();
        atoms
            .insert(
                key("flag"),
                PropertyColumn::Bool(vec![Some(true), Some(false), None]),
            )
            .unwrap();
        atoms
            .insert(
                key("id"),
                PropertyColumn::Int(vec![
                    Some(9_007_199_254_740_992),
                    Some(9_007_199_254_740_993),
                    Some(-3),
                ]),
            )
            .unwrap();
    }
    let indices = |text: &str| {
        parse(text)
            .evaluate(topology.clone(), SelectEvaluateOptions::default())
            .unwrap_or_else(|error| panic!("{text}: {error}"))
            .indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>()
    };
    assert_eq!(indices("property q < 0"), [0, 2]);
    assert_eq!(indices("property q < -0.5 e"), [2]);
    assert_eq!(indices("property radius > 1.6 A"), [0]);
    assert_eq!(indices("property radius 0.15"), [2]);
    assert_eq!(indices("property tag x"), [0]);
    assert_eq!(indices("property flag"), [0]);
    assert_eq!(indices("property flag false"), [1]);
    assert_eq!(indices("property id == 9007199254740993"), [1]);
    assert_eq!(indices("property id 9007199254740992"), [0]);
    assert_eq!(indices("property id < -2.5"), [2]);
    assert_eq!(indices("property id > 9007199254740992.5"), [1]);
    let error = |text: &str| {
        parse(text)
            .evaluate(topology.clone(), SelectEvaluateOptions::default())
            .unwrap_err()
    };
    assert!(matches!(
        error("property missing 1"),
        SelectEvaluateError::UnknownProperty { .. }
    ));
    assert!(matches!(
        error("property q > 1 nm"),
        SelectEvaluateError::PropertyType { .. }
    ));
    assert!(matches!(
        error("property tag 1 nm"),
        SelectEvaluateError::PropertyType { .. }
    ));
}

struct Resolver {
    other: ReferencePoints,
}

impl SelectionResolver for Resolver {
    fn resolve(
        &self,
        term: &HostTerm<'_>,
        request: HostRequest<'_>,
    ) -> Result<HostValue, HostTermError> {
        match (term.kind(), term.text(), request.role()) {
            (HostTermKind::Named, "first", _) => Ok(HostValue::Atoms(
                AtomSelection::from_indices(
                    request.shared_topology(),
                    [kekule::topology::TopologyAtomIndex::new(0)],
                )
                .unwrap(),
            )),
            (HostTermKind::Reference, "2", HostRole::Reference) => {
                Ok(HostValue::Points(self.other.clone()))
            }
            (HostTermKind::Reference, "2", HostRole::Atoms) => Ok(HostValue::Atoms(
                AtomSelection::empty(request.shared_topology()),
            )),
            (HostTermKind::Keyword, "selected", _) => Ok(HostValue::Atoms(AtomSelection::all(
                request.shared_topology(),
            ))),
            (HostTermKind::Reference, "9", _) => Ok(HostValue::Points(self.other.clone())),
            (HostTermKind::Reference, "foreign", _) => Ok(HostValue::Atoms(AtomSelection::all(
                smiles::to_topology("C").unwrap(),
            ))),
            _ => Err(HostTermError::new(format!("unknown {}", term.text()))),
        }
    }
}

#[test]
fn host_terms_resolve_atoms_and_points_across_objects() {
    let topology =
        Arc::new(Topology::from_molecules(smiles::to_molecules("C.C.C").unwrap()).unwrap());
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(20.0, 0.0, 0.0),
    ];
    let model = Model::new(
        topology,
        Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
    )
    .unwrap();
    let resolver = Resolver {
        other: Positions::new(Quantity::new(vec![Point3::new(19.0, 0.0, 0.0)], ANGSTROM))
            .unwrap()
            .into(),
    };
    let options = SelectEvaluateOptions::default().with_resolver(&resolver);
    let indices = |text: &str| {
        select::parse(text, host_options())
            .unwrap()
            .evaluate_structure(&model, options)
            .unwrap_or_else(|error| panic!("{text}: {error}"))
            .indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>()
    };
    assert_eq!(indices("$first"), [0]);
    assert_eq!(indices("within 2 of #2"), [2]);
    assert_eq!(indices("#2"), Vec::<usize>::new());
    assert_eq!(indices("selected and not $first"), [1, 2]);
    assert_eq!(indices("within 11 of $first"), [0, 1]);

    let error = |text: &str, options| {
        select::parse(text, host_options())
            .unwrap()
            .evaluate_structure(&model, options)
            .unwrap_err()
    };
    assert!(matches!(
        error("$first", SelectEvaluateOptions::default()),
        SelectEvaluateError::NoResolver { .. }
    ));
    assert!(matches!(
        error("#9", options),
        SelectEvaluateError::HostPoints { .. }
    ));
    assert!(matches!(
        error("#foreign", options),
        SelectEvaluateError::HostTopologyMismatch { .. }
    ));
    let unknown = error("$nothing", options);
    assert!(matches!(unknown, SelectEvaluateError::Host { .. }));
    assert_eq!(unknown.span(), 0..8);
}

/// Structures kept apart in one frame, as a viewer holds a receptor and its
/// ligand: a `within` reference evaluated on one also measures against the
/// atoms it selects in the others.
struct Scene<'a> {
    models: &'a [Model],
    /// The evaluated model.
    current: usize,
}

impl Scene<'_> {
    fn select(&self, text: &str) -> Result<Vec<usize>, SelectEvaluateError> {
        let options = SelectEvaluateOptions::default().with_resolver(self);
        let selection = parse(text).evaluate_structure(&self.models[self.current], options)?;
        Ok(selection
            .indices()
            .iter()
            .map(|index| index.index())
            .collect())
    }
}

impl SelectionResolver for Scene<'_> {
    fn resolve(&self, term: &HostTerm<'_>, _: HostRequest<'_>) -> Result<HostValue, HostTermError> {
        Err(HostTermError::new(format!("unknown {}", term.text())))
    }

    fn reference_points(
        &self,
        reference: &SelectionExpr,
        request: HostRequest<'_>,
    ) -> Result<ReferencePoints, HostTermError> {
        let evaluated = &self.models[self.current];
        assert_eq!(request.role(), HostRole::Reference);
        assert!(request.topology().shares_layout(evaluated.topology()));
        assert!(request
            .structure()
            .is_some_and(|view| std::ptr::eq(view.positions(), evaluated.positions())));
        let mut points = Vec::new();
        for (index, model) in self.models.iter().enumerate() {
            if index == self.current {
                continue;
            }
            let scene = Scene {
                models: self.models,
                current: index,
            };
            let options = SelectEvaluateOptions::default().with_resolver(&scene);
            let atoms = reference
                .evaluate_structure(model, options)
                .map_err(|error| {
                    HostTermError::new(format!("evaluating on model {index} failed"))
                        .with_source(error)
                })?;
            let positions: &[Point3] = model.positions().values().value();
            points.extend(atoms.indices().iter().map(|atom| positions[atom.index()]));
        }
        Positions::from_vec(Quantity::new(points, CANONICAL_LENGTH_UNIT))
            .map(ReferencePoints::new)
            .map_err(|error| HostTermError::new("invalid points").with_source(error))
    }
}

/// Indices of the points within `cutoff` ångströms of any reference, by brute
/// force over every pair.
fn within_points(points: &[Point3], references: &[Point3], cutoff: f64) -> Vec<usize> {
    let cutoff = cutoff / 10.0;
    points
        .iter()
        .enumerate()
        .filter(|(_, here)| {
            references.iter().any(|there| {
                let (x, y, z) = (here.x - there.x, here.y - there.y, here.z - there.z);
                (x * x + y * y + z * z).sqrt() <= cutoff
            })
        })
        .map(|(index, _)| index)
        .collect()
}

#[test]
fn within_measures_against_structures_the_resolver_holds() {
    // 1AKE split as a docking project keeps it: receptor, ligand, and water
    // as separate models.
    let complex = model("mmcif/1AKE.cif");
    let part = |text: &str| {
        let atoms = parse(text)
            .evaluate_structure(&complex, SelectEvaluateOptions::default())
            .unwrap();
        (complex.subset(&atoms).unwrap(), atoms)
    };
    let (receptor, protein) = part("protein");
    let models = [receptor, part("ligand").0, part("water").0];
    let points = |index: usize| models[index].positions().values().value().to_vec();
    let on_receptor = Scene {
        models: &models,
        current: 0,
    };
    // The subset keeps the protein atoms in dense order.
    let in_complex = |text: &str| {
        parse(text)
            .evaluate_structure(&complex, SelectEvaluateOptions::default())
            .unwrap()
            .indices()
            .iter()
            .map(|index| protein.indices().binary_search(index).unwrap())
            .collect::<Vec<_>>()
    };

    // The pocket: receptor atoms within 5 Å of a ligand atom.
    let pocket = on_receptor.select("within 5 of ligand").unwrap();
    assert!(!pocket.is_empty());
    assert_eq!(pocket, within_points(&points(0), &points(1), 5.0));
    assert_eq!(pocket, in_complex("within 5 of ligand and protein"));

    // Water-mediated contacts: the nested reference reaches the resolver
    // again from the water model, which measures against the ligand model.
    let bridged = on_receptor
        .select("within 3.5 of (water and within 3.5 of ligand)")
        .unwrap();
    let water = points(2);
    let bridging = within_points(&water, &points(1), 3.5)
        .into_iter()
        .map(|index| water[index])
        .collect::<Vec<_>>();
    assert!(!bridging.is_empty() && !bridged.is_empty());
    assert_eq!(bridged, within_points(&points(0), &bridging, 3.5));
    assert_eq!(
        bridged,
        in_complex("within 3.5 of (water and within 3.5 of ligand) and protein")
    );

    // Without reference points, the receptor alone holds no ligand.
    let evaluate = |options| {
        parse("within 5 of ligand")
            .evaluate_structure(&models[0], options)
            .unwrap()
    };
    assert!(evaluate(SelectEvaluateOptions::default()).is_empty());
    let terms_only = Resolver {
        other: ReferencePoints::default(),
    };
    assert!(evaluate(SelectEvaluateOptions::default().with_resolver(&terms_only)).is_empty());
}

/// Structures in one frame answered as a host evaluating one selection over
/// all of them does: each reference's points over every structure, the
/// evaluated one included, built once and shared by every evaluation.
struct SharedScene<'a> {
    models: &'a [Model],
    current: usize,
    /// Reference points by canonical text, and how often they were built.
    built: &'a RefCell<BTreeMap<String, (ReferencePoints, usize)>>,
}

impl SharedScene<'_> {
    fn select(&self, text: &str) -> Vec<usize> {
        let options = SelectEvaluateOptions::default().with_resolver(self);
        parse(text)
            .evaluate_structure(&self.models[self.current], options)
            .unwrap_or_else(|error| panic!("{text}: {error}"))
            .indices()
            .iter()
            .map(|index| index.index())
            .collect()
    }
}

impl SelectionResolver for SharedScene<'_> {
    fn resolve(&self, term: &HostTerm<'_>, _: HostRequest<'_>) -> Result<HostValue, HostTermError> {
        Err(HostTermError::new(format!("unknown {}", term.text())))
    }

    fn reference_points(
        &self,
        reference: &SelectionExpr,
        _: HostRequest<'_>,
    ) -> Result<ReferencePoints, HostTermError> {
        let key = reference.to_string();
        if let Some((points, count)) = self.built.borrow_mut().get_mut(&key) {
            *count += 1;
            return Ok(points.clone());
        }
        let mut points = Vec::new();
        for (index, model) in self.models.iter().enumerate() {
            let scene = SharedScene {
                current: index,
                ..*self
            };
            let options = SelectEvaluateOptions::default().with_resolver(&scene);
            let atoms = reference.evaluate_structure(model, options).unwrap();
            let positions: &[Point3] = model.positions().values().value();
            points.extend(atoms.indices().iter().map(|atom| positions[atom.index()]));
        }
        let points = ReferencePoints::new(
            Positions::from_vec(Quantity::new(points, CANONICAL_LENGTH_UNIT)).unwrap(),
        );
        self.built.borrow_mut().insert(key, (points.clone(), 1));
        Ok(points)
    }
}

#[test]
fn shared_reference_points_select_what_per_structure_points_select() {
    // 1AKE as many structures in one frame: each protein chain cut in two,
    // the ligands, and the water.
    let complex = model("mmcif/1AKE.cif");
    let part = |text: &str| {
        let atoms = parse(text)
            .evaluate_structure(&complex, SelectEvaluateOptions::default())
            .unwrap();
        complex.subset(&atoms).unwrap()
    };
    let models = [
        part("protein and chain A and resid 1-100"),
        part("protein and chain A and resid 101-300"),
        part("protein and chain B and resid 1-100"),
        part("protein and chain B and resid 101-300"),
        part("ligand"),
        part("water"),
    ];
    let expressions = [
        "within 5 of ligand",
        "within 5 of protein",
        "within 4 of protein and not within 2 of protein",
        "within 3.5 of (water and within 3.5 of ligand)",
        "byres (within 5 of ligand) and protein",
        "within 0 of chain A",
    ];
    let built = RefCell::new(BTreeMap::new());
    let mut selected = 0;
    for expression in expressions {
        for current in 0..models.len() {
            let apart = Scene {
                models: &models,
                current,
            }
            .select(expression)
            .unwrap();
            let shared = SharedScene {
                models: &models,
                current,
                built: &built,
            }
            .select(expression);
            assert_eq!(shared, apart, "{expression} on structure {current}");
            selected += shared.len();
        }
    }
    assert!(selected > 0);
    // Each reference was built once and reused by every other structure.
    for (reference, (_, count)) in built.borrow().iter() {
        assert!(*count >= models.len(), "{reference}: {count}");
    }
    assert_eq!(
        built
            .borrow()
            .keys()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        [
            "chain A",
            "ligand",
            "protein",
            "water and within 3.5 of ligand"
        ]
    );
}

#[test]
fn shared_host_points_measure_each_cutoff_exactly() {
    let topology = smiles::to_topology("C.C.C.C").unwrap();
    let points = (0..4)
        .map(|index| Point3::new(f64::from(index), 0.5 * f64::from(index), 0.0))
        .collect::<Vec<_>>();
    let model = Model::new(
        topology,
        Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
    )
    .unwrap();
    let others = vec![Point3::new(3.2, 1.5, 0.0), Point3::new(-1.0, 0.0, 0.0)];
    let resolver = Resolver {
        other: Positions::new(Quantity::new(others, ANGSTROM))
            .unwrap()
            .into(),
    };
    // One set of points, measured at several cutoffs in turn, two of them
    // equal: each answers as a fresh set of points does.
    for cutoff in [1.0, 0.25, 2.5, 1.0, 0.0, 10.0] {
        let text = format!("within {cutoff} of #2");
        let shared = select::parse(&text, host_options())
            .unwrap()
            .evaluate_structure(
                &model,
                SelectEvaluateOptions::default().with_resolver(&resolver),
            )
            .unwrap()
            .indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let expected = within_points(
            model.positions().values().value(),
            resolver.other.positions().values().value(),
            cutoff,
        );
        assert_eq!(shared, expected, "{text}");
    }
    assert_eq!(resolver.other.len(), 2);
}

#[test]
fn resolvers_without_reference_points_keep_todays_results() {
    let model = model("mmcif/1AKE.cif");
    let terms_only = Resolver {
        other: ReferencePoints::default(),
    };
    for (structure, expression, _) in cases() {
        if structure != "1AKE" {
            continue;
        }
        let evaluate = |options| {
            parse(expression)
                .evaluate_structure(&model, options)
                .unwrap()
        };
        assert_eq!(
            evaluate(SelectEvaluateOptions::default().with_resolver(&terms_only)),
            evaluate(SelectEvaluateOptions::default()),
            "{expression}"
        );
    }
}

/// Records what evaluation asks of it.
#[derive(Default)]
struct Recorder {
    calls: RefCell<Vec<String>>,
}

impl SelectionResolver for Recorder {
    fn resolve(
        &self,
        term: &HostTerm<'_>,
        request: HostRequest<'_>,
    ) -> Result<HostValue, HostTermError> {
        self.calls
            .borrow_mut()
            .push(format!("resolve {} as {:?}", term.text(), request.role()));
        Ok(HostValue::Atoms(AtomSelection::empty(
            request.shared_topology(),
        )))
    }

    fn reference_points(
        &self,
        reference: &SelectionExpr,
        request: HostRequest<'_>,
    ) -> Result<ReferencePoints, HostTermError> {
        self.calls
            .borrow_mut()
            .push(format!("points of {reference} as {:?}", request.role()));
        Ok(ReferencePoints::default())
    }
}

#[test]
fn host_term_references_still_resolve_through_resolve() {
    let topology = smiles::to_topology("C.C.C").unwrap();
    let points = (0..3)
        .map(|index| Point3::new(10.0 * f64::from(index), 0.0, 0.0))
        .collect::<Vec<_>>();
    let model = Model::new(
        topology,
        Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
    )
    .unwrap();
    let recorder = Recorder::default();
    let selected = parse("within 2 of #2 or within 3 of (#2) or within 4 of (not $site)")
        .evaluate_structure(
            &model,
            SelectEvaluateOptions::default().with_resolver(&recorder),
        )
        .unwrap();
    assert_eq!(selected.len(), 3);
    assert_eq!(
        *recorder.calls.borrow(),
        [
            "resolve 2 as Reference",
            "resolve 2 as Reference",
            "resolve site as Atoms",
            "points of not $site as Reference",
        ]
    );
}

#[test]
fn reference_point_failures_name_the_reference_and_its_span() {
    // The same molecule twice, perceived only in the evaluated model.
    let aniline = |perceived: bool| {
        let mut topology = smiles::to_topology("c1ccccc1N").unwrap();
        if perceived {
            topology = Arc::new(topology.perceived().unwrap());
        }
        let points = (0..7)
            .map(|index| Point3::new(1.4 * f64::from(index), 0.0, 0.0))
            .collect::<Vec<_>>();
        Model::new(
            topology,
            Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
        )
        .unwrap()
    };
    let models = [aniline(true), aniline(false)];
    let scene = Scene {
        models: &models,
        current: 0,
    };

    // The reference evaluates on the evaluated model and fails where the
    // resolver evaluates it, on the unperceived one.
    let text = "element O or within 5 of (element N and aromatic)";
    let error = scene.select(text).unwrap_err();
    let SelectEvaluateError::Host { term, span, source } = &error else {
        panic!("{error:?}");
    };
    assert_eq!(term, "element N and aromatic");
    assert_eq!(&text[span.clone()], "(element N and aromatic)");
    assert_eq!(source.message(), "evaluating on model 1 failed");
    // The host's evaluation error locates its term in the same text.
    let cause = std::error::Error::source(source)
        .and_then(|cause| cause.downcast_ref::<SelectEvaluateError>())
        .expect("the resolver attaches its evaluation error");
    assert!(matches!(
        cause,
        SelectEvaluateError::MissingPerception { .. }
    ));
    assert_eq!(&text[cause.span()], "aromatic");

    // A nested reference fails first, in the evaluated model's own pass,
    // and names itself.
    let text = "within 5 of (element C and within 3 of aromatic)";
    let error = scene.select(text).unwrap_err();
    assert!(matches!(&error, SelectEvaluateError::Host { term, .. } if term == "aromatic"));
    assert_eq!(&text[error.span()], "aromatic");
}

#[test]
#[ignore = "timing of a one-million-atom spatial query; run in release with --ignored"]
fn within_on_a_million_atoms_takes_under_fifty_milliseconds() {
    let water = smiles::to_molecules("[H]O[H]").unwrap().remove(0);
    let ligand = smiles::to_molecules("c1ccccc1C(=O)O").unwrap().remove(0);
    let mut builder = kekule::topology::TopologyBuilder::new();
    let ligand_instance = builder.add_molecule(ligand).unwrap();
    let definition = builder.add_molecule_definition(water).unwrap();
    let waters = (1_000_000 - builder.atom_count()) / 3;
    builder.reserve_instances(waters).unwrap();
    let instances = (0..waters)
        .map(|_| builder.add_instance(definition).unwrap())
        .collect::<Vec<_>>();
    let hierarchy = builder.hierarchy_mut();
    let chain = hierarchy.add_chain("L", Some("L".into())).unwrap();
    let residue = hierarchy
        .add_residue(chain, "LIG", None, Some("1".into()), None)
        .unwrap();
    for atom in 0..9 {
        hierarchy
            .add_atom_site(
                residue,
                kekule::topology::InstanceAtomId::new(
                    ligand_instance,
                    kekule::molecule::AtomId::new(atom),
                ),
                Default::default(),
            )
            .unwrap();
    }
    let solvent = hierarchy.add_chain("W", Some("W".into())).unwrap();
    for (number, instance) in instances.into_iter().enumerate() {
        let residue = hierarchy
            .add_residue(solvent, "HOH", None, Some(number.to_string()), None)
            .unwrap();
        for atom in 0..3 {
            hierarchy
                .add_atom_site(
                    residue,
                    kekule::topology::InstanceAtomId::new(
                        instance,
                        kekule::molecule::AtomId::new(atom),
                    ),
                    Default::default(),
                )
                .unwrap();
        }
    }
    let topology = builder.build().unwrap();
    // A 21 nm cube at water density.
    let side = 100usize;
    let points = (0..topology.atom_count())
        .map(|index| {
            let cell = index / 3;
            let (x, y, z) = (cell % side, (cell / side) % side, cell / (side * side));
            let offset = (index % 3) as f64 * 0.1;
            Point3::new(x as f64 * 3.1 + offset, y as f64 * 3.1, z as f64 * 3.1)
        })
        .collect::<Vec<_>>();
    let model = Model::new(
        topology,
        Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
    )
    .unwrap();
    let expression = parse("within 5 of resname LIG");
    let mut best = std::time::Duration::MAX;
    let mut count = 0;
    for _ in 0..5 {
        let start = std::time::Instant::now();
        count = expression
            .evaluate_structure(&model, SelectEvaluateOptions::default())
            .unwrap()
            .len();
        best = best.min(start.elapsed());
    }
    eprintln!("{} atoms: {count} selected in {best:?}", model.atom_count());
    assert!(count >= 9);
    if !cfg!(debug_assertions) {
        assert!(best < std::time::Duration::from_millis(50), "took {best:?}");
    }
}

#[test]
#[ignore = "timing of `within` over many structures; run in release with --ignored"]
fn within_over_many_structures_grows_with_atoms_plus_points() {
    // 27 copies of 1AKE without water, about 92k atoms, as a host holds them:
    // apart on a grid, and stacked in one place as an ensemble's members are.
    let complex = model("mmcif/1AKE.cif");
    let dry = complex
        .subset(
            &parse("not water")
                .evaluate_structure(&complex, SelectEvaluateOptions::default())
                .unwrap(),
        )
        .unwrap();
    let points: &[Point3] = dry.positions().values().value();
    let copies = |spacing: f64| {
        (0..27)
            .map(|copy| {
                let (x, y, z) = (copy % 3, (copy / 3) % 3, copy / 9);
                let shift = Point3::new(
                    f64::from(x) * spacing,
                    f64::from(y) * spacing,
                    f64::from(z) * spacing,
                );
                let moved = points
                    .iter()
                    .map(|p| Point3::new(p.x + shift.x, p.y + shift.y, p.z + shift.z))
                    .collect();
                Model::new(
                    dry.shared_topology(),
                    Positions::from_vec(Quantity::new(moved, CANONICAL_LENGTH_UNIT)).unwrap(),
                )
                .unwrap()
            })
            .collect::<Vec<_>>()
    };
    let time = |run: &dyn Fn() -> usize| {
        let mut best = std::time::Duration::MAX;
        let mut count = 0;
        for _ in 0..10 {
            let start = std::time::Instant::now();
            count = run();
            best = best.min(start.elapsed());
        }
        (best, count)
    };
    // 1AKE spans under 7 nm.
    for (layout, spacing) in [("apart", 10.0), ("stacked", 0.0)] {
        let models = copies(spacing);
        let atoms = models.len() * points.len();
        for expression in [
            "within 5 of protein",
            "byres (within 5 of ligand) and protein",
        ] {
            let shared = time(&|| {
                let built = RefCell::new(BTreeMap::new());
                (0..models.len())
                    .map(|current| {
                        SharedScene {
                            models: &models,
                            current,
                            built: &built,
                        }
                        .select(expression)
                        .len()
                    })
                    .sum()
            });
            // Points gathered and gridded anew for each structure.
            let apart = time(&|| {
                (0..models.len())
                    .map(|current| {
                        Scene {
                            models: &models,
                            current,
                        }
                        .select(expression)
                        .unwrap()
                        .len()
                    })
                    .sum()
            });
            eprintln!(
                "{layout}, {atoms} atoms, `{expression}`: shared {:?}, per structure {:?}",
                shared.0, apart.0
            );
            assert_eq!(shared.1, apart.1, "{layout} {expression}");
            assert!(shared.1 > 0);
            if !cfg!(debug_assertions) {
                assert!(
                    shared.0 < std::time::Duration::from_millis(50),
                    "{layout} {expression}: {:?}",
                    shared.0
                );
                // Gridding every other copy's protein for each copy costs
                // copies times points; shared points cost atoms plus points.
                if expression == "within 5 of protein" {
                    assert!(shared.0 * 3 < apart.0, "{layout}: {shared:?} {apart:?}");
                }
            }
        }
    }
}
