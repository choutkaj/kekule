use crate::bio::{MacroMolecule, SmcraAtomSiteMetadata, SmcraHierarchy};
use crate::core::{Atom, BondOrder, Conformer, Element, Molecule, Point3};
use crate::mmcif::{
    self, MmcifAltLocPolicy, MmcifInterpretOptions, MmcifModelSelection, MmcifParseOptions,
    MmcifWriteError, MmcifWriteOptions,
};
use crate::modeling::{Model, ModelBuilder, MoleculeInstanceMetadata, MoleculeRole};
use crate::small::SmallMolecule;

const MIXED: &str = r#"
data_mixed
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 water
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
L 2
W 3
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 N N GLY A 1 1 0.0 0.0 0.0 1
ATOM 2 C CA GLY A 1 1 1.0 0.0 0.0 1
HETATM 3 C C1 LIG L 2 . 2.0 0.0 0.0 1
HETATM 4 O O HOH W 3 . 3.0 0.0 0.0 1
loop_
_audit_author.name
_audit_author.pdbx_ordinal
'Example Author' 1
"#;

fn parse(input: &str) -> mmcif::MmcifDocument {
    mmcif::parse_str(input, MmcifParseOptions::default()).expect("mmCIF parses")
}

#[test]
fn mmcif_parse_preserves_unknown_categories_without_chemistry() {
    let document = parse(MIXED);
    let block = &document.blocks()[0];
    assert!(block.loop_with_tag("_audit_author.name").is_some());
    assert_eq!(
        block
            .loop_with_tag("_atom_site.type_symbol")
            .unwrap()
            .row_count(),
        4
    );
}

#[test]
fn mmcif_parse_preserves_hashes_inside_bare_values() {
    let document = parse("data_hash\nloop_\n_example.id\n_example.label\n1 sample-d2o#1\n#\n");
    let table = document.blocks()[0]
        .loop_with_tag("_example.label")
        .expect("example loop");
    assert_eq!(table.row_count(), 1);
    assert_eq!(
        table.value(0, "_example.label").map(|value| value.text()),
        Some("sample-d2o#1")
    );
}

#[test]
fn interpretation_preserves_first_source_occurrence_for_model_instances() {
    let input = r#"
data_order
loop_
_entity.id
_entity.type
1 polymer
2 polymer
loop_
_struct_asym.id
_struct_asym.entity_id
Z 1
A 2
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.seq_id
Z 1
A 1
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
covale Z C1 1 A C1 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C C1 GLY A 2 1 0.0 0.0 0.0 1
ATOM 2 C C1 GLY Z 1 1 1.0 0.0 0.0 1
"#;
    let result = mmcif::interpret(&parse(input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.model().topology().molecule_count(), 1);
    let instance = result.model().topology().molecules().next().unwrap().1;
    let hierarchy = instance
        .macro_molecule()
        .expect("merged polymer chains remain a macro molecule")
        .hierarchy();
    let hierarchy_model = hierarchy.models().next().unwrap().1;
    let chain_order = hierarchy_model
        .chains()
        .iter()
        .map(|chain| hierarchy.chain(*chain).unwrap().label_id())
        .collect::<Vec<_>>();
    assert_eq!(chain_order, ["Z", "A"]);
}

#[test]
fn interpretation_builds_distinct_typed_instances_and_complete_positions() {
    let interpreted = mmcif::interpret(&parse(MIXED), MmcifInterpretOptions::default()).unwrap();
    let model = interpreted.model();
    assert_eq!(model.topology().molecule_count(), 3);
    assert_eq!(model.atom_count(), 4);
    assert_eq!(model.positions().len(), 4);
    assert!(model.positions().iter().all(|point| point.x.is_finite()));
    let instances = model
        .topology()
        .molecules()
        .map(|(_, molecule)| molecule)
        .collect::<Vec<_>>();
    assert!(instances[0].macro_molecule().is_some());
    assert!(instances[0].has_role(MoleculeRole::Polymer));
    assert!(instances[1].small_molecule().is_some());
    assert!(instances[1].has_role(MoleculeRole::NonPolymer));
    assert!(instances[2].has_role(MoleculeRole::Solvent));
    assert_eq!(interpreted.report().selected_model.as_deref(), Some("1"));
    assert_eq!(interpreted.report().instances.len(), 3);
    assert_eq!(
        interpreted
            .report()
            .instances
            .iter()
            .map(|instance| instance.atoms.len())
            .sum::<usize>(),
        4
    );
    for instance in &instances {
        assert!(instance.graph().props().is_empty());
        assert!(instance
            .graph()
            .atoms()
            .all(|(_, atom)| atom.props.keys().all(|key| !key.starts_with("mmcif."))));
    }
    let first_provenance = &interpreted.report().instances[0];
    assert_eq!(first_provenance.coordinate_model_id, "1");
    assert_eq!(first_provenance.asym_ids, vec!["A"]);
    assert_eq!(first_provenance.entity_ids, vec!["1"]);
    assert_eq!(first_provenance.atoms[0].atom_name, "N");
    assert_eq!(model.topology().bonds().count(), 1);
    assert_eq!(interpreted.report().inferred_bonds(), 1);
    assert!(interpreted.report().issues().iter().any(|issue| matches!(
        issue,
        mmcif::MmcifInterpretIssue::CovalentBondsInferred {
            atom_count: 2,
            bond_count: 1,
        }
    )));
}

#[test]
fn mmcif_covalent_inference_connects_atoms_across_spatial_cells_without_contacts() {
    let input = MIXED.replace(
        "ATOM 2 C CA GLY A 1 1 1.0 0.0 0.0 1",
        "ATOM 2 C CA GLY A 1 1 1.45 0.0 0.0 1\nATOM 5 C C GLY A 1 1 2.90 0.0 0.0 1",
    );
    let interpreted = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    let polymer = interpreted.model().topology().molecules().next().unwrap().1;

    assert_eq!(polymer.graph().atom_count(), 3);
    assert_eq!(polymer.graph().bond_count(), 2);
    assert_eq!(interpreted.report().inferred_bonds(), 2);
}

#[test]
fn multiple_coordinate_models_require_explicit_selection() {
    let input = MIXED.replace(
        "HETATM 4 O O HOH W 3 . 3.0 0.0 0.0 1",
        "HETATM 4 O O HOH W 3 . 3.0 0.0 0.0 1\nHETATM 5 O O HOH W 3 . 8.0 0.0 0.0 2",
    );
    let document = parse(&input);
    let error = mmcif::interpret(&document, MmcifInterpretOptions::default()).unwrap_err();
    assert!(error.message.contains("select one explicitly"));

    let selected = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            model_selection: MmcifModelSelection::Select("2".into()),
            ..MmcifInterpretOptions::default()
        },
    )
    .unwrap();
    assert_eq!(selected.report().selected_model.as_deref(), Some("2"));
    assert_eq!(selected.model().atom_count(), 1);
    assert_eq!(selected.model().positions()[0].x, 8.0);
    assert_eq!(selected.report().ignored_coordinate_models, vec!["1"]);

    let first = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            model_selection: MmcifModelSelection::First,
            ..MmcifInterpretOptions::default()
        },
    )
    .unwrap();
    assert_eq!(first.report().selected_model.as_deref(), Some("1"));
}

#[test]
fn alternate_location_policy_is_explicit_and_reported() {
    let input = MIXED
        .replace(
            "_atom_site.Cartn_x",
            "_atom_site.label_alt_id\n_atom_site.occupancy\n_atom_site.Cartn_x",
        )
        .replace(
            "ATOM 1 N N GLY A 1 1 0.0 0.0 0.0 1",
            "ATOM 1 N N GLY A 1 1 A 0.4 0.0 0.0 0.0 1\nATOM 5 N N GLY A 1 1 B 0.6 5.0 0.0 0.0 1",
        )
        .replace(" 1.0 0.0 0.0 1", " . 1.0 1.0 0.0 0.0 1")
        .replace(" 2.0 0.0 0.0 1", " . 1.0 2.0 0.0 0.0 1")
        .replace(" 3.0 0.0 0.0 1", " . 1.0 3.0 0.0 0.0 1");
    let document = parse(&input);
    let result = mmcif::interpret(&document, MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.model().positions()[0].x, 5.0);
    assert!(result.report().issues.iter().any(|issue| matches!(
        issue,
        mmcif::MmcifInterpretIssue::AlternateLocationOmitted { alt_id: Some(id), .. } if id == "A"
    )));
    assert!(mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            altloc_policy: MmcifAltLocPolicy::ErrorOnAlternateLocations,
            ..MmcifInterpretOptions::default()
        }
    )
    .is_err());
}

#[test]
fn selected_model_requires_complete_positions() {
    let input = MIXED.replace("3.0 0.0 0.0 1", ". . . 1");
    let error = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap_err();
    assert!(error.message.contains("complete position"));
}

#[test]
fn declared_covalent_links_merge_entities_but_noncovalent_links_do_not() {
    let connections = r#"
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
covale A CA 1 L C1 .
hydrog A N 1 W O .
"#;
    let input = format!("{MIXED}\n{connections}");
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.model().topology().molecule_count(), 2);
    let first = result.model().topology().molecules().next().unwrap().1;
    assert!(first.macro_molecule().is_some());
    assert!(first.has_role(MoleculeRole::Polymer));
    assert!(first.has_role(MoleculeRole::NonPolymer));
    assert_eq!(result.report().applied_connections, 1);
}

#[test]
fn symmetry_mate_connections_are_reported_unresolved() {
    let connection = r#"
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_symmetry
disulf A N 1 1_555 A N 1 15_545
"#;
    let input = format!("{MIXED}\n{connection}");
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    assert_eq!(result.report().applied_connections, 0);
    assert!(result.report().issues.iter().any(|issue| matches!(
        issue,
        mmcif::MmcifInterpretIssue::ConnectionUnresolved { connection_type }
            if connection_type == "disulf"
    )));
}

#[test]
fn struct_conn_bond_order_is_interpreted_and_rejected_when_unknown() {
    let connection = r#"
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.pdbx_value_order
covale A N 1 A CA 1 doub
"#;
    let input = format!("{MIXED}\n{connection}");
    let result = mmcif::interpret(&parse(&input), MmcifInterpretOptions::default()).unwrap();
    let first = result.model().topology().molecules().next().unwrap().1;
    assert_eq!(
        first.graph().bonds().next().expect("declared bond").1.order,
        BondOrder::Double
    );

    let error = mmcif::interpret(
        &parse(&input.replace("doub", "arom")),
        MmcifInterpretOptions::default(),
    )
    .unwrap_err();
    assert!(error
        .message
        .contains("unsupported struct_conn bond order `arom`"));
}

#[test]
fn mmcif_writer_round_trips_supported_model_content() {
    let connection = r#"
loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.pdbx_value_order
covale A N 1 A CA 1 doub
"#;
    let original = mmcif::interpret(
        &parse(&format!("{MIXED}\n{connection}")),
        MmcifInterpretOptions::default(),
    )
    .unwrap();
    let written = mmcif::write(
        original.model(),
        MmcifWriteOptions {
            data_block_name: "round_trip".to_owned(),
            coordinate_precision: 4,
        },
    )
    .expect("supported model should write");
    assert!(written.starts_with("data_round_trip\n"));
    assert!(written.contains("_struct_conn.pdbx_value_order"));
    assert!(written.contains("doub"));

    let document = parse(&written);
    let atom_sites = document.blocks()[0]
        .loop_with_tag("_atom_site.type_symbol")
        .expect("writer emits atom-site loop");
    assert_eq!(atom_sites.row_count(), 4);
    let round_trip = mmcif::interpret(&document, MmcifInterpretOptions::default()).unwrap();
    assert_eq!(round_trip.model().topology().molecule_count(), 3);
    assert_eq!(round_trip.model().positions(), original.model().positions());
    let first = round_trip.model().topology().molecules().next().unwrap().1;
    assert!(first.has_role(MoleculeRole::Polymer));
    assert_eq!(
        first
            .graph()
            .bonds()
            .next()
            .expect("round-trip bond")
            .1
            .order,
        BondOrder::Double
    );
}

#[test]
fn mmcif_writer_rejects_unsupported_chemistry_and_incomplete_hierarchy() {
    let aromatic = small_model_with_bond(BondOrder::Aromatic);
    assert!(matches!(
        mmcif::write(&aromatic, MmcifWriteOptions::default()),
        Err(MmcifWriteError::UnsupportedBondOrder {
            order: BondOrder::Aromatic,
            ..
        })
    ));

    let mut graph = Molecule::new();
    let atom = graph.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            atom,
            crate::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    let conformer = graph.add_conformer(conformer).unwrap();
    let mut hierarchy = SmcraHierarchy::new();
    let model = hierarchy.add_model("1");
    let chain = hierarchy.add_chain(model, "A", None).unwrap();
    hierarchy
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    assert_eq!(
        MacroMolecule::try_from_parts(graph, hierarchy)
            .expect_err("incomplete hierarchy must not construct"),
        crate::bio::MacroValidateError::MissingAtomSiteForAtom { atom }
    );
    let _ = conformer;
}

#[test]
fn mmcif_writer_preserves_supported_bond_orders() {
    for order in [
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Quadruple,
    ] {
        let model = small_model_with_bond(order);
        let written = mmcif::write(&model, MmcifWriteOptions::default()).unwrap();
        let interpreted =
            mmcif::interpret(&parse(&written), MmcifInterpretOptions::default()).unwrap();
        let round_trip = interpreted
            .model()
            .topology()
            .molecules()
            .next()
            .unwrap()
            .1
            .graph()
            .bonds()
            .next()
            .unwrap()
            .1
            .order;
        assert_eq!(round_trip, order);
    }
}

#[test]
fn mmcif_writer_rejects_ambiguous_atom_identity_and_unencodable_roles() {
    let carbon = Element::from_symbol("C").unwrap();
    let mut graph = Molecule::new();
    let left = graph.add_atom(Atom::new(carbon));
    let right = graph.add_atom(Atom::new(carbon));
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            left,
            crate::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            right,
            crate::units::Quantity::new(Point3::new(1.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    let conformer = graph.add_conformer(conformer).unwrap();
    let mut hierarchy = SmcraHierarchy::new();
    let model = hierarchy.add_model("1");
    let chain = hierarchy.add_chain(model, "A", None).unwrap();
    let residue = hierarchy
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    for atom in [left, right] {
        hierarchy
            .add_atom_site(
                residue,
                atom,
                SmcraAtomSiteMetadata {
                    label_atom_id: Some("CA".to_owned()),
                    ..SmcraAtomSiteMetadata::default()
                },
            )
            .unwrap();
    }
    let macro_molecule = MacroMolecule::try_from_parts(graph, hierarchy).unwrap();
    let mut builder = ModelBuilder::new();
    builder
        .add_macro_molecule(&macro_molecule, conformer)
        .unwrap();
    assert!(matches!(
        mmcif::write(&builder.build().unwrap(), MmcifWriteOptions::default()),
        Err(MmcifWriteError::DuplicateAtomIdentity(_))
    ));

    let mut metadata = MoleculeInstanceMetadata::default();
    metadata.insert_role(MoleculeRole::Ligand);
    let model = small_model_with_metadata(metadata);
    assert!(matches!(
        mmcif::write(&model, MmcifWriteOptions::default()),
        Err(MmcifWriteError::UnsupportedMoleculeRole {
            role: MoleculeRole::Ligand,
            ..
        })
    ));
}

fn small_model_with_bond(order: BondOrder) -> Model {
    let mut graph = Molecule::new();
    let left = graph.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
    let right = graph.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
    graph.add_bond(left, right, order).unwrap();
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            left,
            crate::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            right,
            crate::units::Quantity::new(Point3::new(1.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    let conformer = graph.add_conformer(conformer).unwrap();
    let molecule = SmallMolecule::from_graph(graph);
    let mut builder = ModelBuilder::new();
    builder.add_small_molecule(&molecule, conformer).unwrap();
    builder.build().unwrap()
}

fn small_model_with_metadata(metadata: MoleculeInstanceMetadata) -> Model {
    let mut graph = Molecule::new();
    let atom = graph.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            atom,
            crate::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    let conformer = graph.add_conformer(conformer).unwrap();
    let molecule = SmallMolecule::from_graph(graph);
    let mut builder = ModelBuilder::new();
    builder
        .add_small_molecule_with_metadata(&molecule, conformer, metadata)
        .unwrap();
    builder.build().unwrap()
}
