const MINIMAL_MMCIF: &str = r#"
data_demo
loop_
_entity.id
_entity.type
1 polymer
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
GLY C1 C2 sing
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
C C1 GLY A 1 1 0.0 0.0 0.0
C C2 GLY A 1 1 1.0 0.0 0.0
"#;

fn with_struct_conn(order: &str) -> String {
    format!(
        r#"{MINIMAL_MMCIF}
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_value_order
duplicate covale A C1 A C2 {order}
"#
    )
}

#[test]
fn mmcif_public_facade_requires_parse_then_interpret() -> Result<(), Box<dyn std::error::Error>> {
    use kekule::mmcif::{
        self, MmcifBlock, MmcifEntityClassifications, MmcifEntityKind, MmcifInterpretOptions,
        MmcifWriteError, MmcifWriteOptions,
    };

    let document = mmcif::parse_str(MINIMAL_MMCIF)?;
    assert_eq!(
        document,
        mmcif::parse_str_with_options(MINIMAL_MMCIF, mmcif::MmcifParseOptions::default())?
    );
    let interpreted = mmcif::interpret(&document, MmcifInterpretOptions::default())?;
    let block: &MmcifBlock = &document.blocks()[0];
    let block_interpreted = mmcif::interpret_block(block, MmcifInterpretOptions::default())?;
    let method_default = block.interpret()?;
    let method_explicit = block.interpret_with_options(MmcifInterpretOptions::default())?;

    assert_eq!(document.blocks().len(), 1);
    assert_eq!(interpreted.report().block_name(), "demo");
    assert_eq!(block_interpreted.report(), interpreted.report());
    assert_eq!(method_default.report(), method_explicit.report());
    assert_eq!(document.interpret()?.report(), method_default.report());
    assert!(block_interpreted
        .topology()
        .same_layout(interpreted.topology()));
    assert_eq!(
        block_interpreted.model().positions(),
        interpreted.model().positions()
    );
    assert_eq!(interpreted.model().topology().instance_count(), 1);
    assert!(!interpreted.model().topology().hierarchy().is_empty());
    assert_eq!(interpreted.model().positions().len(), 2);
    assert_eq!(interpreted.molecules().count(), 1);
    assert!(interpreted
        .molecules()
        .all(|molecule| molecule.perception() == &kekule::core::Perception::default()));
    assert!(interpreted
        .topology()
        .same_layout(interpreted.model().topology()));
    let owned_topology = method_default.clone().to_topology();
    let owned_model = method_default.clone().to_model();
    let owned_molecules = method_default.to_molecules();
    assert!(owned_topology.same_layout(owned_model.topology()));
    assert_eq!(
        owned_topology
            .molecules()
            .map(|occurrence| occurrence.molecule())
            .collect::<Vec<_>>(),
        owned_molecules.iter().collect::<Vec<_>>()
    );
    assert!(block
        .to_topology()?
        .same_layout(block.to_model()?.topology()));
    assert_eq!(block.to_molecules()?.len(), 1);
    assert_eq!(interpreted.report().selected_model(), Some("1"));
    assert_eq!(interpreted.report().instances().len(), 1);
    let provenance = &interpreted.report().instances()[0];
    assert_eq!(provenance.coordinate_model_id(), "1");
    assert_eq!(provenance.atoms().len(), 2);
    assert_eq!(provenance.atoms()[0].atom_name(), "C1");
    assert_eq!(provenance.atoms()[0].label_sequence_id(), None);
    assert_eq!(provenance.atoms()[0].author_sequence_id(), Some("1"));
    assert_eq!(provenance.atoms()[0].insertion_code(), None);
    assert_eq!(provenance.atoms()[0].occurrence(), None);
    assert_eq!(provenance.atoms()[0].selected_alternate_location(), None);
    assert_eq!(
        provenance.atoms()[0].atom().molecule(),
        provenance.molecule()
    );
    let written = mmcif::write_with_report(
        interpreted.model(),
        interpreted.report(),
        MmcifWriteOptions::default(),
    )?;
    assert!(written.starts_with("data_model\n"));
    assert!(mmcif::parse_str(&written).is_ok());
    let custom_written = mmcif::write_with_report(
        interpreted.model(),
        interpreted.report(),
        MmcifWriteOptions {
            block_name: "custom".to_owned(),
            ..MmcifWriteOptions::default()
        },
    )?;
    assert!(custom_written.starts_with("data_custom\n"));
    assert!(matches!(
        mmcif::write_with_report(
            interpreted.model(),
            interpreted.report(),
            MmcifWriteOptions {
                block_name: String::new(),
                ..MmcifWriteOptions::default()
            }
        ),
        Err(MmcifWriteError::InvalidBlockName(name)) if name.is_empty()
    ));
    assert_eq!(
        mmcif::write(interpreted.model(), MmcifWriteOptions::default()),
        Err(MmcifWriteError::MissingEntityClassification(
            provenance.molecule()
        ))
    );
    let mut classifications = MmcifEntityClassifications::new();
    classifications.insert(provenance.molecule(), MmcifEntityKind::Polymer)?;
    let written = mmcif::write_with_classifications(
        interpreted.model(),
        &classifications,
        MmcifWriteOptions::default(),
    )?;
    assert!(written.contains("1 polymer"));
    let model = interpreted.to_model();
    assert_eq!(model.topology().instance_count(), 1);
    assert_eq!(model.positions().len(), 2);

    Ok(())
}

#[test]
fn duplicate_agreeing_authoritative_bond_evidence_is_idempotent() {
    use kekule::mmcif::{self, MmcifInterpretOptions};

    let document = mmcif::parse_str(&with_struct_conn("sing"))
        .expect("duplicate agreeing connectivity parses");
    let interpreted = mmcif::interpret(&document, MmcifInterpretOptions::default())
        .expect("duplicate agreeing connectivity is accepted");
    let graph = interpreted
        .model()
        .topology()
        .definitions()
        .next()
        .expect("one molecule definition")
        .1
        .molecule();

    assert_eq!(graph.bond_count(), 1);
}

#[test]
fn duplicate_conflicting_authoritative_bond_evidence_is_rejected() {
    use kekule::mmcif::{self, MmcifInterpretError, MmcifInterpretOptions};

    let document = mmcif::parse_str(&with_struct_conn("doub"))
        .expect("duplicate conflicting connectivity parses");
    let error: MmcifInterpretError = mmcif::interpret(&document, MmcifInterpretOptions::default())
        .expect_err("conflicting authoritative bond orders must be rejected");

    assert_eq!(error.line(), None);
    assert_eq!(
        error.message(),
        "conflicting authoritative mmCIF bond evidence for one atom pair: Double versus Single"
    );
}
