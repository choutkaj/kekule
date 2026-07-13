const MINIMAL_MMCIF: &str = r#"
data_demo
loop_
_entity.id
_entity.type
1 polymer
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

#[test]
fn mmcif_public_facade_requires_parse_then_interpret() -> Result<(), Box<dyn std::error::Error>> {
    use molecules::mmcif::{self, MmcifInterpretOptions, MmcifParseOptions};

    let document = mmcif::parse_str(MINIMAL_MMCIF, MmcifParseOptions::default())?;
    let interpreted = mmcif::interpret(&document, MmcifInterpretOptions::default())?;

    assert_eq!(document.blocks().len(), 1);
    assert_eq!(interpreted.model().topology().molecule_count(), 1);
    assert!(interpreted
        .model()
        .topology()
        .molecules()
        .next()
        .unwrap()
        .1
        .macro_molecule()
        .is_some());
    assert_eq!(interpreted.model().positions().len(), 2);

    Ok(())
}
