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
C C1 GLY A 1 1
C C2 GLY A 1 1
"#;

#[test]
fn mmcif_read_str_public_facade() -> Result<(), Box<dyn std::error::Error>> {
    use molecules::mmcif::{self, MmcifParseOptions};

    let macro_mol = mmcif::read_str(MINIMAL_MMCIF, MmcifParseOptions::default())?;

    assert_eq!(macro_mol.graph().atom_count(), 2);
    assert_eq!(macro_mol.graph().bond_count(), 0);
    macro_mol.validate()?;

    Ok(())
}

#[test]
fn bio_mmcif_reexport_remains_available() -> Result<(), Box<dyn std::error::Error>> {
    use molecules::bio::{read_mmcif_str, MmcifParseOptions};

    let macro_mol = read_mmcif_str(MINIMAL_MMCIF, MmcifParseOptions::default())?;

    assert_eq!(macro_mol.graph().atom_count(), 2);

    Ok(())
}

#[test]
fn mmcif_document_and_interpretation_public_facade() -> Result<(), Box<dyn std::error::Error>> {
    use molecules::mmcif::{self, MmcifInterpretOptions, MmcifParseOptions};

    let document = mmcif::parse_str(MINIMAL_MMCIF, MmcifParseOptions::default())?;
    let interpreted = mmcif::interpret(&document, MmcifInterpretOptions::default())?;

    assert_eq!(document.blocks().len(), 1);
    assert_eq!(interpreted.contents().macromolecules().count(), 1);
    assert_eq!(interpreted.contents().small_molecules().count(), 0);
    assert!(interpreted.contents().solvent().is_empty());

    Ok(())
}
