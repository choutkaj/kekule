#[test]
fn singleton_polymer_scheme_retains_declared_chain_order() {
    let looped =
        "data_x\nloop_\n_pdbx_poly_seq_scheme.asym_id\n_pdbx_poly_seq_scheme.seq_id\nZ 1\n";
    for input in [
        looped.to_owned(),
        crate::test_fixtures::singleton_loops_as_scalars(looped),
    ] {
        let doc = crate::mmcif::parse_str(&input, Default::default()).unwrap();
        let order = super::polymer_asym_order(&doc.blocks()[0]).unwrap();
        assert_eq!(order.len(), 1);
        assert_eq!(order["Z"], 0);
    }
}
