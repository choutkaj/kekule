//! A decoded hierarchy keeps each distinct row of residue labels and of
//! site metadata once, and every node reads the labels it was written with.

use std::collections::HashSet;

use crate::molecule::AtomId;
use crate::native::{self, NativeReadOptions, NativeWriteOptions};
use crate::smiles;
use crate::topology::{AtomSiteMetadata, InstanceAtomId, Topology, TopologyBuilder};

/// Waters in two chains whose residues and sites repeat labels, with
/// absent, empty and non-ASCII strings side by side.
fn waters() -> Topology {
    let water = smiles::to_molecules("[H]O[H]").unwrap().remove(0);
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(water).unwrap();
    let instances = (0..12)
        .map(|_| builder.add_instance(definition).unwrap())
        .collect::<Vec<_>>();
    let hierarchy = builder.hierarchy_mut();
    let chains = [
        hierarchy.add_chain("A", Some("A".into())).unwrap(),
        hierarchy.add_chain("B", None).unwrap(),
    ];
    for (index, instance) in instances.into_iter().enumerate() {
        let chain = chains[index % 2];
        let residue = hierarchy
            .add_residue(
                chain,
                "HOH",
                (index % 4 < 2).then_some(1),
                Some((index % 4).to_string()),
                (index % 3 == 0).then(|| "Å".to_owned()),
            )
            .unwrap();
        let names = [Some("H1"), Some("Oα"), None];
        for (atom, name) in names.into_iter().enumerate() {
            let symbol = match atom {
                1 => Some("O".to_owned()),
                2 => Some(String::new()),
                _ => None,
            };
            let metadata = AtomSiteMetadata::default()
                .with_type_symbol(symbol)
                .with_label_asym_id(Some(["A", "B"][index % 2].to_owned()))
                .with_auth_asym_id((index % 2 == 0).then(|| "A".to_owned()))
                .with_label_atom_id(name.map(str::to_owned))
                .with_auth_atom_id((index % 5 == 0).then(|| "ö".to_owned()));
            hierarchy
                .add_atom_site(
                    residue,
                    InstanceAtomId::new(instance, AtomId::new(atom as u32)),
                    metadata,
                )
                .unwrap();
        }
    }
    builder.build().unwrap()
}

#[test]
fn decoded_nodes_share_rows_of_equal_labels() {
    let source = waters();
    let bytes = native::write(&source, NativeWriteOptions::default()).unwrap();
    let decoded = native::read(bytes.as_slice(), NativeReadOptions::default())
        .unwrap()
        .into_topology()
        .unwrap();
    assert!(
        native::write(&decoded, NativeWriteOptions::default()).unwrap() == bytes,
        "re-encoding the decoded topology changed its bytes"
    );
    let (hierarchy, written) = (decoded.hierarchy(), source.hierarchy());
    assert_eq!(hierarchy.residue_ids().len(), written.residue_ids().len());
    for residue in written.residue_ids() {
        assert_eq!(
            hierarchy.residue_labels(residue),
            written.residue_labels(residue),
            "labels of {residue:?}"
        );
    }
    assert_eq!(hierarchy.site_atom_ids(), written.site_atom_ids());
    for site in written.atom_site_ids() {
        assert_eq!(
            hierarchy.site_metadata(site),
            written.site_metadata(site),
            "metadata of {site:?}"
        );
    }
    let distinct_labels = written
        .residue_ids()
        .map(|residue| format!("{:?}", written.residue_labels(residue)))
        .collect::<HashSet<_>>();
    let distinct_metadata = written
        .atom_site_ids()
        .map(|site| format!("{:?}", written.site_metadata(site)))
        .collect::<HashSet<_>>();
    assert_eq!(hierarchy.label_rows().len(), distinct_labels.len());
    assert_eq!(hierarchy.metadata_rows().len(), distinct_metadata.len());
    assert!(distinct_labels.len() < written.residue_ids().len());
    assert!(distinct_metadata.len() < written.atom_site_ids().len());
}
