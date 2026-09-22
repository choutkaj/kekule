//! Shared by deterministic regressions and the SMILES fuzz target.
use std::collections::BTreeMap;

use kekule::{core::Molecule, smiles};

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct AtomIdentity {
    element: String,
    isotope: Option<u16>,
    charge: i8,
    radical: Option<(u8, Option<u8>)>,
    map: Option<u32>,
}

// Independent of canonical labeling and atom numbering. Graph and non-graph
// ordinary hydrogens contribute to the same entry; isotope/map/spin information
// remains distinct. Aromatic syntax and fixed/inferred H storage may change.
fn composition(molecule: &Molecule) -> BTreeMap<AtomIdentity, usize> {
    let mut result = BTreeMap::new();
    for (id, atom) in molecule.atoms() {
        *result
            .entry(AtomIdentity {
                element: atom.element.symbol().to_owned(),
                isotope: atom.isotope,
                charge: atom.formal_charge,
                radical: atom
                    .radical
                    .map(|r| (r.electron_count(), r.spin_multiplicity())),
                map: atom.atom_map,
            })
            .or_default() += 1;
        let hydrogens = molecule
            .implicit_hydrogens(id)
            .expect("live atom")
            .expect("perceived H count");
        if hydrogens != 0 {
            *result
                .entry(AtomIdentity {
                    element: "H".to_owned(),
                    isotope: None,
                    charge: 0,
                    radical: None,
                    map: None,
                })
                .or_default() += hydrogens;
        }
    }
    result
}

pub fn assert_output(source: &Molecule, written: &str, canonical: Option<&str>) {
    let document = smiles::parse_str(written).expect("successful writer output must parse");
    let mut restored = document
        .interpret()
        .expect("writer output must interpret")
        .into_molecule()
        .expect("one molecule must remain one component");
    restored
        .perceive()
        .expect("writer output must remain perceivable");
    assert_eq!(
        composition(source),
        composition(&restored),
        "composition changed: {written}"
    );
    if let Some(expected) = canonical {
        assert_eq!(
            smiles::write_canonical(&restored).expect("canonical round trip must write"),
            expected,
            "connectivity/stereo or canonical fixed point changed: {written}"
        );
    }
}

pub fn check_molecule(source: &Molecule) -> usize {
    let canonical = smiles::write_canonical(source);
    let mut checked = 0;
    // Unsupported content and resource exhaustion are valid writer outcomes.
    // Once a writer succeeds, every subsequent round-trip check is mandatory.
    for written in [
        smiles::write(source),
        smiles::write_isomeric(source),
        canonical.clone(),
    ]
    .into_iter()
    .flatten()
    {
        assert_output(source, &written, canonical.as_deref().ok());
        checked += 1;
    }
    checked
}
