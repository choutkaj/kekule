use super::*;

fn code(id: &str) -> Option<(String, CodeOrigin)> {
    component_code(id, &SequenceOptions::default())
        .map(|code| (code.as_str().to_owned(), code.origin()))
}

#[test]
fn standard_and_modified_components() {
    assert_eq!(code("ALA"), Some(("A".into(), CodeOrigin::Component)));
    assert_eq!(code("dc"), Some(("C".into(), CodeOrigin::Component)));
    assert_eq!(code("MSE"), Some(("M".into(), CodeOrigin::Component)));
    assert_eq!(code("CRO"), Some(("TYG".into(), CodeOrigin::Component)));
    // The dictionary states no code for 02O, a modified alanine.
    assert_eq!(code("02O"), Some(("A".into(), CodeOrigin::Parent)));
    // A dinucleotide with no code of its own takes both parents' codes.
    assert_eq!(code("05A"), Some(("TT".into(), CodeOrigin::Parent)));
    // 004 links into peptides but has neither a code nor a parent.
    assert_eq!(code("004"), Some(("X".into(), CodeOrigin::Unknown)));
    assert_eq!(code("HOH"), None);
    assert_eq!(code("NOT-A-COMPONENT"), None);
    // Nucleotide ligands state codes or parents but link into no polymer.
    assert_eq!(code("GTP"), None);
    assert_eq!(code("ADP"), None);
}

#[test]
fn codes_name_their_parents() {
    let parents = |id: &str| {
        component_code(id, &SequenceOptions::default()).map(|code| {
            (
                code.parents().collect::<Vec<_>>().join(","),
                code.is_modified(),
            )
        })
    };
    assert_eq!(parents("ALA"), Some((String::new(), false)));
    assert_eq!(parents("MSE"), Some(("MET".into(), true)));
    assert_eq!(parents("SEP"), Some(("SER".into(), true)));
    assert_eq!(parents("CRO"), Some(("THR,TYR,GLY".into(), true)));
    assert_eq!(parents("02O"), Some(("ALA".into(), true)));
    assert_eq!(parents("004"), Some((String::new(), false)));
}

#[test]
fn stated_components_come_first() {
    let stated = BTreeMap::from([
        (
            "ALA".to_owned(),
            StatedComponent {
                linking: Some(ComponentLinking::Peptide),
                ..StatedComponent::default()
            },
        ),
        (
            "XYZ".to_owned(),
            StatedComponent {
                linking: Some(ComponentLinking::Nucleotide),
                parents: vec!["DA".to_owned()],
                ..StatedComponent::default()
            },
        ),
        (
            "QQQ".to_owned(),
            StatedComponent {
                linking: Some(ComponentLinking::Nucleotide),
                ..StatedComponent::default()
            },
        ),
        (
            "LOOP".to_owned(),
            StatedComponent {
                linking: Some(ComponentLinking::Peptide),
                parents: vec!["LOOP".to_owned()],
                ..StatedComponent::default()
            },
        ),
    ]);
    let sources = CodeSources::new(Some(&stated), None);
    // A stated component without a code still takes the libraries' code.
    assert_eq!(
        sources.code("ALA", None).map(|code| code.code),
        Some("A".into())
    );
    let modified = sources.code("xyz", None).expect("stated parent");
    assert_eq!(
        (modified.as_str(), modified.origin()),
        ("A", CodeOrigin::Parent)
    );
    let unknown = sources.code("QQQ", Some("X")).expect("stated nucleotide");
    assert_eq!((unknown.as_str(), unknown.is_unknown()), ("N", true));
    assert_eq!(
        sources.code("LOOP", None).map(|code| code.code),
        Some("X".into())
    );
}
