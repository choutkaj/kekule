use super::super::altloc::SourceResiduePosition;
use super::*;

fn residue(sequence: i32) -> SourceResidueId {
    SourceResidueId::from_parts("1", "A", SourceResiduePosition::Label(sequence), None)
}

fn names(names: &[&str]) -> Vec<String> {
    names.iter().map(|name| (*name).to_owned()).collect()
}

#[test]
fn issues_read_as_sentences() {
    use StructureInterpretIssue as Issue;
    let cases = [
        (
            Issue::EntityTypeInferred {
                asym_id: "B".into(),
                kind: EntityKind::NonPolymer,
            },
            r#"Chain "B" has no entity type. Its residues give it type "non-polymer"."#,
        ),
        (
            Issue::ElementsInferredFromAtomNames {
                source_lines: vec![12],
            },
            "The atom record at line 12 has no element symbol. Its element comes from its \
                 atom name.",
        ),
        (
            Issue::ElementsInferredFromAtomNames {
                source_lines: vec![12, 13, 1500],
            },
            "3 atom records from line 12 have no element symbol. Their elements come from \
                 their atom names.",
        ),
        (
            Issue::ConnectionsIgnored {
                connection_type: "hydrog".into(),
                count: 1_204,
            },
            r#"1,204 connection records of type "hydrog" are not read as bonds."#,
        ),
        (
            Issue::ConnectionsIgnored {
                connection_type: "hydrog".into(),
                count: 1,
            },
            r#"1 connection record of type "hydrog" is not read as a bond."#,
        ),
        (
            Issue::ConnectionUnresolved {
                connection_id: Some("covale1".into()),
                connection_type: "covale".into(),
                partner: 2,
                source_line: Some(40),
                reason: ConnectionResolutionReason::NoMatchingAtom,
            },
            r#"Partner 2 of connection "covale1" of type "covale" at line 40 is not found. No atom matches the partner selector."#,
        ),
        (
            Issue::ConnectionAmbiguous {
                connection_id: None,
                connection_type: "LINK".into(),
                partner: 1,
                source_line: None,
                candidates: 2,
                reason: ConnectionResolutionReason::MultipleMatchingAtoms,
            },
            r#"Partner 1 of a connection of type "LINK" matches 2 atoms, not one."#,
        ),
        (
            Issue::CoordinateModelIgnored {
                model_id: "2".into(),
                atom_site_rows: 1_632,
            },
            r#"Coordinate model "2" with 1,632 atom sites is not read."#,
        ),
        (
            Issue::AlternateLocationOmitted {
                atom_name: "CG".into(),
                alt_id: Some("B".into()),
                residue: residue(12),
                source_line: 300,
            },
            r#"Alternate location "B" of atom "CG" at line 300 is not read. The atom is in residue 12 of chain "A" in model "1"."#,
        ),
        (
            Issue::AlternateLocationIncomplete {
                residue: residue(12),
                alternative: "B".into(),
                missing_atoms: names(&["CD", "NE"]),
            },
            r#"Alternate location "B" of residue 12 of chain "A" in model "1" lacks atoms "CD" and "NE". It is not read."#,
        ),
        (
            Issue::AlternateLocationOccupancyMissing {
                residue: residue(12),
                source_lines: vec![300],
            },
            r#"Residue 12 of chain "A" in model "1" has 1 alternate-location atom site without occupancy at line 300."#,
        ),
        (
            Issue::AlternateLocationCombinationRejected {
                residue: residue(12),
                alternative: "1".into(),
                reason: "The combination has no alternate location of this residue.".into(),
            },
            r#"The alternate-location combination "1" of residue 12 of chain "A" in model "1" is not read. The combination has no alternate location of this residue."#,
        ),
        (
            Issue::AlternateLocationOccupancyTie {
                residues: vec![residue(12), residue(13)],
                alternatives: names(&["A", "B"]),
            },
            r#"Alternate locations "A" and "B" of residue 12 of chain "A" in model "1" and residue 13 of chain "A" in model "1" have equal occupancy. The first by label is kept."#,
        ),
        (
            Issue::ConnectivityCandidatesInferred {
                atom_count: 2_410,
                candidate_count: 3,
            },
            "A molecule of 2,410 atoms has 3 unbonded atom pairs within bonding distance.",
        ),
        (
            Issue::ResidueBondsInferred {
                residue: residue(401),
                component_id: "LIG".into(),
                bonds: 24,
            },
            r#"Residue 401 of chain "A" in model "1" (component "LIG") has no template or bond table. Its bonds come from atom distances: 24 single bonds."#,
        ),
        (
            Issue::HydrogenAttachedByDistance {
                residue: residue(401),
                atom_name: "H12".into(),
                heavy_atom_name: "C5".into(),
            },
            r#"Hydrogen "H12" of residue 401 of chain "A" in model "1" has no stated bond. It is bonded to "C5", the nearest heavy atom."#,
        ),
        (
            Issue::ChainBreak {
                residue: residue(12),
                next: residue(13),
            },
            r#"Residue 12 of chain "A" in model "1" and residue 13 of chain "A" in model "1" are not bonded. Their linking atoms are too far apart."#,
        ),
        (
            Issue::ChargesFromBonds {
                residue: residue(5),
                component_id: "7MG".into(),
                charged_atoms: vec![("N7".into(), 1)],
            },
            r#"Residue 5 of chain "A" in model "1" (component "7MG") has no template charges. Atom "N7" (+1) takes charges from its bonds."#,
        ),
        (
            Issue::ChargesFromBonds {
                residue: residue(5),
                component_id: "XYZ".into(),
                charged_atoms: vec![("N1".into(), 1), ("B2".into(), -1)],
            },
            r#"Residue 5 of chain "A" in model "1" (component "XYZ") has no template charges. Atoms "N1" (+1) and "B2" (-1) take charges from their bonds."#,
        ),
        (
            Issue::ProtonationFromHydrogens {
                residue: residue(57),
                component_id: "HIS".into(),
                variant: Some("HIS_LSN".into()),
                charged_atoms: Vec::new(),
            },
            r#"The hydrogens of residue 57 of chain "A" in model "1" (component "HIS") differ from its template. It is read as variant "HIS_LSN"."#,
        ),
        (
            Issue::ProtonationFromHydrogens {
                residue: residue(15),
                component_id: "LYS".into(),
                variant: None,
                charged_atoms: vec![("NZ".into(), 0)],
            },
            r#"The hydrogens of residue 15 of chain "A" in model "1" (component "LYS") differ from its template. Atom "NZ" (0) takes its charge from the hydrogens."#,
        ),
        (
            Issue::ResidueFragmented {
                residue: residue(7),
                fragments: 2,
            },
            r#"The atoms of residue 7 of chain "A" in model "1" lie in 2 molecules."#,
        ),
        (
            Issue::ForceFieldTemplate {
                residue: residue(22),
                component_id: "CYX".into(),
                force_field: "AMBER".into(),
            },
            r#"Residue 22 of chain "A" in model "1" (component "CYX") has no dictionary template. It is read with the AMBER template of that name."#,
        ),
        (
            Issue::DisulfideInferred {
                residues: [residue(22), residue(157)],
            },
            r#"Residue 22 of chain "A" in model "1" and residue 157 of chain "A" in model "1" have no stated disulfide. Their SG atoms are within bonding distance, so they are bonded."#,
        ),
        (
            Issue::CellOmitted {
                reason: CellOmissionReason::SpaceGroup {
                    name: "P 21 21 21".into(),
                },
                location: SourceLocation::at_line(8),
            },
            r#"The unit cell at line 8 is not read. The coordinates are the asymmetric unit of space group "P 21 21 21"."#,
        ),
    ];
    for (issue, text) in cases {
        assert_eq!(issue.to_string(), text, "{issue:?}");
    }
}
