//! Disulfides that no record states, found from the positions of cysteine
//! `SG` atoms.
use std::collections::{BTreeSet, HashMap};

use crate::diagnostics::{Cancellation, FormatError};
use crate::molecule::{BondOrder, Element};

use super::atom_site::AtomRow;
use super::connections::{DeclaredConnection, DeclaredCoordination, InstanceUnion};
use super::connectivity::within_bonding_distance;
use super::{
    poll, poll_every, BondEvidence, SourceResidueId, StructureInterpretIssue,
    StructureInterpretationReport,
};

/// The residue names of cysteines that can form a disulfide: the dictionary's
/// and the AMBER name of a cysteine in one. A thiolate, such as AMBER `CYM`,
/// is not one.
const CYSTEINES: [&str; 2] = ["CYS", "CYX"];

/// Bonds the `SG` atoms of two cysteines that lie within bonding distance,
/// when no record bonds either atom to another residue and neither carries a
/// hydrogen. Each `SG` takes its nearest partner, closest pairs first. The
/// bonds join the cysteines' instances and are reported.
pub(super) fn infer_disulfides(
    rows: &[&AtomRow],
    bonds: &mut Vec<DeclaredConnection>,
    coordination: &[DeclaredCoordination],
    union: &mut InstanceUnion,
    report: &mut StructureInterpretationReport,
    cancellation: Option<&Cancellation>,
) -> Result<(), FormatError> {
    let residue = |row: &AtomRow| (row.asym_id.clone(), row.residue_key.clone());
    // An SG that a record links or coordinates outside its residue already
    // has its partner.
    let linked = bonds
        .iter()
        .filter(|bond| !bond.supplements_template)
        .flat_map(|bond| bond.rows)
        .collect::<BTreeSet<_>>();
    let coordinated = coordination
        .iter()
        .flat_map(|record| record.atoms.iter().cloned())
        .collect::<BTreeSet<_>>();
    let mut sulfurs = Vec::new();
    let mut hydrogens = HashMap::<_, Vec<_>>::new();
    for (step, &row) in rows.iter().enumerate() {
        poll_every(cancellation, step)?;
        let Some(point) = row.point else {
            continue;
        };
        if !CYSTEINES.contains(&row.comp_id.to_ascii_uppercase().as_str()) {
            continue;
        }
        if row.element.atomic_number() == 1 {
            hydrogens.entry(residue(row)).or_default().push(point);
        } else if row.atom_name == "SG"
            && row.element == Element::from_symbol("S").expect("sulfur")
            && !linked.contains(&row.row_index)
            && !coordinated.contains(&row.atom_key())
        {
            sulfurs.push((row, point));
        }
    }
    let hydrogen = Element::from_symbol("H").expect("hydrogen");
    sulfurs.retain(|&(row, point)| {
        !hydrogens.get(&residue(row)).is_some_and(|points| {
            points
                .iter()
                .any(|&other| within_bonding_distance((row.element, point), (hydrogen, other)))
        })
    });
    let mut pairs = Vec::new();
    for (index, &(left, left_point)) in sulfurs.iter().enumerate() {
        poll(cancellation)?;
        for &(right, right_point) in &sulfurs[index + 1..] {
            if residue(left) != residue(right)
                && within_bonding_distance((left.element, left_point), (right.element, right_point))
            {
                pairs.push(((left_point - right_point).norm(), left, right));
            }
        }
    }
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut bonded = BTreeSet::new();
    for (_, left, right) in pairs {
        if bonded.contains(&left.row_index) || bonded.contains(&right.row_index) {
            continue;
        }
        bonded.extend([left.row_index, right.row_index]);
        union.union(&left.instance_key, &right.instance_key);
        bonds.push(DeclaredConnection {
            rows: [left.row_index, right.row_index],
            instance_key: left.instance_key.clone(),
            order: BondOrder::Single,
            evidence: BondEvidence::Disulfide,
            supplements_template: false,
        });
        report
            .issues
            .push(StructureInterpretIssue::DisulfideInferred {
                residues: [SourceResidueId::of(left), SourceResidueId::of(right)],
            });
    }
    Ok(())
}
