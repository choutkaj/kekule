//! Disulfides that no record states, found from the positions of cysteine
//! `SG` atoms.
use std::collections::{BTreeSet, HashMap};

use crate::diagnostics::{Cancellation, FormatError};
use crate::molecule::{BondOrder, Element};

use super::atom_site::AtomRow;
use super::connections::{DeclaredConnection, DeclaredCoordination, InstanceUnion};
use super::connectivity::{
    covalent_bonding_radius, within_bonding_distance, OrderSource, COVALENT_BOND_TOLERANCE_ANGSTROM,
};
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
    let residue_of = rows
        .iter()
        .map(|&row| (row.row_index, residue(row)))
        .collect::<HashMap<_, _>>();
    let linked = bonds
        .iter()
        .filter(|bond| {
            let [left, right] = bond.rows.map(|row| residue_of.get(&row));
            left.is_none() || left != right
        })
        .flat_map(|bond| bond.rows)
        .collect::<BTreeSet<_>>();
    let coordinated = coordination
        .iter()
        .flat_map(|record| record.atoms.iter().cloned())
        .collect::<BTreeSet<_>>();
    let sulfur = Element::from_symbol("S").expect("sulfur");
    let mut sulfurs = Vec::new();
    let mut hydrogens = HashMap::<_, Vec<_>>::new();
    for (step, &row) in rows.iter().enumerate() {
        poll_every(cancellation, step)?;
        let Some(point) = row.point else {
            continue;
        };
        if !CYSTEINES
            .iter()
            .any(|name| row.comp_id.eq_ignore_ascii_case(name))
        {
            continue;
        }
        if row.element.atomic_number() == 1 {
            hydrogens.entry(residue(row)).or_default().push(point);
        } else if row.atom_name == "SG"
            && row.element == sulfur
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
    // Sorted along x, each SG meets only those within bonding reach on x.
    let reach = 2.0 * covalent_bonding_radius(sulfur).0 + COVALENT_BOND_TOLERANCE_ANGSTROM;
    sulfurs.sort_by(|a, b| a.1.x.total_cmp(&b.1.x));
    let mut pairs = Vec::new();
    for (index, &(left, left_point)) in sulfurs.iter().enumerate() {
        poll(cancellation)?;
        for &(right, right_point) in &sulfurs[index + 1..] {
            if right_point.x - left_point.x > reach {
                break;
            }
            if residue(left) != residue(right)
                && within_bonding_distance((left.element, left_point), (right.element, right_point))
            {
                let (first, second) = if left.row_index < right.row_index {
                    (left, right)
                } else {
                    (right, left)
                };
                pairs.push(((left_point - right_point).norm(), first, second));
            }
        }
    }
    // Ties fall to source order, so the result does not depend on the sweep.
    pairs.sort_by(|a, b| {
        a.0.total_cmp(&b.0)
            .then_with(|| a.1.row_index.cmp(&b.1.row_index))
            .then_with(|| a.2.row_index.cmp(&b.2.row_index))
    });
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
            order_source: OrderSource::Stated,
            evidence: BondEvidence::Disulfide,
            within_residue: false,
            unstated_metal_bond: false,
        });
        report
            .issues
            .push(StructureInterpretIssue::DisulfideInferred {
                residues: [SourceResidueId::of(left), SourceResidueId::of(right)],
            });
    }
    Ok(())
}
