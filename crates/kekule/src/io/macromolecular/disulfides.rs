//! Disulfides that no record states, found from the positions of cysteine
//! `SG` atoms.
use std::collections::{BTreeSet, HashMap};

use crate::diagnostics::{Cancellation, FormatError};
use crate::math::Point3;
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
pub(super) const CYSTEINES: [&str; 2] = ["CYS", "CYX"];

/// Whether a residue named `component` is a cysteine that can form a
/// disulfide ([`CYSTEINES`]), in any case.
pub(crate) fn is_cysteine(component: &str) -> bool {
    CYSTEINES
        .iter()
        .any(|name| component.eq_ignore_ascii_case(name))
}

/// A cysteine `SG` atom that no record links or coordinates outside its
/// residue and no hydrogen of its residue lies within bonding distance of.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Sulfur<R> {
    /// Its residue.
    pub(crate) residue: R,
    /// Its place in source order.
    pub(crate) order: usize,
    /// Its position in ångströms.
    pub(crate) point: Point3,
}

/// The disulfides readers infer among `sulfurs`, as pairs of their
/// [`Sulfur::order`]s, the earlier first, in the order readers bond them.
/// Each `SG` takes its nearest partner of another residue within bonding
/// distance, closest pairs first; ties fall to source order, so the result
/// does not depend on the order of `sulfurs`. Writers share the rule to
/// find the disulfides a reader adds. `poll` runs once per sulfur.
pub(crate) fn pair_sulfurs<R: PartialEq, E>(
    sulfurs: &mut [Sulfur<R>],
    mut poll: impl FnMut() -> Result<(), E>,
) -> Result<Vec<[usize; 2]>, E> {
    // Sorted along x, each SG meets only those within bonding reach on x.
    let reach = 2.0 * covalent_bonding_radius(Element::S).0 + COVALENT_BOND_TOLERANCE_ANGSTROM;
    sulfurs.sort_by(|a, b| a.point.x.total_cmp(&b.point.x));
    let mut pairs = Vec::new();
    for (index, left) in sulfurs.iter().enumerate() {
        poll()?;
        for right in &sulfurs[index + 1..] {
            if right.point.x - left.point.x > reach {
                break;
            }
            if left.residue != right.residue
                && within_bonding_distance((Element::S, left.point), (Element::S, right.point))
            {
                pairs.push((
                    (left.point - right.point).norm(),
                    left.order.min(right.order),
                    left.order.max(right.order),
                ));
            }
        }
    }
    pairs.sort_by(|a, b| {
        a.0.total_cmp(&b.0)
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.cmp(&b.2))
    });
    let mut bonded = BTreeSet::new();
    let mut disulfides = Vec::new();
    for (_, first, second) in pairs {
        if bonded.contains(&first) || bonded.contains(&second) {
            continue;
        }
        bonded.extend([first, second]);
        disulfides.push([first, second]);
    }
    Ok(disulfides)
}

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
    let mut candidates = HashMap::new();
    let mut hydrogens = HashMap::<_, Vec<_>>::new();
    for (step, &row) in rows.iter().enumerate() {
        poll_every(cancellation, step)?;
        let Some(point) = row.point else {
            continue;
        };
        if !is_cysteine(&row.comp_id) {
            continue;
        }
        if row.element.atomic_number() == 1 {
            hydrogens.entry(residue(row)).or_default().push(point);
        } else if row.atom_name == "SG"
            && row.element == Element::S
            && !linked.contains(&row.row_index)
            && !coordinated.contains(&row.atom_key())
        {
            candidates.insert(row.row_index, (row, point));
        }
    }
    let mut sulfurs = candidates
        .values()
        .filter(|&&(row, point)| {
            !hydrogens.get(&residue(row)).is_some_and(|points| {
                points.iter().any(|&other| {
                    within_bonding_distance((row.element, point), (Element::H, other))
                })
            })
        })
        .map(|&(row, point)| Sulfur {
            residue: residue(row),
            order: row.row_index,
            point,
        })
        .collect::<Vec<_>>();
    for [left, right] in pair_sulfurs(&mut sulfurs, || poll(cancellation))? {
        let (left, right) = (candidates[&left].0, candidates[&right].0);
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
