//! Bonds an mmCIF reader derives without `_struct_conn` records: the bonds
//! of residues their bundled component template rebuilds exactly, and the
//! peptide and phosphodiester links between residues adjacent in sequence.
//!
//! Each test mirrors the reader's connectivity resolution, so a bond is
//! left unstated only when reading the block gives it back with its order.
use super::prepare::{AtomRow, ConnectionRow, Implied};
use crate::components::{Component, ComponentLibrary};
use crate::geometry::Point3;
use crate::io::structure_connectivity::{
    hydrogen_differences, template_bond_order, within_bonding_distance, PolymerLink,
};
use crate::molecule::{BondOrder, Element};
use std::collections::{BTreeMap, BTreeSet};

/// The reader's identity of a residue: rows that share it are one residue.
type ResidueKey<'a> = (
    &'a str,
    &'a str,
    Option<i32>,
    Option<&'a str>,
    Option<&'a str>,
    &'a str,
);

fn residue_key(row: &AtomRow) -> ResidueKey<'_> {
    (
        &row.asym_id,
        &row.entity_id,
        row.label_seq_id,
        row.auth_seq_id.as_deref(),
        row.insertion_code.as_deref(),
        &row.label_comp_id,
    )
}

/// Marks each connection the reader derives, given the rows' elements in
/// row order and coordinates written to `precision` decimal places.
/// Connections reference rows through `row_of`.
pub(super) fn mark_implied(
    rows: &[AtomRow],
    elements: &[Element],
    connections: &mut [ConnectionRow],
    row_of: &BTreeMap<crate::topology::InstanceAtomId, usize>,
    precision: usize,
) {
    let mut residues = BTreeMap::<ResidueKey<'_>, Vec<usize>>::new();
    for (index, row) in rows.iter().enumerate() {
        residues.entry(residue_key(row)).or_default().push(index);
    }
    let mut residue_of = vec![0; rows.len()];
    let residues = residues.into_values().collect::<Vec<_>>();
    for (residue, members) in residues.iter().enumerate() {
        for &row in members {
            residue_of[row] = residue;
        }
    }
    let ends = connections
        .iter()
        .map(|connection| [row_of[&connection.left], row_of[&connection.right]])
        .collect::<Vec<_>>();
    let mut neighbours = vec![Vec::new(); rows.len()];
    for (index, &[left, right]) in ends.iter().enumerate() {
        neighbours[left].push((right, index));
        neighbours[right].push((left, index));
    }

    let library = ComponentLibrary::bundled();
    for members in &residues {
        let Some(component) = library.get(&rows[members[0]].label_comp_id) else {
            continue;
        };
        let internal = |&(other, _): &(usize, usize)| residue_of[other] == residue_of[members[0]];
        if !rebuilt_by_template(component, members, rows, elements, connections, &neighbours) {
            continue;
        }
        for &row in members {
            for &(_, connection) in neighbours[row].iter().filter(|link| internal(link)) {
                connections[connection].implied = Some(Implied::Template);
            }
        }
    }

    // The reader links only residues alone at their position along a chain.
    let mut positions = BTreeMap::<(&str, &str, i32), usize>::new();
    for members in &residues {
        let row = &rows[members[0]];
        if let Some(position) = row.label_seq_id {
            *positions
                .entry((&row.asym_id, &row.entity_id, position))
                .or_default() += 1;
        }
    }
    let alone = |row: &AtomRow, position: i32| {
        positions.get(&(row.asym_id.as_str(), row.entity_id.as_str(), position)) == Some(&1)
    };
    for (index, &[first, second]) in ends.iter().enumerate() {
        if residue_of[first] == residue_of[second] || connections[index].order != BondOrder::Single
        {
            continue;
        }
        let (Some(a), Some(b)) = (rows[first].label_seq_id, rows[second].label_seq_id) else {
            continue;
        };
        let (left, right) = if a.checked_add(1) == Some(b) {
            (first, second)
        } else if b.checked_add(1) == Some(a) {
            (second, first)
        } else {
            continue;
        };
        let (left_row, right_row) = (&rows[left], &rows[right]);
        if left_row.asym_id != right_row.asym_id
            || left_row.entity_id != right_row.entity_id
            || !alone(left_row, a.min(b))
            || !alone(right_row, a.max(b))
        {
            continue;
        }
        let linking = |row: &AtomRow| {
            library.get(&row.label_comp_id).map_or(
                crate::components::ComponentLinking::None,
                Component::linking,
            )
        };
        let residues = &residues;
        let named = |residue: usize| {
            move |name: &str| {
                residues[residue]
                    .iter()
                    .copied()
                    .find(|&row| rows[row].label_atom_id == name)
            }
        };
        let linked = PolymerLink::undeclared(linking(left_row), linking(right_row))
            .and_then(|link| link.atoms(named(residue_of[left]), named(residue_of[right])));
        if linked == Some((left, right))
            && within_bonding_distance(
                (elements[left], written(left_row.position, precision)),
                (elements[right], written(right_row.position, precision)),
            )
        {
            connections[index].implied = Some(Implied::PolymerLink);
        }
    }
}

/// Whether reading the residue with its template gives it exactly its
/// bonds: every atom matches a template atom by name and element, and its
/// bonds are the template's among those atoms with their localized orders,
/// after the protonation variant the reader chooses for its hydrogens, if
/// any, replaces them.
fn rebuilt_by_template(
    component: &Component,
    members: &[usize],
    rows: &[AtomRow],
    elements: &[Element],
    connections: &[ConnectionRow],
    neighbours: &[Vec<(usize, usize)>],
) -> bool {
    let mut matched = BTreeMap::new();
    for &row in members {
        let Some((index, atom)) = component.atom(&rows[row].label_atom_id) else {
            return false;
        };
        if atom.element() != elements[row] || matched.insert(index, row).is_some() {
            return false;
        }
    }
    let mut expected = BTreeMap::new();
    for bond in component.bonds() {
        let [a, b] = bond.atoms();
        let (Some(&left), Some(&right)) = (matched.get(&a), matched.get(&b)) else {
            continue;
        };
        let Some(order) = template_bond_order(bond) else {
            return false;
        };
        expected.insert((left.min(right), left.max(right)), order);
    }
    let in_residue = members.iter().copied().collect::<BTreeSet<_>>();
    let hydrogen = |row: usize| elements[row].atomic_number() == 1;
    let named = matched
        .iter()
        .map(|(&index, &row)| (component.atoms()[index].name(), row))
        .collect::<Vec<_>>();
    let has_hydrogens = members.iter().any(|&row| hydrogen(row));
    let all_zero = |differences: &[(usize, usize, isize)]| {
        differences
            .iter()
            .all(|&(_, _, difference)| difference == 0)
    };
    // As the reader does: the first candidate whose hydrogens match exactly,
    // else the first the hydrogens do not contradict.
    let mut evaluated = Vec::new();
    for candidate in
        std::iter::once(component).chain(ComponentLibrary::bundled().variants(component.id()))
    {
        let differences = hydrogen_differences(candidate, &named, has_hydrogens, hydrogen, |row| {
            let present = neighbours[row]
                .iter()
                .filter(|&&(other, _)| hydrogen(other))
                .count();
            let linked = neighbours[row]
                .iter()
                .any(|&(other, _)| !in_residue.contains(&other));
            Ok::<_, std::convert::Infallible>((present, linked))
        })
        .unwrap_or_else(|never| match never {});
        if let Some(differences) = differences {
            // The reader fails on a difference no charge can hold.
            if differences
                .iter()
                .any(|&(_, _, difference)| i8::try_from(difference).is_err())
            {
                return false;
            }
            evaluated.push((candidate, differences));
        }
    }
    let Some((chosen, _)) = evaluated
        .iter()
        .find(|(_, differences)| all_zero(differences))
        .or_else(|| evaluated.first())
    else {
        return false;
    };
    if chosen.variant().is_some() {
        let row_of = |index: usize| {
            let name = chosen.atoms()[index].name();
            named
                .iter()
                .find(|(candidate, _)| *candidate == name)
                .map(|&(_, row)| row)
        };
        for bond in chosen.bonds() {
            let [a, b] = bond.atoms();
            let (Some(left), Some(right)) = (row_of(a), row_of(b)) else {
                continue;
            };
            let Some(order) = expected.get_mut(&(left.min(right), left.max(right))) else {
                continue;
            };
            let Some(variant_order) = template_bond_order(bond) else {
                return false;
            };
            *order = variant_order;
        }
    }
    let mut actual = BTreeMap::new();
    for &row in members {
        for &(other, connection) in &neighbours[row] {
            if row < other && in_residue.contains(&other) {
                actual.insert((row, other), connections[connection].order);
            }
        }
    }
    actual == expected
}

/// A position as the block writes it, to `precision` decimal places.
fn written(position: Point3, precision: usize) -> Point3 {
    let round = |value: f64| {
        format!("{value:.precision$}")
            .parse()
            .expect("formatted coordinates parse")
    };
    Point3::new(round(position.x), round(position.y), round(position.z))
}
