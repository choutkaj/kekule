//! Bonds an mmCIF reader derives without `_struct_conn` records: the bonds
//! of residues their default template, bundled or force-field, rebuilds
//! exactly, the bonds of components no template covers, which the block's
//! own `_chem_comp_bond` table states, and the peptide and phosphodiester
//! links between residues adjacent in sequence.
//!
//! Each test mirrors the reader's connectivity resolution, so a bond is
//! left unstated only when reading the block gives it back with its order.
use super::prepare::{AtomRow, ConnectionRow, Implied};
use crate::components::{Component, ComponentLibrary};
use crate::io::macromolecular::connectivity::{
    hydrogen_differences, link_valence, template_bond_order, within_bonding_distance, PolymerLink,
};
use crate::io::macromolecular::is_metal;
use crate::math::Point3;
use crate::molecule::{BondOrder, Element};
use std::collections::{BTreeMap, BTreeSet};

/// The reader's identity of a residue: rows that share it are one residue.
pub(super) type ResidueKey<'a> = (
    &'a str,
    &'a str,
    Option<i32>,
    Option<&'a str>,
    Option<&'a str>,
    &'a str,
);

pub(super) fn residue_key(row: &AtomRow) -> ResidueKey<'_> {
    (
        &row.asym_id,
        &row.entity_id,
        row.label_seq_id,
        row.auth_seq_id.as_deref(),
        row.insertion_code.as_deref(),
        &row.label_comp_id,
    )
}

/// The `_chem_comp` row and `_chem_comp_bond` rows of one component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct ComponentTable {
    pub(super) id: String,
    /// `_chem_comp.type`: `non-polymer`, or unknown for a component of a
    /// polymer, whose linking no template states.
    pub(super) polymer: bool,
    pub(super) bonds: Vec<ComponentTableBond>,
}

/// One `_chem_comp_bond` row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct ComponentTableBond {
    pub(super) atoms: [String; 2],
    /// Dative for metal coordination.
    pub(super) order: BondOrder,
    pub(super) aromatic: Option<bool>,
}

/// Marks each connection the reader derives, given the rows' elements in
/// row order and coordinates written to `precision` decimal places, and
/// returns the component tables the block states. Connections reference
/// rows through `row_of`, and `polymer` tells the rows of polymer entities.
pub(super) fn mark_implied(
    rows: &[AtomRow],
    elements: &[Element],
    connections: &mut [ConnectionRow],
    row_of: &BTreeMap<crate::topology::InstanceAtomId, usize>,
    precision: usize,
    polymer: impl Fn(&AtomRow) -> bool,
) -> Vec<ComponentTable> {
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
    // The reader keys component tables by upper-case identifier.
    let mut untemplated = BTreeMap::<String, Vec<usize>>::new();
    let mut templated = BTreeSet::new();
    for (residue, members) in residues.iter().enumerate() {
        let id = &rows[members[0]].label_comp_id;
        let Some(component) = library.get(id).or_else(|| {
            let atoms = members
                .iter()
                .map(|&row| (rows[row].label_atom_id.as_str(), elements[row]));
            ComponentLibrary::force_field_template(id, atoms).map(|(component, _)| component)
        }) else {
            untemplated
                .entry(id.to_ascii_uppercase())
                .or_default()
                .push(residue);
            continue;
        };
        templated.insert(id.to_ascii_uppercase());
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

    let mut tables = Vec::new();
    for (id, of_component) in untemplated {
        // A template covers some residues of the component.
        if templated.contains(&id) {
            continue;
        }
        let Some(table) = component_table(
            &of_component,
            &residues,
            &residue_of,
            rows,
            elements,
            connections,
            &neighbours,
        ) else {
            continue;
        };
        for &residue in &of_component {
            for &row in &residues[residue] {
                for &(other, connection) in &neighbours[row] {
                    if residue_of[other] == residue {
                        connections[connection].implied = Some(Implied::ComponentTable);
                    }
                }
            }
        }
        tables.push(ComponentTable {
            id: rows[residues[of_component[0]][0]].label_comp_id.clone(),
            polymer: of_component
                .iter()
                .any(|&residue| polymer(&rows[residues[residue][0]])),
            bonds: table,
        });
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
    tables
}

/// The `_chem_comp_bond` rows that give every residue of one component its
/// bonds, or `None` when no table can: residues that bond two names
/// differently, or one that leaves unbonded two names another bonds,
/// alternate locations, which repeat names, or a dative bond without a
/// metal, which the coordination flag cannot state. A component without
/// internal bonds needs no table.
fn component_table(
    of_component: &[usize],
    residues: &[Vec<usize>],
    residue_of: &[usize],
    rows: &[AtomRow],
    elements: &[Element],
    connections: &[ConnectionRow],
    neighbours: &[Vec<(usize, usize)>],
) -> Option<Vec<ComponentTableBond>> {
    let mut table = BTreeMap::<[&str; 2], (BondOrder, Option<bool>)>::new();
    for &residue in of_component {
        let mut names = BTreeSet::new();
        for &row in &residues[residue] {
            if rows[row].label_alt_id.is_some() || !names.insert(&rows[row].label_atom_id) {
                return None;
            }
            for &(other, connection) in &neighbours[row] {
                if residue_of[other] != residue || other < row {
                    continue;
                }
                let mut pair = [
                    rows[row].label_atom_id.as_str(),
                    rows[other].label_atom_id.as_str(),
                ];
                pair.sort_unstable();
                let connection = &connections[connection];
                if connection.order == BondOrder::Dative
                    && !is_metal(elements[row])
                    && !is_metal(elements[other])
                {
                    return None;
                }
                let bond = (connection.order, connection.aromatic);
                if *table.entry(pair).or_insert(bond) != bond {
                    return None;
                }
            }
        }
    }
    for &residue in of_component {
        let named = residues[residue]
            .iter()
            .map(|&row| (rows[row].label_atom_id.as_str(), row))
            .collect::<BTreeMap<_, _>>();
        for [left, right] in table.keys() {
            let (Some(&left), Some(&right)) = (named.get(left), named.get(right)) else {
                continue;
            };
            if !neighbours[left].iter().any(|&(other, _)| other == right) {
                return None;
            }
        }
    }
    (!table.is_empty()).then(|| {
        table
            .into_iter()
            .map(|([left, right], (order, aromatic))| ComponentTableBond {
                atoms: [left.to_owned(), right.to_owned()],
                order,
                aromatic,
            })
            .collect()
    })
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
            let links = neighbours[row]
                .iter()
                .filter(|&&(other, _)| !in_residue.contains(&other))
                .map(|&(_, connection)| link_valence(connections[connection].order))
                .sum();
            Ok::<_, std::convert::Infallible>((present, links))
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
