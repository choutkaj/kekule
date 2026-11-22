//! Bonds an mmCIF reader derives without `_struct_conn` records: the bonds
//! of residues their default template, bundled or force-field, rebuilds
//! exactly, and the peptide and phosphodiester links between residues
//! adjacent in sequence.
//!
//! Each test mirrors the reader's connectivity resolution, so a bond is
//! left unstated only when reading the block gives it back with its order.
use super::prepare::{AtomRow, ConnectionRow, Implied};
use crate::components::{Component, ComponentLibrary};
use crate::io::macromolecular::connectivity::{
    dictionary_template, hydrogen_sites, is_coordination, link_valence, reader_protonation,
    reader_template, standard_link, within_bonding_distance, AtomEvidence, PolymerLink,
    Protonation, ResidueEvidence, StandardLink,
};
use crate::math::Point3;
use crate::molecule::{BondOrder, Element};
use crate::topology::InstanceAtomId;
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

/// What a reader makes of one residue of the block: the rows of its
/// template atoms, and the protonation state it gives the residue, or
/// `None` when the residue takes no template or its hydrogens fit none.
pub(super) struct ReadResidue<'c> {
    members: Vec<usize>,
    /// The template the reader applies, if any.
    component: Option<&'c Component>,
    /// Whether the residue's atoms fit the template, so that the reader
    /// adds its bonds.
    fits: bool,
    /// Rows by template atom name; `None` when a row names no template
    /// atom of its element, which the reader rejects.
    matched: Option<BTreeMap<&'c str, usize>>,
    state: Option<Protonation<'c>>,
}

/// The block's rows and connections as a reader sees them: residues,
/// neighbours, and coordination.
pub(super) struct ReadBlock<'c> {
    residues: Vec<ReadResidue<'c>>,
    residue_of: Vec<usize>,
    /// Each row's partners, with the connection's index or `None` for
    /// report coordination.
    neighbours: Vec<Vec<(usize, Option<usize>)>>,
}

impl<'c> ReadBlock<'c> {
    /// Mirrors the reader's templates and protonation states for the rows,
    /// their elements in row order, the connections, and the report's
    /// metal coordination, which the block states as `metalc` rows.
    /// `library` is the caller's component library, which readers prefer.
    pub(super) fn new(
        library: Option<&'c ComponentLibrary>,
        rows: &[AtomRow],
        elements: &[Element],
        connections: &[ConnectionRow],
        coordination: &[[InstanceAtomId; 2]],
        row_of: &BTreeMap<InstanceAtomId, usize>,
    ) -> Self {
        let mut grouped = BTreeMap::<ResidueKey<'_>, Vec<usize>>::new();
        for (index, row) in rows.iter().enumerate() {
            grouped.entry(residue_key(row)).or_default().push(index);
        }
        let grouped = grouped.into_values().collect::<Vec<_>>();
        let mut residue_of = vec![0; rows.len()];
        for (residue, members) in grouped.iter().enumerate() {
            for &row in members {
                residue_of[row] = residue;
            }
        }
        let mut neighbours = vec![Vec::new(); rows.len()];
        for (index, connection) in connections.iter().enumerate() {
            let (left, right) = (row_of[&connection.left], row_of[&connection.right]);
            neighbours[left].push((right, Some(index)));
            neighbours[right].push((left, Some(index)));
        }
        for [left, right] in coordination {
            let (left, right) = (row_of[left], row_of[right]);
            neighbours[left].push((right, None));
            neighbours[right].push((left, None));
        }
        let ends = chain_ends(rows, &grouped);
        let hydrogen = |row: usize| elements[row].atomic_number() == 1;
        let residues = grouped
            .into_iter()
            .zip(ends)
            .map(|(members, ends)| {
                let id = &rows[members[0]].label_comp_id;
                let atoms = members
                    .iter()
                    .map(|&row| (rows[row].label_atom_id.as_str(), elements[row]));
                let Some(template) = reader_template(library, id, atoms) else {
                    return ReadResidue {
                        members,
                        component: None,
                        fits: false,
                        matched: None,
                        state: None,
                    };
                };
                let (component, fits) = (template.component, template.fits);
                let mut matched = Some(BTreeMap::new());
                for &row in &members {
                    let fits = component
                        .atom(&rows[row].label_atom_id)
                        .filter(|(_, atom)| atom.element() == elements[row]);
                    let inserted = match (fits, &mut matched) {
                        (Some((_, atom)), Some(map)) => map.insert(atom.name(), row).is_none(),
                        _ => false,
                    };
                    if !inserted {
                        matched = None;
                    }
                }
                let Some(map) = &matched else {
                    return ReadResidue {
                        members,
                        component: Some(component),
                        fits,
                        matched,
                        state: None,
                    };
                };
                let in_residue = members.iter().copied().collect::<BTreeSet<_>>();
                let heavy = map
                    .iter()
                    .filter(|&(_, &row)| !hydrogen(row))
                    .map(|(&name, &row)| {
                        let mut evidence = AtomEvidence {
                            stated: Some(rows[row].formal_charge),
                            ..AtomEvidence::default()
                        };
                        let mut sites = Vec::new();
                        for &(other, connection) in &neighbours[row] {
                            match connection {
                                None => evidence.coordination += 1,
                                Some(_) if hydrogen(other) => sites.push((
                                    rows[other].type_symbol == "D",
                                    rows[other].occupancy,
                                    rows[other].label_alt_id.as_deref(),
                                )),
                                Some(connection) if !in_residue.contains(&other) => {
                                    evidence.links += link_valence(connections[connection].order);
                                }
                                Some(_) => {}
                            }
                        }
                        evidence.hydrogens = hydrogen_sites(sites);
                        (name, evidence)
                    })
                    .collect();
                let evidence = ResidueEvidence {
                    heavy,
                    has_hydrogens: members.iter().any(|&row| hydrogen(row)),
                    dictionary: component.force_field().is_none(),
                    polymer: rows[members[0]].label_seq_id.is_some(),
                };
                let state = reader_protonation(component, template.variants, evidence, ends);
                ReadResidue {
                    members,
                    component: Some(component),
                    fits,
                    matched,
                    state,
                }
            })
            .collect();
        Self {
            residues,
            residue_of,
            neighbours,
        }
    }

    /// The fixed implicit hydrogen count a reader gives each row, by row
    /// index: atoms that bond across a gap or anchor a leaving atom their
    /// residue lacks, in residues it reads with their template.
    pub(super) fn fixed_hydrogens(
        &self,
        elements: &[Element],
        connections: &[ConnectionRow],
    ) -> BTreeMap<usize, u8> {
        let mut fixed = BTreeMap::new();
        for residue in &self.residues {
            let (Some(matched), Some(state)) = (&residue.matched, &residue.state) else {
                continue;
            };
            if !self.rebuilt(residue, connections) {
                continue;
            }
            for atom in state.atoms.iter().filter(|atom| atom.fixed) {
                let Some(&row) = matched.get(atom.name.as_str()) else {
                    continue;
                };
                let explicit = self.neighbours[row]
                    .iter()
                    .filter(|&&(other, connection)| {
                        connection.is_some() && elements[other].atomic_number() == 1
                    })
                    .count();
                fixed.insert(
                    row,
                    u8::try_from(atom.hydrogens.saturating_sub(explicit)).unwrap_or(u8::MAX),
                );
            }
        }
        fixed
    }

    /// Whether reading the residue with its template gives it exactly its
    /// bonds: every atom matches a template atom by name and element, and
    /// its bonds are the template's among those atoms with the localized
    /// orders of the protonation state the reader chooses.
    fn rebuilt(&self, residue: &ReadResidue<'_>, connections: &[ConnectionRow]) -> bool {
        let (Some(matched), Some(state)) = (&residue.matched, &residue.state) else {
            return false;
        };
        let candidate = state.candidate;
        let row_of = |index: usize| matched.get(candidate.atoms()[index].name()).copied();
        let mut expected = BTreeMap::new();
        for bond in candidate.bonds() {
            let [a, b] = bond.atoms();
            let (Some(left), Some(right)) = (row_of(a), row_of(b)) else {
                continue;
            };
            let Some(order) = state.bond_order(bond) else {
                return false;
            };
            expected.insert((left.min(right), left.max(right)), order);
        }
        let in_residue = residue.members.iter().copied().collect::<BTreeSet<_>>();
        let mut actual = BTreeMap::new();
        for &row in &residue.members {
            for &(other, connection) in &self.neighbours[row] {
                if let Some(connection) = connection {
                    if row < other && in_residue.contains(&other) {
                        actual.insert((row, other), connections[connection].order);
                    }
                }
            }
        }
        actual == expected
    }
}

/// Whether each residue has residues before and after it in its polymer
/// chain, as a reader of the block finds: rows of one asym and entity with
/// a `label_seq_id` form a chain, and a position after the first or a
/// modelled residue on that side counts. The block states no sequence
/// beyond its last modelled residue.
fn chain_ends(rows: &[AtomRow], residues: &[Vec<usize>]) -> Vec<[bool; 2]> {
    let mut bounds = BTreeMap::<(&str, &str), (i32, i32)>::new();
    for members in residues {
        let row = &rows[members[0]];
        if let Some(position) = row.label_seq_id {
            let bound = bounds
                .entry((&row.asym_id, &row.entity_id))
                .or_insert((position, position));
            *bound = (bound.0.min(position), bound.1.max(position));
        }
    }
    residues
        .iter()
        .map(|members| {
            let row = &rows[members[0]];
            match (
                row.label_seq_id,
                bounds.get(&(&*row.asym_id, &*row.entity_id)),
            ) {
                (Some(position), Some(&(first, last))) => {
                    [position > 1 || position > first, position < last]
                }
                _ => [false, false],
            }
        })
        .collect()
}

/// Marks each connection the reader derives, given the rows' elements in
/// row order and coordinates written to `precision` decimal places.
/// Connections reference rows through `row_of`. Returns, for each residue
/// whose template the reader applies but whose atoms lack one of its bonds
/// between heavy atoms, the rows of the first such pair: the reader adds
/// the bond, so no `_struct_conn` row can keep it out.
pub(super) fn mark_implied(
    block: &ReadBlock<'_>,
    rows: &[AtomRow],
    elements: &[Element],
    connections: &mut [ConnectionRow],
    row_of: &BTreeMap<InstanceAtomId, usize>,
    precision: usize,
    library: Option<&ComponentLibrary>,
) -> Vec<[usize; 2]> {
    let ends = connections
        .iter()
        .map(|connection| [row_of[&connection.left], row_of[&connection.right]])
        .collect::<Vec<_>>();
    let residue_of = &block.residue_of;
    let mut missing = Vec::new();
    for residue in &block.residues {
        if !block.rebuilt(residue, connections) {
            if let Some(component) = residue.component.filter(|_| residue.fits) {
                missing.extend(missing_template_bond(
                    component,
                    &residue.members,
                    rows,
                    &block.neighbours,
                ));
            }

            continue;
        }
        for &row in &residue.members {
            for &(other, connection) in &block.neighbours[row] {
                if let Some(connection) = connection {
                    if residue_of[other] == residue_of[row] {
                        connections[connection].implied = Some(Implied::Template);
                    }
                }
            }
        }
    }
    let residues = block
        .residues
        .iter()
        .map(|residue| &residue.members)
        .collect::<Vec<_>>();

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
            dictionary_template(library, &row.label_comp_id).map_or(
                crate::components::ComponentLinking::None,
                |(component, _)| component.linking(),
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
        let connections_ref = &*connections;
        let residue_of = &residue_of;
        // The reader's own rule: a linking atom bonded to another residue
        // keeps the link out, and the block states the bond.
        let outcome = PolymerLink::undeclared(linking(left_row), linking(right_row)).map(|link| {
            standard_link(
                link,
                named(residue_of[left]),
                named(residue_of[right]),
                |atom: usize, partner: usize| {
                    Ok::<_, std::convert::Infallible>(block.neighbours[atom].iter().any(
                        |&(other, connection)| {
                            other != partner
                                && residue_of[other] != residue_of[atom]
                                && connection.is_some_and(|connection| {
                                    !is_coordination(
                                        connections_ref[connection].order,
                                        [elements[atom], elements[other]],
                                    )
                                })
                        },
                    ))
                },
                |earlier, later| {
                    within_bonding_distance(
                        (
                            elements[earlier],
                            written(rows[earlier].position, precision),
                        ),
                        (elements[later], written(rows[later].position, precision)),
                    )
                },
            )
        });
        if let Some(Ok(StandardLink::Linked(pair))) = outcome {
            if pair == [left, right] {
                connections[index].implied = Some(Implied::PolymerLink);
            }
        }
    }
    missing
}

/// The rows of the first bond between two heavy atoms of `component`, whose
/// atoms the residue fits, that the residue lacks.
fn missing_template_bond(
    component: &Component,
    members: &[usize],
    rows: &[AtomRow],
    neighbours: &[Vec<(usize, Option<usize>)>],
) -> Option<[usize; 2]> {
    let row_of = |name: &str| {
        members
            .iter()
            .copied()
            .find(|&row| rows[row].label_atom_id == name)
    };
    component.bonds().iter().find_map(|bond| {
        let [a, b] = bond.atoms().map(|atom| &component.atoms()[atom]);
        if a.element().atomic_number() == 1 || b.element().atomic_number() == 1 {
            return None;
        }
        let (left, right) = (row_of(a.name())?, row_of(b.name())?);
        (!neighbours[left]
            .iter()
            .any(|&(other, connection)| other == right && connection.is_some()))
        .then_some([left, right])
    })
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
