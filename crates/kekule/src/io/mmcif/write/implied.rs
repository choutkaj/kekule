//! Bonds an mmCIF reader derives without `_struct_conn` records: the bonds
//! of residues their default template, bundled or force-field, rebuilds
//! exactly, the bonds of components no template covers, which the block's
//! own `_chem_comp_bond` table states, and the peptide and phosphodiester
//! links between residues adjacent in sequence.
//!
//! Each test mirrors the reader's connectivity resolution, so a bond is
//! left unstated only when reading the block gives it back with its order.
use super::prepare::{AtomRow, ConnectionRow, Implied};
use crate::components::ComponentLinking;
use crate::components::{Component, ComponentLibrary};
use crate::io::macromolecular::connectivity::{
    dictionary_template, hydrogen_sites, is_coordination, link_valence, reader_protonation,
    reader_state, reader_template, standard_link, within_bonding_distance, AtomEvidence,
    PolymerLink, Protonation, ReaderState, ResidueEvidence, StandardLink,
};
use crate::io::macromolecular::disulfides::{is_cysteine, pair_sulfurs, Sulfur};
use crate::io::macromolecular::is_metal;
use crate::math::Point3;
use crate::molecule::{BondOrder, Element};
use crate::topology::InstanceAtomId;
use std::collections::{BTreeMap, BTreeSet};
use std::convert::Infallible;

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
    /// What the state follows from: the template's variants, the
    /// evidence, and the residue's chain ends.
    inputs: Option<(Vec<&'c Component>, ResidueEvidence<'c>, [bool; 2])>,
}

/// The block's rows and connections as a reader sees them: residues,
/// neighbours, and coordination.
pub(super) struct ReadBlock<'c> {
    /// The caller's component library, which readers prefer.
    library: Option<&'c ComponentLibrary>,
    residues: Vec<ReadResidue<'c>>,
    residue_of: Vec<usize>,
    /// Each row's partners, with the connection's index or `None` for
    /// report coordination.
    neighbours: Vec<Vec<(usize, Option<usize>)>>,
}

impl<'c> ReadBlock<'c> {
    /// Mirrors the reader's templates and protonation states for the rows,
    /// their elements in row order, the connections, and the report's
    /// metal coordination, which the block states as `metalc` rows, with
    /// coordinates written to `precision` decimal places and polymer
    /// entities of `lengths` monomers, by entity ID. `library` is the
    /// caller's component library, which readers prefer.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        library: Option<&'c ComponentLibrary>,
        rows: &[AtomRow],
        elements: &[Element],
        connections: &[ConnectionRow],
        coordination: &[[InstanceAtomId; 2]],
        row_of: &BTreeMap<InstanceAtomId, usize>,
        precision: usize,
        lengths: &BTreeMap<String, i32>,
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
        let ends = chain_ends(rows, &grouped, lengths);
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
                        inputs: None,
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
                        inputs: None,
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
                                    Some(written(rows[other].position, precision)),
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
                let state = reader_protonation(
                    component,
                    template.variants.iter().copied(),
                    evidence.clone(),
                    ends,
                );
                ReadResidue {
                    members,
                    component: Some(component),
                    fits,
                    matched,
                    state,
                    inputs: Some((template.variants, evidence, ends)),
                }
            })
            .collect();
        Self {
            library,
            residues,
            residue_of,
            neighbours,
        }
    }

    /// The fixed implicit hydrogen count a reader gives each row, by row
    /// index: atoms that bond across a gap or anchor a leaving atom their
    /// residue lacks, in residues it reads with their template.
    pub(super) fn fixed_hydrogens(&self, elements: &[Element]) -> BTreeMap<usize, u8> {
        let mut fixed = BTreeMap::new();
        for residue in &self.residues {
            let (Some(matched), Some(state)) = (&residue.matched, &residue.state) else {
                continue;
            };
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

    /// The rows whose charge a reader derives differently, with that
    /// charge, in row order: atoms of residues it reads with their template,
    /// given the bond orders `connections` the block states
    /// ([`reader_state`]).
    pub(super) fn unread_charges(
        &self,
        rows: &[AtomRow],
        connections: &[ConnectionRow],
    ) -> BTreeMap<usize, i8> {
        let mut unread = BTreeMap::new();
        for (index, residue) in self.residues.iter().enumerate() {
            let (Some(component), Some(matched), Some((variants, evidence, ends))) =
                (residue.component, &residue.matched, &residue.inputs)
            else {
                continue;
            };
            let name_of = matched
                .iter()
                .map(|(&name, &row)| (row, name))
                .collect::<BTreeMap<_, _>>();
            let mut stated = Vec::new();
            for &row in &residue.members {
                for &(other, connection) in &self.neighbours[row] {
                    let Some(connection) = connection else {
                        continue;
                    };
                    if other <= row
                        || self.residue_of[other] != index
                        || connections[connection].implied.is_some()
                    {
                        continue;
                    }
                    let (Some(&a), Some(&b)) = (name_of.get(&row), name_of.get(&other)) else {
                        continue;
                    };
                    let in_template =
                        component
                            .atom(a)
                            .zip(component.atom(b))
                            .is_some_and(|((a, _), (b, _))| {
                                component.bonds().iter().any(|bond| {
                                    let [x, y] = bond.atoms();
                                    (x, y) == (a, b) || (x, y) == (b, a)
                                })
                            });
                    if in_template {
                        stated.push(([a, b], connections[connection].order));
                    }
                }
            }
            let chosen;
            let state = if stated.is_empty() {
                residue.state.as_ref()
            } else {
                chosen = reader_state(component, variants, &stated, evidence.clone(), *ends);
                match &chosen {
                    ReaderState::Chosen { state, .. } => Some(state),
                    // Readers keep the stated charges, or reject the
                    // residue.
                    ReaderState::OrdersExcludeAll | ReaderState::Unfit => None,
                }
            };
            let Some(state) = state else {
                continue;
            };
            for atom in &state.atoms {
                let Some(&row) = matched.get(atom.name.as_str()) else {
                    continue;
                };
                if rows[row].formal_charge != atom.charge {
                    unread.insert(row, atom.charge);
                }
            }
        }
        unread
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

/// Whether each residue has residues before and after it in its entity's
/// sequence of `lengths` monomers, by entity ID, as a reader of the block
/// finds from `_entity_poly_seq`.
fn chain_ends(
    rows: &[AtomRow],
    residues: &[Vec<usize>],
    lengths: &BTreeMap<String, i32>,
) -> Vec<[bool; 2]> {
    residues
        .iter()
        .map(|members| {
            let row = &rows[members[0]];
            match (row.label_seq_id, lengths.get(&row.entity_id)) {
                (Some(position), Some(&length)) => [position > 1, position < length],
                _ => [false, false],
            }
        })
        .collect()
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

/// What [`mark_implied`] finds beyond the connections the reader derives.
pub(super) struct Implications {
    /// The component tables the block states.
    pub(super) tables: Vec<ComponentTable>,
    /// For each residue whose template the reader applies but whose atoms
    /// lack one of its bonds between heavy atoms, the rows of the first such
    /// pair: the reader adds the bond, so no `_struct_conn` row can keep it
    /// out.
    pub(super) missing: Vec<[usize; 2]>,
    /// The rows of the polymer links the reader makes that no connection
    /// states: bonds the structure lacks.
    pub(super) unstated_links: Vec<[usize; 2]>,
}

/// Marks each connection the reader derives, given the rows' elements in
/// row order and coordinates written to `precision` decimal places.
/// `polymer` tells the rows of polymer entities.
pub(super) fn mark_implied(
    block: &ReadBlock<'_>,
    rows: &[AtomRow],
    elements: &[Element],
    connections: &mut [ConnectionRow],
    precision: usize,
    polymer: impl Fn(&AtomRow) -> bool,
) -> Implications {
    let residue_of = &block.residue_of;
    // The reader keys component tables by upper-case identifier.
    let mut untemplated = BTreeMap::<String, Vec<usize>>::new();
    let mut templated = BTreeSet::new();
    let mut missing = Vec::new();
    for (index, residue) in block.residues.iter().enumerate() {
        let id = rows[residue.members[0]].label_comp_id.to_ascii_uppercase();
        if residue.component.is_none() {
            untemplated.entry(id).or_default().push(index);
            continue;
        }
        templated.insert(id);
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
        .map(|residue| residue.members.as_slice())
        .collect::<Vec<_>>();

    let mut tables = Vec::new();
    for (id, of_component) in untemplated {
        // A template covers some residues of the component.
        if templated.contains(&id) {
            continue;
        }
        let Some(table) = component_table(
            &of_component,
            &residues,
            residue_of,
            rows,
            elements,
            connections,
            &block.neighbours,
        ) else {
            continue;
        };
        for &residue in &of_component {
            for &row in residues[residue] {
                for &(other, connection) in &block.neighbours[row] {
                    if let (Some(connection), true) = (connection, residue_of[other] == residue) {
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

    // A connection between the atoms the reader links is derived; a link
    // no connection states is a bond the structure lacks.
    let mut unstated_links = Vec::new();
    for [left, right] in polymer_links(block, rows, elements, connections, precision, polymer) {
        let stated = block.neighbours[left]
            .iter()
            .find_map(|&(other, connection)| connection.filter(|_| other == right));
        match stated {
            Some(index) if connections[index].order == BondOrder::Single => {
                connections[index].implied = Some(Implied::PolymerLink);
            }
            Some(_) => {}
            None => unstated_links.push([left, right]),
        }
    }
    Implications {
        tables,
        missing,
        unstated_links,
    }
}

/// The rows the reader links between residues alone at consecutive
/// positions along the chain of a polymer entity, the earlier residue's
/// first, by the reader's own rule: a linking atom bonded to another
/// residue keeps the link out, and the atoms lie within bonding distance at
/// coordinates written to `precision` decimal places.
fn polymer_links(
    block: &ReadBlock<'_>,
    rows: &[AtomRow],
    elements: &[Element],
    connections: &[ConnectionRow],
    precision: usize,
    polymer: impl Fn(&AtomRow) -> bool,
) -> Vec<[usize; 2]> {
    let residues = block
        .residues
        .iter()
        .map(|residue| residue.members.as_slice())
        .collect::<Vec<_>>();
    let mut placed = BTreeMap::<(&str, &str, i32), Vec<usize>>::new();
    for (residue, members) in residues.iter().enumerate() {
        let row = &rows[members[0]];
        if let Some(position) = row.label_seq_id {
            placed
                .entry((&row.asym_id, &row.entity_id, position))
                .or_default()
                .push(residue);
        }
    }
    let linking = |row: &AtomRow| {
        dictionary_template(block.library, &row.label_comp_id)
            .map_or(ComponentLinking::None, |(component, _)| component.linking())
    };
    let named = |residue: usize| {
        let members = residues[residue];
        move |name: &str| {
            members
                .iter()
                .copied()
                .find(|&row| rows[row].label_atom_id == name)
        }
    };
    let residue_of = &block.residue_of;
    let mut links = Vec::new();
    for (&(asym, entity, position), earlier) in &placed {
        let &[earlier] = earlier.as_slice() else {
            continue;
        };
        let Some(&[later]) = position
            .checked_add(1)
            .and_then(|next| placed.get(&(asym, entity, next)))
            .map(Vec::as_slice)
        else {
            continue;
        };
        let (left_row, right_row) = (&rows[residues[earlier][0]], &rows[residues[later][0]]);
        if !polymer(left_row) {
            continue;
        }
        let Some(link) = PolymerLink::undeclared(linking(left_row), linking(right_row)) else {
            continue;
        };
        let outcome = standard_link(
            link,
            named(earlier),
            named(later),
            |atom: usize, partner: usize| {
                Ok::<_, Infallible>(block.neighbours[atom].iter().any(|&(other, connection)| {
                    other != partner
                        && residue_of[other] != residue_of[atom]
                        && connection.is_some_and(|connection| {
                            !is_coordination(
                                connections[connection].order,
                                [elements[atom], elements[other]],
                            )
                        })
                }))
            },
            |left, right| {
                within_bonding_distance(
                    (elements[left], written(rows[left].position, precision)),
                    (elements[right], written(rows[right].position, precision)),
                )
            },
        );
        if let Ok(StandardLink::Linked(pair)) = outcome {
            links.push(pair);
        }
    }
    links
}

/// The rows of the disulfides the reader infers that no connection
/// states: between cysteine `SG` atoms that no `_struct_conn` or `metalc`
/// row joins to another residue and no hydrogen of their residue lies
/// within bonding distance of, at coordinates written to `precision`
/// decimal places.
pub(super) fn unstated_disulfides(
    block: &ReadBlock<'_>,
    rows: &[AtomRow],
    elements: &[Element],
    connections: &[ConnectionRow],
    precision: usize,
) -> Vec<[usize; 2]> {
    let residue_of = &block.residue_of;
    let joined = |row: usize| {
        block.neighbours[row].iter().any(|&(other, connection)| {
            residue_of[other] != residue_of[row]
                && connection.is_none_or(|connection| connections[connection].implied.is_none())
        })
    };
    let mut hydrogens = BTreeMap::<usize, Vec<Point3>>::new();
    let mut candidates = Vec::new();
    for (index, row) in rows.iter().enumerate() {
        if !is_cysteine(&row.label_comp_id) {
            continue;
        }
        let point = written(row.position, precision);
        if elements[index].atomic_number() == 1 {
            hydrogens.entry(residue_of[index]).or_default().push(point);
        } else if row.label_atom_id == "SG" && elements[index] == Element::S && !joined(index) {
            candidates.push((index, point));
        }
    }
    let mut sulfurs = candidates
        .into_iter()
        .filter(|&(index, point)| {
            !hydrogens.get(&residue_of[index]).is_some_and(|points| {
                points.iter().any(|&hydrogen| {
                    within_bonding_distance((Element::S, point), (Element::H, hydrogen))
                })
            })
        })
        .map(|(index, point)| Sulfur {
            residue: residue_of[index],
            order: index,
            point,
        })
        .collect::<Vec<_>>();
    let Ok(pairs) = pair_sulfurs(&mut sulfurs, || Ok::<_, Infallible>(()));
    pairs
        .into_iter()
        .filter(|&[left, right]| {
            !block.neighbours[left]
                .iter()
                .any(|&(other, connection)| other == right && connection.is_some())
        })
        .collect()
}

/// The `_chem_comp_bond` rows that give every residue of one component its
/// bonds, or `None` when no table can: residues that bond two names
/// differently, or one that leaves unbonded two names another bonds,
/// alternate locations, which repeat names, or a dative bond without a
/// metal, which the coordination flag cannot state. A component without
/// internal bonds needs no table.
fn component_table(
    of_component: &[usize],
    residues: &[&[usize]],
    residue_of: &[usize],
    rows: &[AtomRow],
    elements: &[Element],
    connections: &[ConnectionRow],
    neighbours: &[Vec<(usize, Option<usize>)>],
) -> Option<Vec<ComponentTableBond>> {
    let mut table = BTreeMap::<[&str; 2], (BondOrder, Option<bool>)>::new();
    for &residue in of_component {
        let mut names = BTreeSet::new();
        for &row in residues[residue] {
            if rows[row].label_alt_id.is_some() || !names.insert(&rows[row].label_atom_id) {
                return None;
            }
            for &(other, connection) in &neighbours[row] {
                let Some(connection) = connection else {
                    continue;
                };
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
            if !neighbours[left]
                .iter()
                .any(|&(other, connection)| other == right && connection.is_some())
            {
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
