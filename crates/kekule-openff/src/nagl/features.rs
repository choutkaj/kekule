//! Configurable atom features using the versioned NAGL normalization/resonance profile.
use super::config::Feature;
use crate::{Error, ErrorKind, Result};

use kekule::{
    core::{AtomId, BondOrder, Molecule, RingBasisModel},
    query::{parse_smarts, QueryGraph},
    substructure::{PreparedTarget, TaggedQuery},
};
use serde::Deserialize;
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    sync::OnceLock,
};

fn limit(detail: impl std::fmt::Display) -> Error {
    Error::new(ErrorKind::ResourceLimit, detail)
}

fn invariant(detail: impl std::fmt::Display) -> Error {
    Error::new(ErrorKind::Chemistry, detail)
}

#[derive(Deserialize)]
struct NormalizationData {
    query: String,
    charges: Vec<(u32, i8)>,
    bonds: Vec<(u32, u32, Option<u8>)>,
}
struct Normalization {
    query: QueryGraph,
    data: NormalizationData,
}
static RULES: OnceLock<Result<Vec<Normalization>>> = OnceLock::new();
pub(crate) fn normalize(input: &Molecule) -> Result<Molecule> {
    let rules = RULES
        .get_or_init(|| {
            serde_json::from_str::<Vec<NormalizationData>>(include_str!(
                "../../data/normalizations.json"
            ))
            .map_err(Error::chemistry)?
            .into_iter()
            .map(|data| {
                Ok(Normalization {
                    query: parse_smarts(&data.query).map_err(Error::chemistry)?,
                    data,
                })
            })
            .collect()
        })
        .as_ref()
        .map_err(Clone::clone)?;
    let mut molecule = input.clone();
    molecule.perceive().map_err(Error::chemistry)?;
    for rule in rules {
        let query = TaggedQuery::new(&rule.query).map_err(Error::chemistry)?;
        let tags: Vec<_> = query.tags().iter().map(|(t, _)| *t).collect();
        // Upstream NAGL stops after 200 applications per rule, which large
        // molecules with many matching groups (arginines, for example) exceed.
        // Each application normalizes one group, so the molecule size bounds
        // every convergent sequence; the guard only catches a cycling rule.
        let applications = molecule.atom_count().max(200);
        let mut converged = false;
        for _ in 0..applications {
            let matches = query
                .find_matches(
                    &PreparedTarget::new(&molecule),
                    crate::assignment::options(molecule.atom_count()),
                    true,
                )
                .map_err(crate::assignment::matching)?;
            let Some(first) = matches.first() else {
                converged = true;
                break;
            };
            let mapping: BTreeMap<_, _> = tags.iter().copied().zip(first.iter().copied()).collect();
            let mut editor = molecule.edit();
            let mut changed = false;
            for (tag, charge) in &rule.data.charges {
                let id = mapping[tag];
                if editor.atom(id).map_err(Error::chemistry)?.formal_charge != *charge {
                    editor.atom_mut(id).map_err(Error::chemistry)?.formal_charge = *charge;
                    changed = true;
                }
            }
            for &(a, b, order) in &rule.data.bonds {
                let Some(order) = order else {
                    continue;
                };
                let id = editor
                    .bond_between(mapping[&a], mapping[&b])
                    .map_err(Error::chemistry)?
                    .ok_or_else(|| invariant("normalization changed connectivity"))?;
                let order = match order {
                    1 => BondOrder::Single,
                    2 => BondOrder::Double,
                    3 => BondOrder::Triple,
                    _ => return Err(invariant("unsupported normalization bond order")),
                };
                if editor.bond(id).map_err(Error::chemistry)?.order != order {
                    editor
                        .bond_mut(id)
                        .map_err(Error::chemistry)?
                        .set_order(order);
                    changed = true;
                }
            }
            if !changed {
                converged = true;
                break;
            }
            molecule = editor.finish().map_err(Error::chemistry)?;
            molecule.perceive().map_err(Error::chemistry)?;
        }
        if !converged {
            return Err(limit(format!(
                "NAGL normalization did not converge within {applications} applications"
            )));
        }
    }
    Ok(molecule)
}

#[derive(Clone)]
struct State {
    /// Formal charges of the fragment's atoms.
    charges: Vec<i8>,
    /// Orders of the fragment's internal bonds.
    orders: Vec<u8>,
}

/// One conjugated fragment with local atom and bond indices.
///
/// Resonance transfers only change charges of fragment atoms and orders of
/// bonds between them, so every state stays local and the work per fragment is
/// independent of the size of the rest of the molecule.
struct ResonanceGraph {
    elements: Vec<u8>,
    /// Internal neighbors and bond indices, in the molecule's bond order.
    neighbors: Vec<Vec<(usize, usize)>>,
    /// Orders of bonds to atoms outside the fragment, which never change.
    fixed: Vec<Vec<u8>>,
}
impl ResonanceGraph {
    // Positive => acceptor, negative => donor; the magnitude identifies a
    // conjugate pair in the upstream ResonanceType registry.
    fn kind(&self, state: &State, i: usize) -> i8 {
        let mut bonds = self.neighbors[i]
            .iter()
            .map(|(_, b)| state.orders[*b])
            .chain(self.fixed[i].iter().copied())
            .collect::<Vec<_>>();
        bonds.sort();
        match (self.elements[i], state.charges[i], bonds.as_slice()) {
            (8, 0, [2]) => 1,
            (8, -1, [1]) => -1,
            (16, 0, [2]) => 2,
            (16, -1, [1]) => -2,
            (7, 1, [1, 1, 2]) => 3,
            (7, 0, [1, 1, 1]) => -3,
            (7, 0, [1, 2]) => 4,
            (7, -1, [1, 1]) => -4,
            (7, 0, [3]) => 5,
            (7, -1, [2]) => -5,
            _ => 0,
        }
    }
    fn energy(&self, state: &State) -> u32 {
        (0..self.elements.len())
            .map(|i| match self.kind(state, i) {
                -1 | -2 | 3 | -4 | -5 => 5,
                _ => 0,
            })
            .sum()
    }
    /// Work per stored state: one charge per atom and one order per bond.
    fn state_size(&self, state: &State) -> usize {
        state.charges.len() + state.orders.len()
    }
    // Separate node and bond stacks avoid re-searching edge IDs at every transfer.
    fn paths(
        &self,
        state: &State,
        path: &mut Vec<usize>,
        bonds: &mut Vec<usize>,
        acceptor: usize,
        output: &mut Vec<State>,
        work: &mut usize,
    ) -> Result<()> {
        *work += 1;
        if *work > 1_000_000 {
            return Err(limit("resonance transfer path limit exceeded"));
        }
        let at = *path.last().unwrap();
        if at == acceptor {
            if bonds.len() >= 2 && bonds.len().is_multiple_of(2) {
                if output.len() >= 10_000
                    || output.len().saturating_mul(self.state_size(state)) > 20_000_000
                {
                    return Err(limit("resonance product limit exceeded"));
                }
                let mut next = state.clone();
                next.charges[path[0]] += 1;
                next.charges[acceptor] -= 1;
                for (i, &b) in bonds.iter().enumerate() {
                    if i % 2 == 0 {
                        next.orders[b] += 1;
                    } else {
                        next.orders[b] -= 1;
                    }
                }
                output.push(next);
            }
            return Ok(());
        }
        for &(next, bond) in &self.neighbors[at] {
            if path.contains(&next) {
                continue;
            }
            if let Some(&previous) = bonds.last() {
                let delta = i16::from(state.orders[bond]) - i16::from(state.orders[previous]);
                if delta != if bonds.len() % 2 == 1 { 1 } else { -1 } {
                    continue;
                }
            }
            path.push(next);
            bonds.push(bond);
            self.paths(state, path, bonds, acceptor, output, work)?;
            bonds.pop();
            path.pop();
        }
        Ok(())
    }
}

/// Average formal charges over the lowest-energy resonance forms of every
/// conjugated fragment, as in the NAGL 0.6.1 preparation profile.
///
/// There is no molecule size limit: fragments are processed independently,
/// and the bounds guard only the combinatorial enumeration within one
/// fragment.
pub(crate) fn average_charges(input: &Molecule) -> Result<Vec<f32>> {
    let molecule = normalize(input)?;
    let ids: Vec<_> = molecule.atom_ids().collect();
    let indices: BTreeMap<_, _> = ids.iter().enumerate().map(|(i, &a)| (a, i)).collect();
    let mut elements = Vec::with_capacity(ids.len());
    let mut charges = Vec::with_capacity(ids.len());
    for &id in &ids {
        let atom = molecule.atom(id).map_err(Error::chemistry)?;
        elements.push(atom.element.atomic_number());
        charges.push(atom.formal_charge);
    }
    let mut neighbors = vec![vec![]; ids.len()];
    let mut orders = Vec::new();
    for (_, b) in molecule.bonds() {
        let order = match b.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => {
                return Err(Error::new(
                    ErrorKind::UnsupportedMolecule,
                    "NAGL features require localized single, double, or triple bonds",
                ))
            }
        };
        let bond = orders.len();
        let a = indices[&b.a()];
        let c = indices[&b.b()];
        neighbors[a].push((c, bond));
        neighbors[c].push((a, bond));
        orders.push(order);
    }
    let mut remaining: BTreeSet<_> = (0..ids.len())
        .filter(|&i| {
            elements[i] != 1
                && !(elements[i] == 6
                    && charges[i] == 0
                    && neighbors[i].len() == 4
                    && neighbors[i].iter().all(|(_, b)| orders[*b] == 1))
        })
        .collect();
    let mut average: Vec<_> = charges.iter().map(|&q| f32::from(q)).collect();
    while let Some(&first) = remaining.first() {
        remaining.remove(&first);
        let mut queue = VecDeque::from([first]);
        let mut fragment = BTreeSet::from([first]);
        while let Some(at) = queue.pop_front() {
            for &(n, _) in &neighbors[at] {
                if remaining.remove(&n) {
                    fragment.insert(n);
                    queue.push_back(n);
                }
            }
        }
        let nodes = fragment.into_iter().collect::<Vec<_>>();
        let local: BTreeMap<usize, usize> =
            nodes.iter().enumerate().map(|(l, &g)| (g, l)).collect();
        // Internal bonds keep the molecule's bond order, and each atom keeps
        // its neighbor order, so enumeration visits states as before.
        let mut internal = BTreeMap::new();
        for &g in &nodes {
            for &(other, bond) in &neighbors[g] {
                if local.contains_key(&other) && !internal.contains_key(&bond) {
                    internal.insert(bond, ());
                }
            }
        }
        let bond_index: BTreeMap<usize, usize> = internal
            .keys()
            .enumerate()
            .map(|(local, &bond)| (bond, local))
            .collect();
        let mut graph = ResonanceGraph {
            elements: nodes.iter().map(|&g| elements[g]).collect(),
            neighbors: vec![vec![]; nodes.len()],
            fixed: vec![vec![]; nodes.len()],
        };
        for (l, &g) in nodes.iter().enumerate() {
            for &(other, bond) in &neighbors[g] {
                match local.get(&other) {
                    Some(&other) => graph.neighbors[l].push((other, bond_index[&bond])),
                    None => graph.fixed[l].push(orders[bond]),
                }
            }
        }
        let initial = State {
            charges: nodes.iter().map(|&g| charges[g]).collect(),
            orders: bond_index.keys().map(|&bond| orders[bond]).collect(),
        };
        let atoms = 0..nodes.len();
        if !atoms.clone().any(|i| graph.kind(&initial, i) > 0)
            || !atoms.clone().any(|i| graph.kind(&initial, i) < 0)
        {
            continue;
        }
        let state_size = graph.state_size(&initial);
        let mut pending = VecDeque::from([initial]);
        let mut visited = BTreeSet::new();
        let mut unique = BTreeMap::new();
        let mut work = 0;
        while let Some(state) = pending.pop_front() {
            let kinds = atoms
                .clone()
                .map(|i| graph.kind(&state, i))
                .collect::<Vec<_>>();
            let key = (kinds.clone(), state.orders.clone());
            if !visited.insert(key) {
                continue;
            }
            if visited.len() > 10_000 || visited.len().saturating_mul(state_size) > 20_000_000 {
                return Err(limit("resonance state limit exceeded"));
            }
            // Upstream deduplicates final forms by donor and acceptor positions.
            let roles = kinds.iter().map(|k| k.signum()).collect::<Vec<_>>();
            unique.insert(roles, state.clone());
            for acceptor in atoms.clone().filter(|&i| kinds[i] > 0) {
                for donor in atoms.clone().filter(|&i| kinds[i] < 0) {
                    let mut products = vec![];
                    graph.paths(
                        &state,
                        &mut vec![donor],
                        &mut vec![],
                        acceptor,
                        &mut products,
                        &mut work,
                    )?;
                    let count = pending.len().saturating_add(products.len());
                    if count > 10_000 || count.saturating_mul(state_size) > 20_000_000 {
                        return Err(limit("resonance queue limit exceeded"));
                    }
                    pending.extend(products);
                }
            }
        }
        let min_energy = unique.values().map(|s| graph.energy(s)).min().unwrap();
        let lowest = unique
            .values()
            .filter(|s| graph.energy(s) == min_energy)
            .collect::<Vec<_>>();
        for (l, &g) in nodes.iter().enumerate() {
            average[g] = lowest.iter().map(|s| f64::from(s.charges[l])).sum::<f64>() as f32
                / lowest.len() as f32;
        }
    }
    Ok(average)
}

pub(super) fn features(
    molecule: &Molecule,
    specifications: &[Feature],
) -> Result<ndarray::Array2<f32>> {
    if molecule.perception().ring_basis_model() != Some(RingBasisModel::FiguerasSssrLike) {
        return Err(invariant("NAGL requires the selected SSSR ring basis"));
    }
    let averages = if specifications
        .iter()
        .any(|f| matches!(f, Feature::AverageFormalCharge))
    {
        average_charges(molecule)?
    } else {
        Vec::new()
    };
    let ring_sizes = if specifications
        .iter()
        .any(|f| matches!(f, Feature::Ring { .. }))
    {
        let mut sizes = BTreeMap::<AtomId, BTreeSet<usize>>::new();
        let rings = molecule
            .ring_set()
            .ok_or_else(|| invariant("NAGL ring features require perceived rings"))?;
        for ring in rings.rings() {
            for &atom in &ring.atoms {
                sizes.entry(atom).or_default().insert(ring.atoms.len());
            }
        }
        sizes
    } else {
        BTreeMap::new()
    };
    let mut result = ndarray::Array2::zeros((
        molecule.atom_count(),
        specifications.iter().map(Feature::width).sum(),
    ));
    for (i, (id, atom)) in molecule.atoms().enumerate() {
        let mut offset = 0;
        for feature in specifications {
            match feature {
                Feature::Element { categories } => {
                    let column = categories
                        .iter()
                        .position(|s| kekule::core::Element::from_symbol(s) == Some(atom.element))
                        .ok_or_else(|| {
                            Error::new(
                                ErrorKind::UnsupportedMolecule,
                                format!(
                                    "element {} absent from NAGL feature categories",
                                    atom.element
                                ),
                            )
                        })?;
                    result[(i, offset + column)] = 1.0;
                }
                Feature::Connectivity { categories } => {
                    let degree = molecule.neighbors(id).map_err(Error::chemistry)?.count();
                    let column = categories
                        .iter()
                        .position(|&n| n == degree)
                        .ok_or_else(|| {
                            Error::new(
                                ErrorKind::UnsupportedMolecule,
                                format!(
                                    "connectivity {degree} absent from NAGL feature categories"
                                ),
                            )
                        })?;
                    result[(i, offset + column)] = 1.0;
                }
                Feature::AverageFormalCharge => result[(i, offset)] = averages[i],
                Feature::Ring { ring_size } => {
                    let member = ring_sizes
                        .get(&id)
                        .is_some_and(|sizes| sizes.contains(ring_size));
                    result[(i, offset)] = if member { 1.0 } else { 0.0 };
                }
            }
            offset += feature.width();
        }
    }
    Ok(result)
}

pub(crate) fn adjacency(molecule: &Molecule) -> Result<Vec<Vec<usize>>> {
    let ids: BTreeMap<AtomId, usize> = molecule
        .atom_ids()
        .enumerate()
        .map(|(i, a)| (a, i))
        .collect();
    molecule
        .atom_ids()
        .map(|a| {
            Ok(molecule
                .neighbors(a)
                .map_err(Error::chemistry)?
                .map(|b| ids[&b])
                .collect())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    fn explicit_hydrogens(smiles: &str) -> kekule::core::Molecule {
        let mut molecule = kekule::smiles::to_molecules(smiles).unwrap().remove(0);
        molecule.perceive().unwrap();
        molecule.add_hydrogens().unwrap();
        molecule
    }

    // Regression: features had a 4096-atom cap, and resonance work grew with
    // the whole molecule for every conjugated fragment.
    #[test]
    fn resonance_averages_scale_past_former_size_cap() {
        // Polyglycine: one amide resonance fragment per residue.
        let molecule = explicit_hydrogens(&format!("N{}CC(=O)O", "CC(=O)N".repeat(800)));
        assert!(molecule.atom_count() > 4096);
        let averages = super::average_charges(&molecule).unwrap();
        assert_eq!(averages.len(), molecule.atom_count());
        // Neutral amides have one dominant form; nothing becomes charged.
        assert!(averages.iter().all(|q| *q == 0.0));
    }

    // Regression: normalization stopped after 200 applications of one rule,
    // as upstream NAGL does, so molecules with more matching groups failed.
    #[test]
    fn normalization_converges_past_two_hundred_matching_groups() {
        // Each zwitterionic imidic group normalizes to a neutral amide.
        let molecule = explicit_hydrogens(&format!("C{}", "C(C([O-])=[NH2+])".repeat(250)));
        assert_eq!(
            molecule
                .atoms()
                .filter(|(_, a)| a.formal_charge != 0)
                .count(),
            500
        );
        let normalized = super::normalize(&molecule).unwrap();
        assert!(normalized.atoms().all(|(_, a)| a.formal_charge == 0));
        assert_eq!(normalized.atom_count(), molecule.atom_count());
    }

    #[test]
    fn complete_ash_features_match_pinned_external_reference() {
        let mut count = 0;
        for r in crate::reference_records() {
            if r["features"]["status"] != "ok" {
                continue;
            }
            let m = kekule::smiles::to_molecules(r["mapped_smiles"].as_str().unwrap())
                .unwrap()
                .remove(0);
            let m = crate::explicit(&m).unwrap();
            let specs = r["features"]["values"]
                .as_array()
                .unwrap()
                .iter()
                .map(|f| serde_json::from_value(f["config"].clone()).unwrap())
                .collect::<Vec<_>>();
            let actual = super::features(&m, &specs).unwrap();
            for (i, (_, atom)) in m.atoms().enumerate() {
                let index = atom.atom_map.unwrap() as usize - 1;
                let expected = r["features"]["values"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .flat_map(|f| f["values"][index].as_array().unwrap())
                    .map(|v| v.as_f64().unwrap())
                    .collect::<Vec<_>>();
                assert_eq!(expected.len(), 22);
                for (j, &expected) in expected.iter().enumerate() {
                    assert!(
                        (f64::from(actual[(i, j)]) - expected).abs() <= 1e-6,
                        "{} atom {} column {j}",
                        r["input"]["id"],
                        index + 1
                    );
                }
            }
            count += 1;
        }
        assert_eq!(count, 22);
    }
}
