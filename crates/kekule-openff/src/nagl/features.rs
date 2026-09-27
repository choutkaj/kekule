//! Configurable atom features using the versioned NAGL normalization/resonance profile.
use super::config::Feature;
use crate::{error, Result};
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
            .map_err(error)?
            .into_iter()
            .map(|data| {
                Ok(Normalization {
                    query: parse_smarts(&data.query).map_err(error)?,
                    data,
                })
            })
            .collect()
        })
        .as_ref()
        .map_err(Clone::clone)?;
    let mut molecule = input.clone();
    molecule.perceive().map_err(error)?;
    for rule in rules {
        let query = TaggedQuery::new(&rule.query).map_err(error)?;
        let tags: Vec<_> = query.tags().iter().map(|(t, _)| *t).collect();
        let mut converged = false;
        for _ in 0..200 {
            let matches = query
                .find_matches(
                    &PreparedTarget::new(&molecule),
                    crate::assignment::options(),
                    true,
                )
                .map_err(error)?;
            let Some(first) = matches.first() else {
                converged = true;
                break;
            };
            let mapping: BTreeMap<_, _> = tags.iter().copied().zip(first.iter().copied()).collect();
            let mut editor = molecule.edit();
            let mut changed = false;
            for (tag, charge) in &rule.data.charges {
                let id = mapping[tag];
                if editor.atom(id).map_err(error)?.formal_charge != *charge {
                    editor.atom_mut(id).map_err(error)?.formal_charge = *charge;
                    changed = true;
                }
            }
            for &(a, b, order) in &rule.data.bonds {
                let Some(order) = order else {
                    continue;
                };
                let id = editor
                    .bond_between(mapping[&a], mapping[&b])
                    .map_err(error)?
                    .ok_or_else(|| error("normalization changed connectivity"))?;
                let order = match order {
                    1 => BondOrder::Single,
                    2 => BondOrder::Double,
                    3 => BondOrder::Triple,
                    _ => return Err(error("unsupported normalization bond order")),
                };
                if editor.bond(id).map_err(error)?.order != order {
                    editor.bond_mut(id).map_err(error)?.set_order(order);
                    changed = true;
                }
            }
            if !changed {
                converged = true;
                break;
            }
            molecule = editor.finish().map_err(error)?;
            molecule.perceive().map_err(error)?;
        }
        if !converged {
            return Err(error(
                "NAGL normalization did not converge within 200 applications",
            ));
        }
    }
    Ok(molecule)
}

#[derive(Clone)]
struct State {
    charges: Vec<i8>,
    orders: Vec<u8>,
}
struct ResonanceGraph {
    elements: Vec<u8>,
    neighbors: Vec<Vec<(usize, usize)>>,
    edges: Vec<(usize, usize)>,
}
impl ResonanceGraph {
    // Positive => acceptor, negative => donor; the magnitude identifies a
    // conjugate pair in the upstream ResonanceType registry.
    fn kind(&self, state: &State, i: usize) -> i8 {
        let mut bonds = self.neighbors[i]
            .iter()
            .map(|(_, b)| state.orders[*b])
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
    fn energy(&self, state: &State, nodes: &[usize]) -> u32 {
        nodes
            .iter()
            .map(|&i| match self.kind(state, i) {
                -1 | -2 | 3 | -4 | -5 => 5,
                _ => 0,
            })
            .sum()
    }
    // Separate node and bond stacks avoid re-searching edge IDs at every transfer.
    #[allow(clippy::too_many_arguments)]
    fn paths(
        &self,
        state: &State,
        allowed: &BTreeSet<usize>,
        path: &mut Vec<usize>,
        bonds: &mut Vec<usize>,
        acceptor: usize,
        output: &mut Vec<State>,
        work: &mut usize,
    ) -> Result<()> {
        *work += 1;
        if *work > 1_000_000 {
            return Err(error("resonance transfer path limit exceeded"));
        }
        let at = *path.last().unwrap();
        if at == acceptor {
            if bonds.len() >= 2 && bonds.len().is_multiple_of(2) {
                if output.len() >= 10_000
                    || output
                        .len()
                        .saturating_mul(state.charges.len() + state.orders.len())
                        > 20_000_000
                {
                    return Err(error("resonance product limit exceeded"));
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
            if !allowed.contains(&next) || path.contains(&next) {
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
            self.paths(state, allowed, path, bonds, acceptor, output, work)?;
            bonds.pop();
            path.pop();
        }
        Ok(())
    }
}
pub(crate) fn average_charges(input: &Molecule) -> Result<Vec<f32>> {
    if input.atom_count() > 4096 {
        return Err(error("NAGL feature atom limit exceeded (4096)"));
    }
    let molecule = normalize(input)?;
    let ids: Vec<_> = molecule.atom_ids().collect();
    let indices: BTreeMap<_, _> = ids.iter().enumerate().map(|(i, &a)| (a, i)).collect();
    let mut graph = ResonanceGraph {
        elements: ids
            .iter()
            .map(|&a| molecule.atom(a).unwrap().element.atomic_number())
            .collect(),
        neighbors: vec![vec![]; ids.len()],
        edges: vec![],
    };
    let mut initial = State {
        charges: ids
            .iter()
            .map(|&a| molecule.atom(a).unwrap().formal_charge)
            .collect(),
        orders: vec![],
    };
    for (_, b) in molecule.bonds() {
        let order = match b.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => return Err(error("unsupported NAGL bond order")),
        };
        let index = graph.edges.len();
        let a = indices[&b.a()];
        let c = indices[&b.b()];
        graph.edges.push((a, c));
        graph.neighbors[a].push((c, index));
        graph.neighbors[c].push((a, index));
        initial.orders.push(order);
    }
    let mut remaining: BTreeSet<_> = (0..ids.len())
        .filter(|&i| {
            graph.elements[i] != 1
                && !(graph.elements[i] == 6
                    && initial.charges[i] == 0
                    && graph.neighbors[i].len() == 4
                    && graph.neighbors[i]
                        .iter()
                        .all(|(_, b)| initial.orders[*b] == 1))
        })
        .collect();
    let mut average: Vec<_> = initial.charges.iter().map(|&q| f32::from(q)).collect();
    while let Some(&first) = remaining.first() {
        remaining.remove(&first);
        let mut queue = VecDeque::from([first]);
        let mut allowed = BTreeSet::from([first]);
        while let Some(at) = queue.pop_front() {
            for &(n, _) in &graph.neighbors[at] {
                if remaining.remove(&n) {
                    allowed.insert(n);
                    queue.push_back(n);
                }
            }
        }
        let nodes = allowed.iter().copied().collect::<Vec<_>>();
        if !nodes.iter().any(|&i| graph.kind(&initial, i) > 0)
            || !nodes.iter().any(|&i| graph.kind(&initial, i) < 0)
        {
            continue;
        }
        let mut pending = VecDeque::from([initial.clone()]);
        let mut visited = BTreeSet::new();
        let mut unique = BTreeMap::new();
        let mut work = 0;
        while let Some(state) = pending.pop_front() {
            let kinds = nodes
                .iter()
                .map(|&i| graph.kind(&state, i))
                .collect::<Vec<_>>();
            let key = (kinds.clone(), state.orders.clone());
            if !visited.insert(key) {
                continue;
            }
            if visited.len() > 10_000
                || visited.len().saturating_mul(ids.len() + graph.edges.len()) > 20_000_000
            {
                return Err(error("resonance state limit exceeded"));
            }
            // Upstream deduplicates final forms by donor and acceptor positions.
            let roles = kinds.iter().map(|k| k.signum()).collect::<Vec<_>>();
            unique.insert(roles, state.clone());
            for &acceptor in nodes.iter().filter(|&&i| graph.kind(&state, i) > 0) {
                for &donor in nodes.iter().filter(|&&i| graph.kind(&state, i) < 0) {
                    let mut products = vec![];
                    graph.paths(
                        &state,
                        &allowed,
                        &mut vec![donor],
                        &mut vec![],
                        acceptor,
                        &mut products,
                        &mut work,
                    )?;
                    let count = pending.len().saturating_add(products.len());
                    if count > 10_000
                        || count.saturating_mul(ids.len() + graph.edges.len()) > 20_000_000
                    {
                        return Err(error("resonance queue limit exceeded"));
                    }
                    pending.extend(products);
                }
            }
        }
        let min_energy = unique
            .values()
            .map(|s| graph.energy(s, &nodes))
            .min()
            .unwrap();
        let lowest = unique
            .values()
            .filter(|s| graph.energy(s, &nodes) == min_energy)
            .collect::<Vec<_>>();
        for &i in &nodes {
            average[i] = lowest.iter().map(|s| f64::from(s.charges[i])).sum::<f64>() as f32
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
        return Err(error("NAGL requires the selected SSSR ring basis"));
    }
    let averages = if specifications
        .iter()
        .any(|f| matches!(f, Feature::AverageFormalCharge))
    {
        average_charges(molecule)?
    } else {
        Vec::new()
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
                            error(format!(
                                "element {} absent from NAGL feature categories",
                                atom.element
                            ))
                        })?;
                    result[(i, offset + column)] = 1.0;
                }
                Feature::Connectivity { categories } => {
                    let degree = molecule.neighbors(id).map_err(error)?.count();
                    let column = categories
                        .iter()
                        .position(|&n| n == degree)
                        .ok_or_else(|| {
                            error(format!(
                                "connectivity {degree} absent from NAGL feature categories"
                            ))
                        })?;
                    result[(i, offset + column)] = 1.0;
                }
                Feature::AverageFormalCharge => result[(i, offset)] = averages[i],
                Feature::Ring { ring_size } => {
                    result[(i, offset)] = if molecule
                        .ring_set()
                        .unwrap()
                        .rings()
                        .iter()
                        .any(|r| r.atoms.len() == *ring_size && r.atoms.contains(&id))
                    {
                        1.0
                    } else {
                        0.0
                    };
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
                .map_err(error)?
                .map(|b| ids[&b])
                .collect())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    #[test]
    fn atom_limit_precedes_normalization_and_feature_allocation() {
        let molecule = kekule::smiles::to_molecules(&"C".repeat(4097))
            .unwrap()
            .remove(0);
        assert_eq!(
            super::average_charges(&molecule).unwrap_err().to_string(),
            "NAGL feature atom limit exceeded (4096)"
        );
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
