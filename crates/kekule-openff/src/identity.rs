//! Lookup selection and whole-molecule lookup correspondence.
//!
//! Upstream keys its lookup table by fixed-H InChI. Native selection matches
//! an input only when it is exactly an entry's own molecule: the same
//! elements, connectivity, formal charges, bond orders and stereo. Ordinary
//! inputs are written like their entries; other resonance or charge-placement
//! forms that upstream's identifier merges fall back to inference.
use crate::{Error, ErrorKind, Result};

use kekule::core::*;
use std::collections::{BTreeMap, HashMap};

fn error(detail: impl std::fmt::Display) -> Error {
    Error::new(ErrorKind::Identity, detail)
}

/// Net charge and every atom's (element, degree), sorted: equal for any two
/// identical graphs.
type Signature = (i64, Vec<(u8, usize)>);

fn signature(molecule: &Molecule) -> Result<Signature> {
    let mut atoms = molecule
        .atoms()
        .map(|(id, atom)| {
            Ok((
                atom.element.atomic_number(),
                molecule.neighbors(id).map_err(Error::chemistry)?.count(),
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    atoms.sort_unstable();
    Ok((molecule.formal_charge(), atoms))
}

/// Prepared lookup entry graphs, indexed for exact selection.
#[derive(Debug, Default)]
pub(crate) struct LookupIndex {
    entries: HashMap<Signature, Vec<(Molecule, usize)>>,
}

impl LookupIndex {
    /// Entry graphs with their table positions; the first of identical graphs wins.
    pub(crate) fn new(entries: impl IntoIterator<Item = (Molecule, usize)>) -> Result<Self> {
        let mut index = Self::default();
        for (graph, position) in entries {
            index
                .entries
                .entry(signature(&graph)?)
                .or_default()
                .push((graph, position));
        }
        Ok(index)
    }

    /// The table position of the entry whose molecule the prepared input is.
    pub(crate) fn select(&self, molecule: &Molecule) -> Result<Option<usize>> {
        let Some(bucket) = self.entries.get(&signature(molecule)?) else {
            return Ok(None);
        };
        for (graph, position) in bucket {
            if find(molecule, graph, 0)?.is_some() {
                return Ok(Some(*position));
            }
        }
        Ok(None)
    }
}

/// Whether two prepared molecules are the same graph, stereo included.
#[cfg(test)]
pub(crate) fn same_graph(a: &Molecule, b: &Molecule) -> Result<bool> {
    Ok(find(a, b, 0)?.is_some())
}

/// The input-to-entry atom correspondence for a lookup hit, with upstream's
/// relaxations: exact, then without charges and bond orders, then also
/// without stereo.
pub(crate) fn mapping(query: &Molecule, entry: &Molecule) -> Result<Vec<AtomId>> {
    if query.atom_count() != entry.atom_count() || query.bond_count() != entry.bond_count() {
        return Err(error("lookup hit has incompatible graph size"));
    }
    for mode in 0..3 {
        if let Some(map) = find(query, entry, mode)? {
            return Ok(map);
        }
    }
    Err(error(
        "lookup hit could not be mapped to the input molecule",
    ))
}

/// The first correspondence for one relaxation mode, in query atom order.
fn find(query: &Molecule, entry: &Molecule, mode: u8) -> Result<Option<Vec<AtomId>>> {
    if query.atom_count() != entry.atom_count() || query.bond_count() != entry.bond_count() {
        return Ok(None);
    }
    let ids: Vec<_> = query.atom_ids().collect();
    let targets: Vec<_> = entry.atom_ids().collect();
    let mut candidates = Vec::new();
    for &a in &ids {
        let qa = query.atom(a).map_err(Error::chemistry)?;
        let degree = query.neighbors(a).map_err(Error::chemistry)?.count();
        let mut row = Vec::new();
        for &b in &targets {
            let eb = entry.atom(b).map_err(Error::chemistry)?;
            if qa.element == eb.element
                && qa.isotope == eb.isotope
                && degree == entry.neighbors(b).map_err(Error::chemistry)?.count()
                && (mode > 0 || qa.formal_charge == eb.formal_charge)
                && query.atom_is_aromatic(a).map_err(Error::chemistry)?
                    == entry.atom_is_aromatic(b).map_err(Error::chemistry)?
            {
                row.push(b);
            }
        }
        candidates.push(row);
    }
    let mut order: Vec<_> = (0..ids.len()).collect();
    order.sort_by_key(|&i| {
        (
            candidates[i].len(),
            std::cmp::Reverse(query.neighbors(ids[i]).unwrap().count()),
            i,
        )
    });
    let mut assigned = BTreeMap::new();
    let mut states = 0;
    let search = Search {
        query,
        entry,
        ids: &ids,
        order: &order,
        candidates: &candidates,
        mode,
        state_limit: 1_000_000_usize.max(ids.len().saturating_mul(1_000)),
    };
    Ok(search
        .visit(0, &mut assigned, &mut states)?
        .then(|| ids.iter().map(|id| assigned[id]).collect()))
}
struct Search<'a> {
    query: &'a Molecule,
    entry: &'a Molecule,
    ids: &'a [AtomId],
    order: &'a [usize],
    candidates: &'a [Vec<AtomId>],
    mode: u8,
    state_limit: usize,
}
impl Search<'_> {
    fn visit(
        &self,
        depth: usize,
        map: &mut BTreeMap<AtomId, AtomId>,
        states: &mut usize,
    ) -> Result<bool> {
        *states += 1;
        // Bounds pathological symmetry, not molecule size.
        if *states > self.state_limit {
            return Err(Error::new(
                ErrorKind::ResourceLimit,
                "lookup mapping search limit exceeded",
            ));
        }
        if depth == self.ids.len() {
            return Ok(self.mode == 2 || stereo_matches(self.query, self.entry, map));
        }
        let i = self.order[depth];
        let a = self.ids[i];
        for &b in &self.candidates[i] {
            if map.values().any(|v| *v == b) {
                continue;
            }
            let mut compatible = true;
            for (&x, &y) in map.iter() {
                let qb = self.query.bond_between(a, x).map_err(Error::chemistry)?;
                let eb = self.entry.bond_between(b, y).map_err(Error::chemistry)?;
                match (qb, eb) {
                    (None, None) => {}
                    (Some(q), Some(e)) => {
                        if self.query.bond_is_aromatic(q).map_err(Error::chemistry)?
                            != self.entry.bond_is_aromatic(e).map_err(Error::chemistry)?
                        {
                            compatible = false;
                            break;
                        }
                        if self.mode == 0
                            && self.query.bond(q).map_err(Error::chemistry)?.order
                                != self.entry.bond(e).map_err(Error::chemistry)?.order
                        {
                            compatible = false;
                            break;
                        }
                    }
                    _ => {
                        compatible = false;
                        break;
                    }
                }
            }
            if compatible {
                map.insert(a, b);
                if self.visit(depth + 1, map, states)? {
                    return Ok(true);
                }
                map.remove(&a);
            }
        }
        Ok(false)
    }
}
fn stereo_matches(query: &Molecule, entry: &Molecule, map: &BTreeMap<AtomId, AtomId>) -> bool {
    if query.stereo_elements().count() != entry.stereo_elements().count() {
        return false;
    }
    let carrier = |c: StereoCarrier| match c {
        StereoCarrier::Atom(a) => StereoCarrier::Atom(map[&a]),
        x => x,
    };
    for (_, element) in query.stereo_elements() {
        let matched =
            entry
                .stereo_elements()
                .any(|(_, other)| match (&element.kind, &other.kind) {
                    (StereoElementKind::Tetrahedral(a), StereoElementKind::Tetrahedral(b))
                        if map[&a.center] == b.center =>
                    {
                        let mut permutation = Vec::new();
                        for c in &a.carriers {
                            match b.carriers.iter().position(|x| *x == carrier(*c)) {
                                Some(i) => permutation.push(i),
                                None => return false,
                            }
                        }
                        let odd = (0..permutation.len())
                            .flat_map(|i| (i + 1..permutation.len()).map(move |j| (i, j)))
                            .filter(|&(i, j)| permutation[i] > permutation[j])
                            .count()
                            % 2
                            == 1;
                        match (a.orientation, b.orientation) {
                            (Some(x), Some(y)) => (x == y) != odd,
                            (None, None) => true,
                            _ => false,
                        }
                    }
                    (StereoElementKind::DoubleBond(a), StereoElementKind::DoubleBond(b)) => {
                        let (left, right) = if map[&a.left] == b.left && map[&a.right] == b.right {
                            (b.left_carrier, b.right_carrier)
                        } else if map[&a.left] == b.right && map[&a.right] == b.left {
                            (b.right_carrier, b.left_carrier)
                        } else {
                            return false;
                        };
                        let odd =
                            (carrier(a.left_carrier) != left) ^ (carrier(a.right_carrier) != right);
                        match (a.orientation, b.orientation) {
                            (Some(x), Some(y)) => (x == y) != odd,
                            (None, None) => true,
                            _ => false,
                        }
                    }
                    _ => false,
                });
        if !matched {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prepared(smiles: &str) -> Molecule {
        crate::explicit(&kekule::smiles::to_molecules(smiles).unwrap().remove(0)).unwrap()
    }

    #[test]
    fn selection_requires_the_exact_entry_molecule() {
        let index = LookupIndex::new([
            (prepared("[H:4][N:2]=[N:1][F:3]"), 7),
            (prepared("[H]O[H]"), 3),
        ])
        .unwrap();
        let select = |s: &str| index.select(&prepared(s)).unwrap();
        assert_eq!(select("FN=N[H]"), Some(7));
        assert_eq!(select("[H]O[H]"), Some(3));
        // Another charge placement, stereo, or constitution is not the entry.
        assert_eq!(select("[H][N-]N=[F+]"), None);
        assert_eq!(select("[H]/N=N/F"), None);
        assert_eq!(select("[H]N([H])F"), None);
    }

    #[test]
    fn mapping_relaxes_charges_and_orders_then_stereo() {
        let entry = prepared("[H:4][N:2]=[N:1][F:3]");
        let query = prepared("[F+:3]=[N:1][N-:2][H:4]");
        let map = mapping(&query, &entry).unwrap();
        let elements = |m: &Molecule, ids: Vec<AtomId>| {
            ids.into_iter()
                .map(|a| m.atom(a).unwrap().element)
                .collect::<Vec<_>>()
        };
        assert_eq!(
            elements(&query, query.atom_ids().collect()),
            elements(&entry, map)
        );
        let other = prepared("[H]N([H])F");
        assert_eq!(
            mapping(&other, &entry).unwrap_err().kind(),
            ErrorKind::Identity
        );
    }
}
