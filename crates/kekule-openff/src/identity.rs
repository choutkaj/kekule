//! Narrow InChI boundary and bounded whole-molecule lookup correspondence.
use crate::{Error, ErrorKind, Result};

use kekule::core::*;
use std::collections::BTreeMap;

fn error(detail: impl std::fmt::Display) -> Error {
    Error::new(ErrorKind::Identity, detail)
}

/// The most atoms the InChI library accepts without its `LargeMolecules`
/// option, which the `inchi` crate does not expose. Upstream can identify larger molecules, so bundles with larger
/// lookup entries are rejected rather than silently missing their hits.
pub(crate) const MAX_INCHI_ATOMS: usize = 1023;

pub(crate) fn fixed_h_inchi(molecule: &Molecule) -> Result<String> {
    if molecule
        .stereo_groups()
        .any(|(_, g)| g.kind != StereoGroupKind::Absolute)
    {
        return Err(error(
            "enhanced stereo groups are not supported by the InChI adapter",
        ));
    }
    let mut native = inchi::Molecule::new();
    let mut indices = BTreeMap::new();
    for (id, atom) in molecule.atoms() {
        // Match the RDKit InChI adapter's H policy, including its element
        // distinction. This affects InChI's tautomer normalization even when
        // the input already contains every hydrogen as a graph vertex.
        let hydrogens = if matches!(atom.element.atomic_number(), 6 | 7 | 8 | 9 | 17 | 35 | 53) {
            inchi::ImplicitH::Auto
        } else {
            inchi::ImplicitH::Exactly(0)
        };
        let mut a = inchi::Atom::new(atom.element.symbol())
            .charge(atom.formal_charge)
            .implicit_hydrogens(hydrogens);
        if let Some(mass) = atom.isotope {
            a = a.isotope(mass);
        }
        indices.insert(id, native.add_atom(a));
    }
    for (_, bond) in molecule.bonds() {
        let order = match bond.order {
            BondOrder::Single => inchi::BondOrder::Single,
            BondOrder::Double => inchi::BondOrder::Double,
            BondOrder::Triple => inchi::BondOrder::Triple,
            _ => return Err(error("InChI requires localized single/double/triple bonds")),
        };
        native
            .add_bond(indices[&bond.a()], indices[&bond.b()], order)
            .map_err(Error::chemistry)?;
    }
    for (_, element) in molecule.stereo_elements() {
        match &element.kind {
            StereoElementKind::Tetrahedral(s) => {
                let center = indices[&s.center];
                let mut neighbors = s
                    .carriers
                    .iter()
                    .map(|c| match c {
                        StereoCarrier::Atom(a) => Ok(indices[a]),
                        StereoCarrier::ImplicitLonePair => Ok(center),
                        _ => Err(error("InChI stereo requires explicit hydrogens")),
                    })
                    .collect::<Result<Vec<_>>>()?;
                // InChI even parity is clockwise viewed from carrier 0;
                // Kekule defines orientation by the signed carrier volume.
                let mut even = matches!(
                    s.orientation,
                    Some(TetrahedralOrientation::CounterClockwise)
                );
                if let Some(i) = neighbors.iter().position(|a| *a == center) {
                    if i != 0 {
                        neighbors.swap(0, i);
                        even = !even;
                    }
                }
                let parity = if s.orientation.is_none() {
                    inchi::Parity::Unknown
                } else if even {
                    inchi::Parity::Even
                } else {
                    inchi::Parity::Odd
                };
                native.add_stereo(inchi::Stereo::Tetrahedral {
                    center,
                    neighbors: neighbors
                        .try_into()
                        .map_err(|_| error("invalid tetrahedral carriers"))?,
                    parity,
                });
            }
            StereoElementKind::DoubleBond(s) => {
                let carrier = |c| match c {
                    StereoCarrier::Atom(a) => Ok(indices[&a]),
                    _ => Err(error("double-bond stereo requires explicit atom carriers")),
                };
                native.add_stereo(inchi::Stereo::DoubleBond {
                    ends: [
                        carrier(s.left_carrier)?,
                        indices[&s.left],
                        indices[&s.right],
                        carrier(s.right_carrier)?,
                    ],
                    parity: match s.orientation {
                        Some(DoubleBondOrientation::Together) => inchi::Parity::Odd,
                        Some(DoubleBondOrientation::Opposite) => inchi::Parity::Even,
                        None => inchi::Parity::Unknown,
                    },
                });
            }
            StereoElementKind::Axis(_) => {
                return Err(error("axial stereo is not supported by the InChI adapter"))
            }
        }
    }
    native
        .to_inchi(inchi::Options::new().fixed_h(true))
        .map(|x| x.into_inchi())
        .map_err(|e| Error::wrap(ErrorKind::Identity, e))
}

pub(crate) fn mapping(query: &Molecule, entry: &Molecule) -> Result<Vec<AtomId>> {
    if query.atom_count() != entry.atom_count() || query.bond_count() != entry.bond_count() {
        return Err(error("InChI lookup hit has incompatible graph size"));
    }
    let ids: Vec<_> = query.atom_ids().collect();
    let targets: Vec<_> = entry.atom_ids().collect();
    for mode in 0..3 {
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
        if search.visit(0, &mut assigned, &mut states)? {
            return Ok(ids.iter().map(|id| assigned[id]).collect());
        }
    }
    Err(error(
        "InChI lookup hit could not be mapped to the input molecule",
    ))
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
    #[test]
    fn adapter_atom_limit_matches_the_inchi_library() {
        let inchi = |smiles: &str| {
            let mut m = kekule::smiles::to_molecules(smiles).unwrap().remove(0);
            m.perceive().unwrap();
            m.add_hydrogens().unwrap();
            fixed_h_inchi(&m)
        };
        // C340H682O has 1,023 atoms; C340H683N has 1,024.
        inchi(&format!("{}O", "C".repeat(340))).unwrap();
        assert_eq!(
            inchi(&format!("{}N", "C".repeat(340))).unwrap_err().kind(),
            ErrorKind::Identity
        );
        assert_eq!(MAX_INCHI_ATOMS, 1023);
    }
    #[test]
    fn identifiers_match_all_pinned_reference_molecules() {
        for r in crate::reference_records() {
            let m = kekule::smiles::to_molecules(r["mapped_smiles"].as_str().unwrap())
                .unwrap()
                .remove(0);
            let m = crate::explicit(&m).unwrap();
            assert_eq!(
                fixed_h_inchi(&m).unwrap(),
                r["fixed_h_inchi"].as_str().unwrap(),
                "{}",
                r["input"]["id"]
            );
        }
    }
    #[test]
    fn alanine_carrier_parity_matches_reference_identifier() {
        // PubChem CID 5950, also present in the independently pinned audit corpus.
        let m=kekule::smiles::to_molecules("[C:1]([C@@:2]([C:3](=[O:4])[O:5][H:11])([N:6]([H:12])[H:13])[H:10])([H:7])([H:8])[H:9]").unwrap().remove(0);
        let id = fixed_h_inchi(&crate::explicit(&m).unwrap()).unwrap();
        assert_eq!(
            id,
            "InChI=1/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1/f/h5H"
        );
    }
}
