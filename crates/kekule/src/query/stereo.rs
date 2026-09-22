use crate::core::{DoubleBondOrientation, TetrahedralOrientation};

use super::{QueryAtomId, QueryBond, QueryBondId, QueryGraphError};

/// A correlated configuration predicate over tetrahedral query centers.
/// Members refer to checked local stereo constraints or Boolean carrier frames.
/// Distinct groups are independent; group order has no matching significance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryStereoGroup {
    pub kind: crate::core::StereoGroupKind,
    pub members: Vec<QueryAtomId>,
}

/// A stereochemical predicate evaluated under a complete query-to-target mapping.
///
/// These constraints describe local configuration, not CIP labels or enhanced
/// stereo groups. They require a specified target assertion. Target coordinates
/// are not interpreted and stereo is not perceived during matching.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryStereoConstraint {
    /// Ordered explicit query neighbors, followed by any omitted target carriers.
    /// With fewer than three neighbors, either handedness can complete the query;
    /// the constraint still requires specified tetrahedral stereo on the target.
    /// A query atom's hydrogen predicate is independent of these carrier references.
    Tetrahedral {
        center: QueryAtomId,
        carriers: Vec<QueryAtomId>,
        orientation: TetrahedralOrientation,
    },
    /// Configuration of the selected carriers at `bond.a()` and `bond.b()`.
    /// The query must contain the two carrier-to-endpoint edges.
    DoubleBond {
        bond: QueryBondId,
        left_carrier: QueryAtomId,
        right_carrier: QueryAtomId,
        orientation: DoubleBondOrientation,
    },
}

impl QueryStereoConstraint {
    pub(super) fn focus(&self) -> (u8, u32) {
        match self {
            Self::Tetrahedral { center, .. } => (0, center.raw()),
            Self::DoubleBond { bond, .. } => (1, bond.raw()),
        }
    }

    pub(super) fn validate(
        &self,
        atoms: usize,
        bonds: &[QueryBond],
        adjacency: &[Vec<QueryBondId>],
    ) -> Result<(), QueryGraphError> {
        let atom = |id: QueryAtomId| {
            if id.index() < atoms {
                Ok(())
            } else {
                Err(QueryGraphError::InvalidAtomId(id))
            }
        };
        let adjacent = |a: QueryAtomId, b: QueryAtomId| {
            adjacency[a.index()]
                .iter()
                .any(|id| bonds[id.index()].other_atom(a) == b)
        };
        match self {
            Self::Tetrahedral {
                center, carriers, ..
            } => {
                atom(*center)?;
                if carriers.len() > 4 || carriers.len() != adjacency[center.index()].len() {
                    return Err(QueryGraphError::InvalidStereo(
                        "tetrahedral carriers must cover all query neighbors, at most four",
                    ));
                }
                for (index, carrier) in carriers.iter().enumerate() {
                    atom(*carrier)?;
                    if !adjacent(*center, *carrier) || carriers[..index].contains(carrier) {
                        return Err(QueryGraphError::InvalidStereo(
                            "tetrahedral carriers must be distinct neighbors",
                        ));
                    }
                }
            }
            Self::DoubleBond {
                bond,
                left_carrier,
                right_carrier,
                ..
            } => {
                let bond = bonds
                    .get(bond.index())
                    .ok_or(QueryGraphError::InvalidBondId(*bond))?;
                atom(*left_carrier)?;
                atom(*right_carrier)?;
                if [bond.a(), bond.b()].contains(left_carrier)
                    || [bond.a(), bond.b()].contains(right_carrier)
                    || left_carrier == right_carrier
                    || !adjacent(bond.a(), *left_carrier)
                    || !adjacent(bond.b(), *right_carrier)
                {
                    return Err(QueryGraphError::InvalidStereo("double-bond carriers must be distinct substituents at the corresponding endpoints"));
                }
            }
        }
        Ok(())
    }

    pub(super) fn canonicalize(&mut self) {
        if let Self::Tetrahedral {
            carriers,
            orientation,
            ..
        } = self
        {
            let odd = carriers
                .iter()
                .enumerate()
                .map(|(i, atom)| {
                    carriers[i + 1..]
                        .iter()
                        .filter(|other| *other < atom)
                        .count()
                })
                .sum::<usize>()
                % 2
                != 0;
            carriers.sort_unstable();
            if odd {
                *orientation = orientation.inverted();
            }
        }
    }
}
