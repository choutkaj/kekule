use crate::core::{AtomId, Molecule, StereoCarrier, StereoElementKind};
use crate::query::{QueryGraph, QueryStereoConstraint};

/// Test inclusion of the query's correlated configuration choices in the target's.
/// Each target group supplies one shared inversion variable. Distinct query
/// variables cannot satisfy a target correlation for all query configurations.
#[derive(Default)]
pub(super) struct GroupMatches {
    target_variables:
        std::collections::BTreeMap<(usize, crate::core::StereoGroupId), (Option<usize>, bool)>,
}

impl GroupMatches {
    pub fn add(
        &mut self,
        target: &Molecule,
        center: AtomId,
        occurrence: usize,
        query: &QueryGraph,
        query_center: crate::query::QueryAtomId,
        inverted: Option<bool>,
    ) -> bool {
        use crate::core::StereoGroupKind as K;
        let Some(element) = target.stereo_elements().find_map(|(_, e)| match &e.kind {
            StereoElementKind::Tetrahedral(s) if s.center == center && s.orientation.is_some() => {
                Some(e)
            }
            _ => None,
        }) else {
            return false;
        };
        let qgroup = query
            .stereo_groups()
            .iter()
            .enumerate()
            .find(|(_, g)| g.members.contains(&query_center));
        let qkind = qgroup.map_or(K::Absolute, |(_, g)| g.kind);
        let tgroup = element
            .group
            .map(|id| (id, target.stereo_group(id).expect("validated stereo group")));
        let tkind = tgroup.map_or(K::Absolute, |(_, g)| g.kind);
        let compatible = match qkind {
            K::Absolute => true,
            K::Or => matches!(tkind, K::Or | K::And),
            K::And => tkind == K::And,
            K::Relative => tkind == K::Relative,
            K::Racemic => tkind == K::Racemic,
        };
        if !compatible {
            return false;
        }
        if tkind == K::Absolute {
            return qkind == K::Absolute && inverted != Some(true);
        }
        let Some(inverted) = inverted else {
            return true;
        };
        let qvariable = qgroup.and_then(|(i, g)| (g.kind != K::Absolute).then_some(i));
        let key = (occurrence, tgroup.unwrap().0);
        let value = (qvariable, inverted);
        match self.target_variables.get(&key) {
            Some(previous) => *previous == value,
            None => {
                self.target_variables.insert(key, value);
                true
            }
        }
    }
}

pub(super) fn matches_constraint(
    target: &Molecule,
    query: &QueryGraph,
    mapping: &[AtomId],
    constraint: &QueryStereoConstraint,
) -> bool {
    match constraint {
        QueryStereoConstraint::Tetrahedral {
            center,
            carriers,
            orientation,
        } => {
            let center = mapping[center.index()];
            let stereo = target
                .stereo_elements()
                .find_map(|(_, element)| match &element.kind {
                    StereoElementKind::Tetrahedral(stereo) if stereo.center == center => {
                        Some(stereo)
                    }
                    _ => None,
                });
            let Some(stereo) = stereo else {
                return false;
            };
            let Some(actual) = stereo.orientation else {
                return false;
            };
            let mut permutation = Vec::with_capacity(4);
            for carrier in carriers {
                let Some(index) = stereo.carriers.iter().position(|target_carrier| {
                    *target_carrier == StereoCarrier::Atom(mapping[carrier.index()])
                }) else {
                    return false;
                };
                permutation.push(index);
            }
            // Two or more unconstrained carriers can be exchanged to satisfy
            // either orientation; one omitted carrier has a unique position.
            if permutation.len() < 3 {
                return true;
            }
            for index in 0..4 {
                if !permutation.contains(&index) {
                    permutation.push(index);
                }
            }
            let odd = permutation
                .iter()
                .enumerate()
                .map(|(i, index)| {
                    permutation[i + 1..]
                        .iter()
                        .filter(|other| *other < index)
                        .count()
                })
                .sum::<usize>()
                % 2
                != 0;
            actual
                == if odd {
                    orientation.inverted()
                } else {
                    *orientation
                }
        }
        QueryStereoConstraint::DoubleBond {
            bond,
            left_carrier,
            right_carrier,
            orientation,
        } => {
            let bond = query.bond(*bond).expect("validated query stereo bond");
            let a = mapping[bond.a().index()];
            let b = mapping[bond.b().index()];
            let Ok(Some(target_bond)) = target.bond_between(a, b) else {
                return false;
            };
            let stereo = target
                .stereo_elements()
                .find_map(|(_, element)| match &element.kind {
                    StereoElementKind::DoubleBond(stereo) if stereo.bond == target_bond => {
                        Some(stereo)
                    }
                    _ => None,
                });
            let Some(stereo) = stereo else {
                return false;
            };
            let Some(actual) = stereo.orientation else {
                return false;
            };
            let mut selected = [
                StereoCarrier::Atom(mapping[left_carrier.index()]),
                StereoCarrier::Atom(mapping[right_carrier.index()]),
            ];
            if stereo.left != a {
                selected.swap(0, 1);
            }
            let odd = (selected[0] != stereo.left_carrier) ^ (selected[1] != stereo.right_carrier);
            actual
                == if odd {
                    orientation.inverted()
                } else {
                    *orientation
                }
        }
    }
}
