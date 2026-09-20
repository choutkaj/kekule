use crate::core::{AtomId, Molecule, StereoCarrier, StereoElementKind};
use crate::query::{QueryGraph, QueryStereoConstraint};

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
