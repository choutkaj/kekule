use crate::diagnostics::{LimitExceeded, Resource};
use std::sync::Arc;

use super::*;
use crate::smiles;
use crate::topology::Topology;

fn topology(smiles: &str) -> Arc<Topology> {
    smiles::to_topology(smiles).unwrap()
}

fn symmetric(topology: &Arc<Topology>, indices: &[usize]) -> SymmetricSelection {
    let selection = AtomSelection::from_indices(
        Arc::clone(topology),
        indices
            .iter()
            .map(|index| crate::topology::TopologyAtomIndex::new(u32::try_from(*index).unwrap())),
    )
    .unwrap();
    SymmetricSelection::new(&selection).unwrap()
}

fn all(topology: &Arc<Topology>) -> SymmetricSelection {
    SymmetricSelection::new(&AtomSelection::all(Arc::clone(topology))).unwrap()
}

fn permutations(selection: &SymmetricSelection) -> Vec<Vec<usize>> {
    selection.permutations().map(<[usize]>::to_vec).collect()
}

#[test]
fn permutation_counts_match_graph_automorphism_group_orders() {
    // Orders of the graph automorphism groups, derived by hand: the dihedral
    // group of a six-ring, one ring flip, the four methyls of neopentane, and
    // both ends and hydrogens of ethane (2 * 3! * 3!).
    let cases = [
        ("c1ccccc1", 12),
        ("Cc1ccccc1", 2),
        ("CC(C)(C)C", 24),
        ("[H]C([H])([H])[H]", 24),
        ("[H]C([H])([H])C([H])([H])[H]", 72),
        // Different hydrogen counts break the skeleton's symmetry.
        ("C=CC", 1),
        ("CC(=O)O", 1),
    ];
    for (smiles, expected) in cases {
        assert_eq!(
            all(&topology(smiles)).permutation_count(),
            expected,
            "{smiles}"
        );
    }
}

#[test]
fn bond_orders_and_charges_do_not_break_symmetry() {
    // Atoms: C0, C1, O2 (double bond), O3 (negative); the oxygens exchange.
    let acetate = topology("CC(=O)[O-]");
    assert_eq!(
        permutations(&all(&acetate)),
        [vec![0, 1, 2, 3], vec![0, 1, 3, 2]]
    );
    // A Kekulé ring is as symmetric as an aromatic one.
    assert_eq!(all(&topology("C1=CC=CC=C1")).permutation_count(), 12);
}

#[test]
fn permutations_restrict_symmetries_to_the_selection() {
    let methane = topology("[H]C([H])([H])[H]");
    // Hydrogen permutations that fix the carbon leave it paired with itself.
    assert_eq!(permutations(&symmetric(&methane, &[1])), [vec![0]]);
    // Fixing one hydrogen setwise still lets the other three exchange, but
    // they are not selected.
    assert_eq!(permutations(&symmetric(&methane, &[0, 1])), [vec![0, 1]]);
    // Selected hydrogens exchange among themselves only.
    assert_eq!(
        permutations(&symmetric(&methane, &[0, 2])),
        [vec![0, 1], vec![1, 0]]
    );

    // Selecting only the acetate oxygens keeps their exchange.
    let acetate = topology("CC(=O)[O-]");
    assert_eq!(
        permutations(&symmetric(&acetate, &[2, 3])),
        [vec![0, 1], vec![1, 0]]
    );
}

#[test]
fn identical_molecule_instances_are_never_exchanged() {
    // Each acetate flips independently; the two ions never swap.
    let pair = topology("CC(=O)[O-].CC(=O)[O-]");
    assert_eq!(pair.instance_count(), 2);
    assert_eq!(
        permutations(&all(&pair)),
        [
            vec![0, 1, 2, 3, 4, 5, 6, 7],
            vec![0, 1, 2, 3, 4, 5, 7, 6],
            vec![0, 1, 3, 2, 4, 5, 6, 7],
            vec![0, 1, 3, 2, 4, 5, 7, 6],
        ]
    );
    // Unselected instances are not searched.
    assert_eq!(
        permutations(&symmetric(&pair, &[2, 3])),
        [vec![0, 1], vec![1, 0]]
    );
}

#[test]
fn searches_beyond_their_limits_fail_instead_of_truncating() {
    let neopentane = topology("CC(C)(C)C");
    let selection = AtomSelection::all(Arc::clone(&neopentane));
    let limits = SymmetryLimits::default().with_max_permutations(23);
    assert_eq!(
        SymmetricSelection::with_limits(&selection, limits).unwrap_err(),
        SymmetryError::LimitExceeded(LimitExceeded::new(Resource::Permutations, 23, 24))
    );
    let limits = SymmetryLimits::default().with_max_search_steps(10);
    assert_eq!(
        SymmetricSelection::with_limits(&selection, limits).unwrap_err(),
        SymmetryError::LimitExceeded(LimitExceeded::new(Resource::SearchStates, 10, 11))
    );
    let exact = SymmetryLimits::default().with_max_permutations(24);
    assert_eq!(
        SymmetricSelection::with_limits(&selection, exact)
            .unwrap()
            .permutation_count(),
        24
    );

    let empty = SymmetricSelection::new(&AtomSelection::empty(Arc::clone(&neopentane))).unwrap();
    assert!(empty.is_empty());
    assert_eq!(permutations(&empty), [Vec::<usize>::new()]);
}
