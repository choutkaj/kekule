use super::*;

fn neighbors_for_edges(size: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut neighbors = vec![Vec::new(); size];
    for &(left, right) in edges {
        neighbors[left].push(right);
        neighbors[right].push(left);
    }
    neighbors
}

#[test]
fn connected_subset_traversal_matches_exhaustive_small_graphs_without_duplicates() {
    // Exhaust every undirected graph on five vertices. The oracle enumerates
    // bit masks and checks reachability independently of the production search.
    let size = 5;
    let possible_edges = (0..size)
        .flat_map(|left| (left + 1..size).map(move |right| (left, right)))
        .collect::<Vec<_>>();
    let indexes = (0..size).collect::<Vec<_>>();
    for graph_mask in 0..1usize << possible_edges.len() {
        let edges = possible_edges
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(index, edge)| (graph_mask & (1 << index) != 0).then_some(edge))
            .collect::<Vec<_>>();
        let neighbors = neighbors_for_edges(size, &edges);
        for subset_size in 1..=size {
            let expected = (1..1usize << size)
                .filter(|mask| mask.count_ones() as usize == subset_size)
                .filter(|&mask| {
                    let mut reached = 1 << mask.trailing_zeros();
                    loop {
                        let previous = reached;
                        for &(left, right) in &edges {
                            let pair = (1 << left) | (1 << right);
                            if pair & mask == pair && pair & reached != 0 {
                                reached |= pair;
                            }
                        }
                        if previous == reached {
                            return reached == mask;
                        }
                    }
                })
                .collect::<BTreeSet<_>>();
            let mut actual = BTreeSet::new();
            let result = visit_connected_ring_subsets(
                &neighbors,
                &indexes,
                subset_size,
                &mut AromaticityWork::new(5_000_000),
                &mut |subset, _| {
                    let mask = subset.iter().fold(0, |mask, ring| mask | (1 << ring));
                    assert!(
                        actual.insert(mask),
                        "duplicate subset in graph {graph_mask}"
                    );
                    Ok(ControlFlow::Continue(()))
                },
            );
            assert!(result.unwrap().is_continue());
            assert_eq!(actual, expected, "graph {graph_mask}, size {subset_size}");
        }
    }
}

#[test]
fn connected_subset_search_scales_with_connected_candidates_in_a_long_chain() {
    // There are only 295 connected six-ring subsets, but C(300, 6) exceeds
    // 962 billion. Enumerating all combinations before testing fusion stalls.
    let size = 300;
    let edges = (1..size).map(|ring| (ring - 1, ring)).collect::<Vec<_>>();
    let neighbors = neighbors_for_edges(size, &edges);
    let indexes = (0..size).collect::<Vec<_>>();
    let mut count = 0;
    let result = visit_connected_ring_subsets(
        &neighbors,
        &indexes,
        6,
        &mut AromaticityWork::new(5_000_000),
        &mut |_, _| {
            count += 1;
            Ok(ControlFlow::Continue(()))
        },
    );
    assert!(result.unwrap().is_continue());
    assert_eq!(count, 295);
}

#[test]
fn connected_subset_search_stops_immediately_when_the_caller_is_done() {
    let size = 100;
    let edges = (0..size)
        .flat_map(|left| (left + 1..size).map(move |right| (left, right)))
        .collect::<Vec<_>>();
    let neighbors = neighbors_for_edges(size, &edges);
    let indexes = (0..size).collect::<Vec<_>>();
    let mut count = 0;
    let result = visit_connected_ring_subsets(
        &neighbors,
        &indexes,
        6,
        &mut AromaticityWork::new(5_000_000),
        &mut |_, _| {
            count += 1;
            Ok(ControlFlow::Break(()))
        },
    );
    assert!(result.unwrap().is_break());
    assert_eq!(count, 1);
}

#[test]
fn fused_neighbors_require_exactly_one_shared_bond_and_at_most_24_ring_bonds() {
    let ring = |bonds: Vec<u32>| Ring {
        atoms: Vec::new(),
        bonds: bonds.into_iter().map(BondId::new).collect(),
    };
    let left = ring((0..24).collect());
    let one_shared = ring(std::iter::once(0).chain(24..47).collect());
    let two_shared = ring(vec![0, 1, 24]);
    let disjoint = ring(vec![24, 25, 26]);
    let too_large = ring(std::iter::once(0).chain(24..48).collect());
    for (other, fused) in [
        (one_shared, true),
        (two_shared, false),
        (disjoint, false),
        (too_large, false),
    ] {
        let neighbors = rdkit_fused_ring_neighbors(
            &[left.clone(), other],
            &[0, 1],
            &mut AromaticityWork::new(5_000_000),
        )
        .unwrap();
        assert_eq!(
            neighbors,
            if fused {
                vec![vec![1], vec![0]]
            } else {
                vec![vec![], vec![]]
            }
        );
    }
}

#[test]
fn dense_subset_search_returns_a_work_error_instead_of_partial_success() {
    let size = 30;
    let edges = (0..size)
        .flat_map(|left| (left + 1..size).map(move |right| (left, right)))
        .collect::<Vec<_>>();
    let neighbors = neighbors_for_edges(size, &edges);
    let result = visit_connected_ring_subsets(
        &neighbors,
        &(0..size).collect::<Vec<_>>(),
        6,
        &mut AromaticityWork::new(10_000),
        &mut |_, _| Ok(ControlFlow::Continue(())),
    );
    assert_eq!(
        result,
        Err(AromaticityError::ResourceLimit { limit: 10_000 })
    );
}

#[test]
fn fusion_discovery_scales_with_independent_rings_and_bounds_dense_families() {
    let independent = (0..10_000)
        .map(|index| Ring {
            atoms: Vec::new(),
            bonds: (index * 6..index * 6 + 6).map(BondId::new).collect(),
        })
        .collect::<Vec<_>>();
    let candidates = (0..independent.len()).collect::<Vec<_>>();
    let neighbors = rdkit_fused_ring_neighbors(
        &independent,
        &candidates,
        &mut AromaticityWork::new(120_000),
    )
    .unwrap();
    assert!(neighbors.iter().all(Vec::is_empty));

    let dense = vec![independent[0].clone(); 100];
    assert_eq!(
        rdkit_fused_ring_neighbors(
            &dense,
            &(0..100).collect::<Vec<_>>(),
            &mut AromaticityWork::new(1_000)
        ),
        Err(AromaticityError::ResourceLimit { limit: 1_000 })
    );
}
