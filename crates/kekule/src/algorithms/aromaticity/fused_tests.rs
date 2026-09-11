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
            let result =
                visit_connected_ring_subsets(&neighbors, &indexes, subset_size, &mut |subset| {
                    let mask = subset.iter().fold(0, |mask, ring| mask | (1 << ring));
                    assert!(
                        actual.insert(mask),
                        "duplicate subset in graph {graph_mask}"
                    );
                    ControlFlow::Continue(())
                });
            assert!(result.is_continue());
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
    let result = visit_connected_ring_subsets(&neighbors, &indexes, 6, &mut |_| {
        count += 1;
        ControlFlow::Continue(())
    });
    assert!(result.is_continue());
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
    let result = visit_connected_ring_subsets(&neighbors, &indexes, 6, &mut |_| {
        count += 1;
        ControlFlow::Break(())
    });
    assert!(result.is_break());
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
    assert!(rdkit_rings_are_fused(&left, &one_shared));
    assert!(!rdkit_rings_are_fused(&left, &two_shared));
    assert!(!rdkit_rings_are_fused(&left, &disjoint));
    assert!(!rdkit_rings_are_fused(&left, &too_large));
    assert!(!rdkit_rings_are_fused(&too_large, &left));
}
