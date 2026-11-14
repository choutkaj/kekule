use super::*;

#[test]
fn ranges_cover_every_index_once_in_order() {
    for threads in 1..=9 {
        let executor = ScopedThreads::new(NonZeroUsize::new(threads).unwrap());
        for len in [0, 1, 7, 64, 1000] {
            let ranges = map_ranges(&executor, len, 4, |range| range);
            assert!(ranges.len() <= threads);
            assert!(ranges
                .iter()
                .all(|range| range.len() >= 4 || ranges.len() == 1));
            let indices = ranges.into_iter().flatten().collect::<Vec<_>>();
            assert_eq!(indices, (0..len).collect::<Vec<_>>());
        }
    }
}

#[test]
fn weighted_ranges_cover_every_index_once_and_share_the_weight() {
    let skewed = (0..60)
        .map(|index| if index < 10 { 50 } else { 1 })
        .collect::<Vec<_>>();
    for threads in 1..=9 {
        let executor = ScopedThreads::new(NonZeroUsize::new(threads).unwrap());
        for weights in [vec![], vec![0; 5], vec![7], vec![1; 100], skewed.clone()] {
            let ranges = map_weighted_ranges(&executor, &weights, |range| range);
            assert!(ranges.len() <= threads.max(1));
            let indices = ranges.iter().cloned().flatten().collect::<Vec<_>>();
            assert_eq!(indices, (0..weights.len()).collect::<Vec<_>>());
        }
        let ranges = map_weighted_ranges(&executor, &skewed, |range| range);
        let total = skewed.iter().sum::<usize>();
        let heaviest = skewed.iter().max().unwrap();
        for range in ranges {
            let weight = skewed[range].iter().sum::<usize>();
            assert!(
                weight <= total / threads + heaviest,
                "{threads} threads: {weight}"
            );
        }
    }
}
