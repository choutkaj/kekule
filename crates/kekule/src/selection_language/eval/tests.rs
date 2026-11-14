use super::{compare_integer_real, Bits};
use crate::topology::{AtomSelection, Topology, TopologyAtomIndex};
use std::cmp::Ordering::{Equal, Greater, Less};
use std::sync::Arc;

#[test]
fn host_selections_become_their_bits() {
    let molecules = crate::smiles::to_molecules(&["C"; 300].join(".")).unwrap();
    let topology = Arc::new(Topology::from_molecules(molecules).unwrap());
    let runs: [&[std::ops::Range<usize>]; 6] = [
        &[],
        &[0..150, 150..300],
        &[0..1, 2..3, 299..300],
        &[0..64, 64..65, 127..129],
        &[5..6, 63..200, 201..202, 230..299],
        &[1..64, 64..128, 128..192, 256..257],
    ];
    for ranges in runs {
        let indices = ranges.iter().cloned().flatten().collect::<Vec<_>>();
        let selection = AtomSelection::from_indices(
            Arc::clone(&topology),
            indices
                .iter()
                .map(|&index| TopologyAtomIndex::new(index as u32)),
        )
        .unwrap();
        let bits = Bits::of_selection(&selection);
        assert_eq!(bits.ones().collect::<Vec<_>>(), indices, "{ranges:?}");
    }
}

#[test]
fn integers_order_exactly_against_reals() {
    assert_eq!(compare_integer_real(3, 3.0), Equal);
    assert_eq!(compare_integer_real(3, 3.5), Less);
    assert_eq!(compare_integer_real(4, 3.5), Greater);
    assert_eq!(compare_integer_real(-4, -3.5), Less);
    assert_eq!(compare_integer_real(-3, -3.5), Greater);
    assert_eq!(compare_integer_real(i64::MAX, 9.3e18), Less);
    assert_eq!(compare_integer_real(i64::MIN, -9.3e18), Greater);
    assert_eq!(
        compare_integer_real(i64::MIN, -9_223_372_036_854_775_808.0),
        Equal
    );
    // 2^53 + 1 is not a float; 2^53 is.
    assert_eq!(
        compare_integer_real(9_007_199_254_740_993, 9_007_199_254_740_992.0),
        Greater
    );
}
