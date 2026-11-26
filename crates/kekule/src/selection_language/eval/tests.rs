use super::{compare, compare_integer_real, value_matches, Bits, Value, ValueRanges};
use crate::selection_language::ast::{Literal, ValueMatch};
use crate::topology::{AtomSelection, Topology, TopologyAtomIndex};
use std::cmp::Ordering::{Equal, Greater, Less};
use std::sync::Arc;

/// A deterministic pseudo-random sequence (xorshift).
struct Numbers(u64);

impl Numbers {
    fn below(&mut self, bound: u64) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 % bound
    }

    /// A number in `-10..100`.
    fn number(&mut self) -> i64 {
        self.below(110) as i64 - 10
    }

    /// Exact values and ranges, some reversed, overlapping or out of range.
    fn list(&mut self, literal: impl Fn(&mut Self) -> Literal) -> Vec<ValueMatch> {
        (0..1 + self.below(12))
            .map(|_| {
                if self.below(3) == 0 {
                    ValueMatch::Exact(literal(self))
                } else {
                    ValueMatch::Range(literal(self), literal(self))
                }
            })
            .collect()
    }
}

#[test]
fn value_ranges_match_what_every_value_of_the_list_matches() {
    let mut numbers = Numbers(0x9e37_79b9_7f4a_7c15);
    let codes = [None, Some("A"), Some("B")];
    for _ in 0..500 {
        let integers = numbers.list(|numbers| Literal::Integer(numbers.number()));
        let ranges = ValueRanges::new(&integers).unwrap();
        for value in -20..120 {
            let compare = |literal: &Literal| compare(Value::Integer(value), literal, 1.0);
            let expected = integers
                .iter()
                .any(|candidate| value_matches(candidate, &compare));
            assert_eq!(ranges.contains(&compare), expected, "{integers:?} {value}");
        }
        let len = 64 + numbers.below(64) as usize;
        let indices = ranges.indices(len).ones().collect::<Vec<_>>();
        let expected = (0..len)
            .filter(|&index| {
                let compare =
                    |literal: &Literal| compare(Value::Integer(index as i64), literal, 1.0);
                integers
                    .iter()
                    .any(|candidate| value_matches(candidate, &compare))
            })
            .collect::<Vec<_>>();
        assert_eq!(indices, expected, "{integers:?}");

        let residues = numbers.list(|numbers| Literal::ResId {
            number: numbers.number() / 10,
            insertion_code: codes[numbers.below(3) as usize].map(str::to_owned),
        });
        let ranges = ValueRanges::new(&residues).unwrap();
        for number in -2..12 {
            for code in codes {
                let compare = |literal: &Literal| compare(Value::ResId(number, code), literal, 1.0);
                let expected = residues
                    .iter()
                    .any(|candidate| value_matches(candidate, &compare));
                assert_eq!(
                    ranges.contains(&compare),
                    expected,
                    "{residues:?} {number}{code:?}"
                );
            }
        }
    }
    let names =
        ["CA", "C", "N", "CA", "OXT"].map(|name| ValueMatch::Exact(Literal::Text(name.into())));
    let ranges = ValueRanges::new(&names).unwrap();
    for name in ["C", "CA", "CB", "N", "O", "OXT", ""] {
        let compare = |literal: &Literal| compare(Value::Text(name), literal, 1.0);
        let expected = names
            .iter()
            .any(|candidate| value_matches(candidate, &compare));
        assert_eq!(ranges.contains(&compare), expected, "{name}");
    }
    // Reals compare scaled, and properties by cell: neither takes ranges.
    let real = [ValueMatch::Exact(Literal::Real {
        value: 1.0,
        unit: None,
    })];
    assert!(ValueRanges::new(&real).is_none());
}

#[test]
fn bits_set_ranges_a_word_at_a_time() {
    for len in [0, 1, 63, 64, 65, 130, 200] {
        for start in (0..len + 3).step_by(7) {
            for end in (start..len + 70).step_by(5) {
                let mut bits = Bits::empty(len);
                bits.set_range(start..end);
                let expected = (start..end.min(len)).collect::<Vec<_>>();
                assert_eq!(
                    bits.ones().collect::<Vec<_>>(),
                    expected,
                    "{len} {start}..{end}"
                );
                assert_eq!(bits.is_empty(), expected.is_empty());
                assert_eq!(
                    bits.is_full(),
                    expected.len() == len,
                    "{len} {start}..{end}"
                );
            }
        }
        assert!(Bits::full(len).is_full());
        assert_eq!(Bits::empty(len).is_full(), len == 0);
    }
}

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
