use super::one_based_serial;

#[test]
fn one_based_serials_widen_before_incrementing() {
    assert_eq!(one_based_serial(0), 1);
    assert_eq!(one_based_serial(u32::MAX), u64::from(u32::MAX) + 1);
}
