use super::*;

#[test]
fn unchanged_run_flag_reuses_the_previous_nonzero_run_length() {
    // Five pairs. Only the first pair supplies a run-length update; later
    // zero flags mean "reuse", not "no run" (xdrfile decompression contract).
    // Pack this focused bitstream independently of the production encoder.
    let mut bits = Vec::new();
    let mut write = |value: u32, count: u32| {
        for bit in (0..count).rev() {
            bits.push((value >> bit) & 1);
        }
    };
    for pair in 0..5 {
        // Absolute mixed-radix [pair, 0, 0], radices [8, 8, 8].
        let packed = pair * 64;
        write(packed & 255, 8);
        write(packed >> 8, 2);
        write(u32::from(pair == 0), 1);
        if pair == 0 {
            write(4, 5);
        } // Run 3; small-index change 0.
          // Delta [5, 4, 4] relative to the absolute coordinate, centered at 4.
        write(356 & 255, 8);
        write(356 >> 8, 1);
    }
    let mut payload = vec![0_u8; bits.len().div_ceil(8)];
    for (index, bit) in bits.into_iter().enumerate() {
        payload[index / 8] |= (bit as u8) << (7 - index % 8);
    }
    let mut output = Vec::new();
    decode_compressed_payload(
        &payload,
        CompressedLayout {
            atom_count: 10,
            minimum: [0; 3],
            maximum: [7; 3],
            initial_small_index: XTC_FIRST_SMALL_INDEX,
            precision: 1.0,
        },
        &mut output,
        "repeated-run",
        0,
    )
    .unwrap();
    let expected = (0..5)
        .flat_map(|pair| [(pair + 1) as f32, 0.0, 0.0, pair as f32, 0.0, 0.0])
        .collect::<Vec<_>>();
    assert_eq!(output, expected);
}
