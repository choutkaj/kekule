use super::*;

/// An index of `format` whose frames record a time and step exactly when
/// `stamped`.
fn index_of(
    format: TrajectoryFormat,
    starts: Vec<u64>,
    byte_count: u64,
    stamped: bool,
) -> TrajectoryFrameIndex {
    TrajectoryFrameIndex::new(
        format,
        TrajectoryFileFingerprint {
            byte_count,
            modified: Some(-1_500_000_000_123),
            content_hash: 0x0123_4567_89ab_cdef,
        },
        Arc::new(FrameOffsets {
            times: stamped.then(|| {
                starts
                    .iter()
                    .map(|&start| start as f64 * 0.25 - 1.0)
                    .collect()
            }),
            steps: stamped.then(|| starts.iter().map(|&start| start * 10).collect()),
            starts,
            end: byte_count,
        }),
        true,
    )
}

fn index(starts: Vec<u64>, byte_count: u64) -> TrajectoryFrameIndex {
    index_of(TrajectoryFormat::Trr, starts, byte_count, true)
}

fn reencode(bytes: &mut [u8]) {
    let body = bytes.len() - CHECKSUM_BYTES;
    let mut checksum = Fnv1a::new();
    checksum.write(&bytes[..body]);
    bytes[body..].copy_from_slice(&checksum.finish().to_le_bytes());
}

fn kind(bytes: &[u8], limits: &TrajectoryIoLimits) -> TrajectoryIoErrorKind {
    TrajectoryFrameIndex::from_bytes(bytes, limits)
        .unwrap_err()
        .kind()
}

#[test]
fn encoding_round_trips_exactly() {
    let original = index(vec![0, 100, 250], 400);
    let bytes = original.to_bytes();
    assert_eq!(bytes.len(), HEADER_BYTES + 3 * 3 * 8 + CHECKSUM_BYTES);
    let decoded = TrajectoryFrameIndex::from_bytes(&bytes, &TrajectoryIoLimits::default()).unwrap();
    assert_eq!(decoded, original);
    assert_eq!(decoded.to_bytes(), bytes);
    assert_eq!(decoded.frame_offsets(), [0, 100, 250]);
    assert_eq!(decoded.steps(), Some(&[0, 1000, 2500][..]));
    let times = decoded.times().unwrap();
    assert_eq!(times.unit(), CANONICAL_TIME_UNIT);
    assert_eq!(*times.value(), [-1.0, 24.0, 61.5]);
    assert!(decoded.fingerprint().modified().unwrap() < UNIX_EPOCH);

    // Formats without times or steps, and empty indexes, have neither.
    for (format, starts) in [
        (TrajectoryFormat::Xyz, vec![0, 100]),
        (TrajectoryFormat::Xtc, Vec::new()),
    ] {
        let original = index_of(format, starts, 400, false);
        let bytes = original.to_bytes();
        assert_eq!(bytes[13], 0);
        let decoded =
            TrajectoryFrameIndex::from_bytes(&bytes, &TrajectoryIoLimits::default()).unwrap();
        assert_eq!(decoded, original);
        assert!(decoded.times().is_none() && decoded.steps().is_none());
    }
    // DCD records steps, and times only under a header-delta policy.
    let mut stepped = index_of(TrajectoryFormat::Dcd, vec![0, 100], 400, true);
    Arc::make_mut(&mut stepped.offsets).times = None;
    for original in [
        stepped,
        index_of(TrajectoryFormat::Dcd, vec![0, 100], 400, true),
    ] {
        let decoded =
            TrajectoryFrameIndex::from_bytes(&original.to_bytes(), &TrajectoryIoLimits::default())
                .unwrap();
        assert_eq!(decoded, original);
    }
}

#[test]
fn version_1_encodings_are_refused() {
    // Version 1 had this header and the offsets alone.
    let mut bytes = index_of(TrajectoryFormat::Xtc, vec![0, 100], 400, false).to_bytes();
    bytes[8..12].copy_from_slice(&1_u32.to_le_bytes());
    reencode(&mut bytes);
    let error =
        TrajectoryFrameIndex::from_bytes(&bytes, &TrajectoryIoLimits::default()).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidIndex);
    assert!(error.to_string().contains("version 1"));
}

#[test]
fn steps_and_times_must_match_the_format() {
    let limits = TrajectoryIoLimits::default();
    let stamped = |format, stamped| index_of(format, vec![0, 100], 400, stamped).to_bytes();
    let with_flags = |mut bytes: Vec<u8>, flags| {
        bytes[13] = flags;
        reencode(&mut bytes);
        bytes
    };
    // Tables the format does not record, or missing ones it does.
    let invalid = [
        stamped(TrajectoryFormat::Xyz, true),
        with_flags(stamped(TrajectoryFormat::Xtc, true), FLAG_STEPS),
        with_flags(stamped(TrajectoryFormat::Trr, true), FLAG_TIMES),
        with_flags(stamped(TrajectoryFormat::Dcd, true), FLAG_TIMES),
        stamped(TrajectoryFormat::Xtc, false),
        with_flags(
            index_of(TrajectoryFormat::Xtc, Vec::new(), 400, false).to_bytes(),
            FLAG_STEPS | FLAG_TIMES,
        ),
    ];
    for bytes in invalid {
        assert_eq!(kind(&bytes, &limits), TrajectoryIoErrorKind::InvalidIndex);
    }
    // Non-finite times.
    for time in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut bytes = stamped(TrajectoryFormat::Xtc, true);
        let last_time = bytes.len() - CHECKSUM_BYTES - 8;
        bytes[last_time..last_time + 8].copy_from_slice(&time.to_le_bytes());
        reencode(&mut bytes);
        let error = TrajectoryFrameIndex::from_bytes(&bytes, &limits).unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidIndex);
        assert!(error.to_string().contains("time 1"));
    }
}

#[test]
fn every_frame_must_record_the_stamps_of_frame_0() {
    let limits = TrajectoryIoLimits::default();
    let stamp = FrameStamp {
        time: Some(1.5),
        step: Some(3),
    };
    let mut offsets = FrameOffsets::empty(0);
    offsets
        .push(0, 10, stamp, &limits, TrajectoryFormat::Xtc, "x")
        .unwrap();
    for different in [
        FrameStamp::default(),
        FrameStamp {
            time: None,
            ..stamp
        },
    ] {
        let error = offsets
            .push(10, 20, different, &limits, TrajectoryFormat::Xtc, "x")
            .unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::InconsistentMetadata);
        assert_eq!(error.frame(), Some(1));
    }
    offsets
        .push(10, 20, stamp, &limits, TrajectoryFormat::Xtc, "x")
        .unwrap();
    assert_eq!(
        (offsets.times(), offsets.steps()),
        (Some(&[1.5, 1.5][..]), Some(&[3, 3][..]))
    );

    let verify = |stamp| {
        offsets.verify_stamp(
            1,
            stamp,
            TrajectoryIoOperation::ReadFrame,
            TrajectoryFormat::Xtc,
            "x",
        )
    };
    verify(stamp).unwrap();
    for different in [
        FrameStamp {
            time: Some(1.25),
            ..stamp
        },
        FrameStamp {
            step: Some(4),
            ..stamp
        },
        FrameStamp::default(),
    ] {
        let error = verify(different).unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch);
        assert_eq!(error.frame(), Some(1));
    }
}

#[test]
fn malformed_encodings_are_typed_errors() {
    let limits = TrajectoryIoLimits::default();
    let bytes = index(vec![0, 100], 400).to_bytes();
    for length in 0..bytes.len() {
        assert_eq!(
            kind(&bytes[..length], &limits),
            TrajectoryIoErrorKind::InvalidIndex
        );
    }
    let mut corrupt = bytes.clone();
    corrupt[60] ^= 1;
    assert_eq!(kind(&corrupt, &limits), TrajectoryIoErrorKind::InvalidIndex);

    let mut version = bytes.clone();
    version[8] = 3;
    reencode(&mut version);
    let error = TrajectoryFrameIndex::from_bytes(&version, &limits).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidIndex);
    assert!(error.to_string().contains("version 3"));

    for (position, value) in [(12, 9), (13, 2), (14, 2), (15, 1)] {
        let mut field = bytes.clone();
        field[position] = value;
        reencode(&mut field);
        assert_eq!(kind(&field, &limits), TrajectoryIoErrorKind::InvalidIndex);
    }
    let mut flags_on_xtc = bytes.clone();
    flags_on_xtc[12] = 4;
    reencode(&mut flags_on_xtc);
    assert_eq!(
        kind(&flags_on_xtc, &limits),
        TrajectoryIoErrorKind::InvalidIndex
    );

    for starts in [vec![100, 100], vec![100, 0], vec![0, 400]] {
        let bytes = index(starts, 400).to_bytes();
        assert_eq!(kind(&bytes, &limits), TrajectoryIoErrorKind::InvalidIndex);
    }
}

#[test]
fn declared_frame_counts_are_bounded_before_allocation() {
    let bytes = index(vec![0, 100, 200], 400).to_bytes();
    let error = TrajectoryFrameIndex::from_bytes(
        &bytes,
        &TrajectoryIoLimits::default().with_max_index_entries(2),
    )
    .unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::LimitExceeded);
    assert_eq!(error.limit_exceeded().unwrap().observed(), 3);

    let mut huge = bytes;
    huge[48..56].copy_from_slice(&u64::MAX.to_le_bytes());
    reencode(&mut huge);
    assert_eq!(
        kind(&huge, &TrajectoryIoLimits::default()),
        TrajectoryIoErrorKind::LimitExceeded
    );
    let mut short = index(vec![0, 100, 200], 400).to_bytes();
    short[48..56].copy_from_slice(&2_u64.to_le_bytes());
    reencode(&mut short);
    assert_eq!(
        kind(&short, &TrajectoryIoLimits::default()),
        TrajectoryIoErrorKind::InvalidIndex
    );
}

#[test]
fn fingerprint_encodings_round_trip_and_match_the_index_fields() {
    let index = index(vec![0, 100], 400);
    let fingerprint = index.fingerprint();
    let bytes = fingerprint.to_bytes();
    assert_eq!(bytes.len(), FINGERPRINT_BYTES);
    assert_eq!(bytes[..8], FINGERPRINT_MAGIC);
    assert_eq!(
        bytes[FINGERPRINT_FIELDS],
        index.to_bytes()[FINGERPRINT_FIELDS]
    );
    assert_eq!(
        &TrajectoryFileFingerprint::from_bytes(&bytes).unwrap(),
        fingerprint
    );
    let text = fingerprint.to_string();
    assert_eq!(text.len(), 2 * FINGERPRINT_BYTES);
    assert_eq!(
        &text.parse::<TrajectoryFileFingerprint>().unwrap(),
        fingerprint
    );

    let absent = TrajectoryFileFingerprint {
        modified: None,
        ..fingerprint.clone()
    };
    let decoded = TrajectoryFileFingerprint::from_bytes(&absent.to_bytes()).unwrap();
    assert_eq!(decoded, absent);
    assert!(decoded.same_content(fingerprint));
    assert_ne!(&decoded, fingerprint);
}

#[test]
fn malformed_fingerprints_are_typed_errors() {
    let fingerprint = index(vec![0], 400).fingerprint().clone();
    let bytes = fingerprint.to_bytes();
    let kind = |bytes: &[u8]| {
        TrajectoryFileFingerprint::from_bytes(bytes)
            .unwrap_err()
            .kind()
    };
    for length in 0..bytes.len() {
        assert_eq!(kind(&bytes[..length]), TrajectoryIoErrorKind::InvalidIndex);
    }
    let mut long = bytes.clone();
    long.extend_from_slice(&[0; 8]);
    reencode(&mut long);
    assert_eq!(kind(&long), TrajectoryIoErrorKind::InvalidIndex);
    let mut corrupt = bytes.clone();
    corrupt[20] ^= 1;
    assert_eq!(kind(&corrupt), TrajectoryIoErrorKind::InvalidIndex);
    assert_eq!(
        kind(&index(vec![0], 400).to_bytes()),
        TrajectoryIoErrorKind::InvalidIndex
    );
    for (position, value) in [(8, 2), (12, 1), (13, 1), (14, 2), (15, 1)] {
        let mut field = bytes.clone();
        field[position] = value;
        reencode(&mut field);
        assert_eq!(kind(&field), TrajectoryIoErrorKind::InvalidIndex);
    }

    let text = fingerprint.to_string();
    for bad in [
        String::new(),
        text[1..].to_owned(),
        text.to_uppercase(),
        format!("{}g", &text[1..]),
        format!("{}é", &text[2..]),
        format!(
            "{}{}",
            &text[..text.len() - 1],
            if text.ends_with('0') { '1' } else { '0' }
        ),
    ] {
        assert_eq!(
            bad.parse::<TrajectoryFileFingerprint>().unwrap_err().kind(),
            TrajectoryIoErrorKind::InvalidIndex,
            "{bad}"
        );
    }
}

#[test]
fn sample_ranges_cover_small_files_whole_and_stay_inside_large_ones() {
    assert_eq!(sample_ranges(0).collect::<Vec<_>>(), [(0, 0)]);
    assert_eq!(
        sample_ranges(WHOLE_FILE_HASH_BYTES).collect::<Vec<_>>(),
        [(0, WHOLE_FILE_HASH_BYTES)]
    );
    let large = WHOLE_FILE_HASH_BYTES + 1;
    let ranges = sample_ranges(large).collect::<Vec<_>>();
    assert_eq!(ranges.len(), 18);
    assert!(ranges
        .iter()
        .all(|&(start, length)| start + length <= large));
}

#[test]
fn frame_ends_must_meet_the_next_start() {
    let offsets = FrameOffsets {
        starts: vec![10, 30],
        end: 50,
        times: None,
        steps: None,
    };
    assert_eq!(offsets.frame(0), Some((10, 30)));
    assert_eq!(offsets.frame(1), Some((30, 50)));
    assert_eq!(offsets.frame(2), None);
    offsets
        .verify_frame_end(1, 50, TrajectoryFormat::Xtc, "x")
        .unwrap();
    let error = offsets
        .verify_frame_end(0, 29, TrajectoryFormat::Xtc, "x")
        .unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch);
    assert_eq!((error.expected(), error.actual()), (Some(30), Some(29)));
}
