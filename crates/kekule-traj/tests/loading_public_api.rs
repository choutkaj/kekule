use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use kekule::properties::{PropertyKey, PropertyValue};
use kekule::units::{DIMENSIONLESS, NANOMETER, PICOSECOND};
use kekule_traj::io::trr::{TrrLambdaPolicy, TrrReadOptions, TRR_LAMBDA_PROPERTY};
use kekule_traj::io::xyz::XyzReadOptions;
use kekule_traj::io::{
    read_trajectory, read_trajectory_with_options, TrajectoryFormatHint, TrajectoryIoLimits,
    TrajectoryOpenOptions,
};
use kekule_traj::{TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat};

mod support;
use support::{codec_kind, linear_carbon_topology, topology};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

struct TemporaryInput(PathBuf);

impl TemporaryInput {
    fn new(contents: &[u8]) -> Self {
        static SEQUENCE: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "kekule-loading-{}-{}.xyz",
            std::process::id(),
            SEQUENCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::write(&path, contents).unwrap();
        Self(path)
    }
}

impl Drop for TemporaryInput {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

#[test]
fn loading_detects_every_supported_format_and_retains_frame_order_and_topology() {
    let water = topology(&["O", "H", "H"], &[(0, 1), (0, 2)]);
    let three_atoms = topology(&["C", "H", "O"], &[(0, 1), (0, 2)]);
    for (name, topology, first_xs) in [
        ("ase-3.26.0-water.xyz", water, [0.0, 0.01]),
        (
            "mdanalysis-2.9.0-three-atoms.dcd",
            three_atoms.clone(),
            [0.0, 0.1],
        ),
        ("mdanalysis-2.9.0-three-atoms.trr", three_atoms, [0.0, 0.1]),
        (
            "mdanalysis-2.9.0-twelve-atoms.xtc",
            linear_carbon_topology(12),
            [0.0, 0.001],
        ),
    ] {
        let trajectory = read_trajectory(fixture(name), topology.clone()).unwrap();
        assert!(Arc::ptr_eq(&topology, &trajectory.shared_topology()));
        assert_eq!(trajectory.len(), 2, "{name}");
        for (frame, expected_x) in trajectory.frames().zip(first_xs) {
            assert_eq!(frame.positions().len(), topology.atom_count());
            let actual_x = frame.positions().values().value()[0].x;
            assert!((actual_x - expected_x).abs() < 1.0e-6, "{name}: {actual_x}");
        }
    }
}

#[test]
fn loading_retains_trr_cell_vectors_time_step_and_properties() {
    let topology = topology(&["C", "H", "O"], &[(0, 1), (0, 2)]);
    let trajectory =
        read_trajectory(fixture("mdanalysis-2.9.0-three-atoms.trr"), topology).unwrap();
    let lambda_key = PropertyKey::new(TRR_LAMBDA_PROPERTY).unwrap();
    for (index, frame) in trajectory.frames().enumerate() {
        let index = index as f64;
        assert!(frame.cell().is_some());
        let velocities = frame.velocities().unwrap();
        let forces = frame.forces().unwrap();
        assert_eq!(velocities.value().len(), 3);
        assert_eq!(forces.value().len(), 3);
        assert!((velocities.value()[0].x - (index + 0.5) * 0.1).abs() < 1.0e-6);
        assert!((forces.value()[0].x - (index + 1.0) * 10.0).abs() < 1.0e-5);
        assert_eq!(
            frame.time().unwrap().value_in(PICOSECOND).unwrap(),
            index * 0.25
        );
        assert_eq!(frame.step(), Some(index as u64));
        assert_eq!(
            frame.properties().get(&lambda_key),
            Some(&PropertyValue::Real {
                value: 0.125 + index * 0.125,
                unit: DIMENSIONLESS,
            })
        );
    }
}

#[test]
fn loading_accepts_owned_topology_and_honors_format_units_and_sequential_limits() {
    let mut input = TemporaryInput::new(b"1\nfirst\nC 1 2 3\n1\nsecond\nC 4 5 6\n");
    let mislabeled_path = input.0.with_extension("xtc");
    fs::rename(&input.0, &mislabeled_path).unwrap();
    input.0 = mislabeled_path;
    // An explicit format must override the conflicting file extension.
    assert!(read_trajectory(&input.0, linear_carbon_topology(1)).is_err());
    let options = TrajectoryOpenOptions::default()
        .with_format_hint(TrajectoryFormatHint::Explicit(TrajectoryFormat::Xyz))
        .with_xyz_options(XyzReadOptions::default().with_length_unit(NANOMETER))
        .with_limits(TrajectoryIoLimits {
            max_frames: 2,
            // Whole-file loading must not require an index.
            max_index_entries: 0,
            max_index_bytes: 0,
            ..Default::default()
        });
    let owned = Arc::try_unwrap(linear_carbon_topology(1)).unwrap();
    let trajectory = read_trajectory_with_options(&input.0, owned, options).unwrap();
    assert_eq!(trajectory.len(), 2);
    assert_eq!(
        trajectory.frame(0).unwrap().positions().values().value()[0].x,
        1.0
    );
    assert_eq!(
        trajectory.frame(1).unwrap().positions().values().value()[0].x,
        4.0
    );
    for frame in trajectory.frames() {
        assert!(frame.cell().is_none());
        assert!(frame.velocities().is_none());
        assert!(frame.forces().is_none());
        assert!(frame.time().is_none());
        assert!(frame.step().is_none());
        assert!(frame.properties().is_empty());
    }
}

#[test]
fn loading_reports_a_late_decode_error_instead_of_returning_a_prefix() {
    let input = TemporaryInput::new(b"1\ncomplete\nC 1 2 3\n1\ntruncated\n");
    let error = read_trajectory(&input.0, linear_carbon_topology(1)).unwrap_err();
    let TrajectoryError::Codec(context) = error else {
        panic!("expected a codec error, got {error}");
    };
    assert_eq!(context.kind(), TrajectoryCodecErrorKind::TruncatedRecord);
    assert_eq!(context.frame(), Some(1));
    assert_eq!(context.format(), Some(TrajectoryFormat::Xyz));
    assert_eq!(context.source_label(), Some(input.0.to_str().unwrap()));
}

#[test]
fn loading_preserves_resource_limits_and_topology_validation() {
    let path = fixture("ase-3.26.0-water.xyz");
    let water = topology(&["O", "H", "H"], &[(0, 1), (0, 2)]);
    let options = TrajectoryOpenOptions::default().with_limits(TrajectoryIoLimits {
        max_frames: 1,
        ..Default::default()
    });
    let error = read_trajectory_with_options(&path, water, options).unwrap_err();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
    );
    // Equal counts do not bypass XYZ's element-order check.
    assert!(read_trajectory(&path, linear_carbon_topology(3)).is_err());
    assert!(read_trajectory(&path, linear_carbon_topology(2)).is_err());
}

#[test]
fn loading_preserves_codec_policy_and_open_errors() {
    let topology = topology(&["C", "H", "O"], &[(0, 1), (0, 2)]);
    let options = TrajectoryOpenOptions::default().with_trr_options(
        TrrReadOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
    );
    let error = read_trajectory_with_options(
        fixture("mdanalysis-2.9.0-three-atoms.trr"),
        topology.clone(),
        options,
    )
    .unwrap_err();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::InconsistentMetadata)
    );

    let input = TemporaryInput::new(b"");
    fs::remove_file(&input.0).unwrap();
    assert!(matches!(
        read_trajectory(&input.0, topology),
        Err(TrajectoryError::Io(_))
    ));
}

#[test]
fn loading_rejects_an_empty_file_even_with_an_explicit_format() {
    // Path readers require a first frame, even when format detection is bypassed.
    let input = TemporaryInput::new(b"");
    let topology = linear_carbon_topology(1);
    let options = TrajectoryOpenOptions::default()
        .with_format_hint(TrajectoryFormatHint::Explicit(TrajectoryFormat::Xyz));
    let error = read_trajectory_with_options(&input.0, topology, options).unwrap_err();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::InvalidHeader)
    );
}
