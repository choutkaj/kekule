//! The files of every readable format version cover every kind that
//! version has.

use std::path::PathBuf;

use crate::io::native::{NativeKind, NATIVE_FORMAT_OLDEST, NATIVE_FORMAT_VERSION};

#[test]
fn every_readable_version_has_files_of_every_kind_it_has() {
    let fixtures = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/native");
    for version in NATIVE_FORMAT_OLDEST..=NATIVE_FORMAT_VERSION {
        let directory = fixtures.join(format!("v{version}"));
        // The kind code is the header's thirteenth byte.
        let stored = std::fs::read_dir(&directory)
            .unwrap_or_else(|error| panic!("{}: {error}", directory.display()))
            .map(|entry| std::fs::read(entry.unwrap().path()).unwrap()[12])
            .collect::<Vec<_>>();
        let kinds = (0..=u8::MAX)
            .filter_map(NativeKind::from_code)
            .filter(|kind| kind.since() <= version);
        for kind in kinds {
            assert!(
                stored.contains(&kind.code()),
                "v{version} has no file of kind {}, {kind}",
                kind.code()
            );
        }
    }
}
