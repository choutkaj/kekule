use std::fs;

use super::*;

#[test]
fn forbid_requires_hard_links_and_leaves_no_probe_behind() {
    let directory = std::env::temp_dir().join(format!(
        "kekule-publish-{}-{}",
        std::process::id(),
        TEMPORARY_FILE_SEQUENCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir(&directory).unwrap();
    let temporary = directory.join(".out.xyz.1.2.tmp");
    fs::write(&temporary, b"trajectory").unwrap();

    require_hard_links(&temporary, |from, to| fs::hard_link(from, to)).unwrap();
    let error = require_hard_links(&temporary, |_, _| {
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "filesystem has no hard links",
        ))
    })
    .unwrap_err();
    assert_eq!(
        (error.kind(), error.to_string()),
        (
            io::ErrorKind::Unsupported,
            "filesystem has no hard links".to_owned()
        )
    );
    let mut names = fs::read_dir(&directory)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect::<Vec<_>>();
    names.sort();
    assert_eq!(names, [".out.xyz.1.2.tmp"]);
    fs::remove_dir_all(directory).unwrap();
}
