# Releasing

The workspace publishes five crates to crates.io. Four share the workspace
version; `kekule-openff-ash` carries its own version and is released only when
its model data changes. Publishing is done by a maintainer from a clean
checkout of `main`; nothing in CI publishes.

## 1. Choose the version

Before 1.0, any breaking change bumps the minor version (0.4 → 0.5) and
everything else bumps the patch version. Run the informational semver check to
see what the release breaks:

```bash
cargo semver-checks check-release --workspace
```

## 2. Bump the version sites

On a release branch based on current `main`:

| Site | What changes |
| --- | --- |
| `Cargo.toml` `[workspace.package] version` | The shared version of `kekule`, `kekule-traj`, `kekule-openff`, `kekule-potentials` and the unpublished benchmark observer |
| `Cargo.toml` `[workspace.dependencies]` | The `version` requirement of `kekule`, `kekule-traj` and `kekule-openff`, kept equal to the workspace version |
| `crates/kekule-openff-ash/Cargo.toml` `version` | Only when the Ash model data changes |
| `Cargo.toml` `[workspace.dependencies] kekule-openff-ash` | Only when a new Ash version is released |
| `Cargo.lock` and `fuzz/Cargo.lock` | Regenerated: `cargo update -w` and `cargo update --manifest-path fuzz/Cargo.toml -p kekule -p kekule-traj` |

Members inherit their requirements with `kekule.workspace = true` and so on, so
no member manifest names a version. Check that nothing else still names the
old version:

```bash
git grep -n '<old version>' -- '*.toml' '*.md' ':!CHANGELOG.md'
```

## 3. Update the CHANGELOG

`CHANGELOG.md` follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

1. Rename `## [Unreleased]` to `## [X.Y.Z] - YYYY-MM-DD` and add a new, empty
   `## [Unreleased]` above it.
2. Add a short summary paragraph under the release heading, and name any crate
   that is not released with it (usually `kekule-openff-ash`).
3. Check that every merged pull request since the last tag has an entry, under
   the headings Added, Changed, Deprecated, Removed, Fixed and Security in that
   order. Mark breaking entries **Breaking:**. A fix to a feature that was never
   released belongs in that feature's Added entry, not under Fixed.
4. Update the link references at the bottom:
   `[Unreleased]: https://github.com/choutkaj/kekule/compare/vX.Y.Z...HEAD` and
   `[X.Y.Z]: https://github.com/choutkaj/kekule/compare/vPREVIOUS...vX.Y.Z`.

## 4. Run the gates

Run every gate from `AGENTS.md` on the release branch, plus the dependency
audits from `SECURITY.md`:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-features --locked
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --locked
cargo +1.89 check --workspace --all-targets --all-features --locked
cargo deny check
cargo audit --file Cargo.lock
cargo audit --file fuzz/Cargo.lock --no-fetch
cargo package --workspace --exclude kekule-bench --locked
```

`cargo package --workspace` resolves the unpublished member versions from the
packages it builds, so it verifies each crate as it will be uploaded.

## 5. Check package sizes and contents

crates.io rejects a `.crate` above 10,485,760 bytes. CI fails above that and
warns above 90 %. Check the sizes and review each file list; every crate
declares an explicit `include` list, and integration tests are not shipped
(cargo notes each ignored test target):

```bash
ls -l target/package/*.crate
cargo package -p kekule --list --locked
```

`kekule-openff-ash` is the largest crate, close to the limit. A model that does
not fit must be split or compressed differently before release.

Merge the release branch into `main` through a pull request once CI passes.

## 6. Publish in dependency order

Publish from the merged `main` commit. A crate can be published only after
every crate it depends on, including through dev-dependencies with a version
requirement, is visible on crates.io:

1. `kekule`, and `kekule-openff-ash` only if its version changed. Publishing an
   unchanged Ash version fails, because crates.io does not accept a version
   twice; skip it.
2. `kekule-openff` and `kekule-traj`.
3. `kekule-potentials`, which depends on `kekule-openff` and has a
   dev-dependency on `kekule-traj`.

```bash
cargo publish -p kekule --locked
cargo publish -p kekule-openff --locked
cargo publish -p kekule-traj --locked
cargo publish -p kekule-potentials --locked
```

If a publish fails partway, fix the cause in a new patch version; a published
version cannot be replaced, only yanked.

## 7. Tag and announce

Tag the published commit and push the tag:

```bash
git tag -a vX.Y.Z -m "Kekule X.Y.Z"
git push origin vX.Y.Z
```

Create a GitHub release from the tag with the version's CHANGELOG section as its
notes. When `kekule-openff-ash` is released on its own, tag it
`kekule-openff-ash-vX.Y.Z`.

## 8. After the release

- Check that docs.rs built every published crate, including the feature-gated
  items.
- Install the published crates in a scratch project (`cargo add kekule@X.Y.Z`)
  and build one README example against them.
- Confirm that the compare links in `CHANGELOG.md` resolve.
- Delete the merged release branch.
