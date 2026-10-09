# kekule-openff-ash

The OpenFF **Ash** NAGL partial-charge model (`openff-gnn-am1bcc-1.0.0`),
packaged as data for [`kekule-openff`](https://crates.io/crates/kekule-openff).
It is the model the bundled Rosemary force field requires.

Use it through `kekule-openff`, which depends on this crate by default:

```rust,ignore
let model = kekule_openff::NaglModel::ash()?;
```

The model is licensed under CC BY 4.0; see [ATTRIBUTION.md](ATTRIBUTION.md) for
the source, checksums, and conversion. The Rust code is MIT OR Apache-2.0.
