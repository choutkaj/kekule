# kekule-potentials

`kekule-potentials` provides concrete potentials for topology-bound [`kekule`](https://crates.io/crates/kekule) structural views. Its default `openff` feature evaluates OpenFF energies and gradients from a `kekule-openff` parameterization, and every potential can be minimized with L-BFGS. The default `ash` feature bundles the Ash charge model that the Rosemary preset requires (about 10 MB); enable `openff` without default features to load charge models only from files.

## Installation

```sh
cargo add kekule-potentials
```

Licensed under MIT or Apache-2.0.
