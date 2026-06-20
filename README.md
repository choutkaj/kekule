<img src="./assets/molecules-character-waves-v3.svg" alt="Molecules banner." width="100%">


<div style="text-align: center;">
    <img src="./assets/molecules-logo.svg" alt="Molecules banner." width="25%" margin="auto">
</div>

<p align="center">
  <img src="./assets/molecules-logo.svg" alt="Molecules banner." width="25%" margin="auto" />
</p>



<pre style="line-height: 1; font-family: monospace;">
╭┬╮╭─╮╷  ╭─╴╭─╴╷ ╷╷  ╭─╴╭─╮
││││ ││  ├╴ │  │ ││  ├╴ ╰─╮
╵ ╵╰─╯╰─╴╰─╴╰─╴╰─╯╰─╴╰─╴╰─╯
  cheminformatics in Rust
</pre>

`molecules` is an AI-coded cheminformatics backend for small molecules and macromolecules written in pure Rust. The repo is organized around feature-scoped development: every meaningful capability is treated as a feature and is validated against established codebases (RDKit and Biopython).

## Current scaffold

- Cargo workspace with `molecules` and `xtask` crates.
- Minimal pure Rust molecule data model skeleton.
- Architecture, roadmap, agent rules, and contribution docs.
- Feature registry, generated dashboard, and feature templates.
- Codex skills for feature work and independent feature review.
- Reference-validation directories for RDKit and Biopython.

## Common commands

```bash
cargo test --workspace
cargo xtask dashboard
cargo xtask dashboard --check
cargo xtask skills --check
cargo xtask validate --feature io.smiles.parse --corpus tiny
cargo xtask validate --feature all --corpus all
```

The `cargo xtask` alias is defined in `.cargo/config.toml`.

## Codex co-author attribution

Install the repository-local commit hook once in each clone:

```bash
bash scripts/setup-codex-attribution.sh
```

The hook appends this trailer exactly once to future non-merge commits:

```text
Co-authored-by: Codex <noreply@openai.com>
```

Inspect the local configuration with:

```bash
bash scripts/show-codex-attribution-status.sh
```

Set `SKIP_AI_COAUTHORS=1` for a single commit that must omit configured AI co-authors.

## License

License is intentionally not selected yet. Choose an open-source license before public release.
