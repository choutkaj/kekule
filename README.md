
<p align="center">
  <img src="assets/molecules_logo.png" alt="Molecules!" width="400">
</p>


<pre style="line-height: 1; font-family: monospace;">
╭┬╮╭─╮╷  ╭─╴╭─╴╷ ╷╷  ╭─╴╭─╮
││││ ││  ├╴ │  │ ││  ├╴ ╰─╮
╵ ╵╰─╯╰─╴╰─╴╰─╴╰─╯╰─╴╰─╴╰─╯
cheminformatics in pure Rust
</pre>


<pre style="line-height: 1; font-family: monospace;">
┏┳┓┏━┓╻  ┏━╸┏━╸╻ ╻╻  ┏━╸┏━┓
┃┃┃┃ ┃┃  ┣╸ ┃  ┃ ┃┃  ┣╸ ┗━┓
╹ ╹┗━┛┗━╸┗━╸┗━╸┗━┛┗━╸┗━╸┗━┛
cheminformatics in pure Rust
</pre>

<pre style="line-height: 1; font-family: monospace;">
┳┳┓┏┓┓ ┏┓┏┓┳┳┓ ┏┓┏┓
┃┃┃┃┃┃ ┣ ┃ ┃┃┃ ┣ ┗┓
┛ ┗┗┛┗┛┗┛┗┛┗┛┗┛┗┛┗┛

</pre>

<pre style="line-height: 1; font-family: monospace;">
░█▄█░█▀█░█░░░█▀▀░█▀▀░█░█░█░░░█▀▀░█▀▀
░█░█░█░█░█░░░█▀▀░█░░░█░█░█░░░█▀▀░▀▀█
░▀░▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀
 cheminformatics in pure Rust
 C─C═C─N─O
</pre>

<pre style="line-height: 1; font-family: monospace;">
███╗   ███╗ ██████╗ ██╗     ███████╗ ██████╗██╗   ██╗██╗     ███████╗███████╗
████╗ ████║██╔═══██╗██║     ██╔════╝██╔════╝██║   ██║██║     ██╔════╝██╔════╝
██╔████╔██║██║   ██║██║     █████╗  ██║     ██║   ██║██║     █████╗  ███████╗
██║╚██╔╝██║██║   ██║██║     ██╔══╝  ██║     ██║   ██║██║     ██╔══╝  ╚════██║
██║ ╚═╝ ██║╚██████╔╝███████╗███████╗╚██████╗╚██████╔╝███████╗███████╗███████║
╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚══════╝╚══════╝
 cheminformatics in pure Rust
</pre>

</pre>



<pre style="line-height: 1; font-family: monospace;">



`molecules` is a pure Rust cheminformatics and molecular-structure backend for small molecules and macromolecules.

The repo is organized around feature-scoped development: every meaningful capability has a feature directory with machine-readable `feature.toml` metadata and one canonical human-readable `feature.md`. RDKit and Biopython are reference implementations for validation only; they are not runtime dependencies of the Rust library.

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
cargo xtask validate --feature core.graph
```

The `cargo xtask` alias is defined in `.cargo/config.toml`.

## License

License is intentionally not selected yet. Choose an open-source license before public release.
