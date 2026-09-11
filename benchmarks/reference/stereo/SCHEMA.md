# Stereo benchmark schema 2

The current goldens preserve externally supplied input hashes and complete
assertions. RDKit 2026.03.6 is an optional reference, never a runtime dependency.

- `stereo.representation` compares canonical stereo elements and groups.
  Source text, bond marks, and assertion provenance are retained separately in
  `document`; carrier reordering preserves the represented configuration.
- `stereo.perception` also compares candidates and inference reports using actual
  Model coordinates. Imported source assertions are already represented;
  `created_element_indices` identifies only explicit coordinate additions.
- `io.smiles.isomeric` checks every source record, including non-stereo controls
  and failures. Actual emitted-and-reparsed graph, hydrogen-declaration, valence,
  neighbor, and CIP fields are compared with the independently parsed source.
  If normalization creates charge on an inferred-H atom, mandatory bracket syntax
  fixes its total H count. Optional metal-neighbor brackets are not imposed.
  `reference_evidence` retains complete original source and RDKit-emission
  projections, required declaration changes, and whole-graph stereo checks.

The `stereo.cip` goldens retain their original RDKit 2026.03.3 reference.
