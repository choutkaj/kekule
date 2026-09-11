# Stereo benchmark schema 2

The schema migration is implemented by `migrate_v2.py`. It requires RDKit
2026.03.6 and never invokes the Rust implementation or reads its outputs.
`schema-v1/` retains the original compressed reference files byte for byte.
`migration-v2.json` records their SHA-256 hashes, the replacement hashes, and
record counts. Input hashes remain those of the externally supplied fixtures.

The three affected benchmarks are `stereo.representation`, `stereo.perception`,
and `io.smiles.isomeric`. The `stereo.cip` references are unchanged.

| Schema 1 assertion | Schema 2 assertion |
| --- | --- |
| `stereo_bond_marks` attached to a molecule | The same marks under `document.stereo_bond_marks`, accompanied by the retained source text and format. |
| Element `source` and `specifiedness` | `document.stereo_sources`, linked to the corresponding canonical element index. Configuration remains explicit as `orientation`, including `null` for unknown configurations. |
| Source-order tetrahedral carriers | Canonically ordered carriers with orientation inverted exactly when the carrier permutation is odd. |
| Arbitrary double-bond reference carriers | The first explicit carrier at each ordered endpoint; replacing one reference reverses the relative orientation. |
| Directional/wedge assertions assembled during general perception | Represented elements published by source interpretation and recorded in `report.assembled_elements`. |
| Coordinate-created element indices | Elements added by explicit coordinate inference using the retained Model positions. Imported source assertions are already present and are not created again. |

The previously reviewed perception golden supplies the directional and wedge
assertions that were absent from the old representation golden. This changes
the stage at which an assertion is checked while preserving its focus, carriers,
configuration, and source evidence. Complete old assertions remain available in
the archive for review.

Coordinate additions are derived independently from signed geometry. For four
tetrahedral points the reference uses the scalar triple product of the three
vectors relative to the fourth point. An implicit hydrogen or lone pair uses
the center as the fourth point; this has the same determinant sign as a virtual
ligand opposite the explicit-ligand vector sum. Double-bond relations use the
dot product of the two ligand-plane normals about the bond axis. Degenerate
geometry produces no assigned configuration. Three-coordinate sulfur and
selenium with at most one double bond are included with a lone-pair carrier.
These are spatial-carrier rules, not claims that every candidate has distinct
CIP priorities.

The isomeric writer reference compares the complete decoded roundtrip with the
independently parsed source graph and CIP-bearing stereo. The asserted fields
include `no_implicit_hydrogens`, explicit and implicit H counts, explicit valence,
and neighbors. The Rust adapter derives them from emitted and reparsed output;
it never copies them from the source or normalizes them to fixed values.

The only declaration projection follows mandatory SMILES grammar: a charged
atom requires brackets and a fixed hydrogen count. If independent chemical
normalization introduces charge on a source atom that still permits hydrogen
inference, the expected emission fixes its total hydrogen count. This rule
applies to every element; it does not apply to uncharged metal neighbors.
`reference_evidence.decoded_source` retains every pre-projection source field,
and `mandatory_bracket_declarations` records the original and required emitted
atom fields. The actual Rust output is still decoded and asserted directly.

Two pinned PubChem 1k records require this projection: CID 173868, containing
`[O-]Cl(=O)(=O)=O`, and CID 423442, containing `OCl(=O)(O)O`.
Normalization gives their chlorine atoms charges +3 and +1, respectively.
It also gives three oxo oxygens in the first case and one in the second a
charge of -1. All six affected source atoms have zero explicit and implicit H with
`no_implicit_hydrogens=false`; the required bracketed output retains zero H
and changes that flag to `true`. Their complete source and reference-emission
projections remain in the corresponding golden evidence.

RDKit chooses extra brackets around metal-adjacent atoms when writing SMILES.
That cosmetic choice can move an H count from implicit to explicit and change
`no_implicit_hydrogens` without changing the chemical graph. Schema 2 therefore
targets source declaration preservation instead of RDKit's bracket preference.
`reference_evidence` retains each literal RDKit emission and its full decoded
projection. An independent canonical isomeric graph identity check must prove
that emission chemically and stereochemically equals its source; generation
fails if this check fails. This use of RDKit canonicalization is a reference
integrity check, not a target for the Rust canonical writer.

Every source record is checked, including non-stereo controls and failures.
Old overlapping assertions remain in the archive. New expected values are derived
from the pinned reference graph, retained assertions, and the explicit transformations
described above. Canonical SMILES behavior is outside this migration.

One archived perception error, PubChem CID 80794 (tert-butyl perchlorate),
is independently accepted by RDKit 2026.03.6. Its replacement report is derived
from that sanitized graph. The other nine archived perception errors remain
errors in the pinned reference. No error expectation is promoted from a Rust
result. Isomeric provenance additionally records any changed paths in previously
asserted records; removing an old field aborts migration.

Run from the repository root with the pinned optional reference environment:

```text
uv run --offline --with rdkit==2026.3.6 --python 3.13 python benchmarks/reference/stereo/migrate_v2.py
uv run --offline --with rdkit==2026.3.6 --python 3.13 python benchmarks/reference/stereo/migrate_v2.py --write
```

Without `--write`, the migration validates inputs and computes expectations
without writing archives, active goldens, manifests, or provenance. Scoped
`--corpus` publication merges provenance entries with previous corpus runs.
This is a reviewed schema migration, not an
implementation-golden acceptance command.
