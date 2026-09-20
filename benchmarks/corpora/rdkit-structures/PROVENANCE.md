# External structure coverage

This corpus preserves 50 input files from RDKit release `Release_2026_03_3`,
commit `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`, under
[`Code/GraphMol/FileParsers/test_data`](https://github.com/rdkit/rdkit/tree/e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f/Code/GraphMol/FileParsers/test_data).

The selection was fixed before testing either implementation:

- Every input SDF variant under `atropisomers` for BMS-986142, JDQ443, Mrtx1719,
  RP-6306, Sotorasib and ZM374979: 47 files. Base structures, 3D structures,
  alternative stereo drawings and variants marked `Bad` are all retained.
- Every `chebi_*.mol` file directly under `test_data`: three V3000 inputs.
- Files containing `.expected` are excluded because they are upstream generated
  outputs, not supplied inputs. No other matching input is omitted.

These are multiple representations of a small set of named compounds, not
50 independently sampled molecules. Their purpose is targeted axial stereo,
drawing, coordinate and V3000 coverage. They do not estimate prevalence in a
chemical database. A variant marked `Bad` may intentionally contain conflicting
stereo; neither filename nor agreement determines applicability or acceptance.
Native and reference failures stay explicit benchmark outcomes.
The ChEBI inputs contain R-group placeholders. They exercise format/query
coverage, not fully specified molecular composition; a reference descriptor
computed with zero-mass dummy atoms is not a physical mass for an unspecified
substituent. Their native interpretation errors remain visible.

The 47 upstream `.sdf` files contain single MOL blocks ending at `M  END`, with
no SDF record delimiter or data fields. They are stored with `.mol` extensions
to identify their actual format; original filenames remain in their source IDs
and upstream paths. This routing was established before native comparison.
The original bytes are stored under `data/`; no atom, bond, coordinate, title,
property, line ending or delimiter was rewritten. Each source has its own ID,
path, SHA-256, Git blob SHA-1 and byte count in `sources.lock.json`. Any multiple
records in a source remain individually measured under that source ID.

`upstream/tree.json` preserves the complete recursive GitHub tree response for
the upstream test-data directory. `upstream/commit.json` pins the release commit,
and `upstream/license.txt` retains RDKit's BSD 3-Clause license. All downloaded
files were verified against their Git blob hashes before adoption. The offline
provenance test checks complete selection, exact bytes and absence of extra
fixtures without depending on either chemistry implementation or the network.

Reference observations are generated independently with pinned RDKit 2026.03.3.
Upstream expected-output files are not used as Kekule goldens. The normal
benchmark format rules select applicable features; query, SMILES-text, mmCIF and
DSSP features are not applicable to these MOL inputs. No existing corpus or
reference observation is replaced by this addition.
