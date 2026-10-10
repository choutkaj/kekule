# SMARTS fuzz seeds

The ten seed patterns are copied verbatim from pinned external SMARTS tables.
Matcher seeds pair each pattern with the PubChem ethanol target (CID 702).
Source locks and licenses are in `provenance/`:

- `provenance/rdkit/`: RDKit 2026.03.3 `Functional_Group_Hierarchy.txt`, pinned
  by the conda-forge `librdkit` archive and file SHA-256 in `sources.lock.json`
  (BSD 3-Clause, `license.txt`). Source rows count the table's non-comment
  query rows from 1.
- `provenance/openff/`: OpenFF force-field and toolkit OFFXML files, pinned by
  commit, URL and SHA-256 in `sources.lock.json` (`LICENSE-forcefields`,
  `LICENSE-toolkit`).

The locks are the original corpus locks; only their upstream entries apply.

| Seed | Source | Source row / parameter |
| --- | --- | --- |
| external-00 | Functional_Group_Hierarchy.smarts | 2 |
| external-01 | Functional_Group_Hierarchy.smarts | 3 |
| external-02 | Functional_Group_Hierarchy.smarts | 5 |
| external-03 | Functional_Group_Hierarchy.smarts | 6 |
| external-04 | chargeincrement-test.offxml | 0 |
| external-05 | chargeincrement-test.offxml | 1 |
| external-06 | chargeincrement-test.offxml | 2 |
| external-07 | chargeincrement-test.offxml | 3 |
| external-08 | openff-2.2.0.offxml | 0 |
| external-09 | openff-2.2.0.offxml | 1 |

Other files created by fuzz execution are mutations of these seeds, not scientific comparison fixtures.
