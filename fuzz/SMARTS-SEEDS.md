# SMARTS fuzz seeds

The ten seed patterns are copied from the pinned external SMARTS fixtures.
Matcher seeds pair each pattern with the existing PubChem ethanol target (CID 702).
Source provenance and licenses are in `../benchmarks/smarts-fixtures`.

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
