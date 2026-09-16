# Stored golden coverage

All 125 dataset/feature pairs exposed by `cargo benchmark` have complete
stored files in `goldens/`, generated without case limits.
Every locked source ID is represented, including reference failures and
cases without the required input format. These are reference-availability
counts, not Kekule agreement scores. Normal comparisons reuse these files.

The preparation added 99 files and preserved the previous 26 byte for byte.
Compressed size: 1,777.0 MiB. Reference software:

- biopython: 1.87
- biopython: Biopython 1.87 / mkdssp version 4.6.1
- rdkit: 2026.03.3

| Dataset | Source IDs | Files | Reference values | Reference errors | Missing-format cases |
| --- | ---: | ---: | ---: | ---: | ---: |
| pubchem-100k | 100,000 | 25 | 3,599,739 | 261 | 200,000 |
| enamine-diversity | 50,240 | 25 | 1,808,640 | 0 | 100,480 |
| pl-rex | 164 | 25 | 5,904 | 0 | 1,148 |
| pdb-1000 | 1,000 | 25 | 1,756 | 244 | 23,000 |
| smoke | 20 | 25 | 452 | 0 | 126 |

An input may contribute to several features and, for molecular algorithms,
through both its supplied SDF and SMILES files. Thus case counts exceed
source-ID counts. Missing formats and reference errors remain failures;
neither is silently omitted or converted into agreement.

## Feature coverage

Each cell is **successful reference values / all cases**. The denominator
includes every source ID and all matching supplied inputs.

| Feature | pubchem-100k | enamine-diversity | pl-rex | pdb-1000 | smoke |
| --- | ---: | ---: | ---: | ---: | ---: |
| `io.smiles.parse` | 99,991 / 100,000 | 50,240 / 50,240 | 0 / 164 | 0 / 1,000 | 8 / 20 |
| `io.smiles.write` | 99,991 / 100,000 | 50,240 / 50,240 | 0 / 164 | 0 / 1,000 | 8 / 20 |
| `io.smiles.canonical` | 99,991 / 100,000 | 50,240 / 50,240 | 0 / 164 | 0 / 1,000 | 8 / 20 |
| `io.smiles.isomeric` | 99,991 / 100,000 | 50,240 / 50,240 | 0 / 164 | 0 / 1,000 | 8 / 20 |
| `io.mol.parse` | 99,991 / 100,000 | 50,240 / 50,240 | 328 / 328 | 0 / 1,000 | 17 / 20 |
| `io.mol.v2000.write` | 99,991 / 100,000 | 50,240 / 50,240 | 328 / 328 | 0 / 1,000 | 17 / 20 |
| `io.mol.v3000.write` | 99,991 / 100,000 | 50,240 / 50,240 | 328 / 328 | 0 / 1,000 | 17 / 20 |
| `io.sdf.parse` | 99,991 / 100,000 | 50,240 / 50,240 | 328 / 328 | 0 / 1,000 | 17 / 20 |
| `io.sdf.v2000.write` | 99,991 / 100,000 | 50,240 / 50,240 | 328 / 328 | 0 / 1,000 | 17 / 20 |
| `io.mmcif.parse` | 0 / 100,000 | 0 / 50,240 | 0 / 164 | 1,000 / 1,000 | 1 / 20 |
| `algo.rings.fast` | 200,000 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `algo.rings.sssr` | 200,000 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `algo.valence.rdkit-like` | 200,000 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `algo.aromaticity.rdkit-like` | 199,982 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `algo.canonical-ranking` | 199,982 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `algo.substructure.vf2` | 199,982 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `query.smarts` | 100,000 / 100,000 | 50,240 / 50,240 | 0 / 164 | 0 / 1,000 | 8 / 20 |
| `chem.perception.default` | 199,982 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `chem.hydrogen-transforms` | 199,982 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `descriptor.molecular` | 199,982 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `descriptor.rotatable-bonds.rdkit-strict` | 199,982 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `stereo.representation` | 199,982 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `stereo.perception` | 199,982 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `stereo.cip` | 199,982 / 200,000 | 100,480 / 100,480 | 328 / 328 | 0 / 1,000 | 25 / 26 |
| `bio.secondary-structure.dssp` | 0 / 100,000 | 0 / 50,240 | 0 / 164 | 756 / 1,000 | 1 / 20 |

## Retained reference failures

RDKit failed on nine PubChem source IDs across the affected features, producing
261 retained reference-error cases. Those inputs remain in the full case sets.

PDB/DSSP produced 756 successful reference values and 244 errors.
Of those error inputs, 230 contain only nucleic-acid polymers.
Direct inspection of the other 14 reference evaluations found six
with no analyzable residues, seven where DSSP produced no output,
and one Biopython failure on a missing `_atom_site.pdbx_PDB_ins_code`
column (9A18). These cases remain in the stored file as errors.
The current runner reports these as `missing or empty residues`,
so its stored message does not preserve the underlying cause.

## Verification

Every compressed JSONL file was read completely and checked for exact
source membership, record indices, input SHA-256 values, duplicate cases,
valid outcomes and consistent reference versions. Generation reports
confirm that no Kekule evaluations were used to produce expectations.

Source membership remains the existing locked corpus membership. The
historical PubChem preselection limitation is described in [GUIDE.md](GUIDE.md).

Use the commands in [GUIDE.md](GUIDE.md) to compare against these files.
Writer validation still requires RDKit to read newly emitted Kekule text.
It does not recalculate the stored expectations.
