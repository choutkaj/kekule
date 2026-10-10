# Small-molecule dataset sources and licences

Each file keeps the bytes its source served, apart from the extraction noted
below. `manifest.jsonl` gives every file's upstream URL and SHA-256, and
`selection.lock.json` records the releases and pools selection drew from.
These files are benchmark inputs only. ChEMBL records in particular must not
be copied into `crates/*/tests/fixtures`, whose contents ship under the
crates' MIT/Apache licences.

| Directory | Source | Licence and attribution |
| --- | --- | --- |
| `chembl/` | [ChEMBL](https://www.ebi.ac.uk/chembl/), EMBL-EBI; release in `selection.lock.json` | [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/). Zdrazil et al., *Nucleic Acids Res.* 52, D1180 (2024). The molfile (`.mol`) and canonical SMILES (`.smi`, plus a newline) are the unmodified `molecule_structures` fields of the API record |
| `ccd/` | [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd), ideal coordinates served by RCSB PDB | [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) |
| `rdkit/` | [RDKit](https://github.com/rdkit/rdkit) `Release_2026_03_3`, commit `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f` | BSD 3-Clause; see `LICENSE-rdkit.txt` |
