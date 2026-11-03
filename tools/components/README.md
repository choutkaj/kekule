# Bundled component templates

`crates/kekule/data/components` holds unmodified entries of the wwPDB
[Chemical Component Dictionary](https://www.wwpdb.org/data/ccd) (CCD), one
file per component as served by RCSB at
`https://files.rcsb.org/ligands/view/<ID>.cif`. Kekule embeds them as the
templates of `kekule::components::ComponentLibrary::bundled`. wwPDB archive data,
including the CCD, are available under CC0 1.0.

| File | Purpose |
| --- | --- |
| `fetch.py` | Selection; downloads the entries and writes the lock |
| `sources.lock.json` | URL, byte count, SHA-256 and dictionary modification date of each entry |
| `test_provenance.py` | Offline check that the shipped files, the locks and the Rust bundle list agree |
| `codes.py` | Generates the one-letter code table from the complete dictionary |
| `codes.lock.json` | URL, SHA-256 and latest modification date of the dictionary the table comes from, and the table's SHA-256 |

## Selection

The selection was fixed before use, by role rather than by observed frequency:

- Amino acids: the 20 standard residues, selenocysteine (`SEC`),
  pyrrolysine (`PYL`) and selenomethionine (`MSE`).
- Nucleotides: `A`, `C`, `G`, `U`, `DA`, `DC`, `DG`, `DT`.
- Water: `HOH` and heavy water `DOD`.
- Ions: `ZN`, `MG`, `CA`, `NA`, `K`, `CL`, `MN`, `FE`, `FE2`, `CU`, `CU1`,
  `NI`, `CO`, `CD`, `HG`, `IOD`, `BR`.
- Cofactors: `HEM`, `HEC`, `FAD`, `FMN`, `NAD`, `NAI`, `NAP`, `NDP`, `ATP`,
  `ADP`, `AMP`, `GTP`, `GDP`, `SAM`, `SAH`, `COA`, `ACO`, `PLP`, `TPP`, `BTN`,
  `SF4`, `FES`, `F3S`, `CLA`.
- Crystallization additives: `SO4`, `PO4`, `GOL`, `EDO`, `ACT`, `PEG`, `DMS`,
  `MPD`.

The dictionary defines histidine only as the imidazolium cation. Kekule derives
its two neutral tautomers in code (`crates/kekule/src/io/components/bundled.rs`);
they are not dictionary entries.

## One-letter codes

`crates/kekule/data/components/codes.tsv` lists every dictionary component
whose type links it into a polypeptide or a polynucleotide (about 2,700),
with its one-letter code (`_chem_comp.one_letter_code`) and parents
(`_chem_comp.mon_nstd_parent_comp_id`), so that modified residues the bundle
does not hold still have codes. `codes.py` generates it from the complete
dictionary, `components.cif.gz`, and records the source in
`codes.lock.json`.

## Reproduce

Run from the repository root. Upstream entries change over time, so a fresh
download may differ from the lock; review such changes like any data update.

```sh
python tools/components/fetch.py
python tools/components/codes.py
python -m unittest discover -s tools/components -p "test_*.py"
```
