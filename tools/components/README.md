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
| `test_provenance.py` | Offline check that the shipped files, the lock and the Rust bundle list agree |

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

## Reproduce

Run from the repository root. Upstream entries change over time, so a fresh
download may differ from the lock; review such changes like any data update.

```sh
python tools/components/fetch.py
python -m unittest discover -s tools/components -p "test_*.py"
```
