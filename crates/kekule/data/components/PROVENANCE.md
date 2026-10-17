# Bundled component template provenance

The 82 `*.cif` files in this directory are unmodified entries of the wwPDB
[Chemical Component Dictionary](https://www.wwpdb.org/data/ccd) (CCD), one file
per component as served by RCSB at
`https://files.rcsb.org/ligands/view/<ID>.cif`, retrieved on 2026-10-11. They
are compiled into the library as the templates of
`kekule::components::ComponentLibrary::bundled`.

## License

wwPDB archive data, including the CCD, are available under the
[CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) public
domain dedication (see the [wwPDB usage policy](https://www.wwpdb.org/about/usage-policies)).
No attribution is required; this notice records the source.

## Lock

The selection, URL, byte count, SHA-256 and dictionary modification date of
every entry are recorded in
[`tools/components/sources.lock.json`](https://github.com/choutkaj/kekule/blob/main/tools/components/sources.lock.json),
which is not shipped in the package. `tools/components/fetch.py` downloads the
entries and rewrites the lock, and `tools/components/test_provenance.py` checks
offline that the shipped files, the lock and the Rust bundle list agree; CI runs
that check. The selection rationale is in
[`tools/components/README.md`](https://github.com/choutkaj/kekule/blob/main/tools/components/README.md).

The neutral histidine tautomers are derived in code
(`src/io/components/bundled.rs`) and are not dictionary entries.
