# OpenFF SMARTS validation corpus

These are unmodified external files. `sources.lock.json` pins every URL and
SHA-256, the OpenFF force-fields commit, and the OpenFF Toolkit 0.18.0 commit.
Licenses are included separately for each repository. PubChem records are
identified by CID and content hash.

The selection contains complete query rows from Sage 2.2.0, TIP4P-FB 1.0.1,
TIP5P 1.0.0, and the Toolkit's charge-increment and ion-charge fixtures. It covers
ordinary bonded/nonbonded terms, charge increments, library charges, and virtual
sites. All 521 rows are retained, including the malformed upstream
`[N:1](H:2)(H:3)` query in `chargeincrement-test.offxml`. Neither parser accepts
that pattern; it must remain a reported rejection rather than being repaired or
excluded from the corpus.

The four PubChem targets are water (962), methane (297), sodium ion (923),
and chloride (312). They supplement the repository's eight existing PubChem
smoke targets, providing positive hydrogen, ion-charge, and virtual-site matches.
Only SMARTS matching is tested: legacy OFFXML units or schema attributes are not
being validated or imported into a force-field implementation.
