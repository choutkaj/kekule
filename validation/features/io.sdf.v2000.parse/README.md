# io.sdf.v2000.parse Validation Fixtures

These fixtures are manual parser inputs, not RDKit-generated goldens.

The valid set covers:

- simple atom and bond blocks
- multi-record SDF streams
- data fields
- aromatic bond order code `4`
- V2000 metadata records such as `M  CHG` and `M  ISO`

The invalid set is intended for negative parser checks. Some valid fixtures include constructs the
prototype may not fully preserve yet, especially charge and isotope metadata. Keep those fixtures as
long-term validation pressure rather than pruning them to current behavior.
