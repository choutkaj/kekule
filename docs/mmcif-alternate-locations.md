# mmCIF alternate locations

`MmcifDocument` preserves atom-site rows, labels, occupancies, coordinate-model
numbers and other source categories. Interpretation selects coordinates without
changing this document. Alternate-location choices and source lines belong to the
interpretation report; canonical `Model` and `Ensemble` keep their existing roles.

## Representative models

The default `MmcifAltLocPolicy::HighestOccupancy` chooses a complete alternative
for each residue. Unlabelled atoms are shared. It does not independently choose
each atom's highest-occupancy row.

The score is the arithmetic mean of labelled-atom occupancies. Shared atoms do not
contribute. Missing occupancies contribute zero and produce a report issue; the
stored occupancy remains missing. Explicitly linked residue groups use the mean
of their residue scores. Equal scores choose the lexically smallest alternative
identifier and produce a tie issue. This avoids favouring alternatives merely
because they have more atoms.

Completeness is relative to the source atom names for the same chemical component,
not a promise that all chemically expected atoms were experimentally observed.
An alternative missing sites present in another alternative of that component is
excluded and reported. Exact selection of an incomplete alternative fails. Atoms
are never borrowed from a different label to complete it. Alternative component
identities at a numbered residue position are mutually exclusive; they can produce
separate Models but may be incompatible with a shared-topology Ensemble.

Source occupancies are site fractions, not probabilities of whole structures.
Independent residue choices do not establish global physical compatibility or a
thermodynamic distribution. Labels shared across distant residues do not by
themselves establish correlation. See the [wwPDB occupancy definition](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_atom_site.occupancy.html)
and [alternate-conformation relationships](https://mmcif.wwpdb.org/docs/tutorials/intro/atom.htm).

## Explicit choices

`SelectLabel("A".into())` requires A at every residue with alternatives, while
retaining shared atoms. `PreferLabel("A".into())` falls back to occupancy ranking
where A cannot be used. `ErrorOnAlternateLocations` rejects even a single labelled
site. For more control, use `Configured(MmcifAltLocSelection { ... })`:

- `default` applies to residues without an override.
- `overrides` maps source residue identities to exact labels.
- `groups` links residues that must use one common label. Overlapping groups are
  united; contradictory exact choices fail.
- `source_conformation` selects a supplied `_atom_sites_alt_gen.ens_id` combination.

`block.alternate_locations()` inventories available labels and residue identities
before selection. Identities include the coordinate model, asymmetry identifier,
sequence namespace and insertion code. Without sequence identifiers, the source
occurrence and component identify the residue. Label identifiers take precedence;
author identifiers are the fallback. A misspelled identity or a group spanning
coordinate models fails instead of silently selecting something else.

When `_atom_sites_alt_gen` supplies combinations, interpretation considers those
combinations across each coordinate model, subject to user overrides and groups.
It does not synthesize additional combinations. A declared combination must cover
every disordered residue in the model to be usable. Its labels may cover several
non-overlapping parts of a residue. Without such metadata or explicit groups,
selection uses residue-sized groups and makes no claim about unknown correlations.

`report.alternate_locations()` records the selected labels, available labels,
selection reason and selected/omitted source lines for the selected coordinate
model. `report.issues()` records ties, missing occupancies, incomplete alternatives
and omissions. All source rows remain validated, even in unselected models.

## Ensembles

`interpret_ensemble` creates one member per deposited coordinate model, applying
the altloc policy within each. A file with one coordinate model produces one
member, whether or not it has alternate locations. No A/B expansion occurs.

`interpret_conformations(&[options_a, options_b])` explicitly interprets the
caller-supplied selections as an unweighted ensemble in that order. It permits
repeated selections and selecting several alternatives of the same source model.
Each selection must produce the same atom identities, chemistry, connectivity and
compatible topology layout. Changed chemistry or atom sets produce an error;
atoms are never padded or discarded to force compatibility. Occupancy is retained
on selected atoms and never converted into member weights.

Dense atom and hierarchy order follow the first source occurrence of each site,
so selecting another alternative does not itself reorder identities. True source
model ordering differences remain subject to the existing ensemble checks.

The API contract lives beside [alternate-location types](../crates/kekule/src/io/mmcif_interpret/altloc.rs)
and [ensemble interpretation](../crates/kekule/src/io/mmcif_interpret/ensemble.rs).
