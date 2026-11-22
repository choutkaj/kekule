# Protonation of macromolecular residues

mmCIF and PDB interpretation gives each residue read with a component
template the protonation state that the model and the source support. It
never changes chemistry without a report issue
(`macromolecular::StructureInterpretIssue`). The code is in
[connectivity](../crates/kekule/src/io/macromolecular/connectivity.rs) and
[protonation](../crates/kekule/src/io/macromolecular/connectivity/protonation.rs).

## Candidates

A residue's candidates are its template and the template's protonation
variants. The bundled library derives the neutral histidine tautomers
(`NE2-protonated`, `ND1-protonated`) and the carboxylic acids protonated on
their other oxygen (`ASP` on `OD1` with `HD1`, `GLU` on `OE1` with `HE1`).
The chosen variant's bond orders replace the template's.

Bonds to other residues take the place of hydrogens: a peptide bond replaces
the amine's leaving `H2`, a disulfide cysteine's `HG`. Metal coordination to
another residue takes the place of one hydrogen on a thiol, an acid, or a
nitrogen: zinc-bound cysteine is a thiolate, and histidine bound through
`NE2` carries its ring hydrogen on `ND1`. Alcohols and water keep their
hydrogens. Each such atom is reported (`HydrogenReplacedByCoordination`).

Writers make the same choice for each residue they leave to its template.
The mmCIF writer states the bonds of a residue that readers would protonate
otherwise. PDB states coordination only as `LINK` records from the
interpretation report, so a residue whose variant depends on unstated
coordination fails to write (`DictionaryBondConflict`) unless omitted.

## Residues with hydrogens

Hydrogens decide when the model has them on heteroatoms, or has hydrogens
and the template places none on its heteroatoms. The candidate that fits
best lacks the fewest hydrogens, then changes the fewest charges; ties go
to the template, then the variants in order.

- A hydrogen beyond the template's first fills the place of a template heavy
  atom the residue lacks: `HO5'` of a 5' end without phosphate replaces `P`.
  Any further hydrogen adds +1 on an atom with single bonds, or on a ring
  nitrogen.
- A hydrogen the model lacks deprotonates an acid: a carboxylic, phosphoric,
  or sulfonic acid, a phenol, or a thiol. On a nitrogen, an alcohol, or
  water it is missing from the model: the atom keeps its template charge and
  takes hydrogens by valence (`HydrogensMissing`).
- Changed charges and a chosen variant are reported
  (`ProtonationFromHydrogens`).

Hydrogens with zero occupancy, in a residue whose heavy atoms have
occupancy, are not read (`ZeroOccupancyHydrogens`). A source whose
occupancies are all zero states none and keeps its hydrogens.

Hydrogens no template names bond by name where convention places them, with
`BondEvidence::HydrogenName`: `H1`-`H3`, `HT1`-`HT3`, and `HN1`-`HN3` on a
peptide's `N`; `HO5'` and `H5T` on a nucleotide's `O5'`, `H3T` on its `O3'`.
A hydrogen a variant names bonds as the variant does. Other hydrogens bond
to the nearest heavy atom (`HydrogenAttachedByDistance`).

## Residues without hydrogens

Amino acids and nucleotides read with a dictionary template take the usual
state of their groups at pH 7:

| Group | State |
| --- | --- |
| Carboxylic acid (`ASP`, `GLU`, the C terminus) | Carboxylate, -1 |
| Phosphoric acid (phosphodiester, 5' phosphate) | One -1 per hydroxyl |
| Free amine of the N terminus | Ammonium, +1 |
| Lysine, arginine | As the dictionary: +1 |
| Histidine | Neutral, hydrogen on `NE2` |
| Thiol, phenol (`CYS`, `TYR`), 5' and 3' hydroxyls | Neutral |

One issue per structure counts these residues by component
(`ProtonationAssumed`). Other residues without hydrogens, such as ligands
and force-field residues, whose names state their protonation, take their
template's state. Coordination applies in both cases.

## Atoms the residue lacks

Template heavy atoms a residue lacks are reported (`AtomsMissing`), and
their bonded atoms take hydrogens by valence. A nucleotide without its
whole 5' phosphate is a 5'-hydroxyl terminus and is not reported. A leaving
atom that nothing replaces, such as the `OXT` of a residue at a chain end or
before a gap, or the `OP3` of a 5' phosphate, is reported
(`TerminalAtomMissing`). Its bonded atom takes no hydrogen for it: its
implicit hydrogens are fixed, so a carboxyl terminus without `OXT` stays an
open carbonyl rather than becoming an aldehyde. The mmCIF and PDB writers
write these fixed counts as they are, since readers fix them again.

## Charges

A charge the source states (`_atom_site.pdbx_formal_charge`, PDB columns
79-80) that differs from the derived one is kept when the atom's bonds and
hydrogens can hold it (`StatedChargeKept`): its valence is valid and it needs
no more hydrogens than it has or the state gives it. Otherwise the derived
charge replaces it (`StatedChargeRejected`); a stated +1 on a carboxylate
oxygen would make an oxonium ion. A template charge the dictionary states as
unknown, such as the iron of heme or of iron-sulfur clusters, is read as zero
and reported (`TemplateChargeUnknown`).

A residue without a template takes stated charges, then the charges of the
block's own `_chem_comp_atom` table, then the charges its bonds imply
(`ChargesFromBonds`). An atom with more bonds than its neutral valence takes
+1, and its terminal neighbours that would otherwise take a hydrogen take
-1, as in nitro groups, nitrates, N-oxides, and azides. A residue of one
atom, or one holding a metal, that nothing gives charges is reported
(`TemplateMissing`).
