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
another residue takes the place of one hydrogen on an acid, a thiol, an
aromatic ring nitrogen, or a cationic nitrogen: zinc-bound cysteine is a
thiolate, and histidine bound through `NE2` carries its ring hydrogen on
`ND1`. A carboxylate bound through either oxygen keeps its template.
Amines and amides with a neutral template, alcohols, and water keep their
hydrogens: adenine's `N6` bound to magnesium stays `NH2`. Each atom that
gives a hydrogen is reported (`HydrogenReplacedByCoordination`).

Writers make the same choice for each residue they leave to its template,
with the component library readers use (`components` of the write
options). The mmCIF writer states the bonds of a residue that readers would
protonate otherwise. PDB states coordination only as `LINK` records from the
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
- A hydrogen the model lacks deprotonates a strong acid: a carboxylic,
  phosphoric, phosphonic, or sulfonic acid, or a sulfate.
- A cationic nitrogen with at least its neutral count of hydrogens is the
  neutral base: lysine's `NZ` with two hydrogens, or arginine lacking one
  guanidinium hydrogen, whose double bond moves to that nitrogen. One whose
  neutral form has no hydrogen there is the neutral base without any, as an
  imidazolium's `ND1` without `HD1`.
- Any other hydrogen the model lacks is missing from the model: on a phenol,
  a thiol, an alcohol, water, or a nitrogen without any of its hydrogens.
  The atom keeps its template charge and takes hydrogens by valence
  (`HydrogensMissing`). A stated charge or metal coordination deprotonates a
  phenol or a thiol.
- Changed charges and a chosen variant are reported
  (`ProtonationFromHydrogens`), with the atoms left charged.

Hydrogens with zero occupancy, in a residue whose heavy atoms have
occupancy, are not read (`ZeroOccupancyHydrogens`). A source whose
occupancies are all zero states none and keeps its hydrogens. A hydrogen and
a deuterium on one heavy atom share one site, as in neutron structures, when
they lie within 0.3 Å of each other and their occupancies sum to at most
one or they are alternate locations of each other.

Hydrogens no template names bond by name where convention places them, with
`BondEvidence::HydrogenName`: `H1`-`H3`, `HT1`-`HT3`, and `HN1`-`HN3` on a
peptide's `N`; `HO5'` and `H5T` on a nucleotide's `O5'`, `H3T` and `HO3*` on
its `O3'`. A hydrogen a variant names bonds as the variant does. Other
hydrogens bond to the nearest heavy atom (`HydrogenAttachedByDistance`).

## Residues without hydrogens

Residues without hydrogens on their heteroatoms, read with a dictionary
template, take an assumed state:

| Group | State |
| --- | --- |
| Carboxylic, phosphoric, phosphonic, or sulfonic acid, sulfate | One -1 per acidic hydroxyl |
| Free amine of a polymer's true N terminus | Ammonium, +1 |
| Lysine, arginine, and other bases | As the dictionary |
| Histidine in a polymer | Neutral, hydrogen on `NE2` |
| Thiol, phenol, hydroxyl | Neutral |

Acid groups take their pH 7 state in every residue, ligands included.
Polymer residues also take the ammonium of a true N terminus and their least
charged variant. Ligands take nothing else. One issue per structure counts
the residues with an assumed state by component (`ProtonationAssumed`).
Force-field residues, whose names state their protonation, take their
template's state. Coordination then takes the place of a hydrogen the state
leaves.

## Chain ends and gaps

A polymer residue's linking atoms are a peptide's `N` and `C` and a
nucleotide's `P` and `O3'`. A linking atom without a bond to the next
residue is at a gap only when the polymer's sequence has residues past it.
The mmCIF `_entity_poly_seq` and the PDB `SEQRES` records give the
sequence. Modelled residues in another `TER` segment or another molecule
make no gap, and a chain without a sequence has none.

The model's terminal atoms end a chain whatever the sequence says: an extra
hydrogen on `N` (`H1`-`H3`, `HT1`-`HT3`, `HN1`-`HN3`), `OXT` or `HXT`, and
`H5T`, `H3T`, or `HO5'`. Without them, an atom at a gap without hydrogens
keeps the hydrogens it has when bonded: an `N` after a gap is an amide with
one fixed hydrogen, not an ammonium. With hydrogens, its charge follows
them. Each gap is reported (`PolymerLinkMissing`), and a `C` at a gap lacks
no `OXT`.

PDB readers match each chain's modelled residues to its `SEQRES` names,
keyed by the chain letter. `TER` records and jumps in residue numbers split
them into segments, and each segment matches a run of the sequence: by one
residue-number offset where one fits every segment, else at its earliest
run, with a lone first residue at its last place before the next segment
and the last segment at its last run. Each polymer residue takes its place
in the sequence as its label sequence number, and each such chain is a
polymer entity of the report, whose monomers are the `SEQRES` names.

The PDB writer writes `SEQRES` for a chain whose sequence has residues the
model lacks, with the residues' names where the report's entity gives them
and `UNK` elsewhere. It writes sequenced residues as `ATOM` records, also
one a chain break leaves a molecule of its own. Each hierarchy chain
holding amino acids or nucleotides ends with `TER` after its last polymer
residue, so ligands bonded to the chain, such as retinal, are no part of
it. The mmCIF writer writes `_entity_poly_seq`, and numbers a chain without
sequence numbers from 2 when its first linking atom is at a gap.

## Atoms the residue lacks

Template heavy atoms a residue lacks are reported (`AtomsMissing`), and
their bonded atoms take hydrogens by valence. A nucleotide without its
whole 5' phosphate is a 5'-hydroxyl terminus and is not reported. A leaving
atom that nothing replaces at a chain end, such as the `OXT` of a C terminus
or the `OP3` of a 5' phosphate, is reported (`TerminalAtomMissing`). Its
bonded atom takes no hydrogen for it: its implicit hydrogens are fixed, so a
carboxyl terminus without `OXT` stays an open carbonyl rather than becoming
an aldehyde. The mmCIF and PDB writers write these fixed counts as they are,
since readers fix them again.

## Charges

A charge the source states (`_atom_site.pdbx_formal_charge`, PDB columns
79-80) that differs from the derived one applies after the state is chosen.
It is kept only where the atom holds it at a valence its element allows
with that charge:

- When another resonance form of the atom's conjugated or aromatic group
  gives it, the double bond moves to that form before any hydrogen is added,
  and the charges are kept (`StatedChargeResonance`). Examples are -1 on an
  aspartate's `OD1`, +1 on an arginine's `NH1` or `NE`, +1 on a histidine's
  `NE2`, and -1 on a phosphate's `P=O` oxygen.
- A double bond between a stated cation and a stated anion becomes single,
  as in charge-separated DMSO, phosphate, and sulfate
  (`StatedChargeResonance`).
- Otherwise, in a residue without hydrogens, the charge sets the atom's
  hydrogen count within its valence (`StatedChargeKept`): +1 on a
  carboxylate oxygen is kept as `-OH2+`. In a residue with hydrogens, the
  charge is kept only with the hydrogens the model has, and no hydrogen is
  added for it: a stated phenolate lacks no hydrogen.
- Any other charge, such as -2 on a hydroxyl oxygen or +1 on a serine's
  `OG` with its one hydrogen, is replaced by the derived one
  (`StatedChargeRejected`).

A stated 0 in a residue without hydrogens is no evidence, since files
without charges write 0 for every atom. Nonzero stated charges stay
evidence. The PDB writer states every charge, 0 included.

A template charge the dictionary states as unknown, such as the iron of
heme or of iron-sulfur clusters, is read as zero and reported
(`TemplateChargeUnknown`).

A residue without a template takes stated charges, then the charges of the
block's own `_chem_comp_atom` table, then the charges its bonds imply
(`ChargesFromBonds`). An atom with more bonds than its neutral valence takes
the charge that makes its valence valid, +1 before -1, and an onium centre
takes +1 rather than a hydrogen. One singly bonded terminal chalcogen or
boron neighbour of either takes -1 with it, as the oxygen of DMSO drawn with
a single bond and one oxygen of a nitro group. A residue of one
atom, or one holding a metal, that nothing gives charges is reported
(`TemplateMissing`).
