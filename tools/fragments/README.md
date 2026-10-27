# Bundled builder fragments

`crates/kekule/src/structure/placement/fragments.sdf` holds the fragment
library of `kekule::structure::placement::fragments`: one SDF record per
fragment, the PubChem computed 3D conformer (MMFF94s-optimized, from OpenEye
OMEGA) of the fragment's capped parent compound, with explicit hydrogens.
PubChem data are in the public domain.

`build_fragments.py` writes the file from the records as served by
`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/<name>/SDF?record_type=3d`
(retrieved 2026-10-26):

```
python -I tools/fragments/build_fragments.py <directory of downloaded records> crates/kekule/src/structure/placement/fragments.sdf
```

Each record keeps PubChem's coordinates and bonds, with the atoms reordered:
a ring's atoms first in its conventional numbering (IUPAC locants; the
script's `RING_ORDER` lists PubChem's atom numbers in that order), a group's
attachment atom first and then its other heavy atoms, then each hydrogen in
the order of the atom that carries it. The script
replaces the title with the fragment's name and the comment line with the CID
and parent name, drops PubChem's data fields except `PUBCHEM_COMPOUND_CID`,
and adds:

- `kind`: `ring` or `group`;
- `attachment`: the 1-based atom that bonds when the fragment is attached;
- `leaving`: the 1-based atom that attaching removes with its branch. It is a
  heavy atom only for nitro (nitromethane's carbon), methylsulfonyl and
  sulfonamide (a methyl carbon). Otherwise it is the attachment atom's least
  hindered hydrogen: the one farthest from the atoms not bonded to the
  attachment atom, so a ring attaches through an equatorial position and a
  benzyl's ring stands perpendicular to its new bond.

| Fragment | Parent (CID) | Attachment |
| --- | --- | --- |
| benzene, phenyl | benzene (241) | C |
| pyridine | pyridine (1049) | C4 |
| pyrimidine | pyrimidine (9260) | C2 |
| pyrazine | pyrazine (9261) | C2 |
| pyrrole | pyrrole (8027) | N1 |
| furan | furan (8029) | C2 |
| thiophene | thiophene (8030) | C2 |
| imidazole | imidazole (795) | N1 |
| pyrazole | pyrazole (1048) | N1 |
| oxazole | oxazole (9255) | C2 |
| thiazole | thiazole (9256) | C2 |
| cyclopropane … cyclohexane | the cycloalkanes (6351, 9250, 9253, 8078) | C |
| piperidine | piperidine (8082) | N1 |
| piperazine | piperazine (4837) | N1 |
| morpholine | morpholine (8083) | N4 |
| tetrahydrofuran | tetrahydrofuran (8028) | C2 |
| naphthalene | naphthalene (931) | C2 |
| indole | indole (798) | C3 |
| quinoline | quinoline (7047) | C4 |
| benzimidazole | benzimidazole (5798) | C2 |
| methyl, ethyl, isopropyl, tert-butyl | methane (297), ethane (6324), propane (6334, C2), isobutane (6360, C2) | C |
| vinyl, ethynyl | ethylene (6325), acetylene (6326) | C |
| benzyl | toluene (1140) | methyl C |
| hydroxyl, methoxy | water (962), methanol (887) | O |
| amino, methylamino, dimethylamino | ammonia (222), methylamine (6329), dimethylamine (674) | N |
| formyl, acetyl | formaldehyde (712), acetaldehyde (177) | carbonyl C |
| carboxyl, methoxycarbonyl, carbamoyl | formic acid (284), methyl formate (7865), formamide (713) | carbonyl C |
| acetamido | acetamide (178) | N |
| nitro | nitromethane (6375) | N |
| cyano | hydrogen cyanide (768) | C |
| trifluoromethyl | fluoroform (6373) | C |
| methylsulfonyl | dimethyl sulfone (6213) | S |
| sulfonamide | methanesulfonamide (72879) | S |
| thiol | hydrogen sulfide (402) | S |
