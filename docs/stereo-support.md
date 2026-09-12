# Stereo support contract

`Molecule` stores explicit stereochemical assertions. Source notation belongs to
format documents; CIP descriptors are derived perception; spatial coordinates
belong to `Model`. Candidate detection does not prove that a center is
stereogenic. CIP assignment does not invent an unasserted configuration.

| Operation | Supported contract |
| --- | --- |
| Represented geometry | Tetrahedral centers, single double bonds, and atropisomeric axes on a single bond |
| CIP assignment | Complete supported assertions using ordered sequence rules and local rooted-digraph auxiliary descriptors; R/S, r/s, E/Z, sequential cis/trans, M/P, and m/p as appropriate |
| Unknown configuration | Preserved as unknown; no specified CIP descriptor is inferred from it |
| Enhanced groups | Molecule-local Absolute, OR, and AND membership; relative membership does not become an absolute configuration |
| Isomeric and canonical SMILES | Tetrahedral and ordinary double-bond stereo, isotope labels, and disconnected component projections |
| V2000/V3000 Model output | Coordinate-consistent tetrahedral, double-bond, and atropisomer stereo; validation uses coordinates rounded to the emitted precision |
| Specified E/Z Model output | Supplied coordinates must encode the assertion; conflicting or degenerate drawings fail explicitly |
| Molfile output of an unasserted alkene | Crossed/either syntax prevents a nondegenerate drawing from inventing specified E/Z; rereading may introduce an explicit unknown element |
| Geometry-free Molfile output | Specified E/Z requires a Model; the writer does not synthesize a layout |
| V3000 atom CFG input | Tetrahedral parity 1/2 and unknown 3; CTfile carrier order puts hydrogen last, including explicit hydrogen rows |
| V3000 enhanced atropisomer groups | Axis endpoint `ATOMS` convention, with unambiguous membership; tetrahedral and axial members may share a group |
| Hydrogen collapse | Carrier remapping preserves valid stereo assertions, IDs, and groups |

V2000 cannot preserve enhanced groups and rejects them. Automatic Molfile version
selection promotes group-bearing structures to V3000. V3000 collections cannot
encode double-bond group members or ambiguous endpoint membership. Absolute groups
may be partitioned across connected components. Cross-component OR/AND relations
are rejected because molecule-local ownership cannot represent their relationship.

Plain isomeric and canonical SMILES reject axial stereo, explicit unknown configurations, and
enhanced groups. A successful ordinary isomeric SMILES write is therefore not a
claim that every represented geometry has a SMILES encoding. Chiral SMARTS is
unsupported. Canonical output preserves supported stereo and isotopes and is
invariant under atom numbering. Ordinary SMILES rejects stereo assertions.
Directional bonds must reconstruct exactly the asserted alkene configurations;
unrepresentable partial assignments fail rather than inventing stereo.

Molfile output preserves unknown tetrahedral and double-bond configurations, but
rejects unknown axes. Molfile input can retain an axis declared unknown by a wavy
member. Source support and export support are distinct contracts.

Allene and longer cumulene axes, extended cumulene E/Z, helical and planar
chirality, and square-planar/trigonal-bipyramidal/octahedral coordination stereo
need additional represented geometry types. Unsupported source assertions reject
the whole input when they cannot be preserved, including other supported centers
in the same input. Wildcard/query atoms are a separate syntax/model boundary.

Assignment limits are explicit. The default expansion depth is 32 and the node
bound is 100,000 per ligand or auxiliary graph. Depth and node exhaustion are
errors, never ties or evidence of nonstereogenicity. Failed assignment leaves the
previous installed perception intact. Callers requiring different bounds must
select them explicitly and handle failure.

These contracts are enforced by focused regressions and format round trips. The
optional [external validation procedure](stereo-validation.md) compares complete
outputs and retains unsupported inputs, unavailable reference results, and known
differences. Finite testing does not establish universal chemical completeness.
