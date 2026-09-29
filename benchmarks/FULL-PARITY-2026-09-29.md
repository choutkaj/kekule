# Full benchmark parity audit, 2026-09-29

## Scope and evidence

Baseline: completed `benchmarks/runs/run-1790612276267-32424.json`, recorded
revision `e554f56b9afa2700f62a68327b31bd51c544fe8e`, clean working tree. This
review started from `764961b3` on `codex/full-benchmark-parity`. The intervening
main changes concern dashboard generation, not scientific algorithms.

The streaming audit inspected all 6,514,328 case records: 6,170,041 applicable
feature/corpus cases and 344,287 not-applicable records, covering 29 features
and seven corpora. Every failure bucket reconciles with the original report:
5,322,328 agreements, 799,470 disagreements and 48,243 errors. These are
feature/input observations, not counts of distinct molecules or defects.

Local evidence and scripts are in `target/full-parity-audit/`: `summary.json`
groups every error and difference path, `failures.jsonl` retains failure
identities, and `examples.jsonl` retains complete examples of differing field
patterns. Counts of difference paths overlap; a path is a diagnostic grouping,
not by itself proof that every case has the same cause.

No golden, source lock, comparison tolerance, asserted field or reference adapter
was changed. Earlier parity and general-chemistry reviews are recoverable from
Git history; several of their deferred issues were fixed before this baseline.

## Changes

### Four-coordinate wedge drawings

The source decoder used the volume of a displaced tetrahedron to interpret a
flat four-neighbor drawing. In crowded projections the result depended on drawn
bond lengths and could invert the represented configuration. The decoder now
uses angular ordering of normalized directions and the alternate projected
volume when the first pair is degenerate or represents the opposite carrier.
This follows the drawing conventions in RDKit's pinned
[chirality implementation](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/Chirality.cpp).

True 3D geometry and virtual-carrier interpretation retain their existing paths;
ambiguous marks remain explicit warnings. Regressions use two local geometries
from supplied PubChem 10524 with distinct test ligands, including reflection,
scale, V2000/V3000 and carrier permutations. The new drawing regression failed
before the fix. All 70 formerly opposite tetrahedral CIP observations now agree;
there are no newly failing CIP cases in the complete paired corpus audit.

### Source-order SMILES emission

Source-order modes now give acyclic branches priority over ring continuations,
use bond order and source atom order within those classes, and keep the last
child as the continuation. Ring closures follow DFS discovery, emit their bond
symbol at the closing endpoint, and normally delay reuse of a closed label until
the next atom. These conventions follow the pinned
[RDKit traversal](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/Canon.cpp).
Canonical and source-order modes share the traversal planner; canonical ranking
is unchanged. Kekule-form output also omits unnecessary single-bond hyphens
between uppercase atoms. Aromatic output still requires `-` between aromatic
atoms joined by a nonaromatic single bond.

The existing 100-live-ring-label regression exposed a capacity loss from delayed
reuse. Source-order export now permits immediate reuse when every other label
is occupied; the full capacity and explicit exhaustion contract are preserved.
This fallback is only reached where the preceding implementation returned a
resource error; none of the successful corpus emissions used it.

Ordinary writing still rejects represented stereochemistry. It does not discard
stereo to eliminate the 20,904 intentional stereo-rejection observations.
Isomeric output retains the stored localized bond assignment, so alternative
Kekule forms can still differ textually from RDKit without changing chemistry.

### V2000 coordinate and atom-field parsing

`BMS-986142_3d_chiral.mol` contains negative z coordinates such as
`-15355.5894`, occupying the otherwise separating column after the ten-character
coordinate field. The parser silently dropped the last digit. It now consumes
that column instead of truncating the value. It also tokenizes atom attributes
after the coordinate block: touching negative coordinates can no longer shift
charge, H-count, valence and atom-map fields. Focused regressions assert all of
those fields and coordinate precision. The supplied structure's coordinate
disagreement disappears without loosening the numerical comparison.

## Measured results

Full-corpus reruns use the original reference archives. Paired comparisons check
case counts, golden archive hashes and individual agreement/error transitions.

| Feature and scope | Agreements before | Agreements after |
| --- | ---: | ---: |
| Ordinary SMILES exact text, all supplied SMILES | 29,991 | 128,275 |
| Isomeric SMILES exact text, all supplied SMILES | 12,800 | 104,789 |
| Canonical SMILES exact text, all supplied SMILES | 144,677 | 144,677 |
| PubChem MOL parsing | 94,031 | 94,089 |
| RDKit structures MOL parsing | 11 | 12 |
| PubChem CIP, both input formats | 199,762 | 199,832 |

The exact-text gains are monotonic: 98,284 ordinary and 91,989 isomeric cases
change from disagreement to agreement, with no newly failing exact-text cases.
Ordinary Enamine output agrees for all 41,874 supported records; PubChem retains
1,060 ordinary-text disagreements. Isomeric text retains 45,450 disagreements;
canonical text retains 5,562. Unsupported inputs and errors remain in the counts.

Independent RDKit rereading of every successful final SMILES emission uses the
existing `read_written` adapter and original identity goldens:

- Ordinary: all 129,335 successful outputs agree on identity.
- Isomeric: 150,234 agree; five Enamine canonical identity strings differ.
- Canonical: 150,175 agree; the existing 64 Enamine identity differences remain.

The five new isomeric identity-string differences require explicit qualification.
They are `Z9237738672`, `Z5188334346`, `Z3510155828`, `Z4875774646` and
`Z3510152401`, all enhanced spiro-stereo cases. All five preserve native
canonical stereo on rereading and pass mutual RDKit substructure matching with
both chirality and enhanced stereo enabled. Three emissions are byte-identical
to RDKit's own independent isomeric output. For four cases, rereading RDKit's
own emitted text produces the same identity discrepancy. The fifth is the
previously documented enhanced-spiro canonicalization instability. Thus this
is an identity-string regression, but the independent stereo-aware comparison
does not indicate lost stereochemistry. The failures remain counted; no
case-specific spelling or relaxed identity comparison was introduced.

Evidence: `spiro-audit.json`, `identity-audit.json`, and a native regression
covering all five supplied sources. Of the 64 existing canonical identity
differences, 62 involve legacy `r` semantics. The other two, `Z9082652630` and
`Z9237738672`, are enhanced-spiro cases. For `Z9082652630`, the emitted text
exactly matches the stored RDKit text, but RDKit rereading changes its canonical
stereo spelling. Legacy relative relationships remain distinct in Kekule;
changing their meaning to match RDKit's chosen absolute representative would
not be a spelling-only fix. The full earlier SMILES follow-up remains in Git
history.

## Second round: shared axis endpoint geometry

The next round replaces two inconsistent endpoint heuristics with one shared,
conservative trigonal-geometry predicate. Source normalization previously used
ring membership or any double bond; coordinate inference used aromaticity or
any double bond. Both could mistake a pyramidal sulfoxide for an axis endpoint.
Coordinate inference could consequently propose both tetrahedral and axis
stereo for the same sulfur environment.

The shared rule counts three explicit ligand directions and the nonbonding
electrons implied by localized bond orders and formal charge. It distinguishes
localized lone pairs from second-row lone pairs adjacent to a pi system, so
amide/pyrrole nitrogen remains eligible while saturated rings, phosphines,
sulfoxides, selenoxides and sulfonium centers do not become trigonal merely
because they touch a ring or a multiple bond. This is a conservative local
atropisomer subset, not a new public or universal hybridization model.

Source normalization supplies only declared hydrogen counts and does not read
or install perception; coordinate inference supplies total implicit counts.
Structural validation of explicitly represented axes remains separate from
inference eligibility. No molecular owner, stereo family or public API changes.

Source-drawing and coordinate-inference regressions both failed before this
change. They now cover sulfur and selenium, opposite wedge configurations,
V2000/V3000 writer round trips, tetrahedral lone-pair carriers, and absence of
spurious axes. Additional local eligibility tests cover conjugated nitrogen,
electron-deficient trigonal centers, saturated rings and pyramidal donors.

Full-corpus CIP rerunning raises PubChem agreement from 199,832 to 199,835,
resolving supplied 146091, 461502 and 461520. The remaining CIP disagreements
are 145, with errors unchanged. PubChem MOL agreement rises from 94,089 to
94,092; the other corpora retain their first-round totals. All 29 features in
`rdkit-structures` have identical observations, including the existing
atropisomer coverage. Detailed second-round evidence and validation logs are
in `target/parity-round2/`.

The completed paired audit (`paired.json`) verifies identical case counts and
golden hashes. Exactly those three MOL observations and three CIP observations
change; each changes from disagreement to agreement. Every other actual
observation in these reruns is unchanged, with no new errors or disagreements.

All commands listed under Validation and reproduction below were repeated after
the second-round implementation and passed, including Rust 1.89, full workspace
tests, warnings-as-errors clippy/documentation, and package validation. The
second-round command results are in `target/parity-round2/checks.json` and
`msrv.log`. Scientific reruns cover MOL parsing and CIP over every corpus, plus
all features over `rdkit-structures`; the complete all-feature/all-corpus run
was not repeated. The same platform and companion-package exclusions below
apply. Existing first-round fixes and their tests remain included.

## Remaining issues and decisions

| Area | Current evidence and decision |
| --- | --- |
| Hydrogen representation | SMILES parsing differs in 6,520 PubChem and 10,923 Enamine cases, dominated by declared versus inferred H storage. MOL/SDF additionally differ in inference permission on 8,241 Enamine records. Preserve native declarations and explicit inference policy; do not normalize observation fields just to match. Enamine group differences overlap these counts. |
| Source stereo and cleanup | Redundant/unknown source assertions remain visible. Parsing and default perception do not perform the separate, explicit symmetry-cleanup operation. The second round fixes the three supplied sulfoxides through shared axis endpoint geometry, as described above. |
| CTfile radicals | Five PubChem MOL inputs differ: 24755, 34052, 177071, 498830 and 499023. The differing carbon atoms declare valence 2 without explicit radical declarations on those atoms; RDKit infers two radical electrons. Reading the supplied sources before and after RDKit sanitization confirms that inference. Resolve the format's radical inference policy consistently with writer round trips before changing represented radical state. |
| Potential stereochemistry | Residuals include RDKit candidates on cumulenic N=C/N=N endpoints and asserted atropisomer axes absent from the native candidate family. For example, PubChem 164 has two double bonds at its nitrogen endpoint. Expanding the represented/candidate geometry requires coordinated validation, inference and matching, rather than treating every RDKit candidate as a supported isolated alkene. |
| CIP | After the second round, PubChem has 145 disagreements involving bond descriptors and source-stereo issues. Two observations of PubChem 158374 exhaust the 100,000-node bound. Keep exhaustion as an error; increasing a cap is not an algorithmic correction. |
| Rotatable bonds | 36,810 baseline disagreements remain. Native strict mode explicitly ignores graph H in terminal/symmetric-group decisions; RDKit's degree-based query can depend on explicit H representation. The old local-aromaticity approximation has already been replaced by fallible default perception. Do not reintroduce representation dependence to improve this raw score. |
| Molecular masses | All 300,821 disagreements are in mass fields, not formula/composition. Native CIAAW/AME constants and consistent charge/electron-mass conventions differ from RDKit. Thirty-six additional PubChem descriptor failures report unavailable standard weights. Preserve supported modern data and explicit absence rather than substituting isotope mass numbers as standard weights. |
| Ring basis and mmCIF | One alternative SSSR basis and 79 multiline-whitespace differences remain. Preserve asserted cycle paths and decoded values; no normalization or field removal was added to the comparator. |
| DSSP | 194 disagreements and 250 error cases remain. Reference conformer selection can differ from native coherent whole-residue selection; chain-boundary and nonpolymer eligibility policies also differ. Native errors include 234 inputs without analyzable residues, four unknown-element X inputs, and two inputs without a complete coherent alternate configuration. Reference failures also include missing metadata and mkdssp process errors. An algorithm comparison needs matched selected coordinates before assigning all end-to-end differences to DSSP itself. The previously identified weak-hydrogen-bond case still warrants a focused geometry audit. |
| Unsupported structures/formats | Three RDKit structures contain query symbols R/R1 outside a concrete Molecule. `BMS-986142_atropBad2` asserts tetrahedral CFG on an unsupported nitrogen carrier geometry. V2000 export rejects enhanced groups; V3000 is the supported lossless route. Preserve these errors rather than discarding assertions or adding query atoms to the molecular graph. |
| Shared failures | Nine invalid-valence PubChem molecules and one SMARTS resource-limit case fail in both implementations. Shared failures remain errors, not agreements. |

The most useful remaining algorithm work is candidate-family coverage, canonical tie-breaking
and a DSSP comparison on identical conformers. This audit does not claim those
broader problems are solved or justify changing their contracts by isolated
benchmark-specific patches.

## Validation and reproduction

Regression tests cover wedge scale/reflection/permutations, source-order branch
and ring traversal, aromatic versus Kekule bond symbols, enhanced-spiro native
round trips, all 100 ring labels and touching/wide V2000 coordinates with atom
attributes. Existing graph/stereo assertions remain enabled.

All final commands below passed. Results are recorded in
`target/full-parity-audit/checks.json`, with the separate Rust 1.89 compatibility
check in `target/full-parity-audit/msrv.log`:

```text
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features --locked --offline
cargo +1.89.0 check --workspace --all-targets --all-features --locked --offline --target-dir target/altloc-msrv
cargo clippy --workspace --all-targets --all-features --locked --offline -- -D warnings
cargo test --workspace --all-features --locked --offline
cargo test --workspace --all-features --doc --locked --offline
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --locked --offline
cargo package -p kekule --locked --offline --allow-dirty
cargo package -p PACKAGE --locked --offline --allow-dirty --list
```

The last command covers `kekule-potentials`, `kekule-traj`, `kekule-openff` and
`kekule-bench`. Companion full package builds follow the repository CI exclusion
until the foundational dependency is published. Linux fuzz compilation/smoke
runs were not run on this Windows host. Unchanged Python/dashboard generators
were not retested. No README was modified.

Full scientific reruns: `io.mol.parse`, `stereo.cip`, and all three exact SMILES
text modes over `--dataset all`; all 29 features over `rdkit-structures`.
Identity rereading is an additional independent audit of emitted text, not a
replacement for the exact-text metric. The complete all-feature/all-corpus
benchmark was not repeated; unaffected feature conclusions use the original
complete run and the exhaustive failure audit. Reports correctly exit nonzero
for retained disagreements and errors.
