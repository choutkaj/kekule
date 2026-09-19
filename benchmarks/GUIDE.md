# Correctness benchmarks

The benchmark measures agreement with independent references on externally
supplied, hash-locked inputs. Generation runs only the reference software;
comparison runs Kekule against stored observations. Matching failures never
establish correctness. Scientific runs are optional, outside CI/release gates.

## Run

```text
cargo benchmark --list
cargo benchmark --feature io.smiles.parse --dataset smoke
cargo benchmark --feature io.sdf.parse --dataset enamine-diversity
cargo benchmark --feature all --dataset smoke --writer-python PATH_TO_RDKIT_PYTHON
cargo benchmark generate --feature FEATURE --dataset DATASET --python PATH --goldens NEW_DIRECTORY
```

A fresh checkout includes the smoke inputs and goldens. Full dataset inputs and
golden archives stay local and are excluded from Git and Cargo packages; no
download or upload is automatic. Their source locks and golden manifests remain
tracked. To run a full dataset on another machine, supply the original inputs
and matching golden archives locally at the paths below. Comparison verifies
their hashes against the tracked provenance. Missing files stop the run.

The default selects all locked source IDs. `--limit N` selects a deterministic
subset before either engine runs; molecular subsets are nested and PDB subsets
retain lock-file order. `--jobs N` controls Kekule concurrency. Reports and per-case
JSONL are saved under `benchmarks/runs/`, with unique run names. `--output FILE`
selects another report location; its summary is also archived in `benchmarks/runs/`
without duplicating its case records. Existing reports are never overwritten.

Normal parser and algorithm checks need Rust and the stored goldens. Writers
also need RDKit to read the emitted text. Select it with `--writer-python`,
`KEKULE_WRITER_PYTHON`, or the activated environment's `python`.
Before evaluating the first applicable writer input in each feature/dataset,
the runner checks that this interpreter can load the reference and matches the
stored tool/version. A setup failure stops with an incomplete report; it does
not turn the entire corpus into apparent writer failures.
Create the pinned reference environment from `reference/environment.yml`;
engine-specific environment files are available for narrower installations.
DSSP requires mkdssp, its shared libraries, and its CIF dictionary on the
activated environment's paths.

Generation is explicit and never evaluates Kekule. Write regenerated goldens to
a new directory, investigate changes, and review them before adopting them.
Missing, stale, malformed or duplicated goldens stop comparison; they do not
trigger regeneration. A partial generation does not cover omitted cases.

## Dashboard

Every comparison refreshes `benchmarks/runs/index.html`, including comparisons that
fail or stop with an error. Open that page once: it checks the adjacent local data
file every five seconds and selects the newest run when history changes. Unchanged
history preserves your selected run and search. Reports are ordered by their recorded
start time, not file modification time. Reference generation is excluded.

Python 3.11+ and its standard library are needed to render the dashboard. The command
uses `KEKULE_DASHBOARD_PYTHON`, an executable path in the ignored file
`benchmarks/.dashboard-python`, the selected writer Python, or Python on `PATH`.
If rendering fails, a warning identifies the problem; reports and the scientific
exit status are preserved. `KEKULE_BENCHMARK_RUNS_DIR` can override the local history
directory. There is no server or upload.

To refresh history manually, including checkpoints left by an aborted process:

```text
python benchmarks/dashboard.py
```

Invalid or incompatible reports are skipped with a warning; identical reports are
deduplicated. Incomplete comparison reports remain labeled in provenance. Original
reports are retained. Old runner filenames carry a sortable start time; other
untimestamped imports sort after timestamped runs, without inferring dates from mtime.

You can also build a fixed snapshot from explicit schema-2 comparison reports:

```text
python benchmarks/dashboard.py target/smoke.json target/full-corpus-sample.json
```

Open `target/benchmark-dashboard/index.html` in a browser, or choose another
destination with `--output PATH.html`. Python 3.11+ and its standard library are
sufficient. This explicit-report mode produces a single file without live polling,
a server, network requests or external assets. It can later be hosted unchanged.
Rerun the command to refresh this derived page; input reports are never modified.

The dashboard shows a searchable feature/dataset comparison table and an overview
of the data corpora: source ID counts (N), supplied formats
and selection notes. Plots group exact and within-precision agreement into one
"Agrees" category; report counts and comparison tolerances are unchanged. Runs
remain separate. Cells distinguish full selections, deterministic samples, stale
provenance and absent results; run completion is recorded in the provenance panel.
No overall correctness score or cross-library speed ratio is computed. Counts
repeat inputs across features and formats; unavailable formats stay outside the
measured denominator.

Only aggregate counts, source membership sizes and hashes/tool versions are
embedded. Inputs, golden payloads, per-case observations, local paths and raw
error messages stay local. The page can be generated without the full datasets
or golden archives. Run history and local Python configuration are excluded from
Git and Cargo packages; fixed snapshots under `target/` are likewise local.
Publication is a separate, deliberate step; this command uploads nothing.
"Download plotted data" saves a sanitized dashboard export, not the original
runner report. Rebuild the page from the original comparison summaries.

```text
python -m unittest discover -s benchmarks -p test_dashboard.py
node benchmarks/test_dashboard.cjs
```

## Observations and comparisons

The versioned contract is [contract.json](contract.json). Successful observations
must satisfy the feature's schema: required fields, types and collection presence
are checked, and unknown fields are rejected. Nonfinite measurements fail before
serialization. Nullable measurements distinguish an unavailable value from a
missing required field. Coordinates occur only when a conformer exists.

Every structural difference and numerical difference is retained, with its path
and both raw values. Discrete quantities compare exactly. For floating-point
values, the allowance is `16 * f64::EPSILON * max(abs(expected), abs(actual))`,
covering the fixed decimal parsing/unit conversion chain. Additional allowances
follow output precision, rather than an empirical fit to results:

| Observation | Additional absolute allowance |
| --- | ---: |
| DSSP legacy energies, phi/psi, reported alpha/kappa | 0.05 kcal/mol or degrees |
| DSSP TCO | 0.0005 |
| DSSP omega (single-precision geometry versus Biopython vectors) | `16 * f32::EPSILON * 180` degrees |
| V2000 writer coordinates | 0.00005 angstrom |
| All other numerical observations | None |

Angles use circular distance. Values are never rounded or replaced. Each case
reports both exact equality and agreement within the declared precision; mass
constant differences and other larger numerical differences remain failures.
Biopython calculates omega with double-precision vectors made from parsed
single-precision coordinates; DSSP 4 geometry uses single-precision arithmetic.
The omega allowance accounts for arithmetic roundoff on the angular domain,
including angles near zero. It is not an allowance for different conformers;
larger deviations and every raw value remain visible.

Only representation order is normalized: undirected bond endpoints/list order,
cyclic ring paths (starting atom, direction and ring-list order), and arbitrary
DSSP sheet/strand/ladder labels and residue order. Dative direction, atom
correspondence, stereo parity and group membership
remain significant. Complete indexed molecular graphs include atoms, bonds,
charges, isotopes, radicals, maps, declared and inferred hydrogens, represented
valence, aromaticity, stereo and relation groups. Both engines apply their
ordinary perception/sanitization before observing this state. SDF properties
remain ordered, including source names beginning with an underscore; RDKit's
internal properties are excluded. Duplicate names that RDKit cannot preserve
are reference errors.

Fast ring membership compares every atom's cyclic flag and every indexed bond's
cyclic flag. Its reference invokes RDKit `FastFindRings`, without computing a
selected ring set; `algo.rings.sssr` separately invokes `GetSymmSSSR`. Zero-order
and dative bonds do not close cycles in these default policies. Multiple input
components remain separate observations. A matching number of ring atoms alone
does not establish agreement.

Selected rings retain their cyclic atom order: different paths through the same
atoms remain different rings. This compares Kekule's Figueras-style selection
with RDKit's symmetrized selection, not a unique mathematical minimum cycle
basis. Candidate traversal can depend on bond iteration order; such differences
remain reported rather than being replaced by counts or cycle-space equivalence.

RDKit-like valence compares the atoms after represented-chemistry normalization.
The reference runs RDKit `Cleanup` on a copy before the non-strict property-cache
update; Kekule normalizes its represented graph during molecule publication.
This includes the charge-separated representation of oxohalogens. The reference
uses RDKit's complete cleanup operation, so any additional normalization
differences remain observable. This feature does not run full sanitization,
aromaticity or radical perception. Formal charge, explicit valence, explicit
hydrogen declarations and inferred implicit hydrogens all remain asserted.

The aromaticity feature compares every indexed atom and bond flag after each
engine's default preparation: Kekule's default perception and RDKit's full
sanitization. It therefore includes upstream radical-inference differences and
strict-valence failures, rather than measuring aromaticity on an externally
equalized graph. A bond between aromatic atoms is not necessarily aromatic;
atom and bond flags are asserted independently. Upstream differences remain
visible and must be distinguished from aromaticity-algorithm defects.

`chem.perception.default` compares the full graph and per-atom valence after
each engine's ordinary preparation. Kekule installs derived valence, rings and
aromaticity without rewriting represented chemistry; the RDKit reference runs
full sanitization and stereo cleanup on a copy. Thus this is an end-to-end
prepared-state comparison, not an isolated test on identical represented input
graphs. Source stereo, hydrogen-storage policy and radical-inference differences
remain visible alongside derived-chemistry differences. Coordinates and SDF
properties are outside this feature; component-local atom correspondence and
every graph/valence field remain asserted. Selected ring lists are checked by
the separate ring feature, not included in this observation.

Hydrogen transformations compare the expanded graph, the number of new H atoms
attached to each original parent, and the collapsed graph after addition followed
by removal. Existing graph hydrogens are not counted as newly added. Removal
acts on all eligible hydrogens, so `round_trip` names the collapsed result and
does not imply recovery of the original graph or atom numbering. Source stereo,
radical inference and hydrogen-policy differences remain observable. The native
remover preserves hydrogen counts and explicit information conservatively;
RDKit uses its default `RemoveHs` policy, including retaining some hydrogens that
define double-bond stereo. Reference defaults can also lose metal-bound hydrogen
counts when the metal has no inferred replacement valence. Such reference
limitations remain disagreements with their raw observations intact; they are
not a reason to discard hydrogen information in Kekule.

Canonical ranking compares partitions of source atom indices, rather than the
engines' arbitrary numeric rank labels. The reference disables tie breaking and
chirality while retaining isotope and atom-map distinctions. Both paths prepare
the molecule first; hydrogen counts contribute through their current total,
independently of fixed versus inferred storage. Every class and every member
remains asserted. This feature does not test canonical traversal strings or
stereochemical equivalence; those have separate features.

The SDF parsing feature uses each engine's SDF reader for every selected
structural input, including standalone MOL files. Kekule explicitly permits
EOF termination for `.mol`/`.mdl` inputs; `.sdf` inputs retain its default
delimiter checks. Records and source fields come from the same SDF parse and
interpretation, rather than a separate metadata pass.

In graph observations, aromatic single/double bonds use the reference's
`AROMATIC` type. Higher orders retain their type alongside the aromatic flag:
an aryne bond can be both `TRIPLE` and aromatic. Aromaticity does not erase its
additional bond order.

These are native API observations, not a pure chemical-equivalence score.
RDKit sanitization can move inferred aromatic hydrogens into its explicit count
and assign radicals from electron deficits. Kekule preserves source hydrogen
declarations and does not infer radicals from bracket notation. Differences in
these fields remain visible; matching total hydrogen counts alone does not make
the complete observations agree. Interpret these distinctions before attributing
every parsing disagreement to incorrect chemistry.

SMILES writers compare complete canonical isomeric CXSMILES identity, without
requiring the same traversal string as RDKit. Canonical mode additionally checks
a read/write fixed point and invariance under reversing atom and bond numbering.
These probes do not exhaust all permutations. CXSMILES source extensions are
parsed by RDKit; unsupported extensions cause an explicit Kekule error.

Versioned MOL/SDF writers must emit the requested format. Each MOL file must end
after its single complete record; extra molecules or SDF data are rejected even
when RDKit would silently ignore them. RDKit reads the emitted
bytes independently and rejects query atoms or bonds in molecular output.
For example, V3000 `HCOUNT` is a query constraint, not a molecular fixed-H
declaration; matching the decoded hydrogen count alone is insufficient.
V3000 coordinates retain floating-point precision, without the V2000
four-decimal quantization allowance. The MOL model API has no title/property input, so its
expected serialized title is empty and its property list is empty; the original
source observation remains intact in the report. SDF document writers retain
titles and properties. Coordinates and all chemical fields remain asserted.

SDF record separators must occupy complete `$$$$` lines. Text resembling a
separator or field header inside a title or field value remains data. Output
validation checks every record and rejects missing terminators or unexpected
content outside the CTAB and data fields, even if RDKit would ignore it.

Valence and descriptors read each library's direct results. No adapter adjusts
aromatic nitrogen, masses, formal charges, or hydrogen-removal policy to force
agreement. Molecular descriptors assert the complete isotope-resolved formula,
formal charge, average mass and monoisotopic mass of each prepared component.
Kekule uses CIAAW 2024 abridged standard weights, CIAAW 2024 natural isotope
abundances, AME2020 isotope masses and the CODATA 2022 electron mass. It subtracts
the charge's electron mass in both mass methods. RDKit's direct `MolWt` uses its
own atomic weights without that charge correction; `ExactMolWt` uses its own
isotope masses and electron correction. These different constants and conventions
remain numerical disagreements, even when the complete formulas agree. Abridged
standard weights are conventional point estimates, not exact sample masses.
Missing standard weights or natural isotope abundances cause explicit native
errors; no representative isotope's mass number is substituted as a weight.

Strict rotatable bonds compare every component-local bond endpoint pair as well
as the count. The reference enumerates RDKit's strict descriptor SMARTS without
a match cap and checks that its endpoint count equals `CalcNumRotatableBonds`
with `Strict` selected explicitly. It retains source graph hydrogens: the
query's degree and hydrogen-count predicates can therefore change the result
when hydrogens are materialized. Kekule deliberately uses heavy-atom degree,
so its terminal-group classification is invariant to that representation.
These policy differences remain visible. Kekule's localized five-/six-member
cycle approximation for resonance classification can also differ from RDKit's
aromaticity-based query in unusual or fused systems. This descriptor estimates
rotatable axes; it does not calculate rotational energy barriers.

`stereo.representation` compares the complete graph after ordinary source
interpretation and preparation, including coordinates, source properties and
enhanced stereo groups. It is not a stereo-only score. Tetrahedral carrier order
and parity are normalized together; double-bond and axis orientations use
canonical endpoint carriers. Focus identity, unknown versus specified
orientation and group membership remain significant. RDKit applies its default
stereo cleanup and interprets Molfile coordinates according to its source-reader
policy. Kekule preserves represented assertions and decodes source wedge and
double-bond drawing marks without general coordinate-stereo materialization.
These distinctions, including symmetry-dependent cleanup and supported stereo
families, remain visible. The separate `stereo.perception` feature explicitly
materializes coordinate stereo and compares candidate detection as well.

`stereo.perception` retains that complete-graph comparison. Kekule reports local
tetrahedral and double-bond candidates, excludes repeated hydrogen ligands and
overcoordinated double-bond endpoints, and preserves existing represented stereo
when materializing coordinate proposals. RDKit's `FindPotentialStereo` also
resolves ligand equivalence and stereo dependencies and can report other stereo
families. Its preceding `AssignStereochemistryFrom3D` call uses the default policy
for replacing existing tags on 3D conformers. The comparison keeps these
differences visible; candidate agreement is not inferred from a count alone.

`stereo.cip` compares atom/bond counts and descriptors assigned to represented
stereo after ordinary source preparation. It does not run the separate coordinate
materialization workflow. Descriptor identity, focus identity and lowercase
pseudoasymmetric labels remain significant. The adapter renders native
`SeqCis`/`SeqTrans` as RDKit's lowercase `z`/`e`; ordinary `Z`/`E` remain distinct.
Native assignment uses the public default bounds (depth 32 and 100,000 nodes),
while RDKit `AssignCIPLabels` is bounded at 1,000,000 recursive iterations.
These are different algorithmic resource units, not equivalent work budgets;
exhaustion remains an explicit error and is included in the reported coverage.

Substructure retains every query-to-target mapping without a match
cap. The shared 18-query input is [queries.smarts](queries.smarts).
Mappings retain query-atom order; only the list of complete mappings is sorted.
Query automorphisms are included, and resource exhaustion remains an error
rather than a partial successful match list. SMARTS single/double bond types
exclude aromatic bonds, while aromatic bond types exclude higher represented
orders even when those bonds carry an aromaticity flag. The low-level query
predicates for localized order and aromatic membership retain their independent
meanings; the SMARTS parser combines them to express these bond types.
`query.smarts` reads each text record's first whitespace-delimited token as
SMARTS, with the remainder as its title. The current external corpus supplies
SMILES strings to exercise the overlapping grammar; this does not validate
molecular interpretation or CXSMILES extensions. The comparison checks syntax
acceptance and atom/bond counts, not predicate equivalence. The behavioral
searches above supplement this limited coverage. Unsupported stereochemical
queries remain reported errors; they are neither stripped nor excluded.

mmCIF compares all decoded tags/values, including non-atom categories and distinct
`.`/`?` tokens. It does not claim full biomolecular topology interpretation
coverage or a complete CIF syntax conformance suite. Biopython's `MMCIF2Dict`
also loses the distinction between quoted literal `.`/`?` and bare missing-value
tokens, and is permissive about malformed source syntax; those semantics need
focused runtime tests. Biopython's multi-block limitation remains a reference
limitation. Exact value comparison deliberately exposes multiline trailing-space
differences: CIF 1.1 permits their removal, Biopython strips them, and Kekule
preserves them. This is a formatting-policy difference, not lost chemical data.
DSSP uses
the first model and original source bytes, including archive metadata. It
compares secondary structure, backbone measurements, hydrogen-bond energies and
partner identities, helix/sheet/ladder information, explicit beta partners, and
CA-C-N-CA omega. Alternate-location/model-policy differences remain visible.
The polymer sequence label is required as an observation but may be null for
nonpolymer residues; such residues remain in the structural comparison.
Kekule analyzes the model's selected conformer, whereas mkdssp 4.6.1 overwrites
backbone coordinates as alternate atom rows arrive. Biopython's independently
computed omega uses its selected atoms. Consequently the reference observation
can itself combine measurements from different conformers; this benchmark does
not establish algorithm parity on identical coordinates in those cases.
Use `--jobs 1` for a memory-constrained full PDB run: large structures can exhaust
memory with the default CPU-count concurrency.

## Integrity and reports

Goldens live in `goldens/<dataset>/<feature>.jsonl.gz`, with adjacent
`.jsonl.meta.json` manifests. Manifests bind the compressed bytes, input lock,
contract/schema, reference version and case count. Newly generated files also
record the adapter source digest. Imported older values explicitly mark the
unrecorded original generator fingerprint; it has not been fabricated. The
manifest is published only after complete generation, as the commit marker.
Repository text fingerprints normalize CRLF to LF, so checkout line endings do
not invalidate a run. External input bytes and compressed goldens use exact hashes.

Goldens are streamed with one lookahead record, ordered by fixture, record index
and source ID. Successful records require reference and input identities; selected
cases verify their exact input bytes. The whole file is checked, even with a
case limit, so a small selection can still incur substantial validation I/O.
Input batches flush at 256 cases or 8 MiB of accumulated source text, whichever
comes first; a single large input can exceed that byte threshold. Evaluation
and output memory still depend on structure size. Writer readers must match the stored
reference tool/version.

Reports record their start time in `started_at_unix_ms`, the Git revision/worktree
state, executable hash, current adapter
and contract hashes, and each feature's golden manifest. A summary is written
before evaluation and checkpointed after each feature. Incomplete runs carry
`complete: false` and the error; per-case JSONL preserves completed observations
and emitted writer text. Reference processes have a 300-second deadline; DSSP
invocations have a 120-second deadline. Failed reference batches are isolated
case by case; Kekule panics become errors. These are resource failures, not
agreements. A process abort/OOM can still terminate a run; its incomplete report
must not be interpreted as a complete score.

`cases`, `agrees`, `disagrees` and `errors` cover applicable selected inputs.
`not_applicable` accounts separately for selected IDs with no supplied format.
Input, reference, Kekule, writer-validation and schema errors are reported
separately. Error categories can overlap when both engines fail. An independent
writer-read failure is not automatically attributed to Kekule: per-case error
diagnostics distinguish `kekule_failed` from `writer_validation_failed`, and
validation preserves the original implementation error if no text was emitted.
No absent format is
counted as agreement. A run passes only if it measured at least one applicable
case and all applicable cases agree. Always inspect coverage alongside agreement;
there is no meaningful universal correctness percentage across these features.

`kekule_ms` measures parallel public-API work plus observation construction and
writer contract checks. `reference_ms` measures generation, or writer reading in
comparison mode. Loading, process startup, transport and comparison are excluded.
These are workflow timings, not evidence of a cross-library speed ratio.

## Inputs and maintenance

`corpora/<dataset>/sources.lock.json` defines membership and supplied file hashes.
All matching formats are tested: one ID can contribute both an SDF and a SMILES
case. Source/feature applicability is independent of either engine's success.
All 50,240 Enamine SDF records and all 1,000 PDB entries remain selected. The
historical PubChem sample was preselected for V2000, size and RDKit success; it
cannot establish performance on compounds that were excluded. Do not describe
that sample as unbiased. Current counts are in [GOLDENS.md](GOLDENS.md).

Large inputs belong in `corpora/<dataset>/data/`; full golden archives belong in
`goldens/<dataset>/`. Keep local bundles, alternate generations and reports under
the ignored repository `target/` directory, or outside the checkout. Only smoke
payloads and provenance belong in commits. `data.py verify`, `data.py pack`, and
`data.py unpack` operate on supplied inputs only (see `--help`).
Do not synthesize benchmark molecules. Toy inputs belong in focused regressions.

```text
cargo test -p kekule-bench --locked
python -m unittest discover -s benchmarks/reference -p test_runner.py
python -m unittest discover -s benchmarks/reference/rdkit -p "test_*.py"
python -m unittest discover -s benchmarks/reference/biopython -p "test_*.py"
```

Fix a demonstrated adapter defect or investigate Kekule when values differ.
Never replace a reference with Kekule output or adjust a threshold to conceal a
mismatch. The separate optional trajectory workflow is documented in
[trajectory validation](reference/trajectory/VALIDATION.md); it is not included
in this executable's feature coverage.
