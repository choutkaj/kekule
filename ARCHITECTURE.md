# Architecture

This document defines ownership boundaries and invariants for contributors.
Detailed API contracts, algorithms, numerical policies, and examples belong in
Rustdoc beside their implementations. The links below point to those sources;
`cargo doc --workspace --all-features --no-deps --locked` builds the API reference.

## Ownership map

| Owner | Authoritative responsibility | Must not own |
| --- | --- | --- |
| `Molecule` | One nonempty connected `Graph`, derived `Perception`, definition-scoped `Properties` | Coordinates, hierarchy, system classification |
| `Topology` | Complete molecule instances, reusable definitions, qualified identities, dense layout, one `Hierarchy`, classification, static properties | Geometry or bonds between separate instances |
| `Model` | One shared `Topology`, positions, optional cell, realization properties | A second chemical or hierarchy model |
| `Ensemble` | One shared topology and non-temporal member payloads | An owned `Model` or topology in each member |
| `Trajectory` (`kekule-traj`) | One shared topology and ordered frame payloads, including optional time, step, velocities, and forces | Implicit topology changes between frames |

A salt, solvent box, or protein-ligand complex is a topology containing connected
molecules. Covalent connectedness determines molecule boundaries; hierarchy and
noncovalent interactions do not. Adding a bond between occurrences constructs a
new connected molecule and publishes a new topology.

Dense numerical containers carry values and units, not atom identity or topology
handles. The owning topology or realization translates semantic IDs to dense
indices and validates every associated array and property column.

See the [crate overview](crates/kekule/src/lib.rs),
[molecular owner](crates/kekule/src/core/molecule.rs),
[topology types](crates/kekule/src/topology/mod.rs),
[structure types](crates/kekule/src/structure/mod.rs), and
[trajectory types](crates/kekule-traj/src/trajectory/mod.rs).

## Represented and derived chemistry

`Graph` is authoritative atom, bond, connectivity, and represented stereo state.
Local `AtomId` and `BondId` identify entities within one molecule. Atom chemistry
includes element, isotope, charge, radical, hydrogen declaration, and atom map;
annotations use properties instead of extending every atom or bond with a map.
Bonds carry localized orders. Aromaticity is perceived state, not a bond order.
Source aromatic and stereo syntax is normalized during interpretation.

Represented stereo must refer to valid focuses and adjacent carriers, with the
required focus bond order and consistent assertions and groups. Validate these
conditions before canonicalization. Changing a focus order or deleting a carrier
bond prunes invalid stereo and group membership, even if an alternate graph path
keeps the molecule connected. Unaffected assertions survive.

`Perception` is reconstructible state derived from the exact represented graph
under an explicit model or policy. It contains fundamental chemistry such as
valence, rings, aromaticity, and installed CIP descriptors. Task-specific
descriptors, force-field typing, scoring, and analyses belong in separate result
objects; attaching selected values as properties is deliberate.

Chemical edits invalidate dependent perception. Property changes do not.
Detached perception is checked against graph references and dimensions before
installation; installation must not repair or rewrite represented chemistry.
Default perception does not add atoms, materialize stereo, or assign CIP.
CIP assignment is transactional, respects represented stereo, and reports
unsupported or exhausted ranking rather than treating unfinished work as a tie.

Owning perception runs once per reusable definition and publishes atomically.
A successful model, ensemble, or trajectory perception operation installs a new
shared topology snapshot without copying coordinate payloads. Existing
selections, prepared objects, and streaming bindings remain attached to their
original snapshot. Failure leaves all definitions and realization state intact.

See [graph storage](crates/kekule/src/core/graph.rs),
[represented stereo](crates/kekule/src/core/stereo.rs),
[perception state and validation](crates/kekule/src/core/perception.rs),
[chemical perception](crates/kekule/src/chemistry/perception.rs), and
[stereo algorithms](crates/kekule/src/algorithms/stereo.rs).

## Topology, hierarchy, and classification

A published topology is nonempty and contains only complete connected molecule
instances and used definitions. Definition reuse is explicit; ordinary molecular
construction does not intern chemically equal inputs. Public system traversal is
instance-first. A definition is a storage and reuse mechanism, not another
chemical owner.

System atom and bond identities qualify a molecule-local ID with its instance.
Hierarchy chain, residue, and atom-site IDs are topology-global. Dense indices
are a separate coordinate/property ordering, never interchangeable with semantic
IDs. Publication establishes complete, deterministic mappings between them.

`Hierarchy` is owned exactly once by `Topology`. Atom sites refer to live
`InstanceAtomId` values; they do not copy atoms. One chain or residue may span
molecules, and one molecule may span chains. Molecular and domain-specific
hierarchy views borrow and filter this single hierarchy. Source labels and author
identifiers remain typed metadata, with their namespaces and insertion codes
preserved. Hierarchy membership never fabricates covalent bonds.

Each reusable definition has one `MoleculeClass` shared by all its instances;
each residue has one `ResidueClass`. These are intrinsic classifications, not
contextual roles such as receptor or ligand, and do not affect molecular identity.
Inference occurs at topology publication without running chemical perception.
It combines conservative component recognition, simple graph evidence, and
inter-residue connectivity; conflicting strong evidence yields `Other`.
Explicit assignments retain their intent separately from inferred cached values.

Unchanged complete entities preserve their classifications through subsets and
ordinary appends. Changed chemistry or informative new context triggers
reclassification of affected entities. Editing overrides apply to current
connected components, not historical staging groups. Published hierarchy removal
discards orphaned residue overrides so later ID reuse cannot revive them.
The exact inference and override rules live with
[classification](crates/kekule/src/topology/classification.rs),
[builder publication](crates/kekule/src/topology/builder.rs), and
[component editing](crates/kekule/src/topology/editor.rs).

See also [hierarchy](crates/kekule/src/topology/hierarchy.rs) and
[qualified lookups](crates/kekule/src/topology/lookup.rs).

## Construction, editing, and subsets

Published structural state changes through builders, editors, or explicit
transformations. Drafts may temporarily violate connectedness; publication may
not. Checked boundaries reject invalid references, dimensions, and empty systems.
`validate`, `try_build`, and `try_finish` support recovery without consuming the
draft. Mutable access must preserve safety immediately, without depending on a
mutation guard's destructor being called.

| Surface | Scope of work | Publication |
| --- | --- | --- |
| `MoleculeEditor` | One molecular graph, possibly disconnected while editing | Exactly one nonempty connected molecule |
| `TopologyBuilder` | Complete definitions, instances, hierarchy, classification, properties | Immutable topology |
| `ModelBuilder` | Complete composition with explicit coordinates | Valid model |
| `TopologyEditor` | Occurrence-local atom, bond, component, and hierarchy changes | Repartitioned topology |
| `ModelEditor` | Structural edits coordinated with one realization | Topology and coordinates together |

`MoleculeEditor` is a draft, not a molecule view. Mutable chemistry access
invalidates perception immediately; validated no-op replacements retain state.
System editors clone touched definitions lazily and preserve unaffected
occurrences. Stable draft handles reject deleted or foreign entities. Joining or
splitting components follows asserted connectivity. Surviving model coordinates
are preserved; new atoms require supplied coordinates, never invented geometry.
Geometry-only edits retain the exact shared topology.

Append-oriented construction preserves existing IDs and dense order. General
structural edits publish deterministic new layouts. Hierarchy is filtered and
empty nodes are pruned; its partition is independent of molecule splits/merges.
Per-entity annotations follow the operation's explicit correspondence. Generic
owner annotations are conservatively cleared when their owner changes, including
ambiguous instance splits or merges. No-op edits preserve them.

Complete model append accepts borrowed realizations. It imports definitions,
reuse, perception, classification, hierarchy, coordinates, and entity properties
without interning independent sources or merging equal hierarchy labels.
Coordinates remain in the supplied coordinate system. Incompatible cells,
property types, or units reject the entire append. Owner-property omissions are
reported. Import mappings are scoped to that append and its draft identities.
See the [append contract](crates/kekule/src/structure/model_editor/append.rs) for
cell comparison tolerances, property transfer, and examples.

Selections bind to one exact shared topology, including empty selections and set
operations. They contain unique atoms in authoritative dense order. Whole-residue
expansion retains selected atoms without residue assignments. Structural subsets
may cut molecules and must repartition the induced graph into connected output
definitions. The same operation-specific mapping transfers hierarchy,
classification, properties, and every realization array. Do not add a universal
topology remapping or provenance framework.

See [molecular editing](crates/kekule/src/core/molecule_edit.rs),
[system editing](crates/kekule/src/topology/editor.rs),
[model editing](crates/kekule/src/structure/model_editor.rs),
[selections](crates/kekule/src/topology/selection.rs), and
[subsets](crates/kekule/src/topology/transform.rs).

## Properties and units

Properties are annotations owned at the narrowest scope whose lifetime matches
their validity. One shared `Properties` / `PropertyTable` substrate stores owner
scalars and typed per-entity columns; do not introduce parallel atom/bond data
containers or maps inside every repeated entity.

| Scope | Permitted entity domains |
| --- | --- |
| Molecule definition | Local atoms and bonds |
| Topology | Instances, qualified atoms and bonds, chains, residues, atom sites |
| Model, ensemble member, trajectory frame | Realization atoms and bonds |
| Ensemble or trajectory collection | Owner values |

Keys are validated. A column has one type and, for real values, one physical unit;
missing entries are explicit. Column length matches its owner's domain. Compatible
units convert to the stored unit; incompatible types or dimensions fail. Borrowed
property reads avoid unnecessary string copies. Realization installation rejects
populated properties in unsupported domains even when dimensions happen to match.

Generic properties do not define chemical identity, topology layout, or perception.
Transformations transfer entity values explicitly and do not infer annotation
validity or recompute arbitrary properties. Arbitrary source fields stay in
format sidecars unless canonical semantics, scope, and domain justify promotion.
Occupancy and B factors use reserved realization atom columns with checked
semantics. Positions, cells, weights, time, steps, velocities, and forces retain
their dedicated APIs.

Kekule uses one runtime unit system across all crates. Public boundaries accept
compatible `Quantity` units; internal numerical state uses the library-wide
canonical units. Mass and amount remain distinct dimensions, with molecular mass
and energy conventions chosen coherently. Unit composition and conversion reject
unrepresentable dimensions or scales; checked numerical storage rejects nonfinite
values. Exact units, tolerances, and conversion policies belong with
[units](crates/kekule/src/units.rs) and
[properties](crates/kekule/src/properties.rs), not in a second architecture table.

## Geometry and realization ownership

A model, ensemble, or trajectory owns its shared topology once. Members and frames
store payloads, not nested models. Borrowed `ModelView` access lets coordinate
algorithms and writers share kernels across models, members, frames, and buffers.
Owned projections explicitly materialize a model while sharing the topology.
There is no implicit collection-of-models constructor or special single-member
ownership model.

Dense positions, velocities, and forces validate values and units without carrying
semantic IDs. Consuming vector constructors and projections transfer storage;
borrowed bulk inputs are evaluated once. The owning realization checks lengths,
property domains, and table dimensions on insertion, replacement, or publication.
Empty property tables may acquire their owner's dimensions; populated tables may
not be silently resized or discarded. Stored editing surfaces enforce these
invariants immediately, including when a guard is forgotten.

Measurements and spatial selections operate on borrowed realizations. Cartesian
measurements do not silently apply periodic imaging. A spatial selection is one
realization's result; reevaluating it across frames is explicit. Numerical
policies and supported geometries live with
[positions](crates/kekule/src/structure/positions.rs),
[models](crates/kekule/src/structure/model.rs),
[ensembles](crates/kekule/src/structure/ensemble.rs),
[measurements](crates/kekule/src/structure/measure.rs), and
[alignment](crates/kekule/src/alignment.rs).

## Trajectories and streaming

A trajectory represents one fixed-topology epoch. Changing chemistry or hierarchy
requires a new topology and explicitly constructed geometry for the next epoch.
Frame selection is separate from atom subsetting: requested order, duplicates,
and empty selections are supported without renumbering stored time or step.
Consuming frame projections transfer payloads rather than cloning them.

File readers receive the topology and interpret file atom index as its dense atom
index. They validate counts and available metadata; equal counts alone do not
prove atom identity. Callers may independently validate semantic atom order.
Streaming buffers bind to the exact topology and publish complete frames
transactionally. Eager reads use the same decoder path through clean EOF.
Path writers stage output and publish only a completed nonempty trajectory;
unsupported fields and collection properties are rejected rather than discarded.

Loaded and streaming transformations share kernels. Copy-returning operations
leave sources intact; explicit in-place operations stage all affected state before
publication. Superposition rotates cells, velocities, and forces consistently.
Periodic reconstruction uses asserted bonds and checks ring closure. Imaging
acts on whole molecules. Unwrapping retains temporal state and requires ordered,
sufficiently close samples; it precedes downsampling. Failed streaming operations
change neither the frame buffer nor temporal state. Diagnostics are opt-in.

RMSF and contact-occupancy reductions consume borrowed frames from the exact source
topology with memory bounded by selected atoms or pairs, not frame count. Failed
observations leave accumulators unchanged. Results retain atom/pair associations;
frame indices are diagnostic labels, not statistical weights. Preprocessing,
including alignment or imaging, is explicit. These are separate analysis results,
not new topology state or a generic reduction framework.

For format profiles, limits, periodic conventions, statistical definitions, and
streaming sequence/reset contracts, see
[trajectory storage](crates/kekule-traj/src/trajectory/collection.rs),
[frame buffers](crates/kekule-traj/src/trajectory/buffer.rs),
[file I/O](crates/kekule-traj/src/io/mod.rs),
[periodic operations](crates/kekule-traj/src/periodic.rs),
[streaming periodic operations](crates/kekule-traj/src/periodic/stream.rs), and
[analysis reductions](crates/kekule-traj/src/analysis/reductions.rs).

## Parsing, interpretation, and export

The authoritative input pipeline is:

```text
source -> syntax Document -> independent Record/Block -> richest Interpretation
                                                       -> borrowed/owned projections
```

Parsing preserves syntax and source metadata. Interpretation constructs represented
chemistry, hierarchy, classification, and available geometry once; it does not run
general perception, standardization, tautomerization, or protonation. Geometry may
inform source stereo even when the caller ultimately requests topology only.
Format convenience functions compose this pipeline without reimplementing it.

| Format scope | Interpretation and projection |
| --- | --- |
| SMILES document | All connected components, with topology projection |
| Molfile document | One model with topology and molecule projections |
| SDF record | Model, title, fields, and reports; each record independent |
| mmCIF block | One model or compatible coordinate models as an ensemble |

Sibling SDF records and mmCIF blocks never merge implicitly. Singular projections
check cardinality; plural molecule projections retain every connected component
in source order. MOL/SDF interpretation supplies deterministic synthetic hierarchy
at topology scope. mmCIF preserves label/author hierarchy identities and interprets
all coordinate-model candidates against the same identity, chemistry, and dense
layout. Malformed rows in unselected models still fail validation. Ensemble
assembly releases redundant topologies while transferring realization payloads;
multiple coordinate models do not imply temporal trajectory semantics.

Interpretations own format reports, metadata, and source correspondence. Borrowed
projections reuse canonical owners; consuming projections explicitly discard richer
information. Document-level reports borrow record reports instead of duplicating
them. Format-specific constructors and save methods do not belong on canonical
owners. See the [public format namespaces](crates/kekule/src/lib.rs),
[MOL interpretation](crates/kekule/src/io/structure_documents.rs),
[SDF records](crates/kekule/src/io/sdf_document.rs), and
[mmCIF interpretation](crates/kekule/src/io/mmcif_interpret/mod.rs).

Writers live in format namespaces and accept the richest supported source.
Models share a borrowed realization path; string-returning conveniences wrap sink
writers. Unsupported representational content must fail explicitly. There is no
universal save trait or implicit trajectory-to-structure export.

- SMILES exports represented chemistry and supported stereo. Canonical output
  requires complete labeling under explicit resource bounds; exhaustion is an
  error, never an unproved canonical answer. Hydrogen and isotope projection must
  agree with emitted syntax and preserve stereo.
- Molfile chooses V2000 when it can represent all supported content, otherwise
  promotes to V3000; explicitly requested versions fail on unsupported content.
  Coordinate stereo must agree with the actual rounded emitted geometry. Writers
  must not invent assertions from an unasserted drawing or hide a conflict.
- SDF represents independent records. mmCIF distinguishes independent model blocks
  from one ensemble's multi-model block. Canonical classification informs entity
  kinds; hierarchy alone does not establish polymer status or contextual roles.
  Expert overrides and source provenance handle finer format distinctions.

Exact capability and numerical policies live with
[SMILES writing](crates/kekule/src/io/smiles/write.rs),
[canonical labeling](crates/kekule/src/io/smiles/canonical.rs),
[MOL/SDF writing](crates/kekule/src/io/molfile_write.rs), and
[mmCIF writing](crates/kekule/src/io/mmcif_write.rs).

## Identity and reconstruction

Molecular equality compares authoritative represented chemistry, excluding
perception, annotations, hierarchy, and topology classification. Topology layout
equality is separate: it includes definitions, instances, classifications,
hierarchy, semantic IDs, and dense order, but excludes perception and properties.
Exact shared snapshot identity is stricter than equal layout. Operations requiring
one `Arc<Topology>` must not silently substitute an independently equal topology.

Geometric correspondence may explicitly pair selected atoms from two exact
snapshots, with each side validated and pairs one-to-one. It need not cover equal
system sizes and does not claim chemical equivalence or rebind either owner.
See [alignment correspondence](crates/kekule/src/alignment.rs) and
[trajectory correspondence](crates/kekule-traj/src/analysis/correspondence.rs).

Reconstruction validates represented graph connectedness first, then property
references/dimensions, then installs checked perception. Topology reconstruction
restores definitions, instances, classes, hierarchy, qualified IDs, and dense
layout consistently before validating realization payloads. Disconnected persisted
graphs must be partitioned or rejected. Export to scientific formats is not exact
native persistence and must not weaken these boundaries. See
[reconstruction regressions](crates/kekule/tests/canonical_reconstruction.rs).

## API and maintenance rules

Use `as_*` for cheap borrowed views, `to_*` for owned conversion while retaining
the source, and `into_*` for consuming conversion. Consuming does not promise zero
allocation, but should transfer compatible storage. Keep mutation on appropriate
editors and semantic owner APIs. Use full `properties` names and avoid redundant
owner wrappers, compatibility aliases, and generic target/remapping abstractions.

When adding state, choose its owner from the ownership map. Substantial derived
results deserve their own concrete types. Format metadata stays with the format
unless deliberate promotion defines its canonical meaning and validity scope.
Keep API inventories and algorithm-specific policies beside the implementation;
update this document when ownership or cross-module invariants change.
