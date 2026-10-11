# Architecture

This document defines ownership boundaries and invariants for contributors.
The [current capability index](docs/capabilities.md) links supported features and
their explicit limits; dated audits describe their recorded revisions.
Detailed API contracts, algorithms, numerical policies, and examples belong in
Rustdoc beside their implementations. The links below point to those sources;
`cargo doc --workspace --all-features --no-deps --locked` builds the API reference.

## Ownership map

| Owner | Authoritative responsibility | Must not own |
| --- | --- | --- |
| `Molecule` | One nonempty connected `Graph`, derived `Perception`, `MoleculeProperties` | Coordinates, hierarchy, system classification |
| `Topology` | Reusable definitions and their complete molecule instances, qualified identities, dense atom and bond order, one `Hierarchy`, classification, `TopologyProperties`, layout identity | Geometry or bonds between separate instances |
| `Model` | One shared `Topology` and one structure of it: dense positions, optional cell, occupancies, B factors, and `StructureProperties`, addressed by dense topology index | A second chemical or hierarchy model; temporal state |
| `Ensemble` | One shared topology, collection owner properties, and a weighted, unordered sample of `EnsembleMember` payloads (positions, optional cell, occupancies, B factors, properties, and a finite positive relative weight) | An owned `Model` or topology in each member; temporal meaning |
| `Trajectory` | One shared topology, collection owner properties, and time-ordered `TrajectoryFrame` payloads (positions, optional cell, properties, and optional velocities, forces, time, and step) | Implicit topology changes between frames; statistical weights; crystallographic site parameters |

`Ensemble` and `Trajectory` are distinct scientific concepts and distinct
types. An ensemble is a weighted sample of a distribution: member order carries
no meaning, every member has a weight, and its statistics must be weighted. A
trajectory is a path produced by dynamics: frame order and spacing are
physical. A trajectory samples an ensemble only under ergodicity, so
`Trajectory::into_ensemble` is an explicit, lossy projection with equal
weights, and nothing converts the other way. Each type owns its topology,
collection properties, and items and implements its own operations, including
those whose meaning is currently the same for both, such as item access,
selection, subsetting, perception, superposition, and RMSD; the two types are
expected to diverge. They may share private per-item helpers, such as row
validation and rigid transforms, and the borrowed `ModelView` kernels, but no
collection storage or generic collection type. Weight- or time-dependent
behavior belongs to one type only; do not add a public abstraction that erases
the difference.

Structures are flat. A structure is one set of coordinates and per-atom state
of a topology: a `Model`, an `EnsembleMember`, or a `TrajectoryFrame`. Each
owns exactly the fields that are meaningful for it; there is no shared payload
type. Positions, cell, and structure properties are common to all three.
Occupancies and B factors are crystallographic per-atom fields of a model or
ensemble member and never belong to a trajectory frame. Velocities,
forces, time, and step belong only to frames, and a weight only to an ensemble
member. Coordinate kernels share code through the borrowed `ModelView`, not
through a common owned payload.

Both live in `kekule::structure`. `kekule-traj` owns file codecs, streaming
reader and writer contracts, the reusable `FrameBuffer`, periodic
preprocessing, and trajectory reductions; it owns no second in-memory
trajectory type.

A salt, solvent box, or protein-ligand complex is a topology containing connected
molecules. Covalent connectedness determines molecule boundaries; hierarchy and
noncovalent interactions do not. Adding a bond between occurrences constructs a
new connected molecule and publishes a new topology.

Dense numerical containers carry values and units, not atom identity or topology
handles. The owning topology or structure translates semantic IDs to dense
indices and validates every associated array and property column.

See the [crate overview](crates/kekule/src/lib.rs),
[molecular owner](crates/kekule/src/core/molecule.rs),
[topology types](crates/kekule/src/topology/mod.rs),
[structure types](crates/kekule/src/structure/mod.rs),
[ensembles](crates/kekule/src/structure/ensemble.rs),
[trajectories](crates/kekule/src/structure/trajectory.rs), and
[streaming contracts](crates/kekule-traj/src/trajectory/mod.rs).

## Represented and derived chemistry

`Graph` is authoritative atom, bond, connectivity, and represented stereo state.
Local `AtomId` and `BondId` identify entities within one molecule. Published IDs
are dense in every ID space; drafts may hold deleted slots, and publication
renumbers survivors in draft order and reports the correspondence. Atom chemistry
includes element, isotope, charge, radical, hydrogen declaration, and atom map;
annotations use properties instead of extending every atom or bond with a map.
Bonds carry localized orders. Aromaticity is perceived state, not a bond order.
Source aromatic and stereo syntax is normalized during interpretation.

Represented stereo must refer to valid focuses and adjacent carriers, with the
required focus bond order and consistent assertions and groups. Validate these
conditions before canonicalization. Changing a focus order or deleting a carrier
bond prunes invalid stereo and group membership, even if an alternate graph path
keeps the molecule connected. Unaffected assertions survive.
Enhanced stereo groups belong to a connected molecule. Splitting a correlated
group across components is rejected at publication; absolute memberships may
split independently. Formats that cannot encode a relationship fail explicitly.

`Perception` is reconstructible state derived from the exact represented graph
under an explicit model or policy. It contains fundamental chemistry such as
valence, rings, aromaticity, conjugation, explicitly prepared resonance groups,
and installed CIP descriptors. Resonance contributors are separately requested
results borrowing the exact source molecule; default perception never enumerates
them or rewrites represented bond orders and charges. Task-specific
descriptors, force-field typing, scoring, and analyses belong in separate result
objects; attaching selected values as properties is deliberate.

Chemical edits invalidate dependent perception. Property changes do not.
Detached perception is checked against graph references and dimensions before
installation; installation must not repair or rewrite represented chemistry.
Default perception does not add atoms, materialize stereo, or assign CIP.
Explicit hydrogens are separate graph atoms. Implicit hydrogens are all attached
hydrogens represented without graph atoms; total counts also include explicit
hydrogen neighbors. An atom's implicit count is either fixed or wholly inferred,
never a mix. Fixed counts remain known without perception, and perception never
infers hydrogens on them or overwrites them; inferred counts are unknown after
chemical edits until perception is recomputed. Explicit/implicit conversion
preserves composition and stereo, fixes a collapsed count that inference would
not reproduce, and retains graph hydrogens whose individual information would be
lost.
CIP assignment is transactional, respects represented stereo, and reports
unsupported or exhausted ranking rather than treating unfinished work as a tie.

Owning perception runs once per reusable definition and publishes atomically.
A successful model, ensemble, or trajectory perception operation installs a new
shared topology snapshot without copying coordinate payloads or static topology
state. The snapshot keeps the source's layout identity, so existing selections,
prepared potentials, and streaming bindings stay attached to their original
snapshot and remain usable with the new one. Failure leaves all definitions and
structure state intact.

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

`Topology` is one flat owner. Definitions and instances are peers: a definition
holds reusable chemistry and an instance is one occurrence of it, and neither
nests inside a separate layout object. Dense orders, hierarchy, classification
overrides, and properties sit beside them. Layout identity is an explicit token
minted when a topology is published, not a consequence of how fields are
allocated or shared.

System atom and bond identities qualify a molecule-local ID with its instance.
Hierarchy chain, residue, and atom-site IDs are topology-global. Dense indices
are a separate coordinate/property ordering, never interchangeable with semantic
IDs. Publication establishes complete, deterministic, constant-time mappings
between them. Dense atom order is stored, not derived: it is chosen at
publication and need not keep an instance's atoms contiguous. Format
interpretation keeps source atom-row order, builders default to instance order
or accept an explicit permutation, and editors and subsets keep surviving atoms
in source order. Dense bond order is instance order, then local bond ID.

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
| `ModelEditor` | Structural edits coordinated with one structure | Topology and coordinates together |

`MoleculeEditor` is a draft, not a molecule view. Mutable chemistry access
invalidates perception immediately; validated no-op replacements retain state.
System editors clone touched definitions lazily and preserve unaffected
occurrences. Stable draft handles reject deleted or foreign entities. Joining or
splitting components follows asserted connectivity. Surviving model coordinates
are preserved; new atoms require supplied coordinates, never invented geometry.
Geometry-only edits retain the exact shared topology.

Append-oriented construction preserves existing IDs and dense order. General
structural edits keep surviving atoms in source dense order and append new
atoms in creation order, while repartitioned molecules receive deterministic new
definitions and instances. Hierarchy is filtered and
empty nodes are pruned; its partition is independent of molecule splits/merges.
Per-entity annotations follow the operation's explicit correspondence. Generic
owner annotations are conservatively cleared when their owner changes, including
ambiguous instance splits or merges. No-op edits preserve them.

Complete model append accepts borrowed structures. It imports definitions,
reuse, perception, classification, hierarchy, coordinates, and entity properties
without interning independent sources or merging equal hierarchy labels.
Coordinates remain in the supplied coordinate system. Incompatible cells,
property types, or units reject the entire append. Owner-property omissions are
reported. Import mappings are scoped to that append and its draft identities.
See the [append contract](crates/kekule/src/structure/model_editor/append.rs) for
cell comparison tolerances, property transfer, and examples.

Selections bind to one topology layout, including empty selections and set
operations. Atom and bond selections independently contain unique entities in
authoritative dense order. Atom-to-bond conversion explicitly chooses both,
either, or exactly one selected endpoint; bond-to-atom conversion selects both
endpoints. Whole-residue and whole-chain expansion retain selected atoms without
hierarchy assignments. Structural subsets
may cut molecules and must repartition the induced graph into connected output
definitions. The same operation-specific mapping transfers hierarchy,
classification, properties, and every structure array. Do not add a universal
topology remapping or provenance framework.

See [molecular editing](crates/kekule/src/core/molecule_edit.rs),
[system editing](crates/kekule/src/topology/editor.rs),
[model editing](crates/kekule/src/structure/model_editor.rs),
[selections](crates/kekule/src/topology/selection.rs), and
[subsets](crates/kekule/src/topology/transform.rs).

## Properties and units

Properties are annotations owned at the narrowest scope whose lifetime matches
their validity. Each scope has one typed owner built from `OwnerProperties`
scalars and `PropertyTable<R>` columns, where the row type `R` names the entity
domain; do not introduce parallel atom/bond data containers or maps inside every
repeated entity.

| Scope | Owner | Row domains |
| --- | --- | --- |
| Molecule definition | `MoleculeProperties` | `AtomId`, `BondId` |
| Topology | `TopologyProperties` | Instances, dense atoms and bonds, chains, residues, atom sites |
| Structure (model, ensemble member, trajectory frame) | `StructureProperties` | Dense atoms and bonds |
| Ensemble or trajectory collection | `OwnerProperties` | Owner values only |

Unsupported domains are unrepresentable rather than rejected at runtime. Reads
take the table's row type. Writes go through length-preserving
`*Mut` guards: whole-column replacement or removal, single-cell writes, and
transactional batches with ordered update semantics. A draft editor addresses
allocated slots, and writing a value to a deleted slot fails with `RemovedRow`.

Keys are validated and never reserved. A column has one type and, for real
values, one physical unit; missing entries are explicit. Column length matches
its owner's domain. Compatible units convert to the stored unit; incompatible
types or dimensions fail. Borrowed property reads avoid unnecessary string copies.

Generic properties do not define chemical identity, topology layout, or perception.
Editing a published topology's properties keeps its layout identity, as
perception does; only builders, editors, and transformations publish a new
layout. Transformations transfer entity values explicitly and do not infer annotation
validity or recompute arbitrary properties. Arbitrary source fields stay in
format sidecars unless canonical semantics, scope, and domain justify promotion.
Occupancy and B factors are typed site fields of models and ensemble members,
with checked semantics independent of any generic property with the same name. Positions, cells,
weights, time, steps, velocities, and forces retain their dedicated APIs.

Kekule uses one runtime unit system across all crates. Public boundaries accept
compatible `Quantity` units; internal numerical state uses the library-wide
canonical units. Mass and amount remain distinct dimensions, with molecular mass
and energy conventions chosen coherently. Unit composition and conversion reject
unrepresentable dimensions or scales; checked numerical storage rejects nonfinite
values. Exact units, tolerances, and conversion policies belong with
[units](crates/kekule/src/units.rs) and
[properties](crates/kekule/src/properties.rs), not in a second architecture table.

## Geometry and structure ownership

A model, ensemble, or trajectory owns its shared topology once. Members and frames
store flat payloads, not nested models; a detached payload has
atom rows only and acquires bond rows when its collection binds it. Borrowed
`ModelView` access, obtained through `AsModelView`, lets coordinate algorithms
and writers share kernels across models, members, frames, and buffers. Owned
projections (`to_model`) explicitly materialize a model while sharing the
topology. Collections built from models (`from_models`) consume them and require
one shared layout; there is no special single-member ownership model.

Read access is through views: `Topology` lookups return `Option` views
(`AtomView`, `BondView`, molecule, chain, residue, and atom-site views), and
`Model::atom` joins an `AtomView` with that atom's structure state. Mutation
goes through explicit, dimension-preserving setters on models, detached
payloads, and the collection item guards returned by `get_mut()`, and through
`properties_mut()` on topologies (`topology_properties_mut()` on models and
collections, which replaces a shared topology with a snapshot of it). Guards
have no `DerefMut`; whole payloads are replaced through the owner.

Dense positions, velocities, and forces validate values and units without carrying
semantic IDs. Consuming vector constructors and projections transfer storage;
borrowed bulk inputs are evaluated once. The owning structure checks lengths,
property domains, and table dimensions on insertion, replacement, or publication.
Empty property tables may acquire their owner's dimensions; populated tables may
not be silently resized or discarded. Stored editing surfaces enforce these
invariants immediately, including when a guard is forgotten.

Measurements and spatial selections operate on borrowed structures. Cartesian
measurements do not silently apply periodic imaging. A spatial selection is one
structure's result; reevaluating it across frames is explicit. Prepared bond
dihedrals bind topology-selected references to one topology layout. CIP
priority and atom-ID ties choose the references without coordinates; undefined
geometry yields an absent value rather than substituting references. Numerical
policies and supported geometries live with
[positions](crates/kekule/src/structure/positions.rs),
[models](crates/kekule/src/structure/model.rs),
[ensembles](crates/kekule/src/structure/ensemble.rs),
[measurements](crates/kekule/src/structure/measure.rs), and
[alignment](crates/kekule/src/alignment.rs).

Cartesian geometry edits mutate `Model` coordinates atomically while retaining
the exact topology, cell, and properties. Prepared distance, angle, and dihedral
edits bind references and moving selections to that layout, never to cached
coordinates. Automatic fragments must separate from fixed references; explicit
selections may deliberately deform boundary bonds. Neither path performs ring
closure, relaxation, or periodic imaging. See the
[geometry editing contract](crates/kekule/src/structure/geometry_edit.rs).

## Trajectories and streaming

A trajectory represents one fixed-topology epoch. Changing chemistry or hierarchy
requires a new topology and explicitly constructed geometry for the next epoch.
Frame selection is separate from atom subsetting: requested order, duplicates,
and empty selections are supported without renumbering stored time or step.
Consuming frame projections transfer payloads rather than cloning them.

File readers receive the topology and interpret file atom index as its dense atom
index. They validate counts and available metadata; equal counts alone do not
prove atom identity. Callers may independently validate semantic atom order.
Because interpretation keeps source atom-row order and structure writers emit
dense order, a coordinate file written in a structure file's row order addresses
the same atoms through that file's topology.
Streaming buffers bind to one topology layout and publish complete frames
transactionally. Eager reads use the same decoder path through clean EOF.
Path writers stage output and publish only a completed nonempty trajectory;
unsupported fields and collection properties are rejected rather than discarded.

Loaded and streaming transformations share kernels. Collection operations,
including superposition and periodic preprocessing, change the collection in
place and stage all affected state before publication; clone first to keep the
source. Superposition and RMSD are methods of both trajectories and ensembles, against
one of its items or an independent view; streaming superposition applies the
same kernel to a buffer. Streaming tools take the caller's frame index for
diagnostics. Superposition rotates cells, velocities, and forces consistently.
Periodic reconstruction uses asserted bonds and checks ring closure. Imaging
acts on whole molecules. Unwrapping retains temporal state and requires ordered,
sufficiently close samples; it precedes downsampling. Failed streaming operations
change neither the frame buffer nor temporal state. Diagnostics are opt-in.

RMSF and contact-occupancy reductions are trajectory statistics in which every
frame counts once; they do not accept ensembles, whose statistics must be
weighted. They consume borrowed frames sharing the source
topology's layout with memory bounded by selected atoms or pairs, not frame count. Failed
observations leave accumulators unchanged. Results retain atom/pair associations;
frame indices are diagnostic labels, not statistical weights. Preprocessing,
including alignment or imaging, is explicit. These are separate analysis results,
not new topology state or a generic reduction framework.

For format profiles, limits, periodic conventions, statistical definitions, and
streaming sequence/reset contracts, see
[trajectory payloads](crates/kekule/src/structure/trajectory.rs),
[collection alignment](crates/kekule/src/alignment/collection.rs),
[frame buffers](crates/kekule-traj/src/trajectory/buffer.rs),
[file I/O](crates/kekule-traj/src/io/mod.rs),
[periodic operations](crates/kekule-traj/src/periodic.rs),
[streaming periodic operations](crates/kekule-traj/src/periodic/stream.rs), and
[analysis reductions](crates/kekule-traj/src/analysis/reductions.rs).

## Parsing, interpretation, and export

The authoritative input pipeline is:

```text
reader -> record/block splitter -> format document (typed syntax + locations)
                                            |
                         options + context -+
                                            v
                         interpretation { canonical owner, report } -> projections
```

The unit of parsing and interpretation is the independent record or block: one
Molfile, one SDF record, one SMILES line, one CIF data block, or one PDB entry.
SDF and SMILES readers split input into records by line framing alone (`$$$$`
lines, line ends), without reading record content, so a malformed record fails
alone and these formats stream with memory bounded by one record. CIF block
boundaries are lexical instead: a semicolon-delimited text value may contain a
line that looks like `data_`, so finding blocks requires tokenizing, and
malformed framing such as an unterminated text value fails the rest of the input
rather than one block. A file-level document is a convenience collection of
parsed records; interpretation is always per record and never depends on a
sibling record.

A format document is the typed syntax of one record or block, with source
locations: what the source states, before chemistry. Anything that requires
reading source text, including CXSMILES extension fields and Molfile property
records, belongs to parsing; interpretation consumes typed syntax and never
re-reads text. Documents keep source locations rather than byte-exact formatting.
CIF documents stay generic data trees because most categories lie outside
Kekule's interpretation.

Interpretation takes a document, explicit options, and optional context: external
resources a format needs but does not carry. It constructs represented chemistry,
hierarchy, classification, and available geometry once into the richest canonical
owner plus a report; it does not run general perception, standardization,
tautomerization, or protonation. Geometry may inform source stereo even when the
caller ultimately requests topology only. Format convenience functions compose
this pipeline without reimplementing it.

Interpretation never silently discards source content. Content the interpreter
does not support is rejected by default; an explicit option downgrades it to typed
report entries that name the record and the kind of content omitted. A deliberately
lossy projection, such as base-only SMILES, says so in its name.

Every stage reports through one diagnostics model. Parse and interpretation
errors carry a structured kind and a source location relative to the original
input, so callers can match on failures and map them to the source without
format-specific arithmetic. Non-fatal report issues carry a structured kind and
identify the source records they concern.

Coordinate-only trajectory formats are the deliberate exception: they carry no
chemistry to interpret, so streaming codecs decode frames directly against a
caller-supplied topology (see Trajectories and streaming).

| Format scope | Interpretation and projection |
| --- | --- |
| SMILES record | All connected components, with topology projection |
| Molfile | One model with topology and molecule projections |
| SDF record | Model, title, fields, and reports; each record independent |
| CIF data block | One model or compatible coordinate models as an ensemble |
| PDB entry | One model or compatible coordinate models as an ensemble |

Sibling SDF records and mmCIF blocks never merge implicitly. Singular projections
check cardinality; plural molecule projections retain every connected component
in source order. MOL/SDF interpretation supplies deterministic synthetic hierarchy
at topology scope. mmCIF preserves label/author hierarchy identities and interprets
all coordinate-model candidates against the same identity, chemistry, and dense
layout. PDB carries author identities only: a residue is a run of consecutive
records, and a number that recurs after other residues names a new residue.
Malformed rows in unselected models still fail validation. Ensemble assembly
releases redundant topologies while transferring structure payloads; multiple
coordinate models do not imply temporal trajectory semantics.

mmCIF alternate rows remain in the source document. Interpretation selects whole
residue alternatives or explicitly constrained groups and records its choices in
the format report. Occupancy ranks representatives, not whole-structure
probabilities. Coordinate-model ensemble conversion never expands alternate
labels; caller-supplied conformation selections must pass the same identity and
topology checks. mmCIF carries no model weights, so every coordinate model
becomes an equally weighted member; occupancy is never converted into a weight.

Macromolecular formats (mmCIF and PDB) do not carry complete bonds, and molecule
boundaries follow bonds, so one resolver shared by these formats derives connectivity
from ranked evidence:

1. Explicit source bonds: component bond tables in the source, and covalent,
   disulfide, and metal-coordination links between residues. Records that list
   connectivity without bond types, such as PDB `CONECT`, supplement a
   residue's template and replace distance inference for residues without one.
2. Component templates supplied as interpretation context, keyed by component
   identifier and matched to source atoms by atom name and element.
3. Polymer linkage between residues adjacent in sequence, accepted only within
   bonding distance; a longer gap is a reported chain break.
4. Distance inference from covalent radii for residues that no template or
   explicit bond table covers. It is on by default, can be disabled, and yields
   single bonds that the report identifies as inferred.

Without entity records, as in PDB, chain termination records delimit the
polymer, and a residue without a polymer template takes its linkage from its
neighbour or its atom names.

Explicit bonds outrank the levels below them. A residue whose internal bonds the
source enumerates, as files that write complete connectivity do, takes exactly
those bonds; templates and inference only fill residues the source leaves without
them. Close contacts that no evidence bonds are counted, never bonded.

Coordination is covalent: metal-ligand links become `Dative` bonds, so a metal
complex or a metalloprotein site is one connected molecule. Links whose records
state no type, such as PDB `LINK`, are coordination between a metal and a
non-metal. Hydrogen bonds, salt bridges, and links to symmetry mates are never
bonds; unresolvable links are reported. A heavy atom whose name does not match
its residue template is an error, not a silently separate molecule. The report
records the evidence for every bond and names residues whose resolved atoms
split into several molecules.

Formal charges come from the source when stated, otherwise from the template,
including when the source's own bond table supplies the bonds. Hydrogens present
in the model are atoms; remaining hydrogens follow from valence, so a template's
protonation is never asserted as a fixed hydrogen count. When a residue's model
hydrogens differ from its template on a heteroatom, a variant template that
matches them exactly, such as a histidine tautomer, is selected; otherwise each
added or removed hydrogen on a heteroatom whose template bonds are all single
changes its charge by one, as one added or removed hydrogen does on a ring
nitrogen whatever its bonds (a protonated cytosine, an imidazolate), and any
other difference fails. Hydrogens absent from carbon are left to valence.

Kekule bundles templates for standard polymer residues, their common protonation
variants, water, common ions, and common cofactors, taken from the Chemical
Component Dictionary. Callers supply further components, up to the complete
dictionary, as interpretation context. Force fields name protonation states,
termini, water, and ions with residue and atom names of their own; libraries
for the AMBER and CHARMM names derive their templates from the bundled entries
and are opt-in context, since force-field names collide with dictionary
identifiers.

Interpretations own format reports, metadata, and source correspondence. Borrowed
projections reuse canonical owners; consuming projections explicitly discard richer
information. Document-level reports borrow record reports instead of duplicating
them. Format-specific constructors and save methods do not belong on canonical
owners. See the [public format namespaces](crates/kekule/src/lib.rs),
[MOL interpretation](crates/kekule/src/io/structure_documents.rs),
[SDF records](crates/kekule/src/io/sdf_document.rs),
[PDB records](crates/kekule/src/io/pdb_document.rs), and
[structure interpretation](crates/kekule/src/io/structure_interpret/mod.rs).

Writers live in format namespaces as one `write(source, options)` returning a
string and one `write_to(writer, source, options)` sink writer per format. Each
source type converts from the inputs that format can represent: molecules,
models, collection items, interpreted records, or mmCIF block sources carrying
classifications or reports. Models and collection items share a borrowed
structure path. Coordinate-bearing writers emit atoms in dense order.
Unsupported representational content must fail explicitly. There is no
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

Molecular queries own predicate graphs, recursive subqueries, output tags, and
local stereo carrier frames; they do not reuse represented molecular graphs.
Matching reads explicit target perception without installing or changing it.
Prepared molecular targets borrow their source; topology matches retain the
exact shared snapshot. Connected component-local searches may be reused per
definition and expanded to instances. Disconnected queries use qualified
occurrence identities without creating bonds between molecules. See the
[SMARTS contract](crates/kekule/src/query/dialect.md).

Molecular equality compares authoritative represented chemistry, excluding
perception, annotations, hierarchy, and topology classification. Topology layout
equality is separate: it includes definitions, instances, classifications,
hierarchy, semantic IDs, and dense order, but excludes perception and properties.
Layout identity is stricter than equal layout. Every publication mints a new
layout identity; perception and topology property edits create snapshots that
keep it. Operations that
address atoms and bonds by index accept any snapshot sharing their layout, and
must not silently substitute an independently equal topology. Operations that
read perception, such as prepared substructure targets, bind the exact snapshot.

Geometric correspondence may explicitly pair selected atoms from two exact
snapshots, with each side validated and pairs one-to-one. It need not cover equal
system sizes and does not claim chemical equivalence or rebind either owner.
See [alignment correspondence](crates/kekule/src/alignment/correspondence.rs) and
[collection alignment](crates/kekule/src/alignment/collection.rs).

Reconstruction validates represented graph connectedness first, then property
references/dimensions, then installs checked perception. Topology reconstruction
restores definitions, instances, classes, hierarchy, qualified IDs, and dense
layout consistently before validating structure payloads. Disconnected persisted
graphs must be partitioned or rejected. Export to scientific formats is not exact
native persistence and must not weaken these boundaries. See
[reconstruction regressions](crates/kekule/tests/molecule/identity.rs).

## API and maintenance rules

Use `as_*` for cheap borrowed views, `to_*` for owned conversion while retaining
the source, and `into_*` for consuming conversion. Consuming does not promise zero
allocation, but should transfer compatible storage. Constructors that take
ownership of chemistry accept values (`Molecule`, `Model`), and callers clone
explicitly. Lookups by identity return `Option`; fallible operations return a
typed error. Keep mutation on appropriate editors, `*_mut` guards, and semantic
owner APIs. Each read lives once on its owner and is reached through `Deref`
or a view rather than forwarded. Use full `properties` names and avoid redundant
owner wrappers, compatibility aliases, and generic target/remapping abstractions.

When adding state, choose its owner from the ownership map. Substantial derived
results deserve their own concrete types. Format metadata stays with the format
unless deliberate promotion defines its canonical meaning and validity scope.
Keep API inventories and algorithm-specific policies beside the implementation;
update this document when ownership or cross-module invariants change.

The unpublished `kekule-bench` workspace package (`benchmarks/observer`) reports
what Kekule's public APIs compute, keyed by source positions, for comparison
with pinned reference toolkits. Datasets, reference adapters, comparison and the
list of known differences stay outside the runtime crates; reference toolkits are
never runtime dependencies. See the [benchmark README](benchmarks/README.md).

`kekule-openff` owns compiled SMIRNOFF rules, configured NAGL inference, and typed
`ParameterizedTopology` results retaining the caller's `Arc<Topology>`.
It prepares temporary per-definition chemistry and expands assignments to
instance-qualified atoms without mutating the topology or adding hydrogens.
All numeric parameters carry canonical units. Complete force-field state stays
in this result rather than unstructured molecule properties. The crate is pure
Rust: charge lookup selects an entry only for that entry's exact molecule, so no
InChI implementation is linked. Checkpoint conversion and the provenance of
shipped data live in `tools/openff`.
Model bundles own feature ordering, supported network configuration, domain and
lookup data. OFFXML retains the required checkpoint identity; parameterization
checks it against the supplied model before assigning any molecule. Chemistry
preparation is a versioned model profile, separate from SMIRNOFF MDL perception.
The converted Ash model that the bundled Rosemary preset requires ships as data
in the separate `kekule-openff-ash` crate, a default feature, so the preset works
without external files while the main crate stays within the crates.io package
size limit. Parameterization has no molecule size limit: its work grows linearly
with the molecule, and bounded searches only stop pathological patterns. Every
failure is one typed error with a stable kind and the failing molecule definition.

`kekule-potentials` owns energy evaluation and geometry optimization; the
foundational crate owns no potential. A prepared `Potential` binds one topology
layout, does not change during evaluation, and accepts any `ModelView` sharing
that layout, including perceived snapshots; an independently equal topology is
incompatible. Energies and
gradients are separate result objects in canonical units, not model properties.
Backends lower explicit parameters, currently an OpenFF `ParameterizedTopology`,
into private functional-form kernels. Parameter assignment stays in its own
companion crate. Unsupported structure state, such as a periodic cell, is
rejected rather than ignored. Gradients are exact derivatives and never capped;
coordinates fail only where a requested energy or gradient is undefined.
Minimization reads a borrowed view and returns new positions without changing
chemistry, topology, cell, or properties. See the
[potential contract](crates/kekule-potentials/src/lib.rs) and
[OpenFF backend](crates/kekule-potentials/src/openff.rs).
