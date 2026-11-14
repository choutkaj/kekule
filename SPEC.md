# Specification

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
| `Geometry` | One arrangement of a topology's atoms, in dense atom order: positions, optional cell, occupancies, B factors, velocities, forces, and `GeometryProperties` | A topology; a place in an order; a membership in a distribution |
| `Structure` | One shared `Topology` and one `Geometry` of it | A second chemical or hierarchy model; temporal state |
| `StructureSet` | One shared topology, set owner properties, `Geometry` values in a stable order (each shared and copied on write), and an optional `Order` with one value per geometry (a `TimeAxis` of times and steps, steps, or one coordinate in a unit) | An order value on a geometry; a weight, a state, or a distribution |
| `Ensemble` | One `StructureSet` and one `Distribution` reading its geometries as draws: samples (all of the set or a sorted subset), log weights, one `Partition`, one `State` per state, and a temperature | Geometries of its own; an energy on a geometry; a stored population |

### Vocabulary

Kekule names what a system is made of with these words and uses no others for
these concepts.

A **geometry** is where the atoms of one topology are, plus the per-atom state
recorded at that arrangement: positions in dense atom order, an optional
periodic cell, occupancies, B factors, velocities, forces, and geometry
properties. A geometry carries no topology; it is *bound* to one when it has a
property row per atom and bond of that topology.

A **structure** is a topology and one geometry of it, so `Structure` is
`Topology` plus `Geometry`. The topology stays coordinate-free and may be
shared with other structures, sets, frames, selections, and prepared
potentials.

A **structure set** is one shared topology and many geometries of it, in a
stable order, with an optional `Order` saying what that order means.

An **ensemble** is a structure set read as draws from a distribution: which of
its geometries are samples, what each weighs, and how they group into states.

A **frame** is a geometry as a trajectory file stores it, with its time and
step. Frames belong to `kekule::traj` alone.

**Sample** and **state** belong to the statistical reading alone: a sample is
a geometry an ensemble's distribution covers, and a state is a group of
samples with one representative.

New types take the prefix of what they belong to: `Geometry*`, `Structure*`,
`StructureSet*`, `Ensemble*`, `Frame*`. `Model` never names a structure. It
stays for perception models (`AromaticityModel`, `ValenceModel`), physical and
machine-learning models (`WaterModel`, NAGL bundles), and the `MODEL` records
of PDB and mmCIF, which this document calls coordinate models.

A salt, solvent box, or protein-ligand complex is a topology containing connected
molecules. Covalent connectedness determines molecule boundaries; hierarchy and
noncovalent interactions do not. Adding a bond between occurrences constructs a
new connected molecule and publishes a new topology.

Dense numerical containers carry values and units, not atom identity or topology
handles. The owning topology or structure translates semantic IDs to dense
indices and validates every associated array and property column.

See the [crate overview](crates/kekule/src/lib.rs),
[molecular owner](crates/kekule/src/molecule/owner.rs),
[topology types](crates/kekule/src/topology/mod.rs),
[structure types](crates/kekule/src/structure/mod.rs),
[structure sets](crates/kekule/src/structure/structure_set.rs), and
[streaming contracts](crates/kekule/src/traj/trajectory/mod.rs).

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
required focus bond order and consistent assertions and groups, and an
implicit-hydrogen carrier must agree with a fixed hydrogen declaration. One
validator checks these conditions for explicit stereo edits, publication, and
`MoleculeEditor::validate`, before canonicalization; graph-edit pruning checks
references alone, since hydrogen transforms change a carrier and its
declaration in separate steps. Changing a focus order, adding a bond to a
focus, or deleting a carrier bond prunes the stereo and group membership it
invalidates, even if an alternate graph path keeps the molecule connected.
Pruning is local to the edited atoms and bonds, so unaffected assertions
survive, and the editor reports every pruned element.
Enhanced stereo groups belong to a connected molecule. Splitting a correlated
group across components is rejected at publication; absolute memberships may
split independently. Formats that cannot encode a relationship fail explicitly.

`Perception` is reconstructible state derived from the exact represented graph
under an explicit model or policy. It contains fundamental chemistry such as
valence, rings, aromaticity, conjugation, resonance delocalization, and
installed CIP descriptors. Resonance delocalization marks, per bond, whether
moving a formal charge changes its order, with the bond's connected group; it
reads the graph locally and enumerates no contributor. Resonance contributors
are separately requested results borrowing the exact source molecule; default
perception never enumerates them or rewrites represented bond orders and
charges. Task-specific
descriptors, force-field typing, scoring, and analyses belong in separate result
objects; attaching selected values as properties is deliberate.

`Perception` exposes only its sections (valence, rings, aromaticity,
conjugation, resonance, and CIP state); each section owns its queries and
records mandatory provenance: the model, plus the policy where results depend
on it. Valence records strict or permissive validation; a ring basis records
the requested model separately from whether the depth-first fallback chose
its rings. Consumers reuse installed state only when its provenance matches
what they would perceive themselves, rather than revalidating it.

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
A definition whose installed perception is the default profile's, installed by
perception and unchanged since, is shared rather than perceived again, so
perceiving after an edit costs the definitions the edit changed.
A successful structure or structure set perception operation installs a new
shared topology snapshot without copying coordinate payloads or static topology
state. The snapshot keeps the source's layout identity, so existing selections,
prepared potentials, and streaming bindings stay attached to their original
snapshot and remain usable with the new one. Failure leaves all definitions and
geometry state intact.

See [graph storage](crates/kekule/src/molecule/graph.rs),
[represented stereo](crates/kekule/src/molecule/stereo.rs),
[perception state and validation](crates/kekule/src/molecule/perception.rs),
[chemical perception](crates/kekule/src/chemistry/perception.rs), and
[stereo algorithms](crates/kekule/src/algorithms/stereo.rs).

## Topology, hierarchy, and classification

A published topology is nonempty and contains only complete connected molecule
instances and used definitions. Definition reuse is explicit; ordinary molecular
construction does not intern chemically equal inputs. Structure-file
interpretation reuses definitions by an explicit rule (see
[parsing and interpretation](#parsing-interpretation-and-export)). Public
system traversal is instance-first. A definition is a storage and reuse
mechanism, not another chemical owner.

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

`Hierarchy` is owned exactly once by `Topology`. Its storage is private: the
chain, residue, and atom-site views are the only read path, and builders stage
it through one `HierarchyBuilder`. Atom sites refer to live
`InstanceAtomId` values; they do not copy atoms. One chain or residue may span
molecules, and one molecule may span chains. Molecular and domain-specific
hierarchy views borrow and filter this single hierarchy. Source labels and author
identifiers remain typed metadata, with their namespaces and insertion codes
preserved. Hierarchy membership never fabricates covalent bonds. Publication
derives one private hierarchy index beside it, shared by every snapshot of the
layout: each site's dense atom and each dense atom's site, and, built on first
use, each instance's sites, residues, and chains, so hierarchy walks cost what
they return. The index numbers the sites' atom names against a table of the
distinct names on first use, which name terms of the selection language read.
The hierarchy stores residues and sites by column; residue labels and site
metadata are keys into immutable row tables shared by every snapshot and draft
that keeps them, so republishing an edited hierarchy copies integers.

Each reusable definition has one `MoleculeClass` shared by all its instances;
each residue has one `ResidueClass`. These are intrinsic classifications, not
contextual roles such as receptor or ligand, and do not affect molecular identity.
Inference occurs at topology publication without running chemical perception.
It combines conservative component recognition, simple graph evidence, and
inter-residue connectivity; conflicting strong evidence yields `Other`.
Explicit assignments retain their intent separately from inferred cached values,
and `ClassAssignment` reports which one a published class is. A residue's
explicit class lives on its staged node, so replacing a staged hierarchy
discards it; clearing an assignment returns the entity to inference.

Unchanged complete entities preserve their classifications through subsets and
ordinary appends. Changed chemistry or informative new context triggers
reclassification of affected entities. Editing overrides apply to current
connected components, not historical staging groups. Published hierarchy removal
discards orphaned residue overrides so later ID reuse cannot revive them.
The exact inference and override rules live with
[classification](crates/kekule/src/topology/classification.rs),
[builder publication](crates/kekule/src/topology/builder.rs), and
[component editing](crates/kekule/src/topology/editor.rs).

One-letter codes and sequences are derived, not hierarchy state. A residue's
code follows from its component identity: the code a component library or the
bundled dictionary codes state, else its parents' codes, else the unknown (`X`
or `N`). Only a monomer of a polymer has a code: its molecule holds another
monomer, the file numbers it in a polymer sequence, or it has no bonds, as in
a Cα trace. A free nucleotide or amino acid bound as a ligand has none. A code
names its component's parents, so a modified residue is told from a standard
one. A chain's sequence joins the codes of the monomers it holds. The
sequence a file states for a polymer entity, with monomers the structure does not
observe, is format metadata: the mmCIF interpretation report owns it and
places it on the interpreted structure's chains and residues. See
[sequences](crates/kekule/src/sequence.rs) and
[entity sequences](crates/kekule/src/io/structure_interpret/sequences.rs).

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
| `StructureBuilder` | Complete composition with explicit coordinates | Valid structure |
| `TopologyEditor` | Occurrence-local atom, bond, component, and hierarchy changes | Repartitioned topology |
| `StructureEditor` | Structural edits coordinated with one structure | Topology and coordinates together |

`MoleculeEditor` is a draft, not a molecule view. Mutable chemistry access
invalidates perception immediately; validated no-op replacements retain state.
System editors clone touched definitions lazily and preserve unaffected
occurrences; publication shares every untouched definition, with its
perception, with the source snapshot, whether or not that snapshot is still
shared, so an edit copies what it changes. A draft started from a snapshot
stores nothing per source atom, bond, or hierarchy node: it reads the source
until an edit changes, moves, or deletes one, and a structure draft shares the
source's positions, occupancies, and B factors until one is written. The
source snapshot is never written, so undo histories keep every snapshot
intact. Publication carries untouched occurrences, residues, and sites over
from the source by index, classifies only what the edit cleared, and costs
the edit plus a few passes over plain per-atom arrays. Stable draft handles reject deleted or foreign entities. Joining or
splitting components follows asserted connectivity. Surviving structure coordinates
are preserved; new atoms require supplied coordinates, never invented geometry,
unless the host asks a placement operation on a `StructureEditor` draft
(`structure::placement`) to derive them from ideal local geometry. Placement
is part of the host's edit and its one publication, never moves an existing
atom, works inside the edited occurrences so its cost follows the edited
atoms rather than the system, and relaxes nothing. Geometry-only edits retain
the exact shared topology.

Append-oriented construction preserves existing IDs and dense order. General
structural edits keep surviving atoms in source dense order and append new
atoms in creation order, while repartitioned molecules receive deterministic new
definitions and instances. Hierarchy is filtered and
empty nodes are pruned; its partition is independent of molecule splits/merges.
Per-entity annotations follow the operation's explicit correspondence. An
editor's correspondence maps surviving handles to published identities and,
densely, the source topology's atom and bond indices, residues, and chains to
the published ones and back, so hosts remap their columns in one pass. Generic
owner annotations are conservatively cleared when their owner changes, including
ambiguous instance splits or merges. No-op edits preserve them.

Complete structure append accepts borrowed structures. It imports definitions,
reuse, perception, classification, hierarchy, coordinates, and entity properties
without interning independent sources or merging equal hierarchy labels.
Coordinates remain in the supplied coordinate system. Incompatible cells,
property types, or units reject the entire append. Owner-property omissions are
reported. Import mappings are scoped to that append and its draft identities.
See the [append contract](crates/kekule/src/structure/structure_editor/append.rs) for
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

See [molecular editing](crates/kekule/src/molecule/molecule_edit.rs),
[system editing](crates/kekule/src/topology/editor.rs),
[structure editing](crates/kekule/src/structure/structure_editor.rs),
[selections](crates/kekule/src/topology/selection.rs), and
[subsets](crates/kekule/src/topology/transform.rs).

The selection language (`select`) is text that evaluates to an
`AtomSelection`; it is a query surface over these owners, not a second
selection type. A parsed `SelectionExpr` is independent of any topology and
prints canonically. Evaluation binds the result to the evaluated topology,
reads hierarchy, classification, properties, and installed perception without
changing them or perceiving, and reads coordinates only when evaluated on a
structure, so `within` is one structure's result like other spatial
selections. Identifier terms follow author identifiers, falling back to label
identifiers, and `label_*` terms read label identifiers. Host applications
extend the language through the `SelectionResolver` extension trait, whose
answers must bind to the evaluated layout or be points for a spatial
reference; other structures enter an evaluation only as such points.
Evaluation may run `within` on a host's `parallel::Executor`, which never
changes the selected atoms. See the
[grammar](docs/selection-language.md) and the
[evaluator](crates/kekule/src/selection_language/eval.rs).

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
| Geometry | `GeometryProperties` | Dense atoms and bonds |
| Structure set | `OwnerProperties` | Owner values only |
| State of a distribution | `OwnerProperties` | Owner values only |

Unsupported domains are unrepresentable rather than rejected at runtime. Reads
take the table's row type. Writes go through length-preserving
`*Mut` guards: whole-column replacement or removal, single-cell writes, and
transactional batches with ordered update semantics. A draft editor addresses
allocated slots, and writing a value to a deleted slot fails with `RemovedRow`.
System editors expose the same guard keyed by editing handle
(`TopologyEditor::properties_mut`), and a deleted handle fails with the
editor's typed invalid-handle error. `StructureEditor` does not dereference to its
topology draft: `topology()` reads it, `topology_mut()` reaches the edits that
need no coordination with geometry, and its own `properties_mut()` holds
structure-scoped values, so the two scopes never shadow each other.

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
Occupancy and B factors are typed site fields of structures and geometries,
with checked semantics independent of any generic property with the same name. Positions, cells,
weights, energies, time, steps, velocities, and forces retain their dedicated APIs.

Kekule uses one runtime unit system across all crates. Every public value with a
physical dimension, accepted or returned, is a `Quantity` in a compatible unit;
internal numerical state uses the library-wide canonical units. Mass and amount remain distinct dimensions, with molecular mass
and energy conventions chosen coherently. Unit composition and conversion reject
unrepresentable dimensions or scales; checked numerical storage rejects nonfinite
values. Exact units, tolerances, and conversion policies belong with
[units](crates/kekule/src/units.rs) and
[properties](crates/kekule/src/properties.rs), not in a second architecture table.

## Geometry and structure ownership

Geometries are flat. One `Geometry` holds all the coordinates and per-atom
state of one arrangement, and a `Structure` holds exactly one of them, so a
structure read from a restart file carries its velocities and a single-point
result carries its forces. A place in an order and a membership in a
distribution belong to the set and to the ensemble, never to a geometry; an
energy belongs to a state, not to any geometry, except as a geometry's own
property. Coordinate kernels share code through the borrowed `StructureView`,
not through a common owned payload.

A structure or structure set owns its shared topology once. A detached
geometry has atom rows only and acquires bond rows when a set or a frame
buffer binds it. Borrowed
`StructureView` access, obtained through the sealed `AsStructureView` (or a frame
buffer's `as_structure_view`), lets coordinate algorithms and writers share
kernels across structures, geometries, frames, and buffers. A view's `kind()`
is `StructureKind::Standalone` or `StructureKind::InSet` with a `SetContext`: the
geometry's place in its order and its velocities and forces. A frame of a codec or buffer reports its time and
step as its place, so writers and reductions never lose that state silently.
Structure writers hold a geometry to one policy, `geometry_metadata`, for all of
that state, whether it belongs to a set or a frame. Per-atom, per-bond, and
own properties of a structure, and the properties of its structure set, follow a
second policy, `properties`. Molfile and SDF, which store no occupancies or B
factors, put them under `geometry_metadata`. Each of these fails by default. Owned
projections (`to_structure`) explicitly materialize a structure while sharing the
topology. A set built from structures (`from_structures`) consumes them and
requires equal layouts (`same_layout`), keeping the first structure's
topology; there is no special single-geometry ownership model.

Read access is through views: `Topology` lookups return `Option` views
(`AtomView`, `BondView`, molecule, chain, residue, and atom-site views), and
`Structure::atom` joins an `AtomView` with that atom's geometry state. Mutation
goes through explicit, dimension-preserving setters on structures, detached
payloads, and the collection item guards returned by `get_mut()`, and through
`properties_mut()` on topologies (`topology_properties_mut()` on structures and
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
[structures](crates/kekule/src/structure/owner.rs),
[structure sets](crates/kekule/src/structure/structure_set.rs),
[measurements](crates/kekule/src/structure/measure.rs), and
[alignment](crates/kekule/src/alignment.rs).

Cartesian geometry edits mutate a structure's coordinates atomically while retaining
the exact topology, cell, and properties. A rigid motion of a whole structure
(`apply_rigid_transform` on a structure, geometry, or frame) rotates its cell with
its coordinates; selection edits such as `transform_atoms` keep the cell. Prepared distance, angle, and dihedral
edits bind references and moving selections to that layout, never to cached
coordinates. Automatic fragments must separate from fixed references; explicit
selections may deliberately deform boundary bonds. Neither path performs ring
closure, relaxation, or periodic imaging. See the
[geometry editing contract](crates/kekule/src/structure/geometry_edit.rs).

## Structure sets and ensembles

`StructureSet` is the one container for many geometries of one topology. What
makes a set a trajectory, a scan, or an ensemble is stated beside the
geometries, never carried by its type. Sources differ along two independent
axes and real ones combine them: an MCMC chain is ordered and a distribution,
a metadynamics run is time-ordered and weighted, a clustered trajectory is
both. A reader cannot know which it has — a multi-model PDB may be an NMR
bundle or an MD run, and a multi-frame XYZ may be conformers or an
optimization history — so a reader loads geometries and the caller states what
they mean, which is an annotation rather than a conversion.

The two axes divide on what the source recorded against what an analysis
infers.

An order is recorded, so it belongs to the set. An MD engine writes each
frame's time and step, an optimizer numbers its steps, and a scan driver
records its coordinate; trajectory files carry these values. An order, when
present, places every geometry: a time axis for a path produced by dynamics,
steps for an optimization, or values of one coordinate in a unit for a scan or
a reaction path (values need not increase, as a scan swept back and forth
shows). A time axis has times, steps, both, or neither, and each column covers
every geometry or none. Order values belong to the order, not to the
geometries: selecting keeps each selected geometry's value, and a removed
geometry takes its value with it.

A distribution is inferred, so it belongs to `Ensemble` and never to the set.
Reweighting a biased run, MBAR, clustering into states, and Boltzmann
populations from energies are interpretations, and one set often carries
several at once. An `Ensemble` owns one `StructureSet` and one `Distribution`
over it, validated together once; `Ensemble::into_parts` returns both without
copying. It lends the set but never lends it mutably, so everything that
changes which geometries exist goes through the ensemble, which remaps the
samples. Several readings of one set are several ensembles sharing every
coordinate.

A distribution covers every geometry of its set or a sorted subset, its
samples. Weights are stored as natural logarithms, one per sample, each
finite: weights from bias reweighting or importance sampling span hundreds of
orders of magnitude, and a stored linear weight would underflow to zero and be
rejected as nonpositive. Normalized weights and probabilities are derived on
request. Without weights the samples weigh equally. Log weights and the
partition are indexed by sample position. Every sample is in exactly one state
and every state has a representative sample; a state's population is the total
weight of its samples, derived and never stored. Only states have energies, on
one scale per distribution: a conformer's minimized energy or a computed free
energy, never a snapshot's instantaneous energy, which stays a property of its
geometry. Selecting and removing geometries remap the samples; a selection
keeps the energy and properties of the states it keeps whole.

A weight means something only against a target, so it is never a field of a
geometry or of a set. The same geometries weigh uniformly at 300 K and
differently reweighted to 310 K, and a metadynamics run has one set of weights
under c(t) reweighting and another under MBAR. What the engine recorded per
geometry — a bias potential, collective-variable values, an instantaneous
energy, a repeat count — is a geometry property, and the weights an analysis
uses are derived from those. `Ensemble::uniform` and `Ensemble::weighted` build
the one-state readings, so weighting a trajectory needs no partition.

A geometry is shared through `Arc` and copied on write, so regrouping,
reweighting, selecting, and collapsing copy no coordinates, and editing one
geometry copies that geometry alone. Set and ensemble operations fail with one
`StructureSetError`, whose `Missing` variant names the `Requirement` an
operation lacks, such as a time order or an energy for every state; a
geometry's own edits fail with `GeometryError`, and perception with
`TopologyPerceptionError`. Statistics and clustering over an ensemble observe
its samples with their probabilities and leave the other geometries out; over
a bare set they observe every geometry equally. Superposition and RMSD act on
every geometry.

Two axes are deliberately absent. Kekule records no runs or segments: replica
exchange, metadynamics walkers, umbrella windows and repeated runs are the
concern of the code that produced them, which knows its own boundaries, and a
pooled weighting across them is still one ensemble. `Order::Coordinate` holds
one coordinate, so a multi-dimensional scan uses `Order::Step` with the
coordinate values as geometry properties.

Where the geometries live is a type parameter. `StructureSet`
(`StructureSet<InMemory>`) holds them and lends them by reference, so its
accessors cannot fail. `StructureSet<IndexedFile>` holds none: it reads
geometry `i` from an indexed trajectory file as an owned `Geometry`, and a
read can fail. Both share the topology, properties, order, and the operations
on them; only a set in memory edits geometries or its order. A set of either
storage reads through the object-safe `StructureSetRead` trait (topology,
length, order, `peek`, `read`), which is open to the host. `peek` lends a
geometry the reader holds; `read` returns an owned one. Everything in kekule
that reads many geometries takes a `&dyn StructureSetRead`: the observation of
statistics, RMSD of listed geometries, pairwise fitted RMSD, and lazy
processing. kekule keeps no geometries a reader read, but for the geometry 0
that lazy processing checks its steps on, and starts no threads: caching,
read-ahead, and parallel observation belong to the host.

A reader describes geometries and their order and never a distribution, so
`Processed` and a host's own readers carry none. Processing neither removes
nor reorders geometries, so indices are preserved, and a host that wants a
processed ensemble reads the processed set and pairs it with the distribution
it already had.

A host loads chosen geometries of any reader with `StructureSet::from_reader`,
which takes indices (a range, a stride, or a list) in any order and makes a
faithful copy: each picked geometry whole, with the reader's order values. It
keeps no set properties. Every index is checked against the reader's length
first, and a repeated index is read once and its geometry shared. Geometries
picked out of ascending order keep their times, so `validate_monotonic_time`
refuses the set when a time decreases. The host changes what the geometries
mean with `set_order` and by building an `Ensemble`; `Ensemble::select`
restricts a distribution to the geometries it keeps. A host that needs
progress or cancellation reads geometry by geometry through
`StructureSetRead`.

These live in `kekule::structure`. `kekule::traj` owns file codecs, streaming
reader and writer contracts, the reusable `FrameBuffer`, the `IndexedFile`
storage of structure sets, periodic preprocessing, streaming superposition,
and lazy processing (`Processed`); it owns no second in-memory trajectory
type. A trajectory writer takes a `FrameView`: a borrowed `Geometry` with its
time and step. `FrameView::from_geometry` lends a geometry of a time-ordered
set and `FrameView::new` one held apart from a set, without copying arrays. A
`FrameBuffer` holds one `Geometry` bound to its topology plus the frame's time
and step; it lends its frame as a `FrameView` and becomes a `Geometry` with
its `Place` through `into_geometry`, again without copying. A frame the buffer
publishes replaces every field, occupancies and B factors included, so a
buffer never keeps a previous frame's state. `read_structure_set` and
`write_structure_set` move whole time-ordered sets.

## Trajectories and streaming

A trajectory, a structure set of a time order, represents one fixed-topology
epoch. Changing chemistry or hierarchy requires a new topology and explicitly
constructed geometry for the next epoch. Frame selection is separate from
atom subsetting: requested order, duplicates, and empty selections are
supported without renumbering stored time or step.

File readers receive the topology and interpret file atom index as its dense atom
index. They validate counts and available metadata; equal counts alone do not
prove atom identity. Callers may independently validate semantic atom order.
Because interpretation keeps source atom-row order and structure writers emit
dense order, a coordinate file written in a structure file's row order addresses
the same atoms through that file's topology.
An indexed reader verifies every frame before random access. One opened from
an exported frame index instead checks that the file begins with the bytes
the index's fingerprint covers, then verifies each frame as it decodes it
and checks that the frame ends where the index places the next and records
the time and step the index holds for it; an index never makes a reader
return a frame it has not verified. A file that only grew by appending is
the same file with more frames, and one whose existing bytes changed is not
(a DCD header's frame counts, which writers rewrite as they append, are not
content): refreshing a reader, or opening the grown file from an index,
verifies and indexes the appended frames without scanning the indexed ones
again. A refresh reads whatever file the path names when it runs, checking,
scanning, and fingerprinting that one file, so the index and its
fingerprint always describe the same bytes. A reader of a growing file
takes a final frame its writer has not finished as pending rather than as
corruption, reads it again at every refresh, and still rejects corrupt
data; sequential reads through an indexed reader end with its indexed
frames. A run restarted from a checkpoint continues its file through
`append_trajectory`, which keeps the frames the checkpoint covers, cuts the
file after them, and appends under the writer settings the file records, so
the continued file equals one written straight through; a reader following
the file sees the cut as an index mismatch and opens the file again.
Streaming buffers bind to one topology layout and publish complete frames
transactionally. A buffer hands its frame out as a geometry with its place,
and takes a geometry back, by moving it; the geometry stays bound and validated. Eager reads use the same decoder path through clean EOF.
Path writers stage output and publish only a completed nonempty trajectory;
stream writers write the same bytes into a host's `Write + Seek` stream and
leave placing the file to the host. `write_structure_set` writes a structure set
of a time order, or of none, through both, and
`read_structure_set` reads a file as a structure set of a time order. Codec writers append frames whole and in
order and never rewrite one, so after a flush the stream holds a complete
trajectory, a DCD header counting its frames, that another process can read
while it grows. Frame state a codec cannot carry fails
unless the caller sets `UnrepresentablePolicy::Omit` for its kind (frame
metadata, frame properties, or the cell), and the writer reports what it
left out. Occupancies and B factors, which no trajectory format stores, are
frame metadata of every codec, so `write_structure_set` returns the codec's
report as it is; collection properties are rejected rather than discarded. A codec
never invents frame state, except the times and steps XTC and TRR number
under an explicit `FrameTimingPolicy::Number` and DCD numbers on its header's
step sequence under an explicit `DcdFrameTimingPolicy::Number`, which reports
each frame step and time it replaces.

`open_indexed_structure_set` opens a file as a `StructureSet<IndexedFile>`,
`open_indexed_structure_set_from_index` opens it from a kept frame index, and
`IndexedFileTrajectoryReader::into_structure_set` makes one from a reader the
host opened, without opening or scanning the file again. The
set's order is the index's time axis, and it exposes the index and its
fingerprint. Geometry `i` is frame `i`, decoded on every read into a geometry
that holds the only copy of its coordinates: the codec decodes them once
into the arrays the geometry keeps. XTC decodes to single precision first,
and DCD with fixed atoms fills a frame's array with the fixed atoms before
its free ones. Reads take `&self` and run on any number of threads at once:
each read takes an idle indexed reader from a pool, or forks one from the
set's reader, and returns it, so no lock is held while decoding. Every
reader reads the file object the set opened, through positional reads, and
shares the set's frame index (`Arc`), so forking reads the start of the
file alone and every reader gives the same geometry. A read allocates the
geometry's arrays on the calling thread, and a pooled reader keeps none for
its next read, so the arrays are allocated and freed on one thread however
the host's threads take the readers. On Windows, reads of
the one shared handle run one after another, so threads take turns reading
the file and decode at once. A read decodes whatever bytes the file holds
at the indexed offsets, and kekule does not watch the file. A file
truncated while open fails every read of the frames it lost with one
`IndexMismatch` error naming the frame. Bytes rewritten in place are read
as they now are, without error when a frame keeps its size, until
`refresh` finds the indexed bytes changed. A file replaced at its path
keeps reading through the open file. `refresh` takes the frames a writer
appended, through the set's reader's refresh, and fails without change
otherwise; the pool's readers are kept and read by the new index. It opens
the path again, so a replacement that begins with the indexed bytes is
taken and every reader moves to it. A truncated or rewritten file fails
`refresh`, and the host opens it again. A host that must know of changes
between refreshes compares the fingerprint, or `TrajectoryFrameIndex::prefix_fingerprint`
of the frames it used, which samples the indexed bytes of those frames alone.
A host keeping the frame index, or a live clone, makes each growing refresh
copy the offsets. An `Ensemble` over the set refreshes through
`Ensemble::refresh`: a distribution over every geometry becomes one over the
geometries it covered, so appended frames are not samples;
`Distribution::extended_to` is the same operation,
`StructureSetDescription::extended_to` gives a saved description the longer
set's order, and `EnsembleDescription::extended_to` does both. A
clone of the set reads nothing: it shares the set's reader, file object and
index, starts with an empty pool, and forks its own reader at its first
refresh, so a host hands it to an analysis while the set goes on refreshing.
Refreshing either one changes only that one. A failed
read names the geometry (`StructureSetReadError::Trajectory`). `load` reads every
geometry, or listed ones, into a set in memory that keeps their order values;
an ensemble's `load` also restricts the distribution, as `select` does.

`Processed` reads another `StructureSetRead` and runs `ProcessingSteps` on each
geometry as it is read: make whole or wrap, center, fit onto a geometry of the
same reader or onto given positions, smooth over an odd window, and project
onto an atom subset, in that order. Each step is optional. Wrapping images
the molecules around the anchors and moves the anchors' centroid to the
cell center, so every molecule's center lies within the cell; centering on
the anchors then keeps them there. Geometries without a cell skip making
whole and wrapping, where the in-memory functions fail with `MissingCell`;
`Processed::skipped_periodic_steps` tells a host which geometries did.
Centering moves the centroid to half the sum of all three cell vectors, or
without a cell to where it lies in geometry 0, so a reader mixing geometries
with and without cells centers them on different targets. Construction
checks the steps and processes geometry 0, so steps that cannot process the
geometries fail there. Making whole and wrapping call the same kernels as the
in-memory periodic functions and `MoleculeImager`, and fitting the same as
`FrameSuperposer`, so each step gives the same bits. The steps change the
geometry in place: a geometry the reader below owns alone is not copied by
them, and a shared one is copied once. Smoothing averages the window's
geometries after the steps before it, each molecule imaged onto the geometry, so
it expects whole molecules; it keeps everything else of the geometry.
Analyses read the same steps without smoothing, through a second reader
with `ProcessingSteps::without_smoothing`. `Processed` keeps geometry 0, the
one construction processed, and lends it through `peek`; it keeps nothing
else it read, so a host caches in front of it.

Loaded and streaming transformations share kernels. Collection operations,
including superposition and periodic preprocessing, change the collection in
place and stage all affected state before publication; clone first to keep the
source. Superposition and RMSD are methods of structure sets, against one of
the set's geometries or an independent view; streaming superposition applies the
same kernel to a buffer. Streaming tools take the caller's frame index for
diagnostics. Superposition rotates cells, velocities, and forces consistently.
Periodic reconstruction uses asserted bonds and checks ring closure. Imaging
acts on whole molecules, around anchors or onto the same molecules of a
borrowed reference; reference imaging is per frame, has no temporal state, and
works in any frame order, using the frame's own cell. Unwrapping retains temporal state and requires ordered,
sufficiently close geometries; it precedes downsampling. Failed streaming operations
change neither the frame buffer nor temporal state. Diagnostics are opt-in.

RMSF and contact occupancy, of atom pairs or of pairs of atom groups, are
weighted statistics in `kekule::statistics`.
Every observation is a borrowed structure together with an explicit, finite,
positive weight, so no structure is reduced with an implied weight.
`Ensemble` methods observe each sample of the distribution with its
probability, its weight over the samples' total, and leave the other geometries
out; `StructureSet` methods observe every geometry with weight one,
taking them as equally weighted equilibrium geometries, and other weights are
passed to the accumulators explicitly. Streamed frames reach the same accumulators through their buffer's
`StructureView`; `kekule::traj` adds no statistics of its own. Accumulators consume
structures sharing the source topology's layout with memory bounded by
selected atoms, pairs, or groups, not observation count. Failed observations leave
accumulators unchanged. Accumulators of the same selection, or of the same
pairs of atoms or groups and cutoff, merge: a host observes ranges of
structures into partial accumulators on any threads and merges them in range
order, so a given split returns the same result on every executor, equal to
serial observation up to rounding. Merging refuses mismatched accumulators
with a typed error and, like observation, changes nothing when it fails.
`SampleObservations` lists the geometries of a reader to observe, every one
or listed indices such as a range with a stride, with weight one, or the
samples of an ensemble, every one or listed ones, with the distribution
restricted to them as `select` restricts it. An accumulator's `observe_samples` reads a range of those
observations in order, and `observe_samples_cancellable` stops on a
`Cancellation`; a set's or ensemble's own statistics run the same loop, so
one read from a file gives the results of one in memory. `alignment::geometry_rmsd`
and `alignment::geometry_aligned_rmsd` measure listed geometries of a reader
through the kernels of a set's RMSD.
Results retain atom/pair associations; observation
indices are diagnostic labels. Preprocessing, including alignment or imaging,
is explicit. These are separate analysis results, not new topology state or a
generic reduction framework.

Clustering is weight-aware and groups the samples of an ensemble's
distribution, or every geometry of a structure set, unweighted. It returns a
`Partition`, the grouping type a distribution owns, and never changes the
set. An ensemble's distances and partitions are indexed by sample position,
so `DistributionMut::regroup` installs a partition directly, and a structure
set's partition makes a distribution of every geometry with
`Distribution::new` and `Ensemble::new`; errors name geometry indices. `DistributionMut::regroup` applies a partition without changing any weight or sample
order: a new state that contains every member of an old state keeps the
energy and properties of its representative's old state, and a state split
apart loses them. Collapsing makes a new collection of each state's
representative carrying the state's population. Visiting order and
tie-breaking are deterministic. Pairwise distances
are a separate `DistanceMatrix` result. A hierarchical `Dendrogram` depends on
distances alone, and weights enter only when choosing representatives;
k-medoids instead minimizes the weight-summed distance to medoids, so weights
shape its clusters. Pairwise fits and leader clustering may run on a
`parallel::Executor`, poll a `Cancellation`, and report progress through a
`ClusteringRun`; fits split into contiguous ranges merged in a fixed order,
and leader clustering still places members one at a time, so every executor
and thread count returns the single-threaded result or error.
`clustering::pairwise_fitted_rmsd` fits the geometries of any `StructureSetRead`,
reading each geometry once and keeping only its fitted atoms, through the
loop of a set's pairwise fits; clustering a matrix works for sets of either
storage. See [clustering](crates/kekule/src/clustering.rs),
[hierarchical clustering](crates/kekule/src/clustering/hierarchical.rs), and
[k-medoids](crates/kekule/src/clustering/medoids.rs).

DSSP assignment reads one structure and returns a snapshot result. A
`DsspPlan` prepares what it reads from the topology once and binds that
layout, so the frames of a trajectory reuse it; it returns what `dssp::assign`
returns for the same structure, errors included. Every stage runs on a
`parallel::Executor` in contiguous ranges merged in DSSP order, so every
executor and thread count returns the single-threaded result. See
[DSSP](crates/kekule/src/dssp/mod.rs).

Non-covalent interactions split by what they depend on. `InteractionFeatures`
classifies one topology's donors, acceptors, charged groups, aromatic rings,
halogen-bond donors, metals, metal ligands, and hydrophobic atoms once, per
reusable definition, and is reused for every structure sharing that layout. It
reads installed perception, or perceives a private copy of a definition without
it, and never changes the topology; it is a separate derived result, not
perception state. A search reads one structure, or two, each side with features
of its own topology, and returns a snapshot of interactions addressed by dense
indices of each side's topology. Criteria are a plain options value in
`Quantity` units with literature defaults. Within one structure, bond separation
excludes ends of one molecule; across structures no bond joins the sides, so
the two searches agree on the parts of a structure no bond joins. Both reduce
hydrophobic contacts as PLIP does, after every pair is found: none between
π-stacked rings, reported or not, then the closest per atom and residue in each
direction. The side with fewer atoms lists its features; the other reads only
atoms within reach of them, and pairs are found on the uniform point grid that
group contact occupancy also uses. Hydrogen-bond occupancy is a weighted,
mergeable accumulator of found interactions, with the transactional observation
and merge rules of the weighted statistics. See
[interactions](crates/kekule/src/interactions/mod.rs).

For format profiles, limits, periodic conventions, statistical definitions, and
streaming sequence/reset contracts, see
[structure sets](crates/kekule/src/structure/structure_set.rs),
[collection alignment](crates/kekule/src/alignment/collection.rs),
[frame buffers](crates/kekule/src/traj/trajectory/buffer.rs),
[file I/O](crates/kekule/src/traj/io/mod.rs),
[structure sets read from files](crates/kekule/src/traj/io/indexed_file.rs),
[lazy processing](crates/kekule/src/traj/processed.rs),
[structure set readers](crates/kekule/src/structure/structure_set/read.rs),
[periodic operations](crates/kekule/src/traj/periodic.rs),
[streaming periodic operations](crates/kekule/src/traj/periodic/stream.rs), and
[weighted statistics](crates/kekule/src/statistics.rs).

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

mmCIF, PDB, and XYZ interpretation gives identical molecules one definition, so
instances map atom for atom and perception runs once per distinct molecule: a
crystal's waters share one definition. Two molecules are identical when, in
local atom order, their atoms (element, isotope, formal charge, radical,
implicit hydrogens, atom map), bonds, represented stereo, and properties are
equal, and so are the labels the format gives each atom: component and atom
name for mmCIF and PDB, none for XYZ. Classification reads those labels, so
requiring them equal keeps one class valid for every instance. Hierarchy stays
per instance, and editors clone a shared definition when an edit touches one
instance, so the others keep it unchanged. Classes are per instance too: an
instance whose class a hierarchy edit cleared is published in a definition of
its own, which shares the molecule, and classified alone, so the other
instances keep their class. Appending a structure or topology to a draft shares
the same way: an appended molecule identical to a definition the draft holds,
from its source or an earlier append, becomes an instance of it, when its
perception and class are equal too. Appending is where copies enter a draft;
publication does not search for identical definitions, so it still costs what
an edit changes. See
[definition reuse](crates/kekule/src/io/definition_reuse.rs) and
[sharing on append](crates/kekule/src/topology/editor/sharing.rs).

Interpretation never silently discards source content. Content the interpreter
does not support is rejected by default; an explicit option downgrades it to typed
report entries that name the record and the kind of content omitted. A deliberately
lossy projection, such as base-only SMILES, says so in its name.

Every stage reports through one diagnostics model. Parse, interpretation, and
write failures of every format are one `diagnostics::FormatError` carrying a
structured kind and a source location relative to the original input, so
callers can match on failures and map them to the source without
format-specific arithmetic. Non-fatal report issues carry a structured kind and
identify the source records they concern.

Coordinate-only trajectory formats are the deliberate exception: they carry no
chemistry to interpret, so streaming codecs decode frames directly against a
caller-supplied topology (see Trajectories and streaming). XYZ is read both
ways: `kekule::xyz` interprets its frames into structures through this
pipeline, and `kekule::traj` streams its coordinates against a topology, such
as the one that interpretation produced. Writing XYZ keeps elements and
coordinates only, so the writer reports the topology content it leaves out in
the shared write report. It writes a structure set as one frame per geometry,
with the PDB writer's rules for set state: a place in the order, a time
order without times or steps, velocities, and forces fail or are omitted with
a report under `geometry_metadata`. XYZ interprets frames as a structure set.

| Format scope | Interpretation and projection |
| --- | --- |
| SMILES record | All connected components, with topology projection |
| Molfile | One structure with topology and molecule projections |
| SDF record | Structure, title, fields, and reports; each record independent |
| CIF data block | One structure, or compatible coordinate models as a structure set |
| PDB entry | One structure, or compatible coordinate models as a structure set |
| XYZ file | One frame as a structure, or frames of the same atoms and inferred bonds as a structure set |

Sibling SDF records and mmCIF blocks never merge implicitly. Singular projections
check cardinality; plural molecule projections retain every connected component
in source order. MOL/SDF interpretation supplies deterministic synthetic hierarchy
at topology scope. In a conformation the coordinates are the configuration:
wedges and atom parities that contradict them yield to them with a warning,
while flat drawings take configurations from those annotations. Coordinates,
flat or not, give E/Z only to double bonds that are eligible before symmetry
analysis, so an aromatic bond or one whose endpoint bears two identical
terminal ligands has none. By the same test, SMILES directional marks assert
nothing at such a bond: interpretation omits the configuration and its report
names the bond, so a parsed molecule holds no E/Z that a conformation cannot
carry. mmCIF preserves label/author hierarchy identities and interprets
all coordinate-model candidates against the same identity, chemistry, and dense
layout. PDB carries author identities only: a residue is a run of consecutive
records, and a number that recurs after other residues names a new residue.
Malformed rows in unselected models still fail validation. Structure set assembly
releases redundant topologies while transferring structure payloads; multiple
coordinate models do not imply temporal trajectory semantics.

mmCIF alternate rows remain in the source document. Interpretation selects whole
residue alternatives or explicitly constrained groups and records its choices in
the format report. Occupancy ranks representatives, not whole-structure
probabilities. Assembling coordinate models into a structure set never expands
alternate labels; caller-supplied conformation selections must pass the same identity and
topology checks. Every coordinate model becomes a geometry; the formats state
no distribution, so a host that reads the models as draws builds an
`Ensemble`, and occupancy is never converted into a weight.

Macromolecular formats (mmCIF and PDB) do not carry complete bonds, and molecule
boundaries follow bonds, so one resolver shared by these formats derives connectivity
from ranked evidence:

1. Explicit source bonds: component bond tables in the source, and covalent
   and disulfide links between residues. Records that list
   connectivity without bond types, such as PDB `CONECT`, supplement a
   residue's template and replace distance inference for residues without one.
2. Component templates supplied as interpretation context, keyed by component
   identifier and matched to source atoms by atom name and element.
3. Polymer linkage between residues adjacent in sequence, and the links of a
   branched entity between residues of one copy of it, accepted only within
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

XYZ states elements and coordinates only. Its interpretation applies the
same distance rule to the whole frame: every bond is inferred and single, or
dative between a metal and a non-metal, every atom fixes zero implicit
hydrogens because the file lists all atoms, and nothing assigns bond orders
or charges. Each frame infers its own bonds, and the frames of a structure set must
agree on them. The report lists every bond with its distance and cutoff,
fallback radii, and contacts too close to bond. A host that holds atoms as
elements and positions builds the same document from them and has it
interpreted the same way; its atoms and frames are located by their position
rather than by a line.

Metal coordination between residues is not a bond. A metal ion, the ligand
it holds, and the protein that coordinates it are separate molecules, since
covalent connectedness alone decides molecule boundaries; the report keeps
each such link as source connection metadata, and consumers perceive
coordination geometrically as an interaction. mmCIF `metalc` links are
coordination, and so are links whose records state no type, PDB `LINK` and
`CONECT`, between a metal and a non-metal. Within one residue, coordination
is a `Dative` bond, so a metal complex that is one component, such as heme,
stays one molecule. A source component bond table states coordination with
its own flag; where the source leaves the flag unstated, as wwPDB entry files
do, the flag comes from the component template, as charges do, and without
one a bond between a metal and a non-metal is coordination. An mmCIF block
and its PDB rendering therefore give the same bonds. XYZ, Molfile, SDF, and
SMILES have no residues: their dative bonds are bonds. Hydrogen bonds, salt bridges, and links to symmetry mates are never
bonds; unresolvable links are reported. The hydrogen bonds and salt bridges a
file lists are annotations, often thousands of them, so the report counts them
once per connection type. A heavy atom whose name does not match
its residue template is an error, not a silently separate molecule. The report
records the evidence for every bond and names residues whose resolved atoms
split into several molecules.

Formal charges come from the source when stated, otherwise from the template,
including when the source's own bond table supplies the bonds. Without a
template, an atom with more bonds than a neutral atom allows takes the charge
that makes its valence valid, +1 before -1, as the methylated N7 of
7-methylguanosine takes +1. Hydrogens present
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
string and one `write_to(writer, source, options)` sink writer per format. Sink
writers stream into the caller's writer without buffering it, as trajectory
codecs do; callers wrap unbuffered sinks, and only path writers buffer. The
exception is the Molfile/SDF writer, which renders one CTAB before emitting it,
because automatic version selection learns whether V2000 can hold a record only
by rendering it; SDF still streams record by record. Each source type
is a struct with private fields, constructors for the inputs it accepts, and
the record's title or name (`with_title`, `with_name`), and it converts from
the inputs that format can represent: molecules,
structures, set geometries, interpreted records, or mmCIF block sources carrying
classifications or reports. Structures and set geometries share a borrowed
structure path. Coordinate-bearing writers emit atoms in dense order.
Unsupported representational content must fail explicitly. Writers read the
view's `kind()`: set or frame state the format cannot carry fails unless
the caller sets `UnrepresentablePolicy::Omit` for it, and a format that can
carry it writes it whichever path supplied the view. Each kind of
unrepresentable content has its own policy, so omitting one never omits
another. The PDB, mmCIF, Molfile, and SDF writers return a `WriteReport`
with the text, or from `write_to`, that names every omission and the number
of entities it affects. There is no universal save trait or implicit trajectory-to-structure
export.

- SMILES exports represented chemistry and supported stereo. Canonical output
  requires complete labeling under explicit resource bounds; exhaustion is an
  error, never an unproved canonical answer. Hydrogen and isotope projection must
  agree with emitted syntax and preserve stereo.
- Molfile chooses V2000 when it can represent all supported content, otherwise
  promotes to V3000; explicitly requested versions fail on unsupported content.
  Coordinate stereo must agree with the actual rounded emitted geometry. Writers
  must not invent assertions from an unasserted drawing or hide a conflict. A
  Molfile or SDF record cannot carry a periodic cell, so writing a structure with
  one fails unless the caller chooses to omit it.
- SDF represents independent records. A geometry writes its place in a step
  or coordinate order as data fields. Set properties, geometry properties,
  occupancies, and B factors fail unless omitted. Only the explicit
  `sdf::structure_set` assembles records into a structure set, whose order
  the steps or coordinates give. A coordinate unit travels as the symbol of
  a unit constant, which both the writer and the reader take from one table
  of `units`. A distribution is never a field of a record: an ensemble is
  inferred, and native files keep it. A Molfile, an
  mmCIF block, and a PDB entry cannot carry a place, and PDB writes a step or
  coordinate order as one `MODEL` per geometry once places are omitted. A
  Molfile, an SDF record, and an mmCIF block cannot carry a frame's time,
  step, velocities, or forces. mmCIF distinguishes independent structure
  blocks from one structure set's multi-model block. mmCIF writes a cell periodic along all three
  axes in the standard crystallographic frame as `_cell` in space group P 1,
  and reads `_cell` as the structure's cell only in P 1 (or with no space group),
  where the coordinates are a whole unit cell; PDB reads `CRYST1` by the same
  rule. PDB's 1 Å placeholder cube is no cell, and stated cell parameters that
  are not the structure's cell are reported as omitted. Biological assemblies are
  report data: mmCIF `_pdbx_struct_assembly`, `_pdbx_struct_assembly_gen` and
  `_pdbx_struct_oper_list`, and PDB `REMARK 350` biomolecules, which the PDB
  parser types while keeping their lines verbatim. Each assembly lists its
  chains by label and the rigid operations that place their copies, products
  of operators already composed; one that cannot be read is reported as
  omitted. Interpretation keeps the deposited coordinates, and
  `Assembly::to_structure` builds the copies on request: identity copies keep
  their chain identifiers, others take the operation's name after a dash, as
  wwPDB assembly files name them, and copies share their molecules'
  definitions. See [assemblies](crates/kekule/src/io/structure_interpret/assemblies.rs).
  Canonical classification informs entity
  kinds; hierarchy alone does not establish polymer status or contextual roles.
  Expert overrides and source provenance handle finer format distinctions.
  The mmCIF writer, like the PDB writer, states only the bonds its reader
  cannot derive, as the PDB's own files do: a residue whose bundled
  template, after the protonation variant its hydrogens select, gives
  exactly its bonds, and peptide and phosphodiester links within bonding
  distance in the written coordinates, are left out; links between residues,
  dative bonds, and the bonds of residues without such a template are
  `_struct_conn` rows. A residue is either derived or stated whole, since a
  stated internal bond makes the reader take the residue's bonds from the
  source alone. Bond order within a molecule follows the reader's evidence,
  which mmCIF does not carry. Metal coordination between residues lives in
  the interpretation report, not the structure, so the writer states it as
  `metalc` rows only when it is given the report, which must still match the
  structure atom for atom; an edited structure is written without it, and a record
  whose partner the structure lacks is left out and reported. A structure's own
  dative bond between residues is written as coordination (`metalc`, or
  PDB `CONECT`), which reads back as a record rather than a bond, and the
  write report counts it. Both writers take the coordination records that
  still match the structure from one shared step: each partner must keep the
  atom names, type symbol, and component the report recorded.
  Every entity instance is its own label asym: a chain the source names as
  one asym stays one, while a chain without a source asym that holds a
  polymer, ligands, and water, as a PDB author chain does, is written as
  several asyms that keep its author chain identifier. Chains that share a
  label, as appended copies of one structure do, are distinct instances: later
  ones are written under the label with a copy serial (`A-2`) and their
  author names, and read back with those labels. The hierarchy keeps
  labels as appended, since a label names a source's instance, not one
  of the combined structure's.
- PDB writes one entry per call, with `MODEL` records only for a structure
  set. It writes what its reader reads back and no more: author
  identifiers, a `TER` record after each chain's polymer, and `CONECT` for
  every bond the reader cannot derive. Standard residues covered by the
  bundled component dictionary, polymer linkage, and water are derived and
  not listed. The dictionary applies to its residues whatever `CONECT`
  lists, so a covered residue whose bonds differ from every variant of its
  template is unrepresentable, as are bonds `CONECT` cannot state.
  Identifiers that do not fit their columns, and residues whose atoms are
  not consecutive in dense order, fail, except a chain identifier of up to
  four characters, such as a chain solvation adds, which an explicit policy
  writes in the segment identifier columns with a blank chain column, as
  CHARMM and NAMD name segments; the PDB reader takes such a segment for the
  chain, so the identifier survives. Hybrid-36 extends serial and residue
  numbers. A cell is written as `CRYST1` under the rule its reader applies.
  Given interpretation reports, the writer states their metal coordination
  as `LINK` records before the coordinates, with the atom records'
  identifiers and the symmetry operator `1555`. `LINK` records belong to
  the entry, not to a model, so a structure set states the coordination every
  geometry's report states. The rest are left out and reported, and so are
  records naming a chain written as a segment identifier, records between
  two metals, which readers would take for a bond, and records whose
  partner identifiers also name another atom.
  Solvation names its waters and ions as the Chemical Component Dictionary
  does, which PDB and mmCIF readers check, and leaves their (zero)
  implicit-hydrogen counts to valence inference, as those readers do, so a
  solvated system writes to both formats and reads back as written.

Exact capability and numerical policies live with
[SMILES writing](crates/kekule/src/io/smiles/write.rs),
[canonical labeling](crates/kekule/src/io/smiles/canonical.rs),
[Molfile writing](crates/kekule/src/io/molfile_write.rs),
[SDF writing](crates/kekule/src/io/sdf_write.rs), and
[mmCIF writing](crates/kekule/src/io/mmcif_write.rs), and
[PDB writing](crates/kekule/src/io/pdb_write.rs). Writer failures say why
no output was produced: the format cannot encode the content, the source lacks
required perception, a work limit was reached, the source is inconsistent, or
the sink failed.

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
Pairing by hierarchy labels reads one namespace, author or label identifiers, as
the selection language's terms of that namespace do, reports the selected atoms
it leaves unpaired, and refuses a label naming more than one selected atom on a
side rather than choosing one.
See [alignment correspondence](crates/kekule/src/alignment/correspondence.rs) and
[collection alignment](crates/kekule/src/alignment/collection.rs).

Symmetry-corrected alignment pairs one selection under every permutation of
equivalent atoms and keeps the smallest RMSD. The permutations are prepared once
per selection and layout by an exact, bounded search for graph automorphisms that
fix the selection setwise and never exchange molecule instances. Equivalence uses
element, isotope, hydrogen count, and connectivity, not bond orders or charges,
so localized resonance forms stay equivalent. Exceeding a bound is an error, never
a truncated result. See [symmetry permutations](crates/kekule/src/alignment/symmetry.rs).

Reconstruction validates represented graph connectedness first, then property
references/dimensions, then installs checked perception. Topology reconstruction
restores definitions, instances, classes, hierarchy, qualified IDs, and dense
layout consistently before validating structure payloads. Disconnected persisted
graphs must be partitioned or rejected. Export to scientific formats is not exact
native persistence and must not weaken these boundaries. See
[reconstruction regressions](crates/kekule/tests/molecule/identity.rs).

Exact persistence is the `native` format module, the one owner of Kekule's
binary encoding of molecules, topologies, structures, structure sets, and
ensembles. A structure or structure set may carry its structure interpretation reports, stored
as they are, so that an application that saves structures in native files
keeps what the source stated beyond the structure. A structure set record may
also keep its geometries elsewhere and store only their number, with its
topology and order, as a `StructureSetDescription`, and an ensemble record
likewise with its distribution, as an `EnsembleDescription`. It stores
each owner's authoritative state, including
explicit or inferred class assignments and adjacency order, and reconstructs
it only through those owners' validating paths: molecule publication,
topology publication with stored classes restored instead of inferred, and
the structure constructors. Stored perception is reinstalled through
`Molecule::install_perception` or discarded. A decoded topology keeps every
semantic ID but is a new publication with a fresh layout identity. Decoding is
bounded before allocation and never trusts the input. Semantic IDs have a
stable encoding (`raw`, and `to_bits` for instance-qualified IDs) so
applications can persist references next to native files. See the
[format specification](docs/native-format.md) and
[topology reconstruction](crates/kekule/src/topology/reconstruct.rs).

## API and maintenance rules

The public API of every crate follows the rules below. Code that predates them
is being migrated for 0.4.0 and is not a precedent; new and changed code follows
them.

### Publication and module layout

Implementation modules are private. Each public module is an explicit list of
re-exports: no glob re-exports and no public implementation files. Every public
item has exactly one public path; `prelude` is the only module that repeats
items. An item that no public path reaches is `pub(crate)`. Workspace lints
(`unreachable_pub`, `missing_docs`, `missing_debug_implementations`) enforce
this.

The public modules of `kekule`, by domain:

| Domain | Modules |
| --- | --- |
| Foundations | `units`, `math`, `properties`, `diagnostics`, `parallel` |
| Chemistry | `molecule`, `stereo` (`cip`, `coordinates`, `cleanup`), `perception` (`valence`, `rings`, `aromaticity`, `conjugation`, `resonance`), `hydrogens`, `descriptors` |
| Queries | `query`, `substructure`, `select` |
| Systems and structures | `topology` (`edit`), `structure` (`measure`, `placement`), `solvation` |
| Analysis | `alignment`, `clustering`, `dssp`, `interactions` |
| Formats | `smiles`, `smarts`, `molfile`, `sdf`, `mmcif`, `pdb`, `xyz`, `macromolecular`, `components`, `native` |
| Trajectories | `traj` (`io` with `xyz`, `dcd`, `trr`, `xtc`; `periodic`; `analysis`) |

A type lives in the module of the domain that owns it:

- perception state with its perception step;
- represented stereo in `stereo`;
- descriptors and other derived per-molecule results in `descriptors`;
- interpretation types shared by mmCIF and PDB in `macromolecular` only;
- the shared failure vocabulary in `diagnostics`.

Format modules stay at the crate root, except the trajectory codecs, which
live in `traj::io` beside the streaming contracts they implement. The
foundational crate has no optional features.

### Names, conversions, and identifiers

- **Conversions.** Use `as_*` for cheap borrowed views, `to_*` for owned
  conversion while retaining the source, and `into_*` for consuming
  conversion. Consuming does not promise zero allocation, but should transfer
  compatible storage.
- **Ownership.** Constructors that take ownership of chemistry accept values
  (`Molecule`, `Structure`), and callers clone explicitly.
- **Getters.** No `get_` prefix. Counts are `*_count()`; `len()` belongs only
  to collection types.
- **Molecule vocabulary.** An *instance* is an occurrence of a molecule in a
  topology and a *definition* is its reusable chemistry. A method named
  `molecule()` returns a `Molecule`.
- **Properties.** Use full `properties` names.
- **Units.** No public name carries a unit (`_angstrom`, `_degrees`).
- **Identifiers.** Semantic IDs and dense indices are distinct types, and
  dense atom and bond indices are `TopologyAtomIndex` and `TopologyBondIndex`.
  Frame indices are `usize`; `u64` is used only for byte offsets.
  Identifiers display in one style.

### Entry points

Every operation has exactly one home.

- **Lifecycle operations** of an owner are its methods: perception and
  hydrogen transforms on `Molecule`, drafts on owners, interpretation on
  format documents, projections on interpretations.
- **Expert algorithms** are free functions in their domain module.
- **Formats** provide parsing and writing as free functions of the format
  module.
- **Options** are always an explicit parameter, and `Default` reproduces the
  documented default behaviour. There are no `*_with_options` twins.
- **Shortcuts.** The only one-call shortcuts that compose several steps are
  `smiles::to_molecules`, `smiles::to_topology`, and `Molecule::perceive`.
- **Reusable work** is prepared explicitly: substructure search, for example,
  runs only on prepared targets.
- **Reads** live once on their owner and are reached through `Deref` or a
  view, never through forwarding wrappers.
- **Mutation** stays on editors, `*_mut` guards, and semantic owner APIs.
- **Avoid** redundant owner wrappers, compatibility aliases, and generic
  target or remapping abstractions.

### Options, results, and enums

- **Options structs** have private fields, `Default`, `#[must_use]` chainable
  `with_*` setters that validate their argument, and getters. Named presets
  are constructors such as `SmilesWriteOptions::canonical()`.
- **Results and reports** have private fields with accessors, have no
  `Default`, and return collections as slices or `ExactSizeIterator`s.
- **Public enums** are `#[non_exhaustive]` unless their documentation
  declares them closed, which is reserved for variant sets fixed by
  definition: orientations, sides, byte order, base dimensions.

### Lookups, errors, and panics

**Lookups.**
- Identity lookups return `Option`.
- Searches by name or label return `Result`, with typed not-found and
  ambiguous cases.
- No API returns `Result<Option<_>>`.

**Errors.**
- Each functional area has one error type with typed payloads, never
  `&'static str` reasons or stringified causes.
- `source()` returns the underlying cause, and `Display` does not repeat it.
  A wrapper keeps its cause as a typed value, never as text. A variant that
  adds nothing to its cause is transparent: it displays the cause and returns
  the cause's own source.
- `diagnostics::ErrorChain` writes an error and its source chain as
  consecutive sentences: `Cannot build the solvated structure. Atom 12 overlaps
  atom 3.`
- Every format's parse, interpretation, and write failures, SMARTS included,
  are one `diagnostics::FormatError`. It carries the stage, a
  `DiagnosticKind`, the format, the source location, and the record.
- Resource exhaustion is reported everywhere as
  `diagnostics::LimitExceeded`, naming the resource, the limit, and the
  observed value. The error of the operation names it in its own sentence:
  `Ring perception exceeded its work limit: 69,400,898 steps of 69,248,000.`
  Default limits are exported as constants.
- Cancellation is cooperative. A long-running operation accepts an optional
  `diagnostics::Cancellation`, through its options or a `*_cancellable`
  lifecycle method when it has none. It polls the cancellation between
  stages and inside long loops, and fails with `diagnostics::Cancelled`, or
  `DiagnosticKind::Cancelled` in a `FormatError`. The error keeps `Cancelled`
  as its source: `Resonance perception stopped. The operation was
  cancelled.` It never returns a partial result. Without a cancellation, it
  does not poll.
- A recoverable draft (`try_build`, `try_finish`) fails with
  `DraftError<D, E>`, which returns the draft with the error. It is the only
  generic error wrapper.

**Error messages.** Applications show `Display` texts to their users.
- A message is one or two short plain sentences, each with a capital letter
  and a full stop. The fact comes first: what failed and why, then the one
  action if there is an obvious one. `Atom 7 has no coordinates. Add them
  before measuring.`
- It names its objects: atoms, bonds, residues, records, lines, files, and
  numbers with their units. Text taken from the input is in double quotes.
- It states facts, not mechanisms. No relative clause explains how Kekule
  works, no rule follows the fact, nothing hedges, and no Rust API is named.
  The variant's documentation explains the cause and names the API that
  avoids it.
- A wrapper's sentence says what failed at its level, `Cannot build the
  solvated structure.`, so that its cause reads as the next sentence.
- A `FormatError`'s first sentence names the stage, the format, and the
  place, `Cannot parse the PDB input at line 12.`, and its message follows.
- Identifiers display as lowercase noun phrases, such as `atom 3`, for use
  inside sentences.

**Panics.** Public functions do not panic on any input their types allow.
Internal invariants may use `expect`, and the message names the invariant.

### Units, topology handles, and traits

**Units.** Every public value with a physical dimension, input or output, is a
`Quantity` (see Properties and units).

**Topology handles.**
- Owners and views expose `topology() -> &Topology` and
  `shared_topology() -> Arc<Topology>`.
- Arguments take `&Topology` when only layout identity is checked, and
  `impl Into<Arc<Topology>>` when the topology is stored.
- Layout compatibility uses `shares_layout` or `same_layout`, never pointer
  identity.

**Traits.** Traits that are not extension points are sealed. An extension
trait documents the contract its implementors must keep.

### Maintenance

When adding state, choose its owner from the ownership map. Substantial derived
results deserve their own concrete types. Format metadata stays with the format
unless deliberate promotion defines its canonical meaning and validity scope.
Keep API inventories and algorithm-specific policies beside the implementation;
update this document when ownership or cross-module invariants change.

### Companion packages

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
failure is one typed error with a stable kind, the failing molecule definition,
and every atom it concerns.

`kekule-potentials` owns energy evaluation and geometry optimization; the
foundational crate owns no potential. A prepared `Potential` binds one topology
layout, does not change during evaluation, and accepts any `StructureView` sharing
that layout, including perceived snapshots; an independently equal topology is
incompatible. Energies and
gradients are separate result objects in canonical units, not geometry properties.
Backends lower explicit parameters, currently an OpenFF `ParameterizedTopology`,
into private functional-form kernels. Parameter assignment stays in its own
companion crate. Unsupported geometry state, such as a periodic cell, is
rejected rather than ignored. Gradients are exact derivatives and never capped;
coordinates fail only where a requested energy or gradient is undefined.
Minimization reads a borrowed view and returns new positions without changing
chemistry, topology, cell, or properties. A minimization scope moves a selection
of atoms; fixed atoms keep their input coordinates bit for bit and still enter
every evaluation, and optional position restraints add to the minimized energy
without becoming part of a potential. See the
[potential contract](crates/kekule-potentials/src/lib.rs) and
[OpenFF backend](crates/kekule-potentials/src/openff.rs).
