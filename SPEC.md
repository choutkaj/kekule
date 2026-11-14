# Specification

This document defines Kekule's ownership boundaries, cross-module invariants,
and API rules: what each owner is responsible for and what must hold between
owners. How an operation works — its options, limits, numerical policies, and
per-format details — belongs in Rustdoc beside the implementation; the links
below point to it, and `cargo doc --workspace --all-features --no-deps --locked`
builds it. The [capability index](docs/capabilities.md) lists supported features
and their limits.

Change this document in the same change that alters ownership, vocabulary, or a
cross-module rule.

## Ownership map

| Owner | Authoritative responsibility | Must not own |
| --- | --- | --- |
| `Molecule` | One nonempty connected `Graph`, derived `Perception`, `MoleculeProperties` | Coordinates, hierarchy, system classification |
| `Topology` | Reusable definitions and their complete molecule instances, qualified identities, dense atom and bond order, one `Hierarchy`, classification, `TopologyProperties`, layout identity | Geometry or bonds between separate instances |
| `Hierarchy` | Chains, residues, and atom sites over live instance atoms, with source labels and author identifiers | Atoms, bonds, or molecule boundaries |
| `Geometry` | The one owned coordinate payload, shared and copied on write: one arrangement of a topology's atoms in dense atom order, with positions, optional cell, occupancies, B factors, velocities, forces, and `GeometryProperties` | A topology; a place in an order; a weight |
| `Structure` | One shared `Topology` and one `Geometry` bound to it | A second chemical or hierarchy model; a place in an order |
| `StructureSet` | One shared topology, set properties, `Geometry` values in a stable order, and an optional `Order` with one value per geometry | An order value on a geometry; a weight, a state, or a distribution |
| `Ensemble` | One `StructureSet` and one `Distribution` over its geometries | Geometries of its own; an energy on a geometry; a stored population |
| `Distribution` | `Samples` (every geometry of the set or a sorted subset), one log weight per sample, one `Partition` of the samples into states, one `StateData` per state, and a temperature | Geometries; per-geometry energies |
| `FrameBuffer` | One `Geometry` bound to its topology plus the frame's time and step, reused across the frames of a stream | A second in-memory trajectory type |

### Vocabulary

Kekule names what a system is made of with these words and uses no others for
these concepts.

A **definition** is a molecule's reusable chemistry stored once in a topology,
and an **instance** is one placement of a definition in a topology. A method
named `molecule()` returns a `Molecule`.

A **geometry** is where the atoms of one topology are, plus the per-atom state
recorded at that arrangement: positions in dense atom order, an optional
periodic cell, occupancies, B factors, velocities, forces, and geometry
properties. A geometry carries no topology; it is *bound* to one when it has a
property row per atom and bond of that topology.

A **structure** is a topology and one geometry of it. The topology stays
coordinate-free and may be shared with other structures, sets, frames,
selections, and prepared potentials.

A **structure set** is one shared topology and many geometries of it, in a
stable order, with an optional `Order` saying what that order means.

An **ensemble** is a structure set read as draws from a distribution: which of
its geometries are samples, what each weighs, and how they group into states.

A **frame** is a geometry as a trajectory file stores it, with its time and
step. Frames belong to `kekule::traj` alone.

**Sample**, **partition**, and **state** belong to the statistical reading
alone. A sample is a geometry an ensemble's distribution covers (`Samples`).
A partition groups the samples into states (`Partition`). A state is one such
group with one representative sample; what a distribution stores for it, an
energy and properties, is its `StateData`, and `StateView` reads a state whole:
its samples, representative, data, and population.

`Model` never names a structure. It stays for perception models
(`AromaticityModel`, `ValenceModel`), physical and machine-learning models
(`WaterModel`, NAGL bundles), and the `MODEL` records of PDB and mmCIF, which
this document calls coordinate models.

A salt, solvent box, or protein-ligand complex is a topology containing
connected molecules. Covalent connectedness determines molecule boundaries;
hierarchy and noncovalent interactions do not. Adding a bond between instances
constructs a new connected molecule and publishes a new topology.

Dense numerical containers carry values and units, not atom identity or
topology handles. The owning topology or structure translates semantic IDs to
dense indices and validates every associated array and property column.

See the [crate overview](crates/kekule/src/lib.rs),
[molecular owner](crates/kekule/src/molecule/owner.rs),
[topology types](crates/kekule/src/topology/mod.rs),
[structure types](crates/kekule/src/structure/mod.rs), and
[structure sets](crates/kekule/src/structure/structure_set.rs).

## Represented and derived chemistry

`Graph` is authoritative atom, bond, connectivity, and represented stereo state.
Local `AtomId` and `BondId` identify entities within one molecule. Published IDs
are dense in every ID space; drafts may hold deleted slots, and publication
renumbers survivors in draft order and reports the correspondence. Atom
chemistry includes element, isotope, charge, radical, hydrogen declaration, and
atom map; annotations use properties instead of extending every atom or bond
with a map. Bonds carry localized orders. Aromaticity is perceived state, not a
bond order. Source aromatic and stereo syntax is normalized during
interpretation.

Represented stereo must refer to valid focuses and adjacent carriers, with the
required focus bond order and consistent assertions and groups, and an
implicit-hydrogen carrier must agree with a fixed hydrogen declaration. One
validator checks these conditions for explicit stereo edits, publication, and
`MoleculeEditor::validate`, before canonicalization; graph-edit pruning checks
references alone, since hydrogen transforms change a carrier and its
declaration in separate steps. An edit that invalidates stereo — changing a
focus order, adding a bond to a focus, deleting a carrier bond — prunes the
affected stereo and group memberships locally, and the editor reports every
pruned element. Enhanced stereo groups belong to a connected molecule:
splitting a correlated group across components is rejected at publication,
while absolute memberships may split. Formats that cannot encode a relationship
fail explicitly.

`Perception` is reconstructible state derived from the exact represented graph
under an explicit model or policy: valence, rings, aromaticity, conjugation,
resonance delocalization, and installed CIP descriptors. Resonance
delocalization marks, per bond, whether moving a formal charge changes its
order; it reads the graph locally and enumerates no contributor. Resonance
contributors are separately requested results borrowing the source molecule;
default perception never enumerates them or rewrites represented bond orders
and charges. Task-specific descriptors, force-field typing, scoring, and
analyses belong in separate result objects; attaching selected values as
properties is deliberate.

`Perception` exposes only its sections, and each section owns its queries and
records mandatory provenance: the model, plus the policy where results depend
on it. Consumers reuse installed state only when its provenance matches what
they would perceive themselves, rather than revalidating it.

Chemical edits invalidate dependent perception; property changes do not.
Detached perception is checked against graph references and dimensions before
installation, and installation never repairs or rewrites represented chemistry.
Default perception does not add atoms, materialize stereo, or assign CIP.

Explicit hydrogens are separate graph atoms. Implicit hydrogens are all
attached hydrogens represented without graph atoms; total counts also include
explicit hydrogen neighbors. An atom's implicit count is either fixed or wholly
inferred, never a mix. Fixed counts remain known without perception, and
perception never overwrites them; inferred counts are unknown after chemical
edits until perception is recomputed. Explicit/implicit conversion preserves
composition and stereo. CIP assignment is transactional, respects represented
stereo, and reports unsupported or exhausted ranking rather than treating
unfinished work as a tie.

Perception runs once per reusable definition and publishes atomically: a
definition whose installed default perception is unchanged since perception
installed it is shared rather than perceived again, so perceiving after an edit
costs the definitions the edit changed. Perceiving a structure or structure set
installs a new topology snapshot that keeps the source's layout identity, so
selections, prepared potentials, and streaming bindings stay usable. Failure
leaves all definitions and geometry state intact.

See [graph storage](crates/kekule/src/molecule/graph.rs),
[represented stereo](crates/kekule/src/molecule/stereo.rs),
[perception state](crates/kekule/src/molecule/perception.rs),
[chemical perception](crates/kekule/src/algorithms/perception.rs), and
[stereo algorithms](crates/kekule/src/algorithms/stereo.rs).

## Topology, hierarchy, and classification

A published topology is nonempty and contains only complete connected molecule
instances and used definitions. Definition reuse is explicit; ordinary
construction does not intern chemically equal inputs, and structure-file
interpretation reuses definitions by the rule in
[Parsing, interpretation, and export](#parsing-interpretation-and-export).
Public traversal is instance-first. A definition is a storage and reuse
mechanism, not another chemical owner.

`Topology` is one flat owner. Definitions and instances are peers, and neither
nests inside a separate layout object. Dense orders, hierarchy, classification
overrides, and properties sit beside them. Layout identity is an explicit token
minted when a topology is published, not a consequence of how fields are
allocated or shared.

System atom and bond identities qualify a molecule-local ID with its instance.
Chain, residue, and atom-site IDs are topology-global. Dense indices are a
separate coordinate and property ordering, never interchangeable with semantic
IDs; publication establishes complete, deterministic, constant-time mappings
between them. Dense atom order is stored, not derived: format interpretation
keeps source atom-row order, builders default to instance order or accept an
explicit permutation, and editors and subsets keep surviving atoms in source
order. Dense bond order is instance order, then local bond ID.

`Hierarchy` is owned exactly once by `Topology`. Its storage is private: chain,
residue, and atom-site views are the only read path, and builders stage it
through one `HierarchyBuilder`. Atom sites refer to live `InstanceAtomId`
values and do not copy atoms. One chain or residue may span molecules, and one
molecule may span chains. Source labels and author identifiers remain typed
metadata with their namespaces and insertion codes. Hierarchy membership never
fabricates covalent bonds. Publication derives a private hierarchy index shared
by every snapshot of the layout, so hierarchy walks cost what they return.

Each reusable definition has one `MoleculeClass` shared by all its instances,
and each residue has one `ResidueClass`. These are intrinsic classifications,
not contextual roles such as receptor or ligand, and do not affect molecular
identity. Inference runs at topology publication without chemical perception
and combines conservative component recognition, graph evidence, and
inter-residue connectivity; conflicting strong evidence yields `Other`.
Explicit assignments keep their intent separately from inferred values, and
`ClassAssignment` reports which one a published class is. Unchanged complete
entities keep their classifications through subsets and appends; changed
chemistry or informative new context reclassifies the affected entities. See
[classification](crates/kekule/src/topology/classification.rs),
[builder publication](crates/kekule/src/topology/builder.rs), and
[component editing](crates/kekule/src/topology/editor.rs).

One-letter codes and sequences are derived, not hierarchy state. A residue's
code follows from its component identity, and only a monomer of a polymer has
one; a chain's sequence joins the codes of the monomers it holds. The sequence
a file states for a polymer entity, including unobserved monomers, is format
metadata owned by the interpretation report. See
[sequences](crates/kekule/src/sequence.rs),
[hierarchy](crates/kekule/src/topology/hierarchy.rs), and
[qualified lookups](crates/kekule/src/topology/lookup.rs).

## Construction, editing, and subsets

Published structural state changes through builders, editors, or explicit
transformations. Drafts may temporarily violate connectedness; publication may
not. Checked boundaries reject invalid references, dimensions, and empty
systems. `validate`, `try_build`, and `try_finish` support recovery without
consuming the draft. Mutable access must preserve safety immediately, without
depending on a mutation guard's destructor being called.

| Surface | Scope of work | Publication |
| --- | --- | --- |
| `MoleculeEditor` | One molecular graph, possibly disconnected while editing | Exactly one nonempty connected molecule |
| `TopologyBuilder` | Complete definitions, instances, hierarchy, classification, properties | Immutable topology |
| `StructureBuilder` | Complete composition with explicit coordinates | Valid structure |
| `TopologyEditor` | Instance-local atom, bond, component, and hierarchy changes | Repartitioned topology |
| `StructureEditor` | Structural edits coordinated with one structure | Topology and coordinates together |

`MoleculeEditor` is a draft, not a molecule view. Mutable chemistry access
invalidates perception immediately; validated no-op replacements retain state.
Drafts copy on write: a draft reads its source snapshot until an edit changes,
moves, or deletes an entity, publication shares every untouched definition and
its perception with the source, and the source snapshot is never written, so an
edit costs what it changes and undo histories keep every snapshot intact.
Stable draft handles reject deleted or foreign entities. Joining or splitting
components follows asserted connectivity. Surviving coordinates are preserved;
new atoms require supplied coordinates, unless the host asks a placement
operation on a `StructureEditor` draft (`structure::placement`) to derive them
from ideal local geometry. Placement is part of the host's edit and its one
publication, never moves an existing atom, and relaxes nothing. Geometry-only
edits retain the exact shared topology.

Append-oriented construction preserves existing IDs and dense order. General
structural edits keep surviving atoms in source dense order and append new
atoms in creation order, while repartitioned molecules receive deterministic
new definitions and instances. Hierarchy is filtered and empty nodes are
pruned; its partition is independent of molecule splits and merges. Per-entity
annotations follow the operation's explicit correspondence: an editor maps
surviving handles to published identities and, densely, source atom and bond
indices, residues, and chains to the published ones and back. Generic owner
annotations are conservatively cleared when their owner changes; no-op edits
preserve them.

Appending a complete structure imports its definitions, perception,
classification, hierarchy, coordinates, and entity properties into the draft
without merging equal hierarchy labels. A molecule identical to a definition
the draft already holds, with equal perception and class, becomes an instance
of that definition, as interpretation shares identical molecules; publication
never searches for identical definitions. Coordinates stay in the supplied
coordinate system. Incompatible cells, property types, or units reject the
whole append, and owner-property omissions are reported. See the
[append contract](crates/kekule/src/structure/structure_editor/append.rs).

Selections bind to one topology layout, including empty selections and set
operations. Atom and bond selections hold unique entities in dense order.
Atom-to-bond conversion explicitly chooses both, either, or exactly one selected
endpoint; bond-to-atom conversion selects both endpoints. Structural subsets
may cut molecules and repartition the induced graph into connected output
definitions; one operation-specific mapping transfers hierarchy,
classification, properties, and every structure array. There is no universal
topology remapping or provenance framework.

The selection language (`select`) is text that evaluates to an
`AtomSelection`; it is a query surface over these owners, not a second
selection type. A parsed `SelectionExpr` is independent of any topology and
prints canonically. Evaluation binds the result to the evaluated topology,
reads hierarchy, classification, properties, and installed perception without
changing them, and reads coordinates only when evaluated on a structure, so
`within` is one structure's result like other spatial selections. Hosts extend
the language through the `SelectionResolver` extension trait, whose answers
must bind to the evaluated layout or be points.

See [molecular editing](crates/kekule/src/molecule/molecule_edit.rs),
[system editing](crates/kekule/src/topology/editor.rs),
[structure editing](crates/kekule/src/structure/structure_editor.rs),
[selections](crates/kekule/src/topology/selection.rs),
[subsets](crates/kekule/src/topology/transform.rs), the
[selection grammar](docs/selection-language.md), and the
[evaluator](crates/kekule/src/selection_language/eval.rs).

## Properties and units

Properties are annotations owned at the narrowest scope whose lifetime matches
their validity. Each scope has one typed owner built from `OwnerProperties`
scalars and `PropertyTable<R>` columns, where the row type `R` names the entity
domain. Do not introduce parallel atom or bond data containers, or maps inside
every repeated entity.

| Scope | Owner | Row domains |
| --- | --- | --- |
| Molecule definition | `MoleculeProperties` | `AtomId`, `BondId` |
| Topology | `TopologyProperties` | Instances, dense atoms and bonds, chains, residues, atom sites |
| Geometry | `GeometryProperties` | Dense atoms and bonds |
| Structure set | `OwnerProperties` | Owner values only |
| State of a distribution | `OwnerProperties` | Owner values only |

Unsupported domains are unrepresentable rather than rejected at runtime. Reads
take the table's row type. Writes go through length-preserving `*Mut` guards:
whole-column replacement or removal, single-cell writes, and transactional
batches. A draft editor addresses allocated slots, and writing to a deleted
slot fails. `StructureEditor` does not dereference to its topology draft:
`topology()` reads it, `topology_mut()` reaches edits that need no geometry,
and its own `properties_mut()` holds structure-scoped values, so the two scopes
never shadow each other.

Keys are validated and never reserved. A column has one type and, for real
values, one physical unit; missing entries are explicit; column length matches
the owner's domain. Generic properties do not define chemical identity,
topology layout, or perception. Editing a published topology's properties keeps
its layout identity. Transformations transfer entity values explicitly and do
not infer annotation validity. Arbitrary source fields stay in format sidecars
unless canonical semantics, scope, and domain justify promotion. Occupancy and
B factors are typed site fields, independent of any generic property with the
same name. Positions, cells, weights, energies, time, steps, velocities, and
forces keep their dedicated APIs.

Kekule uses one runtime unit system. Every public value with a physical
dimension, accepted or returned, is a `Quantity` in a compatible unit; internal
state uses the library-wide canonical units. Mass and amount are distinct
dimensions. Unit conversion rejects unrepresentable dimensions or scales, and
checked numerical storage rejects nonfinite values. See
[units](crates/kekule/src/units.rs) and
[properties](crates/kekule/src/properties.rs).

## Geometries and structures

A `Geometry` holds all the coordinates and per-atom state of one arrangement,
and it is the one owned coordinate payload: a structure is a topology plus one
geometry, a structure set is a topology plus many, and a frame buffer holds
one for a stream. A geometry is shared and copied on write, so moving one
between owners (`Structure::into_parts`, `from_parts`, `to_structure`,
`from_structures`) copies no coordinates. A place in an order belongs to the
set and a weight to the ensemble, never to a geometry; an energy belongs to a
state, except as a geometry's own property. Coordinate kernels share code
through the borrowed `StructureView`.

A structure or structure set owns its shared topology once. A detached
geometry has atom rows only and acquires bond rows when a structure, set, or
frame buffer binds it. `StructureView`, obtained through the sealed
`AsStructureView` or a frame buffer, borrows a topology and a geometry, so
coordinate algorithms and writers share kernels across structures,
geometries, frames, and buffers, and read the geometry's velocities and forces
from the view. A view's `kind()` is `StructureKind::Standalone` or
`StructureKind::InSet` with a `SetContext` reporting the geometry's place and
its set's properties, so writers and reductions never lose that state
silently. Structure writers hold a geometry's place, velocities, and forces to
one policy, `geometry_metadata`, and per-atom, per-bond, own, and set
properties to a second, `properties`; each fails by default. Owned projections
(`to_structure`) materialize a structure with the whole geometry while sharing
the topology. A set built from structures (`from_structures`) consumes them
and requires equal layouts (`same_layout`), keeping the first structure's
topology.

Reads go through views: `Topology` lookups return `Option` views, a structure
reads its geometry through `Deref`, and `Structure::atom` joins an `AtomView`
with that atom's geometry state. Mutation goes through explicit,
dimension-preserving setters on detached geometries, the `GeometryMut` guard
an owner lends for a geometry it keeps bound (`Structure::geometry_mut`,
`StructureSet::get_mut`, `FrameBuffer::geometry_mut`), and `properties_mut()`
(`topology_properties_mut()` on structures and sets, which replaces a shared
topology with a snapshot of it). The owner checks lengths, property domains,
and table dimensions on insertion, replacement, or publication. Empty property
tables may acquire their owner's dimensions; populated tables may not be
silently resized or discarded.

Measurements and spatial selections operate on borrowed structures. Cartesian
measurements do not silently apply periodic imaging. A spatial selection is
one structure's result; reevaluating it across frames is explicit. Cartesian
geometry edits change coordinates atomically while keeping the exact topology,
cell, and properties. A rigid motion of a whole structure, geometry, or frame
rotates its cell with its coordinates; selection edits keep the cell. Prepared
distance, angle, and dihedral edits bind references and moving selections to a
layout, never to cached coordinates, and perform no ring closure, relaxation,
or periodic imaging. See
[positions](crates/kekule/src/structure/positions.rs),
[structures](crates/kekule/src/structure/owned.rs),
[measurements](crates/kekule/src/structure/measure.rs), and the
[geometry editing contract](crates/kekule/src/structure/geometry_edit.rs).

## Structure sets and ensembles

`StructureSet` is the one container for many geometries of one topology. What
makes a set a trajectory, a scan, or an ensemble is stated beside the
geometries, never carried by its type. Sources differ along two independent
axes and real ones combine them: an MCMC chain is ordered and a distribution, a
metadynamics run is time-ordered and weighted, a clustered trajectory is both.
A reader cannot know which it has — a multi-model PDB may be an NMR bundle or
an MD run, a multi-frame XYZ conformers or an optimization history — so a
reader loads geometries and the caller states what they mean, which is an
annotation rather than a conversion.

The two axes divide on what the source recorded against what an analysis
infers.

An order is recorded, so it belongs to the set. An MD engine writes each
frame's time and step, an optimizer numbers its steps, and a scan driver
records its coordinate. An order, when present, places every geometry: a time
axis for a path produced by dynamics, steps for an optimization, or values of
one coordinate in a unit for a scan or reaction path (values need not
increase). A time axis has times, steps, both, or neither, and each column
covers every geometry or none. Order values belong to the order, not to the
geometries: selecting keeps each selected geometry's value, and a removed
geometry takes its value with it.

A distribution is inferred, so it belongs to `Ensemble` and never to the set.
Reweighting a biased run, MBAR, clustering into states, and Boltzmann
populations from energies are interpretations, and one set often carries
several at once. An `Ensemble` owns one `StructureSet` and one `Distribution`,
validated together; `Ensemble::into_parts` returns both without copying. It
lends the set but never mutably, so everything that changes which geometries
exist goes through the ensemble, which remaps the samples. Several readings of
one set are several ensembles sharing every coordinate.

A distribution covers every geometry of its set or a sorted subset, its
samples. Weights are stored as natural logarithms, one finite value per sample:
weights from reweighting span hundreds of orders of magnitude, and a linear
weight would underflow. Normalized weights and probabilities are derived on
request; without weights the samples weigh equally. Log weights and the
partition are indexed by sample position. Every sample is in exactly one state
and every state has a representative sample; a state's population is the total
weight of its samples, derived and never stored. Only states have energies, on
one scale per distribution — a conformer's minimized energy or a computed free
energy, never a snapshot's instantaneous energy, which stays a property of its
geometry. Selecting and removing geometries remap the samples; a selection
keeps the energy and properties of the states it keeps whole.

A weight means something only against a target, so it is never a field of a
geometry or of a set. What the engine recorded per geometry — a bias
potential, collective-variable values, an instantaneous energy, a repeat count
— is a geometry property, and the weights an analysis uses are derived from
those. `Ensemble::uniform` and `Ensemble::weighted` build the one-state
readings, so weighting a trajectory needs no partition.

A geometry is shared through `Arc` and copied on write, so regrouping,
reweighting, selecting, and collapsing copy no coordinates, and editing one
geometry copies that geometry alone. Set and ensemble operations fail with one
`StructureSetError`, whose `Missing` variant names the `Requirement` an
operation lacks; a geometry's own edits fail with `GeometryError`.

Kekule records no runs or segments: replica exchange, metadynamics walkers,
umbrella windows, and repeated runs are the concern of the code that produced
them, and a pooled weighting across them is still one ensemble. `Order`
holds one coordinate, so a multi-dimensional scan uses `Order::Step` with the
coordinate values as geometry properties.

Where the geometries live is a type parameter. `StructureSet`
(`StructureSet<InMemory>`) holds them and lends them by reference, so its
accessors cannot fail. `StructureSet<IndexedFile>` holds none: it reads geometry
`i` from an indexed trajectory file as an owned `Geometry`, and a read can
fail. Both share the topology, properties, order, and the operations on them;
only a set in memory edits geometries or its order. A set of either storage
reads through the object-safe `StructureSetRead` trait, which is open to hosts.
Everything in Kekule that reads many geometries takes a `&dyn StructureSetRead`.
Kekule keeps no geometries a reader read and starts no threads: caching,
read-ahead, and parallel observation belong to the host.

A reader describes geometries and their order and never a distribution.
`StructureSet::from_reader` loads chosen geometries of any reader, in any
order, as a faithful copy with the reader's order values; a repeated index is
read once and shared. The host changes what geometries mean with `set_order`
and by building an `Ensemble`; `Ensemble::select` restricts a distribution to
the geometries it keeps. See
[structure sets](crates/kekule/src/structure/structure_set.rs),
[readers](crates/kekule/src/structure/structure_set/read.rs), and
[ensembles](crates/kekule/src/structure/ensemble/mod.rs).

## Trajectories and streaming

`kekule::traj` owns file codecs, streaming reader and writer contracts, the
reusable `FrameBuffer`, the `IndexedFile` storage of structure sets, periodic
preprocessing, streaming superposition, and lazy processing (`Processed`). It
owns no second in-memory trajectory type. A trajectory writer takes a
`FrameView`: a borrowed `Geometry` with its time and step.
`read_structure_set` reads a file as a structure set of a time order, and
`write_structure_set` writes a set of a time order or of none.

A trajectory represents one fixed-topology epoch. Changing chemistry or
hierarchy requires a new topology and explicitly constructed geometry for the
next epoch. Frame selection is separate from atom subsetting: requested order,
duplicates, and empty selections are supported without renumbering stored time
or step.

File readers receive the topology and interpret file atom index as its dense
atom index. They validate counts and available metadata; equal counts alone do
not prove atom identity. Structure writers emit dense order, so a coordinate
file written in a structure file's row order addresses the same atoms through
that file's topology.

An indexed reader never returns a frame it has not verified. A file that only
grew by appending is the same file with more frames, and refreshing verifies
and indexes only the appended frames; a file whose existing bytes changed is
not the same file. A final frame its writer has not finished is pending, not
corrupt. `append_trajectory` continues a checkpointed run's file so that it
equals a file written straight through.

Streaming buffers bind to one topology layout and publish complete frames
transactionally; a published frame replaces every field. Path writers publish
only a completed nonempty trajectory, and codec writers append whole frames in
order and never rewrite one, so a flushed stream can be read while it grows.
Frame state a codec cannot carry fails unless the caller explicitly allows
omitting it, and the writer reports what it left out. A codec never invents
frame state, except times and steps numbered under an explicit timing policy,
which it reports.

`StructureSet<IndexedFile>` decodes geometry `i` from frame `i` on every read,
from any number of threads at once, and reads whatever bytes the file holds at
the indexed offsets; `refresh` takes appended frames and fails without change
when the indexed bytes changed. A clone shares the set's reader and index and
refreshes independently. An `Ensemble` over the set refreshes through
`Ensemble::refresh`, and appended frames are not samples. `load` reads listed
geometries into a set in memory; an ensemble's `load` also restricts the
distribution. See
[structure sets read from files](crates/kekule/src/traj/io/indexed_file.rs).

`Processed` reads another `StructureSetRead` and runs `ProcessingSteps` on each
geometry as it is read, in a fixed order. Processing neither removes nor
reorders geometries, so indices are preserved, and it carries no distribution.
Its steps call the same kernels as the in-memory functions and give the same
bits. See [lazy processing](crates/kekule/src/traj/processed.rs).

Loaded and streaming transformations share kernels. Collection operations,
including superposition and periodic preprocessing, change the collection in
place and stage all affected state before publication; clone first to keep the
source. Superposition rotates cells, velocities, and forces consistently.
Periodic reconstruction uses asserted bonds and checks ring closure. Imaging
acts on whole molecules, per frame and in any frame order. Unwrapping keeps
temporal state, requires ordered and sufficiently close geometries, and
precedes downsampling. Failed streaming operations change neither the frame
buffer nor temporal state. See
[frame buffers](crates/kekule/src/traj/trajectory/buffer.rs),
[file I/O](crates/kekule/src/traj/io/mod.rs),
[periodic operations](crates/kekule/src/traj/periodic.rs), and
[streaming contracts](crates/kekule/src/traj/trajectory/mod.rs).

## Analysis

Analyses are separate result objects, never new topology state or a generic
reduction framework. Preprocessing such as alignment or imaging is explicit.
Operations that may run on a `parallel::Executor` return the single-threaded
result or error for every executor and thread count.

RMSF and contact occupancy are weighted statistics in `kekule::statistics`.
Every observation is a borrowed structure with an explicit, finite, positive
weight, so no structure is reduced with an implied weight. `Ensemble` methods
observe each sample with its probability and leave other geometries out;
`StructureSet` methods observe every geometry with weight one. Streamed frames
reach the same accumulators through their buffer's `StructureView`. Accumulator
memory is bounded by the selected atoms, pairs, or groups, not by observation
count. Failed observations leave accumulators unchanged. Accumulators of the
same selection merge, so a host observes ranges on any threads and merges them
in range order with the same result on every executor. Results keep their atom
and pair associations; observation indices are diagnostic labels. See
[weighted statistics](crates/kekule/src/statistics.rs).

Clustering groups the samples of an ensemble's distribution, weighted, or every
geometry of a structure set, unweighted. It returns a `Partition`, the
grouping type a distribution owns, and never changes the set.
`DistributionMut::regroup` installs a partition without changing any weight or
sample order: a new state containing every sample of an old state keeps that
state's energy and properties, and a state split apart loses them. Collapsing
makes a new collection of each state's representative carrying the state's
population. Pairwise distances are a separate `DistanceMatrix` result. See
[clustering](crates/kekule/src/clustering.rs).

Superposition and RMSD are methods of structure sets, against one of the set's
geometries or an independent view, and act on every geometry. Geometric
correspondence explicitly pairs selected atoms from two exact snapshots, one to
one, without claiming chemical equivalence; pairing by hierarchy labels reads
one namespace and reports unpaired atoms. Symmetry-corrected alignment keeps
the smallest RMSD over the permutations of equivalent atoms, found by an exact,
bounded automorphism search that never exchanges molecule instances; exceeding
a bound is an error, never a truncated result. See
[alignment](crates/kekule/src/alignment.rs),
[correspondence](crates/kekule/src/alignment/correspondence.rs), and
[symmetry permutations](crates/kekule/src/alignment/symmetry.rs).

DSSP assignment reads one structure and returns a snapshot result. A `DsspPlan`
prepares what it reads from the topology once and binds that layout, so the
frames of a trajectory reuse it. See [DSSP](crates/kekule/src/dssp/mod.rs).

Non-covalent interactions split by what they depend on. `InteractionFeatures`
classifies one topology's interaction features once per reusable definition
and is reused for every structure sharing that layout; it never changes the
topology. A search reads one structure or two and returns a snapshot of
interactions addressed by dense indices of each side's topology, with criteria
as a plain options value in `Quantity` units. Hydrogen-bond occupancy is a
weighted, mergeable accumulator with the rules of the weighted statistics. See
[interactions](crates/kekule/src/interactions/mod.rs).

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
SDF and SMILES readers split records by line framing alone, so a malformed
record fails alone and these formats stream with memory bounded by one record.
CIF block boundaries are lexical, so malformed framing fails the rest of the
input. A file-level document is a convenience collection of parsed records;
interpretation is always per record and never depends on a sibling record.

A format document is the typed syntax of one record or block, with source
locations: what the source states, before chemistry. Reading source text,
including CXSMILES extensions and Molfile property records, belongs to parsing;
interpretation consumes typed syntax and never re-reads text. CIF documents
stay generic data trees.

Interpretation takes a document, explicit options, and optional context —
external resources a format needs but does not carry, such as component
templates. It constructs represented chemistry, hierarchy, classification, and
available geometry once into the richest canonical owner plus a report; it does
not run general perception, standardization, tautomerization, or protonation.
Format convenience functions compose this pipeline without reimplementing it.

mmCIF, PDB, and XYZ interpretation gives identical molecules one definition, so
perception runs once per distinct molecule. Two molecules are identical when
their local atoms, bonds, represented stereo, properties, and the atom labels
the format gives (component and atom name for mmCIF and PDB) are equal.
Hierarchy and classes stay per instance, and editors copy a shared definition
when an edit touches one instance. See
[definition reuse](crates/kekule/src/io/definition_reuse.rs) and
[sharing on append](crates/kekule/src/topology/editor/sharing.rs).

Interpretation never silently discards source content. Unsupported content is
rejected by default; an explicit option downgrades it to typed report entries
naming the record and the kind of content omitted. A deliberately lossy
projection, such as base-only SMILES, says so in its name. Every stage reports
through one diagnostics model: parse, interpretation, and write failures of
every format are one `diagnostics::FormatError` with a structured kind and a
location relative to the original input, and non-fatal report issues identify
the records they concern.

Coordinate-only trajectory formats are the deliberate exception: they carry no
chemistry, so codecs decode frames against a caller-supplied topology. XYZ is
read both ways: `kekule::xyz` interprets its frames into structures through
this pipeline, and `kekule::traj` streams its coordinates against a topology.

| Format scope | Interpretation and projection |
| --- | --- |
| SMILES record | All connected components, with topology projection |
| Molfile | One structure with topology and molecule projections |
| SDF record | Structure, title, fields, and reports; each record independent |
| CIF data block | One structure, or compatible coordinate models as a structure set |
| PDB entry | One structure, or compatible coordinate models as a structure set |
| XYZ file | One frame as a structure, or frames of the same atoms and inferred bonds as a structure set |

Sibling SDF records and mmCIF blocks never merge implicitly. Singular
projections check cardinality; plural molecule projections keep every connected
component in source order. In a conformation the coordinates are the
configuration: annotations that contradict them yield to them with a warning.
mmCIF keeps label and author hierarchy identities; PDB carries author
identities only. Multiple coordinate models do not imply a time order, and the
formats state no distribution: every coordinate model becomes a geometry, a
host that reads the models as draws builds an `Ensemble`, and occupancy is
never converted into a weight. Interpretation selects whole alternate-location
alternatives and records its choices in the report; occupancy ranks
representatives, not structure probabilities.

Macromolecular formats do not carry complete bonds, and molecule boundaries
follow bonds, so one resolver shared by mmCIF and PDB derives connectivity
from ranked evidence: explicit source bonds, component templates supplied as
context, polymer linkage within bonding distance, and distance inference for
residues nothing else covers, which is on by default and reported. Explicit
bonds outrank the levels below them, and close contacts no evidence bonds are
counted, never bonded. Metal coordination between residues is not a bond: the
report keeps it as connection metadata, while coordination within one residue
is a dative bond, so a metal complex that is one component stays one molecule.
Hydrogen bonds, salt bridges, and links to symmetry mates are never bonds.
Kekule bundles Chemical Component Dictionary templates for standard residues,
water, common ions, cofactors, and crystallization additives, plus the neutral
histidine tautomers; callers supply more as context. See
[macromolecular connectivity](crates/kekule/src/io/macromolecular/connectivity.rs)
and [structure interpretation](crates/kekule/src/io/macromolecular/mod.rs).

Interpretations own format reports, metadata, and source correspondence.
Borrowed projections reuse canonical owners; consuming projections explicitly
discard richer information. Document-level reports borrow record reports
instead of duplicating them. Format-specific constructors and save methods do
not belong on canonical owners.

Writers live in format namespaces as one `write(source, options)` returning a
string and one `write_to(writer, source, options)` streaming into the caller's
writer without buffering it. Each source type is a struct with private fields
that converts from the inputs the format can represent. Coordinate-bearing
writers emit atoms in dense order. Unsupported representational content fails
explicitly: set or frame state the format cannot carry fails unless the caller
sets `UnrepresentablePolicy::Omit` for that kind of content, each kind with its
own policy, and the writer's `WriteReport` names every omission. Writers state
what their reader cannot derive and no more, never invent assertions, and
never write a distribution as a field of a record. There is no universal save
trait or implicit trajectory-to-structure export. Writer failures say why no
output was produced: the format cannot encode the content, the source lacks
required perception, a work limit was reached, the source is inconsistent, or
the sink failed. Per-format rules live with the
[format namespaces](crates/kekule/src/lib.rs) and their writers.

## Identity and reconstruction

Molecular queries own predicate graphs, recursive subqueries, output tags, and
local stereo carrier frames; they do not reuse represented molecular graphs.
Matching reads explicit target perception without installing or changing it.
Prepared molecular targets borrow their source; topology matches keep the exact
shared snapshot. See the [SMARTS contract](crates/kekule/src/query/dialect.md).

Molecular equality compares authoritative represented chemistry, excluding
perception, annotations, hierarchy, and classification. Topology layout
equality is separate: it includes definitions, instances, classifications,
hierarchy, semantic IDs, and dense order, but excludes perception and
properties. Layout identity is stricter than equal layout: every publication
mints a new layout identity, and perception and property edits create
snapshots that keep it. Operations that address atoms and bonds by index accept
any snapshot sharing their layout and never substitute an independently equal
topology; operations that read perception bind the exact snapshot.

Reconstruction validates graph connectedness first, then property references
and dimensions, then installs checked perception. Disconnected persisted graphs
must be partitioned or rejected. Export to scientific formats is not exact
persistence and must not weaken these boundaries.

Exact persistence is the `native` format module, the one owner of Kekule's
binary encoding of molecules, topologies, structures, structure sets, and
ensembles. It stores each owner's authoritative state and reconstructs it only
through the owners' validating paths; stored perception is reinstalled or
discarded. A decoded topology keeps every semantic ID but is a new publication
with a fresh layout identity. Decoding is bounded before allocation and never
trusts the input. Semantic IDs have a stable encoding so applications can
persist references next to native files. See the
[format specification](docs/native-format.md) and
[topology reconstruction](crates/kekule/src/topology/reconstruct.rs).

## API and maintenance rules

The public API of every crate follows these rules. Code that breaks one is a
defect to fix, not a precedent.

### Publication and module layout

Implementation modules are private. Each public module is an explicit list of
re-exports: no glob re-exports and no public implementation files. Every public
item has exactly one public path; `prelude` is the only module that repeats
items. An item that no public path reaches is `pub(crate)`. Workspace lints
(`unreachable_pub`, `missing_docs`, `missing_debug_implementations`) enforce
this.

A public module's file holds only its docs, `mod` declarations, and re-export
lists. Its implementation, entry-point functions included, lives in private
child or sibling modules named for what they hold. Format implementations live
in the private `io` module, one directory per format holding its parser,
document types, interpretation, writer, and entry points; plumbing shared by
several formats stays at the `io` level.

The public modules of `kekule`, by domain:

| Domain | Modules |
| --- | --- |
| Foundations | `units`, `math`, `properties`, `diagnostics`, `parallel` |
| Chemistry | `molecule`, `stereo` (`cip`, `coordinates`, `cleanup`), `perception` (`valence`, `rings`, `aromaticity`, `conjugation`, `resonance`), `hydrogens`, `descriptors`, `sequence` |
| Queries | `query`, `smarts`, `substructure`, `select` |
| Systems and structures | `topology` (`edit`), `structure` (`measure`, `placement`), `solvation` |
| Analysis | `alignment`, `clustering`, `statistics`, `dssp`, `interactions` |
| Formats | `smiles`, `molfile`, `sdf`, `mmcif`, `pdb`, `xyz`, `macromolecular`, `components`, `native` |
| Trajectories | `traj` (`io` with `xyz`, `dcd`, `trr`, `xtc`; `periodic`; `analysis`) |

A type lives in the module of the domain that owns it: perception state with
its perception step, represented stereo in `stereo`, derived per-molecule
results in `descriptors`, interpretation types shared by mmCIF and PDB in
`macromolecular`, and the shared failure vocabulary in `diagnostics`. Format
modules stay at the crate root, except the trajectory codecs, which live in
`traj::io` beside the streaming contracts they implement. The foundational
crate has no optional features.

### Names, conversions, and identifiers

- **Names.** Every public type name is unique across the crate, so no import
  needs renaming. A type takes its owner's prefix (`GeometryProperties`,
  `StructureSetError`, `FrameView`) when its plain name would be ambiguous.
- **Conversions.** `as_*` for cheap borrowed views, `to_*` for owned conversion
  keeping the source, `into_*` for consuming conversion, which should transfer
  compatible storage.
- **Mutation verbs.** `perceive` changes its receiver in place; `perceived`
  returns a new value and leaves the receiver unchanged.
- **Ownership.** Constructors that take ownership of chemistry accept values
  (`Molecule`, `Structure`), and callers clone explicitly.
- **Getters.** No `get_` prefix. Counts are `*_count()`; `len()` belongs only
  to collection types.
- **Properties.** Use full `properties` names.
- **Units.** No public name carries a unit (`_angstrom`, `_degrees`).
- **Identifiers.** Semantic IDs and dense indices are distinct types; dense atom
  and bond indices are `TopologyAtomIndex` and `TopologyBondIndex`. A method
  taking a dense index is named `*_at` when a twin takes the semantic ID
  (`atom` and `atom_at`). Frame indices are `usize`; `u64` is used only for
  byte offsets. Identifiers display in one style.

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
- **Mutation** stays on editors, `*_mut` guards, and semantic owner APIs. A
  guard over a value that must stay bound to its topology, such as a set's
  geometry, has no `DerefMut`; whole values are replaced through the owner.
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
API inventories, algorithms, and per-format policies live in Rustdoc beside the
implementation, stated once; this document links to them rather than restating
them.

### Companion packages

The unpublished `kekule-bench` workspace package (`benchmarks/observer`)
reports what Kekule's public APIs compute, keyed by source positions, for
comparison with pinned reference toolkits. Datasets, reference adapters,
comparison, and the list of known differences stay outside the runtime crates;
reference toolkits are never runtime dependencies. See the
[benchmark README](benchmarks/README.md).

`kekule-openff` owns compiled SMIRNOFF rules, configured NAGL inference, and
typed `ParameterizedTopology` results keeping the caller's `Arc<Topology>`. It
prepares temporary per-definition chemistry and expands assignments to
instance-qualified atoms without mutating the topology or adding hydrogens. All
numeric parameters carry canonical units, and complete force-field state stays
in this result rather than in molecule properties. The crate is pure Rust.
Model bundles own feature ordering, network configuration, and lookup data;
OFFXML keeps the required checkpoint identity, which parameterization checks
against the supplied model. The converted Ash model the bundled Rosemary preset
requires ships as data in the separate `kekule-openff-ash` crate, a default
feature, so the preset works without external files while the main crate stays
within the crates.io size limit. Every failure is one typed error with a stable
kind, the failing molecule definition, and every atom it concerns. Checkpoint
conversion and the provenance of shipped data live in `tools/openff`.

`kekule-potentials` owns energy evaluation and geometry optimization; the
foundational crate owns no potential. A prepared `Potential` binds one topology
layout, does not change during evaluation, and accepts any `StructureView`
sharing that layout; an independently equal topology is incompatible. Energies
and gradients are separate result objects in canonical units, not geometry
properties. Backends lower explicit parameters, currently an OpenFF
`ParameterizedTopology`, into private kernels. Unsupported geometry state, such
as a periodic cell, is rejected rather than ignored. Gradients are exact and
never capped. Minimization reads a borrowed view and returns new positions
without changing chemistry, topology, cell, or properties; fixed atoms keep
their input coordinates bit for bit, and position restraints add to the
minimized energy without becoming part of a potential. See the
[potential contract](crates/kekule-potentials/src/lib.rs) and
[OpenFF backend](crates/kekule-potentials/src/openff.rs).
