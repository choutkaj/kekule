# Architecture

## Purpose

This document is the normative architecture contract for `kekule`. It defines
the ownership, semantic boundaries, and invariants of the core molecular data
model. Detailed API behavior belongs in Rustdoc and tests.

`kekule` is a pure-Rust foundation for cheminformatics, structural
bioinformatics, molecular structure handling, and molecular modelling.

The foundational chemical object is `Molecule`: one non-empty connected,
geometry-independent molecular entity. The foundational system object is
`Topology`: one coordinate-free molecular system composed from one or more
molecule instances and their system-level hierarchy. Geometry belongs above
topology.

## Canonical object model

```text
source text / bytes
    -> format-specific parsing
    -> FormatDocument
         -> optional format-specific Record values for record-oriented formats
    -> interpretation / canonical publication
         -> connected canonical Molecule components
              -> Vec<Molecule> when geometry/system context is not requested
              -> Topology + Positions -> Model when one record carries geometry
              -> higher geometry objects such as Ensemble when the format
                 semantically contains several realizations of one topology
```

The intended ownership hierarchy is therefore:

```text
Molecule
  one connected geometry-independent molecular entity
  Graph + Perception + definition-scoped Properties

Topology
  one geometry-independent system made from one or more Molecule instances
  reusable MoleculeDefinition values with canonical MoleculeClass
  topology-wide atom/bond identity and dense ordering
  one system-level Hierarchy, which may be empty
  canonical ResidueClass for hierarchy residues
  system-scoped Properties

Model
  one geometry-dependent realization of a Topology
  realization-scoped Properties

Ensemble
  several non-temporal realizations of one Topology

Trajectory
  several temporally ordered realizations of one Topology
```

A salt, noncovalent complex, solvent box, protein-ligand system, or DNA duplex
is not represented by weakening `Molecule` into a disconnected graph. It is
represented by several connected molecules at the `Topology` level.

An asserted topological `Bond` contributes to molecular connectedness. Spatial
association, hydrogen bonding, ionic attraction, contact, or any other
non-topological interaction does not.

Examples:

```text
ethanol                         -> 1 Molecule
ubiquitin                       -> 1 Molecule
Na+ + acetate-                  -> 2 Molecules
protein + noncovalent ligand    -> 2 Molecules
covalent protein-ligand adduct  -> 1 Molecule
DNA duplex                      -> typically 2 Molecules
disulfide-linked protein chains -> 1 Molecule
```

## `Molecule`

There is exactly one foundational molecular type. Kekule does not distinguish
`SmallMolecule` and `MacroMolecule` as owning wrappers.

The intended core shape is:

```rust
pub struct Molecule {
    graph: Graph,
    perception: Perception,
    properties: Properties,
}
```

These fields have deliberately different semantic roles:

```text
Graph
  authoritative represented chemistry
  required
  defines the connected molecular entity

Perception
  derived chemical interpretation
  reconstructible from represented chemistry plus an explicit perception model
  does not define molecular identity

Properties
  extensible geometry-independent annotations valid for this molecular definition
  may target the molecule itself, its atoms, or its bonds
  do not define represented chemistry or molecular identity
```

Residue, chain, polymer, asymmetry, and atom-site organization do not belong to
`Molecule`. They belong to the system-level `Hierarchy` owned by `Topology`.

This boundary is intentional. Covalent connectedness and biological/source
hierarchy are independent partitions of atoms: one hierarchy chain may span
several disconnected `Molecule` instances, and one connected `Molecule` may span
several hierarchy chains.

`Molecule` contains no coordinates, conformers, velocities, periodic cell,
residue/chain hierarchy, or other system/geometry-dependent state. Its
properties must likewise be geometry-independent and valid for every use of that
molecular definition. Properties that differ between instances of the same
molecule definition belong at `Topology` scope instead.

### Connectedness invariant

Every published `Molecule` is non-empty and connected through asserted
topological bonds. A single atom is a valid connected molecule.

A disconnected graph is never a valid `Molecule`.

Temporary disconnectedness is permitted only inside construction/edit staging
such as `MoleculeEditor`.

This is a type-level architectural invariant, not merely a convention. Public
construction and editing APIs must make it impossible to publish an invalid
disconnected molecule.

## `Graph`

`Graph` is the authoritative chemical graph of a molecule.

Conceptually it owns:

```text
Graph
  atoms
  bonds
  adjacency/connectivity
  stable local AtomId / BondId identity
  represented stereochemical elements/groups
```

The exact physical storage is an implementation detail, but it should remain
purpose-built for chemistry, compact, deterministic, and efficient.

Generic annotations do not live inside `Graph`, `Atom`, or `Bond`. Definition-
scoped object, atom, and bond annotations live in the containing `Molecule`'s
`Properties`. This keeps represented graph chemistry separate from extensible
metadata and prevents generic annotations from affecting molecular identity or
perception invalidation.

### Represented atom chemistry

Fundamental asserted atom state belongs in the graph, for example:

```text
element
isotope
formal charge
radical state
represented hydrogen declaration
atom-map identity when retained as canonical represented chemistry
```

### Represented bond chemistry

Fundamental asserted bond state belongs in the graph:

```text
endpoints
localized represented bond order/kind
represented bond stereochemistry where applicable
```

A `Bond` is a topological relation. A noncovalent interaction must not be added
as a `Bond` merely because it relates two atoms spatially or energetically.

### Stereochemistry

Canonical represented stereochemical state belongs to the graph because it is
part of the asserted molecular representation.

Source-format marks such as SMILES directional syntax or molfile wedges are
format/interpreter state. They must be resolved into Kekule's canonical stereo
representation before a molecule is published.

CIP labels are derived and therefore belong to `Perception`, not `Graph`.

### Aromaticity

Aromaticity is perceived chemistry, not a canonical bond order.

Canonical graph bonding is localized. Aromatic source syntax may be accepted by
format readers/interpreters, but a published `Graph` contains ordinary
represented bond orders. Aromatic atom/bond membership belongs to `Perception`.

## `Hierarchy`

`Hierarchy` is authoritative coordinate-independent organization of atoms at
`Topology` scope.

It may be empty. Its presence does not create a different molecular or topology
type.

Typical hierarchy state includes:

```text
chains
residues
residue classification
polymer organization
residue/chain identifiers
component names
atom-site annotations
mappings from hierarchy atom sites to InstanceAtomId
```

Fixed semantic hierarchy fields remain strongly typed hierarchy state. Generic
annotations targeting chains, residues, atom sites, or the system as a whole are
stored in `Topology`'s `Properties`, not as independent property maps embedded in
each hierarchy node.

Hierarchy is orthogonal to molecular connectedness:

```text
Molecule / Graph answers:
  which atoms form one connected molecular entity, and how are they bonded?

Topology / Hierarchy answers:
  how are system atoms organized into residues, chains, polymers, asymmetry
  groups, and related source-level structural organization?
```

Hierarchy does not own independent atoms or bonds. Every hierarchy atom site
must resolve to one live `InstanceAtomId` in the containing `Topology`.

Hierarchy node identities are topology-global. A chain, residue, or atom-site ID
belongs to one `Topology`; it is not a molecule-local ID requiring an additional
`MoleculeInstanceId` qualifier.

Conceptually:

```text
Hierarchy
  ChainId -> Chain
    ResidueId -> Residue + ResidueClass
      AtomSiteId -> AtomSite -> InstanceAtomId
```

These are the canonical public hierarchy names. `Chain`, `Residue`, and
`AtomSite` are topology-owned storage nodes; topology-bound `ChainView`,
`ResidueView`, and `AtomSiteView` values provide borrowed navigation context.

Kekule intentionally does not reproduce a `Structure -> Model -> Chain ->
Residue -> Atom` object hierarchy. In particular, hierarchy has no `Model`
node, and an `AtomSite` remains metadata and organization referring to the
authoritative chemical atom through `InstanceAtomId`.

### Hierarchy may cross molecule boundaries

Hierarchy and connected molecular identity must not be forced to have the same
boundaries.

Examples:

```text
one source chain with an unresolved break
  -> two connected Molecule instances
  -> one Chain containing residues/atom sites from both instances

covalently disulfide-linked source chains
  -> one connected Molecule instance
  -> two Chains inside one Topology hierarchy

many disconnected waters sharing one source asymmetry identifier
  -> many Molecule instances
  -> hierarchy organization may group their residues under one source chain/asym
```

Kekule must never fabricate bonds merely to preserve hierarchy grouping, and it
must never duplicate hierarchy nodes merely because a chain crosses a molecule
boundary.

### Molecule-centric hierarchy views

`Molecule` itself does not own hierarchy. At system scope, an instance-first
view may nevertheless expose convenient filtered hierarchy navigation:

```text
topology.molecule(instance).chains()
topology.molecule(instance).residues()
topology.molecule(instance).atom_sites()
```

Such APIs are borrowed/filtering views over the one topology-owned hierarchy.
They do not create or own per-molecule hierarchy copies.

### Domain-specific APIs

Protein-, nucleic-acid-, or polymer-specific algorithms do not require owning
wrapper types.

They may operate on `Topology`/`Hierarchy`, `Model`, or lightweight borrowed
validated views such as `ProteinView<'_>` or `NucleicAcidView<'_>`.

Such views interpret existing topology/hierarchy state; they do not own another
molecular object.

## `Perception`

`Perception` is the installed derived interpretation of one exact represented
molecular graph.

It replaces the architectural role formerly described by `PerceptionState`.
Within `Molecule`, `perception: Perception` is idiomatic Rust and unambiguous.

Perception is semantically subordinate to represented chemistry:

```text
Graph
  + perception model/policy
        |
        v
Perception
```

Deleting all perception state must never destroy authoritative information
about the molecule.

Two molecules must not become chemically different merely because one has more
derived perception cached than the other.

Fundamental perception may include sections such as:

```text
valence / implicit hydrogens
ring/cycle information
aromaticity
CIP assignments
```

Not every calculated property belongs in `Perception`. Fingerprints,
descriptors, partial charges, force-field types, pharmacophore features,
rotatable-bond classifications, scoring terms, and other task-specific results
remain separate derived objects by default. When a caller deliberately attaches
such a result to a domain object and its validity matches that owner's scope, it
may be stored through the generic `Properties` layer instead. This does not turn
it into fundamental chemical perception.

### Structural graph derivations versus chemical perception

Kekule may internally distinguish mathematically unique graph caches from
model-dependent chemical perception.

For example:

```text
graph-derived:
  degree
  connected traversal data
  generic cycle membership

chemically perceived:
  valence model
  implicit hydrogens
  aromaticity model
  CIP assignment
```

This distinction may be reflected internally if useful, but the public
`Molecule` architecture remains `Graph + Perception + Properties`.

### Perception installation and invalidation

Perception must always correspond to the current authoritative graph.

Graph-changing edits invalidate affected perception. The implementation should
prefer simple, safe invalidation over a complex dependency engine.

Exact reconstruction of externally stored perception may be supported through
checked installation APIs, but installation must validate references and
dimensions and must never rewrite authoritative graph chemistry.

Public convenience APIs may expose perception-backed queries directly through
`Molecule`; callers should not need to duplicate perceived flags into generic
properties.

Mutating generic properties does not change represented graph chemistry and must
not invalidate perception.

### Perception on canonical owning objects

Explicit perception is available after canonical publication at every owning
level. `Molecule::perceive()` installs the default valence, ring-set, and
aromaticity profile transactionally. `Topology::perceived()` returns a new
topology snapshot with that same profile computed once per reusable molecule
definition, in definition order. It recomputes installed perception and clears
dependent CIP state according to the molecular pipeline's invalidation rules.
It does not assign CIP or add explicit hydrogen atoms.

`Model::perceive()`, `Ensemble::perceive()`, and `Trajectory::perceive()` install
that new topology snapshot only after every definition succeeds. Failure
identifies the source molecule definition and leaves the entire receiving owner,
including its topology allocation and previously installed perception, unchanged.
Collection perception is independent of member/frame count; it never runs once
per realization or once per instance of a reused definition.

Perception preserves represented graphs, definition reuse, instances, semantic
IDs, authoritative dense ordering, hierarchy, classifications, and all stored
properties. Realization arrays, periodic cells, weights, time, step, and other
realization/collection state remain intact without copying realization payloads.
The operation changes no represented chemistry and needs no geometry remapping.

Published shared topologies remain immutable. Successful model/collection
perception always installs a new `Arc<Topology>` snapshot, even when the old
allocation has no other owners. Other owners retain their old snapshot.
Selections, prepared calculations, and trajectory readers/buffers retain their
original topology bindings; layout equality does not implicitly transfer them.

## Editing

All structural mutation of `Molecule` happens through `MoleculeEditor` or an
equivalent transactional staging type.

Conceptually:

```text
Molecule
   |
   | edit
   v
MoleculeEditor
   |
   | arbitrary intermediate structural edits
   | temporary disconnection is allowed
   | temporary incomplete chemistry is allowed
   v
finish()
   |
   +-- invalid -> error
   |
   `-- valid -> Molecule
```

`MoleculeEditor` is allowed to violate publication invariants while editing.
The finished `Molecule` is not.

Editing operations should be simple graph operations; validity is enforced at
publication rather than through specialized mutation APIs.

`finish()` must at minimum validate:

```text
non-empty graph
exactly one connected component
valid atom/bond references
valid adjacency
valid stereo references
valid property-table dimensions for retained atom/bond identity spaces
```

Hierarchy validation is not a `MoleculeEditor` responsibility. Hierarchy is
validated when a `Topology` is published.

Chemical perception need not be valid during editing. On successful
publication, stale perception must be discarded, recomputed, or explicitly
reinstalled through a checked path. Per-atom and per-bond properties may be
projected through editor identity when their target entity survives; owner-level
properties must follow the explicit transformation semantics rather than being
blindly assumed valid for a structurally changed molecule.

The same editor concept may be used for construction from scratch; a separate
public `MoleculeBuilder` is not architecturally required unless it provides
clear ergonomic value without duplicating semantics.

`MoleculeEditor` is the single public construction and editing interface. It
exposes live atom/bond and stereo inspection, connectivity queries, represented
state replacement, checked bond rewiring, batch deletion/retention, and owner and
entity property editing. Complete property columns follow live entity order;
deleted storage slots are handled internally. Fragment append returns atom,
bond, stereo-element, and stereo-group ID correspondence and preserves entity
annotations and represented stereo. It does not import fragment owner properties
or perception; conflicting property data fails transactionally.

`Molecule::edit()` clones into detached state; `Molecule::to_editor()` moves it.
`finish()` consumes the draft. `validate()` checks a snapshot, and `try_finish()`
keeps a rollback snapshot so failed publication can return an unchanged editor
for repair. The latter operations explicitly trade a clone for recoverability.
Published-molecule algorithms require a finished `Molecule`: editors do not
dereference to `Molecule`, including in tests. Internal algorithm fixtures that
intentionally use unfinished state access it explicitly through crate-private
test paths. Bond mutation cannot replace endpoints outside checked rewiring.

Mutable atom access invalidates perception and owner annotations immediately;
invalidation never depends on a guard destructor. Explicit identical atom
replacement and unchanged bond-order setters preserve draft state. Checked
perception installation belongs on the published `Molecule`, after `finish()`.

Public unrestricted mutable access to graph internals should not bypass the
editor and thereby bypass publication validation.

## Parsing and interpretation

### Canonical parsing pipeline

Kekule uses one format boundary and one canonical publication path:

```text
source text / bytes
    -> parse
format-specific Document
    -> select the format's independently interpretable scope when nested
       (Record, Block, or equivalent)
    -> interpret / canonicalize
format-specific Interpretation
    -> borrowed views or owned projections
       -> Vec<Molecule>
       -> Topology
       -> Model
       -> Ensemble / Trajectory when the source semantics justify them
```

Parsing is format-specific. Canonical target semantics are not. The interpretation
object should retain the richest canonical Kekule state justified by the selected
source scope; simpler outputs are projections that deliberately discard
information rather than independent reinterpretations of the source.

The canonical target hierarchy is:

| Source scope | `Vec<Molecule>` | `Topology` | `Model` | `Ensemble` / `Trajectory` |
| --- | --- | --- | --- | --- |
| SMILES record/document | natural | natural | not represented | not represented |
| Molfile document | lossy projection | lossy projection | natural | not represented |
| SDF record | lossy projection | lossy projection | natural | not represented |
| mmCIF block, one selected coordinate model | lossy projection | lossy projection | natural | not represented |
| mmCIF block, several compatible coordinate models | projection of one selected interpretation policy | shared topology | selected-model projection | `Ensemble` is natural |
| trajectory/coordinate stream with external topology | usually not reconstructed from the coordinate format alone | supplied externally | one-frame view/model where useful | `Trajectory` is natural |

"Natural" means that the source scope directly carries the information needed for
that canonical object. "Lossy projection" means the source carries richer state,
usually geometry and/or system hierarchy, that is intentionally discarded.

Formats that contain only coordinates are a separate capability class. They must
not invent molecular chemistry merely to satisfy this table. XTC/DCD-like formats
normally require an external `Topology`; XYZ-like formats require an explicit
connectivity interpretation policy before they can become chemically meaningful
`Model` values.

### Independently interpretable source scopes

A `Document` preserves source-format syntax and container organization. It may
retain source locations, metadata, coordinates, unsupported records, data blocks,
or other information required for faithful interpretation and diagnostics. It is
not itself a canonical chemistry object.

The format's native independently interpretable scope must remain explicit:

```text
SMILES
  SmilesDocument
    -> SmilesInterpretation

Molfile
  MolfileDocument
    -> MolfileInterpretation

SDF
  SdfDocument
    -> records: Vec<SdfRecord>
         -> SdfRecordInterpretation

mmCIF
  MmcifDocument
    -> blocks: Vec<MmcifBlock>
         -> MmcifInterpretation         (one selected coordinate model)
         -> MmcifEnsembleInterpretation (several compatible coordinate models)
```

Formats that intrinsically represent one record do not need a synthetic public
`Record` wrapper merely for uniformity. Record-oriented formats should expose
records, and block-oriented formats should expose blocks. Kekule should not force
all formats through one generic `Document`/`Record` trait with unsupported
operations.

`SdfRecord` is the independently interpretable SDF unit. `MmcifBlock` is the
independently interpretable CIF/mmCIF data-block unit. Sibling SDF records and
sibling mmCIF blocks are independent source scopes and must not be silently merged
into one `Topology`, `Model`, or `Ensemble` merely because they occur in the same
file.

### Consistent parse API

Text formats should use the same parse vocabulary:

```rust
smiles::parse_str(input)?
smiles::parse_str_with_options(input, options)?

molfile::parse_str(input)?
molfile::parse_str_with_options(input, options)?

sdf::parse_str(input)?
sdf::parse_str_with_options(input, options)?

mmcif::parse_str(input)?
mmcif::parse_str_with_options(input, options)?
```

The no-options form uses the format's documented defaults. Explicit options are
available through the `_with_options` form. Future byte-oriented formats should
follow the analogous `parse_bytes` / `parse_bytes_with_options` convention where
appropriate.

Parsing only parses. It does not publish canonical molecules, run chemical
perception, choose a main component, or perform a format-independent modelling
workflow.

### Interpretation is the richest result

Interpretation translates source assertions into canonical Kekule state. The
interpretation object should retain all successfully interpreted canonical state
needed for its format scope plus format-specific reports, mappings, provenance,
and metadata sidecars.

For a geometry-free format such as SMILES, the richest canonical state is the
source-ordered connected molecular components plus interpretation diagnostics.

For a one-realization geometry-bearing scope such as a Molfile, SDF record, or one
selected mmCIF coordinate model, the richest canonical state is a `Model` (or
state exactly equivalent to `Topology + Positions` plus realization properties)
alongside the format-specific report/metadata. Geometry-independent outputs are
projections from that same interpreted state.

In particular, `SdfRecordInterpretation` must retain geometry. It must not eagerly
collapse to only `Vec<Molecule>` plus SDF data fields and thereby make
`to_model()` require a second interpretation path. Conceptually its shape is:

```text
SdfRecordInterpretation
  canonical Model                 # topology + matching Positions
  SDF title/data fields           # source metadata sidecar
  interpretation report/mappings
```

The exact physical field layout is an implementation detail, but after one SDF
record has been interpreted the same interpretation value must be sufficient to
inspect or obtain its molecules, topology, model, metadata, and report without
reinterpreting the source.

For a scope containing several compatible realizations of one topology, such as
multiple coordinate models within one mmCIF block, the richest multi-realization
result is an `Ensemble`. Multiple coordinate models are not automatically a
`Trajectory` because the source does not necessarily assign temporal semantics.

### Borrowed accessors and owned projections

Interpretation APIs should consistently distinguish non-consuming access from
owned projection using Kekule's established naming convention:

```text
interpretation.model()        -> borrowed Model access
interpretation.topology()     -> borrowed Topology access
interpretation.molecules()    -> borrowed/iterated Molecule access

interpretation.to_model()     -> consume/project to owned Model
interpretation.to_topology()  -> consume/project to shared-owned/owned Topology
interpretation.to_molecules() -> consume/project to owned Vec<Molecule>
```

The non-`to_` family leaves the interpretation available so callers may continue
to inspect reports, mappings, provenance, and metadata. The consuming `to_*`
family is for callers that are finished with the format-specific interpretation
wrapper and want to retain only a canonical Kekule object.

The exact ownership type of an owned topology projection may follow Kekule's
shared-topology architecture, for example `Arc<Topology>`, rather than forcing an
expensive independent topology clone. The semantic distinction is borrowed versus
owned/shared-owned access, not the spelling of the smart pointer.

### Format-level ergonomic conversions

The common path should be concise while remaining a composition of parse and
interpret rather than a second code path.

For SMILES:

```rust
let molecules = smiles::to_molecules("CCO.[Na+]")?;
let topology = smiles::to_topology("CCO.[Na+]")?;
```

The explicit path remains available:

```rust
let document = smiles::parse_str("CCO.[Na+]")?;
let interpretation = document.interpret()?;
let molecules = interpretation.to_molecules();
```

For SDF:

```rust
let document = sdf::parse_str(text)?;
let record = &document.records()[0];
let interpretation = record.interpret()?;

let molecules = interpretation.molecules();
let topology = interpretation.topology();
let model = interpretation.model();
```

and the independently interpretable record may expose direct conveniences:

```text
record.to_molecules()
record.to_topology()
record.to_model()
```

For mmCIF the same vocabulary applies at block scope, with explicit
interpretation options where model/alternate-location policy is required:

```rust
let document = mmcif::parse_str(text)?;
let block = &document.blocks()[0];
let interpretation = block.interpret_with_options(options)?;

let molecules = interpretation.molecules();
let topology = interpretation.topology();
let model = interpretation.model();
```

and conceptually:

```text
block.to_molecules(...)
block.to_topology(...)
block.to_model(...)
block.interpret_ensemble_with_options(...)
```

A method form on `Document`/`Record`/`Block` is preferred for ordinary navigation
once that source object already exists. Format-namespace free functions may remain
as concise whole-source conveniences or compatibility wrappers, but there must be
one authoritative implementation path beneath them.

### Component output and cardinality

The canonical molecule-producing result for a source scope that may contain
several disconnected molecular components is:

```rust
Result<Vec<Molecule>>
```

Every returned element is one valid connected molecule. Disconnected source
syntax is partitioned rather than represented as a disconnected `Molecule`.

Examples:

```text
"CCO"               -> [ethanol]
"[Na+].[Cl-]"       -> [sodium, chloride]
"CC(=O)[O-].[Na+]"  -> [acetate, sodium]
```

Component order follows deterministic source/interpreter order. Element zero does
not carry a semantic guarantee that it is the chemically "main" component. A
caller may choose the first component if that is its desired policy, or apply an
explicit largest/organic/main-component policy separately.

An owned conversion should therefore be named `to_molecules()` whenever the
source can produce several components. A strict `to_molecule()` convenience is
appropriate only when the operation either guarantees one connected molecule or
fails loudly unless exactly one component exists. No convenience may silently
select the first component.

### Geometry-bearing projections

For a one-realization geometry-bearing source scope, the canonical relationship
is:

```text
Interpretation
  richest state: Model + format report/metadata

  -> molecules() / to_molecules()
       discard geometry and system organization
       retain the same canonical connected chemistry

  -> topology() / to_topology()
       discard realization-dependent state
       retain system molecule instances, hierarchy, static properties, and order

  -> model() / to_model()
       retain the full one-realization canonical state
```

"Geometry is ignored" in a molecule/topology projection means geometry is not
retained in the resulting canonical object. It does not mean coordinates are
forbidden during interpretation. Coordinates may legitimately participate in
source-stereo normalization, alternate-location/model selection, connectivity
resolution, atom correspondence, or other format semantics before being
discarded.

The chemistry and geometry paths must share one publication pipeline. An
implementation must not independently reinterpret chemistry for `to_molecules()`,
`to_topology()`, and `to_model()`.

Detached `Positions` access is not a headline parsing workflow. `Positions` is
deliberately topology-agnostic dense storage whose semantic meaning depends on
the matching topology order. Callers may construct a `Model` explicitly from a
shared `Topology` and matching `Positions` through the canonical model
constructor, but format APIs should normally return/project the complete `Model`
rather than encourage independently detached topology and position extraction.

### Projection invariants

All projections from one interpretation under one interpretation policy must be
mutually consistent.

For any one-realization source interpretation:

```text
interpretation.model().topology()
    has the same complete static layout as
interpretation.topology()

interpretation.to_model().topology()
    has the same complete static layout as
interpretation.to_topology()

interpretation.molecules()
    corresponds exactly to the molecule instances of interpretation.topology()
    in authoritative instance/source order
```

Equivalent consuming forms must preserve the same relationship. No projection
may silently drop a connected component, select a "main" molecule, reorder
components inconsistently, synthesize different chemistry, or run a different
perception policy.

If a format constructs or preserves hierarchy, the topology obtained directly
from the interpretation and the topology inside its model must carry the same
hierarchy. `to_topology()` must not construct a bare topology while `to_model()`
secretly adds hierarchy.

### Multi-record and multi-block containers

Multi-record formats preserve record boundaries rather than flattening all
components from an entire source into one undifferentiated vector.

For SDF:

```text
SdfDocument
  records: Vec<SdfRecord>
```

`SdfRecord` is independently interpretable. `SdfDocument` must not expose a
conversion that interprets all of its records as one `Model`; independent SDF
records are a collection, not molecule instances of one spatial system.
Whole-document conveniences may return one interpretation/result per record, but
their semantics must remain explicitly record-preserving.

An mmCIF source may likewise contain multiple independent data blocks:

```text
MmcifDocument
  blocks: Vec<MmcifBlock>
```

Sibling blocks must not be combined merely because they occur in one file.
Document-level exact-one-structural-block conveniences may remain ergonomic, but
zero or several independently interpretable structural blocks require explicit
selection/iteration by the caller. Block-level interpretation is authoritative.

### Synthetic MOL/SDF hierarchy

Molfile and SDF do not normally provide PDB/mmCIF-style chain/residue hierarchy,
but a geometry-bearing `Model` benefits from uniform hierarchy-aware selection
and slicing.

When a Molfile or SDF record is interpreted into a `Model`/`Topology`, it should
synthesize minimal hierarchy at topology scope:

```text
one deterministic synthetic chain
  one residue per connected source component
    residue name: UNL
    atom sites -> the component's InstanceAtomId values
```

`UNL` is the conventional unknown-ligand residue name. The exact synthetic chain
identifier and residue numbering policy may be implementation-defined, but must
be deterministic and documented.

This synthetic hierarchy belongs only to the assembled `Topology`; the
underlying `Molecule` definitions remain hierarchy-free. The same hierarchy must
be present whether that topology is observed through `to_topology()` or through
the topology owned by `to_model()`.

### mmCIF coordinate models and hierarchy interpretation

Coordinate-model multiplicity lives inside one `MmcifBlock`, independently of
block multiplicity:

```text
MmcifDocument
  MmcifBlock A
    coordinate model 1
    coordinate model 2
      -> Ensemble when interpreted together with verified shared topology

  MmcifBlock B
    coordinate model 1
      -> Model when one realization is selected
```

One selected coordinate model should expose the same molecule/topology/model
projection family as an SDF record. Several selected compatible coordinate models
may be interpreted as one shared-topology `Ensemble`. Model selection and
alternate-location selection are interpretation policy and must be applied once,
then inherited by every projection from that interpretation.

mmCIF hierarchy must be reconstructed as one topology-level hierarchy, not as
independent copies attached to connected molecules. The interpretation order for
one block is conceptually:

```text
_atom_site and related mmCIF categories
    -> parse source atom/residue/chain/asymmetry identity
    -> select coordinate model and alternate locations
    -> construct asserted/inferred molecular connectivity
    -> partition into connected components
    -> publish canonical Molecule values
    -> install Molecule instances into one Topology
    -> establish source atom -> InstanceAtomId correspondence
    -> construct one Hierarchy over those InstanceAtomId values
    -> classify residues and molecule definitions
    -> attach Positions / realization Properties
    -> publish Model or shared-topology Ensemble
```

The hierarchy must preserve distinct mmCIF label and author identity where both
are present, including at least the relevant chain/asymmetry, residue/component,
sequence, insertion-code, and atom-site identifiers.

Hierarchy construction must not be restricted to polymer entities. Polymer,
branched, non-polymer, ligand/ion, and water atom-site records may all carry
hierarchical/source organization and should participate when representable.

mmCIF entity kinds such as polymer, branched, non-polymer, and water are
format-specific source semantics. They are useful evidence during interpretation
and may be retained in mmCIF provenance for faithful diagnostics or round trips,
but they do not form a parallel canonical Kekule entity-role layer. The
published topology instead carries canonical `MoleculeClass` and `ResidueClass`
values under the ordinary topology classification rules.

If one source chain/asymmetry spans several disconnected graph components, the
result must remain one hierarchy chain referencing atoms from several molecule
instances. If several source chains are covalently connected, they remain
several hierarchy chains referencing one molecule instance.

Recognized semantic realization data such as occupancy and B-factor may be
promoted into reserved atom property columns with dedicated APIs. Arbitrary
source fields are not automatically promoted into canonical properties merely
because a generic property layer exists.

### Canonical constructors are orthogonal to parsing

Canonical domain objects should provide ergonomic format-independent construction
without acquiring format-specific constructors.

In particular, `Topology` should support concise construction from one or more
already canonical connected molecules:

```rust
let topology = Topology::from_molecule(&molecule)?;
let topology = Topology::from_molecules(&molecules)?;
```

The simple constructor semantics are one explicit molecule instance per input
molecule, in input order. Definition interning/reuse remains an advanced builder
concern and must not make the ordinary constructor's scientific semantics
surprising.

Likewise the canonical relationship:

```text
Topology + Positions -> Model
```

should remain available through `Model::new(...)` (with the library's shared
ownership type for topology). This is a general domain constructor, not the
standard parsing workflow. Geometry-bearing format interpretation should normally
produce/project a complete `Model` directly.

Canonical domain objects must not accumulate format-specific constructors or
writers such as `Molecule::from_smiles(...)` or `Model::from_mmcif(...)`.
Equivalent functionality belongs in the format namespace or on the
format-specific source/interpretation types.

### Reports, metadata, and source correspondence

Canonical conversion must not require throwing away useful format diagnostics.
Interpretation may return or retain format-specific reports, source mappings,
warnings, data fields, provenance, or other sidecars alongside the canonical
objects.

Source metadata belongs in a canonical domain object only when its semantics are
part of that object's architecture or when interpretation explicitly promotes it
into the generic property layer with a well-defined owner scope and target.
Otherwise it remains attached to the format document, record, block,
interpretation result, or another explicit sidecar.

SDF data fields therefore remain SDF interpretation/source metadata unless
explicitly promoted. Known mmCIF atom-site quantities such as occupancy and
B-factor may be promoted because their canonical semantics and realization scope
are defined.

### Interpretation and perception

Parsing recognizes source syntax. Interpretation translates source assertions
into canonical Kekule graph state and, where a system object is constructed,
canonical topology hierarchy and realization state.

Interpretation may perform deterministic representation rewrites required to
publish a canonical molecule, such as localization of aromatic source bonding
and conversion of source stereo notation into canonical stereo elements.

Interpretation does not run arbitrary chemical standardization, choose a
tautomer/protonation state, or invent bonds merely to force connectedness or
preserve hierarchy.

If interpretation yields multiple disconnected components, each component is
published independently as a valid `Molecule`.

Chemical perception remains a separate explicit operation. Neither molecule,
topology, model, ensemble, nor trajectory projection implicitly runs default
perception merely because a canonical object is being constructed. Requesting
geometry must not silently change the installed chemical perception relative to
a geometry-independent projection from the same interpretation.

Topology classification is separate from `Molecule::Perception`. Publishing a
`Topology` does assign the lightweight canonical `MoleculeClass` and
`ResidueClass` values defined below; it must not use classification as a reason
to run the molecule's general perception pipeline.

Perception of molecule definitions already installed in a `Topology`, `Model`,
`Ensemble`, or `Trajectory` uses the explicit owning-object operations above.
It is not part of parsing or interpretation semantics. For example, an SDF
workflow may call `let mut model = record.to_model()?; model.perceive()?;`.

## `Topology`

`Topology` is the immutable, geometry-independent representation of one
molecular system.

Its fundamental responsibility is to answer:

> Which molecular entities exist in this system, how are all of their identities
> laid out at system scope, how are their atoms organized hierarchically, and
> what broad molecular/residue classes do those entities belong to?

A topology contains one or more explicit `Molecule` instances. Because every
`Molecule` is connected and topology introduces no bonds between different
instances, the connected components of the topology's asserted covalent graph
are exactly its molecule instances.

Conceptually:

```text
Topology
  molecule definitions + MoleculeClass
  molecule instances
  topology-wide atom/bond identity
  canonical dense atom/bond ordering
  identity <-> dense-index mappings
  Hierarchy + ResidueClass
  Properties
```

Topology contains no positions, velocities, forces, periodic cell, conformers,
frame ordering, or other geometry-dependent state. Its properties are likewise
coordinate-independent and valid for this exact system layout.

### Molecule instances are the public system concept

Scientifically, a topology is a system containing molecules. A caller should
normally think in terms of explicit molecule instances, not storage
normalization.

The primary public abstraction is therefore one instance-qualified molecule.
A borrowed view such as:

```rust
pub struct MoleculeInstanceView<'a> {
    topology: &'a Topology,
    id: MoleculeInstanceId,
}
```

may provide ergonomic access to:

```text
instance identity
underlying Molecule
canonical MoleculeClass through the referenced definition
qualified atoms
qualified bonds
filtered hierarchy nodes touching this instance
```

The exact type name is not normative, but instance-first navigation is.

Typical APIs should make ordinary iteration natural:

```rust
for molecule in topology.molecules() {
    for (atom_id, atom) in molecule.atoms() {
        // atom_id is instance-qualified
    }
}
```

Lower-level definition/instance APIs remain useful for explicit reuse and
advanced system construction.

### Definitions are a storage/reuse mechanism

Repeated identical molecule instances should not require repeated storage of the
same geometry-independent molecular definition.

Topology may therefore intern reusable `Molecule` values as definitions:

```text
MoleculeDefinition
  owns one Molecule
  owns one MoleculeClass

MoleculeInstance
  has one MoleculeInstanceId
  references one MoleculeDefinitionId
```

For example, a box containing many water molecules may store one water
`MoleculeDefinition` and many `MoleculeInstance`s. All instances referencing one
definition necessarily expose the same `MoleculeClass`.

This definition/instance split is part of the storage architecture, but it is
not the primary scientific mental model presented to ordinary callers.

A published topology must not contain unused molecule definitions. Every
`MoleculeDefinition` must be referenced by at least one `MoleculeInstance`.

Instance-specific generic annotations do not require a parallel
`MoleculeInstanceMetadata` object. They belong in the molecule-instance
`PropertyTable` of the containing `Topology`'s `Properties`. Definition-invariant
annotations remain on the reusable `Molecule` definition instead.

### Canonical molecule and residue classification

Kekule has a small canonical topology classification layer for broad structural
navigation and format projection. It is intentionally not a general biological
role ontology. In particular, contextual concepts such as receptor, ligand,
cofactor, substrate, counterion, or structural water are not part of this layer.

The canonical molecule vocabulary is:

```rust
pub enum MoleculeClass {
    Protein,
    Dna,
    Rna,
    Carbohydrate,
    Water,
    Ion,
    SmallMolecule,
    Other,
}
```

`SmallMolecule` is the ordinary class for connected non-polymeric molecular
compounds that do not match one of the more specific classes. It is deliberately
preferred over a contextual name such as `Ligand`. `Other` is reserved for
entities for which Kekule has positive reason not to use one of the named
classes, including irreducibly mixed/conflicting cases.

`MoleculeClass` is definition-scoped state. One `MoleculeDefinition` stores one
class and every instance of that definition shares it. The class is not stored
inside the foundational `Molecule`, because classification may legitimately use
system hierarchy and source component information that exist only when the
molecule is installed in a `Topology`.

Hierarchy residues have a parallel small canonical vocabulary:

```rust
pub enum ResidueClass {
    AminoAcid,
    DnaNucleotide,
    RnaNucleotide,
    Carbohydrate,
    Water,
    Ion,
    Other,
}
```

`ResidueClass` describes broad residue/component identity, not canonicality.
Kekule does not introduce `NoncanonicalAminoAcid` or
`NoncanonicalNucleotide` variants in this foundational taxonomy. A modified
component may be recognized as its broad structural class when the evidence is
strong; otherwise conservative initial recognition may leave it as `Other`.
That does not prevent the enclosing molecule from being recognized from its
polymer connectivity.

Classification is automatic during topology publication, with explicit builder
assignment available as an override. Classification must remain lightweight and
deterministic; it is not a reason to run generic chemical perception, expensive
graph isomorphism, or a large substructure-search suite while loading a
`Topology`.

The initial residue classifier should use a conservative cascade:

```text
explicit user override
  -> exact known component/residue identity
  -> trivial unambiguous water/monoatomic-ion recognition where applicable
  -> Other
```

Known-component tables are fast paths, not the definition of the chemistry. The
initial tables should be small and conservative. Non-canonical or modified
residues that are not explicitly recognized may therefore remain `Other` in v1.

Molecule classification may use the whole connected molecular graph, residue
classes, and inter-residue connectivity. The initial precedence is:

```text
explicit user override
  -> Water
  -> Ion
  -> peptide-connected polymer       -> Protein
  -> DNA phosphodiester polymer      -> Dna
  -> RNA phosphodiester polymer      -> Rna
  -> recognized carbohydrate entity  -> Carbohydrate
  -> ordinary fallback               -> SmallMolecule

conflicting strong polymer identities -> Other
```

Water and monoatomic-ion recognition should use trivial graph/source evidence.
Protein and nucleic-acid recognition should be based on actual inter-residue
covalent linkage patterns plus recognized residues rather than require every
residue name to be canonical. For example, an `Other` residue embedded between
recognized amino-acid residues by peptide bonds does not break an otherwise
unambiguous `Protein` classification. The analogous rule applies to modified
nucleotides embedded in an otherwise unambiguous DNA or RNA backbone.

The initial implementation should prefer a single pass over residue identities
plus a linear scan of bonds between known residue atom groups. Local linkage
predicates such as peptide-bond and phosphodiester-link recognition are preferred
over unrestricted whole-molecule substructure matching. Classification should
therefore remain approximately linear in topology size and negligible compared
with structural parsing/construction.

Because inference may use topology-owned hierarchy while storage is definition-
scoped, `TopologyBuilder::build()` (or the equivalent publication boundary) is
the natural point at which unresolved classifications are finalized. Evidence
from every informative instance of one definition must resolve to one definition
class. If strong instance-derived evidence conflicts and there is no explicit
override, the conservative result is `MoleculeClass::Other`.

Source formats may provide useful classification hints or exact source-level
categories, but source-specific enums are not canonical replacements for
`MoleculeClass` or `ResidueClass`. Format adapters translate between their source
semantics and this canonical layer where appropriate and may retain exact source
provenance separately for round-trip fidelity.

### Instance-qualified molecular identity

`AtomId` and `BondId` are local to one `Molecule` definition. Once a molecule
appears in a topology, system-level molecular identity must qualify the local ID
by its molecule instance.

Conceptually:

```text
InstanceAtomId = (MoleculeInstanceId, AtomId)
InstanceBondId = (MoleculeInstanceId, BondId)
```

This prevents identity collisions when one definition is instantiated multiple
times.

Hierarchy identities are different: chains, residues, and atom sites are owned
by the topology itself and therefore already have topology scope. They must not
be represented as `(MoleculeInstanceId, local hierarchy ID)` merely because an
atom site eventually points into a molecule instance.

### Dense topology ordering

Numerical and property-column state requires a deterministic dense ordering over
the complete system.

Topology therefore owns an authoritative dense atom order and, where useful, a
dense bond order:

```text
InstanceAtomId <-> TopologyAtomIndex
InstanceBondId <-> TopologyBondIndex
```

These concepts have deliberately different roles:

```text
InstanceAtomId / InstanceBondId
  semantic system identity

TopologyAtomIndex / TopologyBondIndex
  dense storage position
```

Dense ordering tells `Properties`, `Model`, `Ensemble`, and `Trajectory` how to
interpret topology-wide atom and bond columns and numerical arrays. Detached
tables/arrays themselves do not own topology identity.

### Hierarchy ownership and invariants

`Topology` is the sole authoritative owner of `Hierarchy`.

A published topology must guarantee at least:

```text
every hierarchy chain/residue/site ID is valid within that Topology
every residue references a live chain
every residue has one canonical ResidueClass
every atom site references a live residue
every atom site resolves to one live InstanceAtomId
atom-site lookup mappings are internally consistent
hierarchy nodes may reference atoms from any molecule instance in the Topology
```

There must not be a second authoritative hierarchy stored inside molecule
definitions or instances. Molecule-centric hierarchy APIs are projections of the
topology hierarchy.

Generic chain/residue/atom-site annotations are likewise stored once at topology
scope in property tables keyed to the topology-owned hierarchy identities.

### No topology-level covalent bonds

Topology must not introduce asserted covalent/topological bonds between
molecule instances.

If atoms from two current molecule instances become connected by an asserted
bond, those atoms belong to one connected `Molecule`. The resulting system must
therefore be represented by a new topology containing the newly connected
molecule rather than by adding an inter-instance bond.

Hydrogen bonds, salt bridges, contacts, coordination hypotheses, force-field
interactions, and other spatial or energetic relations are not topology bonds.

### No connected-components API per instance

Because every published `Molecule` is connected, asking for connected components
inside one molecule instance is redundant: the answer is always exactly that
instance.

Topology should therefore expose molecule-instance membership directly rather
than retain an API such as `connected_components(instance)` whose result is
architecturally predetermined.

### Hierarchy-aware selections and subsets

Hierarchy is a primary system-navigation mechanism. `Topology` should support
selections over chains, residues, atom sites, molecule classes, residue classes,
and their identifiers/labels, with results represented as topology-bound atom
selections.

A hierarchy selection and a structural subset are distinct operations:

```text
selection
  identifies atoms in the existing Topology

subset/slice
  constructs a new Topology containing the selected atoms
```

A structural subset may cut through existing molecule instances. For each
source `Molecule`, the induced selected graph is partitioned into connected
components, and every non-empty component is published as a new valid
`Molecule` instance in the target topology. The target hierarchy is filtered and
remapped onto the resulting target `InstanceAtomId` values; empty residues and
chains are omitted.

The subset operation should return a narrow, operation-specific source-to-target
correspondence sufficient to transfer dense state such as positions,
per-entity property columns, velocities, and forces. This is not a resurrection
of a generic foundational `TopologyMapping` abstraction.

With that primitive, higher-level objects may expose ergonomic operations such
as:

```text
model.slice(selection)
ensemble.slice(selection)
trajectory.slice(selection)
```

`Model` transfers one realization; `Ensemble` and `Trajectory` construct the
subset topology once and apply the same dense-index correspondence to every
member/frame.

Property propagation through structural transformations is conservative.
Per-entity columns may be projected when a valid source-to-target entity
correspondence exists. Owner-level scalar properties are not automatically
copied to a structurally changed owner unless the operation explicitly defines
that their semantics remain valid.

A subset that creates new molecular definitions or changes residue composition
must obtain classifications valid for the resulting topology rather than blindly
copy a source definition class onto a chemically different target definition.

### Construction and invariants

A topology builder may stage definitions, instances, hierarchy, classifications,
and properties and publish an immutable `Topology` only after validation.

A published topology must satisfy at least:

```text
at least one molecule instance
every instance references a live definition
every definition is referenced by at least one instance
every definition has one canonical MoleculeClass
every referenced Molecule satisfies Molecule invariants
instance-qualified atom/bond identities are valid
dense atom/bond ordering is complete and deterministic
identity/index mappings are mutually consistent
hierarchy references and lookups are valid and complete for every stored site
every hierarchy residue has one canonical ResidueClass
all topology property tables match their target-domain cardinalities/orderings
```

Builder property mutation uses fixed-length table views. Appending identities
extends tables with missing values; replacing hierarchy cannot silently truncate
populated columns. Publication rejects remaining dimension mismatches.

High-level format-independent construction should include:

```rust
Topology::from_molecule(&molecule)?
Topology::from_molecules(&molecules)?
```

These constructors create one explicit molecule instance per input molecule in
input order. Explicit builder APIs remain available when callers want reusable
definitions and repeated instances; definition interning must not make the
ordinary constructor's scientific semantics surprising. Classification is
automatic for these ordinary constructors; callers should not have to supply a
parallel classification sidecar merely to build or serialize a normal topology.

### Immutability and topology changes

Published `Topology` is structurally immutable.

Shared exact ownership should use `Arc<Topology>` rather than cloning independent
copies of topology state.

Topology properties are annotations, not structural layout. They may be staged
through `TopologyBuilder` and exposed read-only through a shared published
`Topology`. The architecture does not require interior mutability merely to add
annotations after an `Arc<Topology>` is already shared; workflows that require a
modified static property set may construct a new topology value with the same
layout.

A chemical or structural transformation that changes molecule membership,
connectivity, atom count, bond count, hierarchy identity, or dense layout
produces a new `Topology` rather than mutating an existing topology underneath
geometry-bearing state.

Ordinary append-style extension follows the same publication rule. The primary
API should stage from the existing published topology and publish a new value:

```rust
let mut builder = topology.into_builder();
builder.add_molecule(&ligand);
let topology = builder.build()?;
```

`into_builder()` is a topology transformation boundary, not hidden mutation. For
append-only extension it should preserve the existing definitions, instances,
semantic IDs, authoritative dense order, hierarchy IDs, retained
classifications, and retained static properties, then append new identities
deterministically. A non-consuming clone-based convenience may be added later if
justified, but direct structural mutation such as `topology.add_molecule(...)`
is not the canonical API.

The core architecture does not provide a generic topology-remapping framework.
If a workflow changes topology, geometry or other dense state for the new system
must be constructed explicitly according to that workflow's own semantics.
Operation-specific correspondence, such as the mapping returned by a subset
operation, is allowed when required by that operation.

### Scope of the current Topology design

The current core intentionally remains minimal.

It does not introduce:

```text
generic provenance framework
contextual molecule-role ontology
geometry-dependent interactions
inter-molecule topology bonds
generic topology remapping
```

A unified generic property layer is part of the core because annotations already
exist across molecule, hierarchy, model, ensemble, and trajectory surfaces and
need one coherent ownership/storage model. It is not a generic provenance
framework and does not imply automatic ingestion of arbitrary source metadata.

System-level chain/residue/atom-site hierarchy and broad canonical
molecule/residue classification are not speculative metadata; they are part of
the core Topology architecture because they support structural navigation,
selection, slicing, and format projection without depending on a source-specific
sidecar.

Other concerns should not be added speculatively. They may be introduced later
only as separate concepts when concrete requirements establish their semantics.

## Properties

Kekule has one generic property architecture for extensible annotations across
its canonical molecular and structural objects. The canonical public vocabulary
is:

```text
PropertyKey
PropertyValue
PropertyColumn
PropertyTable
Properties
```

The old parallel concepts `PropMap`, `AtomData`, and `BondData` are not separate
architectural layers. Their useful semantics are folded into this one property
system.

Properties are annotations. They do not replace strongly typed core chemistry,
hierarchy, classification, geometry, or trajectory fields, and they do not
define molecular or topology identity.

### `PropertyKey`

`PropertyKey` identifies one property by name. It is a validated key type rather
than an unconstrained `String` scattered throughout the API.

The exact validation grammar is an implementation detail, but it should be
stable, deterministic, and shared by object values and columns. A conservative
ASCII identifier grammar and bounded key length are preferred.

Conceptually:

```text
PropertyKey("energy")
PropertyKey("partial_charge")
PropertyKey("source_id")
```

### `PropertyValue`

`PropertyValue` is one scalar value of one property attached directly to an
owner.

The initial value domain is deliberately small:

```rust
pub enum PropertyValue {
    Bool(bool),
    Int(i64),
    Real { value: f64, unit: Unit },
    String(String),
}
```

Real-valued properties are always unit-aware. A dimensionless real uses
`DIMENSIONLESS`; there is no parallel untyped floating-point property concept.
Stored real values must be finite.

`PropertyValue` is intended for scalar/object-level annotations such as a model
energy, method label, boolean status, or integer generation number. Large arrays
or structured analysis results should not be forced into scalar property values.

### `PropertyColumn`

`PropertyColumn` represents one property repeated over a homogeneous entity
domain in authoritative owner order.

Conceptually, columns mirror the scalar property value domain:

```rust
pub enum PropertyColumn {
    Bool(Vec<Option<bool>>),
    Int(Vec<Option<i64>>),
    Real { unit: Unit, values: Vec<Option<f64>> },
    String(Vec<Option<String>>),
}
```

`Option` permits a property to be absent for individual entities. One column has
one type, one semantic key, and for real-valued data one stored unit. Compatible
real-valued updates are converted into the stored unit. An entirely absent
column may be normalized away.

For dense topology domains, column position follows the topology's authoritative
dense order. For definition-local molecule atom/bond domains, owner APIs map
stable local identities to the corresponding property-table positions; detached
columns themselves do not resolve semantic IDs.

### `PropertyTable`

`PropertyTable` is the columnar property store for one homogeneous entity
domain.

Conceptually:

```rust
pub struct PropertyTable {
    len: usize,
    columns: BTreeMap<PropertyKey, PropertyColumn>,
}
```

Every column in one table has the table's logical length. Examples include:

```text
Molecule atom PropertyTable
Topology molecule-instance PropertyTable
Topology atom PropertyTable
Topology residue PropertyTable
Model atom PropertyTable
TrajectoryFrame bond PropertyTable
```

Columnar storage is the canonical representation for per-entity properties.
Kekule should not embed an independent map in every `Atom`, `Bond`, `Residue`,
or other repeated entity merely to attach generic annotations. Strongly typed
canonical fields such as `ResidueClass` are not generic properties and remain on
the corresponding domain object.

Stable deterministic iteration is preferred. `BTreeMap` is therefore a suitable
initial implementation unless a demonstrated performance requirement justifies a
different internal map.

### `Properties`

`Properties` is the unified storage concept for generic annotations owned at one
scope. It contains scalar values describing the owner itself and may internally
aggregate zero or more `PropertyTable`s for repeated entity domains addressed by
that owner.

Conceptually:

```text
Properties
  owner values
    PropertyKey -> PropertyValue

  entity-domain tables
    atoms              -> PropertyTable
    bonds              -> PropertyTable
    molecule_instances -> PropertyTable
    chains             -> PropertyTable
    residues           -> PropertyTable
    atom_sites         -> PropertyTable
    ... only where meaningful for that owner
```

The exact physical nesting of those tables inside `Properties` is an
implementation detail. The public API is owner-centric: `properties()` exposes
the owner-level property namespace, while repeated entity domains are exposed
directly by the owning domain object through accessors such as
`atom_properties()` and `bond_properties()`. Callers should not have to navigate
through `properties().atoms()` or drive a generic public target enum. This keeps
the valid property domains of each owner explicit and lets the owner enforce its
identity and mutation invariants.

For example:

```text
molecule.properties()
molecule.atom_properties()
molecule.bond_properties()

model.properties()
model.atom_properties()
model.bond_properties()

topology.properties()
topology.molecule_instance_properties()
topology.atom_properties()
topology.bond_properties()
topology.chain_properties()
topology.residue_properties()
topology.atom_site_properties()

ensemble.properties()
ensemble_member.properties()
ensemble_member.atom_properties()
ensemble_member.bond_properties()

trajectory.properties()
trajectory_frame.properties()
trajectory_frame.atom_properties()
trajectory_frame.bond_properties()
```

Thus `properties()` has one consistent public meaning: properties attached
directly to that owner. Per-entity property tables are reached through the
owner's semantic domain accessors, even if the implementation stores everything
inside one `Properties` value.

The full word `properties` is preferred in the public API over the abbreviation
`props`.

### Ownership and validity

Property scope is determined by the narrowest owner whose lifetime exactly
matches the property's validity.

```text
Molecule
  geometry-independent and invariant across every instance of this definition

Topology
  geometry-independent but specific to this exact system / instance layout

Model
  specific to one concrete realization

EnsembleMember
  specific to one non-temporal realization

TrajectoryFrame
  specific to one temporal realization

Ensemble / Trajectory
  collection-level annotations that apply to the collection itself
```

Target and owner are independent dimensions. An atom property may legitimately
exist at Molecule, Topology, Model, EnsembleMember, or TrajectoryFrame scope, but
those values have different validity semantics.

For example:

```text
Molecule atom property
  definition-invariant atom class

Topology atom property
  static force-field assignment specific to this assembled system

Model atom property
  geometry-dependent SASA or charge for one realization

TrajectoryFrame atom property
  instantaneous per-frame analysis value
```

A property that differs between two instances of one reusable molecular
definition cannot live on the definition's `Molecule`; it belongs at topology or
realization scope.

### Owner shapes

The intended conceptual property targets are:

```text
Molecule Properties
  owner
  atoms
  bonds

Topology Properties
  owner
  molecule_instances
  atoms
  bonds
  chains
  residues
  atom_sites

Model Properties
  owner
  atoms
  bonds

Ensemble Properties
  owner

EnsembleMember Properties
  owner
  atoms
  bonds

Trajectory Properties
  owner

TrajectoryFrame Properties
  owner
  atoms
  bonds
```

Additional target domains should be added only when a concrete owning identity
space exists and the semantics justify them.

Realization installation rejects populated molecule-instance and hierarchy
domains, even when atom and bond dimensions happen to match. Cloning a generic
property container does not implicitly promote its domains to a different scope.

### Canonical scientific state versus generic properties

The property system is an extensibility mechanism, not a replacement for the
type system.

Strongly defined core state remains available through dedicated fields/APIs, for
example:

```text
Atom.element / formal_charge / represented hydrogens
Bond.order
MoleculeDefinition.class
Residue.class
Hierarchy labels and identifiers
Positions
PeriodicCell
EnsembleMember.weight
TrajectoryFrame.time / step
TrajectoryFrame.velocities / forces
```

Canonical per-entity scientific annotations may reuse `PropertyTable` as their
physical storage when that eliminates a redundant data container. Occupancy and
B-factor are the primary initial examples: they belong to realization atom
properties, but retain dedicated semantic APIs, canonical units, validation,
and reserved names rather than becoming arbitrary user-defined strings.

Thus `AtomData` and `BondData` disappear as public architectural concepts while
their strongest implementation idea -- validated, unit-aware columnar storage --
becomes the generic `PropertyTable` / `PropertyColumn` substrate.

### Identity, equality, and perception

Generic properties do not define represented molecular chemistry.

Consequently:

```text
changing Molecule properties
  -> does not change Graph identity
  -> does not invalidate Perception
  -> does not make an otherwise identical represented Molecule chemically unequal

changing Topology properties
  -> does not change topology layout identity
  -> does not change same-layout compatibility
```

APIs that intentionally compare complete annotations may be added separately if
needed, but structural/chemical equality must not accidentally include generic
properties merely because of derived `PartialEq` on storage structs.

### Transformations and propagation

Property propagation is explicit and conservative.

When a transformation provides a valid source-to-target correspondence for a
repeated entity domain, its `PropertyTable` columns may be projected through that
correspondence. This is the natural behavior for retained atoms, bonds, residues,
or frames.

Owner-level `PropertyValue`s are not automatically valid after a structural
transformation. For example, a total energy, system label, or score attached to
a full model may no longer describe a sliced model. Such values are copied only
when the operation explicitly guarantees or defines that behavior.

The generic property layer must not invent semantic recomputation rules.

### Source metadata boundary

The existence of a generic property system does not mean arbitrary parser fields
are automatically copied into canonical objects.

```text
arbitrary source metadata
  -> stays on format Document / Record / Block / interpretation sidecar

recognized canonical or explicitly promoted data
  -> may populate a typed core field or a property with defined owner/target
```

SDF data fields therefore remain SDF record metadata unless explicitly promoted.
Known mmCIF atom-site quantities such as occupancy and B-factor may be promoted
because their canonical semantics and realization scope are defined.

## Physical quantities and units

Kekule has one runtime physical-unit architecture for the entire library. It
must not define independent "model", "trajectory", "QM", or other subsystem
unit conventions. Public APIs may accept values expressed in any compatible
`Unit`; when numerical domain state is normalized for internal storage, it uses
one library-wide canonical unit system.

The intended core remains deliberately small:

```text
BaseDimension
    -> Dimension
         -> Unit

value + Unit
    -> Quantity<T>
```

`Dimension` represents integer powers of the independent physical dimensions.
`Unit` represents one linear unit by its dimension, conversion scale, and
optional symbol. `Quantity<T>` pairs an arbitrary supported value/container with
one runtime `Unit`.

Runtime units are intentional. Kekule should not replace this architecture with
a compile-time type-level quantity system merely to encode dimensions in Rust
types. Runtime units fit parsing, serialization, dynamic APIs, foreign-language
interfaces, trajectories, potentials, and user-selected units while keeping the
core representation straightforward.

Unit conversion is explicit at boundaries. Canonicalization does not mean that
non-canonical units are second-class: `ANGSTROM`, `BOHR`, `KILOCALORIE_PER_MOLE`,
and other compatible units remain valid public inputs and outputs where useful.
The canonical system defines only the preferred numerical representation used
inside Kekule.

### Dimensional coherence at molecular scale

The canonical unit system must be both dimensionally and numerically coherent
for molecular mechanics and dynamics.

`Mass` and `Amount` remain independent base dimensions. Atomic and molecular
masses expressed in daltons are treated dimensionally as molar mass
(`Mass / Amount`), consistent with `1 Da = 1 g/mol` for unit bookkeeping. This
allows molecular-scale mechanics to compose naturally with molar energies such
as `kJ/mol`.

In particular, the canonical basis must satisfy the identity:

```text
CANONICAL_MASS_UNIT
    * CANONICAL_LENGTH_UNIT^2
    / CANONICAL_TIME_UNIT^2

== dimensionally and numerically ==

CANONICAL_ENERGY_UNIT
```

and therefore likewise produce coherent velocity, force, gradient, and force
constant units without hidden subsystem-specific conversion factors.

### Library-wide canonical unit system

The preferred canonical basis is:

```text
CANONICAL_LENGTH_UNIT          = NANOMETER
CANONICAL_MASS_UNIT            = DALTON
CANONICAL_TIME_UNIT            = PICOSECOND
CANONICAL_ENERGY_UNIT          = KILOJOULE_PER_MOLE
CANONICAL_CHARGE_UNIT          = ELEMENTARY_CHARGE
CANONICAL_TEMPERATURE_UNIT     = KELVIN
CANONICAL_ANGLE_UNIT           = RADIAN

CANONICAL_VELOCITY_UNIT        = NANOMETER / PICOSECOND
CANONICAL_FORCE_UNIT           = KILOJOULE_PER_MOLE / NANOMETER
CANONICAL_GRADIENT_UNIT        = KILOJOULE_PER_MOLE / NANOMETER
CANONICAL_FORCE_CONSTANT_UNIT  = KILOJOULE_PER_MOLE / NANOMETER^2
```

These names are library-wide. They must not use a `MODEL_` prefix, because the
canonical convention is not owned by or restricted to the `Model` type. The
same canonical units apply wherever Kekule stores the corresponding physical
quantity internally, including real-valued generic properties and property
columns when normalized storage is appropriate.

Subsystem-specific canonical unit families should not be introduced unless a
future requirement demonstrates that one shared convention is technically
insufficient. Interfaces to external tools may of course convert to whatever
units those tools require without changing Kekule's canonical system.

## Geometry boundary

Geometry starts above `Topology`.

The fundamental relationship is:

```text
Topology + Positions -> Model
```

`Positions` and the other dense realization arrays are numerical storage. They
do not own or retain `Arc<Topology>` and they do not resolve semantic atom or
bond identities themselves.

The intended separation is:

```text
Topology
  defines semantic identities, Hierarchy, classification, static Properties,
  and dense atom/bond ordering

Positions / Velocities / Forces
  store dense numerical geometry/dynamics state
  know their own shape/length and numerical units as appropriate
  do not own Topology
  do not carry topology identity

Properties / PropertyTable / PropertyColumn
  store extensible annotations at the owner scope where they are valid
  per-entity tables know their logical shape and units/types
  detached tables/columns do not own Topology identity

Model / Ensemble / Trajectory
  own the Topology context exactly once
  own realization/collection Properties as appropriate
  validate dense state and property-table dimensions against that Topology
  provide semantic atom/bond/hierarchy access when topology context is required
```

This means operations such as resolving an `InstanceAtomId` to a coordinate or a
realization atom property are operations on `Model`, `Ensemble` member/frame
views, or `Trajectory` frame views, not on a detached `Positions` array or
`PropertyTable`.

Similarly, topology compatibility is an invariant of the owning aggregate. It
is not established by storing repeated `Arc<Topology>` handles inside every
numerical or property subobject.

Geometry-dependent quantities such as positions, velocities, forces, periodic
cell, occupancies, B-factors, and other model/frame state do not belong in
`Molecule` or static `Topology` properties.

### Dense numerical and property containers

`Positions`, `Velocities`, and `Forces` are dense numerical arrays in Kekule's
library-wide canonical units. They validate numerical shape, units, and finite
values as appropriate, but are otherwise topology-agnostic.

`PropertyTable` and `PropertyColumn` provide the corresponding generic columnar
storage for extensible per-entity data. The former `AtomData` and `BondData`
public concepts are folded into this property layer rather than maintained as a
parallel storage architecture.

Primitive dense containers must not expose APIs that require a `Topology`
parameter merely to translate semantic IDs. Semantic navigation belongs to the
higher-level object that owns both topology and dense/property state.

### Canonical construction and realization projection

The canonical format-independent construction graph is:

```text
                          + Positions
                          |
Molecule(s) -> Topology --+----------------------> Model
                          |
                          + EnsembleMember(s) ----> Ensemble
                          |
                          + TrajectoryFrame(s) ---> Trajectory

Ensemble   -- select member --> EnsembleMemberView -- as_model() --> ModelView
                                              `---- to_model() --> Model

Trajectory -- select frame  --> TrajectoryFrameView -- as_model() --> ModelView
                                              `---- to_model() --> Model
```

The native vocabulary of each owning collection remains authoritative. An
`Ensemble` contains `EnsembleMember` payloads, not `Model` values. A `Trajectory`
contains `TrajectoryFrame` payloads, not `Model` values. A selected member or
frame may nevertheless be projected into model semantics because the containing
collection supplies the shared topology context.

The intended ordinary construction surface is:

| Object | Canonical construction |
| --- | --- |
| `Topology` | `Topology::from_molecule(...)`, `Topology::from_molecules(...)`, or `TopologyBuilder` |
| `Model` | `Model::new(topology, positions)` |
| `EnsembleMember` | `EnsembleMember::new(positions)` |
| `TrajectoryFrame` | `TrajectoryFrame::new(positions)` |
| `Ensemble` | `Ensemble::new(topology)` plus `push(member)`, or `Ensemble::from_members(topology, members)` |
| `Trajectory` | `Trajectory::new(topology)` plus `push(frame)`, or `Trajectory::from_frames(topology, frames)` |

The public constructors should make passing either an owned `Topology` or the
library's shared topology handle ergonomic where this does not weaken ownership
semantics. The smart-pointer spelling is an implementation detail; callers
should not need to restructure otherwise natural construction code merely to
wrap an already complete topology.

There is no special singular collection constructor. A one-member ensemble or a
one-frame trajectory uses the same member/frame construction path as any other
cardinality.

Selection from collections should bind the stored topology-free payload to the
collection's topology and return a borrowed topology-aware view:

```text
ensemble.member(i)  -> Option<EnsembleMemberView<'_>>
ensemble.members()  -> Iterator<Item = EnsembleMemberView<'_>>

trajectory.frame(i) -> Option<TrajectoryFrameView<'_>>
trajectory.frames() -> Iterator<Item = TrajectoryFrameView<'_>>
```

The topology-free `EnsembleMember` and `TrajectoryFrame` remain the payload types
used for construction and storage. Their corresponding views are borrowed
navigation/projection helpers, not additional canonical owning objects.

Detached member/frame construction requires only positions. Atom dimensions are
known from positions; bond dimensions come from explicitly supplied columns or
from the containing topology at insertion. Insertion and member replacement
establish dimensions for property tables without retained columns and validate
all populated tables against the topology, without resizing or discarding their
columns. Empty table dimensions are storage metadata, not scientific input.
Stored-member editors and frame buffers preserve their established dimensions.

Projection into model semantics uses an ownership-explicit naming convention:

```text
member.as_model() -> ModelView<'_>
member.to_model() -> Model

frame.as_model()  -> ModelView<'_>
frame.to_model()  -> Model
```

`as_model()` is a zero-copy borrowed projection. `to_model()` explicitly
materializes an owned `Model`, cloning realization state as required while
sharing the immutable topology. `ModelView` should likewise provide the explicit
owned materialization operation needed to support this convention.

The collection itself must not expose convenience methods such as
`ensemble.model(i)` or `trajectory.model(i)`: those names incorrectly suggest
that the collection contains models. The semantic operation is selection of a
member/frame followed by an optional projection to model semantics.

No Model-based `Ensemble`/`Trajectory` construction API is required by this
architecture. If such conveniences are ever justified, they must remain
secondary to the native `EnsembleMember`/`TrajectoryFrame` construction model.

## `Model`

`Model` is one concrete geometry-dependent realization of one topology.

Conceptually:

```text
Model
  shared Topology          <- owned once
  Positions
  optional periodic cell
  Properties
    owner-level values
    atom PropertyTable
    bond PropertyTable
```

A `Model` does not duplicate molecular chemistry or hierarchy. It interprets
dense realization state and property columns against its topology's
authoritative identities, dense layout, and hierarchy.

`Model` is specifically Kekule's geometry-bearing `Topology + Positions`
abstraction. Multiple coordinate models from a source format are several
realizations of the same topology and therefore belong to `Ensemble`, not to a
hierarchy-level node.

Construction validates at least:

```text
Positions length == Topology atom count
Model atom PropertyTable length == Topology atom count
Model bond PropertyTable length == Topology bond count
```

After construction, public mutation APIs must preserve those dimensional
invariants.

Canonical realization atom data such as occupancy and B-factor are represented
inside the atom property table but retain dedicated typed APIs and validation.
There is no parallel model-level `AtomData`/`BondData` architecture.

Semantic operations such as `position(InstanceAtomId)`, atom/bond property
access by semantic ID, and hierarchy-aware selection/slicing belong on `Model`
or a model-level borrowed view because only that layer owns both topology and the
realization state.

## `Ensemble`

`Ensemble` is a finite collection of non-temporal realizations of one topology.

Conceptually:

```text
Ensemble
  shared Topology          <- owned once
  collection Properties
  members[]
    Positions
    optional periodic cell
    Properties
      owner-level values
      atom PropertyTable
      bond PropertyTable
    optional weight
```

Members do not own or repeat the shared topology. Their dense state and
property columns are interpreted in the `Ensemble` topology's authoritative
order.

Insertion/construction validates every member's dimensions against the ensemble
topology. Differences between members are geometric or member-level data, not
molecular or hierarchy identity.

An ensemble weight is contextual to membership in that ensemble and therefore
belongs to the member relation rather than to `Topology`. It remains a dedicated
semantic field/API rather than an arbitrary generic property.

## `Trajectory`

`Trajectory` is an ordered temporal sequence of realizations of one topology.

Conceptually:

```text
Trajectory
  shared Topology          <- owned once
  collection Properties
  ordered frames[]
    Positions
    optional periodic cell
    Properties
      owner-level values
      atom PropertyTable
      bond PropertyTable
    optional Velocities / Forces
    optional time / step
```

Frames do not own or repeat the shared topology. Their dense arrays and property
columns are interpreted in the `Trajectory` topology's authoritative order.

Insertion, decoding, and reusable frame-buffer publication validate frame
shapes and property-table dimensions against the trajectory topology. Streaming
infrastructure may use reusable buffers for allocation efficiency, but those
buffers follow the same ownership rule: topology context is owned once by the
buffer/container rather than repeated inside every numerical/property subobject.

Trajectory readers accept a topology directly, interpreting file coordinate index
`i` as topology dense atom index `i`. They validate counts and available format
metadata (including XYZ element order) automatically. Coordinate-only data cannot
prove atom identity from equal counts. An independently supplied semantic atom
sequence may be checked against topology order without creating a persistent
assertion or binding object. A reader can create a reusable frame buffer sharing
its exact topology; publication rejects buffers belonging to another topology.

Time, step, velocities, and forces remain dedicated semantic fields/APIs rather
than being demoted into arbitrary generic properties.

A trajectory represents one fixed-topology epoch. Topology-changing chemistry or
hierarchy is not represented by silently mutating one shared topology. A workflow
with changing topology should use separate topology epochs/objects and explicitly
construct the geometry belonging to each epoch.

## Molecular identity and equality

Authoritative molecular identity is defined by represented molecular state, not
by derived cache population, generic properties, topology classification, or
system hierarchy.

`Perception` must therefore not make two otherwise identical represented
molecules unequal merely because one has different cache presence. Generic
properties likewise do not participate in represented-molecule equality.

A `Molecule` is independent of the residue/chain context in which one of its
instances appears. The same molecular definition may be instantiated in several
hierarchical contexts without becoming a different `Molecule`. `MoleculeClass`
is stored by the topology's reusable definition wrapper and does not become part
of the foundational `Molecule`'s represented chemical identity.

Topology layout equality is distinct from graph isomorphism or chemical
identity. Full topology layout equality may include molecule definitions,
instances, canonical classifications, hierarchy, semantic IDs, and dense
ordering, but installed perception and generic properties do not alter layout
compatibility. `Topology::same_layout()` therefore remains true after explicit
perception, while exact shared snapshot identity changes. APIs requiring the
same `Arc<Topology>` retain that requirement. Two independently constructed
topologies may represent chemically equivalent systems while still having
different hierarchy IDs or dense layouts.

If complete annotated-state equality is needed, it should be an explicit API
rather than an accidental consequence of deriving `PartialEq` over storage
structs containing properties.

## Persistence and reconstruction

Persistence consumers may store molecular graph/perception, properties, topology
classification, and topology hierarchy separately according to their ownership
boundaries.

Molecule reconstruction order is conceptually:

```text
Graph
  -> validate represented graph
Properties
  -> validate definition-local owner/atom/bond property dimensions
Perception
  -> checked install last
Molecule
```

Persisted disconnected graph data must be partitioned into connected molecules
or rejected before publication.

Loading must never weaken the connectedness invariant.

Topology persistence must reconstruct definitions, definition-scoped
`MoleculeClass`, instances, qualified atom/bond identities, authoritative dense
ordering, hierarchy with `ResidueClass`, and topology property tables
consistently. Hierarchy atom sites are validated against reconstructed
`InstanceAtomId` values. Geometry and realization properties are restored
separately and validated by the owning `Model`, `Ensemble`, or `Trajectory`
against that topology layout.

Runtime domain objects are not required to be generic file-format DTOs. Source
metadata that is not canonical represented molecular/topology state and has not
been explicitly promoted into a property with defined semantics should remain in
format records or other external sidecars.

## Mutation and transformations

A normal molecular edit returns one valid connected `Molecule`.

Operations whose semantic purpose is to split a molecule naturally return more
than one molecule, for example a fragmentation transformation may return
`Vec<Molecule>`.

Topology-changing system operations return a new topology rather than mutating a
published topology in place. They do not automatically remap existing dense
geometry/data unless that operation explicitly defines and returns the necessary
correspondence, as hierarchy-aware slicing does.

Per-entity property columns may follow such an explicit correspondence. Generic
owner-level properties are not blindly copied to structurally changed owners.
Canonical molecule/residue classifications must remain valid for the transformed
topology and are re-inferred when structural changes make direct preservation
unsafe.

Coordinate-only operations never mutate `Graph`, `Perception`, `Hierarchy`, or
`Topology`.

## Naming and module style

The intended molecule field/type naming is idiomatic Rust:

```rust
pub struct Molecule {
    graph: Graph,
    perception: Perception,
    properties: Properties,
}
```

`Topology` owns the system-level `Hierarchy`, canonical classification, and
`Properties`; the exact physical field/module layout is not normative.

The canonical classification vocabulary is:

```text
MoleculeClass
ResidueClass
```

The canonical property vocabulary is:

```text
PropertyKey
PropertyValue
PropertyColumn
PropertyTable
Properties
```

The public API should prefer the full word `properties` rather than `props`.
`PropMap`, `AtomData`, `BondData`, and `RealizationProperties` are not canonical
architectural names.

Field names use `snake_case`; type names use `UpperCamelCase`. Patterns such as
`graph: Graph`, `perception: Perception`, `hierarchy: Hierarchy`, and
`properties: Properties` are normal Rust style and are preferred over redundant
names unless a real ambiguity appears.

The exact file/module layout is not normative; semantic boundaries are.

## Design rules

When deciding where new state belongs:

1. Is it authoritative atom/bond/stereo chemistry of one connected molecule?
   Put it in `Graph`.
2. Is it fundamental chemistry derived from the represented molecular graph?
   Put it in `Perception`.
3. Is it an extensible annotation whose validity exactly matches one molecular
   definition, one system layout, one realization, or one collection? Put it in
   that owner's `Properties`, using an owner-level `PropertyValue` or the
   appropriate per-entity `PropertyTable`.
4. Is it the broad canonical class of a reusable molecule definition or topology
   residue? Use strongly typed `MoleculeClass` / `ResidueClass` at `Topology`
   scope rather than a generic property or source-format entity enum.
5. Does it identify which connected molecules exist in one coordinate-free
   system or define their topology-wide atom/bond layout? Put it in `Topology`.
6. Is it coordinate-independent residue/chain/polymer/atom-site organization of
   system atoms, potentially spanning molecule instances? Put it in the
   `Hierarchy` owned by `Topology`; generic annotations about those hierarchy
   nodes belong in topology property tables.
7. Is it a substantial task-specific analysis, typing, scoring, parameterization,
   or other derived result whose own data model is meaningful? Prefer a separate
   derived object. Attach selected results as properties only deliberately and
   at the scope where their validity is defined.
8. Is it dense coordinate/model/frame data? Store it in a topology-agnostic
   numerical container or the owning realization's property tables above
   `Topology`, according to its semantics.
9. Does an operation need to interpret dense data or property columns by semantic
   atom/bond/hierarchy identity? Perform it at the `Model`, `Ensemble`,
   `Trajectory`, or owning `Topology` level where the identity context exists.
10. Does an asserted new bond connect two current molecule instances? Construct a
    new connected `Molecule` and therefore a new `Topology`.
11. Does a workflow change topology? Construct the new topology and its new dense
    state explicitly; do not rely on a generic remapping layer. Narrow
    operation-specific correspondence is appropriate when required by the
    operation.
12. Is a physical quantity stored numerically inside Kekule, including as a real
    property value/column? Accept compatible units at the boundary and normalize
    consistently rather than creating a subsystem-specific unit convention.
13. Is an operation specific to a file or serialization format? Put it in that
    format namespace or on a format-specific `Document`/`Record`/`Block`, not on
    `Molecule`, `Topology`, `Model`, `Ensemble`, or `Trajectory`. Ergonomic
    helpers may compose the canonical parse/interpret pipeline but must not
    create an independent conversion path.
14. Is an arbitrary source field merely available in an input format? Do not
    automatically turn it into a generic property. Promotion requires defined
    canonical semantics, owner scope, and target domain.

The core invariants are intentionally simple:

> A Kekule `Molecule` is one connected, geometry-independent molecular entity
> represented by authoritative `Graph`, reconstructible `Perception`, and
> definition-scoped generic `Properties` that do not define chemical identity.

> A Kekule `Topology` is one immutable, geometry-independent molecular system
> composed of one or more explicit `Molecule` instances with authoritative
> topology-wide identity, dense layout, system-level `Hierarchy`, canonical
> molecule/residue classification, and system-scoped `Properties` that do not
> define layout identity.

> Every reusable `MoleculeDefinition` has one canonical `MoleculeClass`, shared
> by all of its instances, and every hierarchy `Residue` has one canonical
> `ResidueClass`. Classification is assigned automatically at topology
> publication with explicit builder overrides available for callers.

> `Hierarchy` is owned exactly once by `Topology`; it may span molecule-instance
> boundaries and maps atom sites to topology-qualified `InstanceAtomId` values.

> `Model`, `Ensemble`, and `Trajectory` each own their shared `Topology` exactly
> once. Geometry-bearing realizations use the same generic `Properties` /
> `PropertyTable` substrate instead of parallel `AtomData` and `BondData`
> architectures.

> A property is owned at the narrowest scope whose lifetime matches its validity,
> and per-entity properties are stored column-wise rather than as a separate map
> embedded in every repeated entity.

> Kekule has one library-wide canonical physical-unit system. Runtime
> `Quantity<T>` values may use any compatible unit at interfaces, but internal
> normalized numerical state and real-valued properties do not define
> independent subsystem-specific canonical unit conventions.

## Writing and export

Writing is the format-specific projection of canonical Kekule state into an
external representation. It is related to parsing and interpretation, but it is
not a guaranteed lossless inverse: an export format may be unable to represent
all canonical state owned by `Molecule`, `Topology`, `Model`, or `Ensemble`.

Writing remains format-oriented. Public write/export operations belong in
format namespaces such as `smiles`, `molfile`, `sdf`, and `mmcif`; canonical
objects must not acquire format-specific methods such as `model.write_sdf(...)`,
and Kekule must not introduce a universal `Save`/`Serializable` trait whose
implementations silently discard unsupported state.

The canonical output mapping is:

| Canonical object | SMILES | Molfile / SDF | mmCIF |
| --- | --- | --- | --- |
| `Molecule` | one connected SMILES | coordinate-free Molfile where requested | not a primary target |
| `Topology` | one dot-separated SMILES record | not directly | not directly |
| `Model` | explicit topology projection only | one Molfile / one SDF record | one block, one coordinate model |
| `[Model]` | not direct | one SDF record per model | one block per model |
| `Ensemble` | not direct | one SDF record per member | one block containing multiple coordinate models |

### SMILES projection

A `Molecule` writes as one connected SMILES. A `Topology` writes as one
 dot-separated SMILES record by serializing every explicit molecule instance in
authoritative topology instance order. Reused definitions do not collapse
repeated instances. This projection intentionally discards hierarchy, topology
properties, definition reuse, and all geometry.

`Model` and `Ensemble` should not gain direct SMILES writers merely to discard
geometry implicitly. Callers that want that projection write their topology
explicitly.

### Molfile and SDF projection

A `Model` is the natural geometry-bearing Molfile/SDF structural unit. One model
writes as one CTAB / SDF record, and that record may contain several disconnected
molecule instances from the model topology. Molfile connected-component
boundaries therefore do not imply that a writer must accept exactly one Kekule
`Molecule`.

A sequence of independent models writes to SDF as one record per input model in
input order. The models need not share topology and Kekule must not infer
ensemble semantics merely because two model topologies happen to be compatible.

An `Ensemble` writes to SDF as one record per ensemble member. Members share the
ensemble topology by type semantics, but SDF represents them as separate records.
Format-specific SDF titles and data fields remain SDF-side metadata and are not
invented or copied into canonical `Model` state. Explicit SDF record wrappers may
carry such metadata for round-trip or expert writing.

Molfile version selection is format policy. An automatic policy should prefer
V2000 when the complete record is faithfully representable and promote to V3000
when required by representational limits. An explicitly requested version must
fail rather than silently discard canonical chemistry that it cannot encode.

### mmCIF projection

For mmCIF, data-block multiplicity and coordinate-model multiplicity retain their
separate meanings:

```text
Model
  -> one data block
       one coordinate model

[Model]
  -> one independent data block per model

Ensemble
  -> one data block
       several coordinate models sharing one topology
```

This distinction is architectural. `Vec<Model>` represents independent objects;
`Ensemble` represents several non-temporal realizations of one shared topology.
Even when a sequence of models happens to contain identical topology layouts, a
writer must not reinterpret it as an ensemble.

Consequently:

```text
Vec<Model> + mmCIF -> multiple blocks
Ensemble   + mmCIF -> one multi-model block
```

Generic `Topology` carries canonical `MoleculeClass` and `ResidueClass`, not the
mmCIF-specific polymer/branched/non-polymer/water taxonomy. The mmCIF writer
should normally derive the required mmCIF entity classification automatically
from canonical topology classification plus hierarchy/structural information
where the mmCIF distinction requires it. Typical mappings include protein/DNA/RNA
to polymer, water to water, and ion/small-molecule to non-polymer; carbohydrate
may require hierarchy/connectivity to distinguish polymeric, branched, and
discrete mmCIF representations.

An explicit `MmcifEntityClassifications`-style input may remain as an expert
format-specific override or exact source-preserving aid, but it is not the
ordinary requirement for writing a normal canonical `Model` or `Ensemble`.
Source-preserving round trips may additionally reuse exact mmCIF interpretation
provenance when it carries distinctions not represented by the canonical broad
classification.

The writer must not infer biological/contextual roles such as receptor or ligand,
and it must not use the mere presence of hierarchy as a proxy for polymer status.

### Shared realization-writing path

Where practical, geometry-bearing writers should operate on borrowed model
semantics such as `ModelView` rather than require an owned `Model`. This allows
one structural-writing implementation to serve a `Model`, an
`EnsembleMemberView`, and later a trajectory-frame view without materializing
intermediate owned models.

This is an internal reuse principle, not a requirement for one generic public
writer trait. Clear format-specific implementations are preferred over
abstraction for its own sake.

Trajectory codecs remain a separate specialized capability and are outside this
canonical structural-export contract. Their frame views may reuse compatible
model-view writing machinery in the future where appropriate.

### Output sinks and errors

The foundational writer path should support streaming to an output sink where
practical, especially for multi-record SDF and multi-model mmCIF. String-returning
helpers may wrap the same authoritative implementation for ergonomic use.

Write failures are format-specific and must be explicit. Writers must reject
unsupported selected-version chemistry, unresolved required mmCIF semantics,
invalid format metadata, incompatible realization state, or I/O failures rather
than silently omit canonical state.

### Export is not native persistence

SMILES, Molfile/SDF, and mmCIF are interoperability/export formats. They are not
the versioned native persistence representation of Kekule's complete canonical
object graph. Exact persistence may eventually use a separate Kekule-native
serialization format capable of preserving definitions, canonical
classification, hierarchy, perception, generic properties, collection metadata,
and other canonical state without forcing that state through an external
chemistry format that cannot represent it.
