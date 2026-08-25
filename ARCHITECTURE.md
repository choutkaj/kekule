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
  Graph + Perception

Topology
  one geometry-independent system made from one or more Molecule instances
  topology-wide atom/bond identity and dense ordering
  one system-level Hierarchy, which may be empty

Model
  one geometry-dependent realization of a Topology

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
```

Residue, chain, polymer, asymmetry, and atom-site organization do not belong to
`Molecule`. They belong to the system-level `Hierarchy` owned by `Topology`.

This boundary is intentional. Covalent connectedness and biological/source
hierarchy are independent partitions of atoms: one hierarchy chain may span
several disconnected `Molecule` instances, and one connected `Molecule` may span
several hierarchy chains.

`Molecule` contains no coordinates, conformers, velocities, periodic cell,
residue/chain hierarchy, or other system/geometry-dependent state.

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
polymer organization
residue/chain identifiers
component names
atom-site annotations
mappings from hierarchy atom sites to InstanceAtomId
```

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
    ResidueId -> Residue
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
remain separate derived objects unless explicitly promoted by the architecture.

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
`Molecule` architecture remains `Graph + Perception`.

### Perception installation and invalidation

Perception must always correspond to the current authoritative graph.

Graph-changing edits invalidate affected perception. The implementation should
prefer simple, safe invalidation over a complex dependency engine.

Exact reconstruction of externally stored perception may be supported through
checked installation APIs, but installation must validate references and
dimensions and must never rewrite authoritative graph chemistry.

Public convenience APIs may expose perception-backed queries directly through
`Molecule`; callers should not need to duplicate perceived flags into `Atom` or
`Bond`.

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
```

Hierarchy validation is not a `MoleculeEditor` responsibility. Hierarchy is
validated when a `Topology` is published.

Chemical perception need not be valid during editing. On successful
publication, stale perception must be discarded, recomputed, or explicitly
reinstalled through a checked path.

The same editor concept may be used for construction from scratch; a separate
public `MoleculeBuilder` is not architecturally required unless it provides
clear ergonomic value without duplicating semantics.

Public unrestricted mutable access to graph internals should not bypass the
editor and thereby bypass publication validation.

## Parsing and interpretation

### Two-stage format boundary

Kekule uses a two-stage input architecture:

```text
source text / bytes
    -> parse
format-specific Document
    -> optional format-specific Record selection
    -> interpret / canonicalize
canonical Kekule domain objects
```

Parsing preserves source-format structure and syntax. Interpretation is the
boundary that translates those source assertions into Kekule's canonical domain
model.

A `Document` is format-specific. It may retain syntax, source locations,
metadata, coordinates, unsupported records, data blocks, or other information
required for faithful interpretation and diagnostics. It is not itself a
canonical chemistry object.

Record-oriented formats should expose explicit format-specific record values.
A record is the semantic unit that can be interpreted independently. Formats
that intrinsically represent one record do not need an additional public record
wrapper merely for uniformity.

Conceptually:

```text
SmilesDocument
  one molecular record

MolfileDocument
  one molecular record

SdfDocument
  records: Vec<SdfRecord>

SdfRecord
  one independently interpretable SDF record
```

The exact representation of more complex formats such as mmCIF may follow their
native structure rather than being forced into this simple record shape.

### Format namespaces and domain-object independence

Format-specific syntax belongs to the format boundary, not to canonical domain
objects. Parsing, interpretation, writing, and ergonomic whole-source helpers
should therefore live in a format namespace or on a format-specific
`Document`/`Record` value.

The intended dependency direction is:

```text
format namespace / Document / Record
    -> parse / interpret / write
    -> Molecule / Topology / Model / Ensemble

canonical domain objects
    -> do not depend on a particular source or serialization syntax
```

Consequently, canonical domain objects should not accumulate format-specific
constructors or writers such as:

```text
Molecule::from_smiles(...)
Molecule::to_smiles(...)
Model::from_mmcif(...)
```

Equivalent functionality belongs under the corresponding format surface, for
example conceptually:

```rust
let molecules = smiles::to_molecules("CCO")?;
let text = smiles::write(&molecule)?;
```

Format-specific source objects may still expose natural conversions such as
`SdfRecord::to_molecules()`, `SdfRecord::to_model()`, or an mmCIF interpretation
result's `to_model()`, because those types already represent the format boundary.

High-level conveniences are encouraged when they materially improve ordinary
usage, but they must compose the same authoritative parse -> interpret -> publish
pipeline rather than implement a second interpretation path. Expert callers
must remain able to access the explicit `Document`/`Record`, interpretation
options, reports, mappings, and diagnostics beneath the convenience.

Convenience naming must expose semantic cardinality. If one source record can
produce several connected molecular components, the convenience should return
and be named for `Vec<Molecule>`/`to_molecules()` rather than present itself as a
singular `Molecule` constructor. Likewise, a convenience that promises one
`Model` is appropriate only where its format and selection policy actually
define one model.

### Component output and `to_molecules`

The canonical molecule-producing result for one molecular record is:

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

Component order follows deterministic source/interpreter order. Element zero
does not carry a semantic guarantee that it is the chemically "main" component.
A caller may choose the first component if that is its desired policy, or apply
an explicit largest/organic/main-component policy separately.

At the parsing/interpretation boundary, an owned conversion that returns
`Vec<Molecule>` should be named `to_molecules()`, not `to_molecule()`. A
`to_molecule()` API is appropriate only where the input type or operation
actually guarantees exactly one connected molecule.

Because `Molecule` no longer owns hierarchy, `to_molecules()` must not synthesize
residue or chain organization merely to imitate a structure file. Hierarchy is
constructed only when a system-level `Topology` is being assembled.

### Multi-record formats

Multi-record formats must preserve record boundaries rather than flattening all
components from an entire source into one undifferentiated vector.

For SDF, the intended public structure is:

```text
SdfDocument
  records: Vec<SdfRecord>

SdfRecord
  molfile representation
  SDF data fields
  source metadata / diagnostics
```

`SdfRecord` is the independently interpretable unit. It may expose
`to_molecules()` and, when its coordinates are usable, `to_model()`.

An `SdfDocument` must not expose a conversion that interprets all of its records
as one `Model`. Independent SDF records are a collection, not molecule instances
of one spatial system. Likewise, the primary document API should preserve the
record boundary rather than flatten every record's molecules. Whole-document
conveniences may return one result per record if useful, but their semantics must
remain explicitly record-preserving.

The parsed source record should simply be called `SdfRecord`; a parallel
`SdfRecordDocument` naming layer is unnecessary.

### Coordinate-bearing records and models

A coordinate-bearing record may support two distinct canonical outputs:

```text
record.to_molecules()
  -> canonical connected Molecule values
  -> source geometry is not retained in the returned domain object
  -> no synthetic hierarchy is attached to the Molecules

record.to_model()
  -> the same canonical connected Molecule values
  -> assembled as molecule instances in one Topology
  -> hierarchy assembled or synthesized at Topology scope
  -> source geometry transferred into Positions in matching canonical order
  -> one Model
```

A single coordinate-bearing record may contain several disconnected molecular
components. `to_model()` represents them as several connected `Molecule`
instances in one `Topology`, with the record's coordinates retained as the one
geometric realization.

The chemistry and geometry paths must share one interpretation/publication
pipeline. An implementation must not independently reinterpret chemistry for
`to_molecules()` and `to_model()`. Conceptually, interpretation should first
produce matched canonical component state:

```text
source record
    -> staging chemistry + optional source geometry
    -> partition connected components
    -> canonical publication / atom remapping
    -> canonical Molecule + matching optional component geometry
    -> optional Topology assembly + hierarchy construction
```

`to_molecules()` discards the published geometry and does not create hierarchy.
`to_model()` assembles the published `(Molecule, geometry)` components into
`Topology + Positions` and may construct hierarchy appropriate to the format.

Coordinates may still be consulted during canonical publication when a format's
stereochemical interpretation legitimately requires geometry; "geometry is
ignored" means it is not retained in the resulting `Molecule`, not that the
interpreter is forbidden from consulting it.

### Synthetic MOL/SDF hierarchy

Molfile and SDF do not normally provide PDB/mmCIF-style chain/residue hierarchy,
but a geometry-bearing `Model` benefits from uniform hierarchy-aware selection
and slicing.

When `MolfileDocument::to_model()` or `SdfRecord::to_model()` assembles a
`Topology`, it should synthesize minimal hierarchy at topology scope:

```text
one deterministic synthetic chain
  one residue per connected source component
    residue name: UNL
    atom sites -> the component's InstanceAtomId values
```

`UNL` is the conventional unknown-ligand residue name. The exact synthetic chain
identifier and residue numbering policy may be implementation-defined, but must
be deterministic and documented.

If later format-specific evidence supports a more specific residue/component
identity, it may replace the synthetic default. The initial architecture should
not attempt speculative ligand/ion/water classification beyond information
actually present in the source.

This synthetic hierarchy belongs only to the assembled `Topology`; the
underlying `Molecule` definitions remain hierarchy-free.

### mmCIF hierarchy interpretation

mmCIF hierarchy must be reconstructed as one topology-level hierarchy, not as
independent copies attached to connected molecules.

The interpretation order is conceptually:

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
    -> attach Positions / AtomData and publish Model or Ensemble
```

The hierarchy must preserve distinct mmCIF label and author identity where both
are present, including at least the relevant chain/asymmetry, residue/component,
sequence, insertion-code, and atom-site identifiers.

Hierarchy construction must not be restricted to polymer entities. Polymer,
branched, non-polymer, ligand/ion, and water atom-site records may all carry
hierarchical/source organization and should participate when representable.

Entity classification such as polymer/non-polymer/water must be derived from
mmCIF entity/source semantics, not inferred from whether hierarchy happens to be
present.

If one source chain/asymmetry spans several disconnected graph components, the
result must remain one hierarchy chain referencing atoms from several molecule
instances. If several source chains are covalently connected, they remain
several hierarchy chains referencing one molecule instance.

### Format-specific conversion capabilities

Concrete format document and record types should expose only conversions that
make semantic sense for that format. Kekule does not need one universal public
`Document` trait with unsupported operations.

Typical capabilities are conceptually:

```text
SmilesDocument
  -> to_molecules()

MolfileDocument
  -> to_molecules()
  -> to_model()

SdfDocument
  -> records()

SdfRecord
  -> to_molecules()
  -> to_model()

MmcifDocument
  -> model/ensemble interpretation according to explicit selection policy
```

Higher-level coordinate formats may naturally expose `Model`, `Ensemble`, or
trajectory-oriented interpretation rather than pretending that all formats have
the same conversion surface.

### Reports, metadata, and source correspondence

Canonical conversion must not require throwing away useful format diagnostics.
Interpretation may return or retain format-specific reports, source mappings,
warnings, data fields, provenance, or other sidecars alongside the canonical
objects.

Source metadata belongs in a canonical domain object only when its semantics are
part of that object's architecture. Otherwise it remains attached to the format
record, interpretation result, or another explicit sidecar.

Convenience methods such as `to_molecules()` or `to_model()` may provide the
common owned result, while lower-level interpretation APIs may continue to
expose richer reports and mappings.

### Interpretation and perception

Parsing recognizes source syntax. Interpretation translates source assertions
into canonical Kekule graph state and, where a system object is constructed,
canonical topology hierarchy.

Interpretation may perform deterministic representation rewrites required to
publish a canonical molecule, such as localization of aromatic source bonding
and conversion of source stereo notation into canonical stereo elements.

Interpretation does not run arbitrary chemical standardization, choose a
tautomer/protonation state, or invent bonds merely to force connectedness or
preserve hierarchy.

If interpretation yields multiple disconnected components, each component is
published independently as a valid `Molecule`.

Chemical perception remains a separate explicit operation. Neither
`to_molecules()` nor `to_model()` implicitly runs default perception merely
because a canonical `Molecule`, `Topology`, or `Model` is being constructed.
Requesting geometry must not silently change the installed chemical perception
relative to requesting the geometry-independent molecules from the same record.

If future workflows need an ergonomic way to perceive molecule definitions
already installed in an immutable `Topology` or `Model`, that should be designed
as an explicit perception/topology transformation. It is not part of parsing or
interpretation semantics.

## `Topology`

`Topology` is the immutable, geometry-independent representation of one
molecular system.

Its fundamental responsibility is to answer:

> Which molecular entities exist in this system, how are all of their identities
> laid out at system scope, and how are their atoms organized hierarchically?

A topology contains one or more explicit `Molecule` instances. Because every
`Molecule` is connected and topology introduces no bonds between different
instances, the connected components of the topology's asserted covalent graph
are exactly its molecule instances.

Conceptually:

```text
Topology
  molecule definitions
  molecule instances
  topology-wide atom/bond identity
  canonical dense atom/bond ordering
  identity <-> dense-index mappings
  Hierarchy
```

Topology contains no positions, velocities, forces, periodic cell, conformers,
frame ordering, or other geometry-dependent state.

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

MoleculeInstance
  has one MoleculeInstanceId
  references one MoleculeDefinitionId
```

For example, a box containing many water molecules may store one water
`MoleculeDefinition` and many `MoleculeInstance`s.

This definition/instance split is part of the storage architecture, but it is
not the primary scientific mental model presented to ordinary callers.

A published topology must not contain unused molecule definitions. Every
`MoleculeDefinition` must be referenced by at least one `MoleculeInstance`.

The minimal core does not attach a generic `MoleculeInstanceMetadata` object to
each instance. Contextual roles, annotations, or source metadata can be added
later only when concrete use cases justify their semantics.

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

Numerical and geometry-bearing state requires a deterministic dense ordering over
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

Dense ordering tells `Model`, `Ensemble`, and `Trajectory` how to interpret their
numerical arrays. The arrays themselves do not own topology identity.

### Hierarchy ownership and invariants

`Topology` is the sole authoritative owner of `Hierarchy`.

A published topology must guarantee at least:

```text
every hierarchy chain/residue/site ID is valid within that Topology
every residue references a live chain
every atom site references a live residue
every atom site resolves to one live InstanceAtomId
atom-site lookup mappings are internally consistent
hierarchy nodes may reference atoms from any molecule instance in the Topology
```

There must not be a second authoritative hierarchy stored inside molecule
definitions or instances. Molecule-centric hierarchy APIs are projections of the
topology hierarchy.

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
selections over chains, residues, atom sites, and their identifiers/labels, with
results represented as topology-bound atom selections.

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
correspondence sufficient to transfer dense state such as positions, atom data,
velocities, and forces. This is not a resurrection of a generic foundational
`TopologyMapping` abstraction.

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

### Construction and invariants

A topology builder may stage definitions, instances, and hierarchy and publish
an immutable `Topology` only after validation.

A published topology must satisfy at least:

```text
at least one molecule instance
every instance references a live definition
every definition is referenced by at least one instance
every referenced Molecule satisfies Molecule invariants
instance-qualified atom/bond identities are valid
dense atom/bond ordering is complete and deterministic
identity/index mappings are mutually consistent
hierarchy references and lookups are valid and complete for every stored site
```

Convenience construction may add one fresh definition and one instance in a
single operation. Explicit APIs may separately add a reusable definition and
then instantiate it many times.

### Immutability and topology changes

Published `Topology` is immutable.

Shared exact ownership should use `Arc<Topology>` rather than cloning independent
copies of topology state.

A chemical or structural transformation that changes molecule membership,
connectivity, atom count, bond count, hierarchy identity, or dense layout
produces a new `Topology` rather than mutating an existing topology underneath
geometry-bearing state.

The core architecture does not provide a generic topology-remapping framework.
If a workflow changes topology, geometry or other dense state for the new system
must be constructed explicitly according to that workflow's own semantics.
Operation-specific correspondence, such as the mapping returned by a subset
operation, is allowed when required by that operation.

### Scope of the current Topology design

The current core intentionally remains minimal.

It does not introduce:

```text
generic molecule-instance metadata
generic provenance framework
geometry-dependent interactions
inter-molecule topology bonds
generic topology remapping
```

System-level chain/residue/atom-site hierarchy is not speculative metadata; it
is part of the core Topology architecture because hierarchy can span molecule
boundaries and is required for structural navigation and slicing.

Other concerns should not be added speculatively. They may be introduced later
only as separate concepts when concrete requirements establish their semantics.

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
quantity internally.

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
  defines semantic identities, Hierarchy, and dense atom/bond ordering

Positions / AtomData / BondData / Velocities / Forces
  store dense numerical state
  know their own shape/length and numerical units as appropriate
  do not own Topology
  do not carry topology identity

Model / Ensemble / Trajectory
  own the Topology context exactly once
  validate dense state against that Topology
  provide semantic atom/bond/hierarchy access when topology context is required
```

This means operations such as resolving an `InstanceAtomId` to a coordinate are
operations on `Model`, `Ensemble` member/frame views, or `Trajectory` frame
views, not operations on a detached `Positions` array.

Similarly, topology compatibility is an invariant of the owning aggregate. It
is not established by storing repeated `Arc<Topology>` handles inside every
numerical subobject.

Geometry-dependent quantities such as positions, velocities, forces, periodic
cell, occupancies, B-factors, and other model/frame state do not belong in
`Molecule` or `Topology`.

### Dense numerical containers

`Positions`, `Velocities`, and `Forces` are dense numerical arrays in Kekule's
library-wide canonical units. They validate numerical shape, units, and finite
values as appropriate, but are otherwise topology-agnostic.

`AtomData` and `BondData` are likewise dense data containers rather than
Topology owners. They may retain a logical item count even when all optional
columns are absent. Whether particular atom/bond data belongs at model, member,
frame, or another level is intentionally not settled by this storage rule and
may be refined separately.

The primitive dense containers must not expose APIs that require a `Topology`
parameter merely to translate semantic IDs. Semantic navigation belongs to the
higher-level object that owns both topology and dense state.

## `Model`

`Model` is one concrete geometry-dependent realization of one topology.

Conceptually:

```text
Model
  shared Topology          <- owned once
  Positions
  optional periodic cell
  model-level AtomData
  model-level BondData
```

A `Model` does not duplicate molecular chemistry or hierarchy. It interprets
dense realization state against its topology's authoritative identities, dense
layout, and hierarchy.

`Model` is specifically Kekule's geometry-bearing `Topology + Positions`
abstraction. Multiple coordinate models from a source format are several
realizations of the same topology and therefore belong to `Ensemble`, not to a
hierarchy-level node.

Construction validates at least:

```text
Positions length == Topology atom count
AtomData logical length == Topology atom count
BondData logical length == Topology bond count
```

After construction, public mutation APIs must preserve those dimensional
invariants.

Semantic operations such as `position(InstanceAtomId)` and hierarchy-aware
selection/slicing belong on `Model` or a model-level borrowed view because only
that layer owns both the topology and the numerical state.

## `Ensemble`

`Ensemble` is a finite collection of non-temporal realizations of one topology.

Conceptually:

```text
Ensemble
  shared Topology          <- owned once
  members[]
    Positions
    optional member-level geometry/data
    optional weight
```

Members do not own or repeat the shared topology. Their dense state is interpreted
in the `Ensemble` topology's authoritative order.

Insertion/construction validates every member's dimensions against the ensemble
topology. Differences between members are geometric or member-level data, not
molecular or hierarchy identity.

An ensemble weight is contextual to membership in that ensemble and therefore
belongs to the member relation rather than to `Topology`.

## `Trajectory`

`Trajectory` is an ordered temporal sequence of realizations of one topology.

Conceptually:

```text
Trajectory
  shared Topology          <- owned once
  ordered frames[]
    Positions
    optional periodic cell
    optional AtomData / BondData
    optional Velocities / Forces
    optional time / step
    optional frame properties
```

Frames do not own or repeat the shared topology. Their dense arrays are
interpreted in the `Trajectory` topology's authoritative order.

Insertion, decoding, and reusable frame-buffer publication validate frame
shapes against the trajectory topology. Streaming infrastructure may use
reusable buffers for allocation efficiency, but those buffers follow the same
ownership rule: topology context is owned once by the buffer/container rather
than repeated inside every numerical array.

A trajectory represents one fixed-topology epoch. Topology-changing chemistry or
hierarchy is not represented by silently mutating one shared topology. A workflow
with changing topology should use separate topology epochs/objects and explicitly
construct the geometry belonging to each epoch.

## Molecular identity and equality

Authoritative molecular identity is defined by represented molecular state, not
by derived cache population or system hierarchy.

`Perception` must therefore not make two otherwise identical represented
molecules unequal merely because one has different cache presence.

A `Molecule` is independent of the residue/chain context in which one of its
instances appears. The same molecular definition may be instantiated in several
hierarchical contexts without becoming a different `Molecule`.

Topology layout equality is distinct from graph isomorphism or chemical
identity. Full topology layout equality may include molecule definitions,
instances, hierarchy, semantic IDs, and dense ordering. Two independently
constructed topologies may represent chemically equivalent systems while still
having different hierarchy IDs or dense layouts.

## Persistence and reconstruction

Persistence consumers may store molecular graph/perception and topology
hierarchy separately according to their ownership boundaries.

Molecule reconstruction order is:

```text
Graph
  -> validate represented graph
Perception
  -> checked install last
Molecule
```

Persisted disconnected graph data must be partitioned into connected molecules
or rejected before publication.

Loading must never weaken the connectedness invariant.

Topology persistence must reconstruct definitions, instances, qualified
atom/bond identities, authoritative dense ordering, and hierarchy consistently.
Hierarchy atom sites are validated against reconstructed `InstanceAtomId` values.
Geometry is restored separately and validated by the owning `Model`, `Ensemble`,
or `Trajectory` against that topology layout.

Runtime domain objects are not required to be generic file-format DTOs. Source
metadata that is not canonical represented molecular or topology state should
remain in format records or other external sidecars.

## Mutation and transformations

A normal molecular edit returns one valid connected `Molecule`.

Operations whose semantic purpose is to split a molecule naturally return more
than one molecule, for example a fragmentation transformation may return
`Vec<Molecule>`.

Topology-changing system operations return a new topology rather than mutating a
published topology in place. They do not automatically remap existing dense
geometry/data unless that operation explicitly defines and returns the necessary
correspondence, as hierarchy-aware slicing does.

Coordinate-only operations never mutate `Graph`, `Perception`, `Hierarchy`, or
`Topology`.

## Naming and module style

The intended molecule field/type naming is idiomatic Rust:

```rust
pub struct Molecule {
    graph: Graph,
    perception: Perception,
}
```

`Topology` owns the system-level `Hierarchy`; the exact physical field/module
layout is not normative.

Field names use `snake_case`; type names use `UpperCamelCase`. Patterns such as
`graph: Graph`, `perception: Perception`, and `hierarchy: Hierarchy` are normal
Rust style and are preferred over redundant names unless a real ambiguity
appears.

The exact file/module layout is not normative; semantic boundaries are.

## Design rules

When deciding where new state belongs:

1. Is it authoritative atom/bond/stereo chemistry of one connected molecule?
   Put it in `Graph`.
2. Is it fundamental chemistry derived from the represented molecular graph?
   Put it in `Perception`.
3. Does it identify which connected molecules exist in one coordinate-free
   system or define their topology-wide atom/bond layout? Put it in `Topology`.
4. Is it coordinate-independent residue/chain/polymer/atom-site organization of
   system atoms, potentially spanning molecule instances? Put it in the
   `Hierarchy` owned by `Topology`.
5. Is it task-specific analysis, typing, scoring, or parameterization? Keep it
   in a separate derived object.
6. Is it dense coordinate/model/frame data? Store it in a topology-agnostic
   numerical container above `Topology`.
7. Does an operation need to interpret dense data by semantic atom/bond/hierarchy
   identity? Perform it at the `Model`, `Ensemble`, or `Trajectory` level where
   topology and numerical state meet.
8. Does an asserted new bond connect two current molecule instances? Construct a
   new connected `Molecule` and therefore a new `Topology`.
9. Does a workflow change topology? Construct the new topology and its new dense
   state explicitly; do not rely on a generic remapping layer. Narrow
   operation-specific correspondence is appropriate when required by the
   operation.
10. Is a physical quantity stored numerically inside Kekule? Accept compatible
    units at the boundary and normalize it to the one library-wide canonical
    unit for that quantity rather than creating a subsystem-specific unit
    convention.
11. Is an operation specific to a file or serialization format? Put it in that
    format namespace or on a format-specific `Document`/`Record`, not on
    `Molecule`, `Topology`, `Model`, `Ensemble`, or `Trajectory`. Ergonomic
    helpers may compose the canonical parse/interpret pipeline but must not
    create an independent conversion path.

The core invariants are intentionally simple:

> A Kekule `Molecule` is one connected, geometry-independent molecular entity
> represented by `Graph + Perception`.

> A Kekule `Topology` is one immutable, geometry-independent molecular system
> composed of one or more explicit `Molecule` instances with authoritative
> topology-wide identity, dense layout, and system-level `Hierarchy`.

> `Hierarchy` is owned exactly once by `Topology`; it may span molecule-instance
> boundaries and maps atom sites to topology-qualified `InstanceAtomId` values.

> `Model`, `Ensemble`, and `Trajectory` each own their shared `Topology` exactly
> once. Their dense numerical subobjects do not own topology identity.

> Kekule has one library-wide canonical physical-unit system. Runtime
> `Quantity<T>` values may use any compatible unit at interfaces, but internal
> normalized numerical state does not define independent subsystem-specific
> canonical unit conventions.
