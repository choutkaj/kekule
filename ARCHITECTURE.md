# Architecture

## Purpose

`molecules` is a pure-Rust foundation for cheminformatics, structural
bioinformatics, and molecular modelling. It deliberately serves both small
molecules and biomolecules without making either file-format records or
simulation-engine particles the universal data model.

The foundational chemical concept is `Molecule`: one asserted chemical entity.
The assertion is supplied by a caller, an input record, or an interpretation
policy. Graph connectedness is not an entity invariant. A salt such as
`[Na+].[Cl-]`, a coordination compound such as ferrocene, and an ordinary
connected organic compound can each be one `Molecule`. Connected-component
algorithms still report the graph's actual components.

## Canonical layers

```text
format text
    -> format-specific Document
    -> explicit interpretation + report
    -> SmallMolecule / MacroMolecule / MolecularModel
    -> explicit perception, validation, or sanitization
    -> downstream prepared System or backend object
```

These boundaries are architectural, not merely API conventions:

- Parsing recognizes and preserves format syntax. It does not construct a
  canonical molecule, sanitize chemistry, or run perception.
- Interpretation applies format semantics and documented heuristics. It creates
  canonical objects, source-to-canonical mappings, imported annotations, and a
  report. It still does not sanitize.
- Perception derives chemical state from asserted topology.
- Sanitization is an explicit transactional workflow over canonical objects.
- Modelling preparation may add force-field parameters or mechanical particles,
  but those are downstream representations and never mutate a `MolecularModel`.

## Molecules and wrappers

### `Molecule`

`Molecule` is the raw graph kernel and the asserted entity boundary. It owns:

- stable typed `AtomId` and `BondId` values; deletion leaves tombstones and IDs
  are never reused;
- atoms, bonds, adjacency, graph-adjacent stereo elements, stereo groups, and
  source stereo marks;
- optional source conformers and arbitrary properties;
- one internally consistent `PerceptionState`.

Topology facts stored directly on `Atom`, `Bond`, and stereo elements include
element, isotope, formal charge, radical state, explicit-hydrogen declarations,
atom maps, bond endpoints, `BondOrder` (including `Aromatic`), local stereo, and
source marks. `Molecule` may be disconnected, unsanitized, incomplete, or
chemically invalid.

Implicit hydrogens, ring membership and ring sets, aromatic atom/bond
membership, aromaticity model/provenance, and derived CIP descriptors are not
payload fields. They live in private optional sections of `PerceptionState` and
are available only through read-only queries. There is one installed perception
profile at a time; an alternative-model calculation remains a standalone result
until an explicit perception operation installs it.

Topology or chemistry-relevant mutation clears affected perception immediately.
There is no stale state and no public freshness flag. Property and coordinate
edits are perception-neutral. Failed transactional operations leave their input
unchanged. Imported aromatic annotations may be installed with input provenance;
full aromaticity perception replaces them with model provenance.

### `SmallMolecule`

`SmallMolecule` is the ordinary chemistry wrapper around one `Molecule`. It owns
the ergonomic small-molecule workflows while retaining `graph()` and controlled
`graph_mut()` access. Obtaining mutable access is state-neutral; the selected
`Molecule` mutation operation performs the required targeted invalidation.
`SmallMolecule::from_smiles` is an intentional
parse-then-interpret convenience and does not sanitize;
`from_smiles_sanitized` names the additional operation explicitly.

### `MacroMolecule` and `BioHierarchy`

`MacroMolecule` is one `Molecule` plus `BioHierarchy<AtomId>`. `BioHierarchy`
stores structure labels and metadata—models, chains, residues, atom sites,
author/label identifiers, alternate-location metadata, occupancy, and
B-factors. It maps structural labels to local `AtomId`s but never determines
molecule boundaries.

Small- and macromolecule sanitization/validation remain separate workflows with
separate options, reports, and errors. Chemically general algorithms should
operate on `Molecule` where practical.

## Format documents

There is no generic `Document` trait. Each format exposes a loss-preserving type
appropriate to its grammar:

- `SmilesDocument` preserves source text, tokens and spans, branches, ring and
  stereo marks, and dot-component boundaries. One SMILES record asserts one
  `SmallMolecule`; dots create disconnected graph components.
- `MolfileDocument` auto-detects V2000/V3000 and preserves headers, atom/bond and
  property records, unsupported records, and source lines. Parse and chemical
  interpretation errors are distinct.
- `SdfDocument` owns ordered `SdfRecordDocument`s, each with a raw
  `MolfileDocument` and raw data fields. Interpretation returns canonical
  `SdfRecord`s. Headers and data fields are record metadata, never injected into
  `Molecule::props`.
- `MmcifDocument` preserves blocks, scalar items, loops, missing-value markers,
  unknown categories, and source locations. Interpretation returns
  `MmcifInterpretation { model, report }`.

mmCIF interpretation selects exactly one coordinate-model ID. The default,
`RequireSingle`, rejects ambiguous multi-model input; `Select(id)` and `First`
are explicit alternatives. The selected model must provide one complete finite
position for every interpreted atom. Alternate-location selection is explicit.
Only declared covalent links merge structural instances. Inferred boundaries,
model and altloc selection, ignored models, unresolved connections, and omitted
records belong in the report/provenance.

Interpretation constructs distinct Small/Macro molecule instances and assigns
only conservative, evidence-backed roles. Exact source classification remains
available in report/provenance data.

## `MolecularModel`

`MolecularModel` is the canonical start of modelling workflows:

```text
MolecularModel
  immutable ModelTopology
    ordered MoleculeInstance values
    InstanceAtomId <-> dense ModelAtomIndex
  one complete mutable Point3 array in ModelAtomIndex order
```

Public identifiers are:

```rust
MoleculeInstanceId
InstanceAtomId { molecule: MoleculeInstanceId, atom: AtomId }
InstanceBondId { molecule: MoleculeInstanceId, bond: BondId }
ModelAtomIndex
```

`ModelTopology` owns distinct molecule instances; it never flattens multiple
entities into a synthetic `Molecule`. An instance owns either a conformer-free
`SmallMolecule` or `MacroMolecule`, plus multi-valued `MoleculeRole`s and
properties. Local `AtomId` and `BondId` values—including tombstones—survive
insertion. Qualification adds ownership without remapping local IDs. Dense model
indices exist solely for complete coordinate and gradient arrays.

The supported roles are `Polymer`, `Branched`, `NonPolymer`, `Solvent`, `Ion`,
`Ligand`, and `Cofactor`. A qualified hierarchy view maps atom-site lookup to
`InstanceAtomId`, while a standalone `MacroMolecule` continues to use local IDs.

Model construction copies one selected source conformer into authoritative model
positions and strips conformers from the stored instance payload. Source objects
remain unchanged. Empty models, empty molecules, missing positions, and
non-finite positions are rejected transactionally. Once built, topology and
instance ownership are immutable; only the complete finite position set may
change.

Potentials address topology through `InstanceAtomId`/`InstanceBondId`; gradients
are dense arrays in `ModelAtomIndex` order. Prepared topology signatures include
molecule-instance boundaries. Topology-changing future operations should return
a new model plus explicit lineage mappings rather than mutate an existing model.

Periodic cells, velocities, trajectories, reactions, and model merging are not
part of the current contract.

## Downstream prepared systems

Force-field parameters, virtual sites, Drude particles, constraints, backend
particles, electronic state, and execution-engine objects do not belong in
`MolecularModel`. A future `System`, `MMSystem`, or backend-specific prepared
object may own them, but it must provide explicit mappings between its particles
and `ModelAtomIndex`/`InstanceAtomId` and remain bound to the model topology it
was prepared from.

`molecules-dreiding` demonstrates this boundary: preparation iterates molecule
instances, performs QEq per instance, accepts eligible Small or Macro instances,
and produces a topology-bound potential without sanitizing or mutating the model.

## Public API policy

The public facade is intentionally focused:

```text
core        Molecule graph kernel and local IDs
small       SmallMolecule
bio         MacroMolecule and BioHierarchy
smiles      SmilesDocument parse/interpret and writers
molfile     MolfileDocument parse/interpret and writers
sdf         SdfDocument parse/interpret and record writers
mmcif       MmcifDocument parse/interpret
perception  explicit chemical perception and sanitization
canon       canonicalization algorithms
modeling    ModelTopology, MolecularModel, potentials, minimization
```

The prelude contains only common domain and graph types. Format internals,
specialized reports, modelling types, and expert algorithms remain in focused
namespaces. Parsing, interpretation, sanitization, preparation, and writing must
stay visibly separate in names and documentation; none may be hidden inside a
default parser.
