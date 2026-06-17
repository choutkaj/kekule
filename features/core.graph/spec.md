# Specification

## Scope

`core.graph` defines the shared molecular graph used by `Molecule`, `SmallMolecule`, and `MacroMolecule`.

The feature must provide:

- Typed `AtomId` and `BondId` handles that are not interchangeable with raw indices.
- Stable ID behavior: deleting an atom or bond leaves a tombstone, so IDs for remaining atoms and bonds do not change.
- Storage for atoms, bonds, molecule-level properties, atom properties, bond properties, and computed perception state.
- Public APIs for adding, deleting, fetching, mutating, and iterating atoms and bonds.
- Public APIs for iterating atom neighbors and incident bonds.
- Topology mutation invalidation for computed chemistry state.
- Error reporting for invalid IDs, deleted IDs, duplicate bonds, self-bonds, and bonds whose endpoints are missing.

The same graph must be the single topology source for both small molecules and macromolecules. Biomolecular model, chain, residue, alternate-location, occupancy, B-factor, author, and label identifiers remain in `BioHierarchy`, not in core `Atom`.

## Non-goals

This feature does not implement:

- SDF, Molfile, PDB, or mmCIF parsing.
- Sanitization, valence perception, ring detection, aromaticity perception, stereochemistry assignment, canonicalization, substructure search, or molecule validation.
- Reference-data generation from RDKit or Biopython.
- Conformer coordinate storage beyond preserving any existing scaffold fields if present.
- Biomolecular hierarchy labels in core atoms.
- Runtime dependencies on RDKit, Biopython, or other reference tools.

## Public API

The implementation should expose APIs equivalent to:

- `AtomId` and `BondId` newtypes with explicit raw/index accessors for diagnostics and internal storage lookup.
- `Molecule::new()`.
- `Molecule::add_atom(atom: Atom) -> AtomId`.
- `Molecule::delete_atom(id: AtomId) -> Result<Atom>`.
- `Molecule::atom(id: AtomId) -> Result<&Atom>`.
- `Molecule::atom_mut(id: AtomId) -> Result<&mut Atom>`.
- `Molecule::atoms() -> impl Iterator<Item = (AtomId, &Atom)>`.
- `Molecule::atom_ids() -> impl Iterator<Item = AtomId>`.
- `Molecule::add_bond(a: AtomId, b: AtomId, order: BondOrder) -> Result<BondId>`.
- `Molecule::delete_bond(id: BondId) -> Result<Bond>`.
- `Molecule::bond(id: BondId) -> Result<&Bond>`.
- `Molecule::bond_mut(id: BondId) -> Result<&mut Bond>`.
- `Molecule::bonds() -> impl Iterator<Item = (BondId, &Bond)>`.
- `Molecule::bond_ids() -> impl Iterator<Item = BondId>`.
- `Molecule::neighbors(id: AtomId) -> Result<impl Iterator<Item = AtomId>>`.
- `Molecule::incident_bonds(id: AtomId) -> Result<impl Iterator<Item = (BondId, &Bond)>>`.
- `Molecule::bond_between(a: AtomId, b: AtomId) -> Result<Option<BondId>>`.
- `Molecule::atom_count()`, `Molecule::bond_count()`, and capacity/tombstone-aware helpers if needed internally.
- `props()` and `props_mut()` accessors for molecule, atom, and bond property maps.
- `perception()` accessors for computed state, with mutation invalidation performed by topology-changing operations.

Exact iterator return types may use named iterator structs if `impl Trait` is not sufficient for ergonomics or stability.

## Acceptance criteria

- `AtomId` and `BondId` remain typed newtypes and cannot be mixed accidentally.
- Deleting an atom removes incident bonds or documents and tests equivalent cleanup behavior.
- Deleting an atom or bond does not renumber any remaining atom or bond ID.
- Adding after deletion allocates a new ID or reuses tombstoned slots only if the reuse policy is explicitly documented and tested. The initial implementation should prefer monotonic ID allocation.
- All topology mutations invalidate fresh computed perception state.
- Atom and bond property maps survive ordinary access and mutation, and deleted entities do not appear in iterators.
- Neighbor and incident-bond iteration returns graph-consistent results for isolated atoms, deleted atoms, deleted bonds, and multi-atom molecules.
- Duplicate undirected bonds and self-bonds are rejected with typed errors.
- Missing or deleted endpoints are rejected when adding bonds.
- `SmallMolecule` and `MacroMolecule` continue to wrap or own the same `Molecule` graph rather than duplicate topology logic.
- The feature remains marked unimplemented and unvalidated until code and validation evidence exist.
