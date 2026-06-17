# Algorithm

## Data Model

Represent the molecular graph as slot vectors:

- `Vec<Option<Atom>>` for atom storage.
- `Vec<Option<Bond>>` for bond storage.
- `AtomId` and `BondId` wrap the slot index.
- `None` is a tombstone left by deletion.

The initial implementation should allocate IDs monotonically by appending new slots. This keeps ID lifetime easy to reason about and avoids generation counters for the first graph feature. Reusing tombstones can be reconsidered later only with documented generation or lifetime semantics.

Maintain an adjacency index for efficient neighbor and incident-bond iteration:

- `Vec<Vec<BondId>>` keyed by `AtomId` slot index is the preferred initial representation.
- Adding a bond appends the bond ID to both endpoint adjacency lists.
- Deleting a bond removes the bond ID from both endpoint adjacency lists.
- Deleting an atom deletes all incident bonds before tombstoning the atom.

The adjacency index is derived from the authoritative bond storage. If the implementation ever detects an inconsistency, tests should treat it as a bug rather than a recoverable chemistry error.

## Mutation Rules

Topology mutations are:

- Add atom.
- Delete atom.
- Add bond.
- Delete bond.
- Mutating an existing bond endpoint, if such an API is ever exposed.

Each topology mutation must call the perception invalidation hook so `Fresh` computed states become `Stale`. `Absent` remains `Absent`, and `Stale` remains `Stale`.

Atom property changes, bond property changes, and molecule property changes do not change topology. They should not invalidate topology-derived state unless a future property is documented as chemically interpreted input.

## Add Atom

1. Assign the next `AtomId` from the current atom slot length.
2. Push `Some(atom)` into atom storage.
3. Push an empty incident-bond list into adjacency storage.
4. Invalidate computed chemistry state.
5. Return the assigned ID.

## Delete Atom

1. Resolve the `AtomId`; reject missing or already-deleted IDs.
2. Collect incident `BondId`s before mutating adjacency.
3. Delete each incident bond through the same bond deletion path used by public bond deletion.
4. Replace the atom slot with `None`.
5. Clear the atom adjacency list.
6. Invalidate computed chemistry state once for the operation or allow the nested bond deletions to invalidate as well.
7. Return the removed atom.

Remaining atom IDs and bond IDs must not change.

## Add Bond

1. Resolve both atom endpoints; reject missing or deleted IDs.
2. Reject self-bonds.
3. Check existing live incident bonds for a bond connecting the same unordered endpoint pair.
4. Assign the next `BondId` from the current bond slot length.
5. Store `Bond { a, b, order, ... }`.
6. Append the new bond ID to both endpoint adjacency lists.
7. Invalidate computed chemistry state.
8. Return the assigned ID.

The graph is undirected for connectivity and duplicate detection. Bond direction, wedge information, or stereochemical annotations belong in bond fields and later perception features.

## Delete Bond

1. Resolve the `BondId`; reject missing or already-deleted IDs.
2. Remove the bond ID from both endpoint adjacency lists.
3. Replace the bond slot with `None`.
4. Invalidate computed chemistry state.
5. Return the removed bond.

## Iteration

Atom and bond iterators skip tombstones and yield IDs with references to live records.

Neighbor iteration for an atom:

1. Resolve the atom ID.
2. Traverse its adjacency list.
3. For each live bond, yield the opposite endpoint.
4. Skip no live bond silently; a missing bond referenced by adjacency is an internal invariant failure to cover with tests.

Incident-bond iteration resolves the atom and yields live `(BondId, &Bond)` pairs from the adjacency list.

## Assumptions and Edge Cases

- Molecules can contain zero atoms.
- Isolated atoms are valid.
- Parallel bonds between the same unordered atom pair are rejected in this first feature.
- Self-bonds are rejected.
- Bond order is stored but not chemically validated.
- Atom deletion cascades to incident bond deletion.
- Deleted IDs remain invalid forever under the initial monotonic allocation policy.
- Parsing APIs must not depend on this feature to perform sanitization or perception.
- RDKit and Biopython may be used later to compare behavior in validation infrastructure, but not in Rust runtime code.
