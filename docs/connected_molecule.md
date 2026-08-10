# Connected `Molecule` transition

`Molecule` is moving from a permissive graph container to the canonical boundary for one connected chemical entity.

## Contract

- A published `Molecule` is empty, a singleton, or one connected atom/bond graph.
- `MoleculeBuilder` may be temporarily disconnected while atoms and bonds are assembled; `build()` checks connectedness.
- `MoleculeEditor` works on a private clone and may be temporarily disconnected; `commit()` is transactional and succeeds only for a connected final graph.
- `Topology` remains the system-level container for multiple molecule instances.
- Format `Document` types preserve raw file/record structure even when it contains multiple disconnected chemical components.

## SMILES

A dot-separated SMILES remains one `SmilesDocument`, but interpretation is component-aware:

```text
CC(=O)[O-].[Na+]
    -> SmilesDocument
    -> SmilesInterpretation
       |- acetate SmallMolecule
       `- sodium SmallMolecule
```

Each `SmilesComponentInterpretation` owns one connected `SmallMolecule` and component-local atom/bond mappings back to the original document offsets. `SmallMolecule::from_smiles` is intentionally a single-component convenience and rejects dot-separated input.

## Structural formats

Molfile/SDF and mmCIF need to respect the same published-domain invariant. mmCIF currently has an important transitional case: polymer imports can carry incomplete asserted bond topology while residue-template bonds remain pending. This PR must not invent those bonds from coordinates merely to satisfy connectedness; the interpretation boundary needs to keep incomplete construction state separate from a finalized `Molecule`.
