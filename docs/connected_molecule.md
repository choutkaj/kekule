# Connected `Molecule` transition

`Molecule` is moving from a permissive graph container to the canonical boundary for one connected chemical entity.

## Contract

- A published `Molecule` is intended to be empty, a singleton, or one connected atom/bond graph.
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

## mmCIF connectivity

The public mmCIF interpretation path now completes covalent connectivity from semantic format data before returning the final model:

1. existing `_struct_conn` handling supplies explicit and special covalent links;
2. `_chem_comp_bond` supplies intra-component bonds and bond orders for observed atoms;
3. ordinary polymer links are added only when two observed residues have consecutive `label_seq_id` values and the polymer type defines an unambiguous standard linkage;
4. `_pdbx_entity_branch_link`, resolved through `_pdbx_branch_scheme`, supplies covalent links for branched entities such as oligosaccharides.

Coordinate-distance inference never becomes asserted connectivity. Missing atoms are simply omitted from component-template bonds. Missing sequence positions are not bridged by a fabricated peptide or phosphodiester bond. Polymer entities declaring nonstandard linkage are treated conservatively rather than receiving guessed standard inter-residue bonds.

`MmcifInterpretationReport::template_bonds_pending()` now counts multi-atom molecule instances that remain disconnected after this authoritative completion step. Such a result means the observed structure does not contain enough authoritative connectivity to produce one connected represented graph.

The older coordinate-distance candidate diagnostic is produced by the lower interpretation stage before semantic template completion; it remains diagnostic only and must not be interpreted as the final set of missing bonds.

## Remaining macromolecular design decision

A structure can legitimately omit residues or atoms from a chemically continuous macromolecule. The represented atom graph is then disconnected even though the physical polymer is one molecule. This is the remaining blocker to making connectedness an unconditional `MacroMolecule`/`Molecule` invariant.

The transition should therefore not invent bonds across unresolved gaps. Before the final invariant is enforced, Kekule needs one explicit policy for incomplete macromolecular observations: either represent disconnected observed fragments as separate molecule instances while preserving their common source/entity identity, or introduce a dedicated incomplete-structure boundary distinct from a finalized connected `Molecule`.
