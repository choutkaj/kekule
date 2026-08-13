# Connected `Molecule` transition

`Molecule` is moving from a permissive graph container to the canonical boundary for one connected represented chemical entity.

## Contract

- A published `Molecule` is empty, a singleton, or one connected atom/bond graph.
- `MoleculeBuilder` may be temporarily disconnected while atoms and bonds are assembled; `build()` checks connectedness.
- `MoleculeEditor` works on a private clone and may be temporarily disconnected; `commit()` is transactional and succeeds only for a connected final graph.
- `SmallMolecule` and `MacroMolecule` are connected molecular objects, not containers for mixtures or unresolved graph fragments.
- `Topology` remains the system-level container for multiple molecule instances.
- Source-level identity may span several represented molecules. File/entity/chain provenance is not the same thing as graph connectedness.
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

The public mmCIF interpretation path completes covalent connectivity from semantic format data before producing the final topology:

1. existing `_struct_conn` handling supplies explicit and special covalent links;
2. `_chem_comp_bond` supplies intra-component bonds and bond orders for observed atoms;
3. ordinary polymer links are added only when two observed residues have consecutive `label_seq_id` values and the polymer type defines an unambiguous standard linkage;
4. `_pdbx_entity_branch_link`, resolved through `_pdbx_branch_scheme`, supplies covalent links for branched entities such as oligosaccharides.

Coordinate-distance inference never becomes asserted connectivity. Missing atoms are simply omitted from component-template bonds. Missing sequence positions are not bridged by a fabricated peptide or phosphodiester bond. Polymer entities declaring nonstandard linkage are treated conservatively rather than receiving guessed standard inter-residue bonds.

After authoritative connectivity is applied, every residual connected component becomes its own molecule instance. For example:

```text
source mmCIF entity 1 / asym A
    observed residues 1-50
    missing residues 51-69
    observed residues 70-100

        -> Topology
           |- MacroMolecule fragment 1  (residues 1-50)
           `- MacroMolecule fragment 2  (residues 70-100)
```

Both fragments retain the same source `_entity` / asymmetry / coordinate-model provenance. No direct bond is invented across the unresolved gap. The partition stage carries an explicit source-atom to target-atom mapping so positions, atom data, atom-site IDs, hierarchy metadata, and local identifiers are remapped deliberately rather than relying on dense-order coincidence.

`MmcifInterpretationReport::instances()` therefore describes represented connected molecule instances, while repeated source `entity_ids()` / `asym_ids()` preserve the fact that several fragments came from one source biological entity or chain. After final partitioning, `template_bonds_pending()` is zero: unresolved connectivity is represented structurally as separate molecule instances rather than as a disconnected molecule waiting for guessed chemistry.

The older coordinate-distance candidate diagnostic is produced by the lower interpretation stage before semantic template completion; it remains diagnostic only and must not be interpreted as the final set of missing bonds.

## Complete macromolecule reconstruction

If Kekule later supports reconstructing an unobserved section from sequence/templates, that should be a separate explicit operation. It may create the missing atoms and bonds and then produce a different connected `MacroMolecule` plus partial/derived coordinate state. The basic `Molecule` abstraction should not be weakened to encode atoms that are absent from the represented structure.
