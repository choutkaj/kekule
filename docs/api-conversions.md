# Conversion API audit and migration

The workspace uses `into_editor(self)`. `to_editor()` aliases have been removed
from both `Molecule` and `MoleculeFinishError`.

The general rule is defined in `ARCHITECTURE.md`: `into_*` transfers ownership,
`to_*` retains the source while producing a converted result, and `as_*` returns
a cheap borrowed view. Ordinary getters and operations keep their semantic names.
Copyable views may take `self` without consuming their underlying owner.

## Hard breaks

The audit covers public conversion methods in all workspace crates, plus their
examples, tests, benchmarks, documentation, and private supporting conversions.
No compatibility forwarding methods or deprecated aliases remain for these APIs.

| Owner | Previous API | Current API |
| --- | --- | --- |
| `Molecule`, `MoleculeFinishError` | `to_editor()` | `into_editor()` |
| `SmilesComponentInterpretation` | `to_molecule()`, `to_parts()` | `into_molecule()`, `into_parts()` |
| `SmilesInterpretation` | `to_molecule()`, `to_molecules()`, `to_topology()`, `to_parts()` | Corresponding `into_*` methods |
| `MolfileInterpretation`, `SdfRecordInterpretation`, `MmcifInterpretation` | `to_molecules()`, `to_topology()`, `to_model()`, `to_parts()` | Corresponding `into_*` methods |
| `SdfInterpretation` | `to_records()` | `into_records()` |
| `MmcifEnsembleInterpretation` | `to_ensemble()`, `to_parts()` | `into_ensemble()`, `into_parts()` |
| `Quantity<T>` | `to_value()` | `into_value()` |
| `Quantity<T>` | `to(unit)` borrowing the source | `to_unit(unit)` borrowing the source |
| `Quantity<T>` | `to_unit(unit)` consuming the source | `into_unit(unit)` consuming the source |
| `DcdReader`, `TrrReader`, `XtcReader`, `XyzReader` | `to_indexed()` | `into_indexed()` |
| `MemoryTrajectoryWriter` | `to_trajectory()` | `into_trajectory()` |

The quantity change needs particular care: the `to_unit` spelling now denotes
the borrowed operation and requires a cloneable payload. Existing consuming
calls must migrate to `into_unit` to preserve their ownership and allocation
behavior. The two methods have different contracts; neither is an alias.

`into_*` describes ownership rather than cost. For example, indexing a reader
still scans input, and extracting molecules from a shared topology may clone
definitions. Renaming does not change parsing, chemistry, geometry, or ordering.

## Names retained after review

| API | Reason |
| --- | --- |
| Document, SDF record, and mmCIF block `to_*` conveniences | Borrow the input, interpret it, and return a new canonical result. |
| `smiles::to_molecules(&str)`, `smiles::to_topology(&str)` | Free functions with borrowed input. |
| `ModelView::to_model()`, `EnsembleMemberView::to_model()`, `TrajectoryFrameView::to_model()` | Copyable borrowed views materialize owned geometry and retain the source owner. |
| `TrajectoryFrameView::to_frame()` | Copies the complete frame payload from a borrowed view. |
| `Model::to_builder()` | Clones into staging while retaining the original model. |
| `Model::into_builder()`, `Topology::into_builder()` | Consume their owners to resume assembly. |
| `as_model()`, `Quantity::as_ref()`, `PropertyKey::as_str()` | Borrowed projections. |
| `into_parts()`, editor/builder error recovery, and existing `into_editor()` | Already express ownership transfer correctly. |

For example, `document.to_model()` retains the document, while
`document.interpret()?.into_model()` consumes only the temporary interpretation.
Similarly, `model.edit()` creates a detached draft; `model.into_editor()` moves
the model into its draft. These are distinct operations rather than aliases.

The crate-wide `clippy::wrong_self_convention` exemptions were removed. Existing
regressions exercise the renamed interpretation and reader APIs; focused public
API tests verify borrowed-source retention, consuming conversions of non-cloneable
values, storage transfer, view materialization, and editor recovery.
