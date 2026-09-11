# Conversion APIs

As defined in `ARCHITECTURE.md`, `into_*` transfers ownership, `to_*` retains the
source while producing a converted result, and `as_*` returns a cheap borrowed
view. Ordinary getters and operations use semantic names. Copyable views may take
`self` without consuming their underlying owner.

| Operation | Contract |
| --- | --- |
| `Molecule::into_editor()`, editor error recovery | Transfer ownership into a mutable draft. |
| Interpretation `into_molecule()`, `into_molecules()`, `into_model()`, `into_parts()` | Consume the interpretation and extract canonical results. |
| Document/record/block `to_*` methods | Retain source syntax, interpret it, and produce a canonical result. |
| `smiles::to_molecules(&str)`, `smiles::to_topology(&str)` | Interpret borrowed input. |
| View `to_model()` and `TrajectoryFrameView::to_frame()` | Materialize owned geometry or frame data while retaining the owner. |
| `Model::to_builder()` | Clone into staging. |
| `Model::into_builder()`, `Topology::into_builder()` | Consume the owner to resume assembly. |
| `Quantity::to_unit(unit)` | Convert a borrowed, cloneable payload. |
| `Quantity::into_unit(unit)`, `Quantity::into_value()` | Consume the quantity; the payload need not be cloneable. |
| Reader `into_indexed()`, writer `into_trajectory()` | Transfer reader/writer ownership. |
| `as_model()`, `Quantity::as_ref()`, `PropertyKey::as_str()` | Return borrowed projections. |

Ownership does not imply cost: indexing a reader scans input, and extracting
molecules from shared topology may clone definitions. Likewise,
`document.to_model()` retains the document, while
`document.interpret()?.into_model()` consumes the temporary interpretation.
`model.edit()` creates a detached draft; `model.into_editor()` moves the model
into its draft.

Public API regressions verify source retention, conversions of non-cloneable
values, storage transfer, view materialization, and editor recovery.
