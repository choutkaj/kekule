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
| Collection item view `to_model()`, `payload().clone()` | Materialize an owned model or payload while retaining the collection. |
| `Model::to_builder()` | Clone into staging. |
| `Model::into_builder()`, `Topology::into_builder()` | Consume the owner to resume assembly. |
| `Quantity::to_unit(unit)` | Convert a borrowed, cloneable payload. |
| `Quantity::into_unit(unit)`, `Quantity::into_value()` | Consume the quantity; the payload need not be cloneable. |
| `Ensemble`/`Trajectory` `into_parts()`, `into_items()`; `Trajectory::into_ensemble()` | Consume a collection; members or frames move without copying. `into_ensemble` assigns equal weights. |
| Writer `into_trajectory()` | Transfer writer ownership. |
| `as_model_view()` (`AsModelView`), `Quantity::as_ref()`, `PropertyKey::as_str()` | Return borrowed projections. |

Ownership does not imply cost: opening an indexed reader scans input, and extracting
molecules from shared topology may clone definitions. Likewise,
`document.to_model()` retains the document, while
`document.interpret()?.into_model()` consumes the temporary interpretation.
`model.to_editor()` creates a detached draft; `model.into_editor()` moves the model
into its draft.

Public API regressions verify source retention, conversions of non-cloneable
values, storage transfer, view materialization, and editor recovery.
