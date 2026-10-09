//! Cross-cutting invariants over externally supplied molecules and structures.

#[path = "../support/smiles_contract.rs"]
mod contract;
#[path = "../support/mod.rs"]
mod support;

mod corpus;
mod molfile;
mod perception;
mod smiles;
mod smiles_contract;
mod structures;
