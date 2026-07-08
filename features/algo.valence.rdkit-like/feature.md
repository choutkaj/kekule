# RDKit-like Valence Perception

## Summary

Provide conservative valence perception for common organic molecules.

## Behavior/API

- Exposes `perception::valence::{ValenceModel, ValenceReport, ValenceIssue, perceive_valence}`.
- Computes explicit valence from bond order and explicit hydrogens.
- Assigns implicit hydrogens when a common allowed valence can be selected.
- Handles imported aromatic atoms without counting each aromatic bond as a localized double bond.
- Reports unsupported elements or valence excesses instead of silently accepting them.

## Implementation Notes

- The current model covers common organic elements plus selected charged salts, boranes, simple group-1/group-2 counterions, group-13 through group-16 main-group ions and covalent centers, isolated zero-valence unsupported spectator atoms, hypervalent halogens, flexible transition-metal coordination centers, radicals, and aromatic donor/acceptor cases exposed by external PubChem validation.
- Perception state is marked fresh only after the pass completes.

## Validation

- Unit tests cover neutral organics, charged species, and valence error reporting.
- RDKit-generated goldens compare valence status, explicit valence, and implicit hydrogen assignments for external PubChem fixtures.

## Out Of Scope

- Full RDKit valence tables, organometallics, query atoms, valence tautomer handling, and sanitization orchestration.

## Revision Notes

- v1: Conservative valence perception.
- v2: Validation contract narrowed to valence-specific outputs and passes the RDKit-backed `smoke` corpus.
- v3: Add corpus-driven RDKit-compatible valence cases for charged halides, boron anions, alkali counterions, hypervalent halogens, and simple mercury salts.
- v4: Expand corpus-driven RDKit-compatible valence cases for PubChem-100 salts, silicon, phosphonium, and selected metal centers.
- v5: Generalize pubchem-1k-driven valence handling for transition-metal coordination, group-14/group-15 heavy elements, oxonium centers, chalcogens, and radicals; pubchem-1k still requires further table coverage.
- v6: Add aromatic imported-SMILES valence targets so lowercase aromatic systems sanitize with RDKit-like hydrogen counts.
- v7: Move the public expert API under the `perception::valence` facade.
- v8: Add PubChem-100k as required broad-corpus validation evidence.
- v9: Expand RDKit-like simple-ion and main-group valence support for PubChem salts while leaving actinide and coordination-heavy cases as structured unsupported chemistry.
- v10: Allow isolated unsupported atoms as zero-valence spectators so disconnected PubChem salt fragments do not block sanitization of descriptor-bearing organic components.
