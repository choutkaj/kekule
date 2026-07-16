# Agent rules

These rules apply to AI agents working in this repository.

## Workflow

1. Start from one canonical feature ID under `features/`; list every affected feature ID if the work spans direct dependencies.
2. Read `ARCHITECTURE.md`, this file, and the selected feature directories before editing code.
3. Keep the change scoped. Update feature metadata or docs only when behavior, public API, or validation contracts change.
4. Add or update a regression test for every defect fix or API/behavior contract change.
5. Run applicable checks before handoff and report every command that was not run, with the reason.
6. Commit logical chunks. End every commit message with:

   ```text
   Co-authored-by: codex <codex@openai.com>
   ```

## Branching and PR strategy

Use trunk-based development with short-lived feature branches.

- `main` is the clean integration branch. Keep it buildable, tested, and aligned with `ARCHITECTURE.md`.
- Do not push feature work directly to `main`.
- Start each work branch from current `main`.
- Use one branch per feature ID or tightly scoped set of directly dependent feature IDs.
- Prefer branch names such as `codex/<feature-id>-<short-topic>`, `docs/<short-topic>`, or `hotfix/<short-topic>`.
- For broad changes and/or refactors, namethe branch accordingly.
- Open a PR for every nontrivial change.
- Use draft PRs for larger or staged work.
- Keep PRs small enough to review. Prefer stacked or sequential PRs over one large branch.
- Rebase or merge from `main` before review if the branch has drifted.
- Delete branches after merge.
- Prefer squash merge for normal feature PRs so `main` remains readable.


## Architecture guardrails

- Treat `Molecule` as the raw graph kernel; domain meaning belongs in `SmallMolecule`, `MacroMolecule`, and focused modules.
- Follow the public API shape in `ARCHITECTURE.md`. Do not add broad root re-exports or bloat the prelude casually.
- Keep parsing separate from sanitization, validation, and preparation. Never hide preparation inside parsing or default sanitization.
- Keep small-molecule chemical sanitization separate from macromolecule validation/sanitization; use separate options, reports, and errors.
- Keep biomolecular labels and structure metadata in `SmcraHierarchy`, not core `Atom` or `Bond`, unless chemically general.
- Topology or chemistry-relevant mutation must invalidate affected computed state. Failed transactional operations must leave inputs unchanged.
- Parsers must return structured errors for malformed input. Writers must reject unsupported chemistry rather than silently coercing it.
- RDKit and Biopython are validation/reference tools only, not Rust runtime dependencies.
- Algorithms must document assumptions, edge cases, and resource limits.

## Validation guardrails

- Every tracked feature must have canonical `feature.toml` and `feature.md`.
- Do not claim feature/corpus parity without current generated evidence or documented manual evidence accepted by the validation harness.
- Molecular validation fixtures must be externally supplied and provenance-pinned; toy molecules are allowed only for focused unit regressions.
- Do not weaken comparisons, remove asserted fields, delete regression tests, or regenerate goldens merely to make failures disappear.
- The dashboard is generated from feature metadata; do not hand-edit `features/DASHBOARD.html`.
- Do not claim a check, workflow, branch-protection rule, corpus result, or repository setting was verified unless it was actually inspected or run.
