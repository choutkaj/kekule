"""The checked-in list of known differences, and the rules that keep it honest.

Each entry names a task, signature patterns (`fnmatch` globs over
`path|context|transition`), a verdict, a reason, example records, and
one-way bounds per tier. A case is known only when every one of its
differences matches some entry. A run fails when:

- a case has a difference no entry matches (new);
- an entry is stale: one of its examples reached comparison and the entry
  matched nothing in that example's cases;
- an entry matched more cases than its bound for the tier, or a numeric
  difference exceeded its `max_abs_delta`.

Bounds only move down: `--update-known` lowers them to the observed counts.
"""
from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass, field
from fnmatch import fnmatchcase
from pathlib import Path

VERDICTS = {"kekule-bug", "kekule-gap", "reference-limitation", "intended-policy"}


@dataclass
class Entry:
    id: str
    task: str
    match: list[str]
    verdict: str
    reason: str
    examples: list[str]
    max_cases: dict[str, int]
    max_abs_delta: float | None = None
    matched: dict[str, int] = field(default_factory=dict)

    def applies(self, task_id: str, signature: str) -> bool:
        return fnmatchcase(task_id, self.task) and any(fnmatchcase(signature, p) for p in self.match)

    def covers(self, task_ids) -> bool:
        """Whether a run of these tasks can observe this entry at all."""
        return any(fnmatchcase(task_id, self.task) for task_id in task_ids)


def load(path: Path) -> list[Entry]:
    if not path.exists():
        return []
    with open(path, "rb") as handle:
        document = tomllib.load(handle)
    entries = []
    for value in document.get("known", []):
        entry = Entry(
            id=value["id"],
            task=value["task"],
            match=list(value["match"]),
            verdict=value["verdict"],
            reason=value["reason"],
            examples=list(value.get("examples", [])),
            max_cases=dict(value.get("max_cases", {})),
            max_abs_delta=value.get("max_abs_delta"),
        )
        if entry.verdict not in VERDICTS:
            raise ValueError(f"{entry.id}: unknown verdict {entry.verdict!r}")
        if not entry.examples:
            raise ValueError(f"{entry.id}: needs at least one example record")
        entries.append(entry)
    ids = [entry.id for entry in entries]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate known-difference ids")
    return entries


def lower_bounds(path: Path, tier: str, counts: dict[str, int]) -> list[str]:
    """Rewrite `max_cases.<tier>` downward in place; return the changed ids."""
    text = path.read_text(encoding="utf-8")
    changed = []
    blocks = re.split(r"(?m)^(?=\[\[known\]\])", text)
    for number, block in enumerate(blocks):
        match = re.search(r'(?m)^id\s*=\s*"([^"]+)"', block)
        if not match or match.group(1) not in counts:
            continue
        entry_id, observed = match.group(1), counts[match.group(1)]
        bound = re.search(rf"(?m)^(max_cases\s*=\s*\{{[^}}]*\b{tier}\s*=\s*)(\d+)", block)
        if bound and observed < int(bound.group(2)):
            blocks[number] = block[: bound.start(2)] + str(observed) + block[bound.end(2):]
            changed.append(entry_id)
    path.write_text("".join(blocks), encoding="utf-8", newline="\n")
    return changed
