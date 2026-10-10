"""Write a run's triage report: report.md, report.json and cases.jsonl."""
from __future__ import annotations

import json
import time
from pathlib import Path

from .datasets import cache_root

OUTCOMES = ("agree", "known", "new", "blocked", "both-error")


def write(summary: dict, cases: list, entries: list) -> Path:
    stamp = time.strftime("%Y%m%d-%H%M%S")
    out = cache_root() / "reports" / f"{summary['dataset']}-{summary['tier']}-{stamp}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with open(out / "cases.jsonl", "w", encoding="utf-8", newline="\n") as handle:
        for case in cases:
            if case.outcome == "agree":
                continue
            handle.write(json.dumps({
                "id": case.id, "task": case.task.id, "outcome": case.outcome, "reason": case.reason,
                "entries": case.entries, "differences": [d.to_json() for d in case.differences],
                "errors": {side: value.get("error") for side, value in
                           (("kekule", case.kekule), ("reference", case.reference)) if value.get("status") != "ok"},
            }, ensure_ascii=False, default=str) + "\n")
    (out / "report.md").write_text(markdown(summary, cases), encoding="utf-8", newline="\n")
    return out


def markdown(summary: dict, cases: list) -> str:
    by_id = {case.id: case for case in cases}
    lines = [
        f"# Benchmark: {summary['dataset']} v{summary['version']}, {summary['tier']} tier",
        "",
        f"**{'FAILED' if summary['failed'] else 'PASSED'}** in {summary['seconds']} s. "
        f"Stored references from {', '.join(summary['tools'].values()) or 'unknown tools'}.",
        "",
        "| Task | " + " | ".join(OUTCOMES) + " |",
        "| --- | " + " | ".join("---:" for _ in OUTCOMES) + " |",
    ]
    for task, counts in summary["tasks"].items():
        lines.append(f"| `{task}` | " + " | ".join(str(counts.get(o, 0)) for o in OUTCOMES) + " |")
    if summary["new"]:
        lines += ["", "## New differences", ""]
        for item in summary["new"]:
            lines.append(f"- **{item['cases']}** × `{item['task']}` `{item['signature']}`")
            for case_id in item["examples"]:
                case = by_id[case_id]
                shown = [d for d in case.differences if d.signature == item["signature"]][:3]
                values = "; ".join(f"`{list(d.key)}` reference {short(d.reference)} kekule {short(d.kekule)}"
                                   for d in shown)
                lines.append(f"  - `{case_id}`: {values}")
    for title, key in (("Stale known differences", "stale"), ("Bound violations", "violations")):
        if summary[key]:
            lines += ["", f"## {title}", ""] + [f"- {item}" for item in summary[key]]
    if summary["known"]:
        lines += ["", "## Known differences", "", "| Entry | Verdict | Cases | Bound |", "| --- | --- | ---: | ---: |"]
        for entry, value in summary["known"].items():
            lines.append(f"| `{entry}` | {value['verdict']} | {value['cases']} | {value['bound']} |")
    if summary.get("tightened"):
        lines += ["", "Tightened bounds: " + ", ".join(summary["tightened"])]
    return "\n".join(lines) + "\n"


def short(value) -> str:
    text = json.dumps(value, ensure_ascii=False, default=lambda v: "<missing>")
    return f"`{text if len(text) <= 80 else text[:77] + '...'}`"
