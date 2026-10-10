"""Run one dataset tier: observe with Kekule, compare with stored references, triage.

Needs only the standard library and the Kekule observer. Reference
observations come from the stored file that `bench.py reference` generated.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from . import compare, errors, known, report
from .datasets import Dataset, DatasetError, Record, SourceFile, read_query
from .references import observers
from .tasks import SMALL_FORMATS, SMARTS_TARGET_STRATA, WRITTEN_FORMATS, Task, select

ROOT = Path(__file__).resolve().parents[2]
KNOWN = Path(__file__).resolve().parents[1] / "known-differences.toml"
CURATED_REFERENCES = "references-curated.jsonl.gz"
FULL_REFERENCES = "references.jsonl.gz"


@dataclass
class Case:
    id: str
    task: Task
    record: Record
    file: SourceFile
    text: str
    options: dict
    kekule: dict = field(default_factory=dict)
    reference: dict = field(default_factory=dict)
    outcome: str = ""
    reason: str = ""
    differences: list = field(default_factory=list)
    entries: list = field(default_factory=list)


def read_text(path: Path) -> str:
    data = path.read_bytes()
    if path.suffix == ".gz":
        data = gzip.decompress(data)
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode("latin-1")


def prepared_key(kind: str, path: str) -> str:
    return f"{kind}|{path}"


def reference_key(observer: str, path: str, options: dict) -> str:
    digest = hashlib.sha256(json.dumps(options, sort_keys=True).encode()).hexdigest()[:16]
    return f"{observer}|{path}|{digest}"


def build_cases(dataset: Dataset, tier: str, tasks: list[Task], directory: Path, limit: int | None,
                prepare: Callable[[str, str, str], str]) -> list[Case]:
    records = [r for r in dataset.tier(tier) if not any(f.format == "smarts" for f in r.files)]
    if limit:
        records = records[:limit]
    cases = []
    texts: dict[str, str] = {}

    def text(file: SourceFile) -> str:
        if file.path not in texts:
            texts[file.path] = read_text(directory / file.path)
        return texts[file.path]

    for t in tasks:
        if t.id == "small.smarts":
            continue
        for record in records:
            if t.strata and record.stratum not in t.strata:
                continue
            for file in record.files:
                if file.format in t.formats:
                    source = text(file)
                    if t.prepare:
                        source = prepare(t.prepare, file.path, source)
                    cases.append(Case(f"{record.id}|{file.path}|{t.id}", t, record, file, source, dict(t.options)))
    smarts = [t for t in tasks if t.id == "small.smarts"]
    if smarts:
        curated = dataset.files("curated")
        queries = [
            read_query(curated, record.files[0])[0]
            for record in dataset.tier("curated")
            if record.files[0].format == "smarts"
        ]
        for record in dataset.tier("curated"):
            if record.stratum in SMARTS_TARGET_STRATA and record.files[0].format != "smarts":
                file = record.files[0]
                source = read_text(curated / file.path)
                cases.append(Case(f"{record.id}|{file.path}|small.smarts", smarts[0], record, file, source,
                                  {"queries": queries}))
    return cases


def reference_jobs(case: Case) -> list[tuple[str, tuple]]:
    """The reference observations a case is compared against, by stored key."""
    jobs = []
    if case.file.format in SMALL_FORMATS:
        jobs.append((reference_key("rdkit:parse", case.file.path, {}),
                     ("rdkit:parse", case.file.format, case.text, {})))
    options = case.task.reference_options if case.task.reference_options is not None else case.options
    jobs.append((reference_key(case.task.reference, case.file.path, options),
                 (case.task.reference, case.file.format, case.text, options)))
    return jobs


class StoredReferences:
    """Reference observations generated once by `bench.py reference`."""

    def __init__(self, path: Path, dataset: str, version: int, inputs: str):
        if not path.exists():
            raise DatasetError(f"no stored references at {path}; run `bench.py reference` in the reference environment")
        lines = gzip.decompress(path.read_bytes()).decode("utf-8").splitlines()
        self.header = json.loads(lines[0])
        # Keys name input paths, not contents: observations of another
        # dataset revision or of replaced files would silently describe
        # different inputs.
        made_for = (self.header.get("dataset"), self.header.get("version"))
        if made_for != (dataset, version):
            raise DatasetError(f"{path} holds references for {made_for[0]} v{made_for[1]}, not {dataset} v{version}; "
                               "regenerate it")
        if self.header.get("inputs") != inputs:
            raise DatasetError(f"{path} was generated from other input files than this tier's manifest lists; "
                               "regenerate it")
        if self.header["observers"] != observers.fingerprint(self.header["dataset"]):
            raise DatasetError(f"{path} was generated by different reference observers; regenerate it")
        self.prepared: dict[str, str] = {}
        self.values: dict[str, dict] = {}
        for line in lines[1:]:
            entry = json.loads(line)
            if "prepared" in entry:
                self.prepared[entry["prepared"]] = entry["text"]
            else:
                self.values[entry["key"]] = entry["value"]

    def get(self, key: str) -> dict:
        return self.values.get(key) or {"status": "error", "error": {
            "kind": "no-reference", "message": f"no stored reference {key}; regenerate the references"}}


def progress(started: float, message: str) -> None:
    print(f"[{time.time() - started:7.1f} s] {message}", file=sys.stderr, flush=True)


def observer_binary(build: bool) -> Path:
    name = "kekule-observe.exe" if os.name == "nt" else "kekule-observe"
    binary = ROOT / "target" / "release" / name
    if build:
        subprocess.run(["cargo", "build", "--release", "--locked", "-p", "kekule-bench"], cwd=ROOT, check=True)
    return binary


def run_kekule(binary: Path, requests: list[dict], jobs: int, batch: int = 256) -> dict[str, dict]:
    """Observe every request; a batch whose process dies is bisected down to
    the crashing inputs, which are reported as `crash` errors."""
    responses: dict[str, dict] = {}
    pending = [requests[i:i + batch] for i in range(0, len(requests), batch)]
    while pending:
        chunk = pending.pop()
        payload = "".join(json.dumps(r, separators=(",", ":")) + "\n" for r in chunk)
        result = subprocess.run([str(binary), "--jobs", str(jobs)], input=payload.encode(), capture_output=True)
        if result.returncode == 0:
            for line in result.stdout.decode().splitlines():
                response = json.loads(line)
                responses[response["id"]] = response
        elif len(chunk) == 1:
            message = result.stderr.decode(errors="replace").strip().splitlines()
            responses[chunk[0]["id"]] = {"status": "error", "error": {
                "kind": "crash", "message": f"exit {result.returncode & 0xFFFFFFFF:#x}: {' '.join(message[:1])}"}}
        else:
            middle = len(chunk) // 2
            pending += [chunk[:middle], chunk[middle:]]
    return responses


def written_version(fmt: str, text: str) -> str:
    lines = text.splitlines()
    version = "V3000" if len(lines) > 3 and lines[3].rstrip().endswith("V3000") else "V2000"
    if fmt == "sdf" and not text.rstrip().endswith("$$$$"):
        return f"{version}-unterminated"
    return version


def error_difference(side: str, error: dict) -> compare.Difference:
    """A one-sided failure as a difference: `error|<kind>|ok→<kind>` when Kekule fails."""
    kind = errors.refine(side, error.get("kind", "?"), error.get("message", ""))
    reference, kekule = ("ok", kind) if side == "kekule" else (kind, "ok")
    return compare.Difference(("error",), "error", kind, reference, kekule)


def run(name: str, tier: str, patterns: list[str] | None, jobs: int, limit: int | None,
        build: bool = True, update_known: bool = False) -> int:
    started = time.time()
    dataset = Dataset(name)
    directory = dataset.files(tier)
    stored = StoredReferences(dataset.directory / CURATED_REFERENCES if tier == "curated"
                              else directory / FULL_REFERENCES, name, dataset.version, dataset.inputs_digest(tier))
    tasks = select(name, patterns)
    if any(t.requires for t in tasks):
        parse_task = select(name, ["small.parse"])[0]
        if parse_task not in tasks:
            tasks = [parse_task, *tasks]

    def prepare(kind: str, path: str, text: str) -> str:
        return stored.prepared.get(prepared_key(kind, path), "# no stored prepared input\n")

    cases = build_cases(dataset, tier, tasks, directory, limit, prepare)
    progress(started, f"{len(cases)} cases over {len({c.record.id for c in cases})} records")

    binary = observer_binary(build)
    responses = run_kekule(binary, [
        {"id": c.id, "task": c.task.observer, "format": c.file.format, "text": c.text, "options": c.options}
        for c in cases
    ], jobs)
    for case in cases:
        case.kekule = responses.get(case.id) or {"status": "error", "error": {"kind": "missing", "message": ""}}
    # Writers: Kekule reads its own output back, compared with the reference's
    # reading of the original input, coordinates included.
    written = [c for c in cases if c.task.writes and c.kekule.get("status") == "ok"]
    readback = run_kekule(binary, [
        {"id": c.id, "task": "small.parse", "format": WRITTEN_FORMATS[c.task.writes],
         "text": c.kekule["value"][0][2], "options": {"coordinates": True}}
        for c in written
    ], jobs)
    for case in written:
        text = case.kekule["value"][0][2]
        value = readback.get(case.id) or {"status": "error", "error": {"kind": "missing", "message": ""}}
        if value["status"] != "ok":
            case.kekule = {"status": "error", "error": {"kind": f"unreadable-output:{value['error']['kind']}",
                                                        "message": value["error"]["message"]}}
            continue
        version = written_version(WRITTEN_FORMATS[case.task.writes], text)
        case.kekule = {"status": "ok", "value": [*value["value"], ["version", version]], "written": text}
    progress(started, "Kekule observed")

    contexts = {}
    for case in cases:
        jobs_for = reference_jobs(case)
        if case.file.format in SMALL_FORMATS:
            contexts[case.file.path] = stored.get(jobs_for[0][0])
        case.reference = stored.get(jobs_for[-1][0])
        if case.task.writes and case.reference.get("status") == "ok":
            expected = "V3000" if case.task.writes == "molfile-v3000" else "V2000"
            case.reference = {**case.reference, "facts": [*case.reference["facts"], ["version", expected]]}

    compare_cases(cases, contexts)
    progress(started, "compared")
    # Only entries this run could observe are reported and tightened.
    task_ids = [t.id for t in tasks]
    entries = [e for e in known.load(KNOWN) if e.covers(task_ids)]
    summary = triage(cases, entries, tier)
    if update_known and KNOWN.exists():
        summary["tightened"] = known.lower_bounds(KNOWN, tier, {e.id: e.matched.get(tier, 0) for e in entries})
    summary.update({
        "dataset": name, "tier": tier, "version": dataset.version, "seconds": round(time.time() - started, 1),
        "tools": stored.header.get("tools", {}),
    })
    out = report.write(summary, cases, entries)
    print(out / "report.md", file=sys.stderr)
    return 1 if summary["failed"] else 0


def compare_cases(cases: list[Case], contexts: dict[str, dict]) -> None:
    parse_cases = {(c.record.id, c.file.path): c for c in cases if c.task.id == "small.parse"}
    for case in sorted(cases, key=lambda c: c.task.id != "small.parse"):
        context = contexts.get(case.file.path, {})
        ctx = context.get("context", {}) if context.get("status") == "ok" else {}
        kekule_ok, reference_ok = case.kekule.get("status") == "ok", case.reference.get("status") == "ok"
        parse = parse_cases.get((case.record.id, case.file.path))
        if case.task.requires and parse is not None and parse is not case:
            if parse.outcome == "both-error" or parse.reason.startswith("error"):
                case.outcome, case.reason = "blocked", "parse failed"
                continue
            touched = {d.key[0] for d in parse.differences}
            if touched & case.task.requires:
                case.outcome, case.reason = "blocked", "parse differs"
                continue
            if ctx.get("normalized") and not case.task.writes:
                case.outcome, case.reason = "blocked", "reference normalized the input"
                continue
        if not kekule_ok and not reference_ok:
            case.outcome = "both-error"
            case.reason = f"{case.kekule['error']['kind']} / {case.reference['error']['kind']}"
            continue
        if not kekule_ok:
            case.differences = [error_difference("kekule", case.kekule["error"])]
            case.reason = "error: kekule"
        elif not reference_ok:
            case.differences = [error_difference("reference", case.reference["error"])]
            case.reason = "error: reference"
        else:
            case.differences = compare.diff(case.task.id, case.kekule["value"], case.reference["facts"], ctx,
                                            case.options)
        case.outcome = "differs" if case.differences else "agree"


def triage(cases: list[Case], entries: list[known.Entry], tier: str) -> dict:
    new = defaultdict(list)
    violations = []
    for case in cases:
        if case.outcome != "differs":
            continue
        unmatched = []
        matched = set()
        for difference in case.differences:
            applying = [e for e in entries if e.applies(case.task.id, difference.signature)]
            if not applying:
                unmatched.append(difference)
            for entry in applying:
                matched.add(entry.id)
                delta = difference.delta()
                if entry.max_abs_delta is not None and delta is not None and delta > entry.max_abs_delta:
                    violations.append(f"{entry.id}: |Δ| {delta:g} exceeds {entry.max_abs_delta:g} in {case.id}")
        case.entries = sorted(matched)
        if unmatched:
            case.outcome = "new"
            for signature in sorted({d.signature for d in unmatched}):
                new[(case.task.id, signature)].append(case)
        else:
            case.outcome = "known"
        for entry in entries:
            if entry.id in matched:
                entry.matched[tier] = entry.matched.get(tier, 0) + 1
    stale = []
    by_record = defaultdict(list)
    for case in cases:
        by_record[case.record.id].append(case)
    for entry in entries:
        for example in entry.examples:
            reached = [c for c in by_record.get(example, []) if known.fnmatchcase(c.task.id, entry.task)
                       and c.outcome in ("agree", "known", "new")]
            if reached and not any(entry.id in c.entries for c in reached):
                stale.append(f"{entry.id}: example {example} no longer shows the difference")
        bound = entry.max_cases.get(tier)
        count = entry.matched.get(tier, 0)
        if bound is None:
            if count:
                violations.append(f"{entry.id}: no max_cases.{tier} bound ({count} cases)")
        elif count > bound:
            violations.append(f"{entry.id}: {count} cases exceed max_cases.{tier} = {bound}")
    tasks = defaultdict(Counter)
    for case in cases:
        tasks[case.task.id][case.outcome] += 1
    return {
        "tasks": {task: dict(counts) for task, counts in sorted(tasks.items())},
        "new": [{"task": task, "signature": signature, "cases": len(group),
                 "examples": [c.id for c in sorted(group, key=lambda c: len(c.text))[:3]]}
                for (task, signature), group in sorted(new.items(), key=lambda item: -len(item[1]))],
        "known": {e.id: {"verdict": e.verdict, "cases": e.matched.get(tier, 0), "bound": e.max_cases.get(tier)}
                  for e in entries},
        "stale": stale,
        "violations": violations,
        "failed": bool(new or stale or violations),
    }
