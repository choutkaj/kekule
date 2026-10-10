"""Generate the stored reference observations of a dataset (maintainers only).

Runs inside the environment from `environment.yml`. Every reference
observation the benchmark compares against is computed once here and stored;
`bench.py run` then needs only Kekule and the standard library. Regenerate
when the dataset or a reference observer changes: the run refuses stored
references made by a different observer source.
"""
from __future__ import annotations

import gzip
import hashlib
import importlib
import json
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

from .datasets import Dataset, cache_root
from .tasks import TASKS
from .references import observers

CURATED = "references-curated.jsonl.gz"


def reference_job(job: tuple) -> dict:
    observer, fmt, text, options = job
    tool, task = observer.split(":", 1)
    module = importlib.import_module(observers.MODULES[tool])
    try:
        return {"status": "ok", **module.observe(task, fmt, text, options)}
    except module.Failure as failure:
        return {"status": "error", "error": {"kind": failure.kind, "message": str(failure)}}
    except Exception as error:  # a reference exception is reported per input, never fatal
        return {"status": "error", "error": {"kind": "exception", "message": f"{type(error).__name__}: {error}"}}


def reference_group(group: tuple) -> list[tuple[str, dict]]:
    """Every task for one input text, so the reference reads it once."""
    fmt, text, tasks = group
    return [(key, reference_job((observer, fmt, text, options))) for key, observer, options in tasks]


def compute(groups: list[tuple], workers: int):
    """Observe groups in one pool. When a reference crashes natively and takes
    its worker down, completed groups are kept and the rest are rerun, bisected
    down to the crashing input, whose tasks become `crash` errors."""
    pending = [groups]
    while pending:
        part = pending.pop()
        unfinished = []
        with ProcessPoolExecutor(workers) as pool:
            futures = {pool.submit(reference_group, group): group for group in part}
            for future in as_completed(futures):
                try:
                    yield from future.result()
                except BrokenProcessPool:
                    unfinished.append(futures[future])
        if not unfinished:
            continue
        if len(unfinished) == 1 and len(part) == 1:
            fmt, _, tasks = unfinished[0]
            for key, observer, _ in tasks:
                yield key, {"status": "error", "error": {
                    "kind": "crash", "message": f"{observer} terminated abruptly on this {fmt} input"}}
            continue
        middle = max(1, len(unfinished) // 2)
        pending += [unfinished[:middle], unfinished[middle:]] if len(unfinished) > 1 else [unfinished]


def generate(name: str, workers: int) -> Path:
    from . import run

    started = time.time()
    dataset = Dataset(name)
    directory = dataset.files("full")
    tasks = [t for t in TASKS if t.dataset == name]
    prepared: dict[str, str] = {}

    def prepare(kind: str, path: str, text: str) -> str:
        key = run.prepared_key(kind, path)
        if key not in prepared:
            prepared[key] = observers.prepare(kind, text)
        return prepared[key]

    cases = run.build_cases(dataset, "full", tasks, directory, None, prepare)
    jobs: dict[str, tuple] = {}
    for case in cases:
        for key, job in run.reference_jobs(case):
            jobs[key] = job
    run.progress(started, f"{len(jobs)} reference observations for {len(cases)} cases")
    groups: dict[tuple, list] = defaultdict(list)
    for key, (observer, fmt, text, options) in jobs.items():
        groups[(fmt, text)].append((key, observer, options))
    values = dict(compute([(fmt, text, items) for (fmt, text), items in groups.items()], workers))
    run.progress(started, "references observed")

    header = {"dataset": name, "version": dataset.version, "observers": observers.fingerprint(name),
              "tools": observers.tools(name)}
    curated = {record.id for record in dataset.tier("curated")}
    curated_paths = {file.path for record in dataset.tier("curated") for file in record.files}
    full_path = cache_root() / "references" / f"{name}-v{dataset.version}.jsonl.gz"
    write(full_path, {**header, "inputs": dataset.inputs_digest("full")}, values, prepared)
    write(dataset.directory / CURATED, {**header, "inputs": dataset.inputs_digest("curated")},
          {key: value for key, value in values.items() if key.split("|")[1] in curated_paths},
          {key: text for key, text in prepared.items() if key.split("|", 1)[1] in curated_paths})
    run.progress(started, f"wrote {full_path} and {dataset.directory / CURATED} ({len(curated)} curated records)")
    return full_path


def write(path: Path, header: dict, values: dict[str, dict], prepared: dict[str, str]) -> None:
    """Deterministic gzip JSON lines: header, prepared inputs, then observations."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [json.dumps(header, sort_keys=True)]
    lines += [json.dumps({"prepared": key, "text": prepared[key]}, sort_keys=True) for key in sorted(prepared)]
    lines += [json.dumps({"key": key, "value": values[key]}, sort_keys=True, separators=(",", ":"))
              for key in sorted(values)]
    data = ("\n".join(lines) + "\n").encode("utf-8")
    path.write_bytes(gzip.compress(data, compresslevel=9, mtime=0))
    print(f"{path}: {len(values)} observations, {len(data) / 1e6:.1f} MB uncompressed, "
          f"{path.stat().st_size / 1e6:.1f} MB stored, sha256 {hashlib.sha256(path.read_bytes()).hexdigest()}",
          file=sys.stderr)
