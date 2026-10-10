#!/usr/bin/env python3
"""Kekule scientific benchmark.

Comparison (standard library and the Kekule observer only):
  run DATASET --tier T     observe with Kekule, compare with stored references, triage
Datasets:
  fetch [DATASET]          download and verify the published snapshot archives
  verify DATASET --tier T  check every file of a tier against the manifest
Maintainers only (in the environment from environment.yml):
  select DATASET           resolve selection.toml against upstream databases
  reference DATASET        compute and store every reference observation once
  pack DATASET             build the snapshot archive (inputs and references)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from kbench import datasets, selection  # noqa: E402

RELEASES = "https://github.com/choutkaj/kekule/releases/download"


def dataset_names() -> list[str]:
    return sorted(path.parent.name for path in datasets.DATASETS.glob("*/selection.toml"))


def command_select(arguments) -> None:
    records = selection.run(arguments.dataset)
    counts: dict[str, int] = {}
    for record in records:
        counts[record.stratum] = counts.get(record.stratum, 0) + 1
    print(json.dumps(counts, indent=2))


def command_pack(arguments) -> None:
    dataset = datasets.Dataset(arguments.dataset)
    tag = f"benchmark-data-{dataset.name}-v{dataset.version}"
    name = f"kekule-bench-{dataset.name}-v{dataset.version}.tar.xz"
    staging = datasets.cache_root() / "staging" / dataset.name / f"v{dataset.version}"
    archive = datasets.cache_root() / "archives" / name
    extras = {
        f"NOTICES/{notice}": (dataset.directory / notice).read_bytes()
        for notice in dataset.selection["dataset"]["notices"]
    }
    references = datasets.cache_root() / "references" / f"{dataset.name}-v{dataset.version}.jsonl.gz"
    if not references.exists():
        raise datasets.DatasetError(f"generate the references first: bench.py reference {dataset.name}")
    extras["references.jsonl.gz"] = references.read_bytes()
    digest = datasets.pack(staging, dataset.records.values(), archive, extras)
    lock_path = dataset.directory / "selection.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    lock["snapshot"] = {
        "name": name,
        "url": f"{RELEASES}/{tag}/{name}",
        "sha256": digest,
        "bytes": archive.stat().st_size,
        "extras": {member: datasets.sha256_bytes(data) for member, data in extras.items()},
    }
    lock_path.write_text(json.dumps(lock, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(f"{archive}\nsha256 {digest}\nupload as release asset of tag {tag}")


def command_reference(arguments) -> None:
    from kbench import generate

    generate.generate(arguments.dataset, arguments.jobs)


def command_run(arguments) -> int:
    from kbench import run

    return run.run(arguments.dataset, arguments.tier, arguments.task, arguments.jobs, arguments.limit,
                   build=not arguments.no_build, update_known=arguments.update_known)


def command_fetch(arguments) -> None:
    for name in [arguments.dataset] if arguments.dataset else dataset_names():
        directory = datasets.Dataset(name).files("full")
        print(f"{name}: {directory}")


def command_verify(arguments) -> None:
    dataset = datasets.Dataset(arguments.dataset)
    directory = dataset.files(arguments.tier)
    print(f"{dataset.name} {arguments.tier}: {len(dataset.tier(arguments.tier))} records verified in {directory}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="observe, compare and triage one dataset tier")
    run.add_argument("dataset", choices=dataset_names())
    run.add_argument("--tier", choices=["curated", "full"], default="curated")
    run.add_argument("--task", action="append", help="task id or glob; repeatable (default: all)")
    run.add_argument("--jobs", type=int, default=os.cpu_count() or 1, help="observer and reference workers")
    run.add_argument("--limit", type=int, help="first N records only, for development")
    run.add_argument("--no-build", action="store_true", help="use the existing release observer")
    run.add_argument("--update-known", action="store_true", help="lower known-difference bounds to observed counts")
    run.set_defaults(handler=command_run)
    fetch = commands.add_parser("fetch", help="download and verify snapshot archives")
    fetch.add_argument("dataset", nargs="?", choices=dataset_names())
    fetch.set_defaults(handler=command_fetch)
    verify = commands.add_parser("verify", help="verify a tier against its manifest")
    verify.add_argument("dataset", choices=dataset_names())
    verify.add_argument("--tier", choices=["curated", "full"], default="curated")
    verify.set_defaults(handler=command_verify)
    select = commands.add_parser("select", help="maintainers: resolve the selection upstream")
    select.add_argument("dataset", choices=dataset_names())
    select.set_defaults(handler=command_select)
    reference = commands.add_parser("reference", help="maintainers: store the reference observations")
    reference.add_argument("dataset", choices=dataset_names())
    reference.add_argument("--jobs", type=int, default=os.cpu_count() or 1, help="reference workers")
    reference.set_defaults(handler=command_reference)
    pack = commands.add_parser("pack", help="maintainers: build the snapshot archive")
    pack.add_argument("dataset", choices=dataset_names())
    pack.set_defaults(handler=command_pack)
    arguments = parser.parse_args(argv)
    try:
        return arguments.handler(arguments) or 0
    except datasets.DatasetError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
