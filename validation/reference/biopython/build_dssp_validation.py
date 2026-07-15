#!/usr/bin/env python3
"""Create DSSP manifests and Biopython/mkdssp goldens for PDB corpora."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import Bio


FEATURE_ID = "bio.secondary-structure.dssp"
REFERENCE_VERSION = "Biopython 1.87 / mkdssp version 4.6.1"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--corpus",
        action="append",
        choices=("smoke", "pdb-10", "pdb-100"),
        help="Corpus to generate. May be repeated; defaults to all three.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
    )
    args = parser.parse_args()
    if Bio.__version__ != "1.87":
        raise SystemExit(f"expected Biopython 1.87, found {Bio.__version__}")
    executable = shutil.which("mkdssp")
    if executable is None:
        raise SystemExit("mkdssp is not available in the active reference environment")
    dssp_version = subprocess.run(
        [executable, "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if dssp_version != "mkdssp version 4.6.1":
        raise SystemExit(f"expected mkdssp version 4.6.1, found {dssp_version}")
    executable_sha256 = sha256_file(Path(executable))

    repo = args.repo_root.resolve()
    script = repo / "validation" / "reference" / "biopython" / "run_feature.py"
    for corpus in args.corpus or ["smoke", "pdb-10", "pdb-100"]:
        root = repo / "validation" / "corpora" / corpus
        fixtures = corpus_fixtures(root, corpus)
        manifest = root / "features" / f"{FEATURE_ID}.toml"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.write_text(
            manifest_text(corpus, fixtures, executable_sha256), encoding="utf-8"
        )
        subprocess.run(
            [
                sys.executable,
                str(script),
                "--feature",
                FEATURE_ID,
                "--corpus",
                corpus,
                "--repo-root",
                str(repo),
            ],
            cwd=repo,
            check=True,
        )
    return 0


def corpus_fixtures(root: Path, corpus: str) -> list[str]:
    if corpus == "smoke":
        return ["data/rcsb/1CRN.cif"]
    lock = json.loads((root / "sources.lock.json").read_text(encoding="utf-8"))
    return [
        file["path"]
        for entry in lock["entries"]
        for file in entry["files"]
        if file["record_type"] == "pdbx-mmcif"
    ]


def manifest_text(corpus: str, fixtures: list[str], executable_sha256: str) -> str:
    fixture_lines = "\n".join(f'  "{fixture}",' for fixture in fixtures)
    return (
        f'feature_id = "{FEATURE_ID}"\n'
        f'corpus_id = "{corpus}"\n'
        'reference_tool = "biopython"\n'
        f'reference_version = "{REFERENCE_VERSION}"\n'
        'comparison_mode = "implementation-golden"\n'
        'notes = [\n'
        f'  "DSSP executable SHA256: {executable_sha256}",\n'
        '  "Command: Bio.PDB.DSSP.DSSP(model, highest_occupancy_snapshot, dssp=mkdssp, file_type=MMCIF)",\n'
        '  "Extended command: mkdssp --output-format=mmcif --quiet highest_occupancy_snapshot annotated.cif",\n'
        ']\n\n'
        f"fixtures = [\n{fixture_lines}\n]\n"
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
