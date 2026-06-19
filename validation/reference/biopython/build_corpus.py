#!/usr/bin/env python3
"""Build the deterministic PDB validation corpora."""

from __future__ import annotations

import concurrent.futures
import hashlib
import json
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import Bio
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser

SEED = "molecules-pdb-v1"
CATEGORIES = (
    "multi-model",
    "protein-nucleic-complex",
    "nucleic-only",
    "protein-heterogen",
    "remaining-protein",
)
TARGET_PER_CATEGORY = 20
FEATURES = ("bio.hierarchy.smcra", "io.mmcif.parse")
HOLDINGS_URL = "https://data.rcsb.org/rest/v1/holdings/current/entry_ids"


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    corpus_root = repo / "validation" / "corpora"
    holdings_bytes = fetch(HOLDINGS_URL)
    ids = json.loads(holdings_bytes)
    ranked = sorted(ids, key=lambda pdb_id: hashlib.sha256(f"{SEED}:{pdb_id}".encode()).digest())
    selected = {category: [] for category in CATEGORIES}
    examined = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        for start in range(0, len(ranked), 80):
            batch = ranked[start : start + 80]
            for result in executor.map(fetch_candidate, batch):
                examined += 1
                if result is None:
                    continue
                category = result["category"]
                if len(selected[category]) < TARGET_PER_CATEGORY:
                    selected[category].append(result)
            counts = " ".join(f"{name}={len(selected[name])}" for name in CATEGORIES)
            print(f"examined={examined} {counts}", flush=True)
            if min(len(values) for values in selected.values()) >= TARGET_PER_CATEGORY:
                break
    if min(len(values) for values in selected.values()) < TARGET_PER_CATEGORY:
        raise SystemExit("RCSB holdings exhausted before category quotas were filled")
    ordered = [
        selected[category][index]
        for index in range(TARGET_PER_CATEGORY)
        for category in CATEGORIES
    ]
    holdings_hash = hashlib.sha256(holdings_bytes).hexdigest()
    build_tier(corpus_root / "pdb-100", ordered, holdings_hash)
    build_tier(corpus_root / "pdb-10", ordered[:10], holdings_hash)
    generate_goldens(repo, "pdb-10")
    generate_goldens(repo, "pdb-100")
    return 0


def fetch_candidate(pdb_id: str) -> dict | None:
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        payload = fetch(url)
    except Exception:
        return None
    if not payload or len(payload) > 2 * 1024 * 1024:
        return None
    temp = Path.cwd() / "target" / "validation-pdb" / f"{pdb_id}.cif"
    temp.parent.mkdir(parents=True, exist_ok=True)
    temp.write_bytes(payload)
    try:
        raw = MMCIF2Dict(str(temp))
        MMCIFParser(QUIET=True).get_structure(pdb_id, str(temp))
    except Exception:
        return None
    atom_ids = column(raw, "_atom_site.id")
    if not 1 <= len(atom_ids) <= 20_000:
        return None
    models = set(column(raw, "_atom_site.pdbx_PDB_model_num"))
    polymer_types = " ".join(column(raw, "_entity_poly.type")).lower()
    protein = "polypeptide" in polymer_types
    nucleic = "polyribonucleotide" in polymer_types or "polydeoxyribonucleotide" in polymer_types
    groups = column(raw, "_atom_site.group_PDB")
    components = column(raw, "_atom_site.label_comp_id")
    heterogen = any(
        group == "HETATM" and component.upper() not in {"HOH", "WAT", "DOD"}
        for group, component in zip(groups, components)
    )
    if len(models) > 1:
        category = "multi-model"
    elif protein and nucleic:
        category = "protein-nucleic-complex"
    elif nucleic and not protein:
        category = "nucleic-only"
    elif protein and heterogen:
        category = "protein-heterogen"
    elif protein:
        category = "remaining-protein"
    else:
        return None
    return {"id": pdb_id, "category": category, "url": url, "payload": payload}


def column(raw: dict, name: str) -> list[str]:
    value = raw.get(name, [])
    return [str(item) for item in value] if isinstance(value, list) else [str(value)]


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "molecules-validation/1"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                return response.read()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            if attempt == 4:
                raise
            time.sleep(0.5 * (2**attempt))
    raise RuntimeError("unreachable")


def build_tier(root: Path, entries: list[dict], holdings_hash: str) -> None:
    data = root / "data"
    if data.exists():
        shutil.rmtree(data)
    data.mkdir(parents=True)
    (root / "features").mkdir(parents=True, exist_ok=True)
    lock_entries = []
    fixtures = []
    for item in entries:
        path = data / f"{item['id']}.cif"
        path.write_bytes(item["payload"])
        relative = path.relative_to(root).as_posix()
        fixtures.append(relative)
        lock_entries.append(
            {
                "id": item["id"],
                "category": item["category"],
                "files": [
                    {
                        "path": relative,
                        "url": item["url"],
                        "sha256": hashlib.sha256(item["payload"]).hexdigest(),
                        "record_type": "pdbx-mmcif",
                    }
                ],
            }
        )
    lock = {
        "schema_version": 1,
        "corpus_id": root.name,
        "source": f"RCSB PDB holdings sha256:{holdings_hash}",
        "selection_seed": SEED,
        "entries": lock_entries,
        "packs": [],
    }
    (root / "sources.lock.json").write_text(
        json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    array = ",\n".join(f'  "{fixture}"' for fixture in fixtures)
    for feature in FEATURES:
        text = (
            f'feature_id = "{feature}"\n'
            f'corpus_id = "{root.name}"\n'
            'reference_tool = "biopython"\n'
            f'reference_version = "Biopython {Bio.__version__}"\n'
            'comparison_mode = "implementation-golden"\n\n'
            f"fixtures = [\n{array},\n]\n"
        )
        (root / "features" / f"{feature}.toml").write_text(text, encoding="utf-8")
    descriptor = root / "corpus.toml"
    descriptor.write_text(
        descriptor.read_text(encoding="utf-8").replace("ready = false", "ready = true"),
        encoding="utf-8",
    )


def generate_goldens(repo: Path, corpus: str) -> None:
    script = repo / "validation" / "reference" / "biopython" / "run_feature.py"
    for feature in FEATURES:
        subprocess.run(
            [sys.executable, str(script), "--feature", feature, "--corpus", corpus],
            cwd=repo,
            check=True,
        )


if __name__ == "__main__":
    raise SystemExit(main())
