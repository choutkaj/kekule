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

SEED = "molecules-pdb-v2"
CATEGORIES = (
    "multi-model",
    "protein-nucleic-complex",
    "nucleic-only",
    "protein-heterogen",
    "remaining-protein",
)
BASE_PER_CATEGORY = 20
TARGET_PER_CATEGORY = 200
FEATURE_ID = "io.mmcif.parse"
HOLDINGS_URL = "https://data.rcsb.org/rest/v1/holdings/current/entry_ids"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    corpus_root = repo / "validation" / "corpora"
    base_root = corpus_root / "pdb-100"
    target_root = corpus_root / "pdb-1000"
    base_entries, base_lock_hash = load_base_entries(base_root)
    selected = {
        category: [item for item in base_entries if item["category"] == category]
        for category in CATEGORIES
    }
    if any(len(items) != BASE_PER_CATEGORY for items in selected.values()):
        raise SystemExit("pdb-100 no longer has exactly 20 entries in every category")

    holdings_bytes = fetch(HOLDINGS_URL)
    holdings_ids = set(json.loads(holdings_bytes))
    base_ids = {item["id"] for item in base_entries}
    candidate_pools = {
        category: sorted(
            (
                pdb_id
                for pdb_id in search_entries(query)
                if pdb_id in holdings_ids and pdb_id not in base_ids
            ),
            key=lambda pdb_id: hashlib.sha256(
                f"{SEED}:{pdb_id}".encode()
            ).digest(),
        )
        for category, query in selection_queries().items()
    }
    used_ids = set(base_ids)
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        for category in CATEGORIES:
            examined = 0
            candidates = [
                pdb_id
                for pdb_id in candidate_pools[category]
                if pdb_id not in used_ids
            ]
            for start in range(0, len(candidates), 80):
                batch = candidates[start : start + 80]
                for result in executor.map(fetch_candidate, batch):
                    examined += 1
                    if result is None or result["category"] != category:
                        continue
                    selected[category].append(result)
                    used_ids.add(result["id"])
                print(
                    f"category={category} examined={examined} "
                    f"selected={len(selected[category])}",
                    flush=True,
                )
                if len(selected[category]) >= TARGET_PER_CATEGORY:
                    selected[category] = selected[category][:TARGET_PER_CATEGORY]
                    break
            if len(selected[category]) < TARGET_PER_CATEGORY:
                raise SystemExit(
                    f"RCSB candidate pool exhausted before {category} quota was filled"
                )
    ordered = base_entries + [
        selected[category][index]
        for index in range(BASE_PER_CATEGORY, TARGET_PER_CATEGORY)
        for category in CATEGORIES
    ]
    holdings_hash = hashlib.sha256(holdings_bytes).hexdigest()
    candidate_pool_hash = hashlib.sha256(
        json.dumps(
            candidate_pools, sort_keys=True, separators=(",", ":")
        ).encode()
    ).hexdigest()
    source = (
        f"RCSB PDB holdings sha256:{holdings_hash}; "
        f"RCSB Search API candidate pools sha256:{candidate_pool_hash}; "
        f"pdb-100 lock sha256:{base_lock_hash}"
    )
    build_tier(target_root, ordered, source)
    del ordered, selected, candidate_pools
    write_manifest(base_root, base_entries)
    del base_entries
    generate_golden(repo, "pdb-100")
    generate_golden(repo, "pdb-1000")
    dssp_builder = repo / "validation" / "reference" / "biopython" / "build_dssp_validation.py"
    subprocess.run(
        [
            sys.executable,
            str(dssp_builder),
            "--corpus",
            "pdb-100",
            "--corpus",
            "pdb-1000",
        ],
        cwd=repo,
        check=True,
    )
    mark_ready(target_root)
    return 0


def load_base_entries(root: Path) -> tuple[list[dict], str]:
    lock_path = root / "sources.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    entries = []
    for entry in lock["entries"]:
        files = [
            file for file in entry["files"] if file["record_type"] == "pdbx-mmcif"
        ]
        if len(files) != 1:
            raise SystemExit(f"{entry['id']} does not have exactly one mmCIF file")
        file = files[0]
        path = root / file["path"]
        payload = path.read_bytes() if path.exists() else fetch(file["url"])
        digest = hashlib.sha256(payload).hexdigest()
        if digest != file["sha256"]:
            raise SystemExit(
                f"{entry['id']} SHA256 mismatch: expected {file['sha256']}, found {digest}"
            )
        entries.append(
            {
                "id": entry["id"],
                "category": entry["category"],
                "url": file["url"],
                "payload": payload,
            }
        )
    canonical_lock = json.dumps(lock, sort_keys=True, separators=(",", ":")).encode()
    return entries, hashlib.sha256(canonical_lock).hexdigest()


def selection_queries() -> dict[str, dict]:
    return {
        "multi-model": text_query(
            "rcsb_entry_info.deposited_model_count", "greater", 1
        ),
        "protein-nucleic-complex": text_query(
            "rcsb_entry_info.selected_polymer_entity_types",
            "exact_match",
            "Protein/NA",
        ),
        "nucleic-only": text_query(
            "rcsb_entry_info.selected_polymer_entity_types",
            "exact_match",
            "Nucleic acid (only)",
        ),
        "protein-heterogen": and_query(
            text_query(
                "rcsb_entry_info.selected_polymer_entity_types",
                "exact_match",
                "Protein (only)",
            ),
            text_query(
                "rcsb_entry_info.nonpolymer_entity_count", "greater", 0
            ),
        ),
        "remaining-protein": and_query(
            text_query(
                "rcsb_entry_info.selected_polymer_entity_types",
                "exact_match",
                "Protein (only)",
            ),
            text_query(
                "rcsb_entry_info.nonpolymer_entity_count", "equals", 0
            ),
        ),
    }


def text_query(attribute: str, operator: str, value: object) -> dict:
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": attribute,
            "operator": operator,
            "value": value,
        },
    }


def and_query(*nodes: dict) -> dict:
    return {"type": "group", "logical_operator": "and", "nodes": list(nodes)}


def search_entries(query: dict) -> list[str]:
    payload = {
        "query": query,
        "request_options": {
            "return_all_hits": True,
            "results_verbosity": "compact",
        },
        "return_type": "entry",
    }
    response = post_json(SEARCH_URL, payload)
    results = response.get("result_set")
    if not isinstance(results, list) or not all(
        isinstance(item, str) for item in results
    ):
        raise SystemExit("RCSB Search API returned an unexpected result set")
    return results


def fetch_candidate(pdb_id: str) -> dict | None:
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    temp = Path.cwd() / "target" / "validation-pdb" / f"{pdb_id}.cif"
    if temp.exists():
        payload = temp.read_bytes()
    else:
        try:
            payload = fetch(url, timeout=15, attempts=3)
        except Exception:
            return None
        temp.parent.mkdir(parents=True, exist_ok=True)
        temp.write_bytes(payload)
    if not payload or len(payload) > 2 * 1024 * 1024:
        return None
    try:
        raw = MMCIF2Dict(str(temp))
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


def fetch(url: str, *, timeout: int = 45, attempts: int = 5) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "molecules-validation/1"})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if (
                error.code not in {408, 429, 500, 502, 503, 504}
                or attempt == attempts - 1
            ):
                raise
            time.sleep(0.5 * (2**attempt))
        except (urllib.error.URLError, TimeoutError):
            if attempt == attempts - 1:
                raise
            time.sleep(0.5 * (2**attempt))
    raise RuntimeError("unreachable")


def post_json(url: str, payload: dict) -> dict:
    request = urllib.request.Request(
        url,
        data=json.dumps(payload, separators=(",", ":")).encode(),
        headers={
            "Content-Type": "application/json",
            "User-Agent": "molecules-validation/1",
        },
    )
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                value = json.load(response)
            if not isinstance(value, dict):
                raise SystemExit("RCSB Search API returned a non-object response")
            return value
        except urllib.error.HTTPError as error:
            if error.code not in {408, 429, 500, 502, 503, 504} or attempt == 4:
                raise
            time.sleep(0.5 * (2**attempt))
        except (urllib.error.URLError, TimeoutError):
            if attempt == 4:
                raise
            time.sleep(0.5 * (2**attempt))
    raise RuntimeError("unreachable")


def build_tier(root: Path, entries: list[dict], source: str) -> None:
    data = root / "data"
    if data.exists():
        shutil.rmtree(data)
    data.mkdir(parents=True)
    (root / "features").mkdir(parents=True, exist_ok=True)
    lock_entries = []
    for item in entries:
        path = data / f"{item['id']}.cif"
        path.write_bytes(item["payload"])
        relative = path.relative_to(root).as_posix()
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
        "source": source,
        "selection_seed": SEED,
        "entries": lock_entries,
        "packs": [],
    }
    (root / "sources.lock.json").write_text(
        json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_manifest(root, entries)


def write_manifest(root: Path, entries: list[dict]) -> None:
    fixtures = [f"data/{item['id']}.cif" for item in entries]
    fixture_lines = "\n".join(f'  "{fixture}",' for fixture in fixtures)
    text = (
        f'feature_id = "{FEATURE_ID}"\n'
        f'corpus_id = "{root.name}"\n'
        'reference_tool = "biopython"\n'
        f'reference_version = "Biopython {Bio.__version__}"\n'
        'comparison_mode = "implementation-golden"\n\n'
        f"fixtures = [\n{fixture_lines}\n]\n"
    )
    feature_dir = root / "features"
    feature_dir.mkdir(parents=True, exist_ok=True)
    (feature_dir / f"{FEATURE_ID}.toml").write_text(text, encoding="utf-8")


def generate_golden(repo: Path, corpus: str) -> None:
    output_dir = repo / "validation" / "corpora" / corpus / "golden" / FEATURE_ID
    if output_dir.exists():
        shutil.rmtree(output_dir)
    script = repo / "validation" / "reference" / "biopython" / "run_feature.py"
    subprocess.run(
        [sys.executable, str(script), "--feature", FEATURE_ID, "--corpus", corpus],
        cwd=repo,
        check=True,
    )


def mark_ready(root: Path) -> None:
    descriptor = root / "corpus.toml"
    descriptor.write_text(
        descriptor.read_text(encoding="utf-8").replace("ready = false", "ready = true"),
        encoding="utf-8",
    )

if __name__ == "__main__":
    raise SystemExit(main())
