#!/usr/bin/env python3
"""Build the deterministic PubChem validation corpora."""

from __future__ import annotations

import concurrent.futures
import base64
import hashlib
import json
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

from rdkit import Chem, rdBase

SEED = "molecules-pubchem-v1"
CID_MAX = 200_000_000
CATEGORIES = (
    "disconnected",
    "charged",
    "multi-ring",
    "heteroatom-rich",
    "remaining",
)
TARGET_PER_CATEGORY = 200
SDF_FEATURES = (
    "algo.aromaticity.rdkit-like",
    "algo.rings.fast",
    "algo.rings.sssr",
    "algo.valence.rdkit-like",
    "chem.sanitize.rdkit-like",
    "core.conformers",
    "io.mol.v2000.parse",
    "io.mol.v2000.write",
    "io.sdf.v2000.parse",
    "io.sdf.v2000.write",
)
SMILES_FEATURES = ("io.smiles.parse", "io.smiles.write")


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    corpus_root = repo / "validation" / "corpora"
    checkpoint_path = repo / "target" / "pubchem-corpus-checkpoint.json"
    selected, seen, counter = load_checkpoint(checkpoint_path)

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        while min(len(values) for values in selected.values()) < TARGET_PER_CATEGORY:
            candidates = []
            while len(candidates) < 40:
                cid = candidate_cid(counter)
                counter += 1
                if cid not in seen:
                    seen.add(cid)
                    candidates.append(cid)
            for result in executor.map(fetch_candidate, candidates):
                if result is None:
                    continue
                category = result["category"]
                if len(selected[category]) < TARGET_PER_CATEGORY:
                    selected[category].append(result)
                    save_checkpoint(checkpoint_path, selected, seen, counter)
            counts = " ".join(f"{name}={len(selected[name])}" for name in CATEGORIES)
            print(f"examined={counter} {counts}", flush=True)

    ordered = [
        selected[category][index]
        for index in range(TARGET_PER_CATEGORY)
        for category in CATEGORIES
    ]
    build_tier(corpus_root / "pubchem-1000", ordered)
    build_tier(corpus_root / "pubchem-100", ordered[:100])
    generate_goldens(repo, "pubchem-100")
    generate_goldens(repo, "pubchem-1000")
    checkpoint_path.unlink(missing_ok=True)
    return 0


def load_checkpoint(path: Path) -> tuple[dict[str, list[dict]], set[int], int]:
    if not path.exists():
        return {category: [] for category in CATEGORIES}, set(), 0
    raw = json.loads(path.read_text(encoding="utf-8"))
    selected = {category: [] for category in CATEGORIES}
    for category, entries in raw["selected"].items():
        for entry in entries:
            entry["sdf"] = base64.b64decode(entry["sdf"])
            entry["smiles"] = base64.b64decode(entry["smiles"])
            selected[category].append(entry)
    print(
        "resuming checkpoint "
        + " ".join(f"{name}={len(selected[name])}" for name in CATEGORIES),
        flush=True,
    )
    return selected, set(raw["seen"]), raw["counter"]


def save_checkpoint(
    path: Path, selected: dict[str, list[dict]], seen: set[int], counter: int
) -> None:
    serializable = {}
    for category, entries in selected.items():
        serializable[category] = []
        for entry in entries:
            encoded = dict(entry)
            encoded["sdf"] = base64.b64encode(entry["sdf"]).decode("ascii")
            encoded["smiles"] = base64.b64encode(entry["smiles"]).decode("ascii")
            serializable[category].append(encoded)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(
        json.dumps(
            {"counter": counter, "seen": sorted(seen), "selected": serializable},
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    temporary.replace(path)


def candidate_cid(counter: int) -> int:
    digest = hashlib.sha256(f"{SEED}:{counter}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % CID_MAX + 1


def fetch_candidate(cid: int) -> dict | None:
    smiles_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/"
        "property/IsomericSMILES/TXT"
    )
    try:
        smiles_bytes = fetch(smiles_url)
    except Exception:
        return None
    smiles = smiles_bytes.decode("utf-8").strip()
    record_type = "2d" if "." in smiles else "3d"
    sdf_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
        f"?record_type={record_type}"
    )
    try:
        sdf_bytes = fetch(sdf_url)
    except Exception:
        return None
    text = sdf_bytes.decode("utf-8", "replace")
    if text.count("$$$$") != 1 or " V2000" not in text:
        return None
    mol = Chem.MolFromMolBlock(text, sanitize=False, removeHs=False, strictParsing=False)
    if mol is None:
        return None
    heavy = sum(atom.GetAtomicNum() > 1 for atom in mol.GetAtoms())
    total = mol.GetNumAtoms()
    if not 1 <= heavy <= 150 or total > 200:
        return None
    Chem.GetSymmSSSR(mol)
    fragments = len(Chem.GetMolFrags(mol))
    charged = any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms())
    ring_count = mol.GetRingInfo().NumRings()
    heteroatoms = sum(atom.GetAtomicNum() not in (1, 6) for atom in mol.GetAtoms())
    if fragments > 1:
        category = "disconnected"
    elif charged:
        category = "charged"
    elif ring_count >= 2:
        category = "multi-ring"
    elif heteroatoms >= 3:
        category = "heteroatom-rich"
    else:
        category = "remaining"
    return {
        "id": str(cid),
        "category": category,
        "record_type": record_type,
        "sdf_url": sdf_url,
        "smiles_url": smiles_url,
        "sdf": sdf_bytes,
        "smiles": smiles_bytes,
    }


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "molecules-validation/1"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                payload = response.read()
                time.sleep(0.08)
                return payload
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            if attempt == 4:
                raise
            time.sleep(0.5 * (2**attempt))
    raise RuntimeError("unreachable")


def build_tier(root: Path, entries: list[dict]) -> None:
    data = root / "data"
    if data.exists():
        shutil.rmtree(data)
    (data / "raw").mkdir(parents=True)
    (data / "packs").mkdir(parents=True)
    (root / "features").mkdir(parents=True, exist_ok=True)
    lock_entries = []
    for item in entries:
        cid = item["id"]
        sdf_path = data / "raw" / f"cid_{cid}.sdf"
        smiles_path = data / "raw" / f"cid_{cid}.smi"
        sdf_path.write_bytes(item["sdf"])
        smiles_path.write_bytes(item["smiles"])
        lock_entries.append(
            {
                "id": cid,
                "category": item["category"],
                "files": [
                    source_file(root, sdf_path, item["sdf_url"], item["record_type"]),
                    source_file(root, smiles_path, item["smiles_url"], "isomeric-smiles"),
                ],
            }
        )
    packs = []
    sdf_fixtures = []
    smiles_fixtures = []
    for pack_index in range(0, len(entries), 100):
        members = entries[pack_index : pack_index + 100]
        number = pack_index // 100 + 1
        sdf_pack = data / "packs" / f"pack_{number:02}.sdf"
        smi_pack = data / "packs" / f"pack_{number:02}.smi"
        sdf_pack.write_bytes(b"".join(item["sdf"] for item in members))
        smi_pack.write_text(
            "".join(f"{item['smiles'].decode().strip()} CID:{item['id']}\n" for item in members),
            encoding="utf-8",
        )
        ids = [item["id"] for item in members]
        packs.extend(
            [
                pack_record(root, sdf_pack, "sdf-v2000", ids),
                pack_record(root, smi_pack, "smiles", ids),
            ]
        )
        sdf_fixtures.append(sdf_pack.relative_to(root).as_posix())
        smiles_fixtures.append(smi_pack.relative_to(root).as_posix())
    lock = {
        "schema_version": 1,
        "corpus_id": root.name,
        "source": "PubChem PUG REST",
        "selection_seed": SEED,
        "entries": lock_entries,
        "packs": packs,
    }
    write_json(root / "sources.lock.json", lock)
    write_manifests(root, SDF_FEATURES, sdf_fixtures)
    write_manifests(root, SMILES_FEATURES, smiles_fixtures)
    descriptor = root / "corpus.toml"
    descriptor.write_text(
        descriptor.read_text(encoding="utf-8").replace("ready = false", "ready = true"),
        encoding="utf-8",
    )


def source_file(root: Path, path: Path, url: str, record_type: str) -> dict:
    return {
        "path": path.relative_to(root).as_posix(),
        "url": url,
        "sha256": sha256(path),
        "record_type": record_type,
    }


def pack_record(root: Path, path: Path, format_name: str, members: list[str]) -> dict:
    return {
        "path": path.relative_to(root).as_posix(),
        "format": format_name,
        "count": len(members),
        "members": members,
        "sha256": sha256(path),
    }


def write_manifests(root: Path, features: tuple[str, ...], fixtures: list[str]) -> None:
    array = ",\n".join(f'  "{fixture}"' for fixture in fixtures)
    for feature in features:
        text = (
            f'feature_id = "{feature}"\n'
            f'corpus_id = "{root.name}"\n'
            'reference_tool = "rdkit"\n'
            f'reference_version = "RDKit {rdBase.rdkitVersion}"\n'
            'comparison_mode = "implementation-golden"\n\n'
            f"fixtures = [\n{array},\n]\n"
        )
        (root / "features" / f"{feature}.toml").write_text(text, encoding="utf-8")


def generate_goldens(repo: Path, corpus: str) -> None:
    script = repo / "validation" / "reference" / "rdkit" / "run_feature.py"
    for feature in (*SDF_FEATURES, *SMILES_FEATURES):
        subprocess.run(
            [sys.executable, str(script), "--feature", feature, "--corpus", corpus],
            cwd=repo,
            check=True,
        )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
