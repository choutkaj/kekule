#!/usr/bin/env python3
"""Build deterministic PubChem validation corpora from official bulk/PUG data."""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

from rdkit import Chem, rdBase

SEED = "molecules-pubchem-v1"
CID_MAX = 200_000_000
HASH_COUNTER_LIMIT = 5_000_000
SHARD_CID_MAX = 500_000
POOL_SIZE = 2_000
CATEGORIES = (
    "disconnected",
    "charged",
    "multi-ring",
    "heteroatom-rich",
    "remaining",
)
TARGET_PER_CATEGORY = 200
SMILES_SNAPSHOT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz"
)
SDF_SHARD_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
    "Compound_000000001_000500000.sdf.gz"
)
SDF_FEATURES = (
    "algo.aromaticity.rdkit-like",
    "algo.rings.fast",
    "algo.rings.sssr",
    "algo.valence.rdkit-like",
    "chem.sanitize.rdkit-like",
    "core.conformers",
    "io.mol.v2000.parse",
    "io.mol.v2000.write",
    "io.mol.v3000.parse",
    "io.mol.v3000.write",
    "io.sdf.v2000.parse",
    "io.sdf.v2000.write",
)
SMILES_FEATURES = ("io.smiles.parse", "io.smiles.write")


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    corpus_root = repo / "validation" / "corpora"
    cache = repo / "target" / "pubchem-cache"
    cache.mkdir(parents=True, exist_ok=True)
    snapshot = cache / "CID-SMILES.gz"
    sdf_shard = cache / "Compound_000000001_000500000.sdf.gz"
    pools_path = cache / "shard-candidate-pools.json"

    if not snapshot.exists():
        download_file(SMILES_SNAPSHOT_URL, snapshot)
    if not sdf_shard.exists():
        download_file(SDF_SHARD_URL, sdf_shard)
    selected = load_or_build_pools(snapshot, sdf_shard, pools_path)

    ordered = [
        selected[category][index]
        for index in range(TARGET_PER_CATEGORY)
        for category in CATEGORIES
    ]
    build_tier(corpus_root / "pubchem-1000", ordered, snapshot, sdf_shard)
    build_tier(corpus_root / "pubchem-100", ordered[:100], snapshot, sdf_shard)
    generate_goldens(repo, "pubchem-100")
    generate_goldens(repo, "pubchem-1000")
    return 0


def load_or_build_pools(
    snapshot: Path, sdf_shard: Path, path: Path
) -> dict[str, list[dict]]:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    ranks: dict[int, int] = {}
    for counter in range(HASH_COUNTER_LIMIT):
        cid = candidate_cid(counter)
        if cid <= SHARD_CID_MAX:
            ranks.setdefault(cid, counter)
    candidates: dict[int, dict] = {}
    with gzip.open(snapshot, "rt", encoding="utf-8") as handle:
        for line in handle:
            cid_text, smiles = line.rstrip("\n").split("\t", 1)
            cid = int(cid_text)
            if cid > SHARD_CID_MAX:
                break
            rank = ranks.get(cid)
            if rank is None:
                continue
            classified = classify_smiles(smiles)
            if classified is None:
                continue
            category, heavy_atoms = classified
            candidates[cid] = {
                "id": str(cid),
                "rank": rank,
                "category": category,
                "smiles": smiles,
                "heavy_atoms": heavy_atoms,
            }
    pools = {category: [] for category in CATEGORIES}
    with gzip.open(sdf_shard, "rt", encoding="utf-8", errors="replace") as handle:
        record_lines: list[str] = []
        for line in handle:
            record_lines.append(line)
            if line.rstrip("\r\n") != "$$$$":
                continue
            record = "".join(record_lines)
            record_lines.clear()
            cid = record_cid(record)
            candidate = candidates.get(cid)
            if candidate is None or " V2000" not in record:
                continue
            mol = Chem.MolFromMolBlock(
                record, sanitize=False, removeHs=False, strictParsing=False
            )
            if mol is None or mol.GetNumAtoms() > 200:
                continue
            pools[candidate["category"]].append(
                {
                    **candidate,
                    "record_type": "2d-current-full-shard",
                    "sdf_url": SDF_SHARD_URL,
                    "sdf_hex": record.encode("utf-8").hex(),
                }
            )
    for category in CATEGORIES:
        pools[category].sort(key=lambda item: item["rank"])
        pools[category] = pools[category][:POOL_SIZE]
        if len(pools[category]) < TARGET_PER_CATEGORY:
            raise SystemExit(
                f"only {len(pools[category])} deterministic candidates for {category}"
            )
    path.write_text(json.dumps(pools, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return pools


def record_cid(record: str) -> int:
    marker = "> <PUBCHEM_COMPOUND_CID>"
    position = record.find(marker)
    if position < 0:
        return -1
    value = record[position + len(marker) :].lstrip("\r\n").splitlines()[0]
    return int(value)


def candidate_cid(counter: int) -> int:
    digest = hashlib.sha256(f"{SEED}:{counter}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % CID_MAX + 1


def classify_smiles(smiles: str) -> tuple[str, int] | None:
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    heavy = sum(atom.GetAtomicNum() > 1 for atom in mol.GetAtoms())
    if not 1 <= heavy <= 150:
        return None
    fragments = len(Chem.GetMolFrags(mol))
    charged = any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms())
    try:
        Chem.GetSymmSSSR(mol)
        ring_count = mol.GetRingInfo().NumRings()
    except Exception:
        return None
    heteroatoms = sum(atom.GetAtomicNum() not in (1, 6) for atom in mol.GetAtoms())
    if fragments > 1:
        return "disconnected", heavy
    if charged:
        return "charged", heavy
    if ring_count >= 2:
        return "multi-ring", heavy
    if heteroatoms >= 3:
        return "heteroatom-rich", heavy
    return "remaining", heavy


def download_file(url: str, path: Path) -> None:
    temporary = path.with_suffix(".part")
    with urllib.request.urlopen(url, timeout=120) as response, temporary.open("wb") as output:
        shutil.copyfileobj(response, output)
    temporary.replace(path)


def build_tier(root: Path, entries: list[dict], snapshot: Path, sdf_shard: Path) -> None:
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
        sdf_path.write_bytes(bytes.fromhex(item["sdf_hex"]))
        smiles_path.write_text(item["smiles"] + "\n", encoding="utf-8")
        lock_entries.append(
            {
                "id": cid,
                "category": item["category"],
                "files": [
                    source_file(root, sdf_path, item["sdf_url"], item["record_type"]),
                    source_file(
                        root,
                        smiles_path,
                        SMILES_SNAPSHOT_URL,
                        "isomeric-smiles-snapshot",
                    ),
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
        sdf_pack.write_bytes(b"".join(bytes.fromhex(item["sdf_hex"]) for item in members))
        smi_pack.write_text(
            "".join(f"{item['smiles']} CID:{item['id']}\n" for item in members),
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
        "source": (
            "PubChem CID-SMILES snapshot "
            f"sha256:{hashlib.sha256(snapshot.read_bytes()).hexdigest()} and CURRENT-Full SDF "
            f"shard sha256:{hashlib.sha256(sdf_shard.read_bytes()).hexdigest()}"
        ),
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
