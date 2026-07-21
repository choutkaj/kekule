#!/usr/bin/env python3
"""Build large RDKit-backed validation corpora from local source data."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

RDLogger.DisableLog("rdApp.*")

PUBCHEM_SEED = "molecular-pubchem-100k-v1"
PUBCHEM_TARGET = 100_000
PUBCHEM_CID_MAX = 500_000
PUBCHEM_PACK_SIZE = 1_000
SMILES_SNAPSHOT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz"
)
SDF_SHARD_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
    "Compound_000000001_000500000.sdf.gz"
)

PL_REX_SEED = "pl-rex-primary-v1"
PL_REX_SOURCE_URL = "https://doi.org/10.26434/chemrxiv-2023-zh03k"
ENAMINE_SEED = "enamine-diversity-20260524"
ENAMINE_SOURCE_URL = "https://enamine.net/compound-libraries/diversity-libraries"

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
SMILES_FEATURES = ("io.smiles.parse", "io.smiles.write", "io.smiles.canonical")
ENAMINE_SMILES_FEATURES = ("io.smiles.parse", "io.smiles.write")
PUBCHEM_100K_SDF_FEATURES = (*SDF_FEATURES, "algo.canonical-ranking")
PUBCHEM_100K_SMILES_FEATURES = (*SMILES_FEATURES, "stereo.cip")
PUBCHEM_MANUAL_SMILES_FEATURES = ("stereo.perception", "stereo.representation")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--corpus",
        choices=("all", "pl-rex", "enamine-diversity", "pubchem-100k"),
        default="all",
    )
    parser.add_argument("--skip-goldens", action="store_true")
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
        help="Repository root. Defaults to the script's containing checkout.",
    )
    args = parser.parse_args()

    repo = args.repo_root.resolve()
    selected = (
        ("pl-rex", "enamine-diversity", "pubchem-100k")
        if args.corpus == "all"
        else (args.corpus,)
    )
    for corpus in selected:
        if corpus == "pl-rex":
            build_pl_rex(repo)
            if not args.skip_goldens:
                generate_goldens(repo, corpus, SDF_FEATURES)
        elif corpus == "enamine-diversity":
            build_enamine(repo)
            if not args.skip_goldens:
                generate_goldens(repo, corpus, (*SDF_FEATURES, *ENAMINE_SMILES_FEATURES))
        elif corpus == "pubchem-100k":
            build_pubchem_100k(repo)
            if not args.skip_goldens:
                generate_goldens(
                    repo,
                    corpus,
                    (*PUBCHEM_100K_SDF_FEATURES, *PUBCHEM_100K_SMILES_FEATURES),
                )
                generate_manual_goldens(repo, corpus)
    return 0


def build_pl_rex(repo: Path) -> None:
    root = repo / "validation" / "corpora" / "pl-rex"
    source_files = sorted(root.glob("data/*/structures/*/ligand.sdf"))
    if len(source_files) != 164:
        raise SystemExit(f"expected 164 primary PL-REX ligands, found {len(source_files)}")
    pack_dir = reset_dir(root / "data" / "packs")
    fixtures: list[str] = []
    packs: list[dict] = []
    entries: list[dict] = []
    for pack_index, chunk in enumerate(chunks(source_files, 100), start=1):
        pack_path = pack_dir / f"pack_{pack_index:02}.sdf"
        members = []
        with pack_path.open("w", encoding="utf-8", newline="\n") as output:
            for path in chunk:
                target = path.parts[-4]
                pdb = path.parts[-2]
                member_id = f"{target}:{pdb}"
                members.append(member_id)
                output.write(add_sdf_property(read_single_sdf_record(path), "PL_REX_ID", member_id))
                entries.append(
                    {
                        "id": member_id,
                        "category": "primary-ligand",
                        "files": [
                            {
                                "path": path.relative_to(root).as_posix(),
                                "url": PL_REX_SOURCE_URL,
                                "sha256": sha256(path),
                                "record_type": "pl-rex-primary-ligand-sdf",
                            }
                        ],
                    }
                )
        fixtures.append(pack_path.relative_to(root).as_posix())
        packs.append(
            pack_record(
                root,
                pack_path,
                "sdf-v2000",
                members,
                member_id_property="PL_REX_ID",
            )
        )
    write_json(
        root / "sources.lock.json",
        {
            "schema_version": 1,
            "corpus_id": "pl-rex",
            "source": "PL-REX 1.0.1 primary refined ligand SDF files",
            "selection_seed": PL_REX_SEED,
            "entries": entries,
            "packs": packs,
        },
    )
    reset_dir(root / "features")
    reset_dir(root / "golden")
    write_manifests(root, SDF_FEATURES, fixtures)


def build_enamine(repo: Path) -> None:
    root = repo / "validation" / "corpora" / "enamine-diversity"
    data = root / "data"
    sdf_source = one(data.glob("Enamine_Discovery_Diversity_Set_50_plated_*cmpds_*.sdf"))
    smiles_source = one(data.glob("Enamine_Discovery_Diversity_Set_50_plated_*cmpds_*.smiles"))
    csv_source = one(data.glob("Enamine_Discovery_Diversity_Set_50_plated_*cmpds_*.csv"))
    records = read_sdf_records(sdf_source)
    smiles_rows = read_enamine_smiles(smiles_source)
    if len(records) != 50_240 or len(smiles_rows) != 50_240:
        raise SystemExit(
            f"expected 50240 Enamine records, found {len(records)} SDF and {len(smiles_rows)} SMILES"
        )

    pack_dir = reset_dir(data / "packs")
    sdf_rows: list[tuple[str, str]] = []
    sdf_fixtures: list[str] = []
    smiles_fixtures: list[str] = []
    packs: list[dict] = []
    entries: list[dict] = []
    for index, (record, (smiles, catalog_id)) in enumerate(zip(records, smiles_rows), start=1):
        record = normalize_sdf_record_header(record)
        sdf_id = sdf_property(record, "Catalog ID")
        if sdf_id != catalog_id:
            raise SystemExit(f"Enamine SDF/SMILES order differs at record {index}: {sdf_id} != {catalog_id}")
        entries.append({"id": catalog_id, "category": "diversity", "files": []})
        if sdf_counts_line(record).endswith("V2000"):
            sdf_rows.append((record, catalog_id))

    for pack_index, offset in enumerate(range(0, len(sdf_rows), 1000), start=1):
        sdf_chunk = sdf_rows[offset : offset + 1000]
        ids = [catalog_id for _, catalog_id in sdf_chunk]
        sdf_pack = pack_dir / f"pack_{pack_index:03}.sdf"
        sdf_pack.write_text(
            "".join(record for record, _ in sdf_chunk),
            encoding="utf-8",
            newline="\n",
        )
        sdf_fixtures.append(sdf_pack.relative_to(root).as_posix())
        packs.append(
            pack_record(
                root,
                sdf_pack,
                "sdf-v2000",
                ids,
                member_id_property="Catalog ID",
            )
        )

    for pack_index, offset in enumerate(range(0, len(smiles_rows), 1000), start=1):
        smiles_chunk = smiles_rows[offset : offset + 1000]
        ids = [catalog_id for _, catalog_id in smiles_chunk]
        smiles_pack = pack_dir / f"pack_{pack_index:03}.smi"
        smiles_pack.write_text(
            "".join(f"{smiles} ID:{catalog_id}\n" for smiles, catalog_id in smiles_chunk),
            encoding="utf-8",
            newline="\n",
        )
        smiles_fixtures.append(smiles_pack.relative_to(root).as_posix())
        packs.append(
            pack_record(
                root,
                smiles_pack,
                "smiles",
                ids,
                member_title_prefix="ID:",
            )
        )

    write_json(
        root / "sources.lock.json",
        {
            "schema_version": 1,
            "corpus_id": "enamine-diversity",
            "source": (
                "Enamine Discovery Diversity Set 50 plated 50240 compound files "
                f"sdf sha256:{sha256(sdf_source)}, smiles sha256:{sha256(smiles_source)}, "
                f"csv sha256:{sha256(csv_source)}; SDF validation packs include "
                f"{len(sdf_rows)} V2000 records and SMILES packs include all 50240 records"
            ),
            "selection_seed": ENAMINE_SEED,
            "entries": entries,
            "packs": packs,
        },
    )
    reset_dir(root / "features")
    reset_dir(root / "golden")
    write_manifests(root, SDF_FEATURES, sdf_fixtures)
    write_manifests(root, ENAMINE_SMILES_FEATURES, smiles_fixtures)


def build_pubchem_100k(repo: Path) -> None:
    root = repo / "validation" / "corpora" / "pubchem-100k"
    cache = repo / "target" / "pubchem-cache"
    snapshot = cache / "CID-SMILES.gz"
    sdf_shard = cache / "Compound_000000001_000500000.sdf.gz"
    if not snapshot.exists() or not sdf_shard.exists():
        raise SystemExit(
            "PubChem cache is missing. Run validation/reference/rdkit/build_corpus.py first."
        )
    data = root / "data"
    pack_dir = reset_dir(data / "packs")
    smiles_by_cid = load_pubchem_smiles(snapshot)
    candidates = select_pubchem_candidates(sdf_shard, smiles_by_cid)
    selected_ids = {item["id"] for item in candidates}
    entries = [
        {"id": item["id"], "category": item["category"], "files": []}
        for item in candidates
    ]
    sdf_packs: list[dict] = []
    smiles_packs: list[dict] = []
    sdf_fixtures: list[str] = []
    smiles_fixtures: list[str] = []

    pack_members: list[dict] = []
    with gzip.open(sdf_shard, "rt", encoding="utf-8", errors="replace") as handle:
        record_lines: list[str] = []
        for line in handle:
            record_lines.append(line)
            if line.rstrip("\r\n") != "$$$$":
                continue
            record = "".join(record_lines)
            record_lines.clear()
            cid = record_cid(record)
            if cid in selected_ids:
                pack_members.append(
                    {
                        "id": cid,
                        "sdf": record,
                        "smiles": smiles_by_cid[cid],
                    }
                )
                if len(pack_members) == PUBCHEM_PACK_SIZE:
                    write_pubchem_pack(root, pack_dir, pack_members, sdf_fixtures, smiles_fixtures, sdf_packs, smiles_packs)
                    pack_members = []
        if pack_members:
            write_pubchem_pack(root, pack_dir, pack_members, sdf_fixtures, smiles_fixtures, sdf_packs, smiles_packs)

    ordered_ids = [item["id"] for item in candidates]
    packed_ids = [member for pack in sdf_packs for member in pack["members"]]
    if packed_ids != ordered_ids:
        raise SystemExit("PubChem pack order differs from selected lock order")

    write_json(
        root / "sources.lock.json",
        {
            "schema_version": 1,
            "corpus_id": "pubchem-100k",
            "source": (
                "PubChem CID-SMILES snapshot "
                f"sha256:{sha256(snapshot)} and CURRENT-Full SDF shard "
                f"sha256:{sha256(sdf_shard)}"
            ),
            "selection_seed": PUBCHEM_SEED,
            "entries": entries,
            "packs": [*sdf_packs, *smiles_packs],
        },
    )
    reset_dir(root / "features")
    reset_dir(root / "golden")
    write_manifests(root, PUBCHEM_100K_SDF_FEATURES, sdf_fixtures)
    write_manifests(root, PUBCHEM_100K_SMILES_FEATURES, smiles_fixtures)
    write_manual_manifests(root, smiles_fixtures)


def load_pubchem_smiles(snapshot: Path) -> dict[str, str]:
    smiles_by_cid = {}
    with gzip.open(snapshot, "rt", encoding="utf-8") as handle:
        for line in handle:
            cid_text, smiles = line.rstrip("\n").split("\t", 1)
            cid = int(cid_text)
            if cid > PUBCHEM_CID_MAX:
                break
            smiles_by_cid[cid_text] = smiles
    return smiles_by_cid


def select_pubchem_candidates(sdf_shard: Path, smiles_by_cid: dict[str, str]) -> list[dict]:
    candidates: list[dict] = []
    with gzip.open(sdf_shard, "rt", encoding="utf-8", errors="replace") as handle:
        record_lines: list[str] = []
        for line in handle:
            record_lines.append(line)
            if line.rstrip("\r\n") != "$$$$":
                continue
            record = "".join(record_lines)
            record_lines.clear()
            cid = record_cid(record)
            smiles = smiles_by_cid.get(cid)
            if smiles is None or " V2000" not in record:
                continue
            mol = Chem.MolFromMolBlock(record, sanitize=False, removeHs=False, strictParsing=False)
            if mol is None or mol.GetNumAtoms() > 200:
                continue
            classified = classify_smiles(smiles)
            if classified is None:
                continue
            candidates.append(
                {
                    "id": cid,
                    "rank": pubchem_rank(cid),
                    "category": classified,
                }
            )
    candidates.sort(key=lambda item: item["rank"])
    if len(candidates) < PUBCHEM_TARGET:
        raise SystemExit(f"only {len(candidates)} PubChem candidates available")
    return sorted(candidates[:PUBCHEM_TARGET], key=lambda item: int(item["id"]))


def write_pubchem_pack(
    root: Path,
    pack_dir: Path,
    members: list[dict],
    sdf_fixtures: list[str],
    smiles_fixtures: list[str],
    sdf_packs: list[dict],
    smiles_packs: list[dict],
) -> None:
    number = len(sdf_packs) + 1
    sdf_pack = pack_dir / f"pack_{number:03}.sdf"
    smiles_pack = pack_dir / f"pack_{number:03}.smi"
    ids = [member["id"] for member in members]
    sdf_pack.write_text("".join(member["sdf"] for member in members), encoding="utf-8", newline="\n")
    smiles_pack.write_text(
        "".join(f"{member['smiles']} CID:{member['id']}\n" for member in members),
        encoding="utf-8",
        newline="\n",
    )
    sdf_fixtures.append(sdf_pack.relative_to(root).as_posix())
    smiles_fixtures.append(smiles_pack.relative_to(root).as_posix())
    sdf_packs.append(pack_record(root, sdf_pack, "sdf-v2000", ids))
    smiles_packs.append(pack_record(root, smiles_pack, "smiles", ids))


def classify_smiles(smiles: str) -> str | None:
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
        return "disconnected"
    if charged:
        return "charged"
    if ring_count >= 2:
        return "multi-ring"
    if heteroatoms >= 3:
        return "heteroatom-rich"
    return "remaining"


def pubchem_rank(cid: str) -> str:
    return hashlib.sha256(f"{PUBCHEM_SEED}:{cid}".encode()).hexdigest()


def record_cid(record: str) -> str:
    marker = "> <PUBCHEM_COMPOUND_CID>"
    position = record.find(marker)
    if position < 0:
        return ""
    return record[position + len(marker) :].lstrip("\r\n").splitlines()[0].strip()


def read_enamine_smiles(path: Path) -> list[tuple[str, str]]:
    rows = []
    for index, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines()):
        if index == 0:
            continue
        line = raw_line.strip()
        if not line:
            continue
        smiles, catalog_id, *_ = line.split("\t")
        rows.append((smiles, catalog_id))
    return rows


def read_sdf_records(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")
    records = []
    for record in text.split("$$$$"):
        if record.strip():
            records.append(record.rstrip("\n") + "\n$$$$\n")
    return records


def normalize_sdf_record_header(record: str) -> str:
    lines = record.splitlines()
    if len(lines) >= 5 and not lines[0].strip() and not lines[1].strip() and lines[4].endswith(("V2000", "V3000")):
        lines.pop(1)
        return "\n".join(lines) + "\n"
    return record


def sdf_counts_line(record: str) -> str:
    lines = record.splitlines()
    if len(lines) < 4:
        raise SystemExit("SDF record is too short to contain a counts line")
    return lines[3].strip()


def read_single_sdf_record(path: Path) -> str:
    records = read_sdf_records(path)
    if len(records) != 1:
        raise SystemExit(f"{path} should contain exactly one SDF record")
    return records[0]


def sdf_property(record: str, property_name: str) -> str:
    pattern = re.compile(rf">\s*<{re.escape(property_name)}>[^\n]*\n([^\n\r]+)")
    match = pattern.search(record)
    if not match:
        raise SystemExit(f"SDF record is missing `{property_name}`")
    return match.group(1).strip()


def add_sdf_property(record: str, property_name: str, value: str) -> str:
    stripped = record.rstrip()
    if stripped.endswith("$$$$"):
        stripped = stripped[: -4].rstrip()
    return f"{stripped}\n\n> <{property_name}>\n{value}\n\n$$$$\n"


def write_manifests(root: Path, features: tuple[str, ...], fixtures: list[str]) -> None:
    manifest_dir = root / "features"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    array = ",\n".join(f'  "{fixture}"' for fixture in fixtures)
    for feature in features:
        notes = manifest_notes(root.name, feature)
        notes_text = ""
        if notes:
            notes_array = "\n".join(f'  "{note}",' for note in notes)
            notes_text = f"\nnotes = [\n{notes_array}\n]\n"
        text = (
            f'feature_id = "{feature}"\n'
            f'corpus_id = "{root.name}"\n'
            'reference_tool = "rdkit"\n'
            f'reference_version = "RDKit {rdBase.rdkitVersion}"\n'
            'comparison_mode = "implementation-golden"\n\n'
            f"fixtures = [\n{array},\n]\n"
            f"{notes_text}"
        )
        (manifest_dir / f"{feature}.toml").write_text(text, encoding="utf-8", newline="\n")


def manifest_notes(corpus: str, feature: str) -> tuple[str, ...]:
    if corpus == "pubchem-100k" and feature == "algo.canonical-ranking":
        return (
            "Reference RDKit goldens compare non-stereo atom symmetry partitions, not numeric rank labels.",
        )
    if corpus == "pubchem-100k" and feature == "stereo.cip":
        return (
            "Externally supplied PubChem isomeric SMILES packs provide large small-molecule CIP parity coverage.",
            "Goldens compare RDKit-backed CIP atom and bond descriptor maps, not bytewise SMILES spelling or internal stereo element IDs.",
        )
    return ()


def write_manual_manifests(root: Path, fixtures: list[str]) -> None:
    manifest_dir = root / "features"
    array = ",\n".join(f'  "{fixture}"' for fixture in fixtures)
    subjects = {
        "stereo.perception": (
            "stereo perception regression coverage",
            "semantic stereo perception reports and resulting stereo elements",
        ),
        "stereo.representation": (
            "stereo representation regression coverage",
            "semantic stereo elements, stereo groups, and source bond marks",
        ),
    }
    for feature, (coverage, comparison) in subjects.items():
        text = (
            f'feature_id = "{feature}"\n'
            f'corpus_id = "{root.name}"\n'
            'reference_tool = "pubchem-manual-semantic"\n'
            'reference_version = "PubChem PUG REST 2026-07-05"\n'
            'comparison_mode = "implementation-golden"\n\n'
            f"fixtures = [\n{array},\n]\n\n"
            "notes = [\n"
            f'  "Externally supplied PubChem isomeric SMILES packs provide large small-molecule {coverage}.",\n'
            f'  "Goldens compare {comparison}, not bytewise SMILES spelling.",\n'
            "]\n"
        )
        (manifest_dir / f"{feature}.toml").write_text(
            text, encoding="utf-8", newline="\n"
        )


def generate_goldens(repo: Path, corpus: str, features: tuple[str, ...]) -> None:
    script = repo / "validation" / "reference" / "rdkit" / "run_feature.py"
    for feature in features:
        subprocess.run(
            [sys.executable, str(script), "--feature", feature, "--corpus", corpus],
            cwd=repo,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


def generate_manual_goldens(repo: Path, corpus: str) -> None:
    for feature in PUBCHEM_MANUAL_SMILES_FEATURES:
        subprocess.run(
            [
                "cargo",
                "xtask",
                "validate",
                "--feature",
                feature,
                "--corpus",
                corpus,
                "--accept-implementation-goldens",
            ],
            cwd=repo,
            check=True,
        )


def pack_record(
    root: Path,
    path: Path,
    format_name: str,
    members: list[str],
    *,
    member_id_property: str | None = None,
    member_title_prefix: str | None = None,
) -> dict:
    record = {
        "path": path.relative_to(root).as_posix(),
        "format": format_name,
        "count": len(members),
        "members": members,
        "sha256": sha256(path),
    }
    if member_id_property is not None:
        record["member_id_property"] = member_id_property
    if member_title_prefix is not None:
        record["member_title_prefix"] = member_title_prefix
    return record


def chunks(items: list[Path], size: int) -> list[list[Path]]:
    return [items[index : index + size] for index in range(0, len(items), size)]


def one(items) -> Path:
    matches = list(items)
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one match, found {len(matches)}")
    return matches[0]


def reset_dir(path: Path) -> Path:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    raise SystemExit(main())
