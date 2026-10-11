"""Fetches the bundled Chemical Component Dictionary entries and locks them.

Run from the repository root:

    python tools/components/fetch.py            # download and lock
    python tools/components/fetch.py --lock     # lock the files already present

Each entry is the unmodified wwPDB Chemical Component Dictionary definition
served by RCSB. The dictionary is public domain (CC0 1.0).
"""
import argparse
import datetime
import hashlib
import json
import re
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "crates/kekule/data/components"
LOCK = HERE / "sources.lock.json"
URL = "https://files.rcsb.org/ligands/view/{}.cif"

# Grouped as documented in README.md. Changing the selection requires
# updating the Rust bundle list in crates/kekule/src/io/components/bundled.rs.
SELECTION = {
    "amino acids": [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "SEC", "PYL", "MSE",
    ],
    "nucleotides": ["A", "C", "G", "U", "DA", "DC", "DG", "DT"],
    "water": ["HOH", "DOD"],
    "ions": [
        "ZN", "MG", "CA", "NA", "K", "CL", "MN", "FE", "FE2", "CU", "CU1",
        "NI", "CO", "CD", "HG", "IOD", "BR",
    ],
    "cofactors": [
        "HEM", "HEC", "FAD", "FMN", "NAD", "NAI", "NAP", "NDP", "ATP", "ADP",
        "AMP", "GTP", "GDP", "SAM", "SAH", "COA", "ACO", "PLP", "TPP", "BTN",
        "SF4", "FES", "F3S", "CLA",
    ],
    "additives": ["SO4", "PO4", "GOL", "EDO", "ACT", "PEG", "DMS", "MPD"],
}


def identifiers():
    return [identifier for group in SELECTION.values() for identifier in group]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def modified_date(path):
    text = Path(path).read_text(encoding="utf-8")
    match = re.search(r"^_chem_comp\.pdbx_modified_date\s+(\S+)", text, re.MULTILINE)
    return match.group(1) if match else None


def fetch():
    DATA.mkdir(parents=True, exist_ok=True)
    for identifier in identifiers():
        with urllib.request.urlopen(URL.format(identifier), timeout=60) as response:
            (DATA / f"{identifier}.cif").write_bytes(response.read())


def lock(retrieved):
    entries = []
    for identifier in identifiers():
        path = DATA / f"{identifier}.cif"
        entries.append(
            {
                "id": identifier,
                "path": path.relative_to(ROOT).as_posix(),
                "url": URL.format(identifier),
                "bytes": path.stat().st_size,
                "sha256": digest(path),
                "modified": modified_date(path),
            }
        )
    record = {
        "schema": 1,
        "source": "wwPDB Chemical Component Dictionary entries served by RCSB",
        "license": "CC0 1.0 (wwPDB archive data)",
        "retrieved": retrieved,
        "groups": SELECTION,
        "components": entries,
    }
    LOCK.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", action="store_true", help="lock existing files only")
    parser.add_argument("--retrieved", default=datetime.date.today().isoformat())
    arguments = parser.parse_args()
    if not arguments.lock:
        fetch()
    lock(arguments.retrieved)


if __name__ == "__main__":
    main()
