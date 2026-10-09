"""Package an exported Ash bundle as the data of the `kekule-openff-ash` crate.

Run after `export_model.py` on the checksum-pinned Ash checkpoint:

    python benchmarks/openff/scripts/package_ash.py target/ash-export crates/kekule-openff-ash/data

The manifest is copied unchanged. The float32 weights are stored losslessly as
four byte planes (all first bytes, then all second bytes, ...) compressed with
zlib, which keeps the crate below the crates.io package size limit. The script
verifies that decoding reproduces the manifest's `weights_sha256` exactly.
Only the standard library is required.
"""
import argparse
import hashlib
import json
import zlib
from pathlib import Path

ASH_MODEL = "openff-gnn-am1bcc-1.0.0.pt"
ASH_CHECKPOINT = "7981e7f5b0b1e424c9e10a40d9e7606d96dcd3dd2b095cb4eeff6829f92238ee"


def shuffle(weights: bytes) -> bytes:
    if len(weights) % 4:
        raise ValueError("float32 weights must have a multiple of four bytes")
    return b"".join(weights[plane::4] for plane in range(4))


def unshuffle(planes: bytes) -> bytes:
    count = len(planes) // 4
    out = bytearray(len(planes))
    for plane in range(4):
        out[plane::4] = planes[plane * count:(plane + 1) * count]
    return bytes(out)


def package(source: Path, destination: Path) -> None:
    manifest_bytes = (source / "model.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest.get("schema") != 2 or manifest.get("model") != ASH_MODEL \
            or manifest.get("checkpoint_sha256") != ASH_CHECKPOINT:
        raise ValueError("source must be a schema-2 export of the pinned Ash checkpoint")
    weights = (source / "weights.bin").read_bytes()
    if hashlib.sha256(weights).hexdigest() != manifest["weights_sha256"]:
        raise ValueError("weights.bin does not match the manifest checksum")
    encoded = zlib.compress(shuffle(weights), 9)
    if hashlib.sha256(unshuffle(zlib.decompress(encoded))).hexdigest() != manifest["weights_sha256"]:
        raise AssertionError("lossless round trip failed")
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "model.json").write_bytes(manifest_bytes)
    (destination / "weights.planes.zlib").write_bytes(encoded)
    print(json.dumps({
        "model.json": hashlib.sha256(manifest_bytes).hexdigest(),
        "weights.planes.zlib": hashlib.sha256(encoded).hexdigest(),
        "weights_sha256": manifest["weights_sha256"],
        "encoded_bytes": len(encoded),
    }, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source", type=Path, help="export_model.py output directory")
    parser.add_argument("destination", type=Path, help="crate data directory")
    arguments = parser.parse_args()
    package(arguments.source, arguments.destination)
