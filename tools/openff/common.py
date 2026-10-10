"""Paths, pins and checks shared by the OpenFF data exporters."""
import hashlib
import json
from enum import Enum
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MODEL_HASH = "7981e7f5b0b1e424c9e10a40d9e7606d96dcd3dd2b095cb4eeff6829f92238ee"


def describe_configuration(value):
    # NAGL's activation fields contain Python classes, not JSON strings.
    # Preserve their qualified identities; never silently stringify other objects.
    if isinstance(value, type):
        return dict(python_type=value.__module__ + "." + value.__qualname__)
    if isinstance(value, Enum):
        return value.value
    raise TypeError(f"Unsupported report value: {type(value).__name__}")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_sources():
    """Check every locked upstream source against its recorded SHA-256."""
    lock = json.loads((HERE / "sources.lock.json").read_text(encoding="utf-8"))
    for source in lock["sources"]:
        if digest(ROOT / source["path"]) != source["sha256"]:
            raise ValueError(f"Source checksum mismatch: {source['path']}")
    return lock
