"""Locations shared by the benchmark commands and offline provenance checks."""
import gzip
import json
from pathlib import Path
from zipfile import ZipFile

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
ARCHIVE = HERE / 'archive/previous-validation.zip'


def historical_bytes(name):
    with ZipFile(ARCHIVE) as archive:
        return archive.read(name)


def historical(name):
    data = historical_bytes(name)
    return json.loads(gzip.decompress(data) if name.endswith('.gz') else data)
