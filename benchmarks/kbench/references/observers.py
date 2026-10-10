"""Which reference observer serves which tool, and their source fingerprint.

Importable without the reference toolkits: a run only fingerprints the
observer sources to detect stored references made by different code.
"""
from __future__ import annotations

import hashlib
import importlib
from pathlib import Path

MODULES = {"rdkit": "kbench.references.rdkit_small", "bio": "kbench.references.bio"}
SOURCES = {"rdkit": "rdkit_small.py", "bio": "bio.py"}
USED = {"small": ("rdkit",), "bio": ("bio",)}
HERE = Path(__file__).resolve().parent


def fingerprint(dataset: str) -> dict[str, str]:
    """SHA-256 of the source of each observer a dataset uses, normalized to
    LF line endings, so editing one observer leaves other datasets valid."""
    return {
        tool: hashlib.sha256((HERE / SOURCES[tool]).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        for tool in USED[dataset]
    }


def tools(dataset: str) -> dict[str, str]:
    """Versions of the toolkits behind a dataset's references (imports them)."""
    return {tool: importlib.import_module(MODULES[tool]).TOOL for tool in USED[dataset]}


def prepare(kind: str, text: str) -> str:
    """Reference-side preparation of an input that both observers then read."""
    if kind == "single-conformer":
        from . import bio

        try:
            return bio.single_conformer(text)
        except Exception as error:  # an unpreparable input fails on both sides alike
            return f"# preparation failed: {type(error).__name__}: {error}\n"
    raise ValueError(f"unknown preparation {kind!r}")
