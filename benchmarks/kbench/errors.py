"""Stable error kinds from observer messages.

Observers report a coarse stage (`parse`, `sanitize`, `stereo`, ...). These
reviewed patterns split a stage by cause so distinct failures get distinct
signatures. The first matching pattern wins; unmatched messages keep the
stage alone. Patterns match messages, never inputs.
"""
from __future__ import annotations

import re

PATTERNS = (
    ("kekule", "parse", r"R-group|query atom|\bR#|unsupported element symbol `(?:R\d*|X)`", "query-atom"),
    ("kekule", "parse", r"UnassembledTetrahedralBondMark", "unassembled-wedge"),
    ("kekule", "parse", r"unsupported element symbol `[A-Z]{2}`", "uppercase-element"),
    ("kekule", "parse", r"exceeds configured token count", "token-limit"),
    ("kekule", "interpret", r"unknown atom-site element `D`", "deuterium"),
    ("kekule", "interpret", r"no complete alternate configuration", "altloc"),
    ("kekule", "parse", r"carriers", "stereo-carriers"),
    ("kekule", "write", r"wider|field width|fixed-width|too wide", "field-width"),
    ("kekule", "write", r"[Uu]nsupported|cannot (?:encode|represent)", "unsupported"),
    ("kekule", "perception", r"[Vv]alence", "valence"),
    ("kekule", "cip", r"exhaust|limit|bound", "resource-limit"),
    ("kekule", "resonance", r"work limit", "resource-limit"),
    ("reference", "resonance", r"work limit", "resource-limit"),
    ("reference", "sanitize", r"[Vv]alence", "valence"),
    ("reference", "sanitize", r"[Kk]ekulize", "kekulize"),
    ("reference", "parse", r".", "unreadable"),
)


def refine(side: str, kind: str, message: str) -> str:
    for pattern_side, pattern_kind, pattern, label in PATTERNS:
        if side == pattern_side and kind == pattern_kind and re.search(pattern, message or ""):
            return f"{kind}:{label}"
    return kind
