#!/usr/bin/env python3
"""Generate Biopython-backed golden data for macromolecular benchmark features."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
import warnings
import math
from pathlib import Path
from typing import Any


SUPPORTED_FEATURES = {
    "bio.secondary-structure.dssp",
    "io.mmcif.parse",
}

def import_biopython() -> dict[str, Any]:
    try:
        import Bio
        from Bio.PDB.DSSP import DSSP
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        from Bio.PDB.MMCIFParser import MMCIFParser
    except ImportError as error:
        raise SystemExit(
            "Biopython is not importable. Create the environment from "
            "benchmarks/reference/biopython/environment.yml before generating goldens."
        ) from error
    return {
        "version": Bio.__version__,
        "DSSP": DSSP,
        "MMCIF2Dict": MMCIF2Dict,
        "MMCIFParser": MMCIFParser,
    }


def evaluate(feature_id: str, fixture_path: Path, biopython: dict[str, Any]) -> dict[str, Any]:
    if feature_id == "bio.secondary-structure.dssp":
        return dssp_summary(fixture_path, biopython["MMCIFParser"], biopython["MMCIF2Dict"], biopython["DSSP"])
    if feature_id == "io.mmcif.parse":
        values = biopython["MMCIF2Dict"](str(fixture_path))
        name = values.pop('data_')
        # MMCIF2Dict reads one block. Refuse to silently merge multiple blocks.
        if any(key.lower().startswith('data_') for key in values):
            raise ValueError('Biopython MMCIF2Dict cannot represent multiple data blocks')
        return {'blocks':[{'name':name,'values':{key.lower():value for key,value in values.items()}}]}
    raise ValueError(f"unsupported feature: {feature_id}")


def dssp_reference(biopython_version: str) -> dict[str, Any]:
    executable = shutil.which("mkdssp")
    if executable is None:
        raise SystemExit(
            "mkdssp is not available. Recreate the environment from "
            "benchmarks/reference/biopython/environment.yml."
        )
    version = subprocess.run(
        [executable, "--version"],
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    ).stdout.strip()
    return {
        "tool": "biopython",
        "version": f"Biopython {biopython_version} / {version}",

    }


def dssp_summary(
    fixture_path: Path,
    MMCIFParser: Any,
    MMCIF2Dict: Any,
    DSSP: Any,
) -> dict[str, Any]:
    parser = MMCIFParser(QUIET=True)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = parser.get_structure(fixture_path.stem, str(fixture_path))[0]
            with tempfile.TemporaryDirectory(prefix="kekule-dssp-") as temp_dir:
                legacy = subprocess.run(
                    ['mkdssp', '--output-format=dssp', '--quiet', str(fixture_path)],
                    check=True, capture_output=True, text=True, timeout=120,
                ).stdout
                legacy_path = Path(temp_dir) / 'assignments.dssp'
                legacy_path.write_text(legacy, encoding='utf-8')
                partners = beta_partner_indices(legacy)
                assignments = DSSP(
                    model,
                    str(legacy_path),
                    dssp="mkdssp",
                    file_type="DSSP",
                )
                annotated_path = Path(temp_dir) / "annotated.cif"
                subprocess.run(
                    [
                        "mkdssp",
                        "--output-format=mmcif",
                        "--quiet",
                        str(fixture_path),
                        str(annotated_path),
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                extended_rows = dssp_extended_rows(MMCIF2Dict(str(annotated_path)))
    except Exception as error:
        return {
            "status": "reference_error",
            "error": error.__class__.__name__,
            "message": str(error),
            "residues": [],
        }
    if len(assignments) == 0:
        return {"status": "no_analyzable_residues", "residues": []}

    keys = list(assignments.keys())
    if len(keys) != len(extended_rows):
        raise RuntimeError(
            "Biopython legacy DSSP and DSSP4 mmCIF output contain different "
            f"residue counts: {len(keys)} versus {len(extended_rows)}"
        )
    keys_by_dssp_index = {int(assignments[key][0]): key for key in keys}
    residues = []
    previous_dssp_index = None
    previous_label_chain = None
    for ordinal, (key, extended) in enumerate(zip(keys, extended_rows, strict=True)):
        chain_id, residue_id = key
        _, sequence_id, insertion_code = residue_id
        value = assignments[key]
        dssp_index = int(value[0])
        if previous_dssp_index is None:
            chain_break = "new_chain"
        elif dssp_index != previous_dssp_index + 1:
            chain_break = (
                "new_chain"
                if extended["label_chain_id"] != previous_label_chain
                else "gap"
            )
        else:
            chain_break = "none"
        residues.append(
            {
                "chain_id": chain_id,
                "sequence_id": sequence_id,
                "insertion_code": normalize_missing(str(insertion_code).strip()),
                "label_chain_id": extended["label_chain_id"],
                "label_sequence_id": extended["label_sequence_id"],
                "residue_name": extended["residue_name"],
                "residue_one_letter": value[1],
                "secondary_structure": " " if value[2] == "-" else value[2],
                "chain_break": chain_break,
                "phi_degrees": dssp_optional_angle(value[4]),
                "psi_degrees": dssp_optional_angle(value[5]),
                "tco": extended["tco"],
                "kappa_degrees": extended["kappa_degrees"],
                "alpha_degrees": extended["alpha_degrees"],
                "helix_positions": extended["helix_positions"],
                "sheet": extended["sheet"],
                "strand": extended["strand"],
                "ladders": extended["ladders"],
                "beta_parallel": extended["beta_parallel"],
                "beta_partners": [partner_identity(keys_by_dssp_index[index]) if index else None
                                  for index in partners[dssp_index]],
                "omega_degrees": omega_angle(model, keys, assignments, ordinal),
                "acceptors": [
                    dssp_bond(value[6], value[7], dssp_index, keys_by_dssp_index),
                    dssp_bond(value[10], value[11], dssp_index, keys_by_dssp_index),
                ],
                "donors": [
                    dssp_bond(value[8], value[9], dssp_index, keys_by_dssp_index),
                    dssp_bond(value[12], value[13], dssp_index, keys_by_dssp_index),
                ],
            }
        )
        previous_dssp_index = dssp_index
        previous_label_chain = extended["label_chain_id"]
    return {"status": "ok", "residues": residues}


def beta_partner_indices(text):
    """Read the two explicit BP columns; do not reconstruct partners from ladders."""
    rows = {}
    started = False
    for line in text.splitlines():
        if '  #  RESIDUE' in line:
            started = True
            continue
        if not started or len(line) < 34 or line[9] == ' ':
            continue
        index = int(line[:5])
        if index in rows:
            raise ValueError('duplicate DSSP residue index')
        rows[index] = [int(line[25:29]), int(line[29:33])]
    return rows


def partner_identity(key):
    chain, (_, sequence, insertion) = key
    return {'partner_chain_id': chain, 'partner_sequence_id': sequence,
            'partner_insertion_code': normalize_missing(insertion.strip())}


def omega_angle(model, keys, assignments, index):
    """Biopython dihedral of CA(i), C(i), N(i+1), CA(i+1), in degrees."""
    from Bio.PDB.vectors import calc_dihedral
    if index + 1 == len(keys):
        return None
    key, following = keys[index:index + 2]
    if key[0] != following[0] or assignments[following][0] != assignments[key][0] + 1:
        return None
    def residue(key):
        return next(r for r in model[key[0]] if r.id[1:] == key[1][1:] and r.id[0] != 'W')
    a, b = residue(key), residue(following)
    return math.degrees(calc_dihedral(a['CA'].get_vector(), a['C'].get_vector(),
                                    b['N'].get_vector(), b['CA'].get_vector()))


def dssp_extended_rows(document: dict[str, Any]) -> list[dict[str, Any]]:
    prefix = "_dssp_struct_summary."

    def column(name: str) -> list[str | None]:
        return normalize_mmcif_column(document.get(f"{prefix}{name}"))

    fields = (
        "label_asym_id",
        "label_seq_id",
        "label_comp_id",
        "helix_3_10",
        "helix_alpha",
        "helix_pi",
        "helix_pp",
        "sheet",
        "strand",
        "ladder_1",
        "ladder_2",
        "TCO",
        "kappa",
        "alpha",
    )
    columns = {field: column(field) for field in fields}
    row_count = len(columns["label_asym_id"])

    ladder_ids = normalize_mmcif_column(document.get("_dssp_struct_ladder.id"))
    ladder_types = normalize_mmcif_column(document.get("_dssp_struct_ladder.type"))
    type_by_ladder = {
        ladder_id: ladder_type
        for ladder_id, ladder_type in zip(ladder_ids, ladder_types, strict=True)
    }

    rows = []
    for row in range(row_count):
        ladder_tokens = [
            normalize_missing(columns["ladder_1"][row]),
            normalize_missing(columns["ladder_2"][row]),
        ]
        rows.append(
            {
                "label_chain_id": columns["label_asym_id"][row],
                "label_sequence_id": int(columns["label_seq_id"][row]),
                "residue_name": columns["label_comp_id"][row],
                "tco": dssp_optional_number(columns["TCO"][row]),
                "kappa_degrees": dssp_optional_number(columns["kappa"][row]),
                "alpha_degrees": dssp_optional_number(columns["alpha"][row]),
                "helix_positions": [
                    dssp_helix_position(columns["helix_3_10"][row]),
                    dssp_helix_position(columns["helix_alpha"][row]),
                    dssp_helix_position(columns["helix_pi"][row]),
                    dssp_helix_position(columns["helix_pp"][row]),
                ],
                "sheet": dssp_identifier(columns["sheet"][row], one_based=True),
                "strand": dssp_identifier(columns["strand"][row], one_based=True),
                "ladders": [
                    dssp_identifier(token, one_based=False) for token in ladder_tokens
                ],
                "beta_parallel": [
                    None
                    if token is None
                    else type_by_ladder[token].lower() == "parallel"
                    for token in ladder_tokens
                ],
            }
        )
    return rows


def dssp_optional_number(value: Any) -> float | None:
    normalized = normalize_missing(None if value is None else str(value))
    return None if normalized is None else float(normalized)


def dssp_helix_position(value: Any) -> str:
    normalized = normalize_missing(None if value is None else str(value))
    if normalized is None:
        return "none"
    return {">": "start", "<": "end", "X": "start_and_end"}.get(
        normalized, "middle"
    )


def dssp_identifier(value: Any, *, one_based: bool) -> int | None:
    normalized = normalize_missing(None if value is None else str(value))
    if normalized is None:
        return None
    if normalized.isdigit():
        identifier = int(normalized)
    else:
        identifier = 0
        for character in normalized.upper():
            if not "A" <= character <= "Z":
                raise ValueError(f"unsupported DSSP identifier {normalized!r}")
            identifier = identifier * 26 + ord(character) - ord("A") + 1
        identifier -= 1
    return identifier + 1 if one_based else identifier


def dssp_optional_angle(value: Any) -> float | None:
    angle = float(value)
    return None if angle == 360.0 else angle


def dssp_bond(
    relative_index: Any,
    energy: Any,
    dssp_index: int,
    keys_by_dssp_index: dict[int, Any],
) -> dict[str, Any] | None:
    relative_index = int(relative_index)
    energy = float(energy)
    if relative_index == 0 and energy == 0.0:
        return None
    partner_key = keys_by_dssp_index.get(dssp_index + relative_index)
    if partner_key is None:
        raise RuntimeError(
            f"DSSP bond from index {dssp_index} points to missing relative index "
            f"{relative_index}"
        )
    partner_chain, partner_residue = partner_key
    _, partner_sequence, partner_insertion = partner_residue
    return {
        "partner_chain_id": partner_chain,
        "partner_sequence_id": partner_sequence,
        "partner_insertion_code": normalize_missing(str(partner_insertion).strip()),
        "energy_kcal_per_mol": energy,
    }


def normalize_mmcif_column(value: Any) -> list[str | None]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def normalize_missing(value: str | None) -> str | None:
    if value is None or value in {"", ".", "?"}:
        return None
    return value
