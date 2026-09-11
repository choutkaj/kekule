"""Shared transport for optional reference checks; chemical comparisons stay separate."""

import hashlib
import json
import subprocess

RDKIT_VERSION = "2026.03.6"
DESCRIPTORS = {
    "r": "LowerR", "s": "LowerS", "m": "LowerM", "p": "LowerP",
    "z": "SeqCis", "e": "SeqTrans", "seqCis": "SeqCis", "seqTrans": "SeqTrans",
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_report(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as output:
        output.write(json.dumps(report, indent=2) + "\n")


def run_probe(probe, lines, timeout=30):
    """Return one JSON response or explicit failure per input line, in input order."""
    if not lines:
        return []
    try:
        process = subprocess.run(
            [str(probe)], input="\n".join(lines) + "\n", capture_output=True,
            text=True, encoding="utf-8", timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        if len(lines) > 1:
            return [run_probe(probe, [line], timeout)[0] for line in lines]
        return [{"status": "timeout", "seconds": timeout}]
    except OSError as error:
        return [{"status": "process_error", "message": str(error)} for _ in lines]
    responses = process.stdout.splitlines()
    if process.returncode or len(responses) != len(lines):
        if len(lines) > 1:
            return [run_probe(probe, [line], timeout)[0] for line in lines]
        return [{"status": "process_error", "exit_code": process.returncode,
                 "stdout": process.stdout, "stderr": process.stderr}]
    try:
        values = [json.loads(response) for response in responses]
        if not all(isinstance(value, dict) and "status" in value for value in values):
            raise ValueError("probe response must contain an object with status")
        return values
    except (ValueError, TypeError) as error:
        return [{"status": "protocol_error", "message": str(error)} for _ in lines]
