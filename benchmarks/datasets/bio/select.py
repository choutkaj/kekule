"""Select the biomolecule dataset from the wwPDB versioned archive.

Run through `python benchmarks/bench.py select bio`. Every rule lives in
`selection.toml`; this module only applies it to RCSB metadata.
"""
from __future__ import annotations

import sys

from kbench.datasets import Record
from kbench.selection import pool_hash
from kbench.upstream import NotFound


def terminal(parameters: dict) -> dict:
    parameters = dict(parameters)
    service = parameters.pop("service", "text")
    return {"type": "terminal", "service": service, "parameters": parameters}


def query(cap: dict, stratum: dict) -> dict:
    nodes = [terminal(cap)] + [terminal(node) for node in stratum.get("all", [])]
    if "any" in stratum:
        nodes.append({"type": "group", "logical_operator": "or", "nodes": [terminal(n) for n in stratum["any"]]})
    return {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {"return_all_hits": True, "results_content_type": ["experimental"]},
    }


def select(context) -> list[Record]:
    config = context.selection["rcsb"]
    seen: set[str] = set()
    records = []
    for stratum in context.strata("wwpdb"):
        result = context.client.json(config["search"], query(config["cap"], stratum))
        pool = sorted(hit["identifier"] for hit in result["result_set"])
        forced = [entry for entry in stratum.get("include", []) if entry not in seen]
        order = forced + [entry for entry in context.rng(stratum["id"]).sample(pool, len(pool)) if entry not in forced]
        chosen = []
        unavailable = 0
        for entry in order:
            if len(chosen) == stratum["quota"]:
                break
            if entry in seen:
                continue
            record = fetch(context, config, stratum["id"], entry)
            if record is None:
                unavailable += 1
                continue
            seen.add(entry)
            chosen.append(record)
        context.lock["strata"][stratum["id"]] = {
            "source": "wwpdb",
            "pool_size": len(pool),
            "pool_sha256": pool_hash(pool),
            "included": stratum.get("include", []),
            "unavailable": unavailable,
        }
        print(f"{stratum['id']}: {len(chosen)} of pool {len(pool)}", file=sys.stderr, flush=True)
        records += chosen
    return records


def fetch(context, config: dict, stratum: str, entry: str) -> Record | None:
    accession = context.client.json(config["entry"].format(id=entry))["rcsb_accession_info"]
    code = entry.lower()
    url = config["versioned"].format(
        hash=code[1:3], id=code, major=accession["major_revision"], minor=accession["minor_revision"]
    )
    try:
        data = context.client.get(url)
    except NotFound:
        return None
    name = url.rsplit("/", 1)[-1]
    file = context.stage(f"pdb/{name}", data, "mmcif", url)
    return Record(f"pdb:{entry}", stratum, "wwpdb", entry, (file,))
