"""Select the small-molecule dataset from ChEMBL, the CCD and RDKit.

Run through `python benchmarks/bench.py select small`. Every rule lives in
`selection.toml`; this module only applies it to upstream metadata.
"""
from __future__ import annotations

import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

from kbench.datasets import Record, query_rows
from kbench.selection import pool_hash
from kbench.upstream import NotFound

ELEMENT = re.compile(r"[A-Z][a-z]?")
WORKERS = 6


def select(context) -> list[Record]:
    records = []
    for source in (rdkit_structures, rdkit_queries, ccd, chembl):
        added = source(context)
        print(f"{source.__name__}: {len(added)} records", file=sys.stderr, flush=True)
        records += added
    return records


def ordered(function, items, workers: int = WORKERS):
    """Map `function` over `items` concurrently, yielding results in input order."""
    with ThreadPoolExecutor(workers) as pool:
        for start in range(0, len(items), 4 * workers):
            yield from pool.map(function, items[start : start + 4 * workers])


def chembl(context) -> list[Record]:
    seen: set[str] = set()
    config = context.selection["chembl"]
    api = config["api"]
    status = context.client.json(f"{api}/status.json")
    release = status["chembl_db_version"]
    context.lock["chembl"] = {"release": release, "release_date": status["chembl_release_date"]}
    records = []
    for stratum in context.strata("chembl"):
        query = urllib.parse.urlencode({**config["required"], **stratum["filters"]})
        total = context.client.json(f"{api}/molecule.json?{query}&limit=1&only=molecule_chembl_id")[
            "page_meta"
        ]["total_count"]
        offsets = context.rng(stratum["id"]).sample(range(total), min(total, 3 * stratum["quota"]))

        def draw(offset):
            try:
                page = context.client.json(
                    f"{api}/molecule.json?{query}&limit=1&offset={offset}"
                    "&only=molecule_chembl_id,molecule_structures"
                )
            except urllib.error.HTTPError as error:
                # A record the API persistently fails to serve is unavailable,
                # like a missing structure, but only while the service itself
                # answers; during an outage selection stops and resumes later.
                if error.code >= 500 and service_up(api):
                    return "unavailable"
                raise
            return page["molecules"][0] if page["molecules"] else None

        chosen = []
        drawn = 0
        unavailable = 0
        for molecule in ordered(draw, offsets):
            if len(chosen) == stratum["quota"]:
                break
            drawn += 1
            if molecule == "unavailable":
                unavailable += 1
                continue
            if molecule is None:
                continue
            chembl_id = molecule["molecule_chembl_id"]
            structures = molecule.get("molecule_structures") or {}
            if chembl_id in seen or not structures.get("molfile") or not structures.get("canonical_smiles"):
                continue
            seen.add(chembl_id)
            url = f"{api}/molecule/{chembl_id}.json"
            files = (
                context.stage(f"chembl/{chembl_id}.mol", structures["molfile"].encode(), "mol", url, role="2d"),
                context.stage(
                    f"chembl/{chembl_id}.smi", (structures["canonical_smiles"] + "\n").encode(), "smiles", url
                ),
            )
            chosen.append(Record(f"chembl:{chembl_id}", stratum["id"], "chembl", chembl_id, files))
        context.lock["strata"][stratum["id"]] = {
            "source": "chembl",
            "release": release,
            "query": query,
            "pool_size": total,
            "drawn": drawn,
            "unavailable": unavailable,
        }
        records += chosen
    return records


def service_up(api: str) -> bool:
    request = urllib.request.Request(f"{api}/status.json", headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return response.status == 200
    except (urllib.error.URLError, TimeoutError, ConnectionError):
        return False


def ccd(context) -> list[Record]:
    config = context.selection["ccd"]
    query = {
        "query": {
            "type": "terminal",
            "service": "text_chem",
            "parameters": {"attribute": "rcsb_chem_comp_info.initial_release_date", "operator": "exists"},
        },
        "return_type": "mol_definition",
        "request_options": {"return_all_hits": True, "results_content_type": ["experimental"]},
    }
    result = context.client.json(config["search"], query)
    pool = sorted(hit["identifier"] for hit in result["result_set"] if not hit["identifier"].startswith("PRD_"))
    formulas = ccd_formulas(context, pool)
    common = set(config["common_elements"])
    seen: set[str] = set()
    records = []
    for stratum in context.strata("ccd"):
        candidates = pool
        if stratum.get("uncommon_elements"):
            candidates = [c for c in pool if set(ELEMENT.findall(formulas.get(c) or "")) - common]
        candidates = [c for c in candidates if c not in seen]
        order = context.rng(stratum["id"]).sample(candidates, len(candidates))
        missing = 0
        chosen = []

        def fetch(component):
            try:
                return component, context.client.get(config["ideal_sdf"].format(id=component))
            except NotFound:
                return component, None

        for component, data in ordered(fetch, order[: 2 * stratum["quota"]]):
            if len(chosen) == stratum["quota"]:
                break
            if data is None:
                missing += 1
                continue
            url = config["ideal_sdf"].format(id=component)
            files = (context.stage(f"ccd/{component}_ideal.sdf", data, "sdf", url, role="3d"),)
            chosen.append(Record(f"ccd:{component}", stratum["id"], "ccd", component, files))
            seen.add(component)
        records += chosen
        context.lock["strata"][stratum["id"]] = {
            "source": "ccd",
            "pool_size": len(candidates),
            "pool_sha256": pool_hash(candidates),
            "without_ideal_coordinates": missing,
        }
    return records


def ccd_formulas(context, components: list[str], batch: int = 500) -> dict[str, str]:
    """CCD formulas from the RCSB Data API, e.g. `C34 H32 Fe N4 O4` for HEM."""
    formulas = {}
    for start in range(0, len(components), batch):
        ids = ", ".join(f'"{c}"' for c in components[start : start + batch])
        query = {"query": f"{{ chem_comps(comp_ids: [{ids}]) {{ chem_comp {{ id formula }} }} }}"}
        for entry in context.client.json(context.selection["ccd"]["graphql"], query)["data"]["chem_comps"]:
            if entry and entry.get("chem_comp"):
                formulas[entry["chem_comp"]["id"]] = entry["chem_comp"]["formula"]
    return formulas


def rdkit_raw(context, path: str) -> tuple[str, bytes]:
    config = context.selection["rdkit"]
    url = config["raw"].format(commit=config["commit"], path=path)
    return url, context.client.get(url)


def rdkit_structures(context) -> list[Record]:
    config = context.selection["rdkit"]
    (stratum,) = context.strata("rdkit-structures")
    tree = context.client.json(
        f"https://api.github.com/repos/rdkit/rdkit/git/trees/{config['structures_tree']}?recursive=1"
    )
    if tree.get("truncated"):
        raise ValueError("GitHub truncated the RDKit test-data tree listing")
    paths = sorted(
        entry["path"]
        for entry in tree["tree"]
        if entry["type"] == "blob"
        and ".expected" not in entry["path"]
        and (
            (
                entry["path"].startswith("atropisomers/")
                and entry["path"].endswith(".sdf")
                and entry["path"].removeprefix("atropisomers/").startswith(tuple(config["families"]))
            )
            or re.fullmatch(r"chebi_[^/]*\.mol", entry["path"])
        )
    )
    records = []
    for path in paths:
        url, data = rdkit_raw(context, f"{config['structures_dir']}/{path}")
        # The atropisomer .sdf files hold one MOL block without an SDF delimiter.
        file = context.stage(f"rdkit/test_data/{path}", data, "mol", url)
        records.append(Record(f"rdkit:{path}", stratum["id"], "rdkit", path, (file,)))
    context.lock["strata"][stratum["id"]] = {
        "source": "rdkit",
        "commit": config["commit"],
        "tree": config["structures_tree"],
        "files": len(records),
    }
    return records


def rdkit_queries(context) -> list[Record]:
    config = context.selection["rdkit"]
    (stratum,) = context.strata("rdkit-queries")
    records = []
    for table in config["query_tables"]:
        url, data = rdkit_raw(context, table["path"])
        name = table["path"].rsplit("/", 1)[-1]
        staged = context.stage(f"rdkit/{name}", data, "smarts", url)
        if staged.sha256 != table["sha256"]:
            raise ValueError(f"{table['path']}: sha256 {staged.sha256} != pinned {table['sha256']}")
        for line, _query, _label in query_rows(data.decode("utf-8"), table["layout"]):
            file = type(staged)(staged.path, "smarts", staged.sha256, url, row=line, layout=table["layout"])
            source_id = f"{name}:{line}"
            records.append(Record(f"rdkit:{source_id}", stratum["id"], "rdkit", source_id, (file,)))
    context.lock["strata"][stratum["id"]] = {
        "source": "rdkit",
        "commit": config["commit"],
        "queries": len(records),
    }
    return records
