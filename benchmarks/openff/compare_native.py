"""Compare the native parameterizer with the immutable audited OpenFF observations.

Run with standard-library Python; no scientific reference tools are required.
Numerical parameter comparisons and forced-inference references are supplied by
the additional reference generator when requested.
"""
import argparse
import gzip
import json
import subprocess
from pathlib import Path
from audit import HERE, CHARGE_ATOL, digest


def canonical(atoms, handler):
    atoms = tuple(atoms)
    if handler == "ImproperTorsions":
        outer = sorted((atoms[0], atoms[2], atoms[3]))
        return (outer[0], atoms[1], *outer[1:])
    return min(atoms, atoms[::-1])


def compare(reference, actual):
    errors = []
    if actual.get("status") != "ok":
        return [actual.get("message", "missing native result")], {}
    maps = reference["charge_atom_maps"]
    if len(set(actual["maps"])) != len(maps) or set(actual["maps"]) != set(maps):
        return ["atom identity mismatch"], {}
    order = [actual["maps"].index(tag) for tag in maps]
    if actual["fixed_h_inchi"] != reference["fixed_h_inchi"]:
        errors.append("fixed-H InChI mismatch")
    for handler, rows in actual["labels"].items():
        expected = {canonical([maps[i] for i in row["atoms"]], handler): (row["id"], row["smirks"])
                    for row in reference["assigned_parameters"][handler]}
        observed = {canonical(row["maps"], handler): (row["id"], row["smirks"]) for row in rows}
        if expected != observed:
            errors.append(dict(handler=handler, expected_only=[str((k, v)) for k, v in expected.items() if observed.get(k) != v],
                               actual_only=[str((k, v)) for k, v in observed.items() if expected.get(k) != v]))
    measures = {}
    if reference["features"]["status"] == "ok":
        if actual["features"]["status"] != "ok":
            errors.append(actual["features"])
        else:
            expected = [[] for _ in maps]
            for feature in reference["features"]["values"]:
                for row, value in zip(expected, feature["values"]):
                    row.extend(value if isinstance(value, list) else [value])
            observed = [actual["features"]["values"][i] for i in order]
            if len(observed) != len(expected) or any(len(a) != len(b) for a, b in zip(observed, expected)):
                errors.append("feature shape mismatch")
            else:
                delta = max(abs(a-b) for left, right in zip(observed, expected) for a,b in zip(left,right))
                measures["feature_error"] = delta
                if delta > 1e-6:
                    errors.append(dict(feature_mismatches=[(maps[i], left, right) for i,(left,right) in enumerate(zip(observed,expected)) if any(abs(a-b)>1e-6 for a,b in zip(left,right))]))
    elif actual["features"]["status"] != "error":
        errors.append("missing expected feature domain rejection")
    if reference["charges"]["status"] == "ok":
        if actual["charges"]["status"] != "ok":
            errors.append(actual["charges"])
        else:
            values = actual["charges"]["values"]
            delta = max(abs(values[i]-q) for i,q in zip(order,reference["charges"]["toolkit_normalized"]))
            measures["charge_error_e"] = delta
            if delta > CHARGE_ATOL:
                errors.append(dict(charge_error_e=delta))
    elif actual["charges"]["status"] != "error":
        errors.append("missing expected charge domain rejection")
    values = actual["system"]["charges"]
    expected = [row["charge"] for row in sorted(reference["system"]["charges"],key=lambda r:r["atoms"])]
    delta = max(abs(values[i]-q) for i,q in zip(order,expected))
    measures["system_charge_error_e"] = delta
    if delta > CHARGE_ATOL:
        errors.append(dict(system_charge_error_e=delta))
    expected_collections = reference["system"]["collections"]
    for native, handler in [("bonds","Bonds"),("angles","Angles"),("constraints","Constraints"),("propers","ProperTorsions"),("impropers","ImproperTorsions"),("vdw","vdW")]:
        observed = actual["system"][native]
        count = sum(len(r["parameter"]["terms"]) for r in observed) if native in ("propers","impropers") else len(observed)
        if count != expected_collections[handler]["assignments"]:
            errors.append(dict(handler=handler, count=count, expected_count=expected_collections[handler]["assignments"]))
    return errors, measures


def run(binary, model, output):
    reference = json.loads(gzip.decompress((HERE / "reference.json.gz").read_bytes()))
    rows = reference["records"]
    process = subprocess.run([str(binary.resolve()), str(model.resolve())],
        input="".join(json.dumps(dict(smiles=r["mapped_smiles"]))+"\n" for r in rows),
        text=True, capture_output=True, check=True)
    actual = [json.loads(line) for line in process.stdout.splitlines()]
    if len(actual) != len(rows):
        raise ValueError("Native output count mismatch")
    records = []
    for ref, native in zip(rows, actual):
        errors, measures = compare(ref, native)
        records.append(dict(id=ref["input"]["id"], passed=not errors, errors=errors, measures=measures, native=native))
        print(json.dumps({k:v for k,v in records[-1].items() if k!="native"}))
    report = dict(schema=1, charge_atol_e=CHARGE_ATOL, reference_sha256=digest(HERE/"reference.json.gz"),
                  binary_sha256=digest(binary), records=records)
    output.write_text(json.dumps(report,allow_nan=False)+"\n",encoding="utf-8")
    return int(any(not r["passed"] for r in records))


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary",type=Path,required=True)
    p.add_argument("--model",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    raise SystemExit(run(a.binary,a.model,a.output))
