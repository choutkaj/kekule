#!/usr/bin/env python3
"""Independently reconstruct descriptor masses from formulas and pinned sources.

This optional analysis never changes benchmark observations or tolerances.
RDKit is needed only when running the command, not for the dependency-free tests.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path


FIELDS = ("average_mass_da", "monoisotopic_mass_da")
# CODATA 2022 electron rest mass in unified atomic mass units.
# https://physics.nist.gov/cgi-bin/cuu/Value?meu
NATIVE_ELECTRON_MASS = 5.485799090441e-4


def reconstruct(formula, constants, electron_masses):
    """Return independent sums and a summation-roundoff bound, not uncertainty."""
    charge = formula["formal_charge"]
    if type(charge) is not int or not formula["terms"]:
        raise ValueError("formula needs an integer charge and nonempty terms")
    contributions = [[], []]
    count = 0
    seen = set()
    for term in formula["terms"]:
        symbol, isotope, amount = term["element"], term["isotope"], term["count"]
        if not isinstance(symbol, str) or not symbol:
            raise ValueError("invalid element symbol")
        if isotope is not None and (type(isotope) is not int or isotope <= 0):
            raise ValueError("invalid isotope label")
        if type(amount) is not int or not 0 < amount <= 2**53:
            raise ValueError("invalid constituent count")
        if (symbol, isotope) in seen:
            raise ValueError("duplicate formula term")
        seen.add((symbol, isotope))
        count += amount
        for items, mass in zip(contributions, constants(symbol, isotope), strict=True):
            if not math.isfinite(mass) or mass <= 0:
                raise ValueError("missing or invalid atomic mass")
            items.append(amount * mass)
    result = {}
    for field, items, electron_mass in zip(FIELDS, contributions, electron_masses, strict=True):
        correction = -charge * electron_mass
        computed = math.fsum([*items, correction])
        scale = math.fsum([*items, abs(correction)])
        if not math.isfinite(computed) or not math.isfinite(scale):
            raise ValueError("mass sum is not finite")
        # Atom-by-atom accumulation has at most count additions. Grouped
        # products, the final sum and charge correction need additional rounding
        # steps. The absolute constituent sum also covers cancellation in ions.
        bound = (count + 2 * len(items) + 4) * math.ulp(max(scale, 1.0))
        result[field] = {"reconstructed_da": computed, "roundoff_bound_da": bound}
    return result


def analyze_record(expected, actual, native, reference, reference_electron):
    if expected.get("status") != "ok" or actual.get("status") != "ok":
        return {"classification": "invalid_record_status"}
    for field in ("record_index", "title"):
        if expected.get(field) != actual.get(field):
            return {"classification": "record_identity_difference", "field": field}
    if expected["formula"] != actual["formula"]:
        return {"classification": "formula_difference"}
    try:
        native_values = reconstruct(actual["formula"], native, (NATIVE_ELECTRON_MASS,) * 2)
        reference_values = reconstruct(expected["formula"], reference, (0.0, reference_electron))
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        return {"classification": "unavailable_reconstruction", "message": str(error)}
    verified = True
    for record, values in ((actual, native_values), (expected, reference_values)):
        for field, value in values.items():
            observed = record.get(field)
            if type(observed) not in (int, float) or not math.isfinite(observed):
                return {"classification": "invalid_mass_observation", "field": field}
            residual = observed - value["reconstructed_da"]
            value.update(observed_da=observed, residual_da=residual)
            verified &= abs(residual) <= value["roundoff_bound_da"]
    deltas = {}
    charge = actual["formula"]["formal_charge"]
    for field in FIELDS:
        # Convention contributions have a signed native-minus-reference sense.
        reference_correction = reference_electron if field == "monoisotopic_mass_da" else 0.0
        charge_delta = charge * (reference_correction - NATIVE_ELECTRON_MASS)
        model_delta = native_values[field]["reconstructed_da"] - reference_values[field]["reconstructed_da"]
        deltas[field] = {
            "observed_native_minus_reference_da": actual[field] - expected[field],
            "atomic_data_contribution_da": model_delta - charge_delta,
            "charge_convention_contribution_da": charge_delta,
        }
    return {"classification": "mass_models_verified" if verified else "unexplained_mass",
            "native": native_values, "reference": reference_values, "deltas": deltas}


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load_models(cache):
    from rdkit import Chem, rdBase
    from rdkit.Chem import Descriptors

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("atomic_data_generator", root / "tools/atomic-data/generate.py")
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    sources = generator.acquire_sources(cache, offline=True)
    weights = generator.standard_weights(sources["ciaaw-abridged-2024.html"])
    abundant = generator.natural_isotopes(sources["ciaaw-isotopes-2024.html"])
    entries = generator.isotope_masses(sources["mass_1.mas20"])
    # Fail explicitly if the checked-in table has moved beyond this analysis's
    # source model. Never explain a newer table using stale pinned inputs.
    data_path = root / "crates/kekule/src/descriptors/data.rs"
    if generator.render_data(weights, abundant, entries) != data_path.read_text(encoding="utf-8"):
        raise ValueError("native atomic table differs from verified source regeneration")
    masses = {(z, a): mass for z, a, mass in entries}
    periodic = Chem.GetPeriodicTable()
    parser = Chem.SmilesParserParams()
    parser.removeHs = False
    proton = Chem.MolFromSmiles("[H+]", parser)
    reference_electron = periodic.GetMassForIsotope(1, 1) - Descriptors.ExactMolWt(proton)
    if not math.isfinite(reference_electron) or reference_electron <= 0:
        raise ValueError("invalid reference electron-mass probe")
    used = {}

    def constants(symbol, isotope):
        key = (symbol, isotope)
        if key not in used:
            z = periodic.GetAtomicNumber(symbol)
            if isotope is not None:
                native = (masses[z, isotope], masses[z, isotope])
                mass = periodic.GetMassForIsotope(z, isotope)
                reference = (mass, mass)
            else:
                native = (weights[z], masses[z, abundant[z]])
                reference = (periodic.GetAtomicWeight(z), periodic.GetMostCommonIsotopeMass(z))
            used[key] = {"element": symbol, "isotope": isotope, "native": native, "reference": reference}
        return used[key]

    provenance = {
        "rdkit_version": rdBase.rdkitVersion,
        "native_electron_mass_da": NATIVE_ELECTRON_MASS,
        "reference_electron_mass_da": reference_electron,
        "native_atomic_table_sha256": sha256(data_path),
        "source_parser_sha256": sha256(Path(spec.origin)),
        "sources": {name: {"url": generator.SOURCES[name][0], "sha256": sha256(path)}
                    for name, path in sources.items()},
    }
    return (lambda s, i: constants(s, i)["native"],
            lambda s, i: constants(s, i)["reference"], reference_electron, used, provenance)


def analyze_report(report_path, output, native, reference, reference_electron, provenance, used):
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if not report.get("complete") or report.get("error") is not None:
        raise ValueError("analysis requires a complete report without a run-level error")
    selected = [r for r in report["results"] if r["feature"] == "descriptor.molecular"]
    if not selected:
        raise ValueError("report has no molecular descriptor results")
    for result in selected:
        if result["cases"] == 0:
            continue
        ref = result["golden"]["reference"]
        if ref != {"tool": "rdkit", "version": provenance["rdkit_version"]}:
            raise ValueError("installed RDKit differs from the benchmark reference")
    cases = Path(report["cases"])
    # Benchmark case paths are written relative to the execution directory.
    # Require that path to resolve as recorded, rather than guess another file.
    details_path = output.with_suffix(".details.jsonl.gz")
    if output.exists() or details_path.exists():
        raise ValueError("analysis output already exists; choose a new output path")
    outcomes, classifications = Counter(), Counter()
    opener = gzip.open if cases.suffix == ".gz" else open
    output.parent.mkdir(parents=True, exist_ok=True)
    with opener(cases, "rt", encoding="utf-8") as stream, gzip.open(details_path, "xt", encoding="utf-8") as details:
        for line in stream:
            row = json.loads(line)
            if row["feature"] != "descriptor.molecular":
                continue
            outcomes[row["status"]] += 1
            if row["status"] not in ("agrees", "disagrees"):
                continue
            expected = row["expected"]["value"]["records"]
            actual = row["actual"]["value"]["records"]
            if not expected or len(expected) != len(actual):
                results = [{"classification": "record_count_difference"}]
            else:
                results = [analyze_record(e, a, native, reference, reference_electron)
                           for e, a in zip(expected, actual, strict=True)]
            for ordinal, result in enumerate(results):
                classifications[result["classification"]] += 1
                details.write(json.dumps({"dataset": row["dataset"], "id": row["id"],
                    "fixture": row["fixture"], "record_index": row["record_index"],
                    "component_ordinal": ordinal, **result}, allow_nan=False) + "\n")
    expected_rows = sum(r["cases"] + r["not_applicable"] for r in selected)
    if sum(outcomes.values()) != expected_rows:
        raise ValueError("case file row count differs from completed report")
    summary = {
        "complete": True, "report": str(report_path), "report_sha256": sha256(report_path),
        "cases": str(cases), "cases_sha256": sha256(cases),
        "details": str(details_path), "analyzer_sha256": sha256(Path(__file__)),
        "implementation": report["implementation"], "provenance": provenance,
        "case_outcomes": outcomes, "component_classifications": classifications,
        "used_atomic_constants": list(used.values()),
        "note": "Reconstruction roundoff is not experimental uncertainty or benchmark tolerance. Error/not-applicable cases are counted but not reconstructed.",
    }
    with output.open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--source-cache", required=True, type=Path,
                        help="Existing checksum-pinned cache from tools/atomic-data/generate.py.")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    native, reference, electron, used, provenance = load_models(args.source_cache)
    result = analyze_report(args.report, args.output, native, reference, electron, provenance, used)
    print(json.dumps({key: result[key] for key in ("case_outcomes", "component_classifications")}))
    counts = result["component_classifications"]
    return int(not counts or any(key != "mass_models_verified" for key in counts))


if __name__ == "__main__":
    raise SystemExit(main())
