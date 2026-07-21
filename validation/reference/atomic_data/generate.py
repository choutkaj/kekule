#!/usr/bin/env python3
"""Regenerate the molecular-descriptor atomic reference tables."""

from __future__ import annotations

import argparse
import hashlib
import html
import re
import urllib.request
from pathlib import Path


SOURCES = {
    "ciaaw-abridged-2024.html": (
        "https://ciaaw.org/abridged-atomic-weights.htm",
        "f9e9554471749c55a624aec55151922470a7f4104c62811eb194fed9731b907d",
    ),
    "ciaaw-isotopes-2024.html": (
        "https://ciaaw.org/isotopic-abundances.htm",
        "1ab770a48bb539d26ad2992a94b4c60d51743058b67f106796b754cd5035423b",
    ),
    "mass_1.mas20": (
        "https://amdc.impcas.ac.cn/masstables/Ame2020/mass_1.mas20",
        "e8599c6d7f724fac91934e59f1b9de8fb8f63e820f4b39456b790665ed2a3307",
    ),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "cache",
        help="Source download/cache directory (ignored by git).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[3]
        / "crates"
        / "molecular"
        / "src"
        / "descriptors"
        / "data.rs",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Require all checksum-pinned source files to exist in the cache.",
    )
    args = parser.parse_args()

    sources = acquire_sources(args.cache_dir, args.offline)
    output = render_data(
        standard_weights(sources["ciaaw-abridged-2024.html"]),
        natural_isotopes(sources["ciaaw-isotopes-2024.html"]),
        isotope_masses(sources["mass_1.mas20"]),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(args.output.resolve())
    return 0


def acquire_sources(cache_dir: Path, offline: bool) -> dict[str, Path]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for filename, (url, expected_sha256) in SOURCES.items():
        path = cache_dir / filename
        if not path.exists():
            if offline:
                raise SystemExit(f"missing offline source: {path}")
            with urllib.request.urlopen(url) as response:
                path.write_bytes(response.read())
        actual_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual_sha256 != expected_sha256:
            raise SystemExit(
                f"{path} SHA-256 mismatch: expected {expected_sha256}, "
                f"found {actual_sha256}"
            )
        paths[filename] = path
    return paths


def table_cells(row: str) -> list[str]:
    return [
        re.sub(r"<[^>]+>", "", html.unescape(cell))
        .replace("\xa0", " ")
        .strip()
        for cell in re.findall(r"<td[^>]*>(.*?)</td>", row, re.I | re.S)
    ]


def standard_weights(path: Path) -> list[float | None]:
    source = path.read_text(encoding="utf-8")
    values: dict[int, float | None] = {}
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", source, re.I | re.S):
        fields = table_cells(row)
        if not fields or not fields[0].isdigit():
            continue
        atomic_number = int(fields[0])
        if not 1 <= atomic_number <= 118:
            continue
        values[atomic_number] = (
            None
            if "—" in fields[3]
            else float(re.search(r"\d+(?:\.\d+)?", fields[3]).group())
        )
    if len(values) != 118:
        raise SystemExit(f"expected 118 CIAAW weight rows, found {len(values)}")
    return [None, *(values[index] for index in range(1, 119))]


def abundance_score(text: str) -> float | None:
    if "-" in text:
        return None
    values = [
        float(value)
        for value in re.findall(r"\d+(?:\.\d+)?", text.replace(" ", ""))
    ]
    if not values:
        return None
    return sum(values[:2]) / min(2, len(values)) if text.startswith("[") else values[0]


def natural_isotopes(path: Path) -> list[int | None]:
    source = path.read_text(encoding="utf-8")
    # The published page omits a few opening <tr> tags before rowspan rows.
    source = re.sub(r"(</tr>\s*)(<td rowspan)", r"\1<tr>\2", source)
    best: dict[int, tuple[int, float]] = {}
    current: int | None = None
    for match in re.finditer(r"<tr[^>]*>(.*?)</tr>", source, re.I | re.S):
        row = match.group(0)
        fields = table_cells(match.group(1))
        if not fields:
            continue
        if "childRow" not in row and fields[0].isdigit() and len(fields) >= 5:
            current = int(fields[0])
            mass_text, abundance_text = fields[3], fields[4]
        elif "childRow" in row and current is not None and len(fields) >= 2:
            mass_text, abundance_text = fields[0], fields[1]
        else:
            continue
        mass_match = re.search(r"\d+", mass_text)
        score = abundance_score(abundance_text)
        if mass_match is None or score is None:
            continue
        candidate = (int(mass_match.group()), score)
        if score > best.get(current, (0, -1.0))[1]:
            best[current] = candidate
    if len(best) != 84:
        raise SystemExit(f"expected 84 naturally occurring elements, found {len(best)}")
    return [None, *(best.get(index, (None, 0.0))[0] for index in range(1, 119))]


def isotope_masses(path: Path) -> list[tuple[int, int, float]]:
    rows: list[tuple[int, int, float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        try:
            atomic_number = int(line[9:14])
            mass_number = int(line[14:19])
            whole = int(line[106:109])
            micro_u = float(line[110:123].replace("#", "."))
        except (IndexError, ValueError):
            continue
        if 1 <= atomic_number <= 118:
            rows.append((atomic_number, mass_number, whole + micro_u / 1_000_000.0))
    rows.sort(key=lambda row: row[:2])
    if len(rows) != 3557 or len({row[:2] for row in rows}) != len(rows):
        raise SystemExit(f"expected 3557 unique AME rows, found {len(rows)}")
    return rows


def formatted_array(values: list[object], render, per_line: int = 1) -> str:
    rendered = [render(value) for value in values]
    return "\n".join(
        "    " + ", ".join(rendered[index : index + per_line]) + ","
        for index in range(0, len(rendered), per_line)
    )


def render_data(
    weights: list[float | None],
    abundant: list[int | None],
    masses: list[tuple[int, int, float]],
) -> str:
    option_float = lambda value: "None" if value is None else f"Some({value!r})"
    option_int = lambda value: "None" if value is None else f"Some({value})"
    mass_rows = "\n".join(
        f"    ({atomic_number}, {mass_number}, {mass!r}),"
        for atomic_number, mass_number, mass in masses
    )
    return f'''// This file is mechanically generated from the pinned sources below.
// Regenerate with: python validation/reference/atomic_data/generate.py
// CIAAW Abridged Standard Atomic Weights 2024:
//   https://ciaaw.org/abridged-atomic-weights.htm
//   SHA-256 f9e9554471749c55a624aec55151922470a7f4104c62811eb194fed9731b907d
// CIAAW Isotopic Compositions of the Elements 2024:
//   https://ciaaw.org/isotopic-abundances.htm
//   SHA-256 1ab770a48bb539d26ad2992a94b4c60d51743058b67f106796b754cd5035423b
// AME 2020 mass_1.mas20 (3 March 2021 release):
//   https://amdc.impcas.ac.cn/masstables/Ame2020/mass_1.mas20
//   SHA-256 e8599c6d7f724fac91934e59f1b9de8fb8f63e820f4b39456b790665ed2a3307
// AME rows marked as estimates are retained as evaluated atomic-mass values;
// no integer mass-number fallback is synthesized.
// 2022 CODATA electron rest mass in unified atomic mass units:
//   https://physics.nist.gov/cgi-bin/cuu/Value?meu

use crate::core::Element;

pub(super) const ELECTRON_MASS_DA: f64 = 5.485_799_090_441e-4;

const STANDARD_ATOMIC_WEIGHTS: [Option<f64>; 119] = [
{formatted_array(weights, option_float)}
];

const MOST_ABUNDANT_ISOTOPES: [Option<u16>; 119] = [
{formatted_array(abundant, option_int)}
];

const ISOTOPE_MASSES: &[(u8, u16, f64)] = &[
{mass_rows}
];

pub(super) fn standard_atomic_weight(element: Element) -> Option<f64> {{
    STANDARD_ATOMIC_WEIGHTS[usize::from(element.atomic_number())]
}}

pub(super) fn most_abundant_isotope(element: Element) -> Option<u16> {{
    MOST_ABUNDANT_ISOTOPES[usize::from(element.atomic_number())]
}}

pub(super) fn exact_isotope_mass(element: Element, mass_number: u16) -> Option<f64> {{
    let key = (element.atomic_number(), mass_number);
    ISOTOPE_MASSES
        .binary_search_by_key(&key, |&(atomic_number, mass_number, _)| {{
            (atomic_number, mass_number)
        }})
        .ok()
        .map(|index| ISOTOPE_MASSES[index].2)
}}

#[cfg(test)]
mod tests {{
    use super::*;

    #[test]
    fn generated_tables_cover_pinned_reference_rows() {{
        assert_eq!(ISOTOPE_MASSES.len(), 3557);
        assert_eq!(
            exact_isotope_mass(Element::from_symbol("C").unwrap(), 13),
            Some(13.003_354_835_34)
        );
        assert_eq!(
            most_abundant_isotope(Element::from_symbol("U").unwrap()),
            Some(238)
        );
        assert_eq!(
            standard_atomic_weight(Element::from_symbol("Tc").unwrap()),
            None
        );
    }}
}}
'''


if __name__ == "__main__":
    raise SystemExit(main())
