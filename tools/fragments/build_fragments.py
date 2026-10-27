"""Builds kekule's fragment library SDF from PubChem 3D records.

Usage: python -I build_fragments.py <pubchem-dir> <output.sdf>
"""
import os, sys

# name, PubChem record file, attachment atom (1-based), leaving atom (1-based,
# None for the attachment's least hindered hydrogen), kind
FRAGMENTS = [
    ("benzene", "benzene", 1, None, "ring"),
    ("pyridine", "pyridine", 2, None, "ring"),
    ("pyrimidine", "pyrimidine", 6, None, "ring"),
    ("pyrazine", "pyrazine", 3, None, "ring"),
    ("pyrrole", "pyrrole", 1, None, "ring"),
    ("furan", "furan", 4, None, "ring"),
    ("thiophene", "thiophene", 4, None, "ring"),
    ("imidazole", "imidazole", 1, None, "ring"),
    ("pyrazole", "pyrazole", 1, None, "ring"),
    ("oxazole", "oxazole", 5, None, "ring"),
    ("thiazole", "thiazole", 5, None, "ring"),
    ("cyclopropane", "cyclopropane", 1, None, "ring"),
    ("cyclobutane", "cyclobutane", 1, None, "ring"),
    ("cyclopentane", "cyclopentane", 1, None, "ring"),
    ("cyclohexane", "cyclohexane", 1, None, "ring"),
    ("piperidine", "piperidine", 1, None, "ring"),
    ("piperazine", "piperazine", 1, None, "ring"),
    ("morpholine", "morpholine", 2, None, "ring"),
    ("tetrahydrofuran", "tetrahydrofuran", 4, None, "ring"),
    ("naphthalene", "naphthalene", 7, None, "ring"),
    ("indole", "indole", 4, None, "ring"),
    ("quinoline", "quinoline", 5, None, "ring"),
    ("benzimidazole", "benzimidazole", 8, None, "ring"),
    ("methyl", "methane", 1, None, "group"),
    ("ethyl", "ethane", 1, None, "group"),
    ("isopropyl", "propane", 1, None, "group"),
    ("tert-butyl", "isobutane", 1, None, "group"),
    ("vinyl", "ethylene", 1, None, "group"),
    ("ethynyl", "acetylene", 1, None, "group"),
    ("phenyl", "benzene", 1, None, "group"),
    ("benzyl", "toluene", 4, None, "group"),
    ("hydroxyl", "water", 1, None, "group"),
    ("methoxy", "methanol", 1, None, "group"),
    ("amino", "ammonia", 1, None, "group"),
    ("methylamino", "methylamine", 1, None, "group"),
    ("dimethylamino", "dimethylamine", 1, None, "group"),
    ("formyl", "formaldehyde", 2, None, "group"),
    ("acetyl", "acetaldehyde", 3, None, "group"),
    ("carboxyl", "formic_acid", 3, None, "group"),
    ("methoxycarbonyl", "methyl_formate", 4, None, "group"),
    ("carbamoyl", "formamide", 3, None, "group"),
    ("acetamido", "acetamide", 2, None, "group"),
    ("nitro", "nitromethane", 3, 4, "group"),
    ("cyano", "hydrogen_cyanide", 2, None, "group"),
    ("trifluoromethyl", "fluoroform", 4, None, "group"),
    ("methylsulfonyl", "dimethyl_sulfone", 1, 5, "group"),
    ("sulfonamide", "methanesulfonamide", 1, 5, "group"),
    ("thiol", "hydrogen_sulfide", 1, None, "group"),
]


# Ring atoms of each ring parent, as PubChem atom numbers listed in the
# ring's conventional numbering (IUPAC locants: pyridine N1, C2 ... C6;
# naphthalene C1 ... C4, C4a, C5 ... C8, C8a).
RING_ORDER = {
    "benzene": [1, 2, 4, 6, 5, 3],
    "pyridine": [1, 5, 3, 2, 4, 6],
    "pyrimidine": [1, 6, 2, 5, 3, 4],
    "pyrazine": [1, 3, 4, 2, 6, 5],
    "pyrrole": [1, 2, 4, 5, 3],
    "furan": [1, 4, 2, 3, 5],
    "thiophene": [1, 4, 2, 3, 5],
    "imidazole": [1, 4, 2, 5, 3],
    "pyrazole": [1, 2, 5, 4, 3],
    "oxazole": [1, 5, 2, 4, 3],
    "thiazole": [1, 5, 2, 4, 3],
    "cyclopropane": [1, 2, 3],
    "cyclobutane": [1, 2, 4, 3],
    "cyclopentane": [1, 2, 4, 5, 3],
    "cyclohexane": [1, 2, 4, 6, 5, 3],
    "piperidine": [1, 5, 3, 2, 4, 6],
    "piperazine": [1, 4, 3, 2, 5, 6],
    "morpholine": [1, 5, 3, 2, 4, 6],
    "tetrahydrofuran": [1, 4, 2, 3, 5],
    "naphthalene": [3, 7, 8, 4, 2, 6, 10, 9, 5, 1],
    "indole": [1, 7, 4, 3, 6, 9, 8, 5, 2],
    "quinoline": [1, 10, 9, 5, 2, 4, 7, 8, 6, 3],
    "benzimidazole": [1, 8, 2, 4, 6, 9, 7, 5, 3],
}


# The locants of fused ring systems, in RING_ORDER's order; a single ring's
# atoms are numbered 1 to n.
FUSED_LOCANTS = {
    "naphthalene": ["1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a"],
    "indole": ["1", "2", "3", "3a", "4", "5", "6", "7", "7a"],
    "quinoline": ["1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a"],
    "benzimidazole": ["1", "2", "3", "3a", "4", "5", "6", "7", "7a"],
}


def locants(record):
    """The locant of each ring atom of a ring parent, in RING_ORDER's order."""
    ring = RING_ORDER[record]
    return FUSED_LOCANTS.get(record, [str(i) for i in range(1, len(ring) + 1)])


def atom_order(record, attachment, elements, bonds):
    """The new order of the record's atoms, as old 1-based numbers: the ring
    in its numbering (or the attachment atom first, then the other heavy
    atoms), then each heavy atom's hydrogens in that order."""
    heavy = [i for i, element in enumerate(elements, 1) if element != "H"]
    first = RING_ORDER.get(record, [attachment])
    order = first + [i for i in heavy if i not in first]
    for atom in list(order):
        for a, b, _ in bonds:
            other = b if a == atom else a if b == atom else None
            if other is not None and elements[other - 1] == "H" and other not in order:
                order.append(other)
    order += [i for i in range(1, len(elements) + 1) if i not in order]
    return order


def reordered(lines, n, m, order):
    """The connection table lines with atoms in `order`."""
    new = {old: index for index, old in enumerate(order, 1)}
    atoms = [lines[4 + old - 1] for old in order]
    bonds = []
    for line in lines[4 + n:4 + n + m]:
        a, b = new[int(line[0:3])], new[int(line[3:6])]
        bonds.append(f"{a:>3}{b:>3}{line[6:]}")
    properties = []
    for line in lines[4 + n + m:]:
        if line.startswith("M  CHG"):
            fields = line[6:].split()
            pairs = [(new[int(fields[i])], int(fields[i + 1])) for i in range(1, len(fields), 2)]
            line = "M  CHG" + f"{len(pairs):>3}" + "".join(f" {a:>3} {c:>3}" for a, c in pairs)
        properties.append(line)
    return [lines[3]] + atoms + bonds + properties, new


def main(source, output):
    out = []
    for name, record, attachment, leaving, kind in FRAGMENTS:
        text = open(os.path.join(source, record + ".sdf"), encoding="utf-8").read()
        lines = text.split("\n")
        cid = lines[0].strip()
        end = lines.index("M  END")
        counts = lines[3]
        n, m = int(counts[0:3]), int(counts[3:6])
        elements = [line[31:34].strip() for line in lines[4:4 + n]]
        if leaving is None:
            hydrogens = []
            for line in lines[4 + n:4 + n + m]:
                a, b = int(line[0:3]), int(line[3:6])
                if a == attachment and elements[b - 1] == "H":
                    hydrogens.append(b)
                if b == attachment and elements[a - 1] == "H":
                    hydrogens.append(a)
            # The least hindered hydrogen leaves: the one farthest from the
            # atoms not bonded to the attachment atom, as an equatorial
            # hydrogen of a ring or one perpendicular to an aryl plane.
            coordinates = [
                tuple(float(line[i:i + 10]) for i in (0, 10, 20)) for line in lines[4:4 + n]
            ]
            bonded = {attachment}
            for line in lines[4 + n:4 + n + m]:
                a, b = int(line[0:3]), int(line[3:6])
                if attachment in (a, b):
                    bonded.add(a if b == attachment else b)
            others = [i for i in range(1, n + 1) if i not in bonded]

            def clearance(h):
                if not others:
                    return 0.0
                return round(min(
                    sum((x - y) ** 2 for x, y in zip(coordinates[h - 1], coordinates[o - 1]))
                    for o in others
                ), 4)

            leaving = max(sorted(hydrogens), key=clearance)
        bonds = [
            (int(line[0:3]), int(line[3:6]), int(line[6:9])) for line in lines[4 + n:4 + n + m]
        ]
        order = atom_order(record, attachment, elements, bonds)
        table, new = reordered(lines[:end + 1], n, m, order)
        attachment, leaving = new[attachment], new[leaving]
        block = [name, lines[1], f"PubChem CID {cid}, {record.replace('_', ' ')}"] + table
        block += [
            "> <PUBCHEM_COMPOUND_CID>", cid, "",
            "> <kind>", kind, "",
            "> <attachment>", str(attachment), "",
            "> <leaving>", str(leaving), "",
        ]
        if kind == "ring":
            block += ["> <locants>", " ".join(locants(record)), ""]
        block += ["$$$$"]
        out.extend(block)
    open(output, "w", encoding="utf-8", newline="\n").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
