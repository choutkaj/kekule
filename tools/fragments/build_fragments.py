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
        block = [name, lines[1], f"PubChem CID {cid}, {record.replace('_', ' ')}"] + lines[3:end + 1]
        block += [
            "> <PUBCHEM_COMPOUND_CID>", cid, "",
            "> <kind>", kind, "",
            "> <attachment>", str(attachment), "",
            "> <leaving>", str(leaving), "",
            "$$$$",
        ]
        out.extend(block)
    open(output, "w", encoding="utf-8", newline="\n").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
