"""Compare two fact lists and describe every difference as a signature.

A fact is `[key..., value]`. Facts with equal keys are compared; a key on
one side only is a presence difference. Values are never rounded: floats
and float lists differ only beyond the tolerance justified below.

A signature is `path|context|transition`:
- `path` is the key with source positions removed (`atom.hydrogens`);
- `context` comes from the reference parse of the same input, so a change
  in Kekule never shifts it (`N|ar|+1`, `C-O|2`, `size=6`, `query=...`);
- `transition` is `reference→kekule` for scalar values, `fewer`/`more` for
  counts and `≠` for lists and floats.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass

MISSING = object()

# Absolute tolerances by path. Each is half the step of the precision at which
# the compared value is printed, plus a margin for binary floating point.
TOLERANCES = {
    # Molfile coordinates are written with four decimals (0.0001 Å).
    "atom.xyz": 5e-5 + 1e-9,
    # mmCIF atom sites: coordinates with three decimals, occupancy and B with two.
    "site.xyz": 5e-4 + 1e-9,
    "site.occupancy": 5e-3 + 1e-9,
    "site.b": 5e-3 + 1e-9,
    # mkdssp prints angles and energies with one decimal, TCO with three.
    "dssp.phi": 0.05 + 1e-9,
    "dssp.psi": 0.05 + 1e-9,
    "dssp.kappa": 0.05 + 1e-9,
    "dssp.alpha": 0.05 + 1e-9,
    "dssp.tco": 5e-4 + 1e-9,
    "dssp.acceptor_1_energy": 0.05 + 1e-9,
    "dssp.acceptor_2_energy": 0.05 + 1e-9,
    "dssp.donor_1_energy": 0.05 + 1e-9,
    "dssp.donor_2_energy": 0.05 + 1e-9,
    # Both sides compute RMSD in double precision from the same coordinates.
    "model.rmsd": 1e-6,
}
# Angles compare by circular distance.
ANGLES = {"dssp.phi", "dssp.psi", "dssp.kappa", "dssp.alpha"}
# Record-valued facts compare field by field.
FIELDS = {
    "site": ("label_asym", "auth_asym", "label_seq", "auth_seq", "insertion", "label_comp", "label_atom",
             "element", "altloc", "occupancy", "b", "xyz"),
    "dssp": ("ss", "phi", "psi", "kappa", "alpha", "tco", "acceptor_1", "acceptor_1_energy", "acceptor_2",
             "acceptor_2_energy", "donor_1", "donor_1_energy", "donor_2", "donor_2_energy"),
}


def key_of(fact: list) -> tuple:
    return tuple(
        json.dumps(part, separators=(",", ":")) if isinstance(part, (list, dict)) else part for part in fact[:-1]
    )


def index(facts: list) -> dict:
    values = {}
    for fact in facts:
        values[key_of(fact)] = fact[-1]
    return values


@dataclass(frozen=True)
class Difference:
    key: tuple
    path: str
    context: str
    reference: object
    kekule: object

    @property
    def transition(self) -> str:
        if self.path.endswith(".count") and isinstance(self.reference, int) and isinstance(self.kekule, int):
            return "fewer" if self.kekule < self.reference else "more"
        if self.path in PARTNERS:
            return f"{partner(self.reference)}→{partner(self.kekule)}"
        left, right = label(self.reference), label(self.kekule)
        if left is None or right is None:
            if self.reference is MISSING:
                return "missing→present"
            if self.kekule is MISSING:
                return "present→missing"
            return "≠"
        return f"{left}→{right}"

    @property
    def signature(self) -> str:
        return f"{self.path}|{self.context}|{self.transition}"

    def delta(self) -> float | None:
        """Largest absolute numeric difference, for numeric signatures."""
        pairs = zip(as_floats(self.reference), as_floats(self.kekule))
        values = [abs(a - b) for a, b in pairs]
        return max(values) if values else None

    def to_json(self) -> dict:
        return {
            "key": list(self.key),
            "signature": self.signature,
            "reference": printable(self.reference),
            "kekule": printable(self.kekule),
        }


# DSSP hydrogen-bond partners, `chain:number`, compare by kind of residue.
PARTNERS = {"dssp.acceptor_1", "dssp.acceptor_2", "dssp.donor_1", "dssp.donor_2"}


def partner(value) -> str:
    if value is MISSING or value is None:
        return "none"
    return "unnumbered" if str(value).endswith(":?") else "residue"


def label(value) -> str | None:
    if value is MISSING:
        return "missing"
    if value is None:
        return "null"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str) and len(value) <= 40:
        return value
    return None


def printable(value):
    return "<missing>" if value is MISSING else value


def as_floats(value) -> list[float]:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return [float(value)]
    if isinstance(value, list) and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in value):
        return [float(v) for v in value]
    return []


def close(path: str, reference, kekule) -> bool:
    tolerance = TOLERANCES.get(path)
    if tolerance is None:
        return False
    left, right = as_floats(reference), as_floats(kekule)

    def distance(a, b):
        d = abs(a - b) % 360.0 if path in ANGLES else abs(a - b)
        return min(d, 360.0 - d) if path in ANGLES else d

    return bool(left) and len(left) == len(right) and all(
        math.isfinite(a) and math.isfinite(b) and distance(a, b) <= tolerance for a, b in zip(left, right)
    )


def describe(key: tuple, context: dict, options: dict) -> tuple[str, str]:
    """The path and context of a fact key."""
    kind = key[0]
    atoms = context.get("atoms", {})
    bonds = context.get("bonds", {})

    def atom(i) -> str:
        return atoms.get(str(i), "?")

    def bond(a, b) -> str:
        return bonds.get(f"{a}-{b}") or bonds.get(f"{b}-{a}") or "?"

    def focus(text: str) -> str:
        positions = json.loads(text)
        return atom(positions[0]) if len(positions) == 1 else bond(*positions[:2])

    if kind == "atom":
        return f"atom.{key[2]}", atom(key[1])
    if kind == "bond" and isinstance(key[1], str):
        return f"bond.{key[3]}", ""
    if kind == "bond":
        return f"bond.{key[3]}", bond(key[1], key[2])
    if kind in ("item", "column"):
        return f"{kind}.{key[3]}", key[2].split(".")[0]
    if kind == "residue":
        return f"residue.{key[-1]}", str(key[1]).rsplit("/", 1)[-1]
    if kind in ("site", "dssp"):
        return kind, ""
    if kind in ("molecule", "model", "alpha", "block"):
        return f"{kind}.{key[-1]}", ""
    if kind in ("stereo", "candidate"):
        return f"{kind}.{key[1]}.{key[3]}", focus(key[2])
    if kind == "cip":
        return f"cip.{key[1]}.{key[3]}", focus(key[2])
    if kind == "rings":
        return "rings.count", f"size={key[1]}"
    if kind == "formula":
        if key[2] == "count":
            return "formula.count", f"{key[3]}{'' if key[4] is None else key[4]}"
        return "formula.charge", ""
    if kind == "resonance":
        return f"resonance.{key[1]}.{key[3]}", ""
    if kind == "query":
        queries = options.get("queries", [])
        query = queries[key[1]] if key[1] < len(queries) else "?"
        return f"query.{key[2]}", f"query={query}"
    return ".".join(str(part) for part in key if not isinstance(part, int) or isinstance(part, bool)), ""


def electron_distribution(contributor: str) -> str:
    """A contributor `charges|bonds` reduced to each atom's charge and
    multiple-bond excess. Two contributors with equal reductions differ only
    by alternating cycles, that is, by their Kekulé form; any move of charge
    or of a lone pair changes it. Bonds without an integer order are kept."""
    charges, bonds = contributor.split("|")
    atoms: dict[int, list] = {}
    kept = []
    for item in filter(None, charges.split(",")):
        atom, charge = item.split(":")
        atoms.setdefault(int(atom), [0, 0])[0] = int(charge)
    for item in filter(None, bonds.split(",")):
        pair, order = item.split(":")
        if not order.isdigit():
            kept.append(item)
            continue
        for atom in pair.split("-"):
            atoms.setdefault(int(atom), [0, 0])[1] += int(order) - 1
    excess = ",".join(f"{atom}:{charge}:{extra}" for atom, (charge, extra) in sorted(atoms.items()))
    return f"{excess}|{','.join(sorted(kept))}"


def normalize(task_id: str, kekule: dict, reference: dict, context: dict) -> None:
    """Apply the comparison rules that are not plain equality, in place."""
    if task_id == "bio.hierarchy":
        # Rows Kekule did not keep reflect its documented altloc policy, which
        # unit tests own; only kept rows are compared.
        for key in [k for k in reference if k[0] == "site" and k not in kekule]:
            reference.pop(key)
        return
    if task_id == "bio.connectivity":
        # Altloc policies may keep different rows; compare bonds and molecules
        # over the sites both sides kept.
        def kept(side):
            return {site for key, value in side.items() if key[0] == "molecule" for site in value}

        # The molecule partition follows from the bonds, so it only supplies
        # the kept sites and is not compared again.
        common = kept(kekule) & kept(reference)
        for side in (kekule, reference):
            for key in list(side):
                if key[0] == "molecule" or (key[0] == "bond" and not (key[1] in common and key[2] in common)):
                    side.pop(key)
        return
    if task_id == "small.resonance":
        # Contributors that differ only by a Kekulé swap are one structure;
        # which Kekulé form a toolkit lists is not chemistry.
        for side in (kekule, reference):
            for key in [k for k in side if k[0] == "resonance" and k[1] == "contributors" and k[3] == "set"]:
                side[key] = sorted({electron_distribution(c) for c in side[key]})
                side[key[:3] + ("count",)] = len(side[key])
        return
    if task_id != "small.parse" and not task_id.startswith("small.write."):
        return
    # Kekule localizes bonds the source wrote as aromatic; a specific Kekulé
    # assignment is not chemistry, so only the bond's presence is compared.
    for a, b in context.get("source_aromatic", []):
        key = ("bond", a, b, "order")
        if reference.get(key) == "aromatic" and kekule.get(key) in (1, 2):
            kekule[key] = "aromatic"
    # Total hydrogens are compared only where both toolkits define them.
    for key in [k for k in reference if k[0] == "atom" and k[2] == "hydrogens"]:
        if reference[key] is None or kekule.get(key, MISSING) is None:
            reference.pop(key)
            kekule.pop(key, None)


def diff(task_id: str, kekule_facts: list, reference_facts: list, context: dict, options: dict) -> list[Difference]:
    kekule, reference = index(kekule_facts), index(reference_facts)
    normalize(task_id, kekule, reference, context)
    differences = []
    for key in sorted(set(kekule) | set(reference), key=lambda k: json.dumps(k, default=str)):
        left, right = reference.get(key, MISSING), kekule.get(key, MISSING)
        if left is not MISSING and right is not MISSING and left == right and type(left) is type(right):
            continue
        names = FIELDS.get(key[0])
        if names and isinstance(left, list) and isinstance(right, list) and len(left) == len(right) == len(names):
            where = str(left[5]) if key[0] == "site" else ""
            fields = dict(zip(names, zip(left, right)))
            for name, (a, b) in fields.items():
                path = f"{key[0]}.{name}"
                if (a == b and type(a) is type(b)) or close(path, a, b):
                    continue
                # A hydrogen-bond energy is comparable only for the same partner.
                partner = fields.get(name.removesuffix("_energy")) if name.endswith("_energy") else None
                if partner is not None and partner[0] != partner[1]:
                    continue
                differences.append(Difference((*key, name), path, where, a, b))
            continue
        path, where = describe(key, context, options)
        if close(path, left, right):
            continue
        differences.append(Difference(key, path, where, left, right))
    return differences
