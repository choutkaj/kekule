"""Extract native normalization edits from the pinned NAGL reaction definitions."""
import ast
import json
import re
from pathlib import Path
from rdkit import Chem
from common import HERE, ROOT, verify_sources

verify_sources()
tree = ast.parse((HERE / "upstream/openff.py").read_text())
fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "normalize_molecule")
assignment = next(n for n in fn.body if isinstance(n, ast.Assign) and
                  any(isinstance(t, ast.Name) and t.id == "normalizations" for t in n.targets))
rules = []
for reaction in ast.literal_eval(assignment.value):
    reactant, product = reaction.split(">>")
    charges = []
    for content in re.findall(r"\[([^\]]+)\]", product):
        tag = int(content.rsplit(":", 1)[1])
        charge = re.search(r"([+-])(\d*)", content)
        if charge:
            magnitude = int(charge[2]) if charge[2] else 1
            charges.append([tag, magnitude * (1 if charge[1] == "+" else -1)])
    molecule = Chem.MolFromSmarts(product)
    bonds = [[b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum(),
              None if b.GetIsAromatic() else int(b.GetBondTypeAsDouble())] for b in molecule.GetBonds()]
    rules.append(dict(reaction=reaction, query=reactant, charges=charges, bonds=bonds))
(ROOT / "crates/kekule-openff/data/normalizations.json").write_text(
    json.dumps(rules, indent=2) + "\n", encoding="utf-8", newline="\n")
