"""Independent observations from RDKit; no Kekule algorithms or normalization."""
from rdkit import Chem, rdBase


def record(record, feature):
    if rdBase.rdkitVersion != '2026.03.3':
        raise ValueError('resonance reference requires RDKit 2026.03.3')
    mol = Chem.Mol(record['mol'])
    Chem.SanitizeMol(mol)
    bonds = sorted((tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))), b.GetIdx()) for b in mol.GetBonds())
    result = dict(record_index=record['record_index'], status='ok', title=record['title'])
    if feature == 'algo.conjugation.rdkit-like':
        result['bonds'] = [dict(atoms=list(ends), conjugated=mol.GetBondWithIdx(i).GetIsConjugated()) for ends, i in bonds]
        return result
    supplier = Chem.ResonanceMolSupplier(mol)
    groups = []
    for g in range(supplier.GetNumConjGrps()):
        atoms = [a.GetIdx() for a in mol.GetAtoms() if supplier.GetAtomConjGrpIdx(a.GetIdx()) == g]
        edges = [list(ends) for ends, i in bonds if supplier.GetBondConjGrpIdx(i) == g]
        groups.append(dict(atoms=atoms, bonds=edges))
    result['groups'] = sorted(groups, key=lambda g: (g['atoms'], g['bonds']))
    if feature == 'algo.resonance.groups':
        return result
    result['bonds'] = [list(ends) for ends, _ in bonds]
    profiles = []
    for flags in range(32):
        supplier = Chem.ResonanceMolSupplier(mol, flags, 1000)
        structures = []
        for form in supplier:
            charges = [a.GetFormalCharge() for a in form.GetAtoms()]
            orders = [5 if form.GetBondWithIdx(i).GetBondType() == Chem.BondType.DATIVE else int(form.GetBondWithIdx(i).GetBondTypeAsDouble()) for _, i in bonds]
            structures.append((charges, orders))
        structures.sort()
        profiles.append(dict(flags=flags, max_structures=1000, count=len(structures), structures=[dict(charges=q, orders=bs) for q, bs in structures]))
    result['profiles'] = profiles
    return result
