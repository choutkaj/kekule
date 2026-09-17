"""Independent, indexed molecular observations and external writer validation.

All chemistry is computed by RDKit. Atoms keep their source correspondence;
components and stereo carriers are never replaced with local chemical hashes.
"""
from __future__ import annotations

import json
import io
import re
from typing import Any
from rdkit import Chem


def ordered(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def parity(values):
    return sum(a > b for i, a in enumerate(values) for b in values[i+1:]) % 2


def canonical_carrier(mol, center, other):
    return min((a.GetIdx() for a in mol.GetAtomWithIdx(center).GetNeighbors() if a.GetIdx() != other), default=-1)


def stereo(mol):
    elements = []
    for atom in mol.GetAtoms():
        tag = atom.GetChiralTag()
        if tag == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        if tag not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            raise ValueError(f'RDKit stereo type has no common representation: {tag}')
        carriers = [a.GetIdx() for a in atom.GetNeighbors()]
        if len(carriers) == 3:
            carriers.append(-1 if atom.GetTotalNumHs() else -2)
        if len(carriers) != 4:
            raise ValueError('tetrahedral carrier count is not four')
        elements.append({'type': 'tetrahedral', 'focus': [atom.GetIdx()], 'carriers': sorted(carriers),
                         'parity': int(tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW) ^ parity(carriers)})
    for bond in mol.GetBonds():
        tag = bond.GetStereo()
        if tag == Chem.BondStereo.STEREONONE:
            continue
        ends = [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
        canonical = [canonical_carrier(mol, ends[0], ends[1]), canonical_carrier(mol, ends[1], ends[0])]
        if tag in (Chem.BondStereo.STEREOATROPCW, Chem.BondStereo.STEREOATROPCCW):
            # RDKit Atropisomers.cpp defines orientation using the lowest
            # numbered neighbor on each end, independent of bond insertion order.
            # Positive (axis . (left_reference x right_reference)) is CW.
            p = int(tag == Chem.BondStereo.STEREOATROPCCW)
            if ends[0] > ends[1]:
                ends.reverse(); canonical.reverse()
            elements.append({'type': 'axis', 'focus': ends, 'carriers': canonical, 'parity': p})
            continue
        if tag == Chem.BondStereo.STEREOANY:
            p = None
        elif tag in (Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS, Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
            selected = list(bond.GetStereoAtoms())
            if len(selected) != 2:
                raise ValueError('double bond stereo missing carriers')
            p = int(tag in (Chem.BondStereo.STEREOTRANS, Chem.BondStereo.STEREOE)) ^ int(selected[0] != canonical[0]) ^ int(selected[1] != canonical[1])
        else:
            raise ValueError(f'RDKit stereo type has no common representation: {tag}')
        if ends[0] > ends[1]:
            ends.reverse(); canonical.reverse()
        elements.append({'type': 'double_bond', 'focus': ends, 'carriers': canonical, 'parity': p})
    groups = []
    kinds = {Chem.StereoGroupType.STEREO_ABSOLUTE: 'absolute', Chem.StereoGroupType.STEREO_AND: 'and', Chem.StereoGroupType.STEREO_OR: 'or'}
    for group in mol.GetStereoGroups():
        members = [{'type': 'tetrahedral', 'focus': [a.GetIdx()]} for a in group.GetAtoms()]
        members.extend({'type': 'axis', 'focus': sorted([b.GetBeginAtomIdx(), b.GetEndAtomIdx()])} for b in group.GetBonds())
        groups.append({'kind': kinds[group.GetGroupType()], 'members': sorted(members, key=ordered)})
    return sorted(elements, key=ordered), sorted(groups, key=ordered)


def graph(mol):
    atoms = []
    for atom in mol.GetAtoms():
        value = atom_json(atom)
        value['implicit_hydrogens'] = atom.GetNumImplicitHs()
        value['explicit_valence'] = atom.GetValence(Chem.rdchem.ValenceType.EXPLICIT)
        value['no_implicit_hydrogens'] = atom.GetNoImplicit()
        if mol.GetNumConformers():
            point = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            value['coord'] = [point.x, point.y, point.z]
        atoms.append(value)
    bonds = []
    for bond in mol.GetBonds():
        ends = [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
        if bond.GetBondType() != Chem.BondType.DATIVE:
            ends.sort()
        bonds.append({'begin_atom_index': ends[0], 'end_atom_index': ends[1],
                      'bond_type': str(bond.GetBondType()), 'is_aromatic': bond.GetIsAromatic()})
    elements, groups = stereo(mol)
    return {'atom_count': mol.GetNumAtoms(), 'bond_count': mol.GetNumBonds(),
            'atoms': atoms, 'bonds': bonds, 'stereo': elements, 'groups': groups}


def sdf_record(block):
    # RDKit stores SDF properties by name. Repeated names cannot be represented
    # faithfully by that API, so retain this case as a reference failure.
    headers = re.findall(r'^>\s*.*?<([^>]+)>', block, re.MULTILINE)
    if len(headers) != len(set(headers)):
        raise ValueError('RDKit cannot preserve duplicate SDF property names')
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(block.encode('utf-8')),
        sanitize=False, removeHs=False, strictParsing=True)
    molecules = list(supplier)
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError('RDKit SDF parse failed')
    mol = molecules[0]
    fields = [{'name': name, 'value': mol.GetProp(name)} for name in mol.GetPropNames()]
    return mol, fields


def records(source):
    if source.suffix.lower() in ('.smi', '.smiles', '.txt'):
        params = Chem.SmilesParserParams()
        params.sanitize = False
        params.removeHs = False
        values = []
        for line in source.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            # RDKit owns the grammar for optional CX extensions and names.
            mol = Chem.MolFromSmiles(line.strip(), params)
            if mol is None:
                raise ValueError('RDKit SMILES parse failed')
            mol.UpdatePropertyCache(strict=False)
            Chem.SetBondStereoFromDirections(mol)
            values.append((mol.GetProp('_Name') if mol.HasProp('_Name') else '', mol, []))
        return values
    blocks = read_sdf_blocks(source) if source.suffix.lower() == '.sdf' else [source.read_text()]
    result = []
    for block in blocks:
        if source.suffix.lower() == '.sdf':
            mol, fields = sdf_record(block)
        else:
            mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False, strictParsing=True)
            fields = []
        if mol is None:
            raise ValueError('RDKit MOL parse failed')
        mol.UpdatePropertyCache(strict=False)
        result.append((mol.GetProp('_Name'), mol, fields))
    return result


def evaluate(feature, source):
    outputs = []
    for index, (title, mol, fields) in enumerate(records(source)):
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        if feature == 'stereo.perception':
            if mol.GetNumConformers():
                Chem.AssignStereochemistryFrom3D(mol)
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        components = []
        for fragment in fragments:
            value = graph(fragment)
            if feature == 'stereo.perception':
                candidates = []
                for info in Chem.FindPotentialStereo(fragment):
                    kind = str(info.type)
                    atom_types = {'Atom_Tetrahedral':'tetrahedral', 'Atom_SquarePlanar':'square_planar',
                                  'Atom_TrigonalBipyramidal':'trigonal_bipyramidal', 'Atom_Octahedral':'octahedral'}
                    bond_types = {'Bond_Double':'double_bond', 'Bond_Atropisomer':'axis', 'Bond_Cumulene_Even':'cumulene_even'}
                    if kind in atom_types:
                        candidates.append({'type': atom_types[kind], 'focus': [info.centeredOn]})
                    elif kind in bond_types:
                        bond = fragment.GetBondWithIdx(info.centeredOn)
                        candidates.append({'type': bond_types[kind], 'focus': sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])})
                    else:
                        raise ValueError(f'no common candidate representation for {info.type}')
                value['candidates'] = sorted(candidates, key=ordered)
            components.append(value)
        outputs.append({'record_index': index, 'status': 'ok', 'title': title, 'components': components, 'properties': fields})
    return {'records': outputs}


def writer_value(feature, source):
    if feature.startswith('io.smiles.'):
        outputs = []
        for index, (title, mol, fields) in enumerate(records(source)):
            Chem.SanitizeMol(mol)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            # RDKit's complete canonical isomeric CXSMILES encodes connectivity,
            # isotopes, radicals, maps, stereo groups and stereo placement.
            params = Chem.SmilesWriteParams()
            params.canonical = True
            params.doIsomericSmiles = True
            outputs.append({'record_index': index, 'status': 'ok', 'title': title,
                            'identity': Chem.MolToCXSmiles(mol, params, Chem.CXSmilesFields.CX_ALL)})
        return {'records': outputs}
    # Compare all parsed molecular properties and coordinates. Do not synthesize
    # zero coordinates or copy a title from the input into the writer result.
    value = evaluate('molecule', source)
    if feature.startswith('io.mol.'):
        for record in value['records']:
            record['properties'] = []  # The MOL format has no SDF data section.
    return value


def read_written(feature, value, memory_input):
    if set(value) != {'written'} or not value['written']:
        raise ValueError('writer returned no serialized records')
    outputs = []
    for item in value['written']:
        if set(item) != {'path', 'text'}:
            raise ValueError('invalid writer result')
        validate_written_format(feature, item['path'], item['text'])
        outputs.extend(writer_value(feature, memory_input(item['path'], item['text']))['records'])
    for index, output in enumerate(outputs):
        output['record_index'] = index
    return {'records': outputs}


def validate_written_format(feature, path, text):
    from pathlib import Path
    suffix = Path(path).suffix.lower()
    if feature.startswith('io.smiles.'):
        if suffix != '.smi':
            raise ValueError('SMILES writer must emit a .smi record')
        return
    sdf = feature.startswith('io.sdf.')
    if suffix != ('.sdf' if sdf else '.mol'):
        raise ValueError('writer emitted the wrong file format')
    version = 'V3000' if '.v3000.' in feature else 'V2000'
    blocks = text.split('$$$$') if sdf else [text]
    if sdf and (not text.rstrip().endswith('$$$$') or not blocks[:-1]):
        raise ValueError('SDF writer must terminate each record with $$$$')
    for index, block in enumerate(blocks[:-1] if sdf else blocks):
        # Only the separator newline is removed; blank titles are significant.
        if index:
            block = block.removeprefix('\r\n').removeprefix('\n')
        lines = block.splitlines()
        if len(lines) < 4 or not lines[3].rstrip().endswith(version):
            raise ValueError(f'writer must emit {version}')


def atom_json(atom: Any) -> dict[str, Any]:
    radical, unpaired_electrons = radical_json(atom)
    return {
        "index": atom.GetIdx(),
        "atomic_number": atom.GetAtomicNum(),
        "symbol": atom.GetSymbol(),
        "formal_charge": atom.GetFormalCharge(),
        "isotope": atom.GetIsotope() or None,
        "explicit_hydrogens": atom.GetNumExplicitHs(),
        "atom_map": atom.GetAtomMapNum() or None,
        "radical": radical,
        "unpaired_electrons": unpaired_electrons,
        "aromatic": atom.GetIsAromatic(),
    }


def radical_json(atom: Any) -> tuple[str | None, int]:
    unpaired_electrons = atom.GetNumRadicalElectrons()
    if unpaired_electrons == 0:
        return None, 0
    if unpaired_electrons == 1:
        return "DOUBLET", 1
    if unpaired_electrons == 2:
        return "TRIPLET", 2
    return None, unpaired_electrons


def read_sdf_blocks(fixture_path) -> list[str]:
    text = fixture_path.read_text(encoding="utf-8", errors="replace")
    blocks: list[str] = []
    current: list[str] = []
    for line in text.splitlines():
        if line == "$$$$":
            blocks.append("\n".join(current) + "\n")
            current = []
        else:
            current.append(line)
    if current:
        blocks.append("\n".join(current) + "\n")
    return blocks
