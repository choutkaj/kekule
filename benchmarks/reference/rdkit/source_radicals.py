"""Preserve explicit spin assertions that RDKit's atom model does not retain.

RDKit owns molecular parsing and chemical calculations. These readers extract
only radical annotations from a source already accepted by RDKit. They never
derive spin from electron count or modify RDKit's chemical state.

Sources: BIOVIA CTfile Formats 2020, V2000/V3000 atom radical fields;
https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf
Chemaxon CXSMILES, Radical numbers;
https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html
"""
import re


_ELECTRONS = '_benchmarkSourceRadicalElectrons'
_SPIN = '_benchmarkSourceSpinMultiplicity'
_CTAB = {1: (2, 1), 2: (1, 2), 3: (2, 3)}
_CX = {1: (1, None), 2: (2, None), 3: (2, 1), 4: (2, 3),
       5: (3, None), 6: (3, 2), 7: (3, 4)}


def _attach(mol, assertions):
    # Source atomProp values are arbitrary metadata, including private-looking
    # names. Only annotations extracted here may assert chemical spin.
    for atom in mol.GetAtoms():
        for key in (_ELECTRONS, _SPIN):
            if atom.HasProp(key):
                atom.ClearProp(key)
    for index, (electrons, spin) in assertions.items():
        if not 0 <= index < mol.GetNumAtoms():
            raise ValueError('source radical atom index is out of range')
        atom = mol.GetAtomWithIdx(index)
        atom.SetIntProp(_ELECTRONS, electrons)
        if spin is not None:
            atom.SetIntProp(_SPIN, spin)


def observation(atom):
    """Return RDKit electron occupancy and independently retained source spin."""
    electrons = atom.GetNumRadicalElectrons()
    if atom.HasProp(_ELECTRONS) and atom.GetIntProp(_ELECTRONS) != electrons:
        raise ValueError('RDKit changed the explicitly supplied radical electron count')
    spin = atom.GetIntProp(_SPIN) if atom.HasProp(_SPIN) else None
    return electrons, spin


def attach_ctab(mol, block):
    lines = block.splitlines()
    if len(lines) < 4:
        raise ValueError('missing CTAB counts line')
    if lines[3].rstrip().endswith('V3000'):
        codes = _v3000(lines)
    elif lines[3].rstrip().endswith('V2000') or not lines[3][33:].strip():
        codes = _v2000(lines)
    else:
        raise ValueError('unsupported CTAB version for source radical observations')
    if len(codes) != mol.GetNumAtoms():
        raise ValueError('RDKit/source atom correspondence changed during CTAB reading')
    _attach(mol, {i: _CTAB[code] for i, code in enumerate(codes) if code})


def _v2000(lines):
    count, bonds, lists = (int(lines[3][i:i + 3].strip() or '0') for i in (0, 3, 6))
    codes = [2 if int(line[36:39].strip() or '0') == 4 else 0
             for line in lines[4:4 + count]]
    entries = {}
    override = False
    i = 4 + count + bonds + lists
    while i < len(lines):
        line = lines[i]
        i += 1
        fields = line.split()
        if fields == ['M', 'END']:
            break
        if fields[:1] in (['A'], ['G']):
            i += 1  # The following alias/group label is text, not a property.
        elif fields[:2] == ['S', 'SKP']:
            i += int(fields[2])
        elif fields[:2] in (['M', 'CHG'], ['M', 'RAD']):
            override = True
            if fields[1] == 'CHG':
                continue
            n = int(fields[2])
            if len(fields) != 3 + 2 * n:
                raise ValueError('malformed source M RAD entry count')
            for j in range(n):
                index, code = int(fields[3 + 2 * j]) - 1, int(fields[4 + 2 * j])
                if not 0 <= index < count or code not in (0, 1, 2, 3):
                    raise ValueError('invalid source M RAD entry')
                entries[index] = code
    if override:
        # Either property supersedes *all* charge/radical atom-block fields.
        codes = [0] * count
    for index, code in entries.items():
        codes[index] = code
    return codes


def _v3000_tokens(text):
    """Split logical atom records without interpreting quoted/list properties."""
    tokens, start, depth, quoted = [], None, 0, False
    for i, char in enumerate(text):
        if char.isspace() and not quoted and depth == 0:
            if start is not None:
                tokens.append(text[start:i])
                start = None
            continue
        if start is None:
            start = i
        if char == '"':
            quoted = not quoted
        elif not quoted:
            if char == '(':
                depth += 1
            elif char == ')':
                depth -= 1
    if quoted or depth != 0:
        raise ValueError('unbalanced V3000 property syntax')
    if start is not None:
        tokens.append(text[start:])
    return tokens


def _v3000(lines):
    codes, ids, pending, in_atoms = [], set(), '', False
    for line in lines[4:]:
        if line.strip() == 'M  END':
            break
        if not line.startswith('M  V30 '):
            continue
        text = line[7:].rstrip()
        if text.endswith('-'):
            pending += text[:-1]
            continue
        text, pending = pending + text, ''
        if text == 'BEGIN ATOM':
            in_atoms = True
        elif text == 'END ATOM':
            in_atoms = False
        elif in_atoms:
            fields = _v3000_tokens(text)
            if len(fields) < 6 or int(fields[0]) in ids:
                raise ValueError('invalid V3000 atom correspondence')
            ids.add(int(fields[0]))
            radical = [field for field in fields[6:] if field.startswith('RAD=')]
            if len(radical) > 1 or radical and not re.fullmatch(r'RAD=[0-3]', radical[0]):
                raise ValueError('invalid V3000 source radical property')
            codes.append(int(radical[0][4:]) if radical else 0)
    if pending or in_atoms:
        raise ValueError('incomplete V3000 atom block')
    return codes


def attach_cx(mol):
    if not mol.HasProp('_CXSMILES_Data'):
        return
    text = mol.GetProp('_CXSMILES_Data')
    if not text.startswith('|') or not text.endswith('|'):
        raise ValueError('incomplete preserved CXSMILES extension')
    # Labels and coordinates can contain characters that resemble field syntax.
    # Mask their contents without disturbing top-level commas or atom indices.
    visible, label, depth, field_start = [], False, 0, True
    for char in text[1:-1]:
        if label:
            if char == '$':
                label = False
            visible.append(' ')
        elif depth:
            if char == '(':
                depth += 1
            elif char == ')':
                depth -= 1
            visible.append(' ')
        elif field_start and char == '$':
            label, field_start = True, False
            visible.append(' ')
        elif field_start and char == '(':
            depth, field_start = 1, False
            visible.append(' ')
        else:
            visible.append(char)
            if not char.isspace():
                field_start = char == ','
    if label or depth != 0:
        raise ValueError('unbalanced preserved CXSMILES extension')
    assertions = {}
    for match in re.finditer(r'(?:^|,)\s*\^([1-7]):(\d+(?:,\d+)*)(?=\s*(?:,|$))', ''.join(visible)):
        value = _CX[int(match[1])]
        for index in map(int, match[2].split(',')):
            if index in assertions and assertions[index] != value:
                raise ValueError('conflicting CXSMILES source radical assertions')
            assertions[index] = value
    _attach(mol, assertions)
