import unittest
import sys
from pathlib import Path

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from reference.rdkit import source_radicals as source


def mol_block(code, v3000=False):
    block = Chem.MolToMolBlock(Chem.MolFromSmiles('[CH2]'), forceV3000=v3000)
    if v3000:
        return block.replace('RAD=3', f'RAD={code}')
    return block.replace('M  RAD  1   1   3', f'M  RAD  1   1   {code}')


def read_ctab(block):
    mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False, strictParsing=True)
    if mol is None:
        raise ValueError('fixture rejected by RDKit')
    source.attach_ctab(mol, block)
    Chem.SanitizeMol(mol)
    return mol


def read_smiles(text):
    params = Chem.SmilesParserParams()
    params.sanitize = False
    params.removeHs = False
    mol = Chem.MolFromSmiles(text, params)
    if mol is None:
        raise ValueError('fixture rejected by RDKit')
    source.attach_cx(mol)
    Chem.SanitizeMol(mol)
    return mol


class SourceRadicalTests(unittest.TestCase):
    def test_mol_spin_is_independent_of_electron_count(self):
        for v3000 in (False, True):
            for code, expected in [(1, (2, 1)), (2, (1, 2)), (3, (2, 3))]:
                with self.subTest(v3000=v3000, code=code):
                    # Start with the corresponding occupancy so declared CTAB
                    # valence does not force RDKit to change the source state.
                    original = '[CH3]' if code == 2 else '[CH2]'
                    block = Chem.MolToMolBlock(Chem.MolFromSmiles(original), forceV3000=v3000)
                    if code != 2:
                        block = mol_block(code, v3000)
                    mol = read_ctab(block)
                    self.assertEqual(source.observation(mol.GetAtomWithIdx(0)), expected)

    def test_plain_smiles_does_not_assert_spin(self):
        for text, count in [('[C]', 4), ('[CH]', 3), ('[CH2]', 2), ('[CH3]', 1), ('C', 0)]:
            with self.subTest(text=text):
                mol = read_smiles(text)
                self.assertEqual(source.observation(mol.GetAtomWithIdx(0)), (count, None))

    def test_cx_explicit_and_unspecified_spin(self):
        for code, count, spin in [(1, 1, None), (2, 2, None), (3, 2, 1), (4, 2, 3),
                                  (5, 3, None), (6, 3, 2), (7, 3, 4)]:
            text = f'[CH{4-count}] |^' + str(code) + ':0| title'
            with self.subTest(text=text):
                mol = read_smiles(text)
                self.assertEqual(source.observation(mol.GetAtomWithIdx(0)), (count, spin))

    def test_cx_labels_are_not_radical_annotations(self):
        mol = read_smiles('[CH2] |$^4:0$|')
        self.assertEqual(source.observation(mol.GetAtomWithIdx(0)), (2, None))
        mol = read_smiles('[CH2] |^3:0,$^4:0$|')
        self.assertEqual(source.observation(mol.GetAtomWithIdx(0)), (2, 1))

    def test_cx_atom_lists_and_separate_fields(self):
        mol = read_smiles('[CH2].[CH2].[CH] |^3:0,1,^7:2|')
        self.assertEqual([source.observation(a) for a in mol.GetAtoms()], [(2, 1), (2, 1), (3, 4)])

    def test_cx_literal_atom_properties_do_not_open_labels_or_coordinates(self):
        for value in ['quote " unclosed', 'paren ( unclosed', 'dollar $ unclosed']:
            mol = Chem.MolFromSmiles('[CH2]')
            mol.GetAtomWithIdx(0).SetProp('note', value)
            restored = read_smiles(Chem.MolToCXSmiles(mol))
            self.assertEqual(restored.GetAtomWithIdx(0).GetProp('note'), value)
            self.assertEqual(source.observation(restored.GetAtomWithIdx(0)), (2, None))

    def test_arbitrary_atom_properties_cannot_assert_radical_state(self):
        text = ('[CH2] |atomProp:0._benchmarkSourceSpinMultiplicity.3'
                ':0._benchmarkSourceRadicalElectrons.99|')
        mol = read_smiles(text)
        self.assertEqual(source.observation(mol.GetAtomWithIdx(0)), (2, None))

    def test_annotation_follows_fragment_and_hydrogen_transforms(self):
        mol = read_smiles('[CH2].C |^3:0|')
        fragments = Chem.GetMolFrags(mol, asMols=True)
        transformed = Chem.RemoveHs(Chem.AddHs(fragments[0]))
        self.assertEqual(source.observation(transformed.GetAtomWithIdx(0)), (2, 1))
        self.assertEqual(source.observation(fragments[1].GetAtomWithIdx(0)), (0, None))

    def test_source_annotations_do_not_enter_cx_writer_identity(self):
        block = mol_block(1)
        plain = Chem.MolFromMolBlock(block, removeHs=False)
        annotated = read_ctab(block)
        self.assertEqual(Chem.MolToCXSmiles(plain), Chem.MolToCXSmiles(annotated))

    def test_reference_cannot_silently_change_asserted_electron_count(self):
        mol = read_smiles('[CH4] |^3:0|')
        with self.assertRaisesRegex(ValueError, 'changed.*electron count'):
            source.observation(mol.GetAtomWithIdx(0))

    def test_v2000_property_precedence_and_text_scope(self):
        block = Chem.MolToMolBlock(Chem.MolFromSmiles('CC'))
        lines = block.splitlines()
        lines[4] = lines[4][:36] + '  4' + lines[4][39:]
        base = '\n'.join(lines) + '\n'
        mol = read_ctab(base)
        self.assertEqual(source._v2000(base.splitlines()), [2, 0])
        # RDKit 2026.03.3 ignores the legacy atom-block doublet code. Retain
        # the source assertion and report lost chemistry as reference failure.
        with self.assertRaisesRegex(ValueError, 'changed.*electron count'):
            source.observation(mol.GetAtomWithIdx(0))
        for property_line in ['M  CHG  1   2   1', 'M  RAD  1   2   2']:
            changed = base.replace('M  END', property_line + '\nM  END')
            mol = read_ctab(changed)
            self.assertEqual(source.observation(mol.GetAtomWithIdx(0))[1], None)
        for tail in ['A    1\nM  RAD  1   1   3\n', 'G    1\nM  RAD  1   1   3\n']:
            codes = source._v2000(base.replace('M  END', tail + 'M  END').splitlines())
            self.assertEqual(codes, [2, 0])
        self.assertEqual(source._v2000((base + '>  <note>\nM  RAD  1   1   3\n').splitlines()), [2, 0])

    def test_v3000_source_order_and_continuation(self):
        block = Chem.MolToMolBlock(Chem.MolFromSmiles('[CH2].[CH3]'), forceV3000=True)
        block = block.replace('M  V30 1 C ', 'M  V30 17 C ').replace('M  V30 2 C ', 'M  V30 4 C ')
        block = block.replace('RAD=3 VAL=2', 'RAD=1 -\nM  V30 VAL=2')
        mol = read_ctab(block)
        self.assertEqual([source.observation(a) for a in mol.GetAtoms()], [(2, 1), (1, 2)])

    def test_v3000_quoted_and_list_properties_are_not_radical_codes(self):
        fields = source._v3000_tokens('1 C 0 0 0 0 NOTE="RAD=2 with ""quotes""" RGROUPS=(1 7) RAD=1')
        self.assertEqual(fields[6:], ['NOTE="RAD=2 with ""quotes"""', 'RGROUPS=(1 7)', 'RAD=1'])


if __name__ == '__main__':
    unittest.main()
