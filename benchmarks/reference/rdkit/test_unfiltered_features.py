import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location('rdkit_features', Path(__file__).with_name('run_feature.py'))
features = importlib.util.module_from_spec(spec)
spec.loader.exec_module(features)


class UnfilteredFeaturesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rdkit = features.import_rdkit()

    def test_smiles_reader_attempts_stereo_wildcards_and_invalid_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.smi'
            path.write_text('F[C@H](Cl)Br atom\nF/C=C/F bond\nF\\C=C\\F reverse\n* wildcard\nF[C@ invalid\n')
            for sanitize in (False, True):
                records = features.read_smiles_records(path, self.rdkit['Chem'], sanitize)
                self.assertEqual([r['record_index'] for r in records], list(range(5)))
                self.assertEqual([r['status'] for r in records], ['ok'] * 4 + ['parse_error'])
            for feature in ('io.smiles.parse', 'io.smiles.write', 'stereo.cip'):
                records = features.evaluate(feature, path, self.rdkit)['records']
                self.assertEqual([r['record_index'] for r in records], list(range(5)))
                self.assertEqual([r['status'] for r in records], ['ok'] * 4 + ['parse_error'])

    def test_cip_retains_empty_descriptors_and_preparation_failures(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.smi'
            path.write_text('CC plain\nF[C@H](Cl)Br stereo\n')
            records = features.evaluate('stereo.cip', path, self.rdkit)['records']
            self.assertEqual(len(records), 2)
            self.assertEqual(records[0]['atom_count'], 2)
            self.assertEqual(records[0]['atom_descriptors'], [])
            self.assertEqual(records[0]['bond_descriptors'], [])
            self.assertTrue(records[1]['atom_descriptors'])
            source = features.read_smiles_records(path, self.rdkit['Chem'], False)[0]
            with patch.object(features, 'clone_and_sanitize', return_value=None):
                self.assertEqual(features.stereo_cip_record(source, self.rdkit['Chem'])['status'], 'sanitize_error')
            with patch.object(self.rdkit['Chem'], 'AssignCIPLabels', side_effect=RuntimeError('failed')):
                self.assertEqual(features.stereo_cip_record(source, self.rdkit['Chem'])['status'], 'cip_error')


if __name__ == '__main__':
    unittest.main()
