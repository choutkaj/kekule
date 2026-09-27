"""Offline provenance checks for the independently generated two-model panel."""
import gzip
import hashlib
import json
from pathlib import Path
import unittest
from paths import historical_bytes

from paths import HERE
ROOT = HERE.parent.parent


class ModelReferenceTests(unittest.TestCase):
    def test_archived_model_parity_retains_external_inputs_and_fingerprints(self):
        lock = json.loads((HERE / 'data/models.lock.json').read_bytes())
        for artifact in lock['artifacts']:
            self.assertEqual(hashlib.sha256((ROOT / artifact['path']).read_bytes()).hexdigest(), artifact['sha256'])
        reference_path = ROOT / 'crates/kekule-openff/tests/fixtures/models.json.gz'
        reference = json.loads(gzip.decompress(reference_path.read_bytes()))
        native = json.loads(gzip.decompress((HERE / 'results/models.json.gz').read_bytes()))
        self.assertEqual(native['reference_sha256'], hashlib.sha256(reference_path.read_bytes()).hexdigest())
        for path, sha in reference['sources_sha256'].items():
            data = (HERE / path).read_bytes() if path.startswith('fixtures/') else historical_bytes(path)
            self.assertEqual(hashlib.sha256(data).hexdigest(), sha)
        self.assertEqual(len(reference['models']), 2)
        self.assertEqual(len(native['models']), 2)
        for expected, actual in zip(reference['models'], native['models']):
            self.assertEqual(expected['name'], actual['name'])
            self.assertEqual(len(expected['records']), 66)
            self.assertEqual(actual['passed'], 66)
            self.assertEqual(len(actual['records']), 66)
            self.assertEqual([(r['id'], r['reverse']) for r in expected['records']],
                             [(r['id'], r['reverse']) for r in actual['records']])
            self.assertTrue(all(not r['errors'] for r in actual['records']))


if __name__ == '__main__':
    unittest.main()
