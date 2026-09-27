"""Offline provenance check for independently generated OFFXML fixtures."""
import gzip
import hashlib
import json
import unittest
from paths import historical, historical_bytes
from pathlib import Path

from paths import HERE
ROOT = HERE.parent.parent


class OffxmlReferenceTests(unittest.TestCase):
    def test_external_source_and_reference_fingerprints(self):
        lock = historical('offxml-generalization.lock.json')
        for row in lock['artifacts']:
            self.assertEqual(hashlib.sha256(historical_bytes(Path(row['path']).name) if row['path'].startswith('benchmarks/openff/') else (ROOT / row['path']).read_bytes()).hexdigest(), row['sha256'])
        fixture = ROOT / 'crates/kekule-openff/tests/fixtures/offxml-generalization.json.gz'
        report = json.loads(gzip.decompress(fixture.read_bytes()))
        self.assertEqual(report['source_sha256'], lock['artifacts'][0]['sha256'])
        self.assertEqual(report['toolkit_version'], '0.19.0')
        self.assertEqual(len({r['name'] for r in report['records']}), 6)
        for row in report['records']:
            self.assertIn('NAGLCharges', row['xml'])
            self.assertEqual(set(row['parameters']), {'Bonds', 'Angles', 'ProperTorsions',
                             'ImproperTorsions', 'Constraints', 'vdW', 'LibraryCharges'})


if __name__ == '__main__':
    unittest.main()
