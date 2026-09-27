"""Standard-library tests for the optional scientific validation harness."""
import copy
import gzip
import hashlib
import json
import math
import unittest
from pathlib import Path
from zipfile import ZipFile

from robustness import compare_parameters, digest

from paths import ARCHIVE, HERE, historical, historical_bytes


class RobustnessTests(unittest.TestCase):
    def test_frozen_panel_has_external_provenance_and_complete_coordinates(self):
        inputs=json.loads(gzip.decompress((HERE/'data/inputs.json.gz').read_bytes()))
        rows=inputs['records']
        self.assertEqual(len({r['id'] for r in rows}),110)
        self.assertEqual(sum(r['kind']=='small-molecule' for r in rows),100)
        self.assertEqual(sum(r['kind']=='protein-chain' for r in rows),10)
        for r in rows:
            self.assertEqual(len(r['source']['sha256']),64)
            self.assertEqual(len(r['coordinates_nm']),r['atoms'])
            self.assertTrue(all(len(v)==3 and all(math.isfinite(x) for x in v) for v in r['coordinates_nm']))
        self.assertTrue(any('error' in r for r in inputs['protein_preparation_attempts']))

    def test_parameter_comparison_catches_omissions_and_changed_values(self):
        expected={h:[] for h in ['Bonds','Angles','Constraints','ProperTorsions','ImproperTorsions','vdW']}
        native=dict(maps=[1],system={h:[] for h in ['bonds','angles','constraints','propers','impropers','vdw']})
        expected['vdW']=[dict(atoms=[1],mult=None,values=dict(sigma=.3,epsilon=.4))]
        native['system']['vdw']=[dict(id='n1',sigma=.3,epsilon=.4)]
        self.assertEqual(compare_parameters(expected,native)[0],[])
        changed=copy.deepcopy(native);changed['system']['vdw'][0]['epsilon']=.5
        self.assertTrue(compare_parameters(expected,changed)[0])
        changed=copy.deepcopy(native);changed['system']['vdw']=[]
        self.assertTrue(compare_parameters(expected,changed)[0])
        changed=copy.deepcopy(native);del changed['system']['vdw'][0]['epsilon']
        self.assertTrue(compare_parameters(expected,changed)[0])

    def test_parameter_comparison_preserves_improper_multiplicity(self):
        expected={h:[] for h in ['Bonds','Angles','Constraints','ProperTorsions','ImproperTorsions','vdW']}
        native=dict(maps=[1,2,3,4],system={h:[] for h in ['bonds','angles','constraints','propers','impropers','vdw']})
        p=dict(k=2.0,phase=math.pi,periodicity=2,idivf=3.)
        expected['ImproperTorsions']=[dict(atoms=[2,1,3,4],mult=0,values=p)]*3
        native['system']['impropers']=[dict(maps=[2,1,3,4],parameter=dict(id='i1',terms=[p]))]*3
        self.assertEqual(compare_parameters(expected,native)[0],[])
        native['system']['impropers'].pop()
        self.assertTrue(compare_parameters(expected,native)[0])

    def test_archived_fingerprints(self):
        lock=HERE/'results/manifest.json'
        for row in json.loads(lock.read_text())['artifacts']:
            self.assertEqual(digest(HERE/row['path']),row['sha256'],row['path'])

    def test_cleanup_preserves_every_historical_artifact(self):
        with ZipFile(ARCHIVE) as archive:
            index = json.loads(archive.read('index.json'))
            self.assertEqual(set(archive.namelist()), set(index) | {'index.json'})
            for name, sha in index.items():
                self.assertEqual(hashlib.sha256(archive.read(name)).hexdigest(), sha, name)
        for row in historical('robustness-reports.lock.json')['artifacts']:
            self.assertEqual(hashlib.sha256(historical_bytes(row['path'])).hexdigest(), row['sha256'])

    def test_timings_cover_the_entire_external_panel(self):
        inputs = json.loads(gzip.decompress((HERE/'data/inputs.json.gz').read_bytes()))
        timings = json.loads((HERE/'results/timings.json').read_bytes())
        self.assertEqual(timings['inputs_sha256'], digest(HERE/'data/inputs.json.gz'))
        self.assertEqual([(r['id'], r['atoms']) for r in inputs['records']],
                         [(r['id'], r['atoms']) for r in timings['records']])
        self.assertEqual(timings['warmups_per_operation'], 1)
        self.assertEqual(timings['repetitions'], dict(features_ms=5,
            inference_including_features_ms=5, assign_charges_ms=5, full_parameterization_ms=3))
        for row in timings['records']:
            for key, count in timings['repetitions'].items():
                self.assertEqual(len(row[key]), count)
                self.assertTrue(all(math.isfinite(v) and v > 0 for v in row[key]))


if __name__=='__main__':
    unittest.main()
