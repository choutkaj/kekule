"""Offline end-to-end checks of the complete gradient validation observations."""
import math
import unittest

from gradients import COMPONENTS
from paths import HERE
from robustness import read, digest


class GradientValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = read(HERE / 'data/gradients.json.gz')
        cls.native = read(HERE / 'results/gradients.json.gz')

    def test_complete_panel_retains_known_failures_and_provenance(self):
        inputs = read(HERE / 'data/inputs.json.gz')['records']
        expected = [(r['id'], reverse) for r in inputs for reverse in (False, True)]
        self.assertEqual([(r['id'], r['reverse']) for r in self.native['records']], expected)
        self.assertEqual([r['id'] for r in self.reference['records']], [r['id'] for r in inputs])
        self.assertEqual(self.native['reference_sha256'], digest(HERE / 'data/gradients.json.gz'))
        self.assertEqual(self.native['energy_reference_sha256'], digest(HERE / 'data/reference.json.gz'))
        self.assertEqual(self.reference['energy_reference_sha256'], self.native['energy_reference_sha256'])
        self.assertEqual(self.native['inputs_sha256'], digest(HERE / 'data/inputs.json.gz'))
        failures = [r for r in self.native['records'] if not r['passed']]
        self.assertEqual([(r['id'], r['reverse']) for r in failures], [('pubchem-443915', False), ('pubchem-443915', True)])
        for row in failures:
            self.assertTrue(all(isinstance(e, dict) and e['frame'] == 0 and e['check'] in ('OpenMM', 'rotation') for e in row['errors']))

    def test_every_cartesian_reference_component_is_asserted(self):
        reference = {r['id']: r for r in self.reference['records']}
        total = 0
        for row in self.native['records']:
            n = len(row['maps'])
            self.assertEqual(sorted(row['maps']), list(range(1, n + 1)))
            order = sorted(range(n), key=lambda i: row['maps'][i])
            self.assertEqual(len(row['observations']), 3)
            self.assertEqual(len(row['comparisons']), 3 * 2 * len(COMPONENTS))
            self.assertEqual(len(row['invariants']), 3 * 2 * len(COMPONENTS))
            for frame, observation in enumerate(row['observations']):
                for source in ('reference_charges', 'native_charges'):
                    self.assertEqual(set(observation[source]['gradients']), set(COMPONENTS))
                    for component, vectors in observation[source]['gradients'].items():
                        self.assertEqual(len(vectors), n)
                        self.assertTrue(all(len(v) == 3 and all(math.isfinite(x) for x in v) for v in vectors))
                        if source != 'reference_charges': continue
                        ref = reference[row['id']]['gradients'][frame][component]
                        errors = [abs(vectors[j][axis] - ref[i][axis]) for i,j in enumerate(order) for axis in range(3)]
                        summary = next(s for s in row['comparisons'] if s['frame'] == frame and s['source'] == source and s['component'] == component)
                        self.assertEqual(max(errors), summary['max_abs'])
                        self.assertEqual(len(errors), summary['values'])
                        total += len(errors)
        self.assertEqual(total, 1842876)

    def test_energy_derivatives_converge_and_reference_disagreement_is_visible(self):
        finite = [d for row in self.native['records'] for d in row['finite_differences']]
        self.assertEqual(len(finite), 25410)
        self.assertTrue(all(d['failed'] == 0 for d in finite))
        self.assertTrue(all(not r['finite_differences'] if r['reverse'] else len(r['finite_differences']) == 231 for r in self.native['records']))
        case = next(r for r in self.native['records'] if r['id'] == 'pubchem-443915' and not r['reverse'])
        self.assertEqual([len(a) for a in case['near_linear_angles']], [1, 0, 0])
        diagnostic = case['openmm_finite_differences'][0]
        reference_failures = [c for c in diagnostic['checks'] if c['reference']['failed']]
        self.assertTrue(reference_failures)
        self.assertEqual({c['component'] for c in reference_failures}, {'Angles', 'Total'})
        self.assertGreater(max(c['reference']['max_abs'] for c in reference_failures), 200.)
        self.assertLess(max(c['native']['max_abs'] for c in reference_failures), .05)


if __name__ == '__main__':
    unittest.main()
