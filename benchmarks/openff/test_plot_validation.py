"""Check scientific figure pairing using the externally supplied validation panel."""
import copy
import unittest

from plot_validation import HERE, collect, digest, metrics, read


class PlotValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        reference = read(HERE / 'robustness-reference.json.gz')
        native = read(HERE / 'robustness-native.json.gz')
        chosen = reference['records'][0]
        cls.reference = dict(inputs_sha256=reference['inputs_sha256'], records=[chosen])
        cls.native = dict(inputs_sha256=native['inputs_sha256'],
                          records=[r for r in native['records'] if r['id'] == chosen['id']])

    def test_external_atom_orders_and_geometry_counts_are_preserved(self):
        series, energies = collect(self.reference, self.native)
        self.assertEqual(len(series['charge']), 2 * len(self.reference['records'][0]['charges']))
        self.assertEqual(len(energies['native_charges.Total']), 6)
        self.assertAlmostEqual(metrics(series['charge'])['max_abs'],
                               max(r['charge_max_error_e'] for r in self.native['records']), places=15)

    def test_atom_order_pairing_is_invariant_to_storage_order(self):
        changed = copy.deepcopy(self.native)
        for row in changed['records']:
            obs = row['native']
            for values in (obs['maps'], obs['system']['charges'], obs['system']['vdw'],
                           obs['features']['values']):
                values.reverse()
        self.assertEqual(collect(self.reference, self.native), collect(self.reference, changed))

    def test_missing_cases_maps_and_torsion_multiplicities_are_rejected(self):
        for mutate in (
            lambda d: d['records'].pop(),
            lambda d: d['records'][0]['native']['maps'].pop(),
            lambda d: d['records'][0]['native']['system']['impropers'].pop(),
        ):
            changed = copy.deepcopy(self.native)
            mutate(changed)
            with self.assertRaises(ValueError):
                collect(self.reference, changed)

    def test_rendered_summary_matches_frozen_sources_and_reported_maxima(self):
        summary = read(HERE / 'figures/validation-summary.json')
        for path, sha in summary['source_sha256'].items():
            self.assertEqual(digest(HERE / path), sha)
        self.assertEqual(digest(HERE / 'cpu-timings.json'), summary['timing_sha256'])
        frozen = read(HERE / 'robustness-summary.json')
        self.assertEqual(summary['parameters']['charge']['max_abs'], frozen['charge_max_error_e'])
        for component, expected in frozen['energy_shared_charge_max_error_kj_mol'].items():
            self.assertEqual(summary['energies']['reference_charges.' + component]['max_abs'], expected)
        self.assertEqual(summary['cases'], 220)
        self.assertEqual(summary['all_checks_passed'], 214)


if __name__ == '__main__':
    unittest.main()
