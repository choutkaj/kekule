import copy
from contextlib import redirect_stderr
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import dashboard


class DashboardTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'contract.json').write_text('{"schema": 2}\n')
        lock = self.root / 'corpora/example/sources.lock.json'
        lock.parent.mkdir(parents=True)
        lock.write_text(json.dumps({'corpus_id': 'example', 'entries': [{'id': 'a'}, {'id': 'b'}]}))
        self.golden = {'schema': 2, 'dataset': 'example', 'feature': 'io.smiles.parse',
                       'cases': 4, 'sha256': 'a' * 64, 'input_lock_sha256': dashboard.text_hash(lock),
                       'contract_sha256': dashboard.text_hash(self.root / 'contract.json'),
                       'reference_code_sha256': None,
                       'reference': {'tool': 'independent', 'version': '1'}}
        self.manifest = self.root / 'goldens/example/io.smiles.parse.jsonl.meta.json'
        self.manifest.parent.mkdir(parents=True)
        self.manifest.write_text(json.dumps(self.golden))
        row = dict.fromkeys(dashboard.COUNTS, 0)
        row.update(dataset='example', feature='io.smiles.parse', source_ids=1, cases=1,
                   agrees=1, exact_agrees=1, not_applicable=1, golden=self.golden,
                   kekule_ms=0.1, reference_ms=0.0)
        self.report = {'schema': 2, 'mode': 'compare', 'complete': True, 'passed': True,
                       'error': None, 'goldens': 'C:/private/goldens',
                       'cases': 'C:/private/cases.jsonl',
                       'implementation': {'revision': 'b' * 40, 'dirty': False,
                                          'contract_sha256': self.golden['contract_sha256'],
                                          'executable_sha256': 'c' * 64,
                                          'reference_code_sha256': 'd' * 64,
                                          'working_tree_status_sha256': 'e' * 64},
                       'results': [row]}
        self.path = self.root / 'run.json'

    def load(self, report=None):
        self.path.write_text(json.dumps(self.report if report is None else report))
        return dashboard.load_report(self.path, dashboard.catalogue(self.root))

    def test_samples_keep_unavailable_formats_outside_the_denominator(self):
        row = self.load()['results'][0]
        self.assertEqual(row['coverage'], 'sampled')
        self.assertEqual((row['cases'], row['agrees'], row['not_applicable']), (1, 1, 1))

    def test_corpora_sizes_and_formats_come_from_locked_entries_and_packs(self):
        lock_path = self.root / 'corpora/example/sources.lock.json'
        lock = json.loads(lock_path.read_text())
        lock['entries'][0]['files'] = [{'path': 'data/a.sdf'}, {'path': 'data/shared.txt'}]
        lock['entries'][1]['files'] = [{'path': 'data/shared.txt'}, {'path': 'data/c.cif'},
                                      {'path': 'data/d.mol'}]
        lock['packs'] = [{'path': 'data/packs/b.sdf'}, {'path': 'data/packs/b.sdf'}]
        lock_path.write_text(json.dumps(lock))
        corpus = dashboard.catalogue(self.root)['datasets'][0]
        self.assertEqual(corpus['source_ids'], 2)
        self.assertEqual(corpus['formats'], ['MOL', 'SDF', 'SMILES', 'mmCIF'])
        self.assertFalse((self.root / 'corpora/example/data').exists())

    def test_simplified_page_has_corpora_without_removed_panels_or_logo(self):
        self.load()
        html = dashboard.render([self.path], self.root)
        self.assertIn('<h2>Data corpora</h2>', html)
        self.assertIn('id="corpora"', html)
        self.assertIn('<th scope="col">N</th>', html)
        self.assertIn('<label for="run">Benchmark run</label>', html)
        for removed in ('<th scope="col">Files</th>', '<th scope="col">Source IDs</th>',
                        'Counts describe the complete locked collections.',
                        'Files may contain many records.', 'locked source IDs',
                        'Feature / independent observation', 'Comparison run',
                        'Agreeing / applicable cases.'):
            self.assertNotIn(removed, html)
        for removed in ('metrics', 'coverage-chart', 'outcome-bar', 'detail', 'error-chart', 'notice', 'revision'):
            self.assertNotIn(f'id="{removed}"', html)
        self.assertNotIn('Independent reference comparisons', html)
        self.assertNotIn('Keep the scope attached', html)
        self.assertNotIn('<svg', html)
        self.assertNotIn('kekulé', html.lower())
        self.assertIn('kekule · Benchmarks', html)

    def test_single_agreement_category_preserves_report_counts(self):
        self.report['results'][0].update(source_ids=2, cases=3, agrees=3, exact_agrees=1)
        self.load()
        html = dashboard.render([self.path], self.root)
        self.assertIn('<i class="swatch agree"></i>Agrees', html)
        self.assertNotIn('Within precision', html)
        self.assertNotIn('<i class="swatch exact">', html)
        embedded = html.split('<script id="benchmark-data" type="application/json">')[1].split('</script>')[0]
        row = json.loads(embedded)['runs'][0]['results'][0]
        self.assertEqual((row['cases'], row['agrees'], row['exact_agrees']), (3, 3, 1))

    def test_full_selection_requires_all_golden_records_to_be_accounted_for(self):
        row = self.report['results'][0]
        row['source_ids'] = 2
        with self.assertRaisesRegex(ValueError, 'every golden record'):
            self.load()
        row.update(cases=2, agrees=2, not_applicable=2)
        loaded = self.load()['results'][0]
        self.assertEqual(loaded['coverage'], 'full')
        self.assertEqual(loaded['agrees'] - loaded['exact_agrees'], 1)

    def test_empty_and_partial_runs_stay_incomplete(self):
        self.report.update(complete=False, passed=False, error='C:/private/failure.txt')
        partial = self.load()
        self.assertFalse(partial['complete'])
        self.assertEqual(len(partial['results']), 1)
        self.report['results'] = []
        self.assertEqual(self.load()['results'], [])

    def test_outcome_counts_and_pass_status_cannot_be_fabricated(self):
        for changes in ({'cases': 2}, {'exact_agrees': 2}, {'errors': 1},
                        {'source_ids': 3}, {'reference_errors': 1}):
            report = copy.deepcopy(self.report)
            report['results'][0].update(changes)
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.load(report)
        self.report['passed'] = False
        with self.assertRaisesRegex(ValueError, 'pass status'):
            self.load()

    def test_negative_boolean_and_nonfinite_measurements_are_rejected(self):
        for field, value in [('cases', -1), ('cases', True), ('cases', 2**53),
                             ('kekule_ms', float('nan')), ('reference_ms', float('inf'))]:
            report = copy.deepcopy(self.report)
            report['results'][0][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                self.load(report)

    def test_error_origins_may_overlap_but_do_not_inflate_error_cases(self):
        self.report['passed'] = False
        self.report['results'][0].update(agrees=0, exact_agrees=0, errors=1,
                                         reference_errors=1, kekule_errors=1)
        row = self.load()['results'][0]
        self.assertEqual(row['errors'], 1)
        self.assertEqual(row['reference_errors'] + row['kekule_errors'], 2)

    def test_generation_and_duplicate_rows_are_rejected(self):
        self.report['mode'] = 'generate'
        with self.assertRaisesRegex(ValueError, 'comparison reports'):
            self.load()
        self.report['mode'] = 'compare'
        self.report['results'] *= 2
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            self.load()

    def test_changed_reference_is_visible_as_stale(self):
        self.report['results'][0]['golden'] = dict(self.golden, sha256='f' * 64)
        self.assertEqual(self.load()['results'][0]['coverage'], 'stale')

    def test_reference_identity_and_contract_mismatches_are_rejected(self):
        self.report['results'][0]['golden'] = dict(self.golden, dataset='wrong')
        with self.assertRaisesRegex(ValueError, 'identity mismatch'):
            self.load()
        self.report['results'][0]['golden'] = dict(self.golden, contract_sha256='f' * 64)
        with self.assertRaisesRegex(ValueError, 'contract mismatch'):
            self.load()

    def test_duplicate_json_keys_are_rejected(self):
        self.path.write_text('{"schema":1,"schema":2}')
        with self.assertRaisesRegex(ValueError, 'duplicate JSON key'):
            dashboard.load_report(self.path, dashboard.catalogue(self.root))

    def test_export_needs_no_payloads_and_omits_paths_and_escapes_script_content(self):
        hostile = '</script><script>alert("private")</script>'
        self.golden['reference']['version'] = hostile
        self.manifest.write_text(json.dumps(self.golden))
        self.load()
        html = dashboard.render([self.path], self.root)
        self.assertNotIn(hostile, html)
        self.assertNotIn('C:/private', html)
        self.assertNotIn('DASHBOARD_DATA', html)
        embedded = html.split('<script id="benchmark-data" type="application/json">')[1].split('</script>')[0]
        self.assertEqual(json.loads(embedded)['runs'][0]['results'][0]['golden']['reference']['version'], hostile)

    def test_multiple_runs_are_kept_separate_and_identical_imports_rejected(self):
        self.load()
        with self.assertRaisesRegex(ValueError, 'duplicate comparison report'):
            dashboard.render([self.path, self.path], self.root)
        first = self.root / 'first.json'
        first.write_bytes(self.path.read_bytes())
        self.report['complete'] = self.report['passed'] = False
        self.load()
        html = dashboard.render([first, self.path], self.root)
        embedded = html.split('<script id="benchmark-data" type="application/json">')[1].split('</script>')[0]
        runs = json.loads(embedded)['runs']
        self.assertEqual(len(runs), 2)
        self.assertTrue(runs[0]['complete'])
        self.assertFalse(runs[1]['complete'])

    def test_output_cannot_overwrite_a_report_with_an_html_filename(self):
        report = self.root / 'report.html'
        original = json.dumps(self.report)
        report.write_text(original)
        stderr = io.StringIO()
        with patch('sys.argv', ['dashboard.py', str(report), '--output', str(report)]), \
                redirect_stderr(stderr), self.assertRaises(SystemExit) as stopped:
            dashboard.main()
        self.assertEqual(stopped.exception.code, 2)
        self.assertIn('must not replace an input report', stderr.getvalue())
        self.assertEqual(report.read_text(), original)

    def test_discovery_orders_by_recorded_start_and_keeps_failed_runs(self):
        runs_dir = self.root / 'runs'
        runs_dir.mkdir()
        older = dict(self.report, started_at_unix_ms=1000)
        newer = copy.deepcopy(self.report)
        newer.update(started_at_unix_ms=2000, passed=False)
        newer['results'][0].update(agrees=0, exact_agrees=0, disagrees=1)
        (runs_dir / 'older.json').write_text(json.dumps(older))
        (runs_dir / 'newer.json').write_text(json.dumps(newer))
        os.utime(runs_dir / 'older.json', (999999, 999999))
        os.utime(runs_dir / 'newer.json', (1, 1))
        first = dashboard.discovered(runs_dir, self.root)
        self.assertEqual([run['name'] for run in first['runs']], ['newer.json', 'older.json'])
        self.assertFalse(first['runs'][0]['passed'])
        os.utime(runs_dir / 'newer.json', (9999999, 9999999))
        self.assertEqual(dashboard.discovered(runs_dir, self.root)['signature'], first['signature'])

    def test_discovery_deduplicates_and_skips_invalid_or_generation_reports(self):
        self.load()
        runs_dir = self.root / 'runs'
        runs_dir.mkdir()
        for name in ('first.json', 'copy.json'):
            (runs_dir / name).write_bytes(self.path.read_bytes())
        (runs_dir / 'broken.json').write_text('{')
        (runs_dir / 'generate.json').write_text(json.dumps(dict(self.report, mode='generate')))
        (runs_dir / 'cases.jsonl').write_text('case observations are not summaries')
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            data = dashboard.discovered(runs_dir, self.root)
        self.assertEqual(len(data['runs']), 1)
        self.assertIn('skipped broken.json', stderr.getvalue())
        self.assertNotIn('generate.json', stderr.getvalue())

    def test_refresh_creates_live_empty_history_then_includes_an_incomplete_run(self):
        runs_dir = self.root / 'runs'
        output = dashboard.refresh(runs_dir, root=self.root)
        self.assertTrue(output.is_file())
        script = output.with_name('dashboard-data.js')
        self.assertIn('"runs":[]', script.read_text())
        self.report.update(complete=False, passed=False, started_at_unix_ms=1000)
        (runs_dir / 'partial.json').write_text(json.dumps(self.report))
        dashboard.refresh(runs_dir, root=self.root)
        raw = script.read_text()
        self.assertTrue(raw.startswith('window.updateKekuleBenchmarks('))
        self.assertIn('"complete":false', raw)
        self.assertNotIn('C:/private', raw)
        self.assertEqual(list(runs_dir.glob('*.tmp')), [])

    def test_invalid_run_start_times_are_rejected(self):
        for value in (-1, True, 1.5, 8_640_000_000_000_001):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.load(dict(self.report, started_at_unix_ms=value))

    def test_legacy_runner_names_order_without_inventing_recorded_timestamps(self):
        runs_dir = self.root / 'runs'
        runs_dir.mkdir()
        (runs_dir / 'run-99-1789591874492654300.json').write_text(json.dumps(self.report))
        older = dict(self.report, complete=False, passed=False)
        (runs_dir / 'run-999-1789417953494544600.json').write_text(json.dumps(older))
        runs = dashboard.discovered(runs_dir, self.root)['runs']
        self.assertEqual(runs[0]['name'], 'run-99-1789591874492654300.json')
        self.assertIsNone(runs[0]['started_at_unix_ms'])


if __name__ == '__main__':
    unittest.main()
