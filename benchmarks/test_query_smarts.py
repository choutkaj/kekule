"""Integrity of the standard SMARTS corpus and independent target selection."""
import hashlib
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).parent


class QuerySmartsTests(unittest.TestCase):
    def test_standard_corpus_preserves_every_original_query_and_source_byte(self):
        original = ROOT / 'smarts-fixtures/rdkit-queries'
        standard = ROOT / 'corpora/rdkit-queries'
        lock = json.loads((standard / 'sources.lock.json').read_text())
        self.assertEqual(len(lock['entries']), 518)
        for relative in ['sources.lock.json'] + [p['path'] for p in lock['packs'] + lock['upstream']]:
            self.assertEqual((original / relative).read_bytes(), (standard / relative).read_bytes())
        for item in lock['packs'] + lock['upstream']:
            self.assertEqual(hashlib.sha256((standard / item['path']).read_bytes()).hexdigest(), item['sha256'])

    def test_targets_are_pinned_to_external_sources_and_preselected_by_identity(self):
        contract = json.loads((ROOT / 'query-smarts.json').read_text())
        self.assertEqual(len(contract['targets']), 16)
        self.assertEqual(len({t['id'] for t in contract['targets']}), 16)
        lock = json.loads((ROOT / 'corpora/pubchem-100k/sources.lock.json').read_text())
        selected = sorted((entry['id'] for entry in lock['entries']),
                          key=lambda id: hashlib.sha256(('query-smarts-targets-v2:' + id).encode()).hexdigest())[:8]
        self.assertEqual([t['id'] for t in contract['targets'][8:]], ['pubchem:' + id for id in selected])
        packs = {pack['path']: pack for pack in lock['packs']}
        for target in contract['targets']:
            source = ROOT / target['source']
            if target['id'].startswith('pubchem:'):
                pack = packs[target['source'].removeprefix('corpora/pubchem-100k/')]
                self.assertEqual(pack['sha256'], target['source_sha256'])
                self.assertEqual(pack['members'][target['source_record']], target['id'].split(':')[1])
            # Full PubChem packs are optional local data; supplied files must match.
            if source.exists():
                self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), target['source_sha256'])
                self.assertEqual(source.read_text().splitlines()[target['source_record']].split()[0], target['smiles'])


if __name__ == '__main__':
    unittest.main()
