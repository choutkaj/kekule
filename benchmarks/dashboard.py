"""Refresh local benchmark history or build a fixed report snapshot.

Usage: python benchmarks/dashboard.py [REPORT.json ...] [--runs-dir DIRECTORY]
Only summary counts and provenance are embedded. No inputs, golden payloads,
per-case observations, local paths or raw error messages are exported.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
import sys

ROOT = Path(__file__).resolve().parent
COUNTS = (
    'source_ids', 'cases', 'agrees', 'exact_agrees', 'disagrees', 'errors',
    'not_applicable', 'input_errors', 'reference_errors', 'kekule_errors',
    'writer_validation_errors', 'observation_errors', 'structural_differences',
    'numerical_differences',
)
HASHES = ('contract_sha256', 'input_lock_sha256', 'sha256')
CORPORA = {
    'smoke': 'Small mixed set from PubChem, RDKit test data and the RCSB PDB.',
    'rdkit-queries': 'Complete query rows from three pinned RDKit functional-group/reactivity tables; includes unsupported grammar.',
    'rdkit-structures': 'Pinned RDKit input variants for six medicinal compounds with axial stereo and ChEBI V3000 structures; includes invalid stereo variants.',
    'enamine-diversity': 'Enamine Discovery Diversity Set 50. All supplied records, including CXSMILES.',
    'pdb-1000': 'RCSB PDB structures: proteins, nucleic acids, complexes and multi-model entries.',
    'pl-rex': 'Primary refined ligands from PL-REX 1.0.1.',
    'pubchem-100k': 'PubChem sample preselected for V2000, size and RDKit success.',
}
FORMATS = {'.sdf': 'SDF', '.mol': 'MOL', '.cif': 'mmCIF',
           '.txt': 'SMILES', '.smi': 'SMILES', '.smiles': 'SMILES',
           '.smarts': 'SMARTS', '.sma': 'SMARTS'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def unique_object(items):
    result = {}
    for key, value in items:
        require(key not in result, f'duplicate JSON key: {key}')
        result[key] = value
    return result


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'), object_pairs_hook=unique_object)


def text_hash(path):
    return hashlib.sha256(path.read_text(encoding='utf-8').encode()).hexdigest()


def count(value, name):
    require(type(value) is int and 0 <= value <= 2**53 - 1,
            f'{name} must be a nonnegative, exactly representable integer')
    return value


def digest(value, name, optional=False):
    if optional and value is None:
        return None
    require(isinstance(value, str) and re.fullmatch(r'[0-9a-f]{64}', value),
            f'invalid {name}')
    return value


def golden_metadata(raw):
    require(isinstance(raw, dict) and raw.get('schema') == 2,
            'a schema-2 golden manifest is required')
    result = {key: digest(raw.get(key), key) for key in HASHES}
    for key in ('dataset', 'feature'):
        value = raw.get(key)
        require(isinstance(value, str) and re.fullmatch(r'[a-z0-9][a-z0-9.-]*', value),
                f'invalid golden {key}')
        result[key] = value
    result['cases'] = count(raw.get('cases'), 'golden cases')
    result['reference_code_sha256'] = digest(
        raw.get('reference_code_sha256'), 'reference code hash', optional=True)
    reference = raw.get('reference')
    require(reference is None or (isinstance(reference, dict) and all(
        isinstance(reference.get(key), str) and reference[key] for key in ('tool', 'version'))),
        'invalid reference identity')
    result['reference'] = ({key: reference[key] for key in ('tool', 'version')}
                           if reference else None)
    return result


def catalogue(root=ROOT):
    datasets = []
    for path in sorted((root / 'corpora').glob('*/sources.lock.json')):
        lock = read_json(path)
        require(lock['corpus_id'] == path.parent.name, 'source lock identity mismatch')
        ids = [entry['id'] for entry in lock['entries']]
        require(len(ids) == len(set(ids)) and ids, 'empty or duplicate source membership')
        files = {item['path'] for entry in lock['entries'] for item in entry.get('files', [])}
        files.update(pack['path'] for pack in lock.get('packs', []))
        formats = sorted({FORMATS.get(Path(name).suffix.lower(), Path(name).suffix[1:].upper())
                          for name in files})
        datasets.append({'id': path.parent.name, 'source_ids': len(ids), 'formats': formats,
                         'description': CORPORA.get(path.parent.name, ''),
                         'input_lock_sha256': text_hash(path)})
    datasets.sort(key=lambda item: (item['id'] != 'smoke', item['id']))
    goldens = {}
    for path in sorted((root / 'goldens').glob('*/*.jsonl.meta.json')):
        golden = golden_metadata(read_json(path))
        key = f"{golden['dataset']}/{golden['feature']}"
        require(path.parent.name == golden['dataset'] and
                path.name == golden['feature'] + '.jsonl.meta.json',
                f'manifest identity mismatch: {path}')
        require(key not in goldens, f'duplicate manifest: {key}')
        goldens[key] = golden
    features = sorted({item['feature'] for item in goldens.values()})
    expected = {f"{dataset['id']}/{feature}" for dataset in datasets for feature in features}
    require(goldens and set(goldens) == expected, 'incomplete reference catalogue')
    return {'datasets': datasets, 'features': features, 'goldens': goldens,
            'contract_sha256': text_hash(root / 'contract.json')}


def load_report(path, catalog):
    raw_bytes = path.read_bytes()
    raw = json.loads(raw_bytes, object_pairs_hook=unique_object)
    require(isinstance(raw, dict) and raw.get('schema') == 2 and raw.get('mode') == 'compare',
            'only schema-2 comparison reports are supported (not reference generation)')
    require(type(raw.get('complete')) is bool and type(raw.get('passed')) is bool,
            'missing report completion/pass status')
    identity = raw.get('implementation')
    require(isinstance(identity, dict), 'missing implementation provenance')
    implementation = {key: digest(identity.get(key), key, optional=key == 'working_tree_status_sha256')
                      for key in ('contract_sha256', 'executable_sha256',
                                  'reference_code_sha256', 'working_tree_status_sha256')}
    revision = identity.get('revision')
    require(revision is None or (isinstance(revision, str) and
            re.fullmatch(r'[0-9a-f]{40,64}', revision)), 'invalid implementation revision')
    require(identity.get('dirty') is None or type(identity['dirty']) is bool,
            'invalid worktree status')
    implementation.update(revision=revision, dirty=identity.get('dirty'))
    require(isinstance(raw.get('results'), list), 'missing report results')
    rows, seen = [], set()
    datasets = {item['id']: item for item in catalog['datasets']}
    for item in raw['results']:
        require(isinstance(item, dict), 'report results must be objects')
        golden = golden_metadata(item.get('golden'))
        dataset, feature = item.get('dataset'), item.get('feature')
        key = f'{dataset}/{feature}'
        require(key in catalog['goldens'] and key not in seen,
                f'unknown or duplicate dataset/feature row: {key}')
        seen.add(key)
        require((dataset, feature) == (golden['dataset'], golden['feature']),
                f'golden identity mismatch: {key}')
        row = {name: count(item.get(name), name) for name in COUNTS}
        require(row['cases'] == row['agrees'] + row['disagrees'] + row['errors'],
                f'outcome counts do not sum to cases: {key}')
        require(row['exact_agrees'] <= row['agrees'], f'exact agreements exceed agreements: {key}')
        require(row['source_ids'] > 0, f'invalid selection counts: {key}')
        for field in ('input_errors', 'reference_errors', 'kekule_errors',
                      'writer_validation_errors', 'observation_errors'):
            require(row[field] <= row['errors'], f'{field} exceeds error cases: {key}')
        # A selection may extend beyond a sampled reference archive. Missing
        # references and unreadable inputs remain errors; unavailable formats
        # are counted independently of reference lookup. Other applicable cases
        # still require stored reference records. These aggregate error counts
        # include stored reference failures and can overlap, so give a bound.
        unreferenced_bound = min(row['errors'], row['reference_errors'] + row['input_errors'])
        require(row['cases'] - unreferenced_bound <= golden['cases'],
                f'applicable cases exceed reference coverage: {key}')
        for field in ('kekule_ms', 'reference_ms'):
            value = item.get(field)
            require(type(value) in (int, float) and math.isfinite(value) and value >= 0,
                    f'invalid workflow timing: {key}')
            row[field] = value
        require(golden['contract_sha256'] == implementation['contract_sha256'],
                f'report/manifest comparison contract mismatch: {key}')
        current = (golden == catalog['goldens'][key]
                   and golden['input_lock_sha256'] == datasets[dataset]['input_lock_sha256']
                   and implementation['contract_sha256'] == catalog['contract_sha256'])
        if current:
            require(row['source_ids'] <= datasets[dataset]['source_ids'],
                    f'selection exceeds source membership: {key}')
            full = row['source_ids'] == datasets[dataset]['source_ids']
            # Full describes source selection, not reference coverage or run
            # completion. One input failure can replace many unreadable records.
            if full and raw['complete'] and row['input_errors'] == 0:
                require(row['cases'] + row['not_applicable'] >= golden['cases'],
                        f'full selection does not account for every golden record: {key}')
        row.update(dataset=dataset, feature=feature, golden=golden,
                   coverage=('full' if full else 'sampled') if current else 'stale')
        rows.append(row)
    passed = (raw['complete'] and bool(rows) and sum(row['cases'] for row in rows) > 0
              and all(row['errors'] == row['disagrees'] == 0 for row in rows))
    require(raw['passed'] == passed, 'pass status is inconsistent with measured outcomes')
    for field in COUNTS:
        count(sum(row[field] for row in rows), f'total {field}')
    started = raw.get('started_at_unix_ms')
    if started is not None:
        require(count(started, 'run start time') <= 8_640_000_000_000_000,
                'run start time exceeds the supported date range')
    return {'name': path.name, 'sha256': hashlib.sha256(raw_bytes).hexdigest(),
            'started_at_unix_ms': started,
            'complete': raw['complete'], 'passed': raw['passed'],
            'implementation': implementation, 'results': rows}


def payload(reports, root=ROOT):
    catalog = catalogue(root)
    runs = [load_report(Path(path), catalog) for path in reports]
    require(runs, 'select at least one comparison report')
    require(len({run['sha256'] for run in runs}) == len(runs), 'duplicate comparison report')
    return {'catalog': catalog, 'runs': runs}


def encoded(data):
    data = json.dumps(data, separators=(',', ':'), allow_nan=False)
    # A report string must never be able to close the embedded JSON script element.
    return data.replace('&', r'\u0026').replace('<', r'\u003c').replace('>', r'\u003e')


def page(data):
    template = (ROOT / 'dashboard' / 'index.html').read_text(encoding='utf-8')
    return (template.replace('/* DASHBOARD_CSS */', (ROOT / 'dashboard' / 'style.css').read_text(encoding='utf-8'))
            .replace('/* DASHBOARD_JS */', (ROOT / 'dashboard' / 'app.js').read_text(encoding='utf-8'))
            .replace('DASHBOARD_LIVE_SOURCE', ' file:' if data.get('live') else '')
            .replace('DASHBOARD_DATA', encoded(data)))


def render(reports, root=ROOT):
    return page(payload(reports, root))


def discovered(runs_dir, root=ROOT):
    catalog = catalogue(root)
    runs = {}
    for path in sorted(runs_dir.glob('*.json')):
        try:
            if read_json(path).get('mode') == 'generate':
                continue
            run = load_report(path, catalog)
            runs.setdefault(run['sha256'], run)
        except (ValueError, KeyError, TypeError, AttributeError, OSError) as error:
            print(f'dashboard: skipped {path.name}: {error}', file=sys.stderr)
    def newest(run):
        # Older runner filenames recorded start time in Unix nanoseconds.
        legacy = re.fullmatch(r'run-\d+-(\d{19,})\.json', run['name'])
        started = run['started_at_unix_ms']
        if started is None:
            started = int(legacy[1]) // 1_000_000 if legacy else -1
        return started, run['name']
    data = {'catalog': catalog, 'runs': sorted(runs.values(), key=newest, reverse=True)}
    data['signature'] = hashlib.sha256(encoded(data).encode()).hexdigest()
    data['live'] = True
    return data


def write_atomic(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(dir=path.parent, suffix='.tmp')
    try:
        with os.fdopen(handle, 'w', encoding='utf-8', newline='\n') as output:
            output.write(content)
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


@contextmanager
def refresh_lock(runs_dir):
    runs_dir.mkdir(parents=True, exist_ok=True)
    with (runs_dir / '.dashboard.lock').open('a+b') as lock:
        if os.name == 'nt':
            import msvcrt
            if lock.tell() == 0:
                lock.write(b'\0')
                lock.flush()
            lock.seek(0)
            msvcrt.locking(lock.fileno(), msvcrt.LK_LOCK, 1)
        else:
            import fcntl
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        yield  # Closing the file releases the OS lock, including after exceptions.


def refresh(runs_dir, output=None, root=ROOT):
    output = output or runs_dir / 'index.html'
    with refresh_lock(runs_dir):
        data = discovered(runs_dir, root)
        write_atomic(output.with_name('dashboard-data.js'),
                     f'window.updateKekuleBenchmarks({encoded(data)});\n')
        write_atomic(output, page(data))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reports', nargs='*', type=Path)
    parser.add_argument('--runs-dir', type=Path, help='discover local reports and enable automatic page updates')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        require(not (args.reports and args.runs_dir), 'choose explicit reports or --runs-dir')
        runs_dir = args.runs_dir or ROOT / 'runs'
        output = args.output or (ROOT.parent / 'target/benchmark-dashboard/index.html'
                                 if args.reports else runs_dir / 'index.html')
        require(output.suffix.lower() == '.html', 'output must be an .html file')
        require(output.resolve() not in {path.resolve() for path in args.reports},
                'output must not replace an input report')
        if args.reports:
            write_atomic(output, render(args.reports))
        else:
            refresh(runs_dir, output)
        print(output.resolve())
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.exit(2, f'dashboard: {error}\n')


if __name__ == '__main__':
    main()
