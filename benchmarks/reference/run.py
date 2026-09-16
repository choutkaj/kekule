"""Run an independent reference on every input and return its values or errors.

The reference request contains source data only. No Kekule output is ever used
as an expectation. A separate request can ask RDKit to read writer output.
"""
from __future__ import annotations
import argparse
import io
import json
import sys
import tempfile
import time
from pathlib import Path


# Support direct script execution and unittest discovery without dynamic loaders.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


class MemoryInput:
    def __init__(self, path, text):
        self.path = Path(path)
        self.text = text
        self.suffix = self.path.suffix
        self.name = self.path.name

    def read_text(self, *args, **kwargs):
        return self.text

    def read_bytes(self):
        return self.text.encode('utf-8')

    def open(self, mode='r', *args, **kwargs):
        return io.BytesIO(self.read_bytes()) if 'b' in mode else io.StringIO(self.text)


def adapter(feature):
    engine = 'biopython' if feature in {'io.mmcif.parse', 'bio.secondary-structure.dssp'} else 'rdkit'
    if engine == 'biopython':
        from reference.biopython import run_feature as module
    else:
        from reference.rdkit import run_feature as module
    dependencies = module.import_biopython() if engine == 'biopython' else module.import_rdkit()
    version = dependencies['version']
    if feature == 'bio.secondary-structure.dssp':
        version = module.dssp_reference(version)['version']
    return module, dependencies, {'tool': engine, 'version': version}


def failure(value):
    if isinstance(value, dict):
        if 'status' in value and value['status'] != 'ok':
            return str(value.get('message', value['status']))
        if 'records' in value and (not isinstance(value['records'], list) or not value['records']):
            return 'reference returned no records'
        return next((error for child in value.values() if (error := failure(child))), None)
    if isinstance(value, list):
        return next((error for child in value if (error := failure(child))), None)
    return None


def run(request):
    feature = request['feature']
    if set(request) == {'feature', 'describe'} and request['describe'] is True:
        _, _, reference = adapter(feature)
        return {'reference': reference, 'results': [], 'time_ms': 0}
    written = 'written' in request
    if set(request) != {'feature', 'written' if written else 'inputs'}:
        raise ValueError('request must contain only feature and inputs (or written)')
    items = request['written' if written else 'inputs']
    module, dependencies, reference = adapter(feature)
    if reference['tool'] == 'rdkit':
        from reference.rdkit import molecule
    writer = feature.endswith('.write') or feature in ('io.smiles.canonical', 'io.smiles.isomeric')
    results = []
    elapsed = 0
    with tempfile.TemporaryDirectory(prefix='kekule-reference-') as directory:
        for index, item in enumerate(items):
            start = time.perf_counter_ns()
            try:
                if written:
                    if not writer or reference['tool'] != 'rdkit':
                        raise ValueError('no independent writer reader for this feature')
                    if item['status'] != 'ok':
                        raise ValueError(item['message'])
                    value = molecule.read_written(feature, item['value'], MemoryInput)
                else:
                    if set(item) != {'path', 'text'}:
                        raise ValueError('input must contain only path and text')
                    if reference['tool'] == 'biopython':
                        source = Path(directory) / f'{index}.cif'
                        source.write_bytes(item['text'].encode('utf-8'))
                        value = module.evaluate(feature, source, dependencies)
                    else:
                        source = MemoryInput(item['path'], item['text'])
                        if writer:
                            value = molecule.writer_value(feature, source)
                        elif feature.startswith(('io.smiles.', 'io.mol.', 'io.sdf.')) or feature in ('stereo.representation', 'stereo.perception'):
                            value = molecule.evaluate(feature, source)
                        else:
                            value = module.evaluate(feature, source, dependencies)
                error = failure(value)
                if error:
                    raise ValueError(error)
                collection = {'io.mmcif.parse':'blocks', 'bio.secondary-structure.dssp':'residues'}.get(feature, 'records')
                if not isinstance(value, dict) or not isinstance(value.get(collection), list) or not value[collection]:
                    raise ValueError(f'missing or empty {collection}')
                # Reject NaN/infinity at the individual-case boundary.
                json.dumps(value, allow_nan=False)
                results.append({'status': 'ok', 'value': value})
            except Exception as error:
                results.append({'status': 'error', 'message': f'{type(error).__name__}: {error}'})
            elapsed += time.perf_counter_ns() - start
    return {'reference': reference, 'results': results, 'time_ms': elapsed / 1e6}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--feature')
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.feature or args.input or args.output:
        if not all((args.feature, args.input, args.output)):
            parser.error('--feature, --input and --output must be supplied together')
        request = {'feature': args.feature, 'inputs': [{'path': str(args.input), 'text': args.input.read_bytes().decode('utf-8')}]}
        result = run(request)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x', encoding='utf-8') as output:
            json.dump(result, output, indent=2, allow_nan=False)
    else:
        sys.stdout.write(json.dumps(run(json.load(sys.stdin)), allow_nan=False))


if __name__ == '__main__':
    main()
