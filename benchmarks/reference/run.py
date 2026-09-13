"""Reference adapter: JSON on stdin/stdout, or write a reviewable candidate file.

Run inside the pinned RDKit or Biopython/DSSP environment. Nothing updates a
tracked reference. Imports, preparation and JSON transport are outside timing.
"""
from __future__ import annotations

import argparse
import importlib.util
import hashlib
import io
import json
import sys
import tempfile
import time
from pathlib import Path


class MemoryInput:
    """Path-shaped input for existing reference algorithms, with preloaded bytes."""
    def __init__(self, path: str, text: str):
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
    path = Path(__file__).parent / engine / 'run_feature.py'
    spec = importlib.util.spec_from_file_location('reference_feature', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if feature not in module.SUPPORTED_FEATURES:
        raise ValueError(f'no independent adapter for {feature}')
    dependencies = module.import_biopython() if engine == 'biopython' else module.import_rdkit()
    reference = module.dssp_reference(dependencies['version']) if feature == 'bio.secondary-structure.dssp' else {'tool': engine, 'version': dependencies['version']}
    return module, dependencies, reference


def run(request):
    if 'samples' in request:
        raise ValueError('samples is no longer supported; each input is evaluated once')
    feature = request['feature']
    module, dependencies, reference = adapter(feature)
    with tempfile.TemporaryDirectory(prefix='kekule-reference-') as directory:
        inputs = []
        for index, item in enumerate(request['inputs']):
            # Bio.PDB and DSSP require real paths; stage them before timing.
            if feature in {'io.mmcif.parse', 'bio.secondary-structure.dssp'}:
                path = Path(directory) / f'{index}.cif'
                path.write_text(item['text'], encoding='utf-8', newline='')
                inputs.append(path)
            else:
                inputs.append(MemoryInput(item['path'], item['text']))

        evidence = []
        expected = []
        elapsed_ns = 0
        for source in inputs:
            record_evidence = []
            start = time.perf_counter_ns()
            if reference['tool'] == 'rdkit':
                value = module.evaluate(feature, source, dependencies, record_evidence)
            else:
                value = module.evaluate(feature, source, dependencies)
            elapsed_ns += time.perf_counter_ns() - start
            expected.append(value)
            if reference['tool'] == 'rdkit':
                evidence.append(record_evidence)
        return {
            'reference': reference, 'expected': expected,
            'inputs': [{'path': item['path'], 'sha256': hashlib.sha256(item['text'].encode('utf-8')).hexdigest()} for item in request['inputs']],
            'reference_evidence': evidence,
            'time_ms': elapsed_ns / 1e6,
            'scope': 'one serial evaluation per input: parse + feature + result/evidence materialization; excludes result drop, Python startup/imports/JSON transport; Bio.PDB/DSSP includes staged-file I/O and DSSP subprocesses; no warmup or repetitions',
            'comparable_speed_ratio': False,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--feature')
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.feature or args.input or args.output:
        if not all([args.feature, args.input, args.output]):
            parser.error('--feature, --input and --output must be supplied together')
        request = {'feature': args.feature, 'inputs': [{'path': str(args.input), 'text': args.input.read_bytes().decode('utf-8')}]}
        result = run(request)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x', encoding='utf-8') as output:
            json.dump(result, output, indent=2, allow_nan=False)
    else:
        json.dump(run(json.load(sys.stdin)), sys.stdout, allow_nan=False)


if __name__ == '__main__':
    main()
