"""Replay frozen independent model references through the complete Rust parameterizer.

Rust's models integration test additionally checks every valence/vdW parameter.
This observer retains native outputs and reports charge/feature error magnitudes.
"""
import argparse
import gzip
import json
import subprocess
from pathlib import Path
from audit import digest


def compare(binary, bundles, reference, output):
    if output.exists(): raise FileExistsError(output)
    report = json.loads(gzip.decompress(reference.read_bytes()))
    results = []
    for model in report['models']:
        directory = bundles / model['name'].removesuffix('.pt')
        for file, sha in model['bundle_sha256'].items():
            if digest(directory / file) != sha: raise ValueError(f'Bundle fingerprint mismatch: {file}')
        records = []
        maximum = dict(system_charge_e=0.0, assignment_charge_e=0.0, inference_charge_e=0.0, feature=0.0)
        process = subprocess.Popen([str(binary.resolve()), str(directory.resolve()), str((directory / 'force-field.offxml').resolve())], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        try:
            for r in model['records']:
                process.stdin.write(json.dumps(dict(smiles=r['smiles'])) + '\n'); process.stdin.flush()
                native = json.loads(process.stdout.readline())
                errors = []
                if native['status'] != 'ok':
                    errors.append(native)
                else:
                    maps = native['maps']
                    if sorted(maps) != list(range(1, len(r['system']['charges']) + 1)): raise ValueError('Atom map mismatch')
                    order = [maps.index(i) for i in sorted(maps)]
                    def check(label, observed, expected, atol):
                        if len(observed) != len(expected): raise ValueError(f'{label} shape mismatch')
                        delta = max(abs(observed[i] - expected[j]) for j,i in enumerate(order))
                        maximum[label] = max(maximum[label], delta)
                        if delta > atol: errors.append(dict(field=label, max_abs_error=delta))
                    check('system_charge_e', native['system']['charges'], r['system']['charges'], 5e-5)
                    for field, observed in [('assignment', 'charges'), ('inference','inference')]:
                        if native[observed]['status'] != r[field]['status']:
                            errors.append(dict(field=field, expected_status=r[field]['status'], actual=native[observed]))
                        elif r[field]['status'] == 'ok':
                            check(field + '_charge_e', native[observed]['values'], r[field]['charges'], 1e-6 if field == 'inference' else 5e-5)
                    if r['inference']['status'] == 'ok':
                        expected = r['inference']['features']
                        actual = native['features']['values']
                        if any(len(a) != len(b) for a,b in zip(actual,expected)): raise ValueError('Feature width mismatch')
                        for column in range(len(expected[0])):
                            check('feature', [a[column] for a in actual], [b[column] for b in expected], 1e-6)
                records.append(dict(id=r['id'], reverse=r['reverse'], errors=errors, native=native))
        finally:
            process.stdin.close()
            process.wait()
        if process.returncode: raise RuntimeError(f'Native observer exited {process.returncode}')
        result = dict(name=model['name'], passed=sum(not r['errors'] for r in records), cases=len(records), maximum=maximum, records=records)
        print(json.dumps({k:v for k,v in result.items() if k != 'records'}), flush=True)
        results.append(result)
    payload = dict(schema=1, reference_sha256=digest(reference), binary_sha256=digest(binary), models=results)
    output.write_bytes(gzip.compress((json.dumps(payload, allow_nan=False) + '\n').encode(), mtime=0))
    return int(any(m['passed'] != m['cases'] for m in results))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--bundles', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(compare(args.binary, args.bundles, args.reference, args.output))
