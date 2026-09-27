"""Time every molecule in the frozen panel using the release Rust timing binary."""
import argparse
from datetime import datetime, timezone
import json
import hashlib
import platform
import subprocess
from pathlib import Path
from robustness import read, digest

from paths import HERE
ROOT = HERE.parents[1]


def run(args):
    if args.output.exists(): raise FileExistsError(args.output)
    cases = read(args.inputs)['records']
    requests = [dict(id=r['id'], smiles=r['smiles']) for r in cases]
    result = subprocess.run([str(args.binary.resolve()), str(args.model.resolve())],
                            input=''.join(json.dumps(r)+'\n' for r in requests),
                            text=True, capture_output=True, check=True)
    rows = [json.loads(line) for line in result.stdout.splitlines()]
    loading, *records = rows
    if [r['id'] for r in records] != [r['id'] for r in cases]: raise ValueError('Incomplete timing output')
    if [r['atoms'] for r in records] != [r['atoms'] for r in cases]: raise ValueError('Timing atom counts differ from frozen inputs')
    repetitions = dict(features_ms=5, inference_including_features_ms=5, assign_charges_ms=5, full_parameterization_ms=3)
    for r in records:
        for key, count in repetitions.items():
            if len(r[key]) != count or not all(0 < t < float('inf') for t in r[key]): raise ValueError('Invalid timing repetitions')
    sources = [ROOT/'benchmarks/src/bin/openff_timing.rs', ROOT/'Cargo.lock', *sorted((ROOT/'crates/kekule-openff/src').rglob('*.rs'))]
    report = dict(schema=1, recorded=datetime.now(timezone.utc).isoformat(), cpu=args.cpu,
                  platform=platform.platform(), rustc=subprocess.check_output(['rustc','--version'],text=True).strip(),
                  profile='release; default target CPU; default ndarray CPU backend',
                  method='All 110 frozen inputs, parsed before timing; model and force field resident. One warm-up per operation; five feature/inference/assignment repetitions and three complete parameterizations. Operations measured independently. No energy calculation or GPU timing.',
                  warmups_per_operation=1, repetitions=repetitions, model_load_ms=loading['model_load_ms'],
                  binary_sha256=digest(args.binary), inputs_sha256=digest(args.inputs),
                  model_sha256={name:digest(args.model/name) for name in ('model.json','weights.bin')},
                  source_sha256_lf={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_text(encoding='utf-8').encode()).hexdigest() for p in sources},
                  records=records)
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(f"Timed {len(records)} molecules; model load {loading['model_load_ms']:.3f} ms")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, default=HERE/'data/inputs.json.gz')
    parser.add_argument('--binary', type=Path, default=ROOT/'target/release/openff_timing.exe')
    parser.add_argument('--model', type=Path, default=ROOT/'target/openff-models/openff-gnn-am1bcc-1.0.0')
    parser.add_argument('--cpu', required=True)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
