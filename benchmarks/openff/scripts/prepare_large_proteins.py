"""Freeze a panel of large protein chains before observing native parameterization.

Run in environment.yml. The chains come from the same locked `pdb-1000` corpus
and preparation as the robustness panel, selected in the same ranked order, but
with 300-1,000 residues and at least 4,097 atoms including hydrogens, beyond the
former 4,096-atom NAGL limit. Selection never calls kekule.

    python benchmarks/openff/scripts/prepare_large_proteins.py --output benchmarks/openff/data/large-inputs.json.gz
"""
import argparse
import gzip
import importlib.metadata
import json
from pathlib import Path

from prepare_robustness import CORPORA, proteins, sha

COUNT = 4
RESIDUES = (300, 1000)
MIN_ATOMS = 4097


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    records, attempts = proteins(count=COUNT, residue_range=RESIDUES, min_atoms=MIN_ATOMS)
    packages = ['rdkit', 'openff-toolkit', 'openmm', 'numpy']
    report = dict(schema=1, selection='rosemary-large-proteins-v1',
                  criteria=dict(count=COUNT, residues=list(RESIDUES), min_atoms=MIN_ATOMS),
                  versions={p: importlib.metadata.version(p) for p in packages},
                  corpus_locks={'pdb-1000': sha(CORPORA / 'pdb-1000' / 'sources.lock.json')},
                  protein_preparation_attempts=attempts, records=records)
    args.output.write_bytes(gzip.compress((json.dumps(report, allow_nan=False) + '\n').encode(), mtime=0))
    print(json.dumps(dict(output=str(args.output), sha256=sha(args.output), cases=len(records),
                          atoms=[r['atoms'] for r in records])))


if __name__ == '__main__':
    # Large chains exceed the default Windows thread stack during preparation.
    import sys
    import threading
    sys.setrecursionlimit(100_000)
    threading.stack_size(256 * 1024 * 1024 - 4096)
    worker = threading.Thread(target=main)
    worker.start()
    worker.join()
