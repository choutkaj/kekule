"""Verify, bundle and restore the existing pinned benchmark inputs.

No selection, downloads, chemistry conversion or reference regeneration occurs.
Bundles can be kept on ordinary artifact storage; unpack requires their SHA-256.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tarfile
import tempfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent / 'corpora'
DATASETS = ('pubchem-100k', 'enamine-diversity', 'pl-rex', 'pdb-1000', 'smoke')


def sha256(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def safe_path(root, name):
    path = PurePosixPath(name)
    if '\\' in name or ':' in name or path.is_absolute() or '..' in path.parts or not path.parts or path.parts[0] != 'data':
        raise ValueError(f'unsafe input path: {name}')
    result = root.joinpath(*path.parts)
    if not result.resolve().is_relative_to(root.resolve()):
        raise ValueError(f'input escapes dataset: {name}')
    return result


def files(dataset):
    lock = json.loads((ROOT / dataset / 'sources.lock.json').read_text())
    pinned = {}
    for entry in lock['entries']:
        for item in entry['files']:
            pinned[item['path']] = item['sha256']
    for item in lock.get('packs', []):
        pinned[item['path']] = item['sha256']
    return pinned


def verify(root, pinned):
    failures = []
    for name, expected in pinned.items():
        path = safe_path(root, name)
        if not path.is_file():
            failures.append(f'missing: {name}')
        elif sha256(path) != expected:
            failures.append(f'checksum mismatch: {name}')
    if failures:
        raise ValueError('\n'.join(failures))


def unpack(archive, expected, destination, pinned):
    if sha256(archive) != expected:
        raise ValueError('archive checksum mismatch')
    with tempfile.TemporaryDirectory(prefix='kekule-data-') as directory:
        staged = Path(directory)
        with tarfile.open(archive, 'r:*') as bundle:
            names = set()
            for member in bundle:
                if not member.isfile() or member.name not in pinned or member.name in names:
                    raise ValueError(f'unexpected or duplicate archive member: {member.name}')
                names.add(member.name)
                output = safe_path(staged, member.name)
                output.parent.mkdir(parents=True, exist_ok=True)
                with bundle.extractfile(member) as source, output.open('xb') as target:
                    shutil.copyfileobj(source, target)
        verify(staged, pinned)
        # Validate all existing destinations before installing any files.
        for name, checksum in pinned.items():
            target = safe_path(destination, name)
            if target.exists() and (not target.is_file() or sha256(target) != checksum):
                raise ValueError(f'refusing to overwrite different local input: {name}')
        for name in pinned:
            target = safe_path(destination, name)
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                with safe_path(staged, name).open('rb') as source, target.open('xb') as output:
                    shutil.copyfileobj(source, output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('verify', 'pack', 'unpack'))
    parser.add_argument('--dataset', choices=DATASETS, required=True)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--sha256')
    args = parser.parse_args()
    root, pinned = ROOT / args.dataset, files(args.dataset)
    if args.command == 'verify':
        verify(root, pinned)
        print(f'{args.dataset}: verified {len(pinned)} input files')
    elif args.command == 'pack':
        if args.archive is None:
            parser.error('pack requires --archive')
        verify(root, pinned)
        args.archive.parent.mkdir(parents=True, exist_ok=True)
        with args.archive.open('xb') as raw, tarfile.open(fileobj=raw, mode='w:gz') as bundle:
            for name in sorted(pinned):
                bundle.add(safe_path(root, name), arcname=name, recursive=False)
        print(f'{sha256(args.archive)}  {args.archive}')
    else:
        if args.archive is None or args.sha256 is None:
            parser.error('unpack requires --archive and --sha256')
        unpack(args.archive, args.sha256, root, pinned)
        print(f'{args.dataset}: restored and verified {len(pinned)} input files')


if __name__ == '__main__':
    main()
