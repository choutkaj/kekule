import importlib.util
import io
import tarfile
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


data = load(Path(__file__).parents[1] / 'data.py', 'benchmark_data')
runner = load(Path(__file__).parent / 'run.py', 'benchmark_reference')


class RunnerTests(unittest.TestCase):
    def test_reference_evaluates_each_input_once_and_returns_timed_outputs(self):
        for engine, feature in [('rdkit', 'io.smiles.parse'), ('biopython', 'io.mmcif.parse')]:
            with self.subTest(engine=engine):
                calls = []

                def evaluate(feature, source, dependencies, evidence=None):
                    text = source.read_text()
                    calls.append(text)
                    if evidence is not None:
                        evidence.append({'source': text})
                    return {'value': text, 'evaluation': len(calls)}

                adapter = (SimpleNamespace(evaluate=evaluate), {}, {'tool': engine, 'version': 'test'})
                with patch.object(runner, 'adapter', return_value=adapter), patch.object(
                    runner.time, 'perf_counter_ns', side_effect=[100, 300, 500, 800]
                ):
                    result = runner.run({'feature': feature, 'inputs': [
                        {'path': 'first.smi', 'text': 'first'},
                        {'path': 'second.smi', 'text': 'second'},
                    ]})
                self.assertEqual(calls, ['first', 'second'])
                self.assertEqual(result['expected'], [
                    {'value': 'first', 'evaluation': 1},
                    {'value': 'second', 'evaluation': 2},
                ])
                self.assertEqual(result['time_ms'], 0.0005)
                self.assertEqual(result['reference_evidence'],
                    [[{'source': 'first'}], [{'source': 'second'}]] if engine == 'rdkit' else [])

    def test_repetition_option_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'no longer supported'):
            runner.run({'samples': 1})

    def test_cached_input_preserves_bytes_and_open_modes(self):
        value = runner.MemoryInput('input.sdf', '\r\nexample\n')
        self.assertEqual(value.open('rb').read(), b'\r\nexample\n')
        self.assertEqual(value.open().read(), '\r\nexample\n')

    def test_input_paths_cannot_escape(self):
        for name in ('../file', '/data/file', 'data/../../file', 'data\\file', 'C:/file'):
            with self.assertRaises(ValueError):
                data.safe_path(Path('corpus'), name)

    def test_bundle_checks_bytes_and_never_overwrites_different_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / 'bundle.tar.gz'
            with tarfile.open(archive, 'w:gz') as bundle:
                info = tarfile.TarInfo('data/input.smi')
                info.size = 3
                bundle.addfile(info, io.BytesIO(b'CC\n'))
            pinned = {'data/input.smi': data.hashlib.sha256(b'CC\n').hexdigest()}
            with self.assertRaises(ValueError):
                data.unpack(archive, 'wrong', root / 'out', pinned)
            data.unpack(archive, data.sha256(archive), root / 'out', pinned)
            target = root / 'out/data/input.smi'
            self.assertEqual(target.read_bytes(), b'CC\n')
            target.write_bytes(b'user data')
            with self.assertRaises(ValueError):
                data.unpack(archive, data.sha256(archive), root / 'out', pinned)
            self.assertEqual(target.read_bytes(), b'user data')

    def test_duplicate_archive_members_fail_before_installation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / 'bundle.tar'
            with tarfile.open(archive, 'w') as bundle:
                for _ in range(2):
                    info = tarfile.TarInfo('data/input.smi')
                    info.size = 3
                    bundle.addfile(info, io.BytesIO(b'CC\n'))
            with self.assertRaises(ValueError):
                data.unpack(archive, data.sha256(archive), root / 'out', {'data/input.smi': data.hashlib.sha256(b'CC\n').hexdigest()})
            self.assertFalse((root / 'out').exists())


if __name__ == '__main__':
    unittest.main()
