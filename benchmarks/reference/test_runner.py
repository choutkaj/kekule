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
    def test_original_failure_is_not_replaced_by_empty_collection(self):
        adapter = (SimpleNamespace(evaluate=lambda *args: {'status':'reference_error', 'message':'specific DSSP failure', 'residues':[]}), {}, {'tool':'biopython','version':'test'})
        with patch.object(runner, 'adapter', return_value=adapter):
            value = runner.run({'feature':'bio.secondary-structure.dssp','inputs':[{'path':'input.cif','text':'data_x'}]})
        self.assertIn('specific DSSP failure', value['results'][0]['message'])

    def test_reference_errors_are_isolated_and_later_inputs_still_run(self):
        calls = []
        def evaluate(feature, source, dependencies):
            calls.append(source.read_text())
            if calls[-1] == 'bad':
                raise ValueError('injected reference failure')
            return {'blocks':[{'value':calls[-1]}]}
        adapter = (SimpleNamespace(evaluate=evaluate), {}, {'tool':'biopython','version':'test'})
        with patch.object(runner, 'adapter', return_value=adapter):
            result = runner.run({'feature':'io.mmcif.parse','inputs':[
                {'path':'input.cif','text':text} for text in ('first','bad','last')]})
        self.assertEqual(calls,['first','bad','last'])
        self.assertEqual([r['status'] for r in result['results']],['ok','error','ok'])
        self.assertEqual(result['results'][2]['value'],{'blocks':[{'value':'last'}]})
        self.assertEqual(set(result),{'reference','results','time_ms'})

    def test_matching_failures_empty_records_and_nonfinite_values_are_not_goldens(self):
        for value in (None, {}, {'status':'unsupported'},{'blocks':[]},{'blocks':[{'number':float('nan')}]}, {'blocks':[{'status':'unsupported'}]}):
            adapter = (SimpleNamespace(evaluate=lambda *args:value), {}, {'tool':'biopython','version':'test'})
            with patch.object(runner,'adapter',return_value=adapter):
                result=runner.run({'feature':'io.mmcif.parse','inputs':[{'path':'input.cif','text':'source'}]})
            self.assertEqual(result['results'][0]['status'],'error')

    def test_implementation_outputs_cannot_enter_reference_request(self):
        with self.assertRaisesRegex(ValueError,'only feature'):
            runner.run({'feature':'io.smiles.parse','inputs':[],'actual':{'answer':1}})

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
