"""Manifest parsing contracts shared by the optional reference generators."""

import contextlib
import importlib.util
import io
from pathlib import Path
import tempfile
import tomllib
import unittest
from unittest.mock import patch


def load_generator(engine):
    path = Path(__file__).parent / engine / "run_feature.py"
    spec = importlib.util.spec_from_file_location(f"{engine}_reference", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ManifestTests(unittest.TestCase):
    def test_native_toml_preserves_quoted_strings_and_selection_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "example.toml"
            source = r'''feature_id = "example"
fixtures = [
  "quoted \"record\".smi", # a TOML comment
  'second.smi',
]
'''
            for engine in ("rdkit", "biopython"):
                with self.subTest(engine=engine):
                    path.write_text(source, encoding="utf-8")
                    generator = load_generator(engine)
                    manifest = generator.read_manifest(path)
                    self.assertEqual(manifest["feature_id"], "example")
                    self.assertEqual(manifest["fixtures"], ['quoted "record".smi', "second.smi"])
                    self.assertEqual(generator.selected_fixtures(manifest, list(reversed(manifest["fixtures"]))), manifest["fixtures"])
                    with self.assertRaises(SystemExit):
                        generator.selected_fixtures(manifest, ["missing.smi"])
                    path.write_text('fixtures = [1]\n')
                    with self.assertRaises(SystemExit):
                        generator.read_manifest(path)
                    path.write_text('fixtures = ["unterminated]\n')
                    with self.assertRaises(tomllib.TOMLDecodeError):
                        generator.read_manifest(path)

    def test_biopython_discovers_the_current_corpus_directory(self):
        generator = load_generator("biopython")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            corpus = root / "benchmarks" / "corpora" / "pdb-100"
            manifest = corpus / "features" / "io.mmcif.parse.toml"
            manifest.parent.mkdir(parents=True)
            manifest.write_text('feature_id = "io.mmcif.parse"\ncorpus_id = "pdb-100"\nfixtures = ["data/record.cif"]\n')
            arguments = ["run_feature.py", "--feature", "io.mmcif.parse", "--corpus", "pdb-100", "--repo-root", str(root)]
            with patch("sys.argv", arguments), patch.object(generator, "import_biopython", return_value={}), patch.object(generator, "generate_fixture", return_value=Path("output")) as generate, contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(generator.main(), 0)
            self.assertEqual(generate.call_args.args[0][3], str(corpus.resolve()))


if __name__ == "__main__":
    unittest.main()
