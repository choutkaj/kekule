import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch


SCRIPT = Path(__file__).with_name("generate.py").resolve()
REPOSITORY = next(parent for parent in SCRIPT.parents if (parent / "Cargo.toml").is_file())
spec = importlib.util.spec_from_file_location("atomic_data_generator", SCRIPT)
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)


class GeneratorOutputTests(unittest.TestCase):
    def test_default_output_is_the_repository_descriptor_table(self):
        self.assert_output([], REPOSITORY / "crates/kekule/src/descriptors/data.rs")

    def test_explicit_output_is_respected(self):
        output = REPOSITORY / "target/atomic-data-test/output.rs"
        self.assert_output(["--output", str(output)], output)

    def assert_output(self, arguments, expected):
        sources = {name: Path(name) for name in generator.SOURCES}
        with (
            patch("sys.argv", [str(SCRIPT), "--offline", *arguments]),
            patch.object(generator, "acquire_sources", return_value=sources) as acquire,
            patch.object(generator, "standard_weights", return_value=[]),
            patch.object(generator, "natural_isotopes", return_value=[]),
            patch.object(generator, "isotope_masses", return_value=[]),
            patch.object(generator, "render_data", return_value="generated\n"),
            patch.object(Path, "mkdir", autospec=True),
            patch.object(Path, "write_text", autospec=True) as write,
            patch("builtins.print"),
        ):
            self.assertEqual(generator.main(), 0)
        self.assertTrue(acquire.call_args.args[1])
        write.assert_called_once_with(expected, "generated\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    unittest.main()
