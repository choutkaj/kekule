"""OpenFF validation commands. Run with --help for the workflow and subcommands."""
import argparse
from pathlib import Path
import runpy
import sys

COMMANDS = {
    'energy': 'robustness.py',
    'gradients': 'gradients.py',
    'timings': 'timings.py',
    'plot': 'plot_validation.py',
    'export': 'export_model.py',
    'models-reference': 'reference_models.py',
    'models-compare': 'compare_models.py',
    'offxml-reference': 'reference_offxml.py',
    'prepare': 'prepare_robustness.py',
    'audit': 'audit.py',
}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, epilog='See VALIDATION.md for complete reproducible commands. Scientific commands use environment.yml; each command accepts --help.')
    parser.add_argument('command', choices=COMMANDS)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    scripts = Path(__file__).resolve().parent / 'scripts'
    sys.path.insert(0, str(scripts))
    script = scripts / COMMANDS[args.command]
    sys.argv = [str(script), *args.arguments]
    runpy.run_path(str(script), run_name='__main__')
