"""Classify exhaustive lookup disagreements using the pinned official toolkit.

This diagnostic preserves every native failure and the stored key, even when
the current reference also disagrees. It does not change lookup entries.
"""
import argparse
import json
from pathlib import Path
from collections import Counter
from rdkit import RDLogger
from openff.toolkit import Molecule
from openff.toolkit.utils import RDKitToolkitWrapper, ToolkitRegistry
from openff.toolkit.utils.toolkit_registry import toolkit_registry_manager
from openff.nagl import GNNModel
from openff.nagl_models import get_model
from audit import MODEL_HASH, digest


def run(source, output):
    report = json.loads(source.read_text())
    checkpoint = get_model('openff-gnn-am1bcc-1.0.0.pt')
    assert digest(checkpoint) == MODEL_HASH
    model = GNNModel.load(checkpoint)
    counts = Counter()
    RDLogger.DisableLog('rdApp.*')
    with toolkit_registry_manager(ToolkitRegistry([RDKitToolkitWrapper()])):
        for row in report['disagreements']:
            try:
                mol = Molecule.from_mapped_smiles(row['smiles'], allow_undefined_stereo=True)
                row['reference_domain'] = list(model.chemical_domain.check_molecule(mol, return_error_message=True))
                row['reference_identifier'] = mol.to_inchi(fixed_hydrogens=True)
                if row['actual']['status'] != 'ok':
                    category = 'native-rejected-reference-accepted'
                elif row['reference_identifier'] == row['actual']['fixed_h_inchi']:
                    category = 'reference-and-native-agree-stored-key-differs'
                else:
                    category = 'native-reference-identifier-difference'
            except Exception as e:
                row['reference_error'] = dict(type=type(e).__name__, message=str(e))
                category = 'reference-rejected'
            row['category'] = category
            counts[category] += 1
    report['classification'] = dict(counts)
    report['source_report_sha256'] = digest(source)
    output.write_text(json.dumps(report, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    run(args.source, args.output)
