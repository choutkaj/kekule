"""Independent, immutable two-model OpenFF observations; never calls the Rust engine.

Uses the existing externally sourced 33-molecule panel in both atom orders.
Run in the pinned reference environment. Outputs and bundles must not already exist.
"""
import argparse
import gzip
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from openff.nagl import GNNModel
from openff.nagl_models import get_model
from openff.toolkit import ForceField, Molecule
from openff.toolkit.utils import RDKitToolkitWrapper, ToolkitRegistry
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from openff.toolkit.utils.toolkit_registry import toolkit_registry_manager
from openff.units import unit

from paths import ROOT
from audit import HERE, MODEL_HASH, VERSIONS, describe_configuration, digest, mapped_rdkit
from export_model import export
from cases import sources

MODELS = {
    "openff-gnn-am1bcc-1.0.0.pt": MODEL_HASH,
    "openff-gnn-am1bcc-0.1.0-rc.2.pt": "648b2636580f49f882591aedcc5c404a9cbbecb9ca1e082d98ccd71301db917f",
}


def values(collection):
    """Retain all numeric parameters, multiplicities, and atom identities."""
    units = {
        'length': unit.nanometer, 'distance': unit.nanometer,
        'angle': unit.radian, 'phase': unit.radian,
        'sigma': unit.nanometer, 'epsilon': unit.kilojoule_per_mole,
        'periodicity': unit.dimensionless, 'idivf': unit.dimensionless,
    }
    rows = []
    for key, potential in collection.key_map.items():
        parameters = collection.potentials[potential].parameters
        row = {}
        for name, value in parameters.items():
            target = units.get(name)
            if name == 'k':
                target = unit.kilojoule_per_mole
                if collection.type == 'Bonds': target /= unit.nanometer**2
                elif collection.type == 'Angles': target /= unit.radian**2
            if target is None: raise ValueError(f'Unexpected parameter {name}')
            row[name] = float(value.m_as(target))
        rows.append(dict(maps=[i + 1 for i in key.atom_indices], mult=getattr(key, 'mult', None), parameters=row))
    return rows


def generate(bundle_root, output):
    if output.exists(): raise FileExistsError(output)
    versions = {name: importlib.metadata.version(name) for name in VERSIONS}
    if versions != VERSIONS: raise ValueError(f'Reference environment mismatch: {versions}')
    rosemary = (HERE / 'fixtures/rosemary.offxml').read_text(encoding='utf-8')
    models = []
    for name, sha in MODELS.items():
        checkpoint = Path(get_model(name))
        if digest(checkpoint) != sha: raise ValueError(f'Checkpoint mismatch: {name}')
        directory = bundle_root / name.removesuffix('.pt')
        export(directory, checkpoint, sha, HERE / 'fixtures/LICENSE-models')
        model = GNNModel.load(str(checkpoint))
        xml = rosemary.replace('openff-gnn-am1bcc-1.0.0.pt', name).replace(MODEL_HASH, sha)
        (directory / 'force-field.offxml').write_text(xml, encoding='utf-8', newline='\n')
        ff = ForceField(xml)
        registry = ToolkitRegistry([RDKitToolkitWrapper(), NAGLToolkitWrapper()])
        records = []
        with toolkit_registry_manager(registry):
            for case in sources():
                for reverse in (False, True):
                    rd = mapped_rdkit(case['smiles'])
                    if reverse: rd = Chem.RenumberAtoms(rd, list(range(rd.GetNumAtoms()))[::-1])
                    smiles = Chem.MolToSmiles(rd, canonical=False, allHsExplicit=True)
                    off = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
                    record = dict(id=case['id'], reverse=reverse, smiles=smiles)
                    try:
                        if not model.chemical_domain.check_molecule(off): raise ValueError('outside model domain')
                        features = np.hstack([f.encode(off).numpy().reshape(off.n_atoms, -1) for f in model.config.atom_features])
                        raw = model._compute_properties_nagl(off)['am1bcc_charges'].detach().numpy().flatten().astype(float)
                        raw += (off.total_charge.m_as(unit.elementary_charge) - raw.sum()) / len(raw)
                        record['inference'] = dict(status='ok', features=features.tolist(), charges=raw.tolist())
                    except Exception as e:
                        record['inference'] = dict(status='error', error_type=type(e).__name__, message=str(e))
                    try:
                        charges = model.compute_property(off, readout_name='am1bcc_charges', check_domains=True).flatten().astype(float)
                        charges += (off.total_charge.m_as(unit.elementary_charge) - charges.sum()) / len(charges)
                        record['assignment'] = dict(status='ok', charges=charges.tolist())
                    except Exception as e:
                        record['assignment'] = dict(status='error', error_type=type(e).__name__, message=str(e))
                    system = ff.create_interchange(off.to_topology(), toolkit_registry=registry)
                    record['system'] = dict(
                        charges=[float(v.m_as(unit.elementary_charge)) for k,v in sorted(system.collections['Electrostatics'].charges.items(), key=lambda item:item[0].atom_indices)],
                        parameters={handler: values(system.collections[handler]) for handler in ('Bonds', 'Angles', 'Constraints', 'ProperTorsions', 'ImproperTorsions', 'vdW')})
                    records.append(record)
        models.append(dict(name=name, checkpoint_sha256=sha,
                           source='https://github.com/openforcefield/openff-nagl-models',
                           package_path=f'openff/nagl_models/models/am1bcc/{name}',
                           bundle_sha256={p:digest(directory / p) for p in ('model.json', 'weights.bin', 'force-field.offxml')},
                           config=model.config.model_dump(), lookup_entries=sum(len(t.properties) for t in model.lookup_tables.values()), records=records))
        print(name, len(records), 'reference cases complete', flush=True)
    payload = dict(schema=1, versions=versions, sources_sha256={p:digest(ROOT / p) for p in ('crates/kekule-openff/tests/fixtures/audit.json.gz', 'benchmarks/openff/fixtures/supplementary-sources.lock.json', 'benchmarks/openff/fixtures/rosemary.offxml')}, models=models)
    output.write_bytes(gzip.compress((json.dumps(payload, default=describe_configuration, allow_nan=False) + '\n').encode(), mtime=0))
    print('reference sha256', digest(output), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundles', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    generate(args.bundles, args.output)
