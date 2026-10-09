"""Freeze an externally sourced panel before observing native parameterization.

Run in environment.yml. Selection never calls kekule. Preparation failures remain
in the manifest, alongside the first ten independently preparable protein chains.
"""
import argparse
import gzip
import hashlib
import io
import json
import random
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from openmm import Platform, unit as omm_unit
from openmm.app import ForceField as AmberForceField, Modeller, PDBFile, PDBxFile
from openff.toolkit import Topology
from openff.toolkit.utils import RDKitToolkitWrapper, ToolkitRegistry

from paths import HERE
CORPORA = HERE.parent / 'corpora'
AA = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked(root, entry):
    path = root / entry['path']
    if sha(path) != entry['sha256']:
        raise ValueError(f'corpus hash mismatch: {path}')
    return path


def ranked(value):
    return hashlib.sha256(('rosemary-robustness-v1:' + value).encode()).digest()


def small_molecules():
    root = CORPORA / 'pubchem-100k'
    lock = json.loads((root / 'sources.lock.json').read_text())
    categories = {e['id']: e['category'] for e in lock['entries']}
    candidates = []
    # Spread the independent candidate pool across the entire locked snapshot.
    source_rows = []
    for pack in lock['packs']:
        if not pack['path'].endswith('.smi'):
            continue
        for line in checked(root, pack).read_text().splitlines():
            smiles, cid = line.split()
            cid = cid.removeprefix('CID:')
            source_rows.append((ranked(cid), cid, smiles, pack))
    rejected = []
    for _, cid, smiles, pack in sorted(source_rows)[:4000]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or len(Chem.GetMolFrags(mol)) != 1 or not 3 <= mol.GetNumHeavyAtoms() <= 40:
            continue
        if any(a.GetAtomicNum() not in (1,6,7,8,9,15,16,17,35,53) or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
            continue
        elements = {a.GetSymbol() for a in mol.GetAtoms()}
        rings = mol.GetRingInfo().AtomRings()
        descriptors = {('element', z) for z in elements}
        descriptors |= {('category', categories[cid]), ('size', mol.GetNumHeavyAtoms() // 5),
                        ('charge', Chem.GetFormalCharge(mol)), ('rings', min(len(rings), 5))}
        descriptors |= {('ring-size', len(r)) for r in rings}
        descriptors |= {('atom-environment', a.GetAtomicNum(), a.GetFormalCharge(), str(a.GetHybridization()), a.GetIsAromatic()) for a in mol.GetAtoms()}
        candidates.append((cid, smiles, pack, mol, descriptors))
    selected, seen = [], set()
    for _ in range(100):
        best = max(range(len(candidates)), key=lambda i: len(candidates[i][4] - seen))
        row = candidates.pop(best)
        selected.append(row)
        seen |= row[4]
    result = []
    for cid, source, pack, mol, _ in selected:
        record = dict(id='pubchem-' + cid, kind='small-molecule', category=categories[cid],
                      source=dict(corpus='pubchem-100k', path=pack['path'], sha256=pack['sha256'], smiles=source))
        try:
            mol = Chem.AddHs(mol)
            for i, atom in enumerate(mol.GetAtoms()):
                atom.SetAtomMapNum(i + 1)
            options = AllChem.ETKDGv3()
            options.randomSeed = 20260922
            if AllChem.EmbedMolecule(mol, options) != 0:
                raise ValueError('ETKDG embedding failed')
            record.update(smiles=Chem.MolToSmiles(mol, canonical=False, allHsExplicit=True),
                          coordinates_nm=(mol.GetConformer().GetPositions() * 0.1).tolist(),
                          atoms=mol.GetNumAtoms(), formal_charge=Chem.GetFormalCharge(mol))
        except Exception as exc:
            record['preparation_error'] = f'{type(exc).__name__}: {exc}'
        result.append(record)
    return result, dict(candidate_count=4000, eligible_count=len(candidates)+100,
                        descriptors_covered=sorted(map(str, seen)), rejected=rejected)


def proteins(count=10, residue_range=(25, 180), min_atoms=0):
    """Prepare the first `count` eligible chains in the ranked corpus order.

    The defaults reproduce the frozen robustness panel exactly.
    """
    root = CORPORA / 'pdb-1000'
    lock = json.loads((root / 'sources.lock.json').read_text())
    registry = ToolkitRegistry([RDKitToolkitWrapper()])
    amber = AmberForceField('amber14-all.xml')
    result, attempts = [], []
    for entry in sorted(lock['entries'], key=lambda e: ranked(e['id'])):
        if len(result) == count:
            break
        file = next(f for f in entry['files'] if f['path'].endswith('.cif'))
        path = checked(root, file)
        try:
            pdb = PDBxFile(str(path))
        except Exception as exc:
            attempts.append(dict(id=entry['id'], stage='read', error=f'{type(exc).__name__}: {exc}'))
            continue
        for chain in pdb.topology.chains():
            residues = list(chain.residues())
            if not residue_range[0] <= len(residues) <= residue_range[1] or any(r.name not in AA for r in residues):
                continue
            attempt = dict(id=entry['id'], chain=chain.id, residues=len(residues))
            attempts.append(attempt)
            try:
                # Do not cut an inter-chain disulfide or other covalent linkage.
                if any((a.residue.chain == chain) != (b.residue.chain == chain) for a,b in pdb.topology.bonds()):
                    raise ValueError('chain has an inter-chain covalent bond')
                modeller = Modeller(pdb.topology, pdb.positions)
                modeller.delete([a for a in pdb.topology.atoms() if a.residue.chain != chain or a.element.symbol == 'H'])
                random.seed(20260922)
                variants = modeller.addHydrogens(amber, pH=7.0, platform=Platform.getPlatformByName('Reference'))
                stream = io.StringIO()
                PDBFile.writeFile(modeller.topology, modeller.positions, stream, keepIds=True)
                prepared_pdb = stream.getvalue()
                off_top = Topology.from_pdb(io.StringIO(prepared_pdb), toolkit_registry=registry)
                if off_top.n_molecules != 1:
                    raise ValueError(f'prepared chain has {off_top.n_molecules} disconnected molecules')
                off = off_top.molecule(0)
                if off.n_atoms < min_atoms:
                    raise ValueError(f'prepared chain has {off.n_atoms} atoms, below {min_atoms}')
                rd = off.to_rdkit(toolkit_registry=registry)
                for i, atom in enumerate(rd.GetAtoms()):
                    atom.SetAtomMapNum(i+1)
                coords = np.asarray(modeller.positions.value_in_unit(omm_unit.nanometer))
                assert len(coords) == off.n_atoms
                record = dict(id=f'pdb-{entry["id"]}-{chain.id}', kind='protein-chain',
                              source=dict(corpus='pdb-1000', **file), chain=chain.id,
                              residues=len(residues), sequence=[r.name for r in residues],
                              hydrogen_variants=variants, prepared_pdb=prepared_pdb,
                              smiles=Chem.MolToSmiles(rd, canonical=False, allHsExplicit=True),
                              coordinates_nm=coords.tolist(), atoms=off.n_atoms,
                              formal_charge=int(Chem.GetFormalCharge(rd)))
                result.append(record)
                attempt['selected'] = True
                print(json.dumps({k:v for k,v in record.items() if k in ('id','atoms','residues','formal_charge')}), flush=True)
                break
            except Exception as exc:
                attempt['error'] = f'{type(exc).__name__}: {exc}'
                print(json.dumps(attempt), flush=True)
    if len(result) != count:
        raise ValueError(f'could not prepare {count} protein chains')
    return result, attempts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    small, selection = small_molecules()
    protein, attempts = proteins()
    import importlib.metadata
    packages = ['rdkit','openff-toolkit','openmm','numpy']
    report = dict(schema=1, selection='rosemary-robustness-v1',
                  versions={p:importlib.metadata.version(p) for p in packages},
                  corpus_locks={name:sha(CORPORA/name/'sources.lock.json') for name in ('pubchem-100k','pdb-1000')},
                  small_selection=selection, protein_preparation_attempts=attempts,
                  records=small+protein)
    args.output.write_bytes(gzip.compress((json.dumps(report, allow_nan=False)+'\n').encode(), mtime=0))
    print(json.dumps(dict(output=str(args.output), sha256=sha(args.output), cases=len(small+protein))))


if __name__ == '__main__':
    main()
