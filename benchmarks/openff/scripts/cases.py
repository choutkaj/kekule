"""The 33 independently sourced molecules shared by the model validation runs."""
import gzip
import json
from rdkit import Chem
from audit import mapped_rdkit, digest
from paths import HERE, ROOT


def sources():
    rows = json.loads(gzip.decompress((ROOT / 'crates/kekule-openff/tests/fixtures/audit.json.gz').read_bytes()))['records']
    result = [dict(id=r['input']['id'], smiles=r['mapped_smiles']) for r in rows]
    lock = json.loads((HERE / 'fixtures/supplementary-sources.lock.json').read_bytes())
    for source in lock['sources']:
        path = HERE / source['path']
        if digest(path) != source['sha256']: raise ValueError(f'Source checksum mismatch: {path}')
        prop = json.loads(path.read_bytes())['PropertyTable']['Properties'][0]
        rd = mapped_rdkit(prop['SMILES'])
        result.append(dict(id='pubchem-' + str(prop['CID']), smiles=Chem.MolToSmiles(rd, canonical=False, allHsExplicit=True)))
    return result
