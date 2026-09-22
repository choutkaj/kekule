"""Summarize immutable robustness artifacts without running reference tools."""
import argparse
import collections
import json
from pathlib import Path

from robustness import read, digest

HERE=Path(__file__).resolve().parent


def summarize(inputs,reference,native):
    assert reference['inputs_sha256']==native['inputs_sha256']
    rows=native['records']
    return dict(
        molecules=len(inputs['records']), cases=len(rows),
        parameterization_passed=sum(r.get('parameterization_passed',False) for r in rows),
        all_checks_passed=sum(r['passed'] for r in rows),
        energy_evaluations=sum(len(r.get('energy_errors',[])) for r in rows),
        paired_atoms=sum(len(r.get('native',{}).get('maps',[])) for r in rows),
        charge_max_error_e=max(r.get('charge_max_error_e',0) for r in rows),
        charge_sum_max_error_e=max(r.get('charge_sum_error_e',0) for r in rows),
        feature_max_error=max(r.get('feature_max_error',0) for r in rows),
        energy_shared_charge_max_error_kj_mol={h:max(abs(e['reference_charges'][h]) for r in rows for e in r.get('energy_errors',[])) for h in ('Bonds','Angles','ProperTorsions','ImproperTorsions','vdW','Electrostatics','Total')},
        energy_native_charge_max_error_kj_mol={h:max(abs(e['native_charges'][h]) for r in rows for e in r.get('energy_errors',[])) for h in ('Electrostatics','Total')},
        parameter_rules={h:len({p['id'] for r in reference['records'] for p in r.get('labels',{}).get(h,[])}) for h in ('Bonds','Angles','ProperTorsions','ImproperTorsions','vdW','Constraints')},
        charge_sources=dict(collections.Counter(s.split('{')[0].strip() for r in rows for s in r.get('native',{}).get('charge_sources',[]))),
        remaining_failures=[{k:r[k] for k in ('id','permutation','errors')} for r in rows if not r['passed']],
        proteins=[{k:r[k] for k in ('id','residues','atoms','formal_charge')} for r in inputs['records'] if r['kind']=='protein-chain'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    paths={kind:HERE/f'robustness-{kind}.json.gz' for kind in ('inputs','reference','native')}
    report=summarize(*(read(paths[k]) for k in ('inputs','reference','native')))
    report['sha256']={k:digest(p) for k,p in paths.items()}
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
