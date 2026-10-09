"""Independent reference generation and replayable native scientific validation.

Reference writes an append-only journal; compare never changes its reference.
Vacuum energies retain constrained valence terms. Three geometries per molecule
sample embedded/crystal coordinates and seeded 0.002/0.01 nm displacements.
"""
import argparse
import collections
import gzip
import hashlib
import json
import math
import queue
import subprocess
import threading
import time
from pathlib import Path

from paths import HERE
CHARGE_ATOL = 5e-5
ENERGY_ATOL = 1e-7
ENERGY_RTOL = 1e-10
HANDLERS = ('Bonds','Angles','ProperTorsions','ImproperTorsions','vdW','Electrostatics')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes())


def split_system(omm, nproper, nimproper):
    """Separate exported force components without changing their parameters."""
    import openmm as mm
    # Split the exported OpenMM forces, preserving their exact parameters.
    split = mm.System()
    for i in range(omm.getNumParticles()):
        split.addParticle(omm.getParticleMass(i))
    for force in omm.getForces():
        if isinstance(force,mm.NonbondedForce):
            for component in ('vdW','Electrostatics'):
                clone = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(force))
                for i in range(clone.getNumParticles()):
                    q,s,e = clone.getParticleParameters(i)
                    clone.setParticleParameters(i, q if component=='Electrostatics' else 0, s, e if component=='vdW' else 0)
                for i in range(clone.getNumExceptions()):
                    a,b,q,s,e = clone.getExceptionParameters(i)
                    clone.setExceptionParameters(i,a,b,q if component=='Electrostatics' else 0,s,e if component=='vdW' else 0)
                clone.setForceGroup(HANDLERS.index(component));split.addForce(clone)
        elif isinstance(force,mm.PeriodicTorsionForce):
            assert force.getNumTorsions()==nproper+nimproper
            for component,indices in [('ProperTorsions',range(nproper)),('ImproperTorsions',range(nproper,force.getNumTorsions()))]:
                clone = mm.PeriodicTorsionForce()
                for i in indices:
                    clone.addTorsion(*force.getTorsionParameters(i))
                clone.setForceGroup(HANDLERS.index(component));split.addForce(clone)
        elif isinstance(force,(mm.HarmonicBondForce,mm.HarmonicAngleForce)):
            clone = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(force))
            clone.setForceGroup(HANDLERS.index('Bonds' if isinstance(force,mm.HarmonicBondForce) else 'Angles'))
            split.addForce(clone)
        elif not isinstance(force,mm.CMMotionRemover):
            raise ValueError(f'unexpected OpenMM force {type(force)}')
    return split


def reference(args):
    import numpy as np
    import torch
    import openmm as mm
    from openmm import unit as ou
    from openff.toolkit import Molecule, ForceField
    from openff.toolkit.utils import RDKitToolkitWrapper, ToolkitRegistry
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    from openff.toolkit.utils.toolkit_registry import toolkit_registry_manager
    from openff.nagl import GNNModel
    from openff.nagl_models import get_model
    from openff.units import unit
    from audit import MODEL_HASH
    torch.set_num_threads(2)
    model_path = get_model('openff-gnn-am1bcc-1.0.0.pt')
    assert digest(model_path) == MODEL_HASH
    model = GNNModel.load(model_path)
    registry = ToolkitRegistry([RDKitToolkitWrapper(), NAGLToolkitWrapper()])
    ff = ForceField(str(HERE/'fixtures/rosemary.offxml'))
    units = {
        'Bonds':dict(length=unit.nanometer,k=unit.kilojoule_per_mole/unit.nanometer**2),
        'Angles':dict(angle=unit.radian,k=unit.kilojoule_per_mole/unit.radian**2),
        'Constraints':dict(distance=unit.nanometer),
        'ProperTorsions':dict(k=unit.kilojoule_per_mole,phase=unit.radian,periodicity=unit.dimensionless,idivf=unit.dimensionless),
        'ImproperTorsions':dict(k=unit.kilojoule_per_mole,phase=unit.radian,periodicity=unit.dimensionless,idivf=unit.dimensionless),
        'vdW':dict(sigma=unit.nanometer,epsilon=unit.kilojoule_per_mole)}
    existing = []
    if args.output.exists():
        existing = [json.loads(line) for line in args.output.read_text().splitlines()]
        assert all(r['inputs_sha256']==digest(args.inputs) for r in existing)
    done = {r['id'] for r in existing}
    with args.output.open('a', encoding='utf-8') as output, toolkit_registry_manager(registry):
        for case in read(args.inputs)['records']:
            if case['id'] in done:
                continue
            row = dict(id=case['id'], kind=case['kind'], inputs_sha256=digest(args.inputs),
                       checkpoint_sha256=MODEL_HASH, forcefield_sha256=digest(HERE/'fixtures/rosemary.offxml'))
            started = time.monotonic()
            try:
                if 'preparation_error' in case:
                    raise ValueError(case['preparation_error'])
                off = Molecule.from_mapped_smiles(case['smiles'], allow_undefined_stereo=True)
                row['atoms'] = off.n_atoms
                row['inchi'] = off.to_inchi(fixed_hydrogens=True)
                labels = ff.label_molecules(off.to_topology())[0]
                row['labels'] = {h:[dict(atoms=[i+1 for i in key],id=p.id,smirks=p.smirks) for key,p in entries.items()] for h,entries in labels.items() if h in units}
                system = ff.create_interchange(off.to_topology(), toolkit_registry=registry)
                row['parameters'] = {h:[dict(atoms=[i+1 for i in key.atom_indices],mult=getattr(key,'mult',None),
                    values={name:float(collection.potentials[p].parameters[name].m_as(u)) for name,u in units[h].items()})
                    for key,p in collection.key_map.items()] for h,collection in system.collections.items() if h in units}
                row['charges'] = [float(v.m_as(unit.elementary_charge)) for k,v in sorted(system['Electrostatics'].charges.items(),key=lambda kv:kv[0].atom_indices)]
                row['features'] = np.hstack([f.encode(off).numpy().reshape(off.n_atoms,-1) for f in model.config.atom_features]).tolist()
                coords = np.asarray(case['coordinates_nm'])
                rng = np.random.default_rng(20260922)
                frames = [coords] + [coords+rng.normal(0,scale,coords.shape) for scale in (0.002,0.01)]
                row['coordinates_nm'] = [f.tolist() for f in frames]
                # Select a common vacuum Hamiltonian explicitly, without changing assignments.
                system['vdW'].nonperiodic_method = 'no-cutoff'
                omm = system.to_openmm_system(combine_nonbonded_forces=True, add_constrained_forces=True)
                reference_system_xml = mm.XmlSerializer.serialize(omm)
                nb = next(f for f in omm.getForces() if isinstance(f,mm.NonbondedForce))
                row['exceptions'] = [dict(atoms=[int(i)+1,int(j)+1],charge_product=float(q.value_in_unit(ou.elementary_charge**2)),
                    sigma=float(s.value_in_unit(ou.nanometer)),epsilon=float(e.value_in_unit(ou.kilojoule_per_mole)))
                    for i,j,q,s,e in (nb.getExceptionParameters(k) for k in range(nb.getNumExceptions()))]
                split = split_system(omm, len(system['ProperTorsions'].key_map), len(system['ImproperTorsions'].key_map))
                integrator = mm.VerletIntegrator(0.001)
                context = mm.Context(split,integrator,mm.Platform.getPlatformByName('Reference'))
                row['energies'] = []
                for frame in frames:
                    context.setPositions(frame*ou.nanometer)
                    values = {h:float(context.getState(getEnergy=True,groups=1<<i).getPotentialEnergy().value_in_unit(ou.kilojoule_per_mole)) for i,h in enumerate(HANDLERS)}
                    values['Total'] = float(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(ou.kilojoule_per_mole))
                    row['energies'].append(values)
                del context, integrator
                integrator = mm.VerletIntegrator(0.001)
                context = mm.Context(omm,integrator,mm.Platform.getPlatformByName('Reference'))
                context.setPositions(frames[0]*ou.nanometer)
                unsplit = float(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(ou.kilojoule_per_mole))
                # OpenMM self-consistency of the split export. Its summation order
                # differs from the unsplit system, so rounding grows with the number
                # of pair terms; a check stricter than the comparison is meaningless.
                assert math.isclose(unsplit,row['energies'][0]['Total'],rel_tol=ENERGY_RTOL,abs_tol=1e-8)
                row['unsplit_energy'] = unsplit
                row['openmm_system_xml'] = reference_system_xml
                del context,integrator
                row['status'] = 'ok'
            except Exception as exc:
                row.update(status='error',error=f'{type(exc).__name__}: {exc}')
            row['seconds'] = time.monotonic()-started
            output.write(json.dumps(row,allow_nan=False)+'\n');output.flush()
            print(json.dumps({k:v for k,v in row.items() if k in ('id','atoms','status','error','seconds')}),flush=True)


def canonical(atoms,handler):
    if handler=='ImproperTorsions':
        return (atoms[0],*sorted(atoms[1:]))
    return tuple(min(atoms,atoms[::-1]))


def compare_parameters(expected,native):
    errors, summary = [], {}
    fields = dict(Bonds='bonds',Angles='angles',Constraints='constraints',ProperTorsions='propers',ImproperTorsions='impropers',vdW='vdw')
    for handler,field in fields.items():
        wanted,actual = collections.defaultdict(list),collections.defaultdict(list)
        for row in expected[handler]:
            wanted[(canonical(row['atoms'],handler),row['mult'])].append(row['values'])
        for i,row in enumerate(native['system'][field]):
            atoms = [native['maps'][i]] if field=='vdw' else row['maps']
            p = row if field=='vdw' else row['parameter']
            for mult,values in enumerate(p['terms']) if 'terms' in p else [(None,p)]:
                actual[(canonical(atoms,handler),mult)].append({k:v for k,v in values.items() if k!='id'})
        if wanted.keys()!=actual.keys():
            errors.append(dict(handler=handler,missing=list(map(str,wanted.keys()-actual.keys())),extra=list(map(str,actual.keys()-wanted.keys()))))
        maximum=0.0
        for key in wanted.keys() & actual.keys():
            a,b = wanted[key],actual[key]
            if len(a)!=len(b):
                errors.append(dict(handler=handler,key=str(key),multiplicity=[len(a),len(b)]))
            for x,y in zip(sorted(a,key=lambda v:tuple(sorted(v.items()))),sorted(b,key=lambda v:tuple(sorted(v.items())))):
                if x.keys()!=y.keys():
                    errors.append(dict(handler=handler,key=str(key),fields=[list(x),list(y)]));continue
                for name in x:
                    maximum=max(maximum,abs(x[name]-y[name]))
                    if not math.isclose(x[name],y[name],rel_tol=1e-12,abs_tol=1e-10):
                        errors.append(dict(handler=handler,key=str(key),parameter=name,reference=x[name],native=y[name]))
        summary[handler]=dict(expected_count=sum(map(len,wanted.values())),actual_count=sum(map(len,actual.values())),max_abs_error=maximum)
    return errors,summary


class Observer:
    def __init__(self,binary,model):
        self.process=subprocess.Popen([str(binary.resolve()),str(model.resolve())],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
        self.lines=queue.Queue()
        threading.Thread(target=self._read,daemon=True).start()
    def _read(self):
        for line in self.process.stdout:
            self.lines.put(line)
        self.lines.put('')
    def call(self,request,timeout):
        self.process.stdin.write(json.dumps(request)+'\n');self.process.stdin.flush()
        try:
            line=self.lines.get(timeout=timeout)
        except queue.Empty:
            self.process.kill();self.process.wait()
            raise TimeoutError(f'native observer exceeded {timeout} seconds')
        if not line:
            raise RuntimeError(f'native observer exited: {self.process.poll()}')
        return json.loads(line)
    def close(self):
        if self.process.poll() is None:
            self.process.stdin.close();self.process.wait(timeout=10)


def compare(args):
    import numpy as np
    from rdkit import Chem
    inputs=read(args.inputs)
    reference=read(args.reference)
    assert reference['inputs_sha256']==digest(args.inputs)
    tables=json.loads((args.model/'model.json').read_text(encoding='utf-8'))['lookup_tables']
    lookup_keys={e['inchi'] for entries in tables.values() for e in entries}
    rows=[]
    journal=args.output.with_suffix('.jsonl').open('x',encoding='utf-8')
    observer=Observer(args.binary,args.model)
    try:
        for case,ref in zip(inputs['records'],reference['records'],strict=True):
            assert case['id']==ref['id']
            for permutation in ('original','reversed'):
                row=dict(id=case['id'],kind=case['kind'],permutation=permutation,errors=[])
                start=time.monotonic()
                try:
                    if ref['status']!='ok':
                        raise ValueError('reference: '+ref['error'])
                    # Default RDKit parsing removes explicit H, so retain them explicitly.
                    params=Chem.SmilesParserParams();params.removeHs=False
                    rd=Chem.MolFromSmiles(case['smiles'],params)
                    if permutation=='reversed':
                        rd=Chem.RenumberAtoms(rd,list(reversed(range(rd.GetNumAtoms()))))
                    source=Chem.MolToSmiles(rd,canonical=False,allHsExplicit=True)
                    native=observer.call(dict(smiles=source,mode='validate',coordinates_nm=ref['coordinates_nm'],reference_charges=ref['charges']),args.timeout)
                    row['native']=native
                    if native['status']!='ok':
                        raise ValueError('native: '+native['message'])
                    order=[native['maps'].index(i+1) for i in range(len(ref['charges']))]
                    errors,summary=compare_parameters(ref['parameters'],native)
                    row['errors'].extend(errors);row['parameters']=summary
                    # Upstream selects the entry keyed by its fixed-H InChI, if any.
                    expected=ref['inchi'] if ref['inchi'] in lookup_keys else None
                    if 'identity_error' in native or native['lookup_key']!=expected:
                        row['errors'].append(dict(identity=native.get('identity_error','lookup selection mismatch')))
                    q=np.asarray(native['system']['charges'])[order]
                    delta=np.abs(q-ref['charges'])
                    row['charge_max_error_e']=float(delta.max())
                    row['charge_sum_error_e']=float(abs(q.sum()-case['formal_charge']))
                    if delta.max()>CHARGE_ATOL or row['charge_sum_error_e']>1e-6:
                        row['errors'].append(dict(charges=row['charge_max_error_e'],sum=row['charge_sum_error_e']))
                    if native['features']['status']!='ok':
                        raise ValueError('features: '+native['features']['message'])
                    delta=np.abs(np.asarray(native['features']['values'])[order]-ref['features'])
                    row['feature_max_error']=float(delta.max())
                    if delta.max()>1e-6:
                        row['errors'].append(dict(features=float(delta.max()),mismatches=np.argwhere(delta>1e-6).tolist()))
                    exceptions={tuple(sorted(p['maps'])):p for p in native['system']['exceptions']}
                    expected_exceptions={tuple(sorted(p['atoms'])):p for p in ref['exceptions']}
                    if exceptions.keys()!=expected_exceptions.keys():
                        row['errors'].append('pair exception set mismatch')
                    vdw=[native['system']['vdw'][i] for i in order]
                    for pair in exceptions.keys() & expected_exceptions.keys():
                        i,j=(k-1 for k in pair)
                        a,b=exceptions[pair],expected_exceptions[pair]
                        qp=ref['charges'][i]*ref['charges'][j]*a['electrostatics_scale']
                        ep=math.sqrt(vdw[i]['epsilon']*vdw[j]['epsilon'])*a['vdw_scale']
                        if not math.isclose(qp,b['charge_product'],rel_tol=1e-12,abs_tol=1e-12) or not math.isclose(ep,b['epsilon'],rel_tol=1e-12,abs_tol=1e-12):
                            row['errors'].append(dict(pair=pair,reference=b,native_charge_product=qp,native_epsilon=ep))
                    # Exact IDs/SMIRKS, allowing only conventional key orientation.
                    for handler,labels in ref['labels'].items():
                        def labelkey(atoms):
                            if handler=='ImproperTorsions':
                                return (atoms[1],*sorted([atoms[0],atoms[2],atoms[3]]))
                            return canonical(atoms,handler)
                        a={(labelkey(r['atoms']),r['id'],r['smirks']) for r in labels}
                        b={(labelkey(r['maps']),r['id'],r['smirks']) for r in native['labels'][handler]}
                        if a!=b:
                            row['errors'].append(dict(label_handler=handler,missing=list(map(str,a-b)),extra=list(map(str,b-a))))
                    row['energy_errors']=[]
                    for frame,(expected,observed) in enumerate(zip(ref['energies'],native['energies'],strict=True)):
                        residual={h:observed['reference_charges'][h]-v for h,v in expected.items()}
                        end_to_end={h:observed['native_charges'][h]-v for h,v in expected.items()}
                        # A charge-tolerance-derived Coulomb bound, not a relaxed energy tolerance.
                        xyz=np.asarray(ref['coordinates_nm']);bound=0.0
                        scales={tuple(sorted(p['maps'])):p['electrostatics_scale'] for p in native['system']['exceptions']}
                        qr=np.asarray(ref['charges']);dq=q-qr
                        row_scales=collections.defaultdict(list)
                        for (a,b),scale in scales.items():
                            row_scales[a-1].append((b-1,scale))
                        for i in range(len(q)):
                            weights=np.ones(len(q)-i-1)
                            for j,scale in row_scales[i]:
                                weights[j-i-1]=scale
                            distance=np.linalg.norm(xyz[frame,i]-xyz[frame,i+1:],axis=1)
                            bound+=float(np.sum(138.93545764438198*weights*(abs(qr[i]*dq[i+1:])+abs(qr[i+1:]*dq[i])+abs(dq[i]*dq[i+1:]))/distance))
                        row['energy_errors'].append(dict(frame=frame,reference_charges=residual,native_charges=end_to_end,charge_propagation_bound_kj_mol=bound))
                        for h,v in expected.items():
                            if not math.isclose(observed['reference_charges'][h],v,rel_tol=ENERGY_RTOL,abs_tol=ENERGY_ATOL):
                                row['errors'].append(dict(energy=h,frame=frame,reference=v,residual=residual[h]))
                        if abs(end_to_end['Electrostatics'])>bound+ENERGY_ATOL+ENERGY_RTOL*abs(expected['Electrostatics']):
                            row['errors'].append(dict(electrostatic_end_to_end=frame,error=end_to_end['Electrostatics'],bound=bound))
                except Exception as exc:
                    row['errors'].append(f'{type(exc).__name__}: {exc}')
                    if observer.process.poll() is not None:
                        observer=Observer(args.binary,args.model)
                row.update(passed=not row['errors'],
                           parameterization_passed=not any(not isinstance(e,dict) or 'identity' not in e for e in row['errors']),
                           seconds=time.monotonic()-start)
                rows.append(row)
                print(json.dumps({k:v for k,v in row.items() if k not in ('native','parameters','energy_errors')}),flush=True)
                # Append-only progress also avoids racing readers of the final archive.
                journal.write(json.dumps(row,allow_nan=False)+'\n');journal.flush()
    finally:
        observer.close()
        journal.close()
    report=dict(schema=1,inputs_sha256=digest(args.inputs),reference_sha256=digest(args.reference),binary_sha256=digest(args.binary),
        charge_atol_e=CHARGE_ATOL,energy_atol_kj_mol=ENERGY_ATOL,energy_rtol=ENERGY_RTOL,records=rows)
    args.output.write_bytes(gzip.compress((json.dumps(report,allow_nan=False)+'\n').encode(),mtime=0))
    return int(any(not row['passed'] for row in rows))


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



def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['reference','compare','freeze'])
    parser.add_argument('--inputs',type=Path,default=HERE/'data/inputs.json.gz')
    parser.add_argument('--reference',type=Path,default=HERE/'data/reference.json.gz')
    parser.add_argument('--binary',type=Path,default=HERE.parent.parent/'target/release/openff_parameterize.exe')
    parser.add_argument('--model',type=Path,default=HERE.parent.parent/'target/openff-models/openff-gnn-am1bcc-1.0.0')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=float,default=180)
    args=parser.parse_args()
    if args.mode=='reference':
        reference(args)
        return 0
    if args.output.exists():
        raise FileExistsError(args.output)
    if args.mode=='freeze':
        records=[json.loads(line) for line in args.reference.read_text().splitlines()]
        assert [r['id'] for r in records]==[r['id'] for r in read(args.inputs)['records']]
        args.output.write_bytes(gzip.compress((json.dumps(dict(schema=1,inputs_sha256=digest(args.inputs),records=records),allow_nan=False)+'\n').encode(),mtime=0))
        return 0
    return compare(args)


if __name__=='__main__':
    # Toolkit conversion of large protein chains exceeds the default Windows
    # thread stack, so run on a thread with a larger one.
    import sys
    sys.setrecursionlimit(100_000)
    threading.stack_size(256*1024*1024-4096)
    status=[1]
    worker=threading.Thread(target=lambda: status.__setitem__(0, main()))
    worker.start();worker.join()
    raise SystemExit(status[0])
