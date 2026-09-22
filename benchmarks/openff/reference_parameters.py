"""Optional live OpenFF comparison of numbers, forced inference, and atom permutations."""
import gzip,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from openff.toolkit import Molecule, ForceField
from openff.toolkit.utils import RDKitToolkitWrapper, ToolkitRegistry
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from openff.toolkit.utils.toolkit_registry import toolkit_registry_manager
from openff.nagl import GNNModel
from openff.nagl_models import get_model
from openff.units import unit
from audit import HERE,MODEL_HASH,CHARGE_ATOL,mapped_rdkit,describe_configuration,digest
from compare_native import canonical

def sources():
    rows=json.loads(gzip.decompress((HERE/'reference.json.gz').read_bytes()))['records']
    result=[dict(id=r['input']['id'],smiles=r['mapped_smiles']) for r in rows]
    for s in json.loads((HERE/'supplementary-sources.lock.json').read_text())['sources']:
        path=HERE/s['path'];assert digest(path)==s['sha256']
        prop=json.loads(path.read_bytes())['PropertyTable']['Properties'][0]
        rd=mapped_rdkit(prop['SMILES'])
        result.append(dict(id='pubchem-'+str(prop['CID']),smiles=Chem.MolToSmiles(rd,canonical=False,allHsExplicit=True)))
    return result


def compare_numbers(off,system,native):
    errors=[]
    native_maps=native['maps']
    # Molecule.from_mapped_smiles orders by the consecutive map labels supplied here.
    atom_maps=sorted(native_maps)
    specification={
        'Bonds':('bonds',{'length':unit.nanometer,'k':unit.kilojoule_per_mole/unit.nanometer**2}),
        'Angles':('angles',{'angle':unit.radian,'k':unit.kilojoule_per_mole/unit.radian**2}),
        'Constraints':('constraints',{'distance':unit.nanometer}),
        'ProperTorsions':('propers',{'k':unit.kilojoule_per_mole,'phase':unit.radian,'periodicity':unit.dimensionless,'idivf':unit.dimensionless}),
        'ImproperTorsions':('impropers',{'k':unit.kilojoule_per_mole,'phase':unit.radian,'periodicity':unit.dimensionless,'idivf':unit.dimensionless}),
        'vdW':('vdw',{'sigma':unit.nanometer,'epsilon':unit.kilojoule_per_mole})}
    observations={}
    for handler,(field,units) in specification.items():
        expected={};observed={};expected_terms=[];observed_terms=[]
        collection=system.collections[handler]
        for key,potential in collection.key_map.items():
            atoms=[atom_maps[i] for i in key.atom_indices]
            # Interchange emits the center FIRST. Local atom ordering can select
            # either orientation of the three-member trefoil. Retain every term
            # and independently check its energy on identical mapped coordinates.
            if handler=='ImproperTorsions':
                atomkey=(atoms[0],tuple(sorted(atoms[1:])))
            else:atomkey=canonical(atoms,handler)
            mult=getattr(key,'mult',None)
            params=collection.potentials[potential].parameters
            values={name:float(params[name].m_as(u)) for name,u in units.items()}
            expected.setdefault((atomkey,mult),[]).append(values)
            if handler=='ImproperTorsions':expected_terms.append((atoms,values))
        for i,row in enumerate(native['system'][field]):
            if field=='vdw':atoms=[native_maps[i]];parameter=row
            else:atoms=row['maps'];parameter=row['parameter']
            if handler=='ImproperTorsions':
                atomkey=(atoms[0],tuple(sorted(atoms[1:])))
            else:atomkey=canonical(atoms,handler)
            if field in ('propers','impropers'):
                for mult,term in enumerate(parameter['terms']):
                    values={name:term[name] for name in units}
                    observed.setdefault((atomkey,mult),[]).append(values)
                    if handler=='ImproperTorsions':observed_terms.append((atoms,values))
            else:observed.setdefault((atomkey,None),[]).append({name:parameter[name] for name in units})
        if set(expected)!=set(observed):
            errors.append(dict(handler=handler,missing=[str(k) for k in set(expected)-set(observed)],extra=[str(k) for k in set(observed)-set(expected)]))
        numeric_error=0
        for key in set(expected)&set(observed):
            if len(expected[key])!=len(observed[key]):
                errors.append(dict(handler=handler,key=str(key),expected_multiplicity=len(expected[key]),actual_multiplicity=len(observed[key])))
            for arow,brow in zip(sorted(expected[key],key=lambda p:tuple(p.values())),sorted(observed[key],key=lambda p:tuple(p.values()))):
                for name in units:
                    a,b=arow[name],brow[name]
                    numeric_error=max(numeric_error,abs(a-b))
                    if not np.isclose(a,b,rtol=1e-12,atol=1e-10):errors.append(dict(handler=handler,key=str(key),parameter=name,expected=a,actual=b))
        observations[handler]=dict(max_abs_error=numeric_error,expected_count=sum(map(len,expected.values())),actual_count=sum(map(len,observed.values())),reference=[dict(key=str(k),values=v) for k,v in expected.items()])
        if handler=='ImproperTorsions':
            energies=[]
            for seed in range(3):
                xyz=dict(zip(atom_maps,np.random.default_rng(seed).normal(size=(len(atom_maps),3))))
                a,b=(torsion_energy(terms,xyz) for terms in (expected_terms,observed_terms))
                energies.append(dict(seed=seed,reference_kj_mol=a,native_kj_mol=b))
                if not np.isclose(a,b,rtol=1e-12,atol=1e-10):errors.append(dict(improper_energy_seed=seed,reference=a,native=b))
            observations[handler]['energies']=energies
            observations[handler]['reference_terms']=expected_terms
    return errors,observations


def torsion_energy(terms,xyz):
    total=0.0
    for atoms,p in terms:
        a,b,c,d=(xyz[i] for i in atoms)
        axis=c-b;axis/=np.linalg.norm(axis)
        v=a-b;v-=np.dot(v,axis)*axis
        w=d-c;w-=np.dot(w,axis)*axis
        theta=np.arctan2(np.dot(np.cross(axis,v),w),np.dot(v,w))
        total+=p['k']/p['idivf']*(1+np.cos(p['periodicity']*theta-p['phase']))
    return float(total)


def run(binary,model_directory,output):
    model_path=get_model('openff-gnn-am1bcc-1.0.0.pt');assert digest(model_path)==MODEL_HASH
    model=GNNModel.load(model_path)
    ff=ForceField(str(HERE/'fixtures/rosemary.offxml'))
    registry=ToolkitRegistry([RDKitToolkitWrapper()])
    system_registry=ToolkitRegistry([RDKitToolkitWrapper(),NAGLToolkitWrapper()])
    process=subprocess.Popen([str(binary.resolve()),str(model_directory.resolve())],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
    records=[]
    with toolkit_registry_manager(registry):
        for case in sources():
            for reverse in (False,True):
                rd=mapped_rdkit(case['smiles'])
                if reverse:rd=Chem.RenumberAtoms(rd,list(range(rd.GetNumAtoms()))[::-1])
                smiles=Chem.MolToSmiles(rd,canonical=False,allHsExplicit=True)
                off=Molecule.from_mapped_smiles(smiles,allow_undefined_stereo=True)
                record=dict(id=case['id'],reverse=reverse,smiles=smiles,errors=[])
                try:
                    process.stdin.write(json.dumps(dict(smiles=smiles))+'\n');process.stdin.flush()
                    native=json.loads(process.stdout.readline());record['native']=native
                    if native['status']!='ok':raise ValueError(native['message'])
                    system=ff.create_interchange(off.to_topology(),toolkit_registry=system_registry)
                    errors,numbers=compare_numbers(off,system,native);record['errors'].extend(errors);record['parameters']=numbers
                    if off.n_atoms!=len(native['maps']):raise ValueError('atom count mismatch')
                    order=[native['maps'].index(i) for i in sorted(native['maps'])]
                    if off.to_inchi(fixed_hydrogens=True)!=native['fixed_h_inchi']:record['errors'].append('InChI mismatch')
                    if model.chemical_domain.check_molecule(off):
                        raw=model._compute_properties_nagl(off)['am1bcc_charges'].detach().numpy().flatten().astype(float)
                        raw+=(off.total_charge.m_as(unit.elementary_charge)-raw.sum())/len(raw)
                        record['inference_reference']=raw.tolist()
                        if native['inference']['status']!='ok':raise ValueError(native['inference']['message'])
                        delta=max(abs(raw[j]-native['inference']['values'][i]) for j,i in enumerate(order))
                        record['inference_error_e']=float(delta)
                        if delta>CHARGE_ATOL:record['errors'].append(dict(inference_error_e=float(delta)))
                        features=np.hstack([f.encode(off).numpy().reshape(off.n_atoms,-1) for f in model.config.atom_features])
                        record['features_reference']=features.tolist()
                        if native['features']['status']!='ok':raise ValueError(native['features']['message'])
                        observed=np.array(native['features']['values'])[order]
                        delta=float(np.max(np.abs(features-observed)));record['feature_error']=delta
                        if delta>1e-6:record['errors'].append(dict(feature_error=delta,mismatch=np.argwhere(np.abs(features-observed)>1e-6).tolist()))
                    expected=[v.m_as(unit.elementary_charge) for k,v in sorted(system.collections['Electrostatics'].charges.items(),key=lambda x:x[0].atom_indices)]
                    delta=max(abs(q-native['system']['charges'][i]) for q,i in zip(expected,order))
                    record['system_charge_error_e']=float(delta)
                    if delta>CHARGE_ATOL:record['errors'].append(dict(system_charge_error_e=float(delta)))
                except Exception as e:record['errors'].append(dict(exception=type(e).__name__,message=str(e)))
                record['passed']=not record['errors'];records.append(record)
                print(json.dumps({k:v for k,v in record.items() if k not in ('native','features_reference','inference_reference','parameters','smiles')}),flush=True)
    process.stdin.close();assert process.wait()==0
    report=dict(schema=1,records=records,checkpoint_sha256=MODEL_HASH,binary_sha256=digest(binary),charge_atol_e=CHARGE_ATOL)
    output.write_text(json.dumps(report,default=describe_configuration,allow_nan=False)+'\n',encoding='utf-8')
    return int(any(not r['passed'] for r in records))

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--model',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();raise SystemExit(run(a.binary,a.model,a.output))
