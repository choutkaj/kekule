"""Cartesian gradient validation against OpenMM and numerical energy derivatives.

Uses the frozen 110-molecule Hamiltonians and coordinates. Forces are negated
exactly once at the OpenMM boundary: all stored gradients are dE/dx, kJ/mol/nm.
"""
import argparse
import gzip
import importlib.metadata
import json
import math
from pathlib import Path

from paths import HERE, ROOT
from robustness import HANDLERS, Observer, digest, read, split_system

COMPONENTS = (*HANDLERS, 'Total')
ATOL, RTOL = 1e-5, 1e-10
FD_ATOL, FD_RTOL = 2e-3, 2e-6
ROTATION = [[.36, -.48, .8], [.8, .6, 0.], [-.48, .64, .6]]


def freeze(path, report):
    with path.open('xb') as output:
        output.write(gzip.compress((json.dumps(report, allow_nan=False) + '\n').encode(), mtime=0))


def reference(args):
    import numpy as np
    import openmm as mm
    from openmm import unit
    records = []
    for case in read(args.energy_reference)['records'][:args.limit]:
        if args.case and case['id'] != args.case: continue
        omm = mm.XmlSerializer.deserialize(case['openmm_system_xml'])
        split = split_system(omm, len(case['parameters']['ProperTorsions']), len(case['parameters']['ImproperTorsions']))
        integrator = mm.VerletIntegrator(.001)
        context = mm.Context(split, integrator, mm.Platform.getPlatformByName('Reference'))
        frames = []
        for xyz, energies in zip(case['coordinates_nm'], case['energies'], strict=True):
            context.setPositions(np.asarray(xyz) * unit.nanometer)
            gradients = {}
            for i, component in enumerate(COMPONENTS):
                state = context.getState(getEnergy=True, getForces=True, groups=-1 if component == 'Total' else 1 << i)
                energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                if not math.isclose(energy, energies[component], abs_tol=1e-7, rel_tol=1e-10):
                    raise ValueError(f'Frozen energy changed: {case["id"]}, {component}')
                gradients[component] = (-state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)).tolist()
            frames.append(gradients)
        del context, integrator
        records.append(dict(id=case['id'], gradients=frames))
        print(case['id'], 'reference gradients complete', flush=True)
    if not records: raise ValueError('No matching reference cases')
    freeze(args.output, dict(schema=1, convention='dE/dx; kJ/mol/nm',
        energy_reference_sha256=digest(args.energy_reference), openmm=importlib.metadata.version('openmm'), records=records))


def comparison(expected, actual, atol=ATOL, rtol=RTOL, bound=0.):
    """Reject missing/nonfinite components as well as numerical mismatches."""
    import numpy as np
    a, b = np.asarray(expected, dtype=float), np.asarray(actual, dtype=float)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Invalid gradient shape or nonfinite values')
    error = np.abs(b - a)
    threshold = atol + rtol * np.abs(a) + bound
    return dict(values=a.size, max_abs=float(error.max()), sum_squared=float(np.sum(error**2)),
                failed=int(np.sum(error > threshold)))


def directions(n):
    import numpy as np
    rng = np.random.default_rng(20260927 + n)
    # Coordinate derivatives span both terminal and interior atoms; dense
    # directions exercise the full sum, including every protein atom.
    result = []
    for atom in (0, n // 2, n - 1):
        for axis in range(3):
            d = np.zeros((n, 3)); d[atom, axis] = 1.
            result.append(d)
    for _ in range(2):
        d = rng.normal(size=(n, 3)); d /= np.linalg.norm(d)
        result.append(d)
    return result


def charge_bound(xyz, reference, actual, exceptions):
    import numpy as np
    q, r = np.asarray(actual), np.asarray(reference)
    scales = {tuple(sorted(p['maps'])): p['electrostatics_scale'] for p in exceptions}
    result = np.zeros_like(xyz)
    for i in range(len(q)):
        displacement = xyz[i] - xyz[i + 1:]
        distance = np.linalg.norm(displacement, axis=1)
        weights = np.array([scales.get((i + 1, j + 1), 1.) for j in range(i + 1, len(q))])
        product_error = np.abs(q[i] * q[i + 1:] - r[i] * r[i + 1:])
        pair = 138.93545764438198 * (weights * product_error / distance**3)[:, None] * np.abs(displacement)
        result[i] += pair.sum(axis=0)
        result[i + 1:] += pair
    return result


def linear_angles(case, xyz):
    """Record geometry inside OpenMM Reference's angle cross-product floor."""
    import numpy as np
    result = []
    for parameter in case['parameters']['Angles']:
        a, b, c = (i - 1 for i in parameter['atoms'])
        u, v = xyz[a] - xyz[b], xyz[c] - xyz[b]
        cross_norm = float(np.linalg.norm(np.cross(u, v)))
        if cross_norm < 1e-6:
            result.append(dict(atoms=parameter['atoms'], cross_norm_nm2=cross_norm,
                sine=float(cross_norm / np.linalg.norm(u) / np.linalg.norm(v)), parameters=parameter['values']))
    return result


def openmm_differences(case, xyz, probes, steps):
    """Independently test whether reference forces differentiate reference energies."""
    import openmm as mm
    import numpy as np
    from openmm import unit
    omm = mm.XmlSerializer.deserialize(case['openmm_system_xml'])
    split = split_system(omm, len(case['parameters']['ProperTorsions']), len(case['parameters']['ImproperTorsions']))
    integrator = mm.VerletIntegrator(.001)
    context = mm.Context(split, integrator, mm.Platform.getPlatformByName('Reference'))
    result = []
    for direction in probes:
        values = []
        for h in steps:
            pair = []
            for sign in (-1., 1.):
                context.setPositions((xyz + sign * h * direction) * unit.nanometer)
                pair.append({c:float(context.getState(getEnergy=True, groups=-1 if c == 'Total' else 1 << i).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)) for i, c in enumerate(COMPONENTS)})
            values.append(dict(step_nm=h, minus=pair[0], plus=pair[1]))
        result.append(values)
    context.setPositions((xyz @ np.asarray(ROTATION).T + [1.25, -.75, .5]) * unit.nanometer)
    rotated = {c:(-context.getState(getForces=True, groups=-1 if c == 'Total' else 1 << i).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)).tolist() for i, c in enumerate(COMPONENTS)}
    del context, integrator
    return result, rotated


def compare(args):
    import numpy as np
    from rdkit import Chem
    reference = read(args.reference)
    if reference['energy_reference_sha256'] != digest(args.energy_reference):
        raise ValueError('Energy reference fingerprint mismatch')
    energy = read(args.energy_reference)
    inputs = read(args.inputs)
    if energy['inputs_sha256'] != digest(args.inputs): raise ValueError('Input fingerprint mismatch')
    reference_by_id = {r['id']: r for r in reference['records']}
    inputs_by_id = {r['id']: r for r in inputs['records']}
    records = []
    observer = Observer(args.binary, args.model)
    rotation = np.asarray(ROTATION)
    try:
        for case in energy['records'][:args.limit]:
            if args.case and case['id'] != args.case: continue
            expected = reference_by_id[case['id']]['gradients']
            frames = np.asarray(case['coordinates_nm'])
            n = len(case['charges'])
            probes = directions(n)
            geometry = [linear_angles(case, xyz) for xyz in frames]
            steps_nm = [1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8] if any(geometry) else [1e-5, 5e-6]
            original = None
            for reverse in (False, True):
                row = dict(id=case['id'], reverse=reverse, errors=[], comparisons=[], invariants=[], finite_differences=[])
                try:
                    params = Chem.SmilesParserParams(); params.removeHs = False
                    rd = Chem.MolFromSmiles(inputs_by_id[case['id']]['smiles'], params)
                    if reverse: rd = Chem.RenumberAtoms(rd, list(reversed(range(n))))
                    # Rotated/translated frames test the transformation of vectors,
                    # not merely invariant scalar energies or a zero total force.
                    transformed = frames @ rotation.T + [1.25, -.75, .5]
                    request = dict(smiles=Chem.MolToSmiles(rd, canonical=False, allHsExplicit=True), mode='validate',
                        coordinates_nm=np.concatenate((frames, transformed)).tolist(), reference_charges=case['charges'], gradients=True)
                    native = observer.call(request, args.timeout)
                    if native['status'] != 'ok': raise ValueError(native['message'])
                    maps = native['maps']
                    if sorted(maps) != list(range(1, n + 1)): raise ValueError('Invalid atom correspondence')
                    order = np.argsort(maps)
                    row['maps'] = maps
                    row['observations'] = native['energies'][:3]
                    row['near_linear_angles'] = geometry
                    if len(native['energies']) != 6 or len(expected) != 3: raise ValueError('Missing geometries')
                    charges = np.asarray(native['system']['charges'])[order]
                    row['charges'] = charges.tolist()
                    row['exceptions'] = native['system']['exceptions']
                    for frame, (xyz, ref) in enumerate(zip(frames, expected, strict=True)):
                        bound = charge_bound(xyz, case['charges'], charges, native['system']['exceptions'])
                        order_bound = charge_bound(xyz, np.asarray(original['system']['charges'])[np.argsort(original['maps'])], charges, native['system']['exceptions']) if reverse else 0.
                        for source in ('reference_charges', 'native_charges'):
                            obs = native['energies'][frame][source]
                            rotated = native['energies'][frame + 3][source]
                            if set(obs['gradients']) != set(COMPONENTS): raise ValueError('Missing gradient components')
                            for component in COMPONENTS:
                                g = np.asarray(obs['gradients'][component])[order]
                                if source == 'reference_charges' and not math.isclose(obs[component], case['energies'][frame][component], abs_tol=1e-7, rel_tol=1e-10):
                                    row['errors'].append(dict(check='frozen energy', frame=frame, component=component))
                                stats = comparison(ref[component], g, bound=bound if source == 'native_charges' and component in ('Total', 'Electrostatics') else 0.)
                                row['comparisons'].append(dict(frame=frame, source=source, component=component, **stats))
                                if stats['failed']: row['errors'].append(dict(check='OpenMM', frame=frame, source=source, component=component, **stats))
                                stats = comparison(g @ rotation.T, np.asarray(rotated['gradients'][component])[order])
                                if stats['failed']: row['errors'].append(dict(check='rotation', frame=frame, source=source, component=component, **stats))
                                if not math.isclose(obs[component], rotated[component], abs_tol=1e-7, rel_tol=1e-10):
                                    row['errors'].append(dict(check='rigid energy', frame=frame, component=component))
                                centered = xyz - xyz.mean(axis=0)
                                torque_terms = np.cross(centered, g)
                                force_residual, torque_residual = np.abs(g.sum(axis=0)), np.abs(torque_terms.sum(axis=0))
                                row['invariants'].append(dict(frame=frame, source=source, component=component,
                                    net_gradient=float(force_residual.max()), net_torque=float(torque_residual.max()), rotation=stats))
                                if np.any(force_residual > 1e-7 + 1e-12 * np.abs(g).sum(axis=0)) or np.any(torque_residual > 1e-7 + 1e-12 * np.abs(torque_terms).sum(axis=0)):
                                    row['errors'].append(dict(check='net force/torque', frame=frame, component=component))
                                if reverse:
                                    old = np.asarray(original['energies'][frame][source]['gradients'][component])[np.argsort(original['maps'])]
                                    stats = comparison(old, g, bound=order_bound if source == 'native_charges' and component in ('Total', 'Electrostatics') else 0.)
                                    row['invariants'][-1]['atom_order'] = stats
                                    if stats['failed']: row['errors'].append(dict(check='atom order', frame=frame, source=source, component=component, **stats))
                    if not reverse:
                        original = native
                        # Separate request avoids finite differences on transformed
                        # and reversed duplicates, while retaining all 330 base frames.
                        request.update(coordinates_nm=frames.tolist(), directions=[d.tolist() for d in probes], finite_difference_steps_nm=steps_nm)
                        fd = observer.call(request, args.timeout)
                        if fd['status'] != 'ok': raise ValueError(fd['message'])
                        for frame, observation in enumerate(fd['energies']):
                            for direction, steps in zip(probes, observation['finite_differences'], strict=True):
                                for component in COMPONENTS:
                                    estimates = [(s['plus'][component] - s['minus'][component]) / (2 * s['step_nm']) for s in steps]
                                    extrapolated = [(4 * estimates[i + 1] - estimates[i]) / 3 for i in range(0, len(estimates), 2)]
                                    numerical = extrapolated[-1]
                                    analytic = float(np.sum(np.asarray(native['energies'][frame]['reference_charges']['gradients'][component])[order] * direction))
                                    stats = comparison([analytic], [numerical], FD_ATOL, FD_RTOL)
                                    row['finite_differences'].append(dict(frame=frame, component=component, analytic=analytic, numerical=numerical,
                                        central=estimates, richardson=extrapolated, steps_nm=steps_nm, **stats))
                                    if stats['failed']: row['errors'].append(dict(check='finite difference', frame=frame, component=component, **stats))
                                    if len(extrapolated) > 1 and comparison([extrapolated[-1]], [extrapolated[-2]], FD_ATOL, FD_RTOL)['failed']:
                                        row['errors'].append(dict(check='finite difference convergence', frame=frame, component=component))
                        row['finite_difference_observations'] = [r['finite_differences'] for r in fd['energies']]
                        row['openmm_finite_differences'] = []
                        for frame, angles in enumerate(geometry):
                            if not angles: continue
                            values, rotated = openmm_differences(case, frames[frame], probes, steps_nm)
                            checks = []
                            for direction, steps in zip(probes, values, strict=True):
                                for component in COMPONENTS:
                                    estimates = [(s['plus'][component] - s['minus'][component]) / (2 * s['step_nm']) for s in steps]
                                    numerical = (4 * estimates[-1] - estimates[-2]) / 3
                                    analytic = float(np.sum(np.asarray(expected[frame][component]) * direction))
                                    native_analytic = float(np.sum(np.asarray(native['energies'][frame]['reference_charges']['gradients'][component])[order] * direction))
                                    checks.append(dict(component=component, numerical=numerical, openmm_analytic=analytic, native_analytic=native_analytic,
                                        reference=comparison([analytic], [numerical], FD_ATOL, FD_RTOL), native=comparison([native_analytic], [numerical], FD_ATOL, FD_RTOL)))
                            row['openmm_finite_differences'].append(dict(frame=frame, observations=values, checks=checks,
                                rotated_gradients=rotated, rotation={c:comparison(np.asarray(expected[frame][c]) @ rotation.T, rotated[c]) for c in COMPONENTS}))
                except Exception as exc:
                    row['errors'].append(str(exc))
                row['passed'] = not row['errors']
                records.append(row)
                print(json.dumps(dict(id=row['id'], reverse=reverse, passed=row['passed'], errors=row['errors'])), flush=True)
    finally:
        observer.close()
    if not records: raise ValueError('No matching native cases')
    report = dict(schema=1, convention='dE/dx; kJ/mol/nm', reference_sha256=digest(args.reference),
        energy_reference_sha256=digest(args.energy_reference), inputs_sha256=digest(args.inputs), binary_sha256=digest(args.binary),
        model_sha256={p:digest(args.model / p) for p in ('model.json', 'weights.bin')},
        tolerances=dict(gradient_atol=ATOL, gradient_rtol=RTOL, finite_difference_atol=FD_ATOL, finite_difference_rtol=FD_RTOL),
        method='11 directions (9 coordinate axes, 2 seeded dense unit vectors); central differences at 1e-5 and 5e-6 nm, Richardson extrapolation; fixed charges.',
        rotation=ROTATION, records=records)
    report['source_sha256'] = {p: digest(ROOT / p) for p in ('benchmarks/src/bin/openff_parameterize/energy.rs', 'benchmarks/openff/scripts/gradients.py', 'benchmarks/openff/scripts/robustness.py', 'Cargo.lock')}
    report['method'] += ' Near-linear angles (cross norm < 1e-6 nm^2) additionally use 1e-6/5e-7 and 1e-7/5e-8 nm pairs; require the last two extrapolations to agree. OpenMM energies are independently differentiated there.'
    freeze(args.output, report)
    return int(any(not r['passed'] for r in records))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('reference', 'compare'))
    parser.add_argument('--inputs', type=Path, default=HERE/'data/inputs.json.gz')
    parser.add_argument('--energy-reference', type=Path, default=HERE/'data/reference.json.gz')
    parser.add_argument('--reference', type=Path, default=HERE/'data/gradients.json.gz')
    parser.add_argument('--binary', type=Path, default=ROOT/'target/release/openff_parameterize.exe')
    parser.add_argument('--model', type=Path, default=ROOT/'target/openff-models/openff-gnn-am1bcc-1.0.0')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--limit', type=int)
    parser.add_argument('--case', help='Inspect one named input without changing the frozen panel')
    parser.add_argument('--timeout', type=float, default=600)
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1: parser.error('--limit must be positive')
    if args.output.exists(): raise FileExistsError(args.output)
    if args.mode == 'reference': reference(args)
    else: raise SystemExit(compare(args))
