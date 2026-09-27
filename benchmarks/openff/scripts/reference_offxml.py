"""Freeze independent Toolkit interpretation of supported OFFXML generalizations.

Uses the externally supplied, pinned Rosemary rules. Variants only remove default
attributes/optional sections, change equivalent encodings, or make an explicit
parameter override. No Kekule executable participates in reference generation.
"""
import argparse
import copy
import gzip
import hashlib
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from paths import HERE


def variants(source):
    root = ET.fromstring(source)
    yield 'original', root
    defaults = copy.deepcopy(root)
    for name in ('Bonds', 'Angles', 'ProperTorsions', 'ImproperTorsions', 'vdW'):
        node = defaults.find(name)
        for key in ('potential', 'combining_rules', 'default_idivf',
                    'fractional_bondorder_method', 'fractional_bondorder_interpolation'):
            node.attrib.pop(key, None)
    for name in ('vdW', 'Electrostatics'):
        node = defaults.find(name)
        for key in list(node.attrib):
            if key != 'version':
                del node.attrib[key]
    yield 'specification-defaults', defaults
    legacy = copy.deepcopy(defaults)
    for name in ('Bonds', 'ProperTorsions', 'vdW', 'Electrostatics'):
        legacy.find(name).set('version', '0.3')
    legacy.find('vdW').set('method', 'cutoff')
    legacy.find('Electrostatics').set('method', 'Coulomb')
    yield 'legacy-03-coulomb', legacy
    legacy_default = copy.deepcopy(legacy)
    legacy_default.find('Electrostatics').attrib.pop('method')
    yield 'legacy-03-default-pme', legacy_default
    optional = copy.deepcopy(defaults)
    for name in ('Constraints', 'ImproperTorsions', 'LibraryCharges'):
        optional.remove(optional.find(name))
    yield 'optional-sections-absent', optional
    compact = copy.deepcopy(root)
    for parameter in compact.iter():
        for key, value in list(parameter.attrib.items()):
            if ' * ' in value:
                parameter.set(key, value.replace(' ', ''))
    compact.find('Bonds').set('potential', '(k/2)*(r-length)^2')
    bond = compact.find('Bonds/Bond')
    bond.attrib.pop('id')
    bond.set('k', '2e0*kilocalorie/(mole*angstrom**2)')
    bond.set('length', '.15*nanometers')
    yield 'compact-units-and-anonymous-override', compact


def freeze(source):
    import openff.toolkit
    from openff.toolkit import ForceField
    from openff.units import unit
    units = dict(
        Bonds=dict(length=unit.nanometer, k=unit.kilojoule_per_mole / unit.nanometer**2),
        Angles=dict(angle=unit.radian, k=unit.kilojoule_per_mole / unit.radian**2),
        vdW=dict(sigma=unit.nanometer, epsilon=unit.kilojoule_per_mole),
        Constraints=dict(distance=unit.nanometer))
    records = []
    for name, root in variants(source):
        xml = ET.tostring(root, encoding='unicode')
        ff = ForceField(io.StringIO(xml))
        parameters = {}
        for handler in ('Bonds', 'Angles', 'ProperTorsions', 'ImproperTorsions',
                        'Constraints', 'vdW', 'LibraryCharges'):
            parameters[handler] = []
            if handler not in ff.registered_parameter_handlers:
                continue
            for p in ff[handler].parameters:
                item = dict(id=p.id or '', smirks=p.smirks)
                if handler in units:
                    item.update({k: None if getattr(p, k) is None else float(getattr(p, k).m_as(u))
                                 for k, u in units[handler].items()})
                elif handler == 'LibraryCharges':
                    item['charges'] = [float(q.m_as(unit.elementary_charge)) for q in p.charge]
                else:
                    item['terms'] = [dict(k=float(k.m_as(unit.kilojoule_per_mole)),
                                          phase=float(phase.m_as(unit.radian)), periodicity=int(n),
                                          idivf=float(p.idivf[i]) if p.idivf is not None else
                                          (3.0 if handler == 'ImproperTorsions' else 0.0))
                                     for i, (k, phase, n) in enumerate(zip(p.k, p.phase, p.periodicity, strict=True))]
                parameters[handler].append(item)
        vdw, electro = ff['vdW'], ff['Electrostatics']
        settings = {}
        for prefix, node in [('vdw', vdw), ('electrostatics', electro)]:
            for key in ('cutoff', 'switch_width'):
                settings[f'{prefix}_{key}'] = float(getattr(node, key).m_as(unit.nanometer))
            settings[f'{prefix}_scales'] = [float(getattr(node, f'scale1{i}')) for i in range(2, 6)]
        settings.update(vdw_periodic_method=vdw.periodic_method,
                        vdw_nonperiodic_method=vdw.nonperiodic_method,
                        electrostatics_periodic_method=electro.periodic_potential,
                        electrostatics_nonperiodic_method=electro.nonperiodic_potential)
        records.append(dict(name=name, xml=xml, parameters=parameters, settings=settings))
    return dict(schema=1, toolkit_version=openff.toolkit.__version__,
                source_sha256=hashlib.sha256(source).hexdigest(), records=records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    source = (HERE.parent.parent / 'crates/kekule-openff/data/rosemary.offxml').read_bytes()
    report = freeze(source)
    with args.output.open('xb') as output:
        output.write(gzip.compress((json.dumps(report, allow_nan=False) + '\n').encode(), mtime=0))
    print(json.dumps(dict(cases=len(report['records']), source_sha256=report['source_sha256'],
                         sha256=hashlib.sha256(args.output.read_bytes()).hexdigest())))
