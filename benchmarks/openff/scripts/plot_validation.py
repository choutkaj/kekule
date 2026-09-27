"""Render publication figures from frozen observations, without running either engine.

Only rendering needs NumPy/Matplotlib. Extraction and integrity checks use the
standard library, so their regression tests can run without scientific packages.
"""
import argparse
import collections
import hashlib
import json
import math
import statistics
from pathlib import Path

from robustness import canonical, read, summarize

from paths import HERE
HANDLERS = dict(Bonds='bonds', Angles='angles', Constraints='constraints',
                ProperTorsions='propers', ImproperTorsions='impropers', vdW='vdw')
COMPONENTS = ('Bonds', 'Angles', 'ProperTorsions', 'ImproperTorsions',
              'vdW', 'Electrostatics', 'Total')
COLORS = ('#0072B2', '#D55E00')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metrics(pairs):
    errors = [y - x for _, x, y in pairs]
    require(errors and all(math.isfinite(v) for v in errors), 'invalid paired values')
    return dict(count=len(errors), max_abs=max(map(abs, errors)),
                rmse=math.sqrt(math.fsum(v * v for v in errors) / len(errors)),
                mean_signed=math.fsum(errors) / len(errors))


def collect(reference, native):
    """Pair by atom-map identity; retain every atom order, torsion and geometry."""
    refs = {r['id']: r for r in reference['records']}
    require(len(refs) == len(reference['records']), 'duplicate reference ID')
    require(reference['inputs_sha256'] == native['inputs_sha256'], 'input fingerprint mismatch')
    seen = set()
    series = collections.defaultdict(list)
    energies = collections.defaultdict(list)
    for row in native['records']:
        key = row['id'], row['permutation']
        require(key not in seen, 'duplicate native case')
        seen.add(key)
        ref, obs = refs[row['id']], row['native']
        require(ref['status'] == obs['status'] == 'ok', 'missing successful observation')
        group = int(row['kind'] == 'protein-chain')
        require(row['kind'] == ref['kind'], 'mismatched molecular category')
        maps = obs['maps']
        require(len(maps) == len(set(maps)) == len(ref['charges']), 'invalid atom mapping')
        require(set(maps) == set(range(1, len(maps) + 1)), 'incomplete atom mapping')
        index = {atom: i for i, atom in enumerate(maps)}
        for atom, charge in enumerate(ref['charges'], 1):
            i = index[atom]
            series['charge'].append((group, charge, obs['system']['charges'][i]))
            for x, y in zip(ref['features'][atom - 1], obs['features']['values'][i], strict=True):
                series['feature'].append((group, x, y))
        for handler, field in HANDLERS.items():
            expected, actual = collections.defaultdict(list), collections.defaultdict(list)
            for item in ref['parameters'][handler]:
                expected[(canonical(item['atoms'], handler), item['mult'])].append(item['values'])
            for i, item in enumerate(obs['system'][field]):
                atoms = [maps[i]] if field == 'vdw' else item['maps']
                param = item if field == 'vdw' else item['parameter']
                terms = enumerate(param['terms']) if 'terms' in param else [(None, param)]
                for mult, values in terms:
                    actual[(canonical(atoms, handler), mult)].append(
                        {k: v for k, v in values.items() if k != 'id'})
            require(expected.keys() == actual.keys(), f'{handler}: missing parameter identities')
            sortkey = lambda item: tuple(sorted(item.items()))
            for key in expected:
                a, b = sorted(expected[key], key=sortkey), sorted(actual[key], key=sortkey)
                require(len(a) == len(b), f'{handler}: multiplicity mismatch')
                for x, y in zip(a, b, strict=True):
                    require(x.keys() == y.keys(), f'{handler}: missing parameter fields')
                    for name in x:
                        series[f'{handler}.{name}'].append((group, x[name], y[name]))
        for ref_frame, obs_frame in zip(ref['energies'], obs['energies'], strict=True):
            for mode in ('reference_charges', 'native_charges'):
                for component in COMPONENTS:
                    energies[f'{mode}.{component}'].append(
                        (group, ref_frame[component], obs_frame[mode][component]))
    require(seen == {(key, order) for key in refs for order in ('original', 'reversed')},
            'incomplete molecule/atom-order coverage')
    return dict(series), dict(energies)


def save(fig, output, stem):
    for extension in ('png', 'svg', 'pdf'):
        metadata = {'Date': None} if extension == 'svg' else (
            {'CreationDate': None, 'ModDate': None} if extension == 'pdf' else {})
        fig.savefig(output / f'{stem}.{extension}', dpi=300, metadata=metadata,
                    facecolor='white', bbox_inches='tight')


def render(series, energies, timing, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.labelsize': 9,
        'axes.titlesize': 10, 'axes.titleweight': 'bold',
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': .7, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'legend.fontsize': 8, 'legend.frameon': False,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'svg.hashsalt': 'kekule-openff-validation-2026-09-22',
        'axes.formatter.use_mathtext': True,
    })
    output.mkdir(parents=True, exist_ok=True)
    handles = [Line2D([], [], linestyle='', marker='o', color=c, markersize=4, label=l)
               for c, l in zip(COLORS, ('100 small molecules', '10 protein chains'))]

    def points(ax, data, residual=False, scale=1.):
        data = np.asarray(data)
        for group, color in enumerate(COLORS):
            part = data[data[:, 0] == group]
            ax.scatter(part[:, 1], (part[:, 2] - part[:, 1]) / scale if residual else part[:, 2],
                       color=color, s=8, alpha=.38, linewidths=0, rasterized=True)
        ax.grid(alpha=.15, linewidth=.5)
        ax.set_axisbelow(True)

    fig, axes = plt.subplots(2, 3, figsize=(9.4, 5.8), gridspec_kw={'height_ratios': [1.45, 1]})
    specs = [('vdW.sigma', r'Lennard–Jones $\sigma$', r'$\sigma$ (nm)', 1e-17, r'$\Delta\sigma$ ($10^{-17}$ nm)'),
             ('vdW.epsilon', r'Lennard–Jones $\epsilon$', r'$\epsilon$ (kJ mol$^{-1}$)', 1., r'$\Delta\epsilon$ (kJ mol$^{-1}$)'),
             ('charge', 'Partial charge', r'$q$ ($e$)', 1e-7, r'$\Delta q$ ($10^{-7}\ e$)')]
    for col, (key, title, label, scale, delta_label) in enumerate(specs):
        data = series[key]
        top, bottom = axes[:, col]
        points(top, data)
        bounds = [min(min(x, y) for _, x, y in data), max(max(x, y) for _, x, y in data)]
        pad = .08 * (bounds[1] - bounds[0])
        bounds = [bounds[0] - pad, bounds[1] + pad]
        top.plot(bounds, bounds, '--', color='#555555', lw=.8, zorder=0)
        top.set(xlim=bounds, ylim=bounds, xlabel='OpenFF ' + label, ylabel='Kekule ' + label)
        top.set_aspect('equal', adjustable='box')
        top.set_title(title, pad=9)
        stat = metrics(data)
        top.text(.04, .96, r"max $|\Delta|$ = " + f"{stat['max_abs']:.2e}", transform=top.transAxes,
                 va='top', fontsize=8, bbox=dict(facecolor='white', edgecolor='none', alpha=.85))
        points(bottom, data, residual=True, scale=scale)
        bottom.axhline(0, color='#555555', lw=.7, zorder=0)
        bottom.set(xlim=bounds, xlabel='OpenFF ' + label, ylabel=delta_label)
        extent = max(abs(y - x) for _, x, y in data) / scale
        if extent == 0:
            bottom.set_ylim(-.5, .5)
            bottom.set_yticks([0])
            bottom.text(.5, .78, 'Exact agreement', transform=bottom.transAxes, ha='center', fontsize=8)
        else:
            bottom.set_ylim(-extent * 1.2, extent * 1.2)
        for row, ax in enumerate((top, bottom)):
            ax.text(-.17, 1.06, 'abcdef'[row * 3 + col], transform=ax.transAxes,
                    fontweight='bold', fontsize=12)
            ax.xaxis.set_major_locator(MaxNLocator(4))
    fig.legend(handles=handles, loc='lower center', ncol=2, bbox_to_anchor=(.5, -.015))
    fig.subplots_adjust(left=.085, right=.99, bottom=.12, top=.92, hspace=.65, wspace=.50)
    save(fig, output, 'parameter-parity')
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.4), gridspec_kw={'height_ratios': [1.3, 1]})
    for col, (mode, title, scale, unit) in enumerate([
        ('reference_charges', 'Identical charges', 1e-9, r'$10^{-9}$'),
        ('native_charges', 'Independently assigned charges', 1e-4, r'$10^{-4}$')]):
        data = energies[f'{mode}.Total']
        top, bottom = axes[:, col]
        points(top, data)
        bounds = [-1e4, 2e6]
        require(all(bounds[0] < x < bounds[1] and bounds[0] < y < bounds[1]
                    for _, x, y in data), 'energy plot limits would hide observations')
        top.set_xscale('symlog', linthresh=100)
        top.set_yscale('symlog', linthresh=100)
        top.plot(bounds, bounds, '--', color='#555555', lw=.8, zorder=0)
        top.set(xlim=bounds, ylim=bounds, title=title,
                xlabel=r'OpenMM total energy (kJ mol$^{-1}$)', ylabel=r'Kekule total energy (kJ mol$^{-1}$)')
        top.set_aspect('equal', adjustable='box')
        top.text(.04, .96, r"max $|\Delta E|$ = " + f"{metrics(data)['max_abs']:.2e} kJ/mol",
                 transform=top.transAxes, va='top', fontsize=8,
                 bbox=dict(facecolor='white', edgecolor='none', alpha=.85))
        points(bottom, data, residual=True, scale=scale)
        bottom.set_xscale('symlog', linthresh=100)
        bottom.axhline(0, color='#555555', lw=.7, zorder=0)
        bottom.set(xlim=bounds, xlabel=r'OpenMM total energy (kJ mol$^{-1}$)',
                   ylabel=r'$\Delta E$ (' + unit + r' kJ mol$^{-1}$)')
        for row, ax in enumerate((top, bottom)):
            ax.set_xticks([-1e4, -1e2, 0, 1e2, 1e4, 1e6])
            if row == 0:
                ax.set_yticks([-1e4, -1e2, 0, 1e2, 1e4, 1e6])
            ax.text(-.15, 1.06, 'abcd'[row * 2 + col], transform=ax.transAxes,
                    fontweight='bold', fontsize=12)
    fig.legend(handles=handles, loc='lower center', ncol=2, bbox_to_anchor=(.5, -.015))
    fig.subplots_adjust(left=.11, right=.98, bottom=.12, top=.93, hspace=.62, wspace=.42)
    save(fig, output, 'energy-parity')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    rows = sorted(timing['records'], key=lambda r: r['atoms'])
    for key, label, color, marker in [
        ('inference_including_features_ms', 'Forced inference + preparation', COLORS[0], 'o'),
        ('full_parameterization_ms', 'Full parameterization', COLORS[1], 's')]:
        x = [r['atoms'] for r in rows]
        y = [statistics.median(r[key]) for r in rows]
        limits = [[m - min(r[key]) for r, m in zip(rows, y)], [max(r[key]) - m for r, m in zip(rows, y)]]
        ax.errorbar(x, y, yerr=limits, fmt=marker, color=color, capsize=3,
                    markersize=5, elinewidth=.8, label=label)
    ax.set(xscale='log', yscale='log', xlabel='Atoms including hydrogens', ylabel='Warm execution time (ms)')
    ax.grid(alpha=.15, linewidth=.5, which='both')
    ax.set_axisbelow(True)
    ax.legend(loc='upper left')
    ax.annotate('Lookup hit\n(full pipeline)', xy=(8, statistics.median(rows[0]['full_parameterization_ms'])),
                xytext=(17, .65), fontsize=8, arrowprops=dict(arrowstyle='-', lw=.7, color='#555555'))
    fig.tight_layout()
    save(fig, output, 'cpu-timings')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'figures')
    args = parser.parse_args()
    lock = read(HERE / 'results/manifest.json')
    for artifact in lock['artifacts']:
        require(digest(HERE / artifact['path']) == artifact['sha256'],
                'changed frozen artifact: ' + artifact['path'])
    reference, native = read(HERE / 'data/reference.json.gz'), read(HERE / 'results/native.json.gz')
    require(native['reference_sha256'] == digest(HERE / 'data/reference.json.gz'),
            'reference fingerprint mismatch')
    series, energies = collect(reference, native)
    timing = read(HERE / 'results/timings.json')
    render(series, energies, timing, args.output)
    import matplotlib
    import numpy
    summary = dict(
        source_sha256={a['path']: a['sha256'] for a in lock['artifacts']},
        robustness=summarize(read(HERE / 'data/inputs.json.gz'), reference, native),
        timing_sha256=digest(HERE / 'results/timings.json'),
        renderer=dict(matplotlib=matplotlib.__version__, numpy=numpy.__version__),
        parameters={k: metrics(v) for k, v in series.items()},
        energies={k: metrics(v) for k, v in energies.items()},
        cpu_timings=[dict(id=r['id'], atoms=r['atoms'],
                         median_ms={k: statistics.median(r[k]) for k in timing['repetitions']})
                     for r in timing['records']],
        cases=len(native['records']),
        parameterization_passed=sum(r['parameterization_passed'] for r in native['records']),
        all_checks_passed=sum(r['passed'] for r in native['records']))
    (HERE / 'results/summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
