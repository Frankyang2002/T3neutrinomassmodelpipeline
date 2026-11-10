#!/usr/bin/env python3
"""Corrected D(-2) frozen/running comparison with existing plot/physics code.

Compare two independent switches, using the same retained fixed-Yukawa benchmark:
 E0: freeze renormalisable EFT1 couplings at MF, but evolve full Wilson C
 E1: use saved running renormalisable EFT1 couplings and evolve full Wilson C
 S0: freeze SMEFT running (evaluate C5 at MS, carry unchanged to mu_low)
 S1: evolve full-flavour SM+Weinberg from MS to mu_low.

Every scenario integrates out all eight real scalars simultaneously at MS.
UV pre-MF state is the same saved boundary in ALL scenarios; UV is not toggled.

The heatmap uses the EXISTING scripts.PlotT3OscillationHeatmaps.plot_heatmap
renderer. Missing mixing angles are deliberately omitted: they are not inferred.
Conditional numerical candidates only; not certified physical one-loop matching.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from RunDMinus2CorrectedEFT1 import prepare
from RunDMinus2FinalRoughMasses import native_sm_at_ms, spectrum
from Numerical.parallel_t3.wilson.interchange import Background
from Numerical.parallel_t3.wilson.workflow import run_wilson
from Numerical.parallel_t3.threshold.matching import (
    match_scalar_threshold, scalar_gauge_ward_residuals,
)
from Numerical.running.SMWeinbergEvolution import (
    SMWeinbergInitialConditions, evolve_sm_weinberg,
)
from scripts.PlotT3OscillationHeatmaps import format_cell_value
import matplotlib.colors as mcolors

DEFAULT_INPUT = Path('output/parallel_t3/rough_inputs/fit16_cc2dcbec4a')
DEFAULT_OUT = Path('output/parallel_t3/D_m2_corrected_frozen_running')


def frozen_background(basis, bg):
    n = len(bg.t)
    return Background(
        basis,
        bg.t.copy(),
        np.repeat(bg.y[-1:], n, axis=0),
        np.repeat(bg.quartic[-1:], n, axis=0),
        np.repeat(bg.gauge_couplings[-1:], n, axis=0),
        bg.generators,
        {'background': 'MF renormalisable couplings frozen; Wilson ADM active'},
    )


def record_spectrum(label, C, v):
    sp = spectrum(C, v)
    return dict(label=label, C5=C, mass=sp['masses_eV'],
                dm21=sp['dm21_eV2'], dm31=sp['dm31_eV2'])


def match_one(result, basis, qms, mass2, rotation, ms):
    matched = match_scalar_threshold(
        full_wilson_ms=result['C_final'],
        lepton_indices=np.asarray(basis.leptons, int),
        higgs_indices=np.asarray(basis.higgs, int),
        heavy_indices=np.arange(4, basis.ns, dtype=int),
        full_quartic_ms=qms,
        scalar_mass2_background=mass2,
        R_wilson_from_background=rotation,
        mu_s=ms,
        mode='rge_improved',
        separately_supplied_c5_ms=result['C5_final'],
        projection_tolerance=1e-8,
        enforce_pure_weinberg=True,
    )
    if matched['scalar_components_removed'] != 8 or matched['bubble']['pair_count'] != 36:
        raise RuntimeError('Incomplete simultaneous S1/S2 threshold')
    return matched


def smeft_initial(sm, K):
    return SMWeinbergInitialConditions(
        gY=sm['gY'], g2=sm['g2'], g3=sm['g3'], lambdaH=sm['lambdaH'],
        ye=sm['ye'], yu=sm['yu'], yd=sm['yd'], K=K)


def curve_smeft(sm, Kms, ms, mu_low, npoints, vev):
    # Reuse public existing solver; each target is an independent integration
    # from the same boundary, not an interpolated/invented mass trajectory.
    grid = np.geomspace(mu_low, ms, npoints)
    masses = []
    dm21 = []
    dm31 = []
    for scale in grid:
        if np.isclose(scale, ms, rtol=0, atol=1e-10):
            k = Kms
        else:
            k = evolve_sm_weinberg(
                smeft_initial(sm, Kms), ms, float(scale),
                rtol=1e-8, atol=1e-20).K
        obs = spectrum(k, vev)
        masses.append(obs['masses_eV'])
        dm21.append(obs['dm21_eV2'])
        dm31.append(obs['dm31_eV2'])
    return grid, np.asarray(masses), np.asarray(dm21), np.asarray(dm31)


def save_line(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=240, bbox_inches='tight')
    plt.close(fig)


def make_figures(out, rows, trajectories, ms, mu_low):
    # 2x2 scenario mass spectrum (at the SAME endpoint, including S0 as
    # explicitly unevolved at MS). This reuses standard repository Matplotlib
    # plotting conventions and no fabricated intermediate points.
    names = ('E0_S0', 'E0_S1', 'E1_S0', 'E1_S1')
    x = np.arange(1, 4)
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for name in names:
        ax.plot(x, rows[name]['mass'], marker='o', lw=1.8,
                label=name.replace('_', ' / '))
    ax.set_xticks(x, [r'$m_1$', r'$m_2$', r'$m_3$'])
    ax.set_ylabel('Neutrino mass [eV]')
    ax.set_title('D(-2): corrected frozen versus running — low-scale comparison')
    ax.grid(alpha=0.25)
    ax.legend(title='EFT1 / SMEFT')
    save_line(fig, out / 'D_m2_frozen_running_masses.png')

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for name in names:
        values = [rows[name]['dm21'], rows[name]['dm31']]
        ax.plot([1, 2], values, marker='o', lw=1.8,
                label=name.replace('_', ' / '))
    ax.set_xticks([1, 2], [r'$\Delta m^2_{21}$', r'$\Delta m^2_{31}$'])
    ax.set_yscale('log')
    ax.set_ylabel(r'Mass-squared splitting [eV$^2$]')
    ax.set_title('D(-2): corrected mass-squared splittings')
    ax.grid(alpha=0.25, which='both')
    ax.legend(title='EFT1 / SMEFT')
    save_line(fig, out / 'D_m2_frozen_running_splittings.png')

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for label, (mu, masses, _, _) in trajectories.items():
        ax.plot(mu, masses[:, 2], lw=1.8,
                label=label + r' $m_3$')
    # Exact no-SMEFT-running S0 endpoint values are horizontal across
    # (MS,mu_low), but are deliberately not plotted as dynamical trajectories.
    ax.set_xscale('log')
    ax.axvline(ms, ls='--', color='0.5', lw=1.0, label=r'$M_S$')
    ax.set_xlabel(r'Renormalisation scale $\mu$ [GeV]')
    ax.set_ylabel(r'$m_3$ [eV]')
    ax.set_title('D(-2): SMEFT mass running, EFT1 on versus frozen')
    ax.grid(alpha=0.25, which='both')
    ax.legend()
    save_line(fig, out / 'D_m2_smeft_m3_trajectory.png')

    # Use existing repo heatmap renderer; unlike original 16-model heatmap
    # the rows here represent 4 physically evaluated D(-2) scenarios.
    # Reference fully frozen E0_S0; signed % = (scenario/reference - 1)*100.
    matrix = np.zeros((4, 5), dtype=float)
    keys = ('m1', 'm2', 'm3', 'dm21', 'dm31')
    ref = rows['E0_S0']
    ref_values = np.r_[ref['mass'], ref['dm21'], ref['dm31']]
    for i, name in enumerate(names):
        v = np.r_[rows[name]['mass'], rows[name]['dm21'], rows[name]['dm31']]
        if np.any(np.abs(ref_values) < 1e-100):
            raise ValueError('Reference observables vanish; cannot form percentage heatmap')
        matrix[i] = 100.0 * (v - ref_values) / ref_values
    vmax = max(float(np.max(np.abs(matrix))), 0.01)
    fig, ax = plt.subplots(figsize=(9.3, 4.7))
    im = ax.imshow(matrix, aspect='auto', cmap='RdBu_r',
                   norm=mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    ax.set_xticks(range(5), [r'$m_1$', r'$m_2$', r'$m_3$',
                              r'$\Delta m^2_{21}$', r'$\Delta m^2_{31}$'])
    ax.set_yticks(range(4), [x.replace('_', ' / ') for x in names])
    ax.set_title('D(-2): signed % change versus E0 / S0 (fully frozen)')
    for i in range(4):
        for j in range(5):
            val = matrix[i, j]
            ax.text(j,i,format_cell_value(val),ha='center',va='center',
                    color='white' if abs(val)>.55*vmax else 'black',fontsize=10)
    ax.set_xticks(np.arange(-.5,5,1), minor=True)
    ax.set_yticks(np.arange(-.5,4,1), minor=True)
    ax.grid(which='minor',color='black',lw=.5)
    ax.tick_params(which='minor',bottom=False,left=False)
    cb = fig.colorbar(im, ax=ax, fraction=.035, pad=.04)
    cb.set_label(r'$100\,(O_{scenario}-O_{E0/S0})/O_{E0/S0}$ [%]')
    save_line(fig,out / 'D_m2_frozen_running_heatmap_percent.png')
    with (out / 'D_m2_frozen_running_heatmap.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['scenario', *keys])
        for i, name in enumerate(names):
            w.writerow([name, *(f'{x:.12g}' for x in matrix[i])])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--mu-low', type=float, default=173.0)
    parser.add_argument('--vev', type=float, default=246.22)
    parser.add_argument('--solver', choices=['DOP853', 'RK2'], default='DOP853')
    parser.add_argument('--steps', type=int, default=4)
    parser.add_argument('--npoints', type=int, default=12)
    parser.add_argument('--rtol', type=float, default=1e-6)
    parser.add_argument('--atol', type=float, default=1e-21)
    args = parser.parse_args()
    if args.steps < 1 or args.npoints < 3 or args.vev <= 0 or args.mu_low <= 0:
        parser.error('steps>=1, npoints>=3, positive vev and mu_low required')

    folder = args.input_root / 'a' / 'D_m2'
    basis, bg, C0, scales, oldq, newq = prepare(folder)
    mf, ms = float(scales[-1]), float(scales[0])
    if not 0 < args.mu_low < ms < mf:
        parser.error('Require mu_low < MS < MF')
    cfile = args.input_root / 'c' / 'D_m2' / 'background.npz'
    with np.load(cfile, allow_pickle=False) as f:
        mass = np.asarray(f['mass_squared_gev2_ms'], float)
        R = np.asarray(f['R_wilson_from_background'], float)
        q_saved = np.asarray(f['quartic_real_full_ms'], float)
    if not np.allclose(R, np.eye(8), atol=1e-10):
        raise ValueError('Unexpected nonidentity saved C-stage scalar rotation')
    if np.linalg.norm(q_saved-oldq[0])/max(np.linalg.norm(q_saved), 1e-100) > 1e-9:
        raise ValueError('Retained A/C scalar quartic mismatch')
    generators = np.asarray([s for group in ('SU2','U1')
                             for _, s in basis.generators()[group]])
    scalar_ward = scalar_gauge_ward_residuals(newq[0], mass, generators)
    if not scalar_ward['passed']:
        raise ValueError('Scalar mass/quartic Ward failed: ' + str(scalar_ward))
    sm = native_sm_at_ms(args.input_root, ms)

    print(f'D(-2) E0/E1 x S0/S1 at MF={mf:g}, MS={ms:g}, mu_low={args.mu_low:g} GeV', flush=True)
    rows = {}
    traj = {}
    for eflag, background in ((0, frozen_background(basis,bg)), (1, bg)):
        print(f'E{eflag}: evolving full Wilson tensor...', flush=True)
        wilson = run_wilson(
            basis, C0, background, mf, ms,
            adm_variant='fonseca_eq22_six_term', solver=args.solver,
            steps=args.steps, rtol=args.rtol, atol=args.atol,
        )
        # E0 means EFT1 renormalisable trajectory frozen; scalar threshold is
        # evaluated at the SAME endpoint with frozen MF quartic, unlike E1.
        threshold_quartic = background.quartic[0]
        # Ensure scalar quartic remains gauge-invariant for both branches.
        check = scalar_gauge_ward_residuals(threshold_quartic, mass, generators)
        if not check['passed']:
            raise ValueError(f'E{eflag}: frozen quartic not covariant: {check}')
        matched = match_one(wilson, basis, threshold_quartic, mass, R, ms)
        Kms = np.asarray(matched['c5_smeft_at_ms'], complex)
        label0 = f'E{eflag}_S0'
        rows[label0] = record_spectrum(label0, Kms, args.vev)
        label1 = f'E{eflag}_S1'
        low = evolve_sm_weinberg(smeft_initial(sm, Kms), ms, args.mu_low,
                                 rtol=1e-8, atol=1e-20)
        rows[label1] = record_spectrum(label1, low.K, args.vev)
        traj[f'E{eflag} / S1'] = curve_smeft(
            sm, Kms, ms, args.mu_low, args.npoints, args.vev)
        print(f'E{eflag}: ||C5_SMEFT(MS)||={np.linalg.norm(Kms):.10e}; '
              f'||C5_low||={np.linalg.norm(low.K):.10e}; '
              f'EFT1 RHS={wilson["metadata"]["rhs_evaluations"]}', flush=True)

    args.output.mkdir(parents=True, exist_ok=True)
    make_figures(args.output, rows, traj, ms, args.mu_low)
    with (args.output / 'D_m2_frozen_running_observables.csv').open(
            'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['scenario','EFT1_renormalisable_running','SMEFT_running',
                    'm1_eV','m2_eV','m3_eV','delta_m21_sq_eV2','delta_m31_sq_eV2'])
        for label, item in rows.items():
            w.writerow([label, label[1]=='1', label[-1]=='1',
                        *(float(v) for v in item['mass']), item['dm21'], item['dm31']])
    report = dict(
        model='D_-2', MF_GeV=mf, MS_GeV=ms, mu_low_GeV=args.mu_low,
        comparison='EFT1-renormalisable on/off x SMEFT-Weinberg on/off',
        UV_fermion_boundary_held_fixed=True,
    )
    # Named explicit metadata to avoid conflating a frozen EFT1 scalar background
    # with freezing the entire Wilson evolution.
    report['E0']='MF renormalisable EFT1 couplings frozen, Wilson ADM active; MF quartic at scalar threshold'
    report['E1']='Saved EFT1 renormalisable trajectory; MS quartic at scalar threshold'
    report['S0']='No SMEFT running; matched C5(MS) carried unchanged to mu_low'
    report['S1']='Full-flavour SMEFT Weinberg and SM couplings evolved to mu_low'
    report['weinberg_and_physics_certified']=False
    report['scenarios']={k:{'masses_eV':v['mass'].tolist(),
                            'delta_m21_sq_eV2':v['dm21'],
                            'delta_m31_sq_eV2':v['dm31']}
                         for k,v in rows.items()}
    (args.output/'D_m2_frozen_running_summary.json').write_text(
        json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('\nScenario masses and splittings (all eV / eV^2):')
    for name,row in rows.items():
        print(f'{name:6s} masses={np.array2string(row["mass"],precision=9)} '
              f'dm21={row["dm21"]:.12e} dm31={row["dm31"]:.12e}')
    for name in ['D_m2_frozen_running_masses.png',
                 'D_m2_frozen_running_splittings.png',
                 'D_m2_smeft_m3_trajectory.png',
                 'D_m2_frozen_running_heatmap_percent.png',
                 'D_m2_frozen_running_heatmap.csv',
                 'D_m2_frozen_running_observables.csv',
                 'D_m2_frozen_running_summary.json']:
        print('Saved:',args.output/name)
    print('No benchmark or production code changed; generated figures are numerical candidate.')

if __name__=='__main__':
    main()
