#!/usr/bin/env python3
"""D(-2) rough complete numerical candidate: corrected EFT1 -> both scalars out -> SMEFT.

Read-only existing benchmark; unchanged fitted Yukawas. Uses current repository
Wilson evolution, simultaneous eight-real-scalar threshold, final full-flavour
SM+Weinberg RG solver, and physical mass conversion convention.

Output is numerical-candidate only: formal GammaCC-to-Weyl sign, native singlet
phase mapping, full one-loop ADM and full finite matching not certified.
No old direct scalar LLSS leading-log double counting.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from RunDMinus2CorrectedEFT1 import prepare
from Numerical.parallel_t3.wilson.workflow import run_wilson
from Numerical.parallel_t3.threshold.matching import (
    match_scalar_threshold, scalar_gauge_ward_residuals,
)
from Numerical.running.SMWeinbergEvolution import (
    SMWeinbergInitialConditions, evolve_sm_weinberg,
)
from physics.NeutrinoMass import numerical_neutrino_mass_matrix

DEFAULT_ROOT = Path('output/parallel_t3/rough_inputs/fit16_cc2dcbec4a')
DEFAULT_OUT = Path('output/parallel_t3/D_m2_corrected_rough_final')


def native_sm_at_ms(input_root: Path, ms: float) -> dict:
    """Use the retained same-run EFT1 native beta state, never invented SM inputs."""
    path = input_root / 'a' / 'D_m2' / 'saved_native_trajectory.npz'
    with np.load(path, allow_pickle=False) as f:
        fields = ('gY', 'g2', 'g3', 'lambdaH', 'ye', 'yu', 'yd')
        keys = ('mu_eft1_gev',) + tuple('eft1__' + x for x in fields)
        missing = [x for x in keys if x not in f.files]
        if missing:
            raise ValueError(f'Retained native SMEFT boundary missing: {missing}; found {f.files}')
        mu = np.asarray(f['mu_eft1_gev'], float)
        idx = int(np.argmin(mu))
        if not np.isclose(mu[idx], ms, rtol=1e-10, atol=1e-7):
            raise ValueError(f'Native state low scale {mu[idx]} != scalar threshold {ms}')
        state = {x: np.asarray(f['eft1__' + x][idx]).copy() for x in fields}
    for name in ('gY', 'g2', 'g3', 'lambdaH'):
        value = complex(state[name].item())
        if abs(value.imag) > 1e-10 * max(1.0, abs(value.real)):
            raise ValueError(f'{name}(MS) unexpectedly complex: {value}')
        state[name] = float(value.real)
    for name in ('ye', 'yu', 'yd'):
        state[name] = np.asarray(state[name], complex)
        if state[name].shape != (3, 3):
            raise ValueError(f'Unexpected {name} shape {state[name].shape}')
    return state


def encode(x):
    x = np.asarray(x, complex)
    return {'real': x.real.tolist(), 'imag': x.imag.tolist()}


def spectrum(c5, vev):
    m = numerical_neutrino_mass_matrix(c5, vev_gev=vev)
    mass = np.sort(np.linalg.svd(m, compute_uv=False)) * 1e9
    return {'matrix_GeV': m, 'masses_eV': mass,
            'dm21_eV2': float(mass[1]**2 - mass[0]**2),
            'dm31_eV2': float(mass[2]**2 - mass[0]**2)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, default=DEFAULT_ROOT)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUT,
                        help='New output directory for independent rough results')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--mu-low', type=float, default=173.0,
                        help='Final SMEFT scale in GeV, default 173; not an EW threshold prescription')
    parser.add_argument('--vev', type=float, default=246.22)
    parser.add_argument('--solver', choices=('DOP853', 'RK2'), default='DOP853')
    parser.add_argument('--steps', type=int, default=4)
    parser.add_argument('--rtol', type=float, default=1e-6)
    parser.add_argument('--atol', type=float, default=1e-21)
    parser.add_argument('--variant', choices=('fonseca_eq22_six_term', 'ward_covariant_candidate'),
                        default='fonseca_eq22_six_term')
    parser.add_argument('--no-save', action='store_true')
    args = parser.parse_args()
    if args.steps < 1 or args.mu_low <= 0 or args.vev <= 0:
        parser.error('Steps, final scale and vev must be positive')
    if not args.no_save and args.output.exists() and any(args.output.iterdir()) and not args.overwrite:
        parser.error(f'Output directory is nonempty: {args.output}; choose new --output or --overwrite')

    basis, background, C0, scales, q_old, q_new = prepare(args.input_root / 'a' / 'D_m2')
    mf, ms = float(scales[-1]), float(scales[0])
    if not (0 < args.mu_low < ms < mf):
        parser.error(f'Require 0 < mu_low < MS < MF; got {args.mu_low}, {ms}, {mf}')
    c_path = args.input_root / 'c' / 'D_m2' / 'background.npz'
    with np.load(c_path, allow_pickle=False) as f:
        needed = ('mass_squared_gev2_ms', 'R_wilson_from_background', 'quartic_real_full_ms')
        if any(k not in f.files for k in needed):
            raise ValueError('Missing C-stage scalar mass/rotation/quartic in ' + str(c_path))
        mass2 = np.array(f['mass_squared_gev2_ms'], float)
        R = np.array(f['R_wilson_from_background'], float)
        q_saved = np.array(f['quartic_real_full_ms'], float)
    n = basis.ns-4
    if mass2.shape != (n,n) or R.shape != (n,n) or not np.allclose(R, np.eye(n), atol=1e-10):
        raise ValueError('Saved C scalar background needs an independently checked nonidentity rotation')
    q_drift = np.linalg.norm(q_old[0]-q_saved)/max(np.linalg.norm(q_saved),1e-300)
    if q_drift>1e-9:
        raise ValueError(f'A/C background quartic mismatch: {q_drift}')
    generators = np.asarray([s for group in ('SU2','U1')
                             for _,s in basis.generators()[group]])
    ward = scalar_gauge_ward_residuals(q_new[0], mass2, generators)
    if not ward['passed']:
        raise ValueError(f'Corrected quartic/mass Ward failed: {ward}')
    print(f'D(-2) rough numerical sequence: MF={mf:g}, MS={ms:g}, '
          f'SMEFT endpoint={args.mu_low:g} GeV', flush=True)
    print('Running corrected EFT1 full Wilson tensor...', flush=True)
    evolved = run_wilson(basis, C0, background, mf, ms,
                         adm_variant=args.variant, solver=args.solver,steps=args.steps,
                         rtol=args.rtol,atol=args.atol)
    print('EFT1 ||C5(MS)||:',f'{np.linalg.norm(evolved["C5_final"]):.12e}',flush=True)
    print('Matching all eight scalar real states...',flush=True)
    matched = match_scalar_threshold(
        full_wilson_ms=evolved['C_final'],
        lepton_indices=np.array(basis.leptons,dtype=int),
        higgs_indices=np.array(basis.higgs,dtype=int),
        heavy_indices=np.arange(4,basis.ns,dtype=int),
        full_quartic_ms=q_new[0],
        scalar_mass2_background=mass2,
        R_wilson_from_background=R,
        mu_s=ms,mode='rge_improved',
        separately_supplied_c5_ms=evolved['C5_final'],
        projection_tolerance=1e-8,enforce_pure_weinberg=True)
    if matched['scalar_components_removed']!=8 or matched['bubble']['pair_count']!=36:
        raise RuntimeError('Not all scalars removed together')
    Kms=np.asarray(matched['c5_smeft_at_ms'],complex)
    print('SMEFT ||C5(MS)||:',f'{np.linalg.norm(Kms):.12e}',flush=True)

    sm = native_sm_at_ms(args.input_root, ms)
    initial = SMWeinbergInitialConditions(
        gY=sm['gY'],g2=sm['g2'],g3=sm['g3'],lambdaH=sm['lambdaH'],
        ye=sm['ye'],yu=sm['yu'],yd=sm['yd'],K=Kms)
    print('Running full-flavour SM+Weinberg RG with saved MS SM parameters...',flush=True)
    low = evolve_sm_weinberg(initial,ms,args.mu_low,rtol=1e-8,atol=1e-20)
    at_ms=spectrum(Kms,args.vev)
    at_low=spectrum(low.K,args.vev)
    print('SMEFT RGE evaluations:',low.nfev)
    print('SMEFT ||C5(mu_low)||:',f'{np.linalg.norm(low.K):.12e}')
    print('At MS masses [eV]:',np.array2string(at_ms['masses_eV'],precision=10))
    print('At MS Delta m21^2, Delta m31^2 [eV^2]:',at_ms['dm21_eV2'],at_ms['dm31_eV2'])
    print('Low scale masses [eV]:',np.array2string(at_low['masses_eV'],precision=10))
    print('Low scale Delta m21^2, Delta m31^2 [eV^2]:',at_low['dm21_eV2'],at_low['dm31_eV2'])
    print('Low-scale C5 [GeV^-1]:\n',np.array2string(low.K, precision=9))
    if not args.no_save:
        args.output.mkdir(parents=True,exist_ok=True)
        npz = args.output / 'D_m2_corrected_rough_final.npz'
        summary = args.output / 'D_m2_corrected_rough_final.json'
        if not args.overwrite and (npz.exists() or summary.exists()):
            raise FileExistsError('Output exists; use a fresh --output')
        np.savez_compressed(npz,
            C5_EFT1_MS=evolved['C5_final'],
            C5_scalar_bubble_MS=matched['c5_scalar_bubble'],
            C5_SMEFT_MS=Kms,C5_SMEFT_low=low.K,
            mnu_MS_GeV=at_ms['matrix_GeV'],mnu_low_GeV=at_low['matrix_GeV'],
            masses_MS_eV=at_ms['masses_eV'],masses_low_eV=at_low['masses_eV'])
        record={
            'model':'D_-2','status':'rough_numerical_candidate_not_physics_certified',
            'input_root':str(args.input_root),'mu_f_gev':mf,'mu_s_gev':ms,
            'mu_low_gev':args.mu_low,'vev_gev':args.vev,
            'solver':args.solver,'eft1_nfev':evolved['metadata']['rhs_evaluations'],
            'smeft_nfev':low.nfev,'adm_variant':args.variant,
            'canonical_mapping':'triplet_O_only_no_lepton_epsilon',
            'quartic_mapping':'spin_one_J_portal_at_all_EFT1_scales',
            'C5_EFT1_MS':encode(evolved['C5_final']),
            'C5_scalar_bubble_MS':encode(matched['c5_scalar_bubble']),
            'C5_SMEFT_MS':encode(Kms),'C5_SMEFT_low':encode(low.K),
            'masses_MS_eV':at_ms['masses_eV'].tolist(),
            'delta_m21_sq_MS_eV2':at_ms['dm21_eV2'],
            'delta_m31_sq_MS_eV2':at_ms['dm31_eV2'],
            'masses_low_eV':at_low['masses_eV'].tolist(),
            'delta_m21_sq_low_eV2':at_low['dm21_eV2'],
            'delta_m31_sq_low_eV2':at_low['dm31_eV2'],
            'scalar_components_integrated_out':8,'scalar_eigenpairs':36,
            'weinberg_projection_residual_EFT1':matched['rge_input_C5_residual'],
            'weinberg_projection_residual_bubble':matched['scalar_bubble_projection_residual'],
            'scalar_mass_quartic_ward':ward,
            'known_uncertified':['signed GammaCC/Weyl map','complete one-loop ADM',
                                 'singlet scalar phase/source map',
                                 'complete finite scalar matching',
                                 'EW threshold and vev running prescriptions'],
        }
        summary.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        print('Saved:',npz,'and',summary)
    else:
        print('No files written.')
    print('No Yukawa optimisation or production-code changes. Numerical candidate, not certified physics.')


if __name__=='__main__':
    main()
