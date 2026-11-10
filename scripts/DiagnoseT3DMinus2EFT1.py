#!/usr/bin/env python3
"""Two fixed-Yukawa D(-2) EFT1 diagnostics; no fitting or physics changes.

Run from the repository root. Writes new scratch plus a JSON report under --output.
This intentionally reuses the exploratory source and stage-freeze patches.
"""
from __future__ import annotations

import argparse
import json
import sys
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.complete.all16_rough_yukawa_fit import (
    BENCHMARK, FreshAll16Evaluator, read_json, save_json,
)
from Numerical.parallel_t3.rough_inputs import builder as builder_module
from Numerical.parallel_t3.threshold import preflight as preflight_module
from Numerical.parallel_t3.threshold import rough_numerics as numerics_module
from Numerical.parallel_t3.wilson import workflow, rough_stream2
from Numerical.running import IntermediateScalarRunner
from scripts.BuildT3GeneralisedThesisTables import load_reference
from scripts.RunT3All16ExploratoryFit import (
    corrected_quartic_unit_tensors, consistent_project_c5,
    exploratory_preflight, exploratory_rough_case,
)
from scripts.RunT3GeneralisedRunningComparisons import (
    _frozen_eft1, _frozen_wilson_beta,
)


def load_npz(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f'Missing retained tensor: {path}')
    with np.load(path, allow_pickle=False) as d:
        return {k: d[k].copy() for k in d.files}


def norm(a):
    return float(np.linalg.norm(a))


def complex_pair(z):
    return {'real': float(np.real(z)), 'imag': float(np.imag(z))}


def pack_matrix(a):
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def analyse(root: Path, token: str, model: str):
    # The background.contracts.tag() mapping is authoritative, not the spelling of model.
    from Numerical.parallel_t3.background.contracts import tag
    slug = tag(model)
    result = load_npz(root / 'output/parallel_t3/rough_results' / token / slug / 'C5_smeft_at_MS.npz')
    low = load_npz(root / 'output/parallel_t3/rough_results' / token / slug / 'C5_smeft_low.npz')
    initial_path = root / 'output/parallel_t3/rough_inputs' / token / 'a' / slug
    case = rough_stream2.load_canonical(root, initial_path, model,
                                        c5_boundary='require', assume_formal_identity=True)
    f_candidates = list(initial_path.glob('*.npz'))
    a = np.asarray(result['C5_EFT1_MS'], complex)
    b = np.asarray(result['C5_scalar_bubble_MS'], complex)
    h = np.asarray(result['C5_independent_finite_MS'], complex)
    c = np.asarray(result['C5_SMEFT_MS'], complex)
    assert all(x.shape == (3, 3) for x in (a, b, h, c))
    closure = norm(a + b + h - c) / max(norm(c), 1e-300)
    interference = 2.0 * np.real(np.vdot(a, b))
    return {
        'scratch_token': token,
        'source_boundary_npz_paths': [str(p.relative_to(root)) for p in f_candidates],
        'C_full_MF': case.initial,
        'background_Y_MF': case.background.at(np.log(case.mf))[0],
        'background_lambda_MF': case.background.at(np.log(case.mf))[1],
        'background_gauge_MF': np.asarray([g for g, _ in case.background.at(np.log(case.mf))[2]]),
        'norm_C5_EFT1_MS_GeV_inverse': norm(a),
        'norm_C5_scalar_bubble_MS_GeV_inverse': norm(b),
        'norm_C5_independent_finite_MS_GeV_inverse': norm(h),
        'norm_C5_SMEFT_MS_GeV_inverse': norm(c),
        'two_Re_Tr_C5_EFT1_dagger_C5_bubble_GeV_inverse_squared': float(interference),
        'cross_EFT1_hard': float(2*np.real(np.vdot(a, h))),
        'cross_bubble_hard': float(2*np.real(np.vdot(b, h))),
        'relative_closure_residual': closure,
        'C5_EFT1_MS': pack_matrix(a), 'C5_scalar_bubble_MS': pack_matrix(b),
        'C5_independent_finite_MS': pack_matrix(h), 'C5_SMEFT_MS': pack_matrix(c),
        'C5_SMEFT_low': pack_matrix(np.asarray(low['C5_SMEFT_low'], complex)),
        'norm_C5_SMEFT_low_GeV_inverse': norm(low['C5_SMEFT_low']),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tables', type=Path, required=True)
    p.add_argument('--reference-results', nargs='+', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--batch1', default='output/batch1_general_contract/batch1_native_v1')
    args = p.parse_args()
    root = ROOT.resolve()
    out = (root / args.output).resolve()
    if out.exists() or (root/'output/parallel_t3').resolve() not in out.parents:
        p.error('--output must be a NEW directory under output/parallel_t3')
    tables = (root / args.tables).resolve()
    provenance = read_json(tables / 'provenance.json')
    baseline = read_json(tables / 'results.json')
    if baseline['D_-2']['D_-2']['status'] != 'VALID_EXPLORATORY':
        p.error('Missing valid D(-2) original full-running row')
    refdirs = [(root / x).resolve() for x in args.reference_results]
    ref = load_reference('D_-2', refdirs)
    if str(ref['source']) != provenance['reference_benchmarks']['D_-2']['result_json']:
        p.error('Reference result differs from baseline provenance')
    cfg = read_json(root / BENCHMARK)
    out.mkdir(parents=True)
    report = {'model':'D_-2', 'reference':'D_-2',
              'source':str(ref['source']), 'counting_mode':'rge_improved',
              'conditional_not_physical':True, 'cases':{}}
    for name, frozen in [('111', False), ('101', True)]:
        evaluator = FreshAll16Evaluator(
            root, 'D_-2', cfg, batch1=args.batch1,
            gamma_factor=float(provenance['gamma_factor_assumed']),
            mu_low=float(provenance['mu_low_gev']), vev=float(cfg['vev_GeV']),
            reduced_potential=True, mass_bilinears=provenance['mass_bilinears'],
            wilson_solver=provenance['wilson_solver'],
            wilson_steps=int(provenance['wilson_steps']), keep_trial_files=True)
        with ExitStack() as stack:
            stack.enter_context(patch.object(builder_module, 'quartic_unit_tensors', corrected_quartic_unit_tensors))
            stack.enter_context(patch.object(workflow, 'project_c5', consistent_project_c5))
            stack.enter_context(patch.object(rough_stream2, 'project_c5', consistent_project_c5))
            stack.enter_context(patch.object(preflight_module, 'preflight_pair', exploratory_preflight))
            stack.enter_context(patch.object(numerics_module, 'run_rough_case', exploratory_rough_case))
            if frozen:
                stack.enter_context(patch.object(IntermediateScalarRunner, 'run_intermediate_scalar_segment', _frozen_eft1))
                stack.enter_context(patch.object(workflow, 'beta', _frozen_wilson_beta))
            trial = evaluator(ref['x'].copy())
        token = trial.meta['scratch_token']
        rec = analyse(root, token, 'D_-2')
        boundary_arrays = {k:rec.pop(k) for k in ('C_full_MF','background_Y_MF',
            'background_lambda_MF','background_gauge_MF')}
        if name == '111':
            reference_boundary = boundary_arrays
        else:
            report['MF_boundary_check'] = {k:{
                'relative_difference':norm(v-reference_boundary[k])/max(norm(reference_boundary[k]),1e-300),
                'exact_array_equal':bool(np.array_equal(v,reference_boundary[k]))
            } for k,v in boundary_arrays.items()}
        rec['observables'] = trial.observables
        rec['chi2'] = trial.chi2
        report['cases'][name] = rec
        save_json(out/'diagnostic.json', report)
        print(f'{name}: C5_EFT1={rec["norm_C5_EFT1_MS_GeV_inverse"]:.8e}, '
              f'bubble={rec["norm_C5_scalar_bubble_MS_GeV_inverse"]:.8e}, '
              f'total={rec["norm_C5_SMEFT_MS_GeV_inverse"]:.8e}', flush=True)
    f = report['cases']['111']; z = report['cases']['101']
    def matrix(d, key):
        return np.asarray(d[key]['real']) + 1j*np.asarray(d[key]['imag'])
    A = matrix(f, 'C5_SMEFT_MS'); B = matrix(z, 'C5_SMEFT_MS')
    Al = matrix(f, 'C5_SMEFT_low'); Bl = matrix(z, 'C5_SMEFT_low')
    def alignment(x, y):
        k = np.vdot(y,x)/np.vdot(y,y)
        return {'k':complex_pair(k), 'abs_k':float(abs(k)),
                'relative_residual':norm(x-k*y)/max(norm(x),1e-300)}
    report['comparison'] = {
        'MS_total_norm_ratio_111_over_101':norm(A)/max(norm(B),1e-300),
        'low_total_norm_ratio_111_over_101':norm(Al)/max(norm(Bl),1e-300),
        'MS_proportional_alignment':alignment(A,B),
        'low_proportional_alignment':alignment(Al,Bl),
        'dm21_percent':100*(f['observables']['dm21_sq_eV2']/z['observables']['dm21_sq_eV2']-1),
        'dm31_percent':100*(f['observables']['dm31_sq_eV2']/z['observables']['dm31_sq_eV2']-1),
    }
    save_json(out/'diagnostic.json', report)
    print(json.dumps(report['comparison'],indent=2))
    print('Retained per-trial output is in output/parallel_t3/rough_results/<scratch_token>/')
    print('Report:',out/'diagnostic.json')

if __name__ == '__main__':
    main()
