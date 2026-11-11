"""Experimental A(-4)/D(-2) coupled RG with cached non-HHSS quartics.

Reuses the existing 14-nonmixed simultaneous RGBeta + seven-channel master
DOP853 IVP, UV/F matching, source HHSS, both-scalar finite bubble, SMEFT,
Takagi and stage-freeze machinery. ONLY the provisional non-HHSS quartic
context is replaced by a Batch3B cached tensor. This is not fully matched EFT.

Run from repository root:
 python -X utf8 -m scripts.RunT3CoupledScalarA4D2Experimental --self-test
 python -X utf8 -m scripts.RunT3CoupledScalarA4D2Experimental --response-root OUTPUT_B3B --source-run thesis16_merged_v1 --run-id scalar_cached_A4_D2_v1 --models A_-4,D_-2 --check-only
 python -X utf8 -m scripts.RunT3CoupledScalarA4D2Experimental --response-root OUTPUT_B3B --source-run thesis16_merged_v1 --run-id scalar_cached_A4_D2_v1 --models A_-4,D_-2

Writes only fresh output/studies/general_scalar_extended/
experimental_coupled_cached_scalar_a4d2/<run-id>/ . No GitHub writes.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

from Numerical.running.T3CoupledQuarticContext import (
    CandidateQuarticContext, self_test as context_self_test,
)

MODELS = ('A_-4', 'D_-2')
SCHEMA = 't3-experimental-coupled-rgbeta-cached-nonhhss-a4d2-v1'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def _numbers(record, *keys):
    for key in keys:
        if key in record and record[key] is not None:
            value = record[key]
            if isinstance(value, (int, float)) and math.isfinite(float(value)):
                return float(value)
    return None


def _decode_c5(data):
    if not isinstance(data, dict) or 'real' not in data or 'imag' not in data:
        return None
    re_part = np.asarray(data['real'], float)
    im_part = np.asarray(data['imag'], float)
    if re_part.shape != (3,3) or im_part.shape != (3,3):
        return None
    return re_part + 1j*im_part


def _compare_with_saved(project, model, new_cases):
    parents = {
        'A_-4': Path('output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/simultaneous_A_m4_v2/A_m4/cases.json'),
        'D_-2': Path('output/studies/general_scalar_extended/experimental_simultaneous_eft1_general_nonmixed/simultaneous_nonmixed13_v1/D_m2/cases.json'),
    }
    src = Path(project) / parents[model]
    if not src.is_file():
        return {'status':'ReferenceCasesMissing', 'path':str(src)}
    old = json.loads(src.read_text(encoding='utf-8-sig'))
    if old.get('status') != 'Computed' or not isinstance(old.get('cases'), dict):
        return {'status':'ReferenceCasesNotComputed', 'path':str(src)}
    out = {}
    for key in ('1111','0111','1001','1110'):
        if key not in old['cases'] or key not in new_cases:
            return {'status': 'RequiredReferenceCaseMissing', 'case':key}
        a = old['cases'][key]
        b = new_cases[key]
        entry = {}
        for name in ('c5_at_scalar_threshold', 'c5_simultaneous_running', 'c5_joint_scalar_bubble'):
            old_c = _decode_c5(a.get(name))
            new_c = _decode_c5(b.get(name))
            if old_c is not None and new_c is not None:
                entry[name+'_relative_change'] = float(np.linalg.norm(new_c-old_c) /
                                               max(np.linalg.norm(old_c), 1e-100))
        for name in ('chi2', 'chi_square'):
            if name in a and name in b:
                try:
                    entry[name+'_old'] = float(a[name]); entry[name+'_new'] = float(b[name])
                    entry[name+'_difference'] = float(b[name]-a[name])
                except (TypeError, ValueError):
                    pass
        out[key] = entry
    return {'status':'ComparedConditionalStageFreezeCases', 'path':str(src), 'cases':out,
            'interpretation':'Experimental non-HHSS model difference; no certified physical shift'}



def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--response-root', type=Path)
    p.add_argument('--source-run', default='thesis16_merged_v1')
    p.add_argument('--run-id', default='scalar_cached_A4_D2_v1')
    p.add_argument('--models', default='A_-4,D_-2')
    p.add_argument('--check-only', action='store_true')
    p.add_argument('--spinor-factor', type=float, default=1.0)
    p.add_argument('--vev', type=float, default=246.22)
    p.add_argument('--wolframscript', default='wolframscript')
    args = p.parse_args(argv)
    if args.self_test:
        context_self_test()
        print('PASS: experimental runner CLI and candidate context smoke test')
        return 0
    project = args.project.resolve()
    require(project.is_dir(), 'Project root does not exist')
    require(args.response_root is not None, 'Provide Batch3B output path with --response-root')
    selected = tuple(m.strip() for m in args.models.split(',') if m.strip())
    require(selected and len(set(selected)) == len(selected) and set(selected) <= set(MODELS),
            'Only A_-4 and D_-2 supported; no mixed-sector source guessing')
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', args.run_id),
            'Invalid fresh run id')
    require(math.isfinite(args.spinor_factor) and math.isfinite(args.vev) and args.vev > 0,
            'Invalid spinor-factor or Higgs vev')

    from scripts import RunT3SimultaneousEFT1All16 as simultaneous
    from scripts import RunT3GeneralSequentialFitAblation as original

    cached = {model: CandidateQuarticContext.load(project, args.response_root.resolve(), model)
              for model in selected}
    # The context replacement is process-local and affects only the new runner.
    def replacement_quartic(model, state, hhss):
        if model not in cached:
            raise ValueError('Unexpected model in candidate scalar context: '+model)
        return cached[model].tensor(state, hhss)

    simultaneous.general_quartic = replacement_quartic
    original._nonmixed = simultaneous._nonmixed_general
    original.OUT = project / 'output/studies/general_scalar_extended/experimental_coupled_cached_scalar_a4d2'
    def reject_mixed(model, *a, **kw):
        raise ValueError(f'{model}: experimental scalar response not prepared for mixed B/C')
    original._mixed = reject_mixed
    original_save = original.save
    def save_with_metadata(path, data):
        if isinstance(data, dict):
            data = dict(data)
            if data.get('format') == 't3-full-general-16_U_E_W_S_ablation-v1':
                data['format'] = SCHEMA
                data['running_background'] = 'coupled_RGBeta_7term_with_Batch3B_nonHHSS'
                data['one_scale_C5_used_as_W_on_reference'] = False
                data['complete_physical_C5'] = False
                data['hard_matching_residual_independently_calculated'] = False
            if data.get('status') == 'Computed' and 'cases' in data:
                data['matching_approximation'] = SCHEMA
                data['full_sequential_general_C5_certified'] = False
                data['nonHHSS_quartic_replaced_by_candidate_cached_full_tensor'] = True
                data['source_HHSS_unchanged'] = True
        original_save(path, data)
    original.save = save_with_metadata
    callargs = [
        '--source-run', args.source_run, '--run-id', args.run_id,
        '--models', ','.join(selected), '--spinor-factor', str(args.spinor_factor),
        '--vev', str(args.vev), '--wolframscript', args.wolframscript,
        '--release-dir', str(project/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'),
        '--v4-channel-dir', str(project/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'),
        '--quartic-seed-dir', str(project/'output/full16_symbolic'),
        '--registry', str(project/'RGE/general/generated/general_rgbeta_quartic_registry.json'),
        '--reference-running-dir', str(project/'output/general_observable/all16_v1'),
    ]
    if args.check_only:
        callargs.append('--check-only')
    dest = original.OUT/args.run_id
    if not args.check_only:
        require(not dest.exists(), f'Refusing to overwrite prior run: {dest}')
    result = original.main(callargs)
    if args.check_only:
        return result
    comparisons = {}
    for model in selected:
        tag = original.tag(model)
        path = dest/tag/'cases.json'
        if path.is_file():
            record = json.loads(path.read_text(encoding='utf-8-sig'))
            if record.get('status') == 'Computed':
                comparisons[model] = _compare_with_saved(project, model, record['cases'])
            else:
                comparisons[model] = {'status':'Blocked'}
        else:
            comparisons[model] = {'status':'NoCasesWritten'}
    diagnostics = {
        'format': SCHEMA,
        'status': 'ConditionalExperimentCompleted' if result == 0 else 'ConditionalExperimentPartiallyBlocked',
        'source_HHSS_preserved': True,
        'cached_nonHHSS_used': True,
        'native_to_RGBeta_physical_normalisation_certified': False,
        'fermion_one_loop_finite_matching_certified': False,
        'complete_scalar_threshold_certified': False,
        'scalar_response': {m: c.diagnostics() for m, c in cached.items()},
        'historical_comparison': comparisons,
        'saved_result_path': str(dest),
    }
    (dest/'quartic_context_diagnostics.json').write_text(
        json.dumps(diagnostics, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print('Experimental quartic diagnostics:', dest/'quartic_context_diagnostics.json')
    for model in selected:
        print(model, 'quartic context calls:', cached[model].number_of_calls,
              'HHSS candidate/source discrepancy max:',
              cached[model].max_source_hhss_relative, flush=True)
    return result


if __name__ == '__main__':
    raise SystemExit(main())
