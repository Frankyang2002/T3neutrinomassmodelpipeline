#!/usr/bin/env python3
"""Experimental Workstream A: verify six missing neutral-scalar Yukawas in live RGBeta.

Production RGBeta's T3RGBetaBuild has y1,y2 only. This runner creates one
extra interaction *inside an isolated Wolfram process* for each of the six
models and verifies the native $yukawas record and AddYukawa Ward check.

No production source, RGE beta, Matchete matching, or Github state is changed.
The common Matchete<->RGBeta coupling/component normalization is NOT proven.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / 'Lagrangian' / 'RunT3RGBetaExtraYukawaBridge.wl'
OUTPUT = ROOT / 'Reports' / 'output' / 'workstream_a' / 'rgbeta_extra_yukawa_bridge'
SOURCE = ROOT / 'RGE' / 'running' / 'rgbeta' / 'T3RGBetaModel.wl'
CLOSURE = ROOT / 'Reports' / 'output' / 'workstream_a' / 'uv_closure_stage' / 'native_z2'

# Each tuple is (dS1, dS2, dF, alpha, new coupling, scalar field, other F Weyl field).
MODELS = {
    'A_m2': (1, 3, 2, -2, 'y2Conj', 'S2dag', 'Fc'),
    'A_p0': (1, 3, 2, 0, 'y1Conj', 'S1', 'F'),
    'D_m2': (3, 1, 2, -2, 'y2Conj', 'S2dag', 'Fc'),
    'D_p0': (3, 1, 2, 0, 'y1Conj', 'S1', 'F'),
    'E_m2': (3, 3, 2, -2, 'y2Conj', 'S2dag', 'Fc'),
    'E_p0': (3, 3, 2, 0, 'y1Conj', 'S1', 'F'),
}
NATIVE_SUCCESS = 'ExtraYukawaRegisteredNativeWardCheckPassed'
SUCCESS = 'SixExtraRGBetaYukawasNativeRegisteredPhysicalNormalizationUncertified'


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError(msg)


def expected_fields(model: str) -> tuple[str, str]:
    alpha = MODELS[model][3]
    if alpha == 0:
        return '{S1, l, F}', '{Bar[S1], l, F}'
    return '{Bar[S2], l, Fc}', '{S2, l, Fc}'


def compact(text: str) -> str:
    return re.sub(r'\s+', '', text)


def model_physics_test() -> None:
    require(len(MODELS) == 6 and len({tuple(row[:4]) for row in MODELS.values()}) == 6,
            'Six independent model points required')
    for model, (d1, d2, df, alpha, coupling, scalar, f) in MODELS.items():
        require(df == 2 and alpha in (-2, 0), model + ' must have Dirac doublet F')
        ys1, ys2, yf = alpha / 2, (alpha + 2) / 2, (alpha + 1) / 2
        if alpha == 0:
            require(ys1 == 0 and coupling == 'y1Conj' and scalar == 'S1' and f == 'F',
                    model + ': S1 extra assignment')
            require(abs(df - d1) == 1 and (-1 / 2 + yf + ys1) == 0,
                    model + ': SU2 or U1 hypercharge')
        else:
            require(ys2 == 0 and coupling == 'y2Conj' and scalar == 'S2dag' and f == 'Fc',
                    model + ': S2 extra assignment')
            require(abs(df - d2) == 1 and (-1 / 2 - yf - ys2) == 0,
                    model + ': SU2 or U1 hypercharge')
        # The Yukawa fermion pair 2 x 2 contains a scalar d=1 or d=3 exactly once.
        require((d1 if alpha == 0 else d2) in (1, 3), model + ': SU2 multiplicity')


def source_guard(root: Path = ROOT) -> None:
    p = root / 'RGE' / 'running' / 'rgbeta' / 'T3RGBetaModel.wl'
    require(p.is_file(), 'Missing current RGBeta model source: ' + str(p))
    code = p.read_text(encoding='utf-8-sig')
    # These are factual source anchors, not a proof that the whole RGBeta model
    # contains exactly two Yukawa structures outside the ordinary T3 builder.
    require('T3RGBetaAddT3Yukawas[' in code and
            re.search(r'AddYukawa\[\s*y1\s*,', code) and
            re.search(r'AddYukawa\[\s*y2\s*,', code) and
            'secondHeavy = If[yF === 0, F, Fc]' in code,
            'Production RGBeta T3 Yukawa builder changed; inspect latest version')
    require(not re.search(r'\by1Conj\b|\by2Conj\b', code),
            'Production builder already mentions conjugate Yukawas; audit must be updated')
    require((root / 'Lagrangian' / 'RunT3RGBetaExtraYukawaBridge.wl').is_file(),
            'Experimental Wolfram driver missing')


def check_matchete(model: str, root: Path = ROOT) -> str:
    p = root / 'Reports' / 'output' / 'workstream_a' / 'uv_closure_stage' / 'native_z2' / (model + '.json')
    require(p.is_file(), model + ': missing previous live Matchete parity/Yukawa report: ' + str(p))
    obj = json.loads(p.read_text(encoding='utf-8-sig'))
    require(obj.get('model') == model and
            obj.get('status') == 'PassedActualBSMTermParityAndMatchete' and
            obj.get('extra_yukawas_included') == 1 and
            obj.get('matchete_full_check_passed') is True,
            model + ': previous Matchete extra-Yukawa evidence is not passing')
    return str(p.relative_to(root))


def check_native(model: str, obj: dict) -> dict:
    dims = MODELS[model]
    fields, base = expected_fields(model)
    require(obj.get('format') == 't3-workstream-a-rgbeta-extra-yukawa-registration-v1' and
            obj.get('model') == model and obj.get('status') == NATIVE_SUCCESS,
            model + ': RGBeta returned failure or wrong report')
    require(obj.get('dimensions_alpha') == list(dims[:4]) and
            obj.get('registered_yukawas_before') == 5 and
            obj.get('registered_yukawas_after') == 6 and
            obj.get('extra_yukawa') == dims[4] and
            obj.get('baseline_yukawa') == ('y1' if dims[3] == 0 else 'y2'),
            model + ': registry size, coupling or representation mismatch')
    require(compact(obj.get('extra_fields_inputform', '')) == compact(fields) and
            compact(obj.get('baseline_fields_inputform', '')) == compact(base),
            model + ': native Weyl/scalar field orientation changed')
    require(obj.get('same_native_su2_invariant_as_baseline') is True and
            obj.get('native_ward_check_invoked') is True and
            obj.get('native_ward_check_passed') is True and
            obj.get('independent_flavour_matrix_registered') is True and
            obj.get('hypercharge_scalar_extra') == 0 and
            obj.get('hypercharge_fermion') == dims[3] / 2 + .5,
            model + ': native Ward, invariant, flavour, or hypercharge evidence absent')
    require(obj.get('one_loop_beta_evaluated') is False and
            obj.get('physical_matchete_rgbeta_coupling_map_certified') is False and
            obj.get('physical_component_basis_certified') is False and
            obj.get('production_state_modified') is False,
            model + ': report overclaims physical certification or changes')
    invariant = obj.get('extra_invariant_inputform', '')
    other = obj.get('baseline_invariant_inputform', '')
    require(isinstance(invariant, str) and invariant == other and
            (('delS2[' in invariant) if (dims[0] if dims[3] == 0 else dims[1]) == 3
             else ('eps[' in invariant)),
            model + ': live CG invariant missing or different')
    return {
        'model': model, 'status': 'VerifiedLiveNativeRegistrationAndWardCheck',
        'extra_yukawa': dims[4], 'fields': fields,
        'native_invariant': invariant,
        'matchete_extra_yukawa_already_accepted': True,
        'cross_package_yukawa_normalization_certified': False,
        'beta_evaluated': False,
    }


def make_report(rows: list[dict], root: Path = ROOT) -> dict:
    require(len(rows) == 6 and {row['model'] for row in rows} == set(MODELS),
            'Exactly six unique models are required')
    require(all(row['status'] == 'VerifiedLiveNativeRegistrationAndWardCheck' for row in rows),
            'Some model registration did not pass')
    return {
        'format': 't3-workstream-a-rgbeta-extra-yukawa-bridge-v1',
        'status': SUCCESS,
        'verified_models': 6,
        'rgbeta_yukawa_couplings_before': 5,
        'rgbeta_yukawa_couplings_after': 6,
        'new_yukawa_matrix_types_across_models': 6,
        'full_16_model_uv_yukawa_matrix_types_from_previous_audits': 38,
        'native_rgbeta_ward_check_passed_six': True,
        'physical_yukawa_coupling_map_certified': False,
        'global_triplet_component_basis_certified': False,
        'one_loop_betas_recomputed': False,
        'production_model_changed': False,
        'overall_uv_physical_completeness_certified': False,
        'remaining': [
            'Independent Matchete-RGBeta Yukawa tensor orientation and normalization',
            'Actual common triplet component basis for general scalar vertices',
            'Integration into RGBeta beta functions remains opt-in/unperformed',
        ],
        'models': sorted(rows, key=lambda x: x['model']),
    }


def markdown(report: dict) -> str:
    lines = [
        '# Workstream A — RGBeta generalized neutral-scalar Yukawa bridge', '',
        f'**Status:** `{report["status"]}`', '',
        'The live, isolated RGBeta Yukawa registry and native Ward check passed for the six',
        'neutral-scalar models. Matchete previously accepted the corresponding extra Yukawa',
        'interaction. This does **not** establish cross-package coupling normalization,',
        'recompute beta functions, or modify the production RGBeta model.', '',
        '| Model | New RGBeta Yukawa | Weyl field list | Native Ward check |',
        '|---|---|---|---|',
    ]
    for row in report['models']:
        lines.append(f'| {row["model"]} | `{row["extra_yukawa"]}` | '
                     f'`{row["fields"]}` | Passed |')
    lines += ['', '## Remaining limitations', '']
    lines += [f'- {item}' for item in report['remaining']]
    lines += ['', 'The 38 coupling-matrix types are an aggregate count across the 16 models,',
              'not 38 separate couplings in one model. Each special model has three,',
              'not two, independent BSM Yukawa matrix types.', '']
    return '\n'.join(lines)


def selftest() -> None:
    model_physics_test()
    # Self-test the source guard against source drift and the logical gates.
    source_guard()
    a = 'A_m2'
    good = {
        'format': 't3-workstream-a-rgbeta-extra-yukawa-registration-v1',
        'model': a, 'status': NATIVE_SUCCESS, 'dimensions_alpha': [1, 3, 2, -2],
        'registered_yukawas_before': 5, 'registered_yukawas_after': 6,
        'extra_yukawa': 'y2Conj', 'baseline_yukawa': 'y2',
        'extra_fields_inputform': '{Bar[S2], l, Fc}',
        'baseline_fields_inputform': '{S2, l, Fc}',
        'same_native_su2_invariant_as_baseline': True,
        'native_ward_check_invoked': True, 'native_ward_check_passed': True,
        'independent_flavour_matrix_registered': True,
        'hypercharge_scalar_extra': 0, 'hypercharge_fermion': -.5,
        'one_loop_beta_evaluated': False,
        'physical_matchete_rgbeta_coupling_map_certified': False,
        'physical_component_basis_certified': False,
        'production_state_modified': False,
        'extra_invariant_inputform': 'delS2[SU2L, t3ry1, t3ry2, t3ry3]',
        'baseline_invariant_inputform': 'delS2[SU2L, t3ry1, t3ry2, t3ry3]',
    }
    check_native(a, good)
    for key, bad in [('native_ward_check_passed', False),
                     ('physical_component_basis_certified', True),
                     ('extra_fields_inputform', '{S2, l, Fc}'),
                     ('extra_invariant_inputform', 'eps[SU2L, 1, 2]')]:
        altered = dict(good, **{key: bad})
        try:
            check_native(a, altered)
        except ValueError:
            pass
        else:
            raise AssertionError('Missing negative-control rejection: ' + key)
    print('PASS: six-model Weyl/SU2/U1 census, current-source guard, native-report gates')
    print('PASS: deliberate wrong chirality/scalar-orientation, Ward and certification controls')
    print('NOTE: no live Wolfram kernel executed by --self-test')


def run(args: argparse.Namespace) -> int:
    model_physics_test()
    source_guard()
    rows = []
    out = args.output_dir.resolve()
    target_paths = []
    for model in MODELS:
        path = out / 'native_rgbeta' / (model + '.json')
        target_paths.append(path)
    target_paths += [out/'rgbeta_extra_yukawa_bridge.json',
                     out/'rgbeta_extra_yukawa_bridge.md']
    if not args.overwrite:
        existing = [str(p) for p in target_paths if p.exists()]
        require(not existing, 'Refusing to overwrite previous report; use --overwrite after backup: '\
                + ', '.join(existing))
    (out / 'native_rgbeta').mkdir(parents=True, exist_ok=True)
    for model in MODELS:
        matchete_path = check_matchete(model)
        file = out / 'native_rgbeta' / (model + '.json')
        log = file.with_suffix('.log')
        c = model.split('_')
        cmd = [args.wolframscript, '-file', str(DRIVER), str(file), *c]
        print(f'RGBeta native extra Yukawa {model}: {MODELS[model][4]}', flush=True)
        env = os.environ.copy()
        try:
            cp = subprocess.run(cmd, cwd=str(ROOT), env=env, capture_output=True,
                                text=True, encoding='utf-8', errors='replace',
                                timeout=args.timeout_seconds, check=False)
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise RuntimeError(model + ': Wolfram could not execute: ' + str(exc)) from exc
        log.write_text('command: ' + ' '.join(cmd) + '\nexit code: ' + str(cp.returncode) + '\n'
                       + cp.stdout + '\n' + cp.stderr, encoding='utf-8')
        require(file.is_file(), model + ': no native RGBeta JSON. Inspect ' + str(log))
        obj = json.loads(file.read_text(encoding='utf-8-sig'))
        require(cp.returncode == 0, model + ': native RGBeta returned code '
                + str(cp.returncode) + ': ' + str(obj.get('reason', '')) + '. See ' + str(log))
        row = check_native(model, obj)
        row['matchete_evidence_path'] = matchete_path
        rows.append(row)
        print(f'  {model}: native RGBeta registration + Ward PASS; Matchete extra Yukawa PASS',
              flush=True)
    report = make_report(rows)
    json_path = out/'rgbeta_extra_yukawa_bridge.json'
    md_path = out/'rgbeta_extra_yukawa_bridge.md'
    json_path.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    md_path.write_text(markdown(report), encoding='utf-8')
    print('RESULT: ' + report['status'])
    print('JSON: ' + str(json_path))
    print('REPORT: ' + str(md_path))
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group(required=True)
    modes.add_argument('--self-test', action='store_true')
    modes.add_argument('--all', action='store_true')
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('--wolframscript', default='wolframscript')
    p.add_argument('--timeout-seconds', type=int, default=900)
    p.add_argument('--output-dir', type=Path, default=OUTPUT)
    a = p.parse_args()
    try:
        if a.self_test:
            selftest()
            return 0
        return run(a)
    except (ValueError, RuntimeError, OSError, json.JSONDecodeError) as exc:
        print('FAIL: ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
