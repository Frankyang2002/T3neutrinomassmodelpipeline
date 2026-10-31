#!/usr/bin/env python3
"""Read-only Workstream A: six isolated RGBeta extra-Yukawa registrations.

Native RGBeta v1.2.0 CheckInvariance generated unevaluated symbolic Ward
residuals for A_m2, including SU3 for colour singlet fields. The Wolfram
companion intentionally omits that option and does NOT claim Ward success.

The audit verifies exact quantum-number eligibility, the native registered
field/invariant/chirality/flavour record, and prior Matchete acceptance.
Physical Matchete-RGBeta CG normalization and RGEs remain uncertified.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / 'Lagrangian' / 'RunT3RGBetaExtraYukawaBridge.wl'
OUTPUT = ROOT / 'Reports' / 'output' / 'workstream_a' / 'rgbeta_extra_yukawa_bridge'
SOURCE = ROOT / 'RGE' / 'running' / 'rgbeta' / 'T3RGBetaModel.wl'
MODELS = {
    'A_m2': (1, 3, 2, -2, 'y2Conj'),
    'A_p0': (1, 3, 2, 0, 'y1Conj'),
    'D_m2': (3, 1, 2, -2, 'y2Conj'),
    'D_p0': (3, 1, 2, 0, 'y1Conj'),
    'E_m2': (3, 3, 2, -2, 'y2Conj'),
    'E_p0': (3, 3, 2, 0, 'y1Conj'),
}
NATIVE_SUCCESS = 'ExtraYukawaRegisteredNativeWardUncertified'
SUCCESS = 'SixExtraRGBetaYukawasRegisteredIndependentGaugeCensusPassedWardUncertified'
FORMAT = 't3-workstream-a-rgbeta-extra-yukawa-registration-v1'
WARD_STATUS = 'NotInvokedAfterUnresolvedA_m2SymbolicDiagnostic'


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def compact(txt: str) -> str:
    return re.sub(r'\s+', '', txt)


def expected_fields(model: str) -> tuple[str, str]:
    if MODELS[model][3] == 0:
        return '{S1, l, F}', '{Bar[S1], l, F}'
    return '{Bar[S2], l, Fc}', '{S2, l, Fc}'


def quantum_numbers(model: str) -> tuple[Fraction, int]:
    d1, d2, df, alpha, coupling = MODELS[model]
    require(df == 2 and alpha in (-2, 0), model + ': wrong heavy representation')
    if alpha == 0:
        require(coupling == 'y1Conj', model + ': unexpected coupling')
        d_scalar = d1
        y_scalar, y_heavy = Fraction(alpha, 2), Fraction(alpha+1, 2)
    else:
        require(coupling == 'y2Conj', model + ': unexpected coupling')
        d_scalar = d2
        y_scalar, y_heavy = -Fraction(alpha+2, 2), -Fraction(alpha+1, 2)
    # SU2 singlet appears exactly once in 2 x 2 x (1 or 3).
    require(d_scalar in (1, 3), model + ': SU2 singlet not established')
    total = Fraction(-1, 2) + y_scalar + y_heavy
    require(total == 0 and y_scalar == 0, model + ': hypercharge not conserved')
    return total, d_scalar


def source_guard(root: Path = ROOT) -> None:
    p = root / 'RGE' / 'running' / 'rgbeta' / 'T3RGBetaModel.wl'
    require(p.is_file(), 'Missing source: ' + str(p))
    s = p.read_text(encoding='utf-8-sig')
    require('T3RGBetaAddT3Yukawas[' in s and
            re.search(r'AddYukawa\[\s*y1\s*,', s) and
            re.search(r'AddYukawa\[\s*y2\s*,', s) and
            'secondHeavy = If[yF === 0, F, Fc]' in s,
            'Production RGBeta Yukawa builder changed')
    require(not re.search(r'\by1Conj\b|\by2Conj\b', s),
            'Production RGBeta builder already contains extra Yukawas')
    require((root/'Lagrangian'/'RunT3RGBetaExtraYukawaBridge.wl').is_file(),
            'Wolfram driver absent')


def check_matchete(model: str, root: Path = ROOT) -> str:
    p = root/'Reports'/'output'/'workstream_a'/'uv_closure_stage'/'native_z2'/f'{model}.json'
    require(p.is_file(), model + ': missing prior Matchete report: ' + str(p))
    r = json.loads(p.read_text(encoding='utf-8-sig'))
    require(r.get('model') == model and
            r.get('status') == 'PassedActualBSMTermParityAndMatchete' and
            r.get('extra_yukawas_included') == 1 and
            r.get('matchete_full_check_passed') is True,
            model + ': Matchete extra Yukawa was not accepted')
    return str(p.relative_to(root))


def check_native(model: str, r: dict) -> dict:
    require(isinstance(r, dict) and r.get('format') == FORMAT and
            r.get('model') == model and r.get('status') == NATIVE_SUCCESS,
            model + ': unexpected native format/status')
    dims = MODELS[model]
    total, d_scalar = quantum_numbers(model)
    new_fields, old_fields = expected_fields(model)
    require(r.get('dimensions_alpha') == list(dims[:4]) and
            r.get('extra_yukawa') == dims[4] and
            r.get('baseline_yukawa') == ('y1' if dims[3] == 0 else 'y2') and
            r.get('registered_yukawas_before') == 5 and
            r.get('registered_yukawas_after') == 6,
            model + ': unexpected native registry census')
    require(compact(r.get('extra_fields_inputform', '')) == compact(new_fields) and
            compact(r.get('baseline_fields_inputform', '')) == compact(old_fields),
            model + ': wrong scalar or Weyl orientation')
    invariant = r.get('extra_invariant_inputform', '')
    require(isinstance(invariant, str) and invariant and
            invariant == r.get('baseline_invariant_inputform') and
            (('delS2[' in invariant) if d_scalar == 3 else ('eps[' in invariant)) and
            r.get('same_native_su2_invariant_as_baseline') is True,
            model + ': native invariant not the expected one')
    require(r.get('independent_flavour_matrix_registered') is True and
            r.get('hypercharge_scalar_extra') == 0 and
            r.get('hypercharge_fermion') == Fraction(dims[3]+1,2) and
            r.get('hypercharge_sum') == str(total) and
            r.get('independent_u1_charge_selection_passed') is True and
            r.get('independent_su2_singlet_multiplicity_one') is True and
            r.get('bsm_color_singlet_by_model_source') is True,
            model + ': independent representation census or flavour missing')
    require(r.get('native_ward_check_invoked') is False and
            r.get('native_ward_check_passed') is False and
            r.get('native_ward_check_status') == WARD_STATUS and
            r.get('native_ward_check_option_disabled_explicitly') is False,
            model + ': Ward check is misreported as certified')
    for flag in ('physical_matchete_rgbeta_coupling_map_certified',
                 'physical_component_basis_certified',
                 'triplet_cross_package_basis_certified',
                 'one_loop_beta_evaluated', 'production_state_modified'):
        require(r.get(flag) is False, model + ': overclaimed ' + flag)
    return {
        'model': model,
        'status': 'NativeRegistrationAndRepresentationCensusPassedWardUncertified',
        'extra_yukawa': dims[4], 'fields': new_fields,
        'su2_scalar_dimension': d_scalar,
        'invariant': invariant,
        'native_ward_check_passed': False,
        'physical_coupling_mapping_certified': False,
    }


def make_report(rows: list[dict]) -> dict:
    require(len(rows) == 6 and {r['model'] for r in rows} == set(MODELS) and
            all(r['status'] == 'NativeRegistrationAndRepresentationCensusPassedWardUncertified'
                for r in rows), 'Six independent passing registration records required')
    return {
        'format': 't3-workstream-a-rgbeta-extra-yukawa-bridge-v2',
        'status': SUCCESS, 'verified_native_registrations': 6,
        'rgbeta_yukawa_couplings_before_per_model': 5,
        'rgbeta_yukawa_couplings_after_per_model': 6,
        'extra_yukawa_matrix_types_across_models': 6,
        'full_uv_yukawa_matrix_types_across_16_models': 38,
        'independent_u1_and_su2_representation_eligibility_passed': True,
        'matchete_extra_yukawas_previously_accepted': True,
        'native_rgbeta_ward_check_certified': False,
        'physical_matchete_rgbeta_yukawa_normalization_certified': False,
        'triplet_component_basis_certified': False,
        'production_rgbeta_extended': False,
        'beta_functions_recomputed': False,
        'remaining': [
            'Installed native CheckInvariance produces unresolved symbolic Ward residuals on A_m2',
            'SU2 eligibility is not proof of RGBeta tensor covariance in its actual component basis',
            'Cross-package absolute CG/coupling normalization remains conditional',
            'Extra Yukawas are not installed into the production RGBeta model or beta functions',
        ],
        'models': sorted(rows, key=lambda row:row['model']),
    }


def markdown(r: dict) -> str:
    out = ['# Workstream A — generalized RGBeta extra Yukawa registration', '',
           f'**Status:** `{r["status"]}`', '',
           'Six isolated registrations passed, with independent quantum-number',
           'eligibility checks and earlier Matchete-accepted operators. The installed',
           'RGBeta native Ward diagnostic is unresolved and was **not invoked** in',
           'this registration exercise. No physical coupling normalization is certified.', '',
           '| Model | New coupling | Registered fields | Native Ward |',
           '|---|---|---|---|']
    for row in r['models']:
        out.append(f'| {row["model"]} | `{row["extra_yukawa"]}` | '
                   f'`{row["fields"]}` | **Not certified** |')
    out.extend(['', '## Remaining limitations', ''])
    out.extend('- '+x for x in r['remaining'])
    out.append('')
    return '\n'.join(out)


def selftest() -> None:
    source_guard()
    require(len(MODELS) == 6, 'Six models expected')
    for model in MODELS:
        quantum_numbers(model)
    good = {
        'format': FORMAT, 'model': 'A_m2', 'status': NATIVE_SUCCESS,
        'dimensions_alpha':[1,3,2,-2], 'extra_yukawa':'y2Conj',
        'baseline_yukawa':'y2', 'registered_yukawas_before':5,
        'registered_yukawas_after':6,
        'extra_fields_inputform':'{Bar[S2],l,Fc}',
        'baseline_fields_inputform':'{S2,l,Fc}',
        'extra_invariant_inputform':'delS2[SU2L,t3ry1,t3ry2,t3ry3]',
        'baseline_invariant_inputform':'delS2[SU2L,t3ry1,t3ry2,t3ry3]',
        'same_native_su2_invariant_as_baseline':True,
        'independent_flavour_matrix_registered':True,
        'hypercharge_scalar_extra':0,'hypercharge_fermion':Fraction(-1,2),
        'hypercharge_sum':'0', 'independent_u1_charge_selection_passed':True,
        'independent_su2_singlet_multiplicity_one':True,
        'bsm_color_singlet_by_model_source':True,
        'native_ward_check_invoked':False, 'native_ward_check_passed':False,
        'native_ward_check_status':WARD_STATUS,
        'native_ward_check_option_disabled_explicitly':False,
        'physical_matchete_rgbeta_coupling_map_certified':False,
        'physical_component_basis_certified':False,
        'triplet_cross_package_basis_certified':False,
        'one_loop_beta_evaluated':False,'production_state_modified':False,
    }
    check_native('A_m2',good)
    for key, wrong in (
        ('native_ward_check_passed',True),
        ('physical_matchete_rgbeta_coupling_map_certified',True),
        ('native_ward_check_status','Passed'),
        ('extra_fields_inputform','{S2,l,Fc}'),
        ('hypercharge_sum','1'),
        ('extra_invariant_inputform','eps[SU2L,1,2]'),
    ):
        bad = {**good,key:wrong}
        try:
            check_native('A_m2',bad)
        except ValueError:
            continue
        raise AssertionError('Negative control not rejected: '+key)
    print('PASS: six exact U1/SU2 selection-rule censuses and native-registry format')
    print('PASS: negative controls reject false Ward, fields, hypercharge or normalization')
    print('NOTE: self-test does not run live RGBeta or certify a native Ward identity')


def run(args: argparse.Namespace) -> None:
    source_guard()
    rows = []
    out = args.output_dir.resolve()
    paths = [out/'native_rgbeta'/f'{m}.json' for m in MODELS]
    paths.extend([out/'rgbeta_extra_yukawa_bridge.json',out/'rgbeta_extra_yukawa_bridge.md'])
    if not args.overwrite:
        require(not any(p.exists() for p in paths),
                'Existing evidence would be overwritten: back up output directory and use --overwrite')
    (out/'native_rgbeta').mkdir(parents=True,exist_ok=True)
    for model in MODELS:
        quantum_numbers(model)
        matchete = check_matchete(model)
        target=out/'native_rgbeta'/f'{model}.json'
        log=target.with_suffix('.log')
        family,tag=model.split('_')
        cmd=[args.wolframscript,'-file',str(DRIVER),str(target),family,tag]
        print(f'RGBeta native register-only {model}: {MODELS[model][4]}',flush=True)
        try:
            process=subprocess.run(cmd,cwd=str(ROOT),env=os.environ.copy(),
                capture_output=True,text=True,encoding='utf-8',errors='replace',
                timeout=args.timeout_seconds,check=False)
        except (OSError,subprocess.TimeoutExpired) as ex:
            raise RuntimeError(f'{model}: Wolfram launch failed: {ex}') from ex
        log.write_text('command: '+' '.join(cmd)+'\nexit code: '+str(process.returncode)+'\n'
                       +process.stdout+'\n'+process.stderr,encoding='utf-8')
        require(target.is_file(),f'{model}: native JSON missing; inspect {log}')
        native=json.loads(target.read_text(encoding='utf-8-sig'))
        require(process.returncode == 0,
                f'{model}: native registration exit {process.returncode}: '
                f'{native.get("reason","")}; inspect {log}')
        row=check_native(model,native)
        row['prior_matchete_report']=matchete
        rows.append(row)
        print(f'  {model}: registry + independent U1/SU2 census PASS; native Ward UNCERTIFIED',
              flush=True)
    report=make_report(rows)
    jp=out/'rgbeta_extra_yukawa_bridge.json'
    mp=out/'rgbeta_extra_yukawa_bridge.md'
    jp.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    mp.write_text(markdown(report),encoding='utf-8')
    print('RESULT: '+report['status'])
    print('JSON: '+str(jp))
    print('REPORT: '+str(mp))


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--self-test',action='store_true')
    modes.add_argument('--all',action='store_true')
    parser.add_argument('--overwrite',action='store_true')
    parser.add_argument('--wolframscript',default='wolframscript')
    parser.add_argument('--timeout-seconds',type=int,default=900)
    parser.add_argument('--output-dir',type=Path,default=OUTPUT)
    args=parser.parse_args()
    try:
        if args.self_test:
            selftest()
        else:
            run(args)
        return 0
    except (ValueError,RuntimeError,OSError,json.JSONDecodeError) as ex:
        print('FAIL: '+str(ex),file=sys.stderr)
        return 1

if __name__=='__main__':
    raise SystemExit(main())
