#!/usr/bin/env python3
"""Run opted-in RGBeta UV plus complete generated scalar quartic projection for 16 models.

This runs the EXISTING RGE/running/rgbeta/RunGeneralT3RGBeta.wl runner and writes
separate output files; no historical beta reports are overwritten.

python scripts/RunT3ExtraYukawaStack.py --self-test
python scripts/RunT3ExtraYukawaStack.py --smoke  # five representative models
python scripts/RunT3ExtraYukawaStack.py --all    # all sixteen
python scripts/RunT3ExtraYukawaStack.py --check-only
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / 'RGE/running/rgbeta/T3RGBetaModel.wl'
RUNNER = ROOT / 'RGE/running/rgbeta/RunGeneralT3RGBeta.wl'
OUT = ROOT / 'Reports/output/workstream_a/integrated_general_yukawa_uv'

SPECS: dict[str, tuple[int, int, int, int]] = {
    'A_m4': (1, 3, 2, -4), 'A_m2': (1, 3, 2, -2), 'A_p0': (1, 3, 2, 0),
    'B_m3': (2, 2, 1, -3), 'B_m1': (2, 2, 1, -1), 'B_p1': (2, 2, 1, 1),
    'C_m3': (2, 2, 3, -3), 'C_m1': (2, 2, 3, -1), 'C_p1': (2, 2, 3, 1),
    'D_m2': (3, 1, 2, -2), 'D_p0': (3, 1, 2, 0), 'D_p2': (3, 1, 2, 2),
    'E_m4': (3, 3, 2, -4), 'E_m2': (3, 3, 2, -2), 'E_p0': (3, 3, 2, 0),
    'E_p2': (3, 3, 2, 2),
}
EXTRAS = {k: ('y1Conj' if v[3] == 0 else 'y2Conj') for k, v in SPECS.items()
          if k[0] in 'ADE' and v[3] in (0, -2)}
SMOKE = ('A_p0', 'A_m2', 'B_m1', 'D_m2', 'E_p0')


def require(ok: bool, detail: str) -> None:
    if not ok:
        raise ValueError(detail)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_checks(root: Path) -> str:
    p = root / MODEL.relative_to(ROOT)
    require(p.is_file(), 'Missing installed RGBeta model: '+str(p))
    s = p.read_text(encoding='utf-8')
    require('T3RGBetaAddGeneralFermionYukawas[' in s and
            'T3_EXPERIMENTAL_GENERAL_FERMIONS' in s and
            'eps[T3RGBetaRep[dF], #1, #2]' in s and
            '"y1Conj" -> Quiet[BetaTerm[y1Conj, 1]]' in s and
            '"y2Conj" -> Quiet[BetaTerm[y2Conj, 1]]' in s,
            'Integrated model not installed, or charged MF correction missing')
    return digest(p)


def check_payload(model: str, d: dict) -> dict:
    require(isinstance(d, dict) and d.get('status') == 'Success',
            f'{model}: RGBeta did not return Success')
    dims = SPECS[model]
    meta = d.get('metadata')
    require(isinstance(meta, dict) and
            all(meta.get(k) == v for k, v in zip(('dS1', 'dS2', 'dF', 'alpha'), dims)),
            f'{model}: unexpected model metadata')
    require(meta.get('MajoranaFermion') == (dims[3] == -1),
            f'{model}: incorrect Majorana flag')
    report = d.get('report_betas')
    require(isinstance(report, dict), f'{model}: no report_betas dictionary')
    beta = d.get('betas')
    require(isinstance(beta, dict) and set(report) == set(beta),
            f'{model}: raw versus report beta keys mismatch')
    require(meta.get('GeneralScalarStage') == 'UV' and
            meta.get('GeneralScalarModel') == model[0] + '_' +
            ('-' if dims[3] < 0 else '+') + str(abs(dims[3])),
            f'{model}: general scalar UV projection metadata missing')
    required = {'gY', 'g2', 'g3', 'yu', 'yd', 'ye', 'y1', 'y2', 'MF',
                'mS1Sq', 'mS2Sq', 'lambdaT3'}
    require(required.issubset(report), f'{model}: required legacy couplings missing')
    found_extras = {'y1Conj', 'y2Conj'} & set(report)
    want_extras = {EXTRAS[model]} if model in EXTRAS else set()
    require(found_extras == want_extras,
            f'{model}: expected extra Yukawas {sorted(want_extras)}, got {sorted(found_extras)}')
    for k, value in report.items():
        require(isinstance(value, str) and value.strip() and
                '$Failed' not in value and '$Aborted' not in value and
                'BetaTerm[' not in value,
                f'{model}: unevaluated or malformed beta for {k}')
    return {'model': model, 'dimensions_alpha': list(dims),
            'extra_yukawas': sorted(found_extras),
            'beta_count': len(report), 'status': 'VerifiedRGBetaUVReport',
            'physical_matchete_rgbeta_dictionary_certified': False}


def read_report(path: Path) -> dict:
    require(path.is_file(), f'Missing native RGBeta result: {path}')
    with path.open(encoding='utf-8-sig') as handle:
        result = json.load(handle)
    require(isinstance(result, dict), 'RGBeta result is not a dictionary')
    return result


def wolfram_command(root: Path, model: str, result: Path) -> list[str]:
    """Build Wolfram invocation with a positional, non-option alpha token.

    The existing T3RGBetaParseIntegerToken parser accepts m4/m2/m1 for
    negative integers. Passing '-2' directly to wolframscript can be parsed
    as a command-line option and terminate the program before the WL script
    runs. Preserve the runner's A_-2 naming convention separately.
    """
    require(model in SPECS, 'Unknown T3 model: '+model)
    d1, d2, df, alpha = SPECS[model]
    alpha_token = 'm'+str(-alpha) if alpha < 0 else str(alpha)
    native_name = model[0] + '_' + ('-' if alpha < 0 else '+') + str(abs(alpha))
    return ['wolframscript', '-file', str(root / RUNNER.relative_to(ROOT)),
            'UV', str(d1), str(d2), str(df), alpha_token,
            native_name, str(result)]


def run_one(root: Path, model: str, timeout: int) -> dict:
    modeldir = root / OUT.relative_to(ROOT)
    (modeldir / 'native').mkdir(parents=True, exist_ok=True)
    (modeldir / 'logs').mkdir(parents=True, exist_ok=True)
    result = modeldir / 'native' / f'{model}.json'
    log = modeldir / 'logs' / f'{model}.log'
    if result.exists():
        result.unlink()   # never allow stale native outputs to look like a success
    env = os.environ.copy()
    env['T3_EXPERIMENTAL_GENERAL_FERMIONS'] = '1'
    cmd = wolfram_command(root, model, result)
    try:
        p = subprocess.run(cmd, cwd=root, env=env,
                           capture_output=True, text=True, errors='replace',
                           timeout=timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f'{model}: Wolfram timed out; see {log}') from exc
    log.write_text(p.stdout+'\nSTDERR:\n'+p.stderr, encoding='utf-8')
    if p.returncode != 0:
        # Surface the actionable tail; retain complete Wolfram output on disk.
        log_lines = [line.strip() for line in (p.stdout+'\n'+p.stderr).splitlines()
                     if line.strip()]
        tail = '\n'.join(log_lines[-10:])[-1800:]
        raise RuntimeError(f'{model}: Wolfram failed with exit {p.returncode}; '
                           f'see {log}\nWolfram tail:\n{tail}')
    return check_payload(model, read_report(result))


def save_summary(root: Path, rows: list[dict], mode: str, local_sha: str) -> Path:
    outdir = root / OUT.relative_to(ROOT)
    outdir.mkdir(parents=True, exist_ok=True)
    models = {r['model'] for r in rows if r.get('status') == 'VerifiedRGBetaUVReport'}
    expected = set(SPECS) if mode == 'all' else set(SMOKE)
    passed = expected == models and all(r.get('status') == 'VerifiedRGBetaUVReport' for r in rows)
    data = {
        'format': 't3-general-fermion-production-beta-stack-v1',
        'status': ('VerifiedIntegratedRGBetaUVAll16' if passed and mode == 'all' else
                   'VerifiedIntegratedRGBetaUVSmoke' if passed else
                   'IntegratedRGBetaUVPartialOrBlocked'),
        'mode': mode, 'utc_generated_at': datetime.now(timezone.utc).isoformat(),
        'local_model_sha256': local_sha,
        'all_requested_models_succeeded': passed,
        'six_independent_extra_yukawas_present_in_catalogue':
            len([x for x in rows if x.get('status') == 'VerifiedRGBetaUVReport' and x['extra_yukawas']]) == (6 if mode == 'all' else 4),
        'rows': rows,
        'physical_matchete_rgbeta_yukawa_dictionary_certified': False,
        'neutrino_mass_matching_rerun': False,
        'numerical_coupled_general_fermion_uv_evolution_certified': False,
        'note': 'RGBeta full generated scalar quartic projection and opt-in Yukawa betas. Physical Matchete map and coupled extra-Yukawa numerical state remain uncertified.',
    }
    path = outdir / 'summary.json'
    path.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
    return path


def check_only(root: Path, mode: str) -> int:
    h = source_checks(root)
    labels = list(SPECS) if mode == 'all' else list(SMOKE)
    rows = [check_payload(k, read_report(root / OUT.relative_to(ROOT) / 'native' / (k+'.json')))
            for k in labels]
    path = save_summary(root, rows, mode, h)
    print(f'PASS: {len(rows)} pre-existing native reports checked: {path}')
    return 0


def self_test() -> None:
    for name in SPECS:
        d1, d2, df, alpha = SPECS[name]
        p = {'status': 'Success', 'metadata': {
            'dS1': d1, 'dS2': d2, 'dF': df, 'alpha': alpha,
            'MajoranaFermion': alpha == -1,
            'GeneralScalarStage': 'UV',
            'GeneralScalarModel': name[0] + '_' + ('-' if alpha < 0 else '+') + str(abs(alpha))},
            'report_betas': {k:'1' for k in (
                'gY','g2','g3','yu','yd','ye','y1','y2','MF','mS1Sq','mS2Sq','lambdaT3')},
        }
        if name in EXTRAS:
            p['report_betas'][EXTRAS[name]] = '2'
        p['betas'] = dict(p['report_betas'])
        row = check_payload(name, p)
        require(row['extra_yukawas'] == ([EXTRAS[name]] if name in EXTRAS else []),
                'Model catalogue mismatch')
        p['report_betas']['y1Conj' if name not in EXTRAS else
                         ('y2Conj' if EXTRAS[name]=='y1Conj' else 'y1Conj')] = '3'
        p['betas'] = dict(p['report_betas'])
        try:
            check_payload(name, p)
        except ValueError:
            pass
        else:
            raise AssertionError('Unexpected extra Yukawa not detected')
    require(len(EXTRAS)==6 and len(SPECS)==16, 'Wrong T3 catalogue')
    for name, dims in SPECS.items():
        cmd = wolfram_command(ROOT, name, Path('output.json'))
        alpha = dims[3]
        token = cmd[-3]
        expected = 'm'+str(-alpha) if alpha < 0 else str(alpha)
        require(token == expected, f'{name}: invalid Wolfram alpha token {token}')
        require(cmd[-2] == name[0]+'_'+('-' if alpha < 0 else '+')+str(abs(alpha)),
                f'{name}: wrong runner model label')
        require(not token.startswith('-'), f'{name}: raw negative CLI argument')
    print('PASS: 16-model catalogue; 6 extra Yukawa keys; 16 safe Wolfram alpha tokens; negative controls')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    gg = ap.add_mutually_exclusive_group(required=True)
    gg.add_argument('--self-test', action='store_true')
    gg.add_argument('--smoke', action='store_true')
    gg.add_argument('--all', action='store_true')
    gg.add_argument('--check-only', action='store_true')
    ap.add_argument('--timeout', type=int, default=900,
                    help='Maximum Wolfram execution time per model in seconds')
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    mode = 'all' if args.all or args.check_only else 'smoke'
    hash0 = source_checks(ROOT)
    require((ROOT / RUNNER.relative_to(ROOT)).is_file(),
            'Missing RunGeneralT3RGBeta.wl quartic projection runner')
    if args.check_only:
        return check_only(ROOT, mode)
    require(args.timeout > 0, 'timeout must be positive')
    labels = list(SPECS) if mode == 'all' else list(SMOKE)
    rows = []
    for n, label in enumerate(labels, 1):
        try:
            row = run_one(ROOT, label, args.timeout)
            rows.append(row)
            print(f'[{n}/{len(labels)}] PASS {label}: {row["beta_count"]} betas; '
                  f'extras={row["extra_yukawas"]}', flush=True)
        except (ValueError, OSError, RuntimeError) as exc:
            rows.append({'model': label, 'status': 'Blocked',
                         'error': str(exc), 'extra_yukawas': []})
            print(f'[{n}/{len(labels)}] BLOCKED {label}: {exc}', flush=True)
            break  # do not report remaining models as successful
    path = save_summary(ROOT, rows, mode, hash0)
    success = len(rows)==len(labels) and all(x['status']=='VerifiedRGBetaUVReport' for x in rows)
    print('Summary: '+str(path))
    print('Absolute physical Matchete/RGBeta dictionary: NOT CERTIFIED')
    return 0 if success else 3


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print('BLOCKED: '+str(exc), file=sys.stderr)
        raise SystemExit(2)
