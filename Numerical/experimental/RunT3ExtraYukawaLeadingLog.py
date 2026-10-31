#!/usr/bin/env python3
"""RGBeta-only numerical one-step diagnostic for six new T3 Yukawa matrices.

Consumes the *new* experimental RGBeta UV beta exports and uses existing
fitted-reference numbers as a convenient numerical test point, NOT as a
certified Matchete-to-RGBeta physical conversion. Independent extra Yukawa
matrices are seeded in RGBeta convention for this diagnostic only.

No coupled ODE, matching, C5 calculation, or neutrino-mass claims.

python Numerical/experimental/RunT3ExtraYukawaLeadingLog.py --self-test
python Numerical/experimental/RunT3ExtraYukawaLeadingLog.py --all
"""
from __future__ import annotations

import argparse
import json
from math import log, pi
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

EXPORT_DIR = ROOT / 'Reports/output/workstream_a/integrated_general_yukawa_uv'
OUTPUT = ROOT / 'Reports/output/workstream_a/general_yukawa_leading_log'
CONFIG_DIR = ROOT / 'Numerical/benchmarks/fitted_reference'
ROWS = {
    'A_m2': (1,3,2,-2,'y2Conj'), 'A_p0': (1,3,2,0,'y1Conj'),
    'D_m2': (3,1,2,-2,'y2Conj'), 'D_p0': (3,1,2,0,'y1Conj'),
    'E_m2': (3,3,2,-2,'y2Conj'), 'E_p0': (3,3,2,0,'y1Conj'),
}


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def complex_value(x: object) -> complex:
    if isinstance(x, dict):
        require(set(x) == {'real','imag'}, 'Complex scalar expected real and imag keys')
        z = complex(float(x['real']), float(x['imag']))
    else:
        z = complex(float(x))
    require(bool(np.isfinite([z.real, z.imag]).all()), 'Non-finite complex scalar')
    return z


def matrix_value(x: object) -> np.ndarray:
    require(isinstance(x, dict) and set(x) == {'real','imag'}, 'Matrix expected real/imag')
    re = np.asarray(x['real'], dtype=float)
    im = np.asarray(x['imag'], dtype=float)
    require(re.shape==(3,3) and im.shape==(3,3), 'Expected 3x3 Yukawa/mass flavour matrix')
    require(bool(np.isfinite(re).all() and np.isfinite(im).all()), 'Non-finite input matrix')
    return re+1j*im


def encode_matrix(x: np.ndarray) -> dict:
    require(x.shape == (3,3) and np.isfinite(x).all(), 'Non-finite result matrix')
    return {'real': x.real.tolist(), 'imag': x.imag.tolist()}


def base_environment(config: dict, extra: str, magnitude: float) -> dict:
    state = config['base_state']
    sm = state['sm']
    ordinary = state['ordinary']
    env = {key: float(sm[key]) for key in ('gY','g2','g3','lambdaH')}
    env.update({key: matrix_value(sm[key]) for key in ('yu','yd','ye')})
    env.update({key: matrix_value(ordinary[key]) for key in ('y1','y2','MF')})
    for key in ('mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12'):
        env[key] = float(ordinary[key])
    env['lambdaT3'] = complex_value(ordinary['lambdaT3'])
    optional_quartics = ('lambdaH1Adj','lambdaH2Adj','lambdaS1Adj','lambdaS2Adj',
                         'lambda12Adj','lambda12Cross',
                         'lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
                         'lambdaS1bar2S2bar2','lambdaS1barS2S2bar2',
                         'lambdaS1S1bar2S2bar','lambdaHHdagS1barS2barCross')
    for name in optional_quartics:
        if name in ordinary and ordinary[name] is not None:
            env[name] = complex_value(ordinary[name])
    # A nontrivial trial matrix avoids incorrectly classifying the optional
    # coupling as inert when its own beta is homogeneous in the coupling.
    seed = np.array([[1,0.2+0.15j,-0.1],[0.3,0.8,0.2j],[-0.05j,0.1,0.6]],dtype=complex)
    env[extra] = magnitude*seed
    return env


def evaluate_subset(expressions: dict, env: dict) -> dict[str,np.ndarray]:
    from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
    extra = 'y1Conj' if 'y1Conj' in env else 'y2Conj'
    results = {}
    for name in ('y1','y2',extra,'MF'):
        value = evaluate_inputform_expression(expressions[name], env)
        result = np.asarray(value, dtype=complex)
        require(result.shape == (3,3) and np.isfinite(result).all(),
                'Expected finite 3x3 matrix for beta '+name)
        results[name] = result
    return results


def single_model(model: str) -> dict:
    spec = ROWS[model]
    ds1,ds2,df,alpha,extra = spec
    report_path = EXPORT_DIR / 'native' / (model+'.json')
    require(report_path.is_file(), 'First generate RGBeta model report: '+str(report_path))
    data = json.loads(report_path.read_text(encoding='utf-8-sig'))
    require(data.get('status')=='Success', model+': UV beta export not successful')
    meta = data.get('metadata',{})
    for key, number in zip(('dS1','dS2','dF','alpha'),(ds1,ds2,df,alpha)):
        require(meta.get(key)==number, model+': mismatched '+key)
    betas = data.get('report_betas')
    require(isinstance(betas,dict) and extra in betas, model+': extra beta not present')
    cfg_name = f'T3_dS1_{ds1}_dS2_{ds2}_dF_{df}_alpha_'+('m'+str(-alpha) if alpha<0 else 'p'+str(alpha))
    cfg_path = CONFIG_DIR / (cfg_name+'.json')
    require(cfg_path.is_file(), 'Missing existing fitted reference: '+str(cfg_path))
    config = json.loads(cfg_path.read_text(encoding='utf-8-sig'))
    r = config['representation']
    require((r['d_s1'],r['d_s2'],r['d_f'],r['alpha'])==(ds1,ds2,df,alpha),
            model+': benchmark model representation differs')
    original = config['base_state']['ordinary']
    norm = max(float(np.max(np.abs(matrix_value(original['y1'])))),
               float(np.max(np.abs(matrix_value(original['y2'])))),1e-6)
    seed_magnitude = 0.3 * norm
    env = base_environment(config, extra, seed_magnitude)
    zeros = dict(env, **{extra: np.zeros((3,3),dtype=complex)})
    beta_seed = evaluate_subset(betas, env)
    beta_zero = evaluate_subset(betas, zeros)
    mu_uv = float(config['base_state']['mu_gev'])
    mu_ir = float(config['scales']['mu_fermion_threshold_gev']) * 1.1
    require(np.isfinite(mu_uv) and np.isfinite(mu_ir) and mu_uv>mu_ir>0,
            model+': invalid numerical test scales')
    dlog = log(mu_ir/mu_uv)/(16*pi*pi)
    shifts = {k: dlog*v for k,v in beta_seed.items()}
    evolved = {k:env[k]+shifts[k] for k in beta_seed}
    require(all(np.isfinite(v).all() for v in evolved.values()),
            model+': nonfinite one-step diagnostic')
    return {
        'model':model,'status':'RGBetaOnlyFrozenBackgroundOneStepPassed',
        'alpha':alpha,'extra_coupling':extra,
        'uv_scale_gev':mu_uv,'ir_test_scale_gev':mu_ir,
        'one_loop_log_factor':dlog,'extra_yukawa_seed_magnitude':seed_magnitude,
        'yukawa_delta_frobenius_norm':{k:float(np.linalg.norm(v)) for k,v in shifts.items()},
        'beta_difference_from_extra_zero_frobenius_norm':{
            k:float(np.linalg.norm(beta_seed[k]-beta_zero[k])) for k in beta_seed},
        'beta_at_extra_zero_frobenius_norm':{
            k:float(np.linalg.norm(v)) for k,v in beta_zero.items()},
        'new_yukawa_matrix_initial':encode_matrix(env[extra]),
        'new_yukawa_matrix_one_step':encode_matrix(evolved[extra]),
        'physical_coupling_normalization_certified':False,
        'neutrino_mass_prediction_performed':False,
    }


def self_test() -> None:
    m = {'real': [[1,0,0],[0,2,0],[0,0,3]],
         'imag': [[0,0,0],[0,0,0],[0,0,0]]}
    require(matrix_value(m).shape == (3,3), 'Complex matrix fixture failed')
    require(encode_matrix(matrix_value(m)) == m, 'Complex matrix round trip failed')
    require(complex_value({'real':1,'imag':-2}) == 1-2j,
            'Complex scalar conversion failed')
    require(len(ROWS) == 6 and {v[4] for v in ROWS.values()}=={'y1Conj','y2Conj'},
            'Extra-Yukawa model catalogue incorrect')
    for invalid in ({'real':[[0]],'imag':[[0]]},{'real':[[0]],'imag':[]}):
        try:
            matrix_value(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid matrix fixture accepted')
    print('PASS: six-model numerical registry, exact 3x3 state shape, complex data and negative controls')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument('--self-test',action='store_true')
    grp.add_argument('--all',action='store_true')
    grp.add_argument('--model',choices=list(ROWS))
    opts=parser.parse_args()
    if opts.self_test:
        self_test();return 0
    OUTPUT.mkdir(parents=True,exist_ok=True)
    labels = list(ROWS) if opts.all else [opts.model]
    rows=[]
    for model in labels:
        try:
            row=single_model(model)
            rows.append(row)
            print(f'{model}: PASS beta({row["extra_coupling"]}) one-step '
                  f'|delta|={row["yukawa_delta_frobenius_norm"][row["extra_coupling"]]:.5g}',
                  flush=True)
        except (ValueError,KeyError,FileNotFoundError,TypeError) as exc:
            print(f'{model}: BLOCKED: {exc}',file=sys.stderr,flush=True)
            rows.append({'model':model,'status':'Blocked','reason':str(exc)})
            break
    result={'format':'t3-opt-in-yukawa-leading-log-benchmark-v1',
            'status':('RGBetaOnlySixModelOneStepPassed' if len(rows)==6 and
                     all(row['status']=='RGBetaOnlyFrozenBackgroundOneStepPassed'
                         for row in rows) else
                     'RGBetaOnlyOneModelOneStepPassed' if len(rows)==1 and
                     rows[0]['status']=='RGBetaOnlyFrozenBackgroundOneStepPassed' else
                     'RGBetaOnlyLeadingLogBlocked'),
            'rows':rows,
            'input_basis':'RGBeta-only numerical placeholder: legacy fitted-reference magnitudes treated as RGBeta input coordinates; independent extra matrix seeded for test',
            'one_step_only':True,
            'full_coupled_running':False,
            'physical_neutrino_mass_calculation':False,
            'matchete_rgbeta_physical_coupling_dictionary_certified':False}
    path=OUTPUT/'summary.json'
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Report: '+str(path))
    return 0 if all(r.get('status')=='RGBetaOnlyFrozenBackgroundOneStepPassed' for r in rows) and len(rows)==len(labels) else 3

if __name__=='__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, json.JSONDecodeError) as e:
        print('BLOCKED: '+str(e),file=sys.stderr)
        raise SystemExit(2)
