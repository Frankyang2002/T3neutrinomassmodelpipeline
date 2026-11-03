"""Zee-style T1/T2/T3 scans through the EXISTING general-T3 numerical evaluators.

NONMIXED: canonical full-general UV/EFT1 RGBeta and one-scale Matchete C5
           Yukawa-channel evaluator (RunT3GeneralAll16FromReference).
B(-1)/C(-1): existing mixed-scalar, conditional one-scale C5 evaluator.

Each model defaults to its BEST ACCEPTED FIT selected by the saved 16-model
merged manifest (original 11 + continued 5). Every gamma point re-evolves
UV and scalar EFT1. No legacy/restricted beta or C5 fallback, new fit,
Wolfram matching generation, or GitHub modification.
This is candidate general-Lagrangian sensitivity, NOT certified sequential
fermion-first Wilson matching and NOT a fully physical 16-model prediction.

PowerShell from project root:
    python -X utf8 -m scripts.RunT3GeneralZeeSensitivity --models A_-2,B_-1,C_-1 --preflight
    python -X utf8 -m scripts.RunT3GeneralZeeSensitivity --models A_-2 --run-id zee_general_smoke_v1
    python -X utf8 -m scripts.RunT3GeneralZeeSensitivity --models all --run-id zee_general_all_v1
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import csv
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIT_ROOT = ROOT / 'output/studies/general_scalar_extended/thesis_general_all16'
OUT_ROOT = ROOT / 'output/studies/general_scalar_extended/zee_style_general_sensitivity'
CANONICAL_ROOT = ROOT / 'output/studies/general_scalar_extended/complete_scalar_symbolic'
MERGED_FILE = FIT_ROOT / 'continue_general16_v1' / 'merged_all16.json'
MERGED_FORMAT = 't3-general-fit-continuation-18-yukawas-conditional-v1'
MERGED_SOURCES = {
    'OriginalGeneralCampaign': 'thesis16_general_c5_full_v1',
    'GeneralFitContinuation': 'continue_general16_v1',
}
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
       'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
DIRECTIONS = ('T1', 'T2', 'T3')
MODEL_ORDER = tuple(f'{fam}_{a:+d}' for fam, alphas in {
    'A': (-4,-2,0), 'B': (-3,-1,1), 'C': (-3,-1,1),
    'D': (-2,0,2), 'E': (-4,-2,0,2),
}.items() for a in alphas)
MIXED = frozenset(('B_-1','C_-1'))
SCHEMA = 't3-general-zee-compensated-sensitivity-v1'


def tag(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


def load(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(str(path))
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(value, dict):
        raise ValueError('Expected JSON object: ' + str(path))
    return value


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def matrix(spec: dict) -> np.ndarray:
    real = np.asarray(spec['real'], dtype=float)
    imag = np.asarray(spec.get('imag', np.zeros_like(real)), dtype=float)
    if real.shape != (3,3) or imag.shape != real.shape:
        raise ValueError('Expected 3x3 complex Yukawa matrix')
    z = real + 1j * imag
    if not np.isfinite(z).all():
        raise ValueError('Nonfinite Yukawa matrix')
    return z


def scalar(spec: dict) -> complex:
    z = complex(float(spec['real']), float(spec.get('imag',0)))
    if not np.isfinite([z.real,z.imag]).all() or z == 0:
        raise ValueError('lambdaT3 must be finite and nonzero')
    return z


def gamma_grid(lo: float, hi: float, points: int) -> list[float]:
    if not (np.isfinite([lo, hi]).all() and 0 < lo < 1 < hi and
            points >= 3 and points % 2):
        raise ValueError('Require 0 < gamma-min < 1 < gamma-max and odd points >= 3')
    h = points // 2
    grid = [*np.linspace(lo, 1., h+1), *np.linspace(1., hi, h+1)[1:]]
    grid[h] = 1.0
    return [float(x) for x in grid]


def deform(config: dict, direction: str, gamma: float) -> dict:
    """Only lambdaT3/y1/y2 change: all other general UV inputs remain fixed."""
    if direction not in DIRECTIONS or not math.isfinite(gamma) or gamma <= 0:
        raise ValueError('Bad direction/gamma')
    out = deepcopy(config)
    ordinary = out['base_state']['ordinary']
    factors = {'T1': (gamma, gamma**-.5, gamma**-.5),
               'T2': (gamma, 1.0, gamma**-1),
               'T3': (1.0, gamma, gamma**-1)}
    fl, f1, f2 = factors[direction]
    lam = scalar(ordinary['lambdaT3']) * fl
    ordinary['lambdaT3'] = {'real':float(lam.real),'imag':float(lam.imag)}
    for name, factor in (('y1',f1),('y2',f2)):
        arr = matrix(ordinary[name]) * factor
        ordinary[name] = {'real':arr.real.tolist(),'imag':arr.imag.tolist()}
    return out


def leading_tensor(config: dict) -> np.ndarray:
    """Algebraic nominal topology only; extra general Yukawa channels may contribute."""
    o = config['base_state']['ordinary']
    y1, y2 = matrix(o['y1']), matrix(o['y2'])
    return scalar(o['lambdaT3']) * (y1 @ y2.T + y2 @ y1.T)


def verify_nominal_invariance(base: dict, shifted: dict) -> None:
    initial, end = leading_tensor(base), leading_tensor(shifted)
    if not np.allclose(initial, end, rtol=3e-13, atol=3e-13*max(1e-40, np.max(abs(initial)))):
        raise RuntimeError('Nominal lambdaT3*y1*y2 topology changed under T1/T2/T3')


def params_from(config: dict) -> dict[str,float]:
    result = {}
    for name in ('y1','y2'):
        z = matrix(config['base_state']['ordinary'][name])
        if np.max(abs(z.imag)) > 1e-12:
            raise ValueError('Existing general fitter accepts real fitted y1/y2 only')
        result.update({f'{name}_{i}{j}':float(z[i,j].real)
                       for i in range(3) for j in range(3)})
    return result


def fit_selections(models: list[str], fit_run_id: str) -> dict[str,dict]:
    """Select exactly the saved, model-specific best general fit, not a shared seed.

    The merged file lists one *result.json* source per model; the neighbouring
    fitted_config.json and general_uv_overrides.json must belong to that SAME
    source. Never combine continued parameters with original scalar boundaries.
    """
    if fit_run_id != 'merged':
        return {m: {'run_id': fit_run_id, 'source': 'ExplicitCampaign',
                    'manifest': None} for m in models}

    doc = load(MERGED_FILE)
    if (doc.get('format') != MERGED_FORMAT or doc.get('accepted_count') != 16 or
            doc.get('blocked_count') != 0 or not isinstance(doc.get('rows'), list) or
            len(doc['rows']) != len(MODEL_ORDER)):
        raise ValueError('Merged general-fit manifest is not the accepted 16-model campaign')
    indexed = {}
    for row in doc['rows']:
        if not isinstance(row, dict):
            raise ValueError('Malformed merged model row')
        model = row.get('model')
        if model not in MODEL_ORDER or model in indexed:
            raise ValueError('Duplicate or unrecognised model in merged fit manifest: '+str(model))
        source = row.get('source')
        if source not in MERGED_SOURCES:
            raise ValueError('Unrecognised merged fit provenance: '+str(source))
        run_id = MERGED_SOURCES[source]
        expected = (Path('output/studies/general_scalar_extended/thesis_general_all16') /
                    run_id / tag(model) / 'result.json').as_posix()
        declared = str(row.get('result_path','')).replace('\\','/')
        if declared != expected:
            raise ValueError(f'{model}: merged source does not match expected path: {declared}')
        chi2 = float(row.get('chi2_final',float('nan')))
        if (row.get('status') != 'ExploratoryAcceptedFit' or
                not math.isfinite(chi2) or not 0 <= chi2 <= 5.86):
            raise ValueError(f'{model}: merged source does not pass saved fit screen')
        indexed[model] = {'run_id': run_id, 'source': source,
                          'manifest': {'result_path': expected, 'chi2_final': chi2}}
    if set(indexed) != set(MODEL_ORDER):
        raise ValueError('Incomplete merged 16-model fit manifest')
    return {m: indexed[m] for m in models}


def _check_fit_yukawas(config: dict, fit: dict) -> None:
    """Prevent silently using a saved config from a different optimisation."""
    stored = params_from(config)
    optimum = fit.get('best_parameters')
    if not isinstance(optimum, dict) or set(stored) != set(optimum):
        raise ValueError('Fitted Yukawa config and saved optimiser parameters differ')
    for name, value in stored.items():
        if not math.isclose(value,float(optimum[name]),rel_tol=1e-9,abs_tol=1e-13):
            raise ValueError('Fitted Yukawa mismatch: '+name)


def _verify_conditional_matching(fit: dict) -> None:
    """Both historical and continuation result schemas explicitly deny sequential matching."""
    names=('full_general_sequential_matching_performed',
           'full_sequential_general_scalar_matching_performed')
    statuses=[fit[name] for name in names if name in fit]
    if not statuses or any(status is not False for status in statuses):
        raise ValueError('Missing or unexpected general sequential-matching certification')


def _bc_input_fractions(fit: dict) -> tuple[float,float]:
    provenance = fit.get('matching_and_rge_provenance', {})
    if not isinstance(provenance, dict):
        raise ValueError('Malformed continued general-fit provenance')
    q = fit.get('extra_quartic_fraction', provenance.get('quartic_fraction'))
    mu = fit.get('extra_mixed_scalar_bilinear_fraction', provenance.get('mu2_fraction'))
    if q is None or mu is None:
        raise ValueError('Missing B/C saved independent scalar boundary fractions')
    q, mu = float(q),float(mu)
    if not (math.isfinite(q) and 0 < q <= 1 and math.isfinite(mu) and 0 <= mu <= .2):
        raise ValueError('Invalid B/C saved independent scalar boundary fractions')
    return q,mu


def _matrix_from_uv_entries(rows: list) -> np.ndarray:
    values = np.asarray([[complex(float(z['real']),float(z['imag'])) for z in row]
                         for row in rows],dtype=complex)
    if values.shape != (3,3) or not np.isfinite(values).all():
        raise ValueError('Bad saved general UV Yukawa override')
    return values


def source_paths(model: str, fit_run: str, canonical_run: str) -> dict[str,Path]:
    root = FIT_ROOT / fit_run / tag(model)
    canon = CANONICAL_ROOT / canonical_run / tag(model)
    p = {'config':root/'fitted_config.json', 'fit':root/'result.json'}
    if model in MIXED:
        from scripts import RunT3BCGeneralFromReference as br
        from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
        from scripts.RunT3ZeroExtraBoundaryBenchmark import beta_path
        p.update({'reference':br.REF/(br._model_key(model)+'.json'),
                  'uv_beta':beta_path(model,'UV'), 'eft1_beta':beta_path(model,'EFT1'),
                  'matching':Path(DEFAULT_SOURCE)/({'B_-1':'B_m1_final',
                      'C_-1':'C_m1_final'}[model]+'_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')})
    else:
        p.update({'overrides':root/'general_uv_overrides.json',
                  'uv_beta':canon/'uv_rgbeta.json','eft1_beta':canon/'eft1_rgbeta.json',
                  'matching':ROOT/'output/general_matchete_matching'/tag(model)/'c5_coefficient.txt',
                  'matching_status':ROOT/'output/general_matchete_matching'/tag(model)/'general_matching_status.json'})
    return p


def sources(model: str, args, selection: dict | None = None) -> tuple[dict,dict,dict[str,Path]]:
    if selection is None:
        selection = fit_selections([model],args.fit_run_id)[model]
    run_id = selection['run_id']
    paths = source_paths(model,run_id,args.canonical_run_id)
    missing = [str(p) for p in paths.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError('Missing general source(s): '+ '; '.join(missing))
    config, fit = load(paths['config']), load(paths['fit'])
    if fit.get('model') != model or not isinstance(config.get('base_state'),dict):
        raise ValueError('Wrong saved general-fit model/source')
    chi2 = float(fit.get('chi2_final',float('nan')))
    if not math.isfinite(chi2):
        raise ValueError('General source lacks a finite chi2_final')
    _check_fit_yukawas(config,fit)
    _verify_conditional_matching(fit)
    manifest = selection.get('manifest')
    if manifest:
        if not math.isclose(chi2,manifest['chi2_final'],rel_tol=1e-9,abs_tol=1e-13):
            raise ValueError('Merged chi2 differs from saved fitted result')
        if selection['source'] == 'GeneralFitContinuation':
            continuation = config.get('general_fit_continuation',{})
            provenance = fit.get('matching_and_rge_provenance',{})
            expected = ('mixed-scalar-conditional-one-scale' if model in MIXED
                        else 'full-general-matchete-one-scale-channels')
            if (continuation.get('source_model') != model or
                    continuation.get('source_campaign') != 'thesis16_general_c5_full_v1' or
                    continuation.get('source_scalar_boundary_modified') is not False or
                    provenance.get('modified_scalar_boundary') is not False or
                    provenance.get('method') != expected or
                    fit.get('physical_general_c5_certified') is not False):
                raise ValueError('Continued general-fit physics provenance differs')
            original_hash = fit.get('original_result_sha256')
            original_path = FIT_ROOT / 'thesis16_general_c5_full_v1' / tag(model) / 'result.json'
            if original_hash is not None:
                if not original_path.is_file() or digest(original_path) != original_hash:
                    raise ValueError('Continued fit no longer matches original-fit source SHA256')
                paths['original_fit'] = original_path
        elif selection['source'] != 'OriginalGeneralCampaign':
            raise ValueError('Unrecognised general-fit source')
    if model in MIXED:
        _bc_input_fractions(fit)
    else:
        if selection['source'] == 'GeneralFitContinuation':
            general_provenance = fit['matching_and_rge_provenance'].get('general_symbolic_provenance',{})
            if (general_provenance.get('source') != 'canonical_full_general_symbolic' or
                    general_provenance.get('model') != model or
                    general_provenance.get('run_id') != args.canonical_run_id):
                raise ValueError('Continued fit did not use selected canonical general RGEs')
        elif fit.get('general_rges_used') is not True or fit.get('general_matchete_matching_used') is not True:
            raise ValueError('Input was not generated by general RGE + Matchete matching')
        if load(paths['matching_status']).get('status') != 'Matched':
            raise ValueError('Full-general one-scale C5 is not marked Matched')
        overrides = load(paths['overrides'])
        uv = load(paths['uv_beta'])
        if set(overrides) != set(uv.get('report_betas',{})):
            raise ValueError('Saved general UV overrides do not cover entire UV beta inventory')
        for k in ('y1','y2','lambdaT3'):
            if k not in overrides:
                raise ValueError('Missing mandatory general UV override '+k)
        for name in ('y1','y2'):
            if not np.allclose(_matrix_from_uv_entries(overrides[name]),
                               matrix(config['base_state']['ordinary'][name]),
                               rtol=1e-9,atol=1e-13):
                raise ValueError('Saved general UV override/Yukawa config mismatch: '+name)
        lam_override = complex(float(overrides['lambdaT3']['real']),
                               float(overrides['lambdaT3']['imag']))
        if not np.isclose(lam_override,scalar(config['base_state']['ordinary']['lambdaT3']),
                          rtol=1e-12,atol=1e-14):
            raise ValueError('Saved general UV override/scalar config mismatch')
    return config,fit,paths


def target_object():
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET
    return OscillationFitTarget.from_json(TARGET)


def make_runner(model: str, config: dict, fit: dict, paths: dict[str,Path], args):
    target = target_object()
    if model not in MIXED:
        from scripts import RunT3GeneralAll16FromReference as general
        uv, eft, provenance, kernel, matched, scheme = general.nonmixed_inputs(
            model, args.canonical_run_id, generate_general_c5=False)
        if matched.resolve() != paths['matching'].resolve():
            raise ValueError('Unexpected general matching source')
        frozen = load(paths['overrides'])
        def evaluate(shifted: dict, gamma: float):
            override = deepcopy(frozen)
            lam = shifted['base_state']['ordinary']['lambdaT3']
            override['lambdaT3'] = deepcopy(lam)
            for key in ('y1','y2'):
                override[key] = general.entry_matrix(matrix(shifted['base_state']['ordinary'][key]))
            e = general.NonMixedEvaluator(model,uv,eft,kernel,override,target)
            return e.evaluate(params_from(shifted))
        return evaluate,scheme

    from scripts import RunT3BCGeneralFromReference as br
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    # This starts with the tracked reference ONLY to identify the general beta
    # registry; actual numeric inputs are replaced with the saved GENERAL fit.
    reference = br.build_reference_source(model)
    reference['cfg'] = deepcopy(config)
    reference['lambdaT3'] = scalar(config['base_state']['ordinary']['lambdaT3'])
    bundle = bc.beta_bundle(model,reference,reference['registry'])
    compiled = bc.compile_source(paths['matching'])
    initial_lambda_abs = abs(scalar(config['base_state']['ordinary']['lambdaT3']))
    q0, mu0 = _bc_input_fractions(fit)
    if not (0 < q0 <= 1 and 0 <= mu0 <= 0.2):
        raise ValueError('Invalid saved B/C general scalar boundary fractions')
    def evaluate(shifted: dict, gamma: float):
        record = dict(reference)
        record['cfg'] = shifted
        new_lam = scalar(shifted['base_state']['ordinary']['lambdaT3'])
        record['lambdaT3'] = new_lam
        # B/C initial new/reused quartics are defined internally as
        # quartic_fraction * abs(lambdaT3). Counter-scale fraction so the
        # INDEPENDENT physical quartic boundary is unchanged at every gamma.
        q = compensated_quartic_fraction(q0, initial_lambda_abs, new_lam)
        e = bc.CandidateEvaluator('reference_scalar_general', model, record,
                                  reference['registry'],bundle,compiled,target,
                                  quartic_fraction=q,mu2_fraction=mu0)
        return e.evaluate(params_from(shifted))
    return evaluate,'conditional-mixed-scalar-one-scale-signed-quartic-map'


def compensated_quartic_fraction(q0: float, initial_lambda: complex,
                                  shifted_lambda: complex) -> float:
    """Keep B/C independent quartic boundary values fixed if lambdaT3 moves."""
    old, new = abs(initial_lambda), abs(shifted_lambda)
    if not (math.isfinite(q0) and 0 < q0 <= 1 and old > 0 and new > 0):
        raise ValueError('Invalid B/C quartic compensation')
    fraction = q0 * old / new
    if not 0 < fraction <= 1:
        raise ValueError('Cannot keep B/C independent quartics fixed: effective fraction > 1')
    return fraction


def prediction(row) -> dict[str,float]:
    vals = dict(zip(row.observable_names, map(float,row.prediction)))
    if set(OBS)-set(vals) or not all(math.isfinite(vals[k]) for k in OBS):
        raise ValueError('Nonfinite/missing neutrino predictions')
    return {k:vals[k] for k in OBS}


def plot_solar(rows: list[dict], output: Path, model: str) -> None:
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(6.9,5.0))
    for direction in DIRECTIONS:
        points = [r for r in rows if r['direction']==direction and r['status']=='Success']
        points.sort(key=lambda r:r['gamma'])
        if not points: continue
        ax.plot([p['sin2_theta12'] for p in points],
                [p['delta_m21_sq_ev2'] for p in points], '-o',ms=3,label=direction)
    ax.set(xlabel=r'$\sin^2\theta_{12}$',ylabel=r'$\Delta m^2_{21}\ (\mathrm{eV}^2)$',
           title=f'{model}: general-Lagrangian candidate coupling sensitivity')
    ax.grid(alpha=0.2);ax.legend();fig.tight_layout()
    fig.savefig(output/'solar_sensitivity.png',dpi=180)
    plt.close(fig)


def baseline_consistency(gamma_one: list[dict], saved: dict, tolerance: float = 1e-5) -> None:
    if len(gamma_one)!=3 or any(r['status']!='Success' for r in gamma_one):
        raise RuntimeError('At least one gamma=1 baseline failed')
    for k in OBS:
        expected = float(saved['prediction'][k])
        for r in gamma_one:
            if abs(r[k]-expected) > max(1e-10, tolerance*abs(expected)):
                raise RuntimeError(f'gamma=1 {k}: recomputed={r[k]:g}, saved={expected:g}')


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--models', default='A_-2',help='Comma-separated, or all')
    p.add_argument('--fit-run-id',default='merged',
                   help='Default: merged accepted best-fit source per model; alternatively a single campaign run ID')
    p.add_argument('--canonical-run-id',default='full16_general_v1')
    p.add_argument('--run-id', default='zee_general_pilot_v1')
    p.add_argument('--gamma-min',type=float,default=0.99)
    p.add_argument('--gamma-max',type=float,default=1.01)
    p.add_argument('--points',type=int,default=11)
    p.add_argument('--preflight',action='store_true',help='Read-only check, no outputs/physics')
    args = p.parse_args(argv)
    wanted = list(MODEL_ORDER) if args.models == 'all' else [x.strip() for x in args.models.split(',') if x.strip()]
    if not wanted or len(wanted)!=len(set(wanted)) or set(wanted)-set(MODEL_ORDER):
        p.error('Invalid models list')
    if not all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',v)
               for v in (args.run_id,args.fit_run_id,args.canonical_run_id)):
        p.error('Invalid run identifier')
    try: grid = gamma_grid(args.gamma_min,args.gamma_max,args.points)
    except ValueError as exc: p.error(str(exc))
    selected = [m for m in MODEL_ORDER if m in wanted]
    try:
        selections = fit_selections(selected,args.fit_run_id)
    except (ValueError,FileNotFoundError) as exc:
        p.error(str(exc))
    if args.preflight:
        failed=0
        for model in selected:
            try:
                _,fit,paths = sources(model,args,selections[model])
                print(f'READY {model}: {selections[model]["source"]} '+
                      f'({selections[model]["run_id"]}), general fitted input, UV/EFT1 betas, '
                      f'C5 matching, chi2={float(fit["chi2_final"]):.7g}; '
                      f'matching={"conditional B/C" if model in MIXED else "full-general one-scale"}')
            except Exception as exc:
                failed+=1
                print(f'BLOCKED {model}: {type(exc).__name__}: {exc}')
        print(f'General source preflight: ready={len(selected)-failed}, blocked={failed}')
        return 2 if failed else 0

    dest = OUT_ROOT / args.run_id
    if dest.exists():
        p.error(f'Refusing to overwrite {dest}; choose a new --run-id')
    dest.mkdir(parents=True)
    all_rows=[]; summaries=[]
    for model in selected:
        folder = dest / tag(model)
        folder.mkdir()
        rows=[]
        try:
            config,fit,paths=sources(model,args,selections[model])
            evaluator,scheme=make_runner(model,config,fit,paths,args)
            for direction in DIRECTIONS:
                for gamma in grid:
                    r={'model':model,'direction':direction,'gamma':gamma,'status':'Failed',
                       'chi2':None,'error':''}
                    try:
                        altered=deform(config,direction,gamma)
                        verify_nominal_invariance(config,altered)
                        evaluated=evaluator(altered,gamma)
                        r.update(prediction(evaluated))
                        score=float(evaluated.chi2)
                        if not math.isfinite(score):raise ValueError('Nonfinite chi2')
                        r.update(chi2=score,status='Success')
                    except Exception as exc:
                        r['error']=f'{type(exc).__name__}: {exc}'[:750]
                    rows.append(r)
                    print(f'{model} {direction} gamma={gamma:.6g} {r["status"]}',flush=True)
            baseline_consistency([r for r in rows if r['gamma']==1.0],fit)
            failures=sum(r['status']!='Success' for r in rows)
            status='Complete' if failures==0 else 'PartialSuccess'
            plot_solar(rows,folder,model)
            save(folder/'scan.json',{'schema':SCHEMA,'model':model,'matching':scheme,
                 'status':status,'baseline_chi2_saved':float(fit['chi2_final']),
                 'baseline_acceptable_chi2_le_5p86':float(fit['chi2_final'])<=5.86,
                 'source_sha256':{name:digest(path) for name,path in paths.items()},
                 'fit_source': selections[model],
                 'merged_manifest_sha256':digest(MERGED_FILE) if args.fit_run_id=='merged' else None,
                 'gamma_grid':grid,'points':rows,
                 'nominal_invariant':'lambdaT3 * (y1 y2.T + y2 y1.T) at UV',
                 'extra_scalar_quartics_and_bilinears_held_fixed':True,
                 'extra_conjugate_yukawas_held_fixed':True,
                 'full_general_sequential_matching_performed':False,
                 'scope':'conditional general UV/EFT1 RG and one-scale C5; not full sequential EFT'})
            summary={'model':model,'status':status,'matching':scheme,'failed_points':failures,
                     'baseline_chi2_saved':float(fit['chi2_final']),
                     'fit_source':selections[model]}
        except Exception as exc:
            summary={'model':model,'status':'Blocked', 'failed_points':None,
                     'error':f'{type(exc).__name__}: {exc}'[:1000],
                     'fit_source':selections[model]}
            print('BLOCKED',model,summary['error'],flush=True)
        all_rows.extend(rows);summaries.append(summary)
        save(dest/'summary.json',{'schema':SCHEMA,'models':summaries,'gamma_grid':grid,
                                   'fit_run_id_requested':args.fit_run_id,
                                   'merged_manifest_sha256':digest(MERGED_FILE) if args.fit_run_id=='merged' else None})
    with (dest/'points.csv').open('w',newline='',encoding='utf-8') as f:
        fields=('model','direction','gamma','status','chi2',*OBS,'error')
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        w.writerows([{k:r.get(k,'') for k in fields} for r in all_rows])
    complete=sum(r['status']=='Complete' for r in summaries)
    print(f'Completed {complete}/{len(selected)} model scans: {dest}')
    return 0 if complete==len(selected) else 2


if __name__ == '__main__':
    raise SystemExit(main())
