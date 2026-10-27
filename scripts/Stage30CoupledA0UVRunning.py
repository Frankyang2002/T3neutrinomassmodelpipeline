"""Stage 30: isolated, coupled one-loop UV integration for A_+0 + 5 new quartics.

Requires an actual fresh RGBeta export after registering the generated A_+0
quartics. Uses the existing production UV state, parser, beta packer and
benchmark config unchanged. New quartics live only in an additional 10-real
coordinate sidecar, so both standard and new beta functions backreact.

  python -m scripts.Stage30CoupledA0UVRunning
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Mapping

import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
DATA = Path('output/general_t3_stage22')
MODEL = 'A_+0'
CONFIG = ROOT / 'Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_p0.json'
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
BETA_SOURCE = ROOT / 'scripts/ExportGeneralT3Stage30A0OneLoop.wl'
WL_MODEL = ROOT / 'RGE/running/rgbeta/T3RGBetaModel.wl'
WL_EXTRA = ROOT / 'RGE/general/generated/A_+0_extra_quartics.wl'
START = 1e7
STOP = 1e5
BETA_CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'


def fail(message: str) -> None:
    raise ValueError('Stage30 fail-closed: ' + message)


def read_json(path: Path) -> dict:
    if not path.is_file():
        fail('Missing input ' + str(path))
    result = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(result, dict):
        fail('Input is not a JSON object: ' + str(path))
    return result


def file_hash(path: Path) -> str:
    if not path.is_file():
        fail('Missing provenance file: ' + str(path))
    return sha256(path.read_bytes()).hexdigest()



def verify_stage29_registry_provenance(stage29: dict, registry_file: Path = REGISTRY) -> None:
    """Compare a Stage29 source digest using platform-independent path keys.

    Stage29 serializes str(Path.relative_to(root)), which uses backslashes on
    Windows and slashes on POSIX. Require one unambiguous matching key and the
    actual SHA-256 digest; do not accept a missing key as a changed file.
    """
    recorded = stage29.get('input_sha256')
    if not isinstance(recorded, dict):
        fail('Stage29 provenance hashes missing')
    canonical = 'RGE/general/generated/general_rgbeta_quartic_registry.json'
    matching = [key for key in recorded if isinstance(key, str) and
                key.replace('\\', '/') == canonical]
    if len(matching) != 1:
        fail('Stage29 registry provenance key absent or ambiguous (Windows/POSIX path)')
    claimed = recorded[matching[0]]
    actual = file_hash(registry_file)
    if not isinstance(claimed, str) or claimed != actual:
        fail('Stage29 registry source changed: regenerate Stage29 handoff')


def quartic_names(registry: dict) -> tuple[str, ...]:
    totals = registry.get('totals', {})
    if totals.get('rgbeta_new_vertices') != 38 or totals.get('matchete_extra_quartics') != 50:
        fail('Unexpected global quartic registry inventory')
    entry = registry.get('models', {}).get(MODEL)
    if not isinstance(entry, dict) or entry.get('total_additional_quartics') != 5 or \
            entry.get('preexisting_rgbeta_symbols') != []:
        fail('A_+0 registry is not 5 new/0 reused')
    names = tuple(x['name'] for x in entry['new_rgbeta_symbols'])
    if len(names) != 5 or len(set(names)) != 5 or \
            any(not n.startswith('lambdaT3GeneralQ') for n in names):
        fail('Invalid A_+0 new quartic symbols')
    return names


def uv_quartics(stage29: dict, registry: dict) -> dict[str, complex]:
    names = quartic_names(registry)
    if stage29.get('status') != 'NumericalBoundaryCandidatesOnlyNoRGEIntegration' or \
            (stage29.get('models_checked'), stage29.get('slots_checked'),
             stage29.get('blocks_checked')) != (16, 50, 46) or \
            stage29.get('physical_coupling_mapping_certified') is not False or \
            stage29.get('rge_evolution_performed') is not False or \
            stage29.get('benchmark_scale_gev') != START:
        fail('Stage29 numerical boundary evidence missing or unexpectedly certified')
    values = stage29.get('models', {}).get(MODEL, {}).get('rgbeta_quartics')
    if not isinstance(values, dict) or set(values) != set(names):
        fail('A_+0 Stage29 boundary coverage differs from exact registry')
    out = {}
    for name in names:
        value = values[name]
        if value.get('registry_class') != 'new':
            fail('Stage29 A_+0 slot is not a new quartic: ' + name)
        exact = sp.sympify(value['exact'])
        if exact.free_symbols or sp.simplify(sp.im(exact)) != 0:
            fail('Stage29 extra boundary is not an exact real number: ' + name)
        number = float(sp.N(exact, 17))
        if not math.isfinite(number) or not math.isclose(number, value['float'], rel_tol=1e-12, abs_tol=1e-15):
            fail('Stage29 exact and float boundaries disagree: ' + name)
        out[name] = complex(number)
    return out


def validate_export(export: dict, expected_new: tuple[str, ...], base_names: tuple[str, ...]) -> dict[str, str]:
    if export.get('format') != 'stage30-A-p0-coupled-uv-rgbeta-betas-v1' or \
            export.get('status') != 'Success' or export.get('model') != MODEL or \
            export.get('physical_coupling_mapping_certified') is not False:
        fail('Missing/invalid live Stage30 one-loop beta export')
    meta = export.get('metadata', {})
    if not isinstance(meta, dict) or any(meta.get(k) != v for k, v in {
        'SharedScalar': False, 'dS1': 1, 'dS2': 3, 'dF': 2, 'alpha': 0,
        'ReportBetaConvention': BETA_CONVENTION,
    }.items()):
        fail('RGBeta export model or gauge-beta convention changed')
    if not isinstance(export.get('new_quartic_symbols'), list) or \
            set(export['new_quartic_symbols']) != set(expected_new) or \
            len(export['new_quartic_symbols']) != len(expected_new):
        fail('Exported five-coupling inventory differs from Stage29')
    if not isinstance(export.get('base_beta_symbols'), list) or \
            set(export['base_beta_symbols']) != set(base_names) or \
            len(export['base_beta_symbols']) != len(base_names):
        fail('Full production UV beta inventory differs from packed numerical state')
    betas = export.get('report_betas')
    if not isinstance(betas, dict) or set(betas) != set(expected_new) | set(base_names) or \
            len(betas) != len(expected_new) + len(base_names) or \
            any(not isinstance(v, str) or not v.strip() for v in betas.values()):
        fail('Incomplete coupled one-loop beta dictionary')
    if any('BetaTerm[' in v or '$Failed' in v or 'Missing[' in v for v in betas.values()):
        fail('Unresolved RGBeta function or failure leaked into beta expressions')
    return betas


def combined_rhs(t: float, y: np.ndarray, *, layout: Any, names: tuple[str, ...],
                 betas: Mapping[str, str]) -> np.ndarray:
    """Evaluate all coupled one-loop betas using production RGBeta arithmetic."""
    from Numerical.core.StateVector import unpack_uv_state
    from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression, state_environment
    from Numerical.core.BetaVector import beta_values_to_derivative

    base_size = layout.size
    if y.shape != (base_size + 2*len(names),):
        fail('Extended UV solver vector has an unexpected shape')
    if not np.all(np.isfinite(y)):
        fail('Nonfinite extended numerical state')
    mu = float(np.exp(t))
    state = unpack_uv_state(y[:base_size], layout, mu_gev=mu)
    env = state_environment(state)
    env.update({name: complex(y[base_size+2*i], y[base_size+2*i+1])
                for i, name in enumerate(names)})
    evaluated = {name: evaluate_inputform_expression(expr, env)
                 for name, expr in betas.items()}
    base = beta_values_to_derivative(
        {name: evaluated[name] for name in layout.names}, layout,
        divide_by_loop_factor=True)
    extra = np.empty(2*len(names), dtype=float)
    for i, name in enumerate(names):
        value = evaluated[name]
        if isinstance(value, np.ndarray):
            fail('New scalar-quartic beta unexpectedly became matrix-valued: ' + name)
        try:
            z = complex(value)
        except (TypeError, ValueError, OverflowError) as exc:
            fail('New scalar beta not numeric: ' + name + ': ' + str(exc))
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            fail('Nonfinite quartic beta: ' + name)
        extra[2*i] = z.real / (16.0*math.pi**2)
        extra[2*i+1] = z.imag / (16.0*math.pi**2)
    result = np.concatenate((base, extra))
    if not np.all(np.isfinite(result)):
        fail('Nonfinite full coupled beta vector')
    return result


def integrate_vector(y0: np.ndarray, rhs: Callable[[float, np.ndarray], np.ndarray],
                     *, n_points: int = 24, rtol: float = 1e-8, atol: float = 1e-11):
    """Small deterministic ODE kernel; actual betas are passed from RGBeta."""
    if not (np.ndim(y0) == 1 and len(y0) > 0 and np.all(np.isfinite(y0))):
        fail('Initial ODE vector is not finite and one-dimensional')
    if n_points < 3 or n_points > 200 or rtol <= 0 or atol <= 0:
        fail('Invalid ODE sampling or tolerances')
    t0, t1 = math.log(START), math.log(STOP)
    result = solve_ivp(rhs, (t0, t1), np.asarray(y0, dtype=float),
                       method='DOP853', rtol=rtol, atol=atol,
                       t_eval=np.linspace(t0,t1,n_points), max_step=0.35)
    if not result.success or result.y.shape != (len(y0),n_points) or \
            not np.all(np.isfinite(result.y)):
        fail('Full one-loop integration failed: ' + result.message)
    return result


def benchmark_state(config: dict):
    from Numerical.fitting.ScanCLI import build_uv_state_from_config
    from Numerical.core.StateVector import pack_uv_state
    if config.get('representation') != {
            'd_s1': 1, 'd_s2': 3, 'd_f': 2, 'alpha': 0, 'shared_scalar': False}:
        fail('Fitted benchmark representation changed; expected A_+0')
    raw = dict(config)
    raw.setdefault('scan', {'bindings': {}})
    initial = build_uv_state_from_config(SimpleNamespace(raw=raw), {}).validated()
    if initial.mu_gev != START:
        fail('Fitted benchmark UV scale does not match Stage29 handoff')
    vec, layout = pack_uv_state(initial)
    return initial, vec, layout


def run(*, folder: Path = DATA, output: Path | None = None,
        export_file: Path | None = None, config_file: Path = CONFIG,
        n_points: int = 24, rtol: float = 1e-8, atol: float = 1e-11) -> dict:
    stage29_file = folder / 'stage29_convention_uv_handoff.json'
    beta_file = export_file or (folder / 'stage30_A_p0_full_uv_betas.json')
    registry = read_json(REGISTRY)
    stage29 = read_json(stage29_file)
    export = read_json(beta_file)
    config = read_json(config_file)
    names = quartic_names(registry)
    boundary = uv_quartics(stage29, registry)
    verify_stage29_registry_provenance(stage29)
    initial, base_vec, layout = benchmark_state(config)
    betas = validate_export(export, names, layout.names)
    y0 = np.concatenate((base_vec,np.array([v for name in names for v in
                (boundary[name].real,boundary[name].imag)],dtype=float)))
    def rhs(t, y):
        return combined_rhs(t, y, layout=layout, names=names, betas=betas)
    # Independent start-of-interval evaluation catches unsupported expression
    # heads or stale parameter names before invoking the numerical integrator.
    start_derivative = rhs(math.log(START), y0)
    result = integrate_vector(y0, rhs, n_points=n_points, rtol=rtol, atol=atol)
    rows = []
    from Numerical.core.StateVector import unpack_uv_state
    for col in range(result.y.shape[1]):
        vec = result.y[:,col]
        scale = float(math.exp(result.t[col]))
        sm = unpack_uv_state(vec[:layout.size],layout,mu_gev=scale).sm
        values = {name: {'real': float(vec[layout.size+2*i]),
                          'imag': float(vec[layout.size+2*i+1])}
                  for i, name in enumerate(names)}
        rows.append({'mu_gev': scale, 'gY': float(sm.gY), 'g2': float(sm.g2),
                     'g3': float(sm.g3), 'lambdaH': float(sm.lambdaH),
                     'extra_quartics': values})
    initial_values = rows[0]['extra_quartics']
    final_values = rows[-1]['extra_quartics']
    changes = {}
    for name in names:
        u = complex(initial_values[name]['real'],initial_values[name]['imag'])
        v = complex(final_values[name]['real'],final_values[name]['imag'])
        changes[name] = {'start_real':u.real,'start_imag':u.imag,
                         'final_real':v.real,'final_imag':v.imag,
                         'delta_abs':abs(v-u),'final_abs':abs(v)}
    inputs = [stage29_file,beta_file,config_file,REGISTRY,BETA_SOURCE,WL_MODEL,WL_EXTRA]
    report = {
        'format':'stage30-A-p0-coupled-one-loop-running-v1',
        'status':'SuccessOneModelNumericalUVRunning',
        'model':MODEL, 'mu_initial_gev':START,'mu_final_gev':STOP,
        'one_loop_beta_evolution_performed':True,
        'full_coupled_uv_beta_system_evolved':True,
        'new_quartic_couplings_evolved':5,
        'production_uv_state_extended':False,
        'sidecar_extra_real_dimensions':10,
        'baseline_real_dimensions':layout.size,
        'source_one_loop_beta_convention': BETA_CONVENTION,
        'numerical_solver':'scipy.solve_ivp/DOP853',
        'n_saved_points':len(rows), 'nfev':result.nfev,
        'rtol':rtol,'atol':atol,
        'initial_rhs_max_abs':float(np.max(np.abs(start_derivative))),
        'conditional_uv_mapping':True,
        'physical_coupling_mapping_certified':False,
        'barred_triplet_component_intertwiner_certified':False,
        'model_input': 'Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_p0.json',
        'fit_input_used_without_refitting':True,
        'input_sha256':{str(p):file_hash(p) for p in inputs},
        'new_quartic_changes':changes,
        'trajectory':rows,
        'limitations':[ 'A_+0 only, UV interval only; no fermion/scalar threshold matching',
            'Extra quartic Matchete/RGBeta map depends on provisional component convention',
            'Illustrative additional quartics from Stage29, not an independent physical fit',
            'Production ODE state and RGBeta source are unchanged',
            'No neutrino-mass calculation or oscillation validation'],
    }
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,default=DATA)
    p.add_argument('--betas',type=Path,default=None)
    p.add_argument('--config',type=Path,default=CONFIG)
    p.add_argument('--output',type=Path,
        default=DATA/'stage30_A_p0_coupled_uv_running.json')
    p.add_argument('--n-points',type=int,default=24)
    ns=p.parse_args()
    try:
        r=run(folder=ns.input_dir,output=ns.output,export_file=ns.betas,
              config_file=ns.config,n_points=ns.n_points)
    except (ValueError,KeyError,RuntimeError,OverflowError,OSError,
            json.JSONDecodeError,sp.SympifyError) as exc:
        p.exit(2,str(exc)+'\n')
    print('Stage30 A_+0 coupled one-loop running: 5/5 new quartics +',
          r['baseline_real_dimensions'],'baseline real ODE coordinates')
    print('Scale: 1e7 -> 1e5 GeV; DOP853; nfev:',r['nfev'])
    for name,item in r['new_quartic_changes'].items():
        print(name+':',format(item['start_real'],'.8g'),'->',
              format(item['final_real'],'.8g'),
              'Im(final)=',format(item['final_imag'],'.4g'))
    print('Physical coupling map: UNCERTIFIED (convention-labelled)')
    print('Report:',ns.output)


if __name__ == '__main__':
    main()
