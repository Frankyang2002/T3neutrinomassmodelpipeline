"""Batch 20: genuinely coupled native RGBeta EFT1 running for B(+1).

Integrates the same 18 renormalisable Stage34 coordinates plus the complex
lambdaT3GeneralQ022000I1 sidecar from the checked Stage33 B(+1) handoff.
Exports native RGBeta EFT1 one-loop equations in a new Wolfram process.
Also evolves a zero-extra-coupling control; NO Wilson matching, neutrino fit,
or cross-package physical coupling-map certification is claimed.

    python -m scripts.RunT3BPlus1GeneralEFT1 --run-id b_p1_eft1_v1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import subprocess
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'output/general_t3_stage22'
OUTPUT = ROOT / 'output/studies/general_scalar_extended/b_p1_eft1_running'
WOLFRAM = ROOT / 'scripts/ExportT3BPlus1EFT1GeneralOneLoop.wl'
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
EXTRA_WL = ROOT / 'RGE/general/generated/B_+1_extra_quartics.wl'
RGBETA_MODEL = ROOT / 'RGE/running/rgbeta/T3RGBetaModel.wl'
RGBETA_COMMON = ROOT / 'RGE/running/rgbeta/T3RGBetaRunnerCommon.wl'
NEW = 'lambdaT3GeneralQ022000I1'
MODEL = 'B_+1'
MU_F = 1e5
MU_S = 1e3
FORMAT = 'batch20-b-p1-coupled-eft1-rgbeta-betas-v1'
RESULT_FORMAT = 'batch20-b-p1-coupled-eft1-running-v1'
CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'
LOOP = 16.0 * math.pi**2


def fail(message: str) -> None:
    raise ValueError('Batch 20 B(+1) coupled EFT1: ' + message)


def read_json(path: Path) -> dict:
    if not path.is_file():
        fail('Missing file: ' + str(path))
    try:
        record = json.loads(path.read_text(encoding='utf-8-sig'))
    except (OSError, ValueError) as exc:
        fail('Invalid JSON ' + str(path) + ': ' + str(exc))
    if not isinstance(record, dict):
        fail('Not a JSON object: ' + str(path))
    return record


def digest(path: Path) -> str:
    if not path.is_file():
        fail('Missing provenance source: ' + str(path))
    return sha256(path.read_bytes()).hexdigest()


def validate_run_id(run_id: str) -> str:
    if not isinstance(run_id, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', run_id):
        fail('Run ID requires 1–64 ASCII alphanumeric/_/- characters')
    return run_id


def validate_registry(registry: dict) -> None:
    models = registry.get('models')
    if not isinstance(models, dict) or len(models) != 16:
        fail('Expected the complete 16-model registry')
    row = models.get(MODEL)
    if (not isinstance(row, dict) or row.get('dimensions') != [2, 2, 1] or
            row.get('alpha') != 1 or row.get('total_additional_quartics') != 1 or
            row.get('preexisting_rgbeta_symbols') != [] or
            [x.get('name') for x in row.get('new_rgbeta_symbols', [])] != [NEW]):
        fail('B(+1) model registry no longer has exactly one new quartic')


def validate_export(payload: dict, state: Any) -> tuple[dict[str, str], Any]:
    """Require exact native 18+1 beta inventory; no frozen sidecar betas."""
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state

    vector, layout = pack_intermediate_scalar_state(state)
    rep = state.representation
    meta = payload.get('metadata')
    expected_meta = {'SharedScalar': False, 'dS1': 2, 'dS2': 2, 'dF': 1,
                     'alpha': 1, 'IntegratedField': 'F',
                     'ActiveBSMFields': ['S1', 'S2'],
                     'ReportBetaConvention': CONVENTION}
    if (payload.get('format') != FORMAT or payload.get('status') != 'Success' or
            payload.get('model') != MODEL or
            payload.get('physical_coupling_mapping_certified') is not False or
            payload.get('wilson_matching_performed') is not False or
            payload.get('stage29_quartic_evolved_in_this_export') is not False or
            not isinstance(meta, dict) or
            any(meta.get(k) != v for k, v in expected_meta.items()) or
            (rep.d_s1, rep.d_s2, rep.d_f, rep.alpha, rep.shared_scalar) !=
            (2, 2, 1, 1, False)):
        fail('Native export metadata differs from B(+1) F-absent EFT1')
    base = tuple(layout.names)
    if (len(base) != 18 or len(set(base)) != 18 or
            payload.get('base_beta_symbols') is None or
            set(payload['base_beta_symbols']) != set(base) or
            len(payload['base_beta_symbols']) != len(base) or
            payload.get('new_quartic_symbols') != [NEW]):
        fail('18+1 native beta symbol inventory differs from packed EFT1 state')
    betas = payload.get('report_betas')
    if (not isinstance(betas, dict) or set(betas) != set(base) | {NEW} or
            len(betas) != 19 or not np.all(np.isfinite(vector)) or
            any(not isinstance(v, str) or not v.strip() or
                    any(err in v for err in ('BetaTerm[','Missing[','$Failed',
                                              '$Aborted','Indeterminate','ComplexInfinity'))
                    for v in betas.values())):
        fail('Native projected RGBeta beta expressions incomplete or unresolved')
    if (NEW not in betas[NEW] or
            not any(NEW in betas[k] for k in base)):
        fail('Full projected EFT1 beta dictionary has no coupled new-quartic response')
    return betas, layout


def sidecar_complex(raw: Any) -> complex:
    if not isinstance(raw, dict) or set(raw) != {'real', 'imag'}:
        fail('Stage33 new complex quartic sidecar must have real/imag')
    try:
        z = complex(float(raw['real']), float(raw['imag']))
    except (TypeError, ValueError, OverflowError) as exc:
        fail('Invalid Stage33 new quartic value: ' + str(exc))
    if not math.isfinite(z.real) or not math.isfinite(z.imag):
        fail('Nonfinite Stage33 quartic sidecar')
    return z


def derivative(t: float, y: np.ndarray, *, layout: Any,
               betas: Mapping[str, str]) -> np.ndarray:
    """Use production InputForm evaluator and beta-vector normalization."""
    from Numerical.core.BetaVector import beta_values_to_derivative
    from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
    from Numerical.running.IntermediateScalarRGBetaEvaluator import (
        intermediate_scalar_state_environment,
    )
    from Numerical.running.IntermediateScalarStateVector import (
        unpack_intermediate_scalar_state,
    )

    y = np.asarray(y, dtype=float)
    if y.shape != (layout.size + 2,) or not np.all(np.isfinite(y)):
        fail('Expected finite 18-block EFT1 vector + one complex sidecar')
    mu = math.exp(float(t))
    state = unpack_intermediate_scalar_state(y[:layout.size], layout, mu_gev=mu)
    env = intermediate_scalar_state_environment(state)
    env[NEW] = complex(float(y[-2]), float(y[-1]))
    evaluated = {name: evaluate_inputform_expression(expr, env)
                 for name, expr in betas.items()}
    baseline = beta_values_to_derivative(
        {k: evaluated[k] for k in layout.names}, layout,
        divide_by_loop_factor=True)
    if np.asarray(baseline).shape != (layout.size,):
        fail('Production EFT1 derivative has incorrect shape')
    value = evaluated[NEW]
    if isinstance(value, np.ndarray):
        fail('New scalar quartic beta evaluated to a matrix')
    try:
        z = complex(value)
    except (ValueError, TypeError, OverflowError) as exc:
        fail('New quartic beta not a scalar: ' + str(exc))
    result = np.concatenate((baseline, [z.real / LOOP, z.imag / LOOP]))
    if not np.all(np.isfinite(result)):
        fail('Nonfinite coupled EFT1 beta vector')
    return result


def integrate(initial: np.ndarray, layout: Any, betas: Mapping[str, str],
              *, points: int = 24):
    if not isinstance(points, int) or points < 3 or points > 200:
        fail('Trajectory points must be integer 3–200')
    y0 = np.asarray(initial, dtype=float)
    if y0.shape != (layout.size+2,) or not np.all(np.isfinite(y0)):
        fail('Invalid coupled-EFT1 initial vector')
    t0, t1 = math.log(MU_F), math.log(MU_S)
    rhs = lambda t, y: derivative(t, y, layout=layout, betas=betas)
    initial_rhs = rhs(t0, y0)
    solution = solve_ivp(rhs, (t0, t1), y0, method='DOP853',
                         rtol=1e-8, atol=1e-11,
                         t_eval=np.linspace(t0,t1,points))
    if not solution.success or solution.y.shape != (y0.size, points) or \
            not np.all(np.isfinite(solution.y)):
        fail('Native coupled EFT1 DOP853 integration failed: ' + solution.message)
    if not np.allclose(solution.y[:,0], y0, atol=2e-13, rtol=1e-13):
        fail('Integrated EFT1 initial vector differs from checked Stage33 boundary')
    return solution, initial_rhs


def state_rows(solution: Any, layout: Any) -> tuple[list[dict], Any]:
    from Numerical.running.IntermediateScalarStateVector import unpack_intermediate_scalar_state
    from scripts import Stage33FermionThresholdHandoff as s33
    rows = []
    for i, t in enumerate(solution.t):
        st = unpack_intermediate_scalar_state(solution.y[:layout.size, i], layout,
                                              mu_gev=float(math.exp(t)))
        rows.append({'mu_gev': st.mu_gev, 'mS1Sq_gev2': st.mS1Sq,
                     'mS2Sq_gev2': st.mS2Sq,
                     'lambdaH': st.sm.lambdaH,
                     'lambdaT3': s33.json_state(st.lambdaT3),
                     NEW: s33.json_state(complex(solution.y[-2,i],solution.y[-1,i]))})
    last = unpack_intermediate_scalar_state(solution.y[:layout.size,-1], layout,
                                             mu_gev=MU_S)
    return rows, last


def baseline_control_compare(reference: Path, zero_state: Any,
                             source_sha: str) -> dict:
    """Optional comparison to existing Stage34, but only for same Stage33 input."""
    if not reference.is_file():
        return {'status':'ReferenceNotAvailable'}
    record = read_json(reference)
    hashes = record.get('input_sha256', {})
    expected_name = 'stage33_B_p1_fermion_threshold_handoff.json'
    keys = [k for k in hashes if isinstance(k,str) and
            k.replace('\\','/').endswith('/'+expected_name)]
    if (record.get('model') != MODEL or
            record.get('status') != 'SuccessSupportedRenormalisableIntermediateRunning' or
            len(keys) != 1 or hashes[keys[0]] != source_sha):
        return {'status':'ReferenceStaleOrIncompatible'}
    old = record['intermediate_final_state']
    fields = ['mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12']
    differences = {name: abs(float(old[name])-float(getattr(zero_state,name)))
                   for name in fields}
    max_relative = max(differences[k]/max(1.0,abs(float(old[k]))) for k in fields)
    if max_relative > 3e-6:
        fail(f'Zero-extra control differs from saved Stage34 reference, max relative {max_relative:.4e}')
    return {'status':'ZeroExtraAgreesWithSavedStage34',
            'max_relative_scalar_parameter_difference':max_relative}


def execute(run_id: str, *, folder: Path = INPUT, root: Path = ROOT,
            wolframscript: str = 'wolframscript', timeout_sec: int = 3600,
            points: int = 24, output_root: Path | None = None) -> Path:
    from scripts import Stage34IntermediateEFTRunning as s34
    from Numerical.running.IntermediateScalarStateVector import (
        pack_intermediate_scalar_state, unpack_intermediate_scalar_state,
    )
    from scripts import Stage33FermionThresholdHandoff as s33
    from Numerical.running.ScalarThresholdBoundary import scalar_threshold_masses

    run_id = validate_run_id(run_id)
    validate_registry(read_json(Path(root)/'RGE/general/generated/general_rgbeta_quartic_registry.json'))
    single, state, extra, reused, hashes = s34.verify_boundary(MODEL, Path(folder))
    if set(extra) != {NEW} or reused != () or state.mu_gev != MU_F:
        fail('Stage33 B(+1) boundary extra-quartic inventory or scale changed')
    z0 = sidecar_complex(extra[NEW])
    base, _ = pack_intermediate_scalar_state(state)
    if output_root is None:
        output_root = Path(root)/'output/studies/general_scalar_extended/b_p1_eft1_running'
    target = Path(output_root)/run_id
    if target.exists():
        fail('Refusing to overwrite existing run ID: '+str(target))
    if timeout_sec < 1 or points < 3 or points > 200:
        fail('Invalid export timeout or number of saved samples')
    sources = [WOLFRAM if root == ROOT else Path(root)/'scripts/ExportT3BPlus1EFT1GeneralOneLoop.wl',
               Path(root)/'RGE/general/generated/general_rgbeta_quartic_registry.json',
               Path(root)/'RGE/general/generated/B_+1_extra_quartics.wl',
               Path(root)/'RGE/running/rgbeta/T3RGBetaModel.wl',
               Path(root)/'RGE/running/rgbeta/T3RGBetaRunnerCommon.wl']
    source_hashes = {str(p): digest(p) for p in sources}
    target.mkdir(parents=True)
    beta_path = target/'native_eft1_betas.json'
    log_path = target/'wolfram_export.log'
    with log_path.open('w', encoding='utf-8') as log:
        try:
            job = subprocess.run([wolframscript,'-file',str(sources[0]),str(beta_path)],
                                 stdout=log, stderr=subprocess.STDOUT,
                                 timeout=timeout_sec, check=False)
        except (OSError, subprocess.TimeoutExpired) as exc:
            fail('Wolfram exporter failed (see '+str(log_path)+'): '+str(exc))
    if job.returncode != 0:
        fail(f'Wolfram exporter exit {job.returncode}; see {log_path}')
    betas, layout = validate_export(read_json(beta_path), state)
    if not np.allclose(base, pack_intermediate_scalar_state(state)[0],rtol=0,atol=0):
        fail('Stage33 EFT1 input changed after native export')
    initial = np.concatenate((base, [z0.real, z0.imag]))
    live, first_rhs = integrate(initial, layout, betas, points=points)
    control_start = initial.copy()
    control_start[-2:] = 0.0
    control, _ = integrate(control_start, layout, betas, points=points)
    # In B(+1) the holomorphic extra quartic is multiplicatively renormalised.
    # If this ceases to be true, the restricted-limit comparison must be revisited.
    if np.max(np.abs(control.y[-2:])) > 5e-11:
        fail('Zero-extra control generates a nonzero new quartic: no restricted-limit equivalence')
    rows, last = state_rows(live,layout)
    _, zero_last = state_rows(control,layout)
    masses = scalar_threshold_masses(last)
    mass_vals = np.asarray(masses.masses_gev, dtype=float)
    if mass_vals.shape != (2,) or not np.all(np.isfinite(mass_vals)) or np.any(mass_vals<=0):
        fail('Invalid or tachyonic B(+1) scalar threshold masses')
    ref_path = Path(folder)/'stage34_B_p1_intermediate_eft_running.json'
    comparison = baseline_control_compare(ref_path,zero_last,
                                          digest(Path(folder)/'stage33_B_p1_fermion_threshold_handoff.json'))
    differences = {key: float(getattr(last,key)-getattr(zero_last,key)) for key in
                   ('mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1',
                    'lambdaH2','lambda12','lambdaH1Adj','lambdaH2Adj','lambda12Adj')}
    from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
    coupling_names = [name for name in layout.names if name.startswith('lambda')]
    max_quartic = max(abs(complex(intermediate_scalar_state_environment(
                     unpack_intermediate_scalar_state(
                         live.y[:layout.size,i],layout,
                         mu_gev=float(math.exp(live.t[i]))))[k]))
                       for i in range(points) for k in coupling_names)
    max_quartic = max(float(max_quartic),float(np.max(np.hypot(live.y[-2],live.y[-1]))))
    inputs = {**hashes, **source_hashes, str(beta_path):digest(beta_path)}
    report = {'format':RESULT_FORMAT,'status':'SuccessCoupledBPlus1GeneralEFT1Running',
              'model':MODEL,'mu_initial_gev':MU_F,'mu_final_gev':MU_S,
              'new_quartic_symbol':NEW,'source_quartic_at_muF':s33.json_state(z0),
              'quartic_at_muS':s33.json_state(complex(live.y[-2,-1],live.y[-1,-1])),
              'initial_derivative_max_abs':float(np.max(abs(first_rhs))),
              'beta_symbols':list(betas),'base_state_real_dimensions':int(layout.size),
              'complex_sidecar_real_dimensions':2,'live_nfev':int(live.nfev),
              'zero_control_nfev':int(control.nfev),'n_saved_points':points,
              'scalar_masses_at_muS_gev':mass_vals.tolist(),
              'max_sampled_quartic_abs':max_quartic,
              'delta_full_minus_zero_at_muS':differences,
              'zero_control_vs_original_stage34':comparison,
              'eft1_final_state':s33.json_state(last),
              'zero_extra_control_final_state':s33.json_state(zero_last),
              'trajectory':rows,'input_sha256':inputs,
              'native_rgbeta_full_coupled_quartic_eft1':True,
              'new_stage29_quartic_evolved_in_eft1':True,
              'wilson_coefficient_matching_performed':False,
              'neutrino_c5_evaluated':False,
              'physical_coupling_mapping_certified':False,
              'limitations':['Stage33 illustrative UV coupling and renormalisable no-finite threshold boundary',
                             'One-loop renormalisable EFT1 only; no Wilson matching or running',
                             'No oscillation refits in this batch',
                             'General scalar Matchete/RGBeta physical map uncertified']}
    (target/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    summary = (f'# Batch 20 — B(+1) coupled general-scalar EFT1 RG\n\n'
               f'Status: `{report["status"]}`. Native one-loop RGBeta 18+1 beta equations.\n\n'
               f'- New quartic `{NEW}` at muF: `{z0.real:.10g} + {z0.imag:.10g}i`; '
               f'at muS: `{live.y[-2,-1]:.10g} + {live.y[-1,-1]:.10g}i`.\n'
               f'- Scalar masses at 1 TeV: `{mass_vals[0]/1000:.8g}` and '
               f'`{mass_vals[1]/1000:.8g}` TeV (unbroken quadratic running masses).\n'
               f'- Native integration: `{live.nfev}` RHS evaluations; zero-extra control: '
               f'`{control.nfev}`.\n'
               f'- Zero control vs original Stage34: `{comparison["status"]}`.\n'
               f'- Max |scalar quartic coordinate| along saved trajectory: `{max_quartic:.7g}`.\n\n'
               f'## Full minus zero-extra at 1 TeV\n\n'
               '| Coordinate | Difference |\n|---|---:|\n' +
               ''.join(f'| `{key}` | {value:.9g} |\n' for key,value in differences.items())+
               '\nNo Wilson C5 matching, no neutrino masses, and no refits in Batch 20. '
               'Matchete/RGBeta physical normalization remains conditional.\n')
    (target/'summary.md').write_text(summary,encoding='utf-8')
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id',required=True)
    parser.add_argument('--input-dir',type=Path,default=INPUT)
    parser.add_argument('--wolframscript',default='wolframscript')
    parser.add_argument('--timeout-sec',type=int,default=3600)
    parser.add_argument('--points',type=int,default=24)
    args = parser.parse_args()
    try:
        folder = execute(args.run_id,folder=args.input_dir,
                         wolframscript=args.wolframscript,
                         timeout_sec=args.timeout_sec,points=args.points)
    except (ValueError,TypeError,OSError,RuntimeError,KeyError,OverflowError) as exc:
        parser.exit(2,str(exc)+'\n')
    print('B(+1) coupled general-scalar EFT1 running:',folder/'summary.md',flush=True)


if __name__ == '__main__':
    main()
