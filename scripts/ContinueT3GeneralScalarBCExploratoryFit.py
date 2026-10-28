"""Batch 12: continue existing B/C exploratory Yukawa fits, 18 real variables.

Reuses the unchanged Batch 11 candidate general-scalar RG/mixed-one-scale-C5
point evaluator, existing Sobol/least-squares optimiser, and 100 GeV NuFIT
objective. The saved best point in an existing run is the initial point, NOT
historical Yukawas. The local optimisation coordinates use each seed's scale,
not the historical fixed [-3,3] interval for every derivative.

Outputs are new-only in output/studies/general_scalar_extended/<case>/<run-id>.
Neither the historical fits nor any previous exploratory run are modified.
Physical Matchete/RGBeta coupling and sequential-matching status: UNCERTIFIED.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from scripts import RunT3GeneralScalarBCExploratoryFit as old

ROOT = old.ROOT
FORMAT = 't3-bc-general-scalar-exploratory-yukawa-continuation-v1'
NAMES = old.NAMES
FIELDS = old.FIELDS


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError('B/C exploratory continuation: ' + message)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def previous_paths(case: str, model: str, previous_run_id: str,
                   *, root: Path = ROOT) -> tuple[Path, Path]:
    require(case in old.b9.CASES, 'unknown case')
    require(model in old.b9.MODELS, 'only B_-1/C_-1 supported')
    run_id = old.b9.parse_run_id(previous_run_id)
    folder = Path(root) / 'output/studies/general_scalar_extended' / case / run_id
    tag = old.b9.safe_tag(model)
    return (folder / f'{tag}_exploratory_fit.json',
            folder / f'{tag}_exploratory_fitted_config.json')


def values_from_record(result: Mapping[str, Any]) -> dict[str, float]:
    values = result.get('best_parameters')
    require(isinstance(values, dict) and set(values) == set(NAMES),
            'saved best point must contain exactly 18 named Yukawa coordinates')
    numeric = {name: old.real_finite(values[name], name) for name in NAMES}
    return numeric


def load_previous(case: str, model: str, previous_run_id: str, *,
                  root: Path = ROOT) -> tuple[dict, dict, Path, Path]:
    result_path, config_path = previous_paths(case, model, previous_run_id, root=root)
    require(result_path.is_file(), 'missing previous fit: ' + str(result_path))
    require(config_path.is_file(), 'missing previous fitted config: ' + str(config_path))
    prior = json.loads(result_path.read_text(encoding='utf-8'))
    config = json.loads(config_path.read_text(encoding='utf-8'))
    require(prior.get('format') == old.FORMAT,
            'requires a complete Batch 11 exploratory fit, not a probe or another format')
    require(prior.get('case') == case and prior.get('model') == model and
            prior.get('run_id') == previous_run_id, 'previous run metadata mismatch')
    require(prior.get('status') in ('ExploratoryPoorFit', 'ExploratoryAcceptedFit'),
            'previous run is not a completed fitted result')
    require(prior.get('physical_matchete_rgbeta_map_certified') is False and
            prior.get('full_sequential_general_scalar_matching_performed') is False,
            'previous physical certification flags changed unexpectedly')
    require(prior.get('UV_running_enabled') == old.b9.CASES[case][1],
            'historical UV-running mode changed')
    require(config.get('exploratory_general_scalar_fit', {}).get('run_id') == previous_run_id,
            'previous fitted-config provenance mismatch')
    require(config['exploratory_general_scalar_fit'].get('case') == case and
            config['exploratory_general_scalar_fit'].get('model') == model,
            'previous fitted-config case or model mismatch')
    previous_chi2 = float(prior.get('chi2_final', float('nan')))
    require(math.isfinite(previous_chi2) and previous_chi2 >= 0,
            'invalid saved fit chi-square')
    seed = values_from_record(prior)
    ordinary = config.get('base_state', {}).get('ordinary', {})
    for k in ('y1', 'y2'):
        spec = ordinary.get(k)
        require(isinstance(spec, dict), 'missing previous fitted Yukawas: ' + k)
        re = np.asarray(spec.get('real'), dtype=float)
        im = np.asarray(spec.get('imag'), dtype=float)
        require(re.shape == (3, 3) and im.shape == (3, 3) and
                np.isfinite(re).all() and np.isfinite(im).all() and
                float(np.max(np.abs(im))) <= 1e-13,
                'invalid previous fitted 3x3 real Yukawas: ' + k)
        expected = np.array([seed[f'{k}_{i}{j}'] for i in range(3)
                             for j in range(3)], dtype=float).reshape(3, 3)
        require(np.array_equal(re, expected),
                'previous fitted configuration does not exactly match saved best point: ' + k)
    return prior, config, result_path, config_path


def check_previous_source_hashes(prior: Mapping[str, Any], *, root: Path = ROOT) -> None:
    hashes = prior.get('source_manifest', {}).get('files_sha256')
    require(isinstance(hashes, dict) and len(hashes) >= 5,
            'missing original physics-source hash manifest')
    base = Path(root).resolve()
    for relative, expected in hashes.items():
        require(isinstance(relative, str) and isinstance(expected, str),
                'invalid source hash entry')
        full = (base / relative).resolve()
        require(full.is_relative_to(base), 'source hash path escapes project root')
        require(full.is_file(), 'original physics source missing: ' + relative)
        require(digest(full) == expected,
                'physics source changed since previous fit: ' + relative)


def adaptive_definitions(seed: Mapping[str, float], *, bound: float,
                         relative_window: float, floor_ratio: float,
                         absolute_floor: float) -> tuple[Any, ...]:
    """Local bounded coordinates ensure sensible Yukawa derivatives at both scales."""
    from Numerical.fitting.SobolBenchmarkSearch import BenchmarkParameter
    bound = old.real_finite(bound, 'bound')
    relative_window = old.real_finite(relative_window, 'relative_window')
    floor_ratio = old.real_finite(floor_ratio, 'floor_ratio')
    absolute_floor = old.real_finite(absolute_floor, 'absolute_floor')
    require(bound > 0 and 0 < relative_window <= 4 and
            0 < floor_ratio <= 1 and absolute_floor > 0,
            'invalid adaptive bound settings')
    x = values_from_record({'best_parameters': dict(seed)})
    scale = max(abs(v) for v in x.values())
    require(scale > 0 and scale < bound, 'saved Yukawas are zero or exceed the global bound')
    floor = max(absolute_floor, floor_ratio * scale)
    out = []
    for name in NAMES:
        center = x[name]
        half = max(relative_window * abs(center), floor)
        low = max(-bound, center - half)
        high = min(bound, center + half)
        require(low < center < high, 'empty/truncated parameter interval for ' + name)
        field, index = name.split('_')
        binding = f'ordinary.{field}.real.{index[0]}.{index[1]}'
        out.append(BenchmarkParameter(name, binding, low, high))
    return tuple(out)


def evaluate_candidates(evaluator: Any, seed: dict[str, float], definitions: tuple,
                        *, sobol_points: int, sobol_seed: int) -> tuple[list[dict], list[dict]]:
    from Numerical.fitting.SobolBenchmarkSearch import generate_sobol_parameter_points
    fit = evaluator(seed)
    candidates = [old.result_record(fit, seed, 'saved_previous_best_recomputed')]
    failures: list[dict] = []
    if sobol_points:
        require(sobol_points >= 2 and sobol_points & (sobol_points - 1) == 0,
                '--sobol-points must be zero or a power of two >=2')
        points = generate_sobol_parameter_points(definitions, sobol_points,
                                                  seed=sobol_seed)
        for i, point in enumerate(points, 1):
            try:
                candidates.append(old.result_record(evaluator(point), point,
                                                    f'local_sobol_{i}'))
            except Exception as exc:
                failures.append({'point': i, 'reason': str(exc)[:300]})
    return sorted(candidates, key=lambda rec: float(rec['chi2'])), failures


def result_markdown(result: Mapping[str, Any]) -> str:
    lines = [
        '# General-scalar B/C exploratory fit continuation', '',
        f"Historical case: `{result['case']}`; model: `{result['model']}`; run: `{result['run_id']}`", '',
        f"Restarted from `{result['previous_run_id']}`, not the historical Yukawas.",
        'All 18 real Yukawa coordinates enabled in the local least-squares stage.',
        'Candidate RGBeta UV/EFT1 + mass-mixed one-scale Matchete C5 + SM/Weinberg to 100 GeV.',
        '**Conditional** Stage23/26 signed map for two direct quartics; NOT a certified physical or sequential match.', '',
        f"Prior fit chi2: {result['chi2_previous']:.10g}",
        f"Recomputed start chi2: {result['chi2_start_recomputed']:.10g}",
        f"Continued best chi2: {result['chi2_final']:.10g}",
        f"Chi2 <= 5.86: {'YES' if result['passes_approx_chi2_5p86'] else 'NO'}", '',
        '| Observable | Continued result |', '|---|---:|',
        *[f'| {name} | {result["prediction"][name]:.10g} |' for name in FIELDS], '',
        f"Evaluations (successful): {result['candidate_evaluation_count']}; local Sobol failures: {result['local_sobol_failure_count']}",
        'No previously saved results or historical fit files were modified.', ''
    ]
    return '\n'.join(lines)


def continue_run(case: str, model: str, previous_run_id: str, run_id: str, *,
                 sobol_points: int = 0, sobol_seed: int = 20261029,
                 bound: float = 3.0, relative_window: float = 0.5,
                 floor_ratio: float = 0.02, absolute_floor: float = 1e-8,
                 max_nfev: int = 60, root: Path = ROOT) -> dict:
    # All validation precedes any output write.
    dest = old.fit_destination(case, run_id, root=root)
    require(previous_run_id != run_id, 'new run ID must differ from prior run ID')
    require(sobol_points == 0 or (sobol_points >= 2 and
                                 sobol_points & (sobol_points-1) == 0),
            'invalid --sobol-points')
    require(10 <= max_nfev <= 500, '--max-nfev must be 10..500')
    prior, prior_cfg, fit_path, config_path = load_previous(
        case, model, previous_run_id, root=root)
    check_previous_source_hashes(prior, root=root)
    seed = values_from_record(prior)
    definitions = adaptive_definitions(seed, bound=bound,
         relative_window=relative_window, floor_ratio=floor_ratio,
         absolute_floor=absolute_floor)

    from Numerical.experimental.GeneralT3ScalarRGState import load_registry
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.fitting.SobolBenchmarkSearch import refine_sobol_candidates
    from Numerical.orchestration import RunSmallLambdaRefit as base
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE

    registry = load_registry(old.b7.REGISTRY)
    source = old.b9.study_source(case, model, registry, root=root)
    require(source['cfg_path'].relative_to(Path(root)).as_posix() ==
            prior.get('historical_fitted_config'),
            'previous refit uses a different historical fitted benchmark')
    expected_uv = source['uv_on']
    require(expected_uv == prior['UV_running_enabled'],
            'changed UV running mode')
    manifest = prior.get('source_manifest', {})
    old.b10.check_convention_evidence(
        old.b7.read_json(old.b10.EVIDENCE23),
        old.b7.read_json(old.b10.EVIDENCE26), models=(model,))
    coefficient = Path(DEFAULT_SOURCE) / (
        {'B_-1':'B_m1_final','C_-1':'C_m1_final'}[model] +
        '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    compiled = old.compile_source(coefficient)
    beta = old.beta_bundle(model, source, registry)
    target = OscillationFitTarget.from_json(base.TARGET)
    q = old.real_finite(prior['input_extra_quartic_fraction_of_abs_lambdaT3'],
                        'quartic_fraction')
    u = old.real_finite(prior['input_mu2_fraction_of_min_scalar_m2'],
                        'mu2_fraction')
    require(0 < q <= 1 and 0 <= u <= 0.2, 'invalid previous scalar boundary')
    require(manifest.get('fit_options', {}).get('quartic_fraction') == q and
            manifest.get('fit_options', {}).get('mu2_fraction') == u,
            'previous scalar boundary differs from original fit manifest')
    evaluator = old.CandidateEvaluator(case,model,source,registry,beta,compiled,target,
                               quartic_fraction=q,mu2_fraction=u)
    candidates, sobol_failures = evaluate_candidates(evaluator, seed, definitions,
                           sobol_points=sobol_points, sobol_seed=sobol_seed)
    recomputed = candidates[0] if not sobol_points else next(
        c for c in candidates if c['source'] == 'saved_previous_best_recomputed')
    start_chi2 = float(recomputed['chi2'])
    require(math.isclose(start_chi2, float(prior['chi2_final']), rel_tol=1e-5,
                         abs_tol=1e-5),
            'saved best point no longer reproduces the previous chi2; check source')
    previous_candidate = recomputed
    if candidates[0]['chi2'] <= old.b9.FIT_CUTOFF:
        selected = candidates[0]
        history={'enabled':False,'skipped':True,'reason':'local candidate already accepted'}
    else:
        selected, history = refine_sobol_candidates(
            evaluator=evaluator, observable_names=FIELDS, definitions=definitions,
            candidates=candidates[:8], target_chi2=old.b9.FIT_CUTOFF,
            seed_count=1, active_counts=(18,), max_nfev=(max_nfev,),
            sensitivity_step=0.02)
    # Never regress below the previously recomputed point (including on failures).
    if selected['chi2'] > previous_candidate['chi2']:
        selected = previous_candidate
    final, details = evaluator.evaluate(selected['parameters'], detail=True)
    require(float(final.chi2) <= float(prior['chi2_final']) + 1e-5,
            'continuation did not reproduce or improve previous chi2')
    accepted = bool(float(final.chi2) <= old.b9.FIT_CUTOFF)
    fitted_config = old.trial_record(source,np.array(
        [selected['parameters'][name] for name in NAMES],dtype=float))['cfg']
    fitted_config['benchmark_search']={'enabled':False}
    fitted_config['sensitivity']={'enabled':False}
    fitted_config['exploratory_general_scalar_fit']={
        'status':'ExploratoryAcceptedFit' if accepted else 'ExploratoryPoorFit',
        'case':case, 'model':model, 'run_id':run_id,
        'continued_from':previous_run_id,
        'prior_fit_sha256':digest(fit_path),
        'prior_fitted_config_sha256':digest(config_path),
        'chi2_general_scalar':float(final.chi2),
        'boundary_extra_quartic_fraction':q,
        'boundary_quadratic_mixing_fraction':u,
        'conditional_signed_direct_c5_map_only':True,
        'full_sequential_general_scalar_matching_performed':False,
        'physical_general_scalar_fit_certified':False,
        'original_historical_fit_metadata_retained_for_provenance_only':True}
    result = {
        'format':FORMAT, 'case':case, 'model':model, 'run_id':run_id,
        'previous_run_id':previous_run_id,
        'status':'ExploratoryAcceptedFit' if accepted else 'ExploratoryPoorFit',
        'chi2_previous':float(prior['chi2_final']),
        'chi2_start_recomputed':start_chi2,
        'chi2_final':float(final.chi2),
        'passes_approx_chi2_5p86':accepted,
        'historical_fit_chi2':float(prior['historical_fit_chi2']),
        'prediction':dict(zip(final.observable_names,map(float,final.prediction))),
        'whitened_residual':np.asarray(final.whitened_residual,dtype=float).tolist(),
        'initial_parameters':seed,
        'best_parameters':{k:float(selected['parameters'][k]) for k in NAMES},
        'optimiser':history,
        'local_sobol_failure_count':len(sobol_failures),
        'local_sobol_failures':sobol_failures,
        'candidate_evaluation_count':evaluator.calls,
        'final_diagnostics':details,
        'source_manifest':{
            'previous_fit_sha256':digest(fit_path),
            'previous_config_sha256':digest(config_path),
            'prior_input_file_sha256':manifest['files_sha256'],
            'quartic_fraction':q,'mu2_fraction':u,
            'optimizer':{'all_real_yukawa_parameters_active':18,
               'adaptive_bound_global_abs_max':float(bound),
               'relative_window':float(relative_window),
               'floor_ratio':float(floor_ratio),
               'absolute_floor':float(absolute_floor),
               'sobol_points':sobol_points, 'sobol_seed':sobol_seed,
               'max_nfev':max_nfev}},
        'UV_running_enabled':expected_uv,
        'physical_matchete_rgbeta_map_certified':False,
        'one_scale_matchete_expression_reused':True,
        'full_sequential_general_scalar_matching_performed':False,
        'fits_not_independent_predictions':True,
        'historical_outputs_overwritten':False}
    dest.mkdir(parents=True, exist_ok=False)
    tag=old.b9.safe_tag(model)
    old.b7.save_json(dest/f'{tag}_exploratory_fit.json',result)
    old.b7.save_json(dest/f'{tag}_exploratory_fitted_config.json',fitted_config)
    old.b7.save_json(dest/'summary.json',{
        'format':FORMAT,'case':case,'model':model,'run_id':run_id,
        'previous_run_id':previous_run_id,'chi2':float(final.chi2),
        'passes_approx_chi2_5p86':accepted,
        'physical_c5_certified':False,'historical_outputs_overwritten':False})
    (dest/'summary.md').write_text(result_markdown(result),encoding='utf-8')
    return result


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',choices=tuple(old.b9.CASES),required=True)
    p.add_argument('--model',choices=old.b9.MODELS,required=True)
    p.add_argument('--previous-run-id',required=True)
    p.add_argument('--run-id',required=True)
    p.add_argument('--max-nfev',type=int,default=60,
                   help='Maximum local least-squares function evaluations (default 60).')
    p.add_argument('--relative-window',type=float,default=0.5)
    p.add_argument('--floor-ratio',type=float,default=0.02)
    p.add_argument('--absolute-floor',type=float,default=1e-8)
    p.add_argument('--bound',type=float,default=3.0)
    p.add_argument('--sobol-points',type=int,default=0,
                   help='Optional additional local Sobol points (0, 8, 16, ...).')
    p.add_argument('--sobol-seed',type=int,default=20261029)
    a=p.parse_args()
    result=continue_run(a.case,a.model,a.previous_run_id,a.run_id,
        max_nfev=a.max_nfev, relative_window=a.relative_window,
        floor_ratio=a.floor_ratio,absolute_floor=a.absolute_floor,
        bound=a.bound,sobol_points=a.sobol_points,sobol_seed=a.sobol_seed)
    print(f"{a.case}: {a.model} {result['status']} chi2={result['chi2_final']:.8g}")
    print(f"Saved: {old.OUT_ROOT/a.case/a.run_id/'summary.md'}")

if __name__=='__main__':
    main()
