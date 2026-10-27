"""Stage 31: genuine coupled one-loop UV running for 15 ordinary T3 models.

Uses the production numerical state and evaluator WITHOUT modifying either.
New quartic parameters use a complex sidecar, reused quartics overwrite a
copied fitted benchmark UV state. No cross-package physical basis certificate.

    python -m scripts.Stage31CoupledUVRunning --model A_-2
    python -m scripts.Stage31CoupledUVRunning --summary
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import sympy as sp

from scripts.Stage30CoupledA0UVRunning import (
    ROOT, DATA, REGISTRY, START, STOP, BETA_CONVENTION,
    combined_rhs, integrate_vector, file_hash, read_json,
    verify_stage29_registry_provenance,
)

MODELS = ('A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
          'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
          'E_+0', 'E_+2', 'E_-2', 'E_-4')
BASELINE = 'A_+0'
EXPORT_SOURCE = ROOT / 'scripts/ExportGeneralT3Stage31OneLoop.wl'
WL_MODEL = ROOT / 'RGE/running/rgbeta/T3RGBetaModel.wl'
COMPLEX_REUSED = frozenset({
    'lambdaHHdagS2S2', 'lambdaHHdagS1barS1bar',
    'lambdaS1bar2S2bar2', 'lambdaS1barS2S2bar2',
    'lambdaS1S1bar2S2bar', 'lambdaHHdagS1barS2barCross',
})


def fail(message: str) -> None:
    raise ValueError('Stage31 fail-closed: ' + message)


def inventories(registry: dict, model: str) -> tuple[dict, tuple[str, ...], tuple[str, ...]]:
    totals = registry.get('totals', {})
    if any(totals.get(k) != v for k, v in {
            'matchete_extra_quartics': 50,
            'rgbeta_new_vertices': 38,
            'rgbeta_preexisting_structural': 12}.items()):
        fail('Global registry totals changed')
    models = registry.get('models')
    if not isinstance(models, dict) or len(models) != 16 or model not in MODELS:
        fail('Unexpected model selection or registry model count')
    entry = models.get(model)
    if not isinstance(entry, dict) or entry.get('model') != model:
        fail('Registry model entry missing or inconsistent')
    new = tuple(x['name'] for x in entry.get('new_rgbeta_symbols', []))
    reused = tuple(x['existing_name'] for x in entry.get('preexisting_rgbeta_symbols', []))
    if len(new) + len(reused) != entry.get('total_additional_quartics') or \
            len(set(new + reused)) != len(new) + len(reused) or \
            any(not n.startswith('lambdaT3GeneralQ') for n in new) or \
            any(n not in COMPLEX_REUSED for n in reused):
        fail('Unexpected new/reused quartic inventory for ' + model)
    dimensions = entry.get('dimensions')
    if not isinstance(dimensions, list) or len(dimensions) != 3 or \
            any(type(n) is not int or n not in (1, 2, 3) for n in dimensions) or \
            type(entry.get('alpha')) is not int:
        fail('Invalid registry representation')
    return entry, new, reused


def stage29_values(stage29: dict, entry: dict, new: tuple[str, ...],
                   reused: tuple[str, ...]) -> dict[str, complex]:
    if stage29.get('status') != 'NumericalBoundaryCandidatesOnlyNoRGEIntegration' or \
            (stage29.get('models_checked'), stage29.get('slots_checked'),
             stage29.get('blocks_checked')) != (16, 50, 46) or \
            stage29.get('physical_coupling_mapping_certified') is not False or \
            stage29.get('rge_evolution_performed') is not False or \
            stage29.get('benchmark_scale_gev') != START:
        fail('Stage29 conditional benchmark evidence inconsistent')
    row = stage29.get('models', {}).get(entry['model'], {})
    quartics = row.get('rgbeta_quartics')
    if not isinstance(quartics, dict) or set(quartics) != set(new + reused):
        fail('Stage29 quartic boundary coverage inconsistent for ' + entry['model'])
    result = {}
    for symbol in new + reused:
        obj = quartics[symbol]
        if not isinstance(obj, dict) or obj.get('registry_class') != ('new' if symbol in new else 'reused'):
            fail('Stage29 new/reused classification altered: ' + symbol)
        exact = sp.sympify(obj['exact'])
        if exact.free_symbols or sp.simplify(sp.im(exact)) != 0:
            fail('Stage29 UV value not exactly real: ' + symbol)
        v = float(sp.N(exact, 17))
        if not math.isfinite(v) or not math.isclose(v, obj.get('float', float('nan')),
                                                    rel_tol=1e-12, abs_tol=1e-15):
            fail('Stage29 UV exact/float mismatch: ' + symbol)
        result[symbol] = complex(v)
    return result


def config_path(entry: dict) -> Path:
    d1, d2, df = entry['dimensions']
    alpha = entry['alpha']
    token = ('m' if alpha < 0 else 'p') + str(abs(alpha))
    return ROOT / 'Numerical/benchmarks/fitted_reference' / (
        f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{token}.json')


def prepare_state(config: dict, entry: dict, reused: tuple[str, ...],
                  values: dict[str, complex]):
    """Override only reused quartics in an in-memory copy of existing UV input."""
    from Numerical.fitting.ScanCLI import build_uv_state_from_config
    from Numerical.core.StateVector import pack_uv_state
    d1, d2, df = entry['dimensions']
    rep = {'d_s1': d1, 'd_s2': d2, 'd_f': df,
           'alpha': entry['alpha'], 'shared_scalar': False}
    if config.get('representation') != rep:
        fail('Fitted benchmark representation differs from registry')
    copied = deepcopy(config)
    ordinary = copied.get('base_state', {}).get('ordinary')
    if not isinstance(ordinary, dict):
        fail('Missing ordinary base-state configuration')
    for symbol in reused:
        if symbol not in ordinary:
            fail('Reused quartic not in benchmark UV state: ' + symbol)
        value = values[symbol]
        if symbol in COMPLEX_REUSED:
            ordinary[symbol] = {'real': value.real, 'imag': value.imag}
        else:
            if value.imag != 0:
                fail('Complex boundary for real reused quartic: ' + symbol)
            ordinary[symbol] = value.real
    copied.setdefault('scan', {'bindings': {}})
    state = build_uv_state_from_config(SimpleNamespace(raw=copied), {}).validated()
    if state.mu_gev != START:
        fail('Fitted benchmark UV scale changed')
    vec, layout = pack_uv_state(state)
    if not set(reused).issubset(layout.names):
        fail('Reused quartics not represented by production numerical state')
    return vec, layout


def validate_export(export: dict, entry: dict, new: tuple[str, ...],
                    reused: tuple[str, ...], base_names: tuple[str, ...]) -> dict[str, str]:
    d1, d2, df = entry['dimensions']
    if export.get('format') != 'stage31-multi-model-coupled-uv-rgbeta-betas-v1' or \
       export.get('status') != 'Success' or export.get('model') != entry['model'] or \
       export.get('physical_coupling_mapping_certified') is not False:
        fail('Invalid live Stage31 beta export')
    meta = export.get('metadata')
    if not isinstance(meta, dict) or any(meta.get(k) != v for k, v in {
          'SharedScalar': False, 'dS1': d1, 'dS2': d2, 'dF': df,
          'alpha': entry['alpha'], 'ReportBetaConvention': BETA_CONVENTION}.items()):
        fail('Live RGBeta export representation or one-loop convention changed')
    for field, expected in (('new_quartic_symbols', new),
                            ('reused_quartic_symbols', reused),
                            ('base_beta_symbols', base_names)):
        actual = export.get(field)
        if not isinstance(actual, list) or len(actual) != len(expected) or set(actual) != set(expected):
            fail('Live RGBeta beta inventory mismatch: ' + field)
    if not set(reused).issubset(base_names) or set(new) & set(base_names):
        fail('Reused/new couplings overlap wrong numerical blocks')
    betas = export.get('report_betas')
    if not isinstance(betas, dict) or set(betas) != set(new) | set(base_names) or \
            any(not isinstance(v, str) or not v.strip() or
                any(t in v for t in ('BetaTerm[', 'Missing[', '$Failed', '$Aborted'))
                for v in betas.values()):
        fail('Live full one-loop beta equations absent or unresolved')
    return betas


def run_model(model: str, *, folder: Path = DATA, output: Path | None = None,
              n_points: int = 24) -> dict:
    if model not in MODELS:
        fail('Select one of the 15 remaining models')
    registry = read_json(REGISTRY)
    entry, new, reused = inventories(registry, model)
    stage29_file = folder / 'stage29_convention_uv_handoff.json'
    stage29 = read_json(stage29_file)
    verify_stage29_registry_provenance(stage29)
    values = stage29_values(stage29, entry, new, reused)
    conf_path = config_path(entry)
    config = read_json(conf_path)
    base_vec, layout = prepare_state(config, entry, reused, values)
    safe_name = model.replace('+', 'p').replace('-', 'm')
    beta_file = folder / f'stage31_{safe_name}_full_uv_betas.json'
    export = read_json(beta_file)
    betas = validate_export(export, entry, new, reused, layout.names)
    y0 = np.concatenate((base_vec, np.array([
        part for name in new for part in (values[name].real, values[name].imag)], dtype=float)))

    def rhs(t, y):
        return combined_rhs(t, y, layout=layout, names=new, betas=betas)

    initial_rhs = rhs(math.log(START), y0)
    solved = integrate_vector(y0, rhs, n_points=n_points)
    rows = []
    from Numerical.core.StateVector import unpack_uv_state
    for i in range(solved.y.shape[1]):
        vec = solved.y[:, i]
        mu = float(math.exp(solved.t[i]))
        state = unpack_uv_state(vec[:layout.size], layout, mu_gev=mu)
        extra = {name: {'real': float(vec[layout.size+2*j]),
                        'imag': float(vec[layout.size+2*j+1])}
                 for j, name in enumerate(new)}
        reused_rows = {name: {'real': complex(getattr(state, name)).real,
                              'imag': complex(getattr(state, name)).imag}
                       for name in reused}
        rows.append({'mu_gev': mu, 'gY': float(state.sm.gY),
                     'g2': float(state.sm.g2), 'g3': float(state.sm.g3),
                     'lambdaH': float(state.sm.lambdaH),
                     'extra_quartics': extra, 'reused_quartics': reused_rows})
    changes = {}
    for kind, group in (('new', new), ('reused', reused)):
        key = 'extra_quartics' if kind == 'new' else 'reused_quartics'
        for name in group:
            start = complex(**{'real': rows[0][key][name]['real'],
                               'imag': rows[0][key][name]['imag']})
            finish = complex(**{'real': rows[-1][key][name]['real'],
                                'imag': rows[-1][key][name]['imag']})
            changes[name] = {'registry_class': kind,
                             'start_real': start.real, 'start_imag': start.imag,
                             'final_real': finish.real, 'final_imag': finish.imag,
                             'delta_abs': abs(finish-start), 'final_abs': abs(finish)}
    wl_extension = ROOT / 'RGE/general/generated' / f'{model}_extra_quartics.wl'
    inputs = (stage29_file, beta_file, conf_path, REGISTRY,
              EXPORT_SOURCE, WL_MODEL, wl_extension)
    report = {
        'format': 'stage31-coupled-ordinary-t3-uv-running-v1',
        'status': 'SuccessOneModelNumericalUVRunning',
        'model': model, 'mu_initial_gev': START, 'mu_final_gev': STOP,
        'one_loop_beta_evolution_performed': True,
        'full_coupled_uv_beta_system_evolved': True,
        'new_quartic_couplings_evolved': len(new),
        'reused_quartic_couplings_evolved': len(reused),
        'new_quartic_symbols': list(new), 'reused_quartic_symbols': list(reused),
        'sidecar_extra_real_dimensions': 2 * len(new),
        'baseline_real_dimensions': layout.size,
        'production_uv_state_extended': False,
        'source_one_loop_beta_convention': BETA_CONVENTION,
        'numerical_solver': 'scipy.solve_ivp/DOP853',
        'n_saved_points': len(rows), 'nfev': solved.nfev,
        'rtol': 1e-8, 'atol': 1e-11,
        'initial_rhs_max_abs': float(np.max(np.abs(initial_rhs))),
        'conditional_uv_mapping': True,
        'physical_coupling_mapping_certified': False,
        'barred_triplet_component_intertwiner_certified': False,
        'model_input': conf_path.relative_to(ROOT).as_posix(),
        'fit_input_used_without_refitting': True,
        'reused_quartics_set_from_stage29_boundary': len(reused),
        'input_sha256': {str(p): file_hash(p) for p in inputs},
        'quartic_changes': changes, 'trajectory': rows,
        'limitations': ['UV interval only; no threshold matching',
                        'Conventional cross-package quartic map UNCERTIFIED',
                        'Stage29 illustrative extra quartics, no independent fit',
                        'Existing production UV state/source unchanged',
                        'No neutrino mass or oscillation validation'],
    }
    if output is None:
        output = folder / f'stage31_{safe_name}_coupled_uv_running.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return report


def summarize(*, folder: Path = DATA, output: Path | None = None) -> dict:
    registry = read_json(REGISTRY)
    report_paths = {}
    for model in MODELS:
        safe_name = model.replace('+', 'p').replace('-', 'm')
        report_paths[model] = folder / f'stage31_{safe_name}_coupled_uv_running.json'
    baseline_path = folder / 'stage30_A_p0_coupled_uv_running.json'
    baseline = read_json(baseline_path)
    if baseline.get('status') != 'SuccessOneModelNumericalUVRunning' or \
            baseline.get('model') != BASELINE or \
            baseline.get('new_quartic_couplings_evolved') != 5 or \
            baseline.get('physical_coupling_mapping_certified') is not False:
        fail('Stage30 A_+0 baseline missing or inconsistent')
    entries = {}
    count_new = 0
    count_reused = 0
    for model in MODELS:
        entry, new, reused = inventories(registry, model)
        data = read_json(report_paths[model])
        if data.get('status') != 'SuccessOneModelNumericalUVRunning' or \
                data.get('model') != model or \
                data.get('one_loop_beta_evolution_performed') is not True or \
                data.get('physical_coupling_mapping_certified') is not False or \
                data.get('mu_initial_gev') != START or \
                data.get('mu_final_gev') != STOP or \
                data.get('new_quartic_couplings_evolved') != len(new) or \
                data.get('reused_quartic_couplings_evolved') != len(reused) or \
                set(data.get('quartic_changes', {})) != set(new + reused):
            fail('Stage31 summary contains inconsistent model: ' + model)
        count_new += len(new)
        count_reused += len(reused)
        entries[model] = {'new': len(new), 'reused': len(reused),
                          'nfev': data['nfev'],
                          'higgs_quartic_final': data['trajectory'][-1]['lambdaH'],
                          'max_extra_quartic_final_abs': max(
                              [v['final_abs'] for v in data['quartic_changes'].values()],
                              default=None)}
    if (count_new, count_reused) != (33, 12):
        fail('Remaining quartic inventory is not 33 new + 12 reused')
    summary = {'format': 'stage31-16-model-uv-summary-v1',
               'status': 'Success16ModelConditionalUVRunning',
               'stage31_models': 15, 'stage30_baseline_models': 1,
               'total_models': 16, 'new_quartic_slots_total': count_new + 5,
               'reused_quartic_slots_total': count_reused,
               'total_additional_quartic_slots': count_new + 5 + count_reused,
               'mu_initial_gev': START, 'mu_final_gev': STOP,
               'physical_coupling_mapping_certified': False,
               'baseline_stage30_sha256': file_hash(baseline_path),
               'stage31_sha256': {k: file_hash(v) for k, v in report_paths.items()},
               'models': entries,
               'limitations': ['UV running only, no threshold matching',
                               'Physical Matchete/RGBeta component map UNCERTIFIED']}
    if output is None:
        output = folder / 'stage31_16_model_uv_summary.json'
    output.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--model', choices=MODELS)
    group.add_argument('--summary', action='store_true')
    p.add_argument('--input-dir', type=Path, default=DATA)
    p.add_argument('--n-points', type=int, default=24)
    ns = p.parse_args()
    try:
        if ns.summary:
            result = summarize(folder=ns.input_dir)
            print('Stage31 full UV benchmark: 16/16 models, 38 new + 12 reused slots')
            print('Summary:', ns.input_dir / 'stage31_16_model_uv_summary.json')
        else:
            result = run_model(ns.model, folder=ns.input_dir, n_points=ns.n_points)
            print('Stage31', ns.model, ': one-loop UV integration PASSED;',
                  result['new_quartic_couplings_evolved'], 'new,',
                  result['reused_quartic_couplings_evolved'], 'reused;',
                  'nfev=', result['nfev'])
            print('Physical coupling mapping: UNCERTIFIED')
            safe_name = ns.model.replace('+', 'p').replace('-', 'm')
            print('Report:', ns.input_dir / f'stage31_{safe_name}_coupled_uv_running.json')
    except (ValueError, KeyError, TypeError, OSError, OverflowError,
            json.JSONDecodeError, sp.SympifyError) as exc:
        p.exit(2, str(exc) + '\n')


if __name__ == '__main__':
    main()
