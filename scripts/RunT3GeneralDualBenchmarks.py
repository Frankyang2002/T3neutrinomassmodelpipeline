"""Compare D(-2) and A(-4) fitted Yukawa textures across 16 general T3 models.

Read-only source models: general thesis16_general_c5_full_v1 fitted_config.json.
Each source fixes the 18 real UV y1/y2 entries for every target representation.
No fit/refit, no source overwrite, no GitHub access. All general scalar and
fermion interactions are retained *when permitted by target symmetries*.

Uses the existing canonical general UV/EFT1 beta solvers and conditional
one-scale Matchete C5 adapters, NOT a certified fully sequential EFT.
B/C(-1) use the existing separate mixed-scalar candidate C5 adapter.

Run at the project root:
  python -X utf8 -m scripts.RunT3GeneralDualBenchmarks --self-test
  python -X utf8 -m scripts.RunT3GeneralDualBenchmarks --run-id thesis51_dual_v1 --preflight
  python -u -X utf8 -m scripts.RunT3GeneralDualBenchmarks --run-id thesis51_dual_v1
  python -u -X utf8 -m scripts.RunT3GeneralDualBenchmarks --run-id thesis51_dual_v1 --resume
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import traceback
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GENERAL_FITS = (ROOT / 'output' / 'studies' / 'general_scalar_extended'
                / 'thesis_general_all16' / 'thesis16_general_c5_full_v1')
OUTPUT = ROOT / 'output' / 'studies' / 'general_scalar_extended' / 'shared_general_dual'
MODELS = ('A_-4', 'A_-2', 'A_+0', 'B_-3', 'B_-1', 'B_+1',
          'C_-3', 'C_-1', 'C_+1', 'D_-2', 'D_+0', 'D_+2',
          'E_-4', 'E_-2', 'E_+0', 'E_+2')
SOURCE_MODELS = ('D_-2', 'A_-4')
MIXED = frozenset(('B_-1', 'C_-1'))
COMPARISON_PAIRS = (('A_-2', 'D_-2'), ('B_-3', 'C_-3'), ('B_-1', 'C_-1'))
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
       'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
FORMAT = 't3-general-dual-fixed-yukawa-conditional-c5-v1'


def model_tag(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


def json_read(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f'Missing required input: {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise ValueError(f'Expected JSON object: {path}')
    return obj


def json_save(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def file_sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def compact_matrix(raw: dict) -> np.ndarray:
    out = (np.asarray(raw['real'], dtype=float)
           + 1j * np.asarray(raw['imag'], dtype=float))
    if out.shape != (3, 3) or not np.isfinite(out).all():
        raise ValueError('Invalid complex 3x3 matrix')
    return out


def entry_matrix(raw: list[list[dict]]) -> np.ndarray:
    out = np.asarray([[complex(v['real'], v['imag']) for v in row]
                      for row in raw], dtype=complex)
    if out.shape != (3, 3) or not np.isfinite(out).all():
        raise ValueError('Invalid entry-wise complex 3x3 matrix')
    return out


def matrix_json(a: np.ndarray) -> dict:
    arr = np.asarray(a, dtype=complex)
    if arr.shape != (3, 3) or not np.isfinite(arr).all():
        raise ValueError('Cannot save invalid 3x3 matrix')
    return {'real': arr.real.tolist(), 'imag': arr.imag.tolist()}


def source_files(model: str) -> dict[str, Path]:
    base = GENERAL_FITS / model_tag(model)
    return {'fitted': base / 'fitted_config.json',
            'result': base / 'result.json',
            'overrides': base / 'general_uv_overrides.json'}


def load_source(model: str) -> dict:
    paths = source_files(model)
    cfg, result = json_read(paths['fitted']), json_read(paths['result'])
    overrides = json_read(paths['overrides'])
    if result.get('model') != model or not result.get('general_rges_used') or not result.get('general_matchete_matching_used'):
        raise ValueError(f'{model}: selected fit is not a general RGE+matching result')
    expected = {'A_-4': (1, 3, 2, -4), 'D_-2': (3, 1, 2, -2)}[model]
    rep = cfg['representation']
    actual = tuple(int(rep[k]) for k in ('d_s1', 'd_s2', 'd_f', 'alpha'))
    if actual != expected or rep.get('shared_scalar', True):
        raise ValueError(f'{model}: unexpected source representation: {rep}')
    if not result.get('chi2_final', math.inf) <= 5.86:
        raise ValueError(f'{model}: source fit does not pass original chi2 screen')
    if float(cfg['base_state']['mu_gev']) != 1e7:
        raise ValueError('Source must be defined at 1e7 GeV')
    arrays = {name: compact_matrix(cfg['base_state']['ordinary'][name])
              for name in ('y1', 'y2')}
    if any(np.max(np.abs(v.imag)) > 1e-12 for v in arrays.values()):
        raise ValueError('Existing general-fitting adapters only accept 18 real Yukawa entries')
    params = {f'{name}_{i}{j}': float(arrays[name][i, j].real)
              for name in ('y1', 'y2') for i in range(3) for j in range(3)}
    extra = {}
    for key in ('y1Conj', 'y2Conj'):
        if key in overrides:
            extra[key] = entry_matrix(overrides[key])
    return {'model': model, 'files': paths, 'config': cfg, 'result': result,
            'params': params, 'additional_yukawas': extra,
            'chi2': float(result['chi2_final']),
            'hashes': {name: file_sha(p) for name, p in paths.items()}}



def source_preview(source_data: dict, dest: Path) -> None:
    """Immediate independently fitted source-matrix comparison (NOT a sweep)."""
    if not all(m in source_data for m in SOURCE_MODELS):
        return
    d = source_data['D_-2']['result']['numerical_details']['c5_at_low_scale']
    a = source_data['A_-4']['result']['numerical_details']['c5_at_low_scale']
    dm, am = compact_matrix(d), compact_matrix(a)
    k, eps = proportionality(dm, am)
    payload = {
        'status': 'SavedIndependentFittedSourcesOnly_NotSharedYukawaSweep',
        'scale_gev': 100,
        'D_-2_chi2': source_data['D_-2']['chi2'],
        'A_-4_chi2': source_data['A_-4']['chi2'],
        'D_-2_C5_low': matrix_json(dm),
        'A_-4_C5_low': matrix_json(am),
        'k_Aminus4_to_Dminus2': {'real': float(k.real), 'imag': float(k.imag)},
        'epsilon_Aminus4_to_Dminus2': eps,
        'D_-2_C5_norm': float(np.linalg.norm(dm)),
        'A_-4_C5_norm': float(np.linalg.norm(am)),
        'caution': ('These matrices use different fitted Yukawa inputs, so this '
                    'does not test representation dependence at common Yukawas.')}
    json_save(dest / 'source_independent_fit_preview.json', payload)


def coefficients(detail: dict, mixed: bool) -> tuple[np.ndarray, np.ndarray]:
    a = 'C5_muS_gev_inv' if mixed else 'c5_at_scalar_threshold'
    b = 'C5_low_gev_inv' if mixed else 'c5_at_low_scale'
    scalar = compact_matrix(detail[a])
    low = compact_matrix(detail[b])
    for label, mat in (('scalar', scalar), ('low', low)):
        if not np.allclose(mat, mat.T, rtol=1e-7, atol=1e-23):
            raise ValueError(f'Asymmetric {label} C5 matrix')
    return scalar, low


def proportionality(reference: np.ndarray, candidate: np.ndarray) -> tuple[complex, float]:
    """Best complex scalar multiplier and Frobenius residual."""
    r, c = np.asarray(reference, dtype=complex), np.asarray(candidate, dtype=complex)
    if r.shape != (3, 3) or c.shape != (3, 3) or not (np.isfinite(r).all() and np.isfinite(c).all()):
        raise ValueError('Invalid C5 matrices')
    rn, cn = np.linalg.norm(r), np.linalg.norm(c)
    if rn < 1e-30 or cn < 1e-30:
        raise ValueError('Cannot compare vanishing C5 matrix')
    k = np.vdot(r, c) / np.vdot(r, r)
    eps = float(np.linalg.norm(c - k * r) / cn)
    return complex(k), eps


def eval_source_model(source: dict, model: str, args: argparse.Namespace,
                      *, execute: bool) -> tuple[dict, dict | None]:
    """Use the existing general evaluators. Never substitute restricted C5."""
    from scripts import RunT3GeneralAll16FromReference as general
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET

    target = OscillationFitTarget.from_json(TARGET)
    params = source['params']
    extra_source = source['additional_yukawas']
    if model not in MIXED:
        uv, eft, provenance, kernel, c5_path, method = general.nonmixed_inputs(
            model, args.canonical_run_id, generate_general_c5=False)
        ref_path = general.ref_path(model)
        ref = json_read(ref_path)
        fixed = general.fixed_overrides(
            ref, uv, quartic_fraction=args.quartic_fraction,
            mu2_fraction=args.mu2_fraction,
            extra_yukawa_fraction=args.extra_yukawa_fraction)
        extra_policy = {}
        for key in ('y1Conj', 'y2Conj'):
            if key not in fixed:
                continue
            if key in extra_source:
                values = extra_source[key]
                extra_policy[key] = 'exact_saved_source_matrix'
            else:
                parent = 'y1' if key == 'y1Conj' else 'y2'
                values = np.asarray([[params[f'{parent}_{i}{j}'] for j in range(3)]
                                     for i in range(3)], dtype=float) * args.extra_yukawa_fraction
                extra_policy[key] = 'assumed_fraction_times_shared_' + parent
            fixed[key] = general.entry_matrix(values)
        row = {'model': model, 'source_model': source['model'],
               'status': 'Ready', 'method': method,
               'canonical_run_id': provenance.get('run_id'),
               'general_beta_physical_mapping_certified': False,
               'target_reference_sha256': file_sha(ref_path),
               'general_c5_kernel_sha256': file_sha(c5_path),
               'extra_yukawas': extra_policy,
               'additional_quartic_coordinates': sorted(
                   k for k in fixed if k.startswith('lambdaT3General')),
               'additional_bilinear_coordinates': sorted(
                   k for k in fixed if k.startswith('muT3General'))}
        if not execute:
            return row, None
        evaluator = general.NonMixedEvaluator(model, uv, eft, kernel, fixed, target)
        fit, detail = evaluator.evaluate(params, detail=True)
        scalar, low = coefficients(detail, mixed=False)
    else:
        # Exact same conditional B/C(-1) adapter as the existing general fit.
        from scripts import RunT3BCGeneralFromReference as bc_ref
        from scripts import RunT3GeneralScalarBCExploratoryFit as bc
        from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE

        bc_source = bc_ref.build_reference_source(model)
        beta = bc_ref.ensure_betas(model, bc_source, export=False,
                                  wolframscript=args.wolframscript)
        c5_path = DEFAULT_SOURCE / (
            {'B_-1': 'B_m1_final', 'C_-1': 'C_m1_final'}[model]
            + '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
        if not c5_path.is_file():
            raise FileNotFoundError(f'Missing B/C conditional matching: {c5_path}')
        bundle = bc.beta_bundle(model, bc_source, bc_source['registry'])
        row = {'model': model, 'source_model': source['model'],
               'status': 'Ready', 'method': 'BC_mixed_scalar_one_scale_candidate',
               'general_beta_physical_mapping_certified': False,
               'target_reference_sha256': file_sha(bc_source['cfg_path']),
               'general_c5_kernel_sha256': file_sha(c5_path),
               'general_uv_beta_sha256': file_sha(beta['UV']),
               'general_eft1_beta_sha256': file_sha(beta['EFT1']),
               'extra_yukawas': {},
               'additional_quartic_coordinates': list(bc_source['new']) + list(bc_source['reused']),
               'additional_bilinear_coordinates': ['muT3GeneralS1S2']}
        if not execute:
            return row, None
        compiled = bc.compile_source(c5_path)
        evaluator = bc.CandidateEvaluator(
            'reference_scalar_general', model, bc_source, bc_source['registry'],
            bundle, compiled, target, quartic_fraction=args.quartic_fraction,
            mu2_fraction=args.mu2_fraction)
        fit, detail = evaluator.evaluate(params, detail=True)
        scalar, low = coefficients(detail, mixed=True)
    pred = dict(zip(fit.observable_names, map(float, fit.prediction)))
    if set(OBS) - set(pred):
        raise ValueError(f'Missing oscillation observable(s): {sorted(set(OBS)-set(pred))}')
    values = {'model': model, 'source_model': source['model'],
              'C5_at_scalar_gev_inv': matrix_json(scalar),
              'C5_at_low_gev_inv': matrix_json(low),
              'prediction': pred, 'chi2': float(fit.chi2),
              'matching_diagnostics': detail,
              'physical_sequential_eft_certified': False}
    # Cross-check the original fitted source point, where the source is also
    # the target model. Do not silently mislabel a non-reproduced baseline.
    if model == source['model']:
        prior = source['result']['numerical_details']['c5_at_low_scale']
        prior_mat = compact_matrix(prior)
        diff = float(np.linalg.norm(low - prior_mat) / np.linalg.norm(prior_mat))
        values['original_source_C5_relative_difference'] = diff
        row['source_C5_reproduction_relative_difference'] = diff
        row['source_C5_reproduction_check'] = 'Passed' if diff < 1e-5 else 'Mismatch'
        if diff >= 1e-5:
            raise ValueError(f'{model} source C5 was not reproduced (relative difference={diff:.5g})')
    singular = np.linalg.svd(low, compute_uv=False)
    row.update({'status': 'CompletedConditional', 'chi2': float(fit.chi2),
                'fit_screen_pass_chi2_le_5p86': bool(float(fit.chi2) <= 5.86),
                'C5_muS_norm_gev_inv': float(np.linalg.norm(scalar)),
                'C5_low_norm_gev_inv': float(np.linalg.norm(low)),
                'neutrino_masses_low_ev': sorted(float(z*(246.22**2)/2*1e9) for z in singular),
                'prediction': pred})
    return row, values


def gather_results(dest: Path, sources: tuple[str, ...], models: tuple[str, ...]) -> tuple[list[dict], list[dict], list[dict]]:
    """Assemble per-model results and conditional source-relative comparisons."""
    base_rows, matrix_lookup = [], {}
    for source in sources:
        for model in models:
            st = dest / model_tag(source) / model_tag(model) / 'status.json'
            if not st.is_file():
                continue
            row = json_read(st)
            base_rows.append(row)
            if row['status'] != 'CompletedConditional':
                continue
            value = json_read(st.with_name('result.json'))
            matrix_lookup[(source, model)] = compact_matrix(value['C5_at_low_gev_inv'])
    for row in base_rows:
        src, model = row['source_model'], row['model']
        self_ref = matrix_lookup.get((src, src))
        candidate = matrix_lookup.get((src, model))
        if self_ref is None or candidate is None:
            row['k_to_source_real'] = row['k_to_source_imag'] = row['epsilon_to_source'] = None
        else:
            k, eps = proportionality(self_ref, candidate)
            row.update({'k_to_source_real': float(k.real),
                        'k_to_source_imag': float(k.imag),
                        'epsilon_to_source': eps})
    pairs = []
    for source in sources:
        for left, right in COMPARISON_PAIRS:
            if left not in models or right not in models:
                continue
            a, b = matrix_lookup.get((source, left)), matrix_lookup.get((source, right))
            if a is None or b is None:
                pairs.append({'benchmark': source, 'model_a': left, 'model_b': right,
                              'status': 'Unavailable'})
                continue
            k, eps = proportionality(a, b)
            pairs.append({'benchmark': source, 'model_a': left, 'model_b': right,
                          'status': 'CompletedConditional', 'k_real': float(k.real),
                          'k_imag': float(k.imag), 'epsilon_pair': eps,
                          'C5_norm_ratio_b_over_a': float(np.linalg.norm(b)/np.linalg.norm(a))})
    crosses = []
    if len(sources) == 2:
        for model in models:
            a = matrix_lookup.get((sources[0], model))
            b = matrix_lookup.get((sources[1], model))
            if a is None or b is None:
                crosses.append({'model': model, 'status': 'Unavailable'})
                continue
            k, eps = proportionality(a, b)
            crosses.append({'model': model, 'status': 'CompletedConditional',
                            'k_AoverD_real': float(k.real), 'k_AoverD_imag': float(k.imag),
                            'epsilon_AoverD_at_same_model': eps})
    return base_rows, pairs, crosses


def write_csv(dest: Path, rows: list[dict], pairs: list[dict], cross: list[dict]) -> None:
    cols = ('source_model', 'model', 'status', 'method', 'chi2',
            'fit_screen_pass_chi2_le_5p86', 'k_to_source_real', 'k_to_source_imag',
            'epsilon_to_source', 'C5_low_norm_gev_inv', 'C5_muS_norm_gev_inv',
            'm1_ev', 'm2_ev', 'm3_ev', *OBS, 'source_C5_reproduction_check', 'message')
    with (dest / 'summary.csv').open('w', newline='', encoding='utf-8') as fp:
        writer = csv.DictWriter(fp, fieldnames=cols)
        writer.writeheader()
        for r in rows:
            data = {x: r.get(x) for x in cols}
            data.update(zip(('m1_ev', 'm2_ev', 'm3_ev'),
                            r.get('neutrino_masses_low_ev', [None] * 3)))
            data.update(r.get('prediction', {}))
            writer.writerow(data)
    with (dest / 'six_pair_comparison.csv').open('w', newline='', encoding='utf-8') as fp:
        cols = ('benchmark', 'model_a', 'model_b', 'status', 'k_real', 'k_imag',
                'epsilon_pair', 'C5_norm_ratio_b_over_a')
        writer = csv.DictWriter(fp, fieldnames=cols)
        writer.writeheader()
        for x in pairs:
            writer.writerow({k: x.get(k) for k in cols})
    with (dest / 'same_model_texture_comparison.csv').open('w', newline='', encoding='utf-8') as fp:
        cols = ('model', 'status', 'k_AoverD_real', 'k_AoverD_imag',
                'epsilon_AoverD_at_same_model')
        writer = csv.DictWriter(fp, fieldnames=cols)
        writer.writeheader()
        for x in cross:
            writer.writerow({k: x.get(k) for k in cols})


def make_plots(dest: Path, rows: list[dict], pairs: list[dict], crosses: list[dict],
               sources: tuple[str, ...], models: tuple[str, ...]) -> None:
    if not any(r['status'] == 'CompletedConditional' for r in rows):
        return
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return
    look = {(r['source_model'], r['model']): r for r in rows}
    xs = np.arange(len(models))
    labels = [m.replace('_', '(') + ')' for m in models]
    for key, ylabel, title, name, logy in (
        ('chi2', r'Conditional $\chi^2$', 'Two fitted Yukawa textures: no refitting',
         'fig_dual_fixed_yukawa_chi2', True),
        ('epsilon_to_source', r'$\epsilon_{C_5}$ to own source',
         'General C5 shape changes by representation',
         'fig_dual_C5_epsilon', False),
        ('sin2_theta23', r'$\sin^2\theta_{23}$',
         'Mixing-angle response to representation',
         'fig_dual_theta23', False),
    ):
        fig, ax = plt.subplots(figsize=(max(9, 0.64*len(models)), 4.6))
        any_line = False
        for source in sources:
            ys = []
            for m in models:
                r = look.get((source, m), {})
                value = r.get(key, r.get('prediction', {}).get(key))
                ys.append(max(float(value), 1e-14) if logy and value is not None
                          else (float(value) if value is not None else np.nan))
            if any(math.isfinite(v) for v in ys):
                ax.plot(xs, ys, marker='o', label=f'{source} fixed Yukawas')
                any_line = True
        if not any_line:
            plt.close(fig)
            continue
        if logy:
            ax.set_yscale('log')
            ax.axhline(5.86, color='gray', linewidth=1, linestyle='--', label='Fit screen 5.86')
        ax.set_xticks(xs, labels, rotation=60, ha='right')
        ax.set_ylabel(ylabel)
        ax.set_title(title + ' (conditional matching)')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.15)
        fig.tight_layout()
        fig.savefig(dest / (name + '.png'), dpi=190)
        fig.savefig(dest / (name + '.pdf'))
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    lbls = [f'{p["benchmark"]}: {p["model_a"]}/{p["model_b"]}' for p in pairs]
    vals = [p.get('epsilon_pair', np.nan) for p in pairs]
    if lbls and any(math.isfinite(v) for v in vals):
        ax.barh(np.arange(len(lbls)), vals)
        ax.set_yticks(np.arange(len(lbls)), lbls)
        ax.set_xlabel(r'Best-fit C5 pairwise residual $\epsilon$')
        ax.set_title('Original thesis 5.1 representation pairs, two textures')
        fig.tight_layout()
        fig.savefig(dest / 'fig_six_representation_pairs.png', dpi=190)
        fig.savefig(dest / 'fig_six_representation_pairs.pdf')
    plt.close(fig)


def summarise(dest: Path, source_data: dict, models: tuple[str, ...],
              *, preflight: bool, plot: bool = False) -> None:
    sources = tuple(source_data)
    rows, pairs, crosses = gather_results(dest, sources, models)
    status_counts = {key: sum(r.get('status') == key for r in rows)
                     for key in ('Ready', 'CompletedConditional', 'Blocked')}
    result = {'format': FORMAT, 'source_models': sources,
              'source_fitted_chi2': {name: source_data[name]['chi2'] for name in sources},
              'source_input_sha256': {name: source_data[name]['hashes'] for name in sources},
              'fixed_yukawa_policy': 'source-specific fitted y1/y2 held identical for every target model',
              'extra_yukawa_policy': ('copy exact fitted-source auxiliary Yukawa when source permits it; '
                                      'otherwise use 0.05 times corresponding shared ordinary Yukawa '
                                      'when allowed by target symmetry'),
              'matching_scope': ('full-general UV/EFT1 renormalisable running, one-scale conditional C5; '
                                 'B/C(-1) use separate conditional mixed-scalar mass adapter'),
              'physical_sequential_eft_certified': False,
              'status_counts': status_counts,
              'rows': rows, 'six_representation_pairs': pairs,
              'same_model_texture_comparison': crosses}
    json_save(dest / 'summary.json', result)
    write_csv(dest, rows, pairs, crosses)
    if plot and not preflight:
        make_plots(dest, rows, pairs, crosses, sources, models)
    lines = ['# General T3 dual benchmark — thesis Section 5.1', '',
             'Benchmarks: D(-2) and A(-4) independently fitted full-general Yukawa matrices.',
             'For each benchmark the same 3x3 real y1,y2 is used in every target model.',
             'Other fixed UV inputs come from target reference configurations (not necessarily identical across models).',
             'Added quartic boundary = 0.1 * |lambdaT3|; scalar bilinear boundary = 0.01 * min(mS1Sq,mS2Sq).',
             'Additional Yukawas: exact source y2Conj copied for D(-2),',
             'otherwise 0.05 * shared y1/y2 where that extra channel exists.',
             'The exact-source and assumed-fraction cases are not equivalent controls.', '',
             f'Completed conditional evaluations: {status_counts["CompletedConditional"]}/{2*len(models)}.',
             f'Blocked: {status_counts["Blocked"]}; preflight-ready: {status_counts["Ready"]}.', '',
             'summary.csv: C5 norm, source-relative residual, masses and oscillation observables.',
             'six_pair_comparison.csv: A(-2)/D(-2), B(-3)/C(-3), B(-1)/C(-1) under both Yukawa benchmarks.',
             'same_model_texture_comparison.csv: D(-2) vs A(-4) Yukawas in each fixed target model.',
             'See per-model result.json and status.json for the full 3x3 C5 matrices and provenance.', '',
             'STATUS: conditional one-scale general matching, not physically certified sequential EFT.',
             'The physical Matchete-RGBeta coupling dictionary is uncertified.',
             'B/C(-1) uses a separate matching evaluator; cross-adapter normalisation is uncertified.',
             'D(0) can be blocked by an unsupported negative real self-bilinear rotation.',
             'A large fixed-Yukawa chi2 does not exclude a model after independent Yukawa refitting.', '']
    (dest / 'README.md').write_text('\n'.join(lines), encoding='utf-8')


def self_test() -> None:
    z = np.array([[1, 2j, 0], [2j, 3, 1], [0, 1, 4]], dtype=complex) * 1e-15
    k, eps = proportionality(z, (2+1j)*z)
    assert abs(k - (2+1j)) < 1e-12 and eps < 1e-12
    assert np.allclose(compact_matrix(matrix_json(z)), z)
    assert entry_matrix([[{'real': float(i==j), 'imag': 0} for j in range(3)]
                         for i in range(3)]).shape == (3, 3)
    assert len(MODELS) == 16 and len(COMPARISON_PAIRS) == 3
    print('Self-test passed: complex C5 scaling, matrix I/O, model coverage')


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', help='Unique output directory name')
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--preflight', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--models', default=','.join(MODELS),
                   help='Comma-separated subset of models; default all 16')
    p.add_argument('--canonical-run-id', default='full16_general_v1')
    p.add_argument('--quartic-fraction', type=float, default=0.1)
    p.add_argument('--mu2-fraction', type=float, default=0.01)
    p.add_argument('--extra-yukawa-fraction', type=float, default=0.05)
    p.add_argument('--wolframscript', default='wolframscript')
    args = p.parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    if not args.run_id or not all(ch.isalnum() or ch in '-_' for ch in args.run_id):
        p.error('--run-id must contain only letters, numbers, underscore or hyphen')
    if not (0 < args.quartic_fraction <= 1 and 0 <= args.mu2_fraction <= 0.2
            and 0 <= args.extra_yukawa_fraction <= 1):
        p.error('Invalid general-scalar or extra-Yukawa fraction')
    selected = [m.strip() for m in args.models.split(',') if m.strip()]
    if not selected or len(selected) != len(set(selected)) or set(selected) - set(MODELS):
        p.error('--models must select unique valid model identifiers')
    selected = tuple(m for m in MODELS if m in selected)
    source_data = {name: load_source(name) for name in SOURCE_MODELS}
    dest = OUTPUT / args.run_id
    settings_path = dest / 'settings.json'
    signature = {'format': FORMAT, 'source_input_sha256': {m: d['hashes'] for m, d in source_data.items()},
                 'source_models': SOURCE_MODELS, 'models': selected,
                 'canonical_run_id': args.canonical_run_id,
                 'quartic_fraction': args.quartic_fraction,
                 'mu2_fraction': args.mu2_fraction,
                 'extra_yukawa_fraction': args.extra_yukawa_fraction}
    if settings_path.is_file():
        if json_read(settings_path) != signature:
            p.error('Existing run settings changed: use a new --run-id')
        if not args.resume and not args.preflight:
            any_done = any((dest / model_tag(s) / model_tag(m) / 'result.json').is_file()
                           for s in SOURCE_MODELS for m in selected)
            if any_done:
                p.error('Existing results: use --resume or a new --run-id')
    else:
        if dest.exists() and any(dest.iterdir()):
            p.error('Existing nonempty directory without validated settings.json')
        json_save(settings_path, signature)
    source_preview(source_data, dest)
    for src in SOURCE_MODELS:
        source = source_data[src]
        for i, model in enumerate(selected, 1):
            st_dir = dest / model_tag(src) / model_tag(model)
            status_path = st_dir / 'status.json'
            result_path = st_dir / 'result.json'
            print(f'[{src} {i}/{len(selected)}] {model}', flush=True)
            if args.preflight and result_path.is_file():
                # A preflight must never replace a completed numerical status.
                print('    Existing numerical result preserved', flush=True)
                continue
            if args.resume and not args.preflight and result_path.is_file():
                try:
                    row = json_read(status_path)
                    stored = json_read(result_path)
                    if row.get('status') != 'CompletedConditional' or stored.get('source_model') != src or stored.get('model') != model:
                        raise ValueError('Resume identity/status mismatch')
                    print('    Reused completed conditional result', flush=True)
                    continue
                except Exception as exc:
                    print(f'    Invalid previous result: {exc}', flush=True)
                    # Never overwrite unverified prior outputs.
                    raise
            try:
                row, result = eval_source_model(source, model, args,
                                                execute=not args.preflight)
                if result is not None:
                    json_save(result_path, result)
                if row.get('chi2') is not None:
                    print(f"    conditional chi2={row['chi2']:.9g}", flush=True)
                else:
                    print('    Preflight Ready', flush=True)
            except Exception as exc:
                row = {'model': model, 'source_model': src, 'status': 'Blocked',
                       'message': f'{type(exc).__name__}: {exc}',
                       'traceback': traceback.format_exc(limit=6)}
                print('    BLOCKED:', row['message'], flush=True)
            json_save(status_path, row)
            summarise(dest, source_data, selected, preflight=args.preflight)
    summarise(dest, source_data, selected, preflight=args.preflight, plot=not args.preflight)
    print(f'Results: {dest / "summary.csv"}', flush=True)
    return 0 if all(row['status'] in ('Ready', 'CompletedConditional')
                    for row in gather_results(dest, SOURCE_MODELS, selected)[0]) else 2


if __name__ == '__main__':
    raise SystemExit(main())
