"""Fixed-A(-2)-Yukawa comparison across the 16 general ordinary T3 models.

Source: the *general* A_-2 fitted_config.json from thesis16_general_c5_full_v1.
No optimisation or updates to pre-existing outputs. Each model uses its own
symmetry-allowed general scalar/bilinear basis and non-Yukawa benchmark values.
Canonical general RGEs and the pre-existing candidate general C5 evaluators
are used; the result is NOT certified full sequential EFT matching.

Run from project root:
    python -X utf8 -m scripts.RunT3GeneralSharedA2 --self-test
    python -X utf8 -m scripts.RunT3GeneralSharedA2 --run-id thesis51_shared_a2_v1 --preflight
    python -u -X utf8 -m scripts.RunT3GeneralSharedA2 --run-id thesis51_shared_a2_v1
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

from scripts import RunT3GeneralAll16FromReference as general

ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT / 'output' / 'studies' / 'general_scalar_extended'
          / 'thesis_general_all16' / 'thesis16_general_c5_full_v1'
          / 'A_m2' / 'fitted_config.json')
SOURCE_RESULT = SOURCE.with_name('result.json')
SOURCE_OVERRIDES = SOURCE.with_name('general_uv_overrides.json')
OUTPUT = ROOT / 'output' / 'studies' / 'general_scalar_extended' / 'shared_general_A_m2'
FORMAT = 't3-general-shared-a2-all16-conditional-v1'
MODELS = general.MODELS
MIXED = general.MIXED
YUKAWA_NAMES = general.NAMES
OBSERVABLES = general.OBS


def save(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def shared_parameters(source: dict) -> dict[str, float]:
    ordinary = source['base_state']['ordinary']
    result: dict[str, float] = {}
    for name in ('y1', 'y2'):
        matrix = general.matrix(ordinary[name])
        if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
            raise ValueError(f'Invalid A_-2 {name} matrix')
        if np.max(np.abs(matrix.imag)) > 1e-13:
            raise ValueError(f'{name} is complex; existing general pilot assumes real fitted matrices')
        for i in range(3):
            for j in range(3):
                result[f'{name}_{i}{j}'] = float(matrix[i, j].real)
    return result


def make_matrix(spec: dict) -> np.ndarray:
    """Accept saved {real:[[...]],imag:[[...]]} complex-matrix format."""
    return np.asarray(spec['real'], dtype=float) + 1j * np.asarray(spec['imag'], dtype=float)


def c5_matrices(details: dict, *, mixed: bool) -> tuple[np.ndarray, np.ndarray]:
    scalar_key = 'C5_muS_gev_inv' if mixed else 'c5_at_scalar_threshold'
    low_key = 'C5_low_gev_inv' if mixed else 'c5_at_low_scale'
    scalar, low = make_matrix(details[scalar_key]), make_matrix(details[low_key])
    for name, mat in (('scalar', scalar), ('low', low)):
        if mat.shape != (3, 3) or not np.isfinite(mat).all():
            raise ValueError(f'Invalid {name} C5 matrix')
        if not np.allclose(mat, mat.T, rtol=1e-8, atol=1e-20):
            raise ValueError(f'{name} C5 matrix not symmetric')
    return scalar, low


def matrix_json(mat: np.ndarray) -> dict:
    return general.matrix_json(mat)


def proportionality(ref: np.ndarray, test: np.ndarray) -> tuple[complex, float]:
    denom = np.vdot(ref, ref)
    if abs(denom) < 1e-50 or np.linalg.norm(test) < 1e-30:
        raise ValueError('C5 too small for meaningful relative proportionality')
    k = np.vdot(ref, test) / denom
    residual = float(np.linalg.norm(test - k * ref) / np.linalg.norm(test))
    return complex(k), residual


def masses(c5: np.ndarray, *, vev: float = 246.22) -> list[float]:
    """Takagi singular values, with m_nu = -v^2*C5/2 and GeV -> eV."""
    sv = np.linalg.svd(c5, compute_uv=False)
    return sorted(float(x * vev**2 * 0.5 * 1e9) for x in sv)


def status_model(model: str, params: dict, args: argparse.Namespace, *, execute: bool) -> tuple[dict, dict | None]:
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET

    target = OscillationFitTarget.from_json(TARGET)
    if model not in MIXED:
        uv, eft, provenance, kernel, match_path, method = general.nonmixed_inputs(
            model, args.canonical_run_id, generate_general_c5=False)
        reference_path = general.ref_path(model)
        reference = general.read(reference_path)
        fixed = general.fixed_overrides(
            reference, uv,
            quartic_fraction=args.quartic_fraction,
            mu2_fraction=args.mu2_fraction,
            extra_yukawa_fraction=args.extra_yukawa_fraction)
        # The *actual* saved A(-2) general benchmark also contains a fixed
        # y2Conj matrix, which was set from its ORIGINAL reference input
        # before fitting y1/y2. It is NOT 0.05 times the fitted y2 matrix.
        # Transfer that original A(-2) y2Conj unchanged wherever allowed.
        # For alpha=0 the allowed distinct interaction is y1Conj, which
        # A(-2) does not possess. Use the documented candidate convention
        # y1Conj = fraction * shared y1, and label it a new choice.
        if 'y2Conj' in fixed:
            source_y2conj = general.read(SOURCE_OVERRIDES)['y2Conj']
            z = np.array([[complex(v['real'], v['imag']) for v in row]
                          for row in source_y2conj], dtype=complex)
            if z.shape != (3, 3) or not np.isfinite(z).all():
                raise ValueError('A(-2) general y2Conj invalid')
            fixed['y2Conj'] = source_y2conj
        if 'y1Conj' in fixed:
            arr = np.array([[params[f'y1_{i}{j}'] for j in range(3)]
                            for i in range(3)], dtype=float)
            fixed['y1Conj'] = general.entry_matrix(arr * args.extra_yukawa_fraction)
        row = {'model': model, 'method': method, 'scope': 'conditional_one_scale',
               'provenance_run': provenance.get('run_id'),
               'general_matchete_rge_dictionary_certified': False,
               'model_reference_sha256': sha(reference_path),
               'general_c5_kernel_sha256': sha(match_path),
               'extra_yukawas': sorted(x for x in fixed if x in ('y1Conj', 'y2Conj')),
               'y2Conj_source': ('saved_general_A2' if 'y2Conj' in fixed else None),
               'y1Conj_source': ('assumed_fraction_shared_y1' if 'y1Conj' in fixed else None),
               'extra_quartics': sum(x.startswith('lambdaT3General') for x in fixed),
               'extra_bilinears': sum(x.startswith('muT3General') for x in fixed)}
        if not execute:
            return dict(row, status='Ready'), None
        evaluator = general.NonMixedEvaluator(model, uv, eft, kernel, fixed, target)
        fit, detail = evaluator.evaluate(params, detail=True)
        scalar, low = c5_matrices(detail, mixed=False)
    else:
        # B/C(-1) have a mixed scalar mass term. Reuse the pre-existing
        # specialised one-scale general-quartic+mixed-mass numerical adapter.
        # This is not the 14-model generic Matchete adapter.
        from scripts import RunT3BCGeneralFromReference as bc_ref
        from scripts import RunT3GeneralScalarBCExploratoryFit as bc
        from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
        source = bc_ref.build_reference_source(model)
        beta_paths = bc_ref.ensure_betas(model, source, export=False, wolframscript='wolframscript')
        coefficient = DEFAULT_SOURCE / (
            {'B_-1': 'B_m1_final', 'C_-1': 'C_m1_final'}[model]
            + '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
        if not coefficient.is_file():
            raise FileNotFoundError(f'Missing mixed-scalar general one-scale coefficient: {coefficient}')
        bundle = bc.beta_bundle(model, source, source['registry'])
        row = {'model': model, 'method': 'BC_mixed_scalar_conditional_one_scale',
               'scope': 'conditional_one_scale',
               'general_matchete_rge_dictionary_certified': False,
               'model_reference_sha256': sha(source['cfg_path']),
               'general_c5_kernel_sha256': sha(coefficient),
               'general_uv_beta_sha256': sha(beta_paths['UV']),
               'general_eft1_beta_sha256': sha(beta_paths['EFT1']),
               'extra_yukawas': [],
               'note': 'Separate B/C mixed-mass adapter; not cross-certified to generic Matchete matching.'}
        if not execute:
            return dict(row, status='Ready'), None
        compiled = bc.compile_source(coefficient)
        evaluator = bc.CandidateEvaluator(
            'reference_scalar_general', model, source, source['registry'],
            bundle, compiled, target,
            quartic_fraction=args.quartic_fraction, mu2_fraction=args.mu2_fraction)
        fit, detail = evaluator.evaluate(params, detail=True)
        scalar, low = c5_matrices(detail, mixed=True)
    predictions = dict(zip(fit.observable_names, (float(x) for x in fit.prediction)))
    row.update({'status': 'CompletedConditional', 'chi2': float(fit.chi2),
                'fit_screen_pass_chi2_le_5p86': bool(float(fit.chi2) <= 5.86),
                'C5_muS_frobenius_gev_inv': float(np.linalg.norm(scalar)),
                'C5_low_frobenius_gev_inv': float(np.linalg.norm(low)),
                'neutrino_masses_low_ev': masses(low),
                'prediction': predictions})
    detail_out = {'model': model, 'source_model': 'A_-2',
                  'C5_at_scalar_gev_inv': matrix_json(scalar),
                  'C5_at_low_gev_inv': matrix_json(low),
                  'prediction': predictions, 'chi2': row['chi2'],
                  'conditional_matching_details': detail,
                  'physical_sequential_eft_certified': False}
    return row, detail_out


def summarise(dest: Path, rows: list[dict], params: dict, a2_sha: str) -> None:
    completed = [r for r in rows if r['status'] == 'CompletedConditional']
    matrices: dict[str, np.ndarray] = {}
    for row in completed:
        path = dest / general.tag(row['model']) / 'result.json'
        obj = general.read(path)
        matrices[row['model']] = make_matrix(obj['C5_at_low_gev_inv'])
    reference = matrices.get('A_-2')
    for row in rows:
        if row['model'] in matrices and reference is not None:
            k, eps = proportionality(reference, matrices[row['model']])
            row['k_to_A_m2_real'] = float(k.real)
            row['k_to_A_m2_imag'] = float(k.imag)
            row['epsilon_to_A_m2'] = eps
        else:
            row['k_to_A_m2_real'] = None
            row['k_to_A_m2_imag'] = None
            row['epsilon_to_A_m2'] = None
    totals = {'completed_conditional': len(completed),
              'blocked': sum(r['status'] == 'Blocked' for r in rows),
              'ready_not_executed': sum(r['status'] == 'Ready' for r in rows)}
    summary = {'format': FORMAT, 'source_model': 'A_-2',
               'source_fitted_config_sha256': a2_sha,
               'source_general_overrides_sha256': sha(SOURCE_OVERRIDES),
               'source_original_fit_chi2': float(general.read(SOURCE_RESULT)['chi2_final']),
               'physical_sequential_eft_certified': False,
               'common_yukawa_parameters': params,
               'totals': totals, 'rows': rows}
    save(dest / 'summary.json', summary)
    cols = ('model', 'status', 'method', 'chi2', 'fit_screen_pass_chi2_le_5p86',
            'k_to_A_m2_real', 'k_to_A_m2_imag', 'epsilon_to_A_m2',
            'C5_muS_frobenius_gev_inv', 'C5_low_frobenius_gev_inv',
            'm1_ev', 'm2_ev', 'm3_ev', *OBSERVABLES, 'message')
    with (dest / 'summary.csv').open('w', newline='', encoding='utf-8') as handle:
        w = csv.DictWriter(handle, fieldnames=cols)
        w.writeheader()
        for row in rows:
            ms = row.get('neutrino_masses_low_ev', [None] * 3)
            vals = {x: row.get(x) for x in cols}
            vals.update(zip(('m1_ev', 'm2_ev', 'm3_ev'), ms))
            vals.update(row.get('prediction', {}))
            w.writerow(vals)
    lines = ['# General T3 shared-A(-2) Yukawa comparison', '',
             f"Models with conditional C5 results: {totals['completed_conditional']}/16.",
             f"Blocked: {totals['blocked']}; ready-only: {totals['ready_not_executed']}.", '',
             'The same non-diagonal real A(-2) *general-fit* y1/y2 are applied',
             'at the UV boundary to every target model. Other model-specific',
             'SM, fermion mass, scalar mass and scalar-sector inputs are preserved',
             'from their tracked reference configurations, plus fixed general',
             'extra-coupling fractions. Additional conjugate Yukawas exist',
             'only when the target symmetry allows and are proportional to',
             'the shared y1 for neutral S1 (this is an additional convention).',
             'For neutral S2, the exact saved A(-2) y2Conj benchmark matrix',
             'is copied. It is not 0.05 times newly fitted y2.', '',
             'NO OPTIMISATION. NO validated physical Matchete/RGBeta dictionary.',
             'NO complete fermions-first, scalars-second EFT Wilson matching.',
             'B/C(-1) use a distinct exploratory mixed-scalar matching adapter.',
             'Do not interpret a poor chi2 as representation exclusion.', '',
             f"Source original general A(-2) fit chi2 = {float(general.read(SOURCE_RESULT)['chi2_final']):.9g} (check fit-screen status).",
             '', 'See summary.csv and model result.json files.']
    (dest / 'README.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')




def make_figures(dest: Path, rows: list[dict]) -> None:
    """Optional, read-only thesis figures from finished conditional evaluations."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not installed: summary CSV/JSON still available', flush=True)
        return
    valid = [r for r in rows if r['status'] == 'CompletedConditional']
    if not valid:
        return
    labels = [r['model'].replace('_', '(') + ')' for r in valid]
    xs = np.arange(len(valid))
    fig, ax = plt.subplots(figsize=(max(8, len(valid)*0.58), 4.5))
    y = [max(float(r['chi2']), 1e-12) for r in valid]
    ax.bar(xs, y)
    ax.axhline(5.86, color='black', linestyle='--', linewidth=1,
               label='Campaign screening threshold: 5.86')
    ax.set_yscale('log')
    ax.set_xticks(xs, labels, rotation=65)
    ax.set_ylabel(r'Conditional five-observable $\chi^2$')
    ax.set_title('Shared general A(-2) Yukawas: no refitting')
    ax.legend(loc='upper left', fontsize=8)
    fig.tight_layout()
    fig.savefig(dest/'fig_shared_a2_general_chi2.png', dpi=190)
    fig.savefig(dest/'fig_shared_a2_general_chi2.pdf')
    plt.close(fig)
    valid_eps = [r for r in valid if r.get('epsilon_to_A_m2') is not None]
    if valid_eps:
        fig, ax = plt.subplots(figsize=(max(8,len(valid_eps)*0.58),4.5))
        xx = np.arange(len(valid_eps))
        ax.bar(xx, [r['epsilon_to_A_m2'] for r in valid_eps])
        ax.set_xticks(xx, [r['model'].replace('_','(')+')' for r in valid_eps],
                      rotation=65)
        ax.set_ylabel(r'$\epsilon_{C_5}$ relative to general A(-2)')
        ax.set_title('Low-scale C5 matrix proportionality; conditional matching')
        fig.tight_layout()
        fig.savefig(dest/'fig_shared_a2_general_C5_residual.png', dpi=190)
        fig.savefig(dest/'fig_shared_a2_general_C5_residual.pdf')
        plt.close(fig)

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', required=False)
    p.add_argument('--preflight', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--models', default=','.join(MODELS),
                   help='Comma-separated subset of the 16 ordinary T3 models')
    p.add_argument('--canonical-run-id', default='full16_general_v1')
    p.add_argument('--quartic-fraction', type=float, default=0.1)
    p.add_argument('--mu2-fraction', type=float, default=0.01)
    p.add_argument('--extra-yukawa-fraction', type=float, default=0.05)
    args = p.parse_args(argv)
    if args.self_test:
        eye = np.diag([1., 2., 3.]) * 1e-15
        k, eps = proportionality(eye, 2 * eye)
        assert abs(k - 2) < 1e-12 and eps < 1e-12
        assert len(masses(eye)) == 3
        print('Self-test passed')
        return 0
    if not args.run_id or not args.run_id.replace('_', '').replace('-', '').isalnum():
        p.error('--run-id must contain only letters, numbers, underscore or hyphen')
    if not 0 < args.quartic_fraction <= 1 or not 0 <= args.mu2_fraction <= 0.2:
        p.error('Invalid general-scalar quartic/mu2 fractions')
    if not 0 <= args.extra_yukawa_fraction <= 1:
        p.error('Invalid extra-Yukawa fraction')
    selected = [x.strip() for x in args.models.split(',') if x.strip()]
    if len(selected) != len(set(selected)) or not selected or set(selected) - set(MODELS):
        p.error('Invalid or repeated --models')
    selected = [m for m in MODELS if m in selected]
    if not SOURCE.is_file() or not SOURCE_RESULT.is_file() or not SOURCE_OVERRIDES.is_file():
        p.error(f'Missing general A(-2) benchmark: {SOURCE} / {SOURCE_RESULT}')
    source_result = general.read(SOURCE_RESULT)
    if source_result.get('model') != 'A_-2' or not source_result.get('general_rges_used'):
        p.error('A(-2) source is not recorded as a general-RGE fit')
    source_config = general.read(SOURCE)
    if source_config['representation'] != {'d_s1': 1, 'd_s2': 3, 'd_f': 2,
                                           'alpha': -2, 'shared_scalar': False}:
        p.error('Unexpected A(-2) representation in source config')
    params = shared_parameters(source_config)
    a2_sha = sha(SOURCE)
    dest = OUTPUT / args.run_id
    signature = {'format': FORMAT, 'source_sha256': a2_sha,
                 'source_result_sha256': sha(SOURCE_RESULT),
                 'source_general_overrides_sha256': sha(SOURCE_OVERRIDES),
                 'models': selected, 'canonical_run_id': args.canonical_run_id,
                 'quartic_fraction': args.quartic_fraction,
                 'mu2_fraction': args.mu2_fraction,
                 'extra_yukawa_fraction': args.extra_yukawa_fraction}
    settings = dest / 'settings.json'
    if settings.is_file():
        if general.read(settings) != signature:
            p.error('Existing campaign settings differ; use a new --run-id')
        if not args.resume and not args.preflight:
            # A completed preflight writes settings.json but no C5 results.
            # Allow exactly that transition; never silently overwrite physics.
            if any((dest / general.tag(m) / 'result.json').exists() for m in selected):
                p.error('Existing numerical results; specify --resume or a new --run-id')
    else:
        if dest.exists() and any(dest.iterdir()):
            p.error('Nonempty output directory without matching settings')
        save(settings, signature)
    rows: list[dict] = []
    for index, model in enumerate(selected, 1):
        print(f'[{index}/{len(selected)}] {model}', flush=True)
        result_file = dest / general.tag(model) / 'result.json'
        try:
            if args.resume and result_file.is_file() and not args.preflight:
                detail = general.read(result_file)
                if detail.get('model') != model or detail.get('source_model') != 'A_-2':
                    raise ValueError('Existing model result does not match source')
                row = general.read(dest / general.tag(model) / 'status.json')
                if row.get('status') != 'CompletedConditional':
                    raise ValueError('Existing resume status not completed')
            else:
                row, detail = status_model(model, params, args,
                                           execute=not args.preflight)
                if detail is not None:
                    save(result_file, detail)
            if row['status'] == 'CompletedConditional':
                print(f"    conditional chi2={row['chi2']:.7g}", flush=True)
            else:
                print(f"    {row['status']}", flush=True)
        except Exception as exc:
            row = {'model': model, 'status': 'Blocked',
                   'message': f'{type(exc).__name__}: {exc}',
                   'traceback': traceback.format_exc(limit=5)}
            print(f"    BLOCKED: {row['message']}", flush=True)
        save(dest / general.tag(model) / 'status.json', row)
        rows.append(row)
        summarise(dest, rows, params, a2_sha)
    if not args.preflight:
        make_figures(dest, rows)
    print('Results: ', dest / 'summary.csv', flush=True)
    return 0 if all(r['status'] in ('Ready', 'CompletedConditional') for r in rows) else 2


if __name__ == '__main__':
    raise SystemExit(main())
