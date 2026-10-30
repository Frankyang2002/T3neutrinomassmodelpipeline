"""Continue the six-model shared-A(-2) general-scalar benchmark for B(-3)/C(-3).

This read-only continuation resolves a measured factor-of-two convention mismatch
between the existing legacy fast C5 kernel and the existing symbolic reference.
The existing campaign and shared-comparison implementations are not modified.
The original parity check remains enforced at relative tolerance 3e-8.

The calculated predictions are CONDITIONAL: legacy common-scale matching,
incomplete general-scalar quadratic sector, no certified general sequential EFT.

PowerShell (repository root):
    python -X utf8 -m scripts.ContinueSixGeneralScalarLegacyC5

Inputs: output/studies/general_scalar_extended/bulk_48/bulk48_v1/
        output/studies/general_scalar_extended/six_benchmarks/shared_a2_v1/
Output: output/studies/general_scalar_extended/six_benchmarks/
        shared_a2_legacy_normalisation_v1/
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from scripts import CompareT3GeneralScalarSharedA2 as compare
from scripts import RunT3GeneralScalar48Campaign as bulk

ROOT = Path(__file__).resolve().parents[1]
CASE = 'large_uv_frozen_0p5'
MODELS = ('A_-2', 'D_-2', 'B_-1', 'C_-1', 'B_-3', 'C_-3')
LEGACY = ('B_-3', 'C_-3')
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12',
       'sin2_theta13', 'sin2_theta23')
FORMAT = 't3-conditional-shared-a2-legacy-parity-normalisation-v1'


def require(condition: bool, msg: str) -> None:
    if not condition:
        raise ValueError('Six-model legacy C5 continuation: ' + msg)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n',
                    encoding='utf-8')


def finite_array(mat: Any, label: str) -> np.ndarray:
    arr = np.asarray(mat, dtype=complex)
    require(arr.shape == (3, 3), label + ': expected 3x3')
    require(np.isfinite(arr.real).all() and np.isfinite(arr.imag).all(),
            label + ': nonfinite matrix')
    return arr


def encoded_matrix(matrix: np.ndarray) -> dict:
    arr = finite_array(matrix, 'C5')
    return {'real': arr.real.tolist(), 'imag': arr.imag.tolist()}


def evaluate_legacy(model: str, source: Path, point: dict,
                    extra_real: float) -> dict:
    from Numerical.orchestration import RunSmallLambdaRefit as reference
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
    from Numerical.running.SMWeinbergStage import evaluate_symbolic_c5
    from types import SimpleNamespace

    registry = bulk.read(bulk.REGISTRY)
    entry, new, reused = bulk.inventory(registry, model)
    require(not reused, model + ': unexpected reused quartics')
    history = bulk.historical(CASE, model, entry)
    require(not history['uv_on'], 'unexpected UV-running reference')
    beta_path = source / 'betas' / (bulk.safe(model) + '.json')
    payload = bulk.read(beta_path)
    kernel, c5_path, kind = bulk.matching(model)
    require(kind == 'legacy-common-scale' and isinstance(kernel, bulk.LegacyKernel),
            model + ': matching changed, expected LegacyKernel')

    # Source-verified: the symbolic legacy evaluator includes 1/2 in its
    # symmetric Yukawa contraction; the inherited fast evaluator does not.
    # The measured uncorrected ratio fast/reference is 2 for both selected models.
    # Correct only this kernel instance; do NOT alter any other evaluator.
    uncorrected = kernel.fun
    kernel.fun = lambda *args: 0.5 * uncorrected(*args)

    target = OscillationFitTarget.from_json(reference.TARGET)
    ev = bulk.Evaluator(CASE, model, entry, new, history, payload,
                        kernel, target, extra_real)
    x = np.asarray([point[k] for k in bulk.NAMES], dtype=float)
    require(x.shape == (18,) and np.isfinite(x).all(), 'invalid donor Yukawas')
    evaluation = ev.evaluate(x, 'shared_a_m2_legacy_reference_normalised')
    require(evaluation.get('screen') is True,
            model + ': failed physical screen: ' + str(evaluation.get('reason')))
    require(evaluation.get('status') == 'NeutrinoEvaluated',
            model + ': ' + str(evaluation.get('status')) + ': '
            + str(evaluation.get('reason')))
    require(ev.parity is not None and float(ev.parity) < 3e-8,
            model + ': corrected parity not established')

    # Save threshold C5, and independently check the symbolic reference again.
    config = bulk.update_yukawas(history, x)
    vec, layout, unpack = bulk.uv_base(config, entry)
    uv = unpack(vec, layout, mu_gev=1e5)
    eft, masses, _, _ = ev.frozen
    fast = finite_array(kernel.evaluate(uv, eft), 'corrected fast C5')
    view = SimpleNamespace(uv_threshold_state=uv,
                           eft1_threshold_state=eft,
                           mu_scalar_threshold_gev=1e3)
    exact = finite_array(evaluate_symbolic_c5(
        c5_path, extract_final_c5_inputs(view).as_config()), 'reference C5')
    parity = float(np.max(np.abs(fast - exact)) /
                   max(float(np.max(np.abs(exact))), 1e-30))
    require(parity < 3e-8, model + ': independent parity check failed')

    return {
        'model': model,
        'status': 'CandidateSharedNeutrinoEvaluated',
        'chi2': float(evaluation['chi2']),
        'prediction': {k: float(evaluation['prediction'][k]) for k in OBS},
        'screen': True,
        'rhoF': float(evaluation['rhoF']),
        'scalar_masses_muS_tev': evaluation['scalar_masses_muS_tev'],
        'parity_relative_error': parity,
        'matched_c5_at_muS': encoded_matrix(fast),
        'reference_c5_at_muS': encoded_matrix(exact),
        'source_c5_kind': kind,
        'source_matching': c5_path.relative_to(ROOT).as_posix(),
        'source_matching_sha256': digest(c5_path),
        'source_beta': beta_path.relative_to(ROOT).as_posix(),
        'source_beta_sha256': digest(beta_path),
        'legacy_kernel_normalisation': 0.5,
        'matching_physics_certified': False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', default='bulk48_v1')
    parser.add_argument('--prior-run', default='shared_a2_v1')
    parser.add_argument('--run-id', default='shared_a2_legacy_normalisation_v1')
    parser.add_argument('--extra-real', type=float, default=-0.01)
    args = parser.parse_args()
    require(0 < abs(args.extra_real) <= 0.5, 'invalid extra quartic input')
    require(args.run_id not in (args.prior_run, '', '.', '..'), 'invalid run id')
    require(all(c.isalnum() or c in '_-' for c in args.run_id),
            'run id may only contain letters, digits, hyphens, underscores')

    source = compare.BASE / args.campaign
    prior = compare.OUT / args.prior_run / 'shared_results.json'
    require(prior.is_file(), 'missing completed prior shared_results.json: ' + str(prior))
    old = compare.read(prior)
    require(old.get('format') == compare.FORMAT, 'unexpected prior report format')
    require(old.get('source_sha256') == digest(source / 'campaign.json'),
            'prior comparison belongs to different source campaign')
    settings = old.get('settings', {})
    require(settings.get('extra_real') == args.extra_real,
            'extra-real setting differs from original comparison')
    require(not settings.get('allow_poor_donor') and
            not settings.get('dry_run') and
            settings.get('bc_quartic_fraction') == 0.1 and
            settings.get('bc_mu2_fraction') == 0.0,
            'prior comparison settings do not match expected fixed benchmark')
    donor = compare.donor(source, CASE, allow_poor=False)
    require(donor['status'] == 'AcceptedDonor', 'accepted A(-2) donor unavailable')
    old_donor = old.get('donors', {}).get(CASE, {})
    require(old_donor.get('source_fit_sha256') == donor['source_fit_sha256'] and
            old_donor.get('source_fitted_config_sha256') ==
            donor['source_fitted_config_sha256'],
            'prior and current donor hashes disagree')
    previous = {(r['case'], r['model']): r for r in old['rows']}
    require(all((CASE, m) in previous for m in MODELS),
            'prior comparison is missing selected models')

    output = compare.OUT / args.run_id
    require(not output.exists(), 'output already exists: ' + str(output))
    output.mkdir(parents=True)
    # Failures leave a diagnostic directory but never overwrite the original run.
    corrected = {}
    for model in LEGACY:
        result = evaluate_legacy(model, source, donor['parameters'], args.extra_real)
        corrected[model] = result
        print(f'{model}: corrected C5 parity={result["parity_relative_error"]:.3e}, '
              f'chi2={result["chi2"]:.8g}', flush=True)
        write_json(output / (model.replace('-', 'm') + '_corrected.json'), result)

    all_rows = []
    for model in MODELS:
        if model in corrected:
            result = corrected[model]
            kind = result['source_c5_kind']
            parity = result['parity_relative_error']
        else:
            result = previous[(CASE, model)]
            kind = result.get('source_c5_kind',
                              'conditional-one-scale-mixed-scalar' if model in
                              ('B_-1', 'C_-1') else 'hierarchical-hard-plus-direct')
            parity = None
        pred = result.get('prediction') or {}
        all_rows.append({'model': model, 'case': CASE,
                         'status': result.get('status'),
                         'chi2': result.get('chi2'),
                         **{k: pred.get(k) for k in OBS},
                         'source_c5_kind': kind,
                         'legacy_parity_relative_error': parity})
    fields = list(all_rows[0])
    with (output / 'six_model_comparison.csv').open(
            'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_rows)
    write_json(output / 'manifest.json', {
        'format': FORMAT, 'scenario': CASE,
        'six_models': list(MODELS),
        'donor_fit_sha256': donor['source_fit_sha256'],
        'donor_fitted_config_sha256': donor['source_fitted_config_sha256'],
        'prior_shared_results_sha256': digest(prior),
        'legacy_normalisation_multiplier': 0.5,
        'legacy_c5_parity_tolerance': 3e-8,
        'rows': all_rows,
        'limitations': [
            'Conditional exploratory general-quartic model, not certified full general-scalar EFT',
            'B(-3)/C(-3) use legacy common-scale matching, not hierarchical hard-plus-direct',
            'B(-1)/C(-1) use conditional one-scale mixed-scalar matching',
            'General-scalar bilinear evolution and complete matching remain incomplete',
            'A(-2), D(-2), B(-1), C(-1) come from the prior unmodified shared study',
        ],
    })
    print('Combined six-model CSV:', output / 'six_model_comparison.csv')
    print('No original campaign, matching, or saved shared results were modified.')


if __name__ == '__main__':
    main()
