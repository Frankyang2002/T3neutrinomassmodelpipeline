"""Read-only 16-representation shared-A(-2) GENERAL-SCALAR candidate comparison.

Uses existing model-specific source evaluators. The only numerical correction is
0.5 applied to LegacyKernel.fun, identical to the established B(-3)/C(-3)
parity diagnosis; the original parity check is still enforced. No fitting.

Run from project root:
 python -X utf8 -m scripts.AnalyzeT3GeneralScalarSharedA2Sixteen
Re-running with --resume uses saved per-model physics results after verifying
input hashes. The report is a conditional exploratory comparison, NOT certified
full general-scalar sequential EFT matching.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scripts import CompareT3GeneralScalarSharedA2 as comparison
from scripts import RunT3GeneralScalar48Campaign as bulk

ROOT = Path(__file__).resolve().parents[1]
CASE = 'large_uv_frozen_0p5'
MODELS = tuple(comparison.MODELS)
OBS = tuple(comparison.OBS)
OUT_BASE = ROOT / 'output/studies/general_scalar_extended/six_benchmarks'
FORM = 't3-general-scalar-shared-a2-sixteen-matrix-study-v1'


def require(ok, message):
    if not ok:
        raise ValueError('Shared-A2 sixteen analysis: ' + message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temp.replace(path)


def mat(z):
    a = np.asarray(z, dtype=complex)
    require(a.shape == (3, 3) and np.isfinite(a.real).all() and np.isfinite(a.imag).all(),
            'invalid 3x3 C5')
    return a


def encode(a):
    a = mat(a)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def decode(d):
    require(isinstance(d, dict) and 'real' in d and 'imag' in d,
            'unrecognized saved C5 matrix representation')
    return mat(np.asarray(d['real'], dtype=float) + 1j * np.asarray(d['imag'], dtype=float))


def spectrum(c5_low, vev):
    # Singular values of a symmetric Majorana matrix are Takagi masses.
    m = -0.5 * vev**2 * mat(c5_low) * 1e9  # GeV -> eV
    vals = np.sort(np.linalg.svd(m, compute_uv=False))
    return vals.tolist()


def record(model, pred, chi2, c5s, c5low, vev, scheme, parity, sources, old):
    c5s = mat(c5s); c5low = mat(c5low)
    require(np.allclose(c5s, c5s.T, rtol=1e-9, atol=1e-27), model + ': C5 not symmetric')
    require(np.allclose(c5low, c5low.T, rtol=1e-9, atol=1e-27), model + ': low C5 not symmetric')
    pred = {k: float(pred[k]) for k in OBS}
    if old.get('chi2') is not None:
        require(math.isclose(float(chi2), float(old['chi2']), rel_tol=3e-6, abs_tol=5e-8),
                model + ': newly calculated chi2 disagrees with saved shared run')
        for k in OBS:
            require(math.isclose(pred[k], float(old['prediction'][k]), rel_tol=3e-6, abs_tol=1e-9),
                    model + ': observable ' + k + ' disagrees with saved shared run')
    return {'model': model, 'case': CASE, 'status': 'EvaluatedConditional',
            'matching_scheme': scheme, 'chi2': float(chi2), 'prediction': pred,
            'neutrino_masses_eV_ascending': spectrum(c5low, vev),
            'C5_muS_GeV_inv': encode(c5s), 'C5_low_GeV_inv': encode(c5low),
            'parity_relative_error': parity, 'sources': sources,
            'physics_certified_complete_general_scalar': False}


def bulk_model(model, point, source, old, *, extra_real):
    from Numerical.orchestration import RunSmallLambdaRefit as reference
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from scripts.RunT3GeneralScalarBCExploratoryFit import low_energy_fit
    from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
    from Numerical.running.SMWeinbergStage import evaluate_symbolic_c5
    registry = bulk.read(bulk.REGISTRY)
    entry, new, reused = bulk.inventory(registry, model)
    require(not reused, model + ': model has reused special quartics; wrong evaluator')
    history = bulk.historical(CASE, model, entry)
    require(not history['uv_on'], 'benchmark must be frozen UV')
    beta_path = source / 'betas' / (bulk.safe(model) + '.json')
    payload = bulk.read(beta_path)
    kernel, matching_path, kind = bulk.matching(model)
    if kind == 'legacy-common-scale':
        # Independently diagnosed in B(-3)/C(-3): inherited rapid symmetric
        # contraction lacks 1/2 present in evaluate_symbolic_c5.
        require(isinstance(kernel, bulk.LegacyKernel), model + ': unexpected legacy kernel')
        uncorrected = kernel.fun
        kernel.fun = lambda *args: 0.5 * uncorrected(*args)
    target = OscillationFitTarget.from_json(reference.TARGET)
    ev = bulk.Evaluator(CASE, model, entry, new, history, payload, kernel, target, extra_real)
    x = np.asarray([point[k] for k in bulk.NAMES], dtype=float)
    row = ev.evaluate(x, 'shared_a_m2_sixteen')
    require(row.get('screen') is True, model + ': physical screen failed: ' + str(row.get('reason')))
    require(row['status'] == 'NeutrinoEvaluated', model + ': ' + row['status'] + ': ' + str(row.get('reason')))
    require(ev.parity is not None and float(ev.parity) < 3e-8, model + ': parity check not passed')
    cfg = bulk.update_yukawas(history, x)
    vec, layout, unpack = bulk.uv_base(cfg, entry)
    uv = unpack(vec, layout, mu_gev=1e5)
    eft, _, _, _ = ev.frozen
    c5 = mat(kernel.evaluate(uv, eft))
    if kind == 'legacy-common-scale':
        view = SimpleNamespace(uv_threshold_state=uv, eft1_threshold_state=eft,
                               mu_scalar_threshold_gev=1e3)
        exact = mat(evaluate_symbolic_c5(matching_path,
                         extract_final_c5_inputs(view).as_config()))
        err = float(np.max(abs(c5-exact))/max(np.max(abs(exact)), 1e-30))
        require(err < 3e-8, model + ': symbolic legacy parity failed')
    _, low = low_energy_fit(c5, eft, target, vev=float(cfg['numerical']['vev_gev']))
    return record(model, row['prediction'], row['chi2'], c5, low.K,
                  float(cfg['numerical']['vev_gev']), kind, float(ev.parity),
                  {'matching': str(matching_path.relative_to(ROOT)),
                   'matching_sha256': sha(matching_path),
                   'beta_sha256': sha(beta_path)}, old)


def bc_model(model, point, old, *, quartic_fraction, mu2_fraction):
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from scripts import RunT3GeneralScalarBCPilot as pilot
    from scripts import RunGeneralT3NumericalC5 as c5src
    from Numerical.experimental.GeneralT3ScalarRGState import load_registry
    from Numerical.orchestration import RunSmallLambdaRefit as reference
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    registry = load_registry(bc.b7.REGISTRY)
    source = pilot.study_source(CASE, model, registry, root=ROOT)
    bundle = bc.beta_bundle(model, source, registry)
    coefficient = c5src.DEFAULT_SOURCE / ({'B_-1': 'B_m1_final', 'C_-1': 'C_m1_final'}[model] +
                        '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    compiled = bc.compile_source(coefficient)
    target = OscillationFitTarget.from_json(reference.TARGET)
    ev = bc.CandidateEvaluator(CASE, model, source, registry, bundle, compiled, target,
                         quartic_fraction=quartic_fraction, mu2_fraction=mu2_fraction)
    fit, detail = ev.evaluate(point, detail=True)
    return record(model, dict(zip(fit.observable_names, map(float, fit.prediction))),
                  float(fit.chi2), decode(detail['C5_muS_gev_inv']),
                  decode(detail['C5_low_gev_inv']),
                  float(source['cfg']['numerical']['vev_gev']),
                  'conditional-one-scale-mixed-scalar', None,
                  {'matching': str(coefficient.relative_to(ROOT)),
                   'matching_sha256': sha(coefficient)}, old)


def bp1_model(point, old, *, extra_real, beta_id):
    from scripts import FitT3BPlus1GeneralThreeStudies as b21
    from scripts import RunT3BPlus1GeneralEFT1 as b20
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration import RunSmallLambdaRefit as reference
    historical = b21.load_historical(CASE, ROOT)
    registry = b21.read(ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json')
    b20.validate_registry(registry)
    config = b21.candidate_config(historical, historical['seed'])
    uv, _ = b21.prepare_uv(ROOT, config, uv_on=False)
    initial = build_intermediate_scalar_boundary(uv, matching_scale_gev=b20.MU_F)
    payload, layout, beta_path, beta_report = b21.verify_beta_source(ROOT, beta_id, initial)
    kernel = b21.MatchedKernel(ROOT / b21.C5_REL)
    target = OscillationFitTarget.from_json(reference.TARGET)
    evaluator = b21.Evaluator(historical, payload, None, kernel, target, root=ROOT,
                       extra_fraction=None, extra_real=extra_real, rho_max=b21.RHO_MAX,
                       scalar_min_tev=b21.SCALAR_MIN_TEV, scalar_max_tev=b21.SCALAR_MAX_TEV)
    row = evaluator.evaluate(np.asarray([point[k] for k in b21.NAMES], dtype=float), label='shared_a_m2_sixteen')
    require(row.get('screen_pass') is True and row.get('chi2') is not None,
            'B_+1: screen/parity/evaluation failed: ' + str(row.get('reason')))
    return record('B_+1', row['prediction'], row['chi2'], decode(row['c5_at_muS']),
                  decode(row['c5_at_mu100']),
                  float(historical['cfg']['numerical']['vev_gev']),
                  'hierarchical-hard-plus-direct', evaluator.parity,
                  {'matching': str(b21.C5_REL), 'matching_sha256': sha(ROOT / b21.C5_REL),
                   'beta_sha256': sha(beta_path)}, old)


def coeff(a, b):
    denominator = float(np.vdot(a, a).real)
    require(denominator > 1e-65, 'zero reference C5 norm')
    k = np.vdot(a, b) / denominator
    residual = np.linalg.norm(b - k*a) / max(np.linalg.norm(b), 1e-60)
    return {'scale_real': float(k.real), 'scale_imag': float(k.imag),
            'relative_frobenius_residual': float(residual)}


def numerical_table(rows):
    out = []
    for m in MODELS:
        r = rows[m]
        d = {'model': m, 'status': r['status'], 'matching_scheme': r.get('matching_scheme'),
             'chi2': r.get('chi2'), 'mass1_eV': None, 'mass2_eV': None, 'mass3_eV': None,
             **{k: None for k in OBS}, 'C5_scale_vs_A2_real': None,
             'C5_scale_vs_A2_imag': None, 'C5_residual_vs_A2': None}
        if r['status'] == 'EvaluatedConditional':
            d.update({k: r['prediction'][k] for k in OBS})
            d.update(dict(zip(('mass1_eV','mass2_eV','mass3_eV'), r['neutrino_masses_eV_ascending'])))
            d['C5_scale_vs_A2_real'] = r['matrix_comparison_to_A2_muS']['scale_real']
            d['C5_scale_vs_A2_imag'] = r['matrix_comparison_to_A2_muS']['scale_imag']
            d['C5_residual_vs_A2'] = r['matrix_comparison_to_A2_muS']['relative_frobenius_residual']
        out.append(d)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign', default='bulk48_v1')
    p.add_argument('--prior-run', default='shared_a2_v1')
    p.add_argument('--run-id', default='shared_a2_16_matrix_v1')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--extra-real', type=float, default=-0.01)
    a = p.parse_args()
    require(a.run_id.isascii() and a.run_id.replace('_','').replace('-','').isalnum(),
            'unsafe run id')
    source = comparison.BASE / a.campaign
    saved = OUT_BASE / a.prior_run / 'shared_results.json'
    prior = comparison.read(saved)
    require(prior.get('format') == comparison.FORMAT and
            prior.get('source_sha256') == sha(source / 'campaign.json'),
            'incompatible shared campaign')
    settings = prior['settings']
    require(settings.get('extra_real') == a.extra_real and
            not settings.get('dry_run') and not settings.get('allow_poor_donor'),
            'incompatible shared benchmark settings')
    donor = comparison.donor(source, CASE, allow_poor=False)
    require(donor['status'] == 'AcceptedDonor', 'no accepted A(-2) donor')
    require(prior['donors'][CASE]['source_fit_sha256'] == donor['source_fit_sha256'] and
            prior['donors'][CASE]['source_fitted_config_sha256'] == donor['source_fitted_config_sha256'],
            'donor changed')
    prior_rows = {(r['case'], r['model']): r for r in prior['rows']}
    require(all((CASE, m) in prior_rows for m in MODELS), 'shared results lack one of 16 models')
    dest = OUT_BASE / a.run_id
    require(dest != OUT_BASE / a.prior_run, 'refusing to overwrite shared results')
    require(not dest.exists() or a.resume, 'output exists; use --resume or new --run-id')
    dest.mkdir(parents=True, exist_ok=True)
    all_rows = {}
    for m in MODELS:
        path = dest / 'models' / (bulk.safe(m) + '.json')
        if a.resume and path.is_file():
            row = comparison.read(path)
            require(row.get('source_shared_sha256') == sha(saved) and
                    row.get('source_donor_sha256') == donor['source_fit_sha256'],
                    m + ': stale resume checkpoint')
            print(m, 'resume', row['status'], flush=True)
        else:
            print(m, 'evaluating...', flush=True)
            try:
                if m in ('B_-1', 'C_-1'):
                    row = bc_model(m, donor['parameters'], prior_rows[CASE,m],
                                   quartic_fraction=settings['bc_quartic_fraction'],
                                   mu2_fraction=settings['bc_mu2_fraction'])
                elif m == 'B_+1':
                    row = bp1_model(donor['parameters'], prior_rows[CASE,m],
                                    extra_real=a.extra_real, beta_id=settings['bp1_beta_id'])
                else:
                    row = bulk_model(m, donor['parameters'], source,
                                     prior_rows[CASE,m], extra_real=a.extra_real)
            except Exception as e:
                row = {'model': m, 'case': CASE, 'status': 'EvaluationFailed',
                       'reason': f'{type(e).__name__}: {e}'[:850],
                       'previous_status': prior_rows[CASE,m].get('status'),
                       'physics_certified_complete_general_scalar': False}
            row['source_shared_sha256'] = sha(saved)
            row['source_donor_sha256'] = donor['source_fit_sha256']
            save(path, row)
            print(m, row['status'], row.get('reason',''), flush=True)
        all_rows[m] = row
    if all_rows['A_-2']['status'] == 'EvaluatedConditional':
        reference = decode(all_rows['A_-2']['C5_muS_GeV_inv'])
        for m, row in all_rows.items():
            if row['status'] == 'EvaluatedConditional':
                row['matrix_comparison_to_A2_muS'] = coeff(reference, decode(row['C5_muS_GeV_inv']))
                save(dest / 'models' / (bulk.safe(m)+'.json'), row)
    else:
        print('A(-2) C5 unavailable, proportionality comparisons skipped', flush=True)
    pairs = {}
    for b,c in [('B_-3','C_-3'),('B_-1','C_-1'),('A_-2','D_-2')]:
        if all_rows[b]['status'] == all_rows[c]['status'] == 'EvaluatedConditional':
            mb=decode(all_rows[b]['C5_muS_GeV_inv']);mc=decode(all_rows[c]['C5_muS_GeV_inv'])
            p=coeff(mb,mc)
            p['solar_splitting_ratio'] = all_rows[c]['prediction']['delta_m21_sq_ev2'] / all_rows[b]['prediction']['delta_m21_sq_ev2']
            p['atmospheric_splitting_ratio'] = all_rows[c]['prediction']['delta_m3l_sq_ev2'] / all_rows[b]['prediction']['delta_m3l_sq_ev2']
            pairs[b+'_to_'+c] = p
    table = numerical_table(all_rows)
    with (dest / 'sixteen_model_comparison.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(table[0]))
        writer.writeheader(); writer.writerows(table)
    save(dest / 'matrix_analysis.json', {'format': FORM, 'case': CASE,
          'model_count':16, 'evaluated_count':sum(r['status']=='EvaluatedConditional' for r in all_rows.values()),
          'donor_fit_sha256':donor['source_fit_sha256'], 'prior_shared_sha256':sha(saved),
          'pairwise':pairs, 'rows':table, 'limitations':[
              'All results are conditional candidate general-scalar calculations, not certified general EFT matching.',
              'General quadratic sector incomplete in bulk models; model-dependent matching schemes.',
              'B/C(-1) uses one-scale conditional mixed-scalar calculation, distinct from bulk hard+direct.',
              'Extreme-hypercharge legacy C5 uses measured and parity-verified factor 0.5 correction.',
              'One accepted A(-2) Yukawa donor transferred without fitting; this is not a viability exclusion test.',
              'C5 matrix ratios evaluated at muS=1 TeV; low-energy masses at mu=100 GeV.'
          ]})
    lines=['# Sixteen general-scalar representations: shared A(-2) Yukawas', '',
           'Scenario: `large_uv_frozen_0p5`; frozen UV, evolved EFT1 and SM+Weinberg.',
           'These are **conditional exploratory** calculations, not certified full general-scalar matching.', '',
           f"Successful numerical evaluations: **{sum(r['status']=='EvaluatedConditional' for r in all_rows.values())}/16**.", '',
           '| Model | Status | chi2 | Δm21² [eV²] | Δm3l² [eV²] | C5 residual vs A(-2) |',
           '|---|---|---:|---:|---:|---:|']
    for r in table:
        fmt=lambda x: '—' if x is None else f'{x:.6g}'
        lines.append(f"| {r['model']} | {r['status']} | {fmt(r['chi2'])} | {fmt(r['delta_m21_sq_ev2'])} | {fmt(r['delta_m3l_sq_ev2'])} | {fmt(r['C5_residual_vs_A2'])} |")
    lines.extend(['', '## Pairwise C5 scaling at 1 TeV', ''])
    for key,pair in pairs.items():
        lines.append(f"- {key}: k=({pair['scale_real']:.9g})+i({pair['scale_imag']:.9g}); residual={pair['relative_frobenius_residual']:.6g}; solar ratio={pair['solar_splitting_ratio']:.9g}; atmospheric ratio={pair['atmospheric_splitting_ratio']:.9g}.")
    lines.extend(['', 'A small matrix residual establishes approximate threshold C5 proportionality within the specific implemented schemes; it does not certify representation physics.',
                  'Some fits may fail experimental screening even when their forward numerical evaluation succeeds.',
                  'Individual failures and full matrices are saved in `models/` JSON files.', ''])
    (dest / 'analysis.md').write_text('\n'.join(lines), encoding='utf-8')
    print('Saved:', dest / 'sixteen_model_comparison.csv')
    print('Saved:', dest / 'matrix_analysis.json')
    print('Saved:', dest / 'analysis.md')
    print('Evaluated:', sum(r['status']=='EvaluatedConditional' for r in all_rows.values()), '/16')


if __name__=='__main__':
    main()
