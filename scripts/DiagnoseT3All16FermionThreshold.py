"""Reconstruct fitted fermion-threshold formal LLSS tensors for all 16 T3 models.

Project path: scripts/DiagnoseT3All16FermionThreshold.py

Uses the exact UV and (for B_-1/C_-1 only) EFT1 trajectories used in
RunT3GeneralSequentialFitAblation, with the same source construction.
Does not refit benchmarks or alter the existing corrected LLSS campaign.

C_LLSS is a six-index formal Matchete commuting-spinor tensor in GeV^-1.
Its absolute GammaCC-to-Weyl normalization is not certified. The complete
hard one-loop C5(muF) is not calculated. For mixed B/C, the reconstructed
fermion-source tensor uses U(muS) as in the saved corrected campaign and is
NOT a strictly same-scale muF matching coefficient.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import traceback

import numpy as np

from scripts import RunT3GeneralSequentialFitAblation as seq
from scripts import RunT3GeneralAll16FromReference as general

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_BASE = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
F_SCALE = seq.SCALES[1]


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def norm(x):
    return float(np.linalg.norm(np.asarray(x, np.complex128)))


def encode(z):
    z = complex(z)
    return {'real': float(z.real), 'imag': float(z.imag)}


def nonmixed(model, saved, canonical, wolframscript):
    from Numerical.orchestration.GeneralRunningBenchmark import initial_state
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.GeneralT3DynamicFermionWilson import nonmixed_boundary

    beta_uv, _beta_eft, _unused, _kernel, _coefficient, _kind = general.nonmixed_inputs(
        model, canonical, generate_general_c5=False, wolframscript=wolframscript)
    overrides = read(saved / 'general_uv_overrides.json')
    _, initial = initial_state(beta_uv, mu_uv=seq.SCALES[0],
                               mu_f=F_SCALE, mu_s=seq.SCALES[2],
                               overrides=overrides)
    sol = run_uv_segment(initial, beta_uv, F_SCALE, rtol=1e-7, atol=1e-10)
    source, info = nonmixed_boundary(model, sol.final_state, repo_root=ROOT,
                                     mu_f=F_SCALE)
    return source, info, {'uv_nfev': int(sol.nfev),
                          'input_sha256': {'general_uv_overrides': digest(saved / 'general_uv_overrides.json')},
                          'mixed_uses_scalar_U_muS': False}


def mixed(model, saved):
    """Match the original mixed model U=1,E=1 branch, including U(muS)."""
    from scripts import RunT3BCGeneralFromReference as runner
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from scripts import RunT3GeneralScalarBCPilot as pilot
    from scripts import RunT3ZeroExtraBoundaryBenchmark as core
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state, project_general_after_fermion,
        pack_general_scalar_state)
    from Numerical.orchestration.GeneralRunningBenchmark import _state_values
    from Numerical.running.GeneralT3MixedExactFermionWilson import mixed_exact_boundary

    fit = read(saved / 'result.json')
    source_config = runner.build_reference_source(model)
    source_config['cfg'] = read(saved / 'fitted_config.json')
    registry, new, reused = (source_config['registry'], source_config['new'],
                             source_config['reused'])
    q, m = seq.mixed_fit_fractions(fit)
    beta_bundle = bc.beta_bundle(model, source_config, registry)
    initial, uv_layout, _, _ = pilot.prepare_initial(
        source_config, registry,
        'general_quartics_and_mixing' if m > 0 else 'general_quartics', q, m)
    uv_sol = core.integrate_segment(
        initial,
        lambda t, x: core.rhs_at(t, x, layout=uv_layout,
                                 registry=registry, betas=beta_bundle['UV'],
                                 stage='UV'),
        seq.SCALES[0], F_SCALE, n_points=3)
    uv_vec = uv_sol.y[:, -1].copy()
    uv = unpack_general_scalar_state(uv_vec[:uv_layout.size], uv_layout,
                                     registry, uv=True, mu_gev=F_SCALE)
    heavy = uv.baseline
    eft = project_general_after_fermion(uv, registry)
    eft_vec, eft_layout = pack_general_scalar_state(eft, registry)
    eft_vec = np.concatenate((eft_vec, uv_vec[uv_layout.size:]))
    eft_beta = core.check_payload(
        beta_bundle['payloads']['EFT1'], model=model, stage='EFT1',
        rep=eft_layout.representation,
        base_names=tuple(eft_layout.names[:-len(new)]),
        new_names=new, reused_names=reused)
    eft_sol = core.integrate_segment(
        eft_vec,
        lambda t, x: core.rhs_at(t, x, layout=eft_layout,
                                 registry=registry, betas=eft_beta,
                                 stage='EFT1'),
        F_SCALE, seq.SCALES[2], n_points=3)
    end_vec = eft_sol.y[:, -1].copy()
    eft_endpoint = unpack_general_scalar_state(
        end_vec[:eft_layout.size], eft_layout,
        registry, uv=False, mu_gev=seq.SCALES[2]).baseline
    state = _state_values(eft_endpoint)
    endpoint = core.snapshot(end_vec, eft_layout, registry, stage='EFT1',
                             mu_gev=seq.SCALES[2], new=new, reused=reused)
    for name, value in endpoint['new_quartics'].items():
        state.setdefault(name, value)
    for name, value in endpoint['reused_quartics'].items():
        state.setdefault(name, value)
    state['muT3GeneralS1S2'] = endpoint['quadratic_mu2_gev2']
    boundary, info = mixed_exact_boundary(model, heavy, state, repo_root=ROOT,
                                          mu_f=F_SCALE)
    return boundary, info, {
        'uv_nfev': int(uv_sol.nfev), 'eft1_nfev': int(eft_sol.nfev),
        'input_sha256': {'fitted_config': digest(saved / 'fitted_config.json'),
                         'fitted_result': digest(saved / 'result.json')},
        'mixed_uses_scalar_U_muS': True,
        'warning': ('Mixed-source scalar mixing U is evaluated at muS in the '
                    'corrected campaign; the tensor is not a strictly local '
                    'muF-only Wilson boundary.')}


def check_saved(model, case_record, source, info):
    if case_record.get('model') != model or case_record.get('status') != 'Computed':
        raise ValueError('Saved campaign model not computed')
    if '1111' not in case_record['cases']:
        raise ValueError('Saved 1111 case missing')
    saved = case_record.get('provenance', {})
    if model in ('B_-1', 'C_-1'):
        expected = saved.get('mixed_native_CG_exact_majorana_source', {}).get('11')
        if not isinstance(expected, dict):
            raise ValueError('Saved mixed U=1,E=1 source diagnostics missing')
        expected_indices = expected.get('source_indices')
    else:
        expected = saved.get('fermion_tree', {}).get('1')
        if not isinstance(expected, dict):
            raise ValueError('Saved U=1 F-tree source diagnostics missing')
        expected_indices = [r['source_index'] for r in expected.get('sources', [])]
    if list(source.source_indices) != [int(k) for k in expected_indices]:
        raise ValueError('Source indices do not reproduce saved campaign')
    current_masses = np.asarray(info['heavy_mass_eigenvalues_gev'], float)
    expected_masses = np.asarray(expected['heavy_mass_eigenvalues_gev'], float)
    mismatch = float(np.linalg.norm(current_masses-expected_masses) /
                     max(np.linalg.norm(expected_masses), 1e-40))
    if mismatch > 1e-7:
        raise ValueError(f'Heavy-fermion masses differ from saved campaign: {mismatch:.3e}')
    # For B/C, the same scalar mixing matrix must also be reproduced.
    mixing_residual = None
    if model in ('B_-1', 'C_-1'):
        def matrix(raw):
            return np.asarray(raw['real'],float)+1j*np.asarray(raw['imag'],float)
        a, b = matrix(info['mixing_U_muS']), matrix(expected['mixing_U_muS'])
        mixing_residual = norm(a-b)/max(norm(b),1e-40)
        if mixing_residual > 1e-7:
            raise ValueError(f'Mixed scalar U(muS) differs from saved campaign: {mixing_residual:.3e}')
    return mismatch, mixing_residual


def write_csv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def analyze(model, args, source_root, canonical, out):
    saved = source_root / seq.tag(model)
    case_path = CAMPAIGN_BASE / args.reference_run / seq.tag(model) / 'cases.json'
    if not case_path.is_file():
        raise FileNotFoundError(case_path)
    record = read(case_path)
    if model in ('B_-1', 'C_-1'):
        boundary, info, run = mixed(model, saved)
    else:
        boundary, info, run = nonmixed(model, saved, canonical, args.wolframscript)
    masses_residual, mixing_residual = check_saved(model, record, boundary, info)
    tensor = np.asarray(boundary.raw, np.complex128)
    by = np.asarray(boundary.by_source, np.complex128)
    if tensor.ndim != 6 or tensor.shape[:4] != (3,3,2,2) or tensor.shape[-1] != tensor.shape[-2]:
        raise ValueError(f'Unexpected tensor axes: {tensor.shape}')
    if by.shape != (len(boundary.source_indices), *tensor.shape):
        raise ValueError('Source tensor array shape mismatch')
    total = norm(tensor)
    if not np.isfinite(total) or total <= 0:
        raise ValueError('Nonpositive or nonfinite formal LLSS norm')
    residual = norm(tensor-by.sum(axis=0))/total
    if residual > 1e-10:
        raise ValueError(f'LLSS source addition failed: {residual:.3e}')
    case = record['cases']['1111']
    log_norm = float(case['c5_log_frobenius'])
    bubble = case.get('c5_joint_scalar_bubble')
    if bubble is None:
        raise ValueError('Saved finite bubble missing')
    bubble_norm = norm(np.asarray(bubble['real'])+1j*np.asarray(bubble['imag']))
    folder = out / seq.tag(model)
    folder.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(folder / 'LLSS_muF_formal.npz',
                        C_LLSS_formal=tensor,
                        C_LLSS_formal_per_source=by,
                        source_indices=np.asarray(boundary.source_indices,dtype=int),
                        mu_F_gev=np.asarray(F_SCALE))
    source_rows = []
    for i, idx in enumerate(boundary.source_indices):
        part = by[i]
        source_rows.append({
            'model':model, 'source_index':int(idx),
            'source_norm_gev_inverse':norm(part),
            'signed_projection_percent':100*float(np.vdot(tensor,part).real)/total**2,
            'source_norm_over_total':norm(part)/total})
    write_csv(folder / 'LLSS_sources.csv',source_rows,list(source_rows[0]))
    with (folder / 'LLSS_nonzero_components.csv').open('w',newline='',encoding='utf-8') as fh:
        writer=csv.writer(fh)
        writer.writerow(['alpha_0based','beta_0based','su2_i_0based',
                         'su2_j_0based','scalar_a_0based','scalar_b_0based',
                         'real_gev_inverse','imag_gev_inverse'])
        for indices in np.argwhere(tensor != 0):
            coord=tuple(int(i) for i in indices)
            z=complex(tensor[coord]);writer.writerow([*coord,z.real,z.imag])
    summary={
        'model':model, 'mu_F_gev':F_SCALE, 'source_campaign':args.source_run,
        'reference_campaign':args.reference_run, 'reference_case':'1111',
        'tensor_shape':list(tensor.shape),
        'tensor_axes':['lepton_flavour_alpha','lepton_flavour_beta',
                       'lepton_su2_i','lepton_su2_j',
                       'heavy_real_scalar_a','heavy_real_scalar_b'],
        'scalar_labels':list(boundary.labels),
        'C_LLSS_formal_frobenius_gev_inverse':total,
        'C5_LLSS_log_frobenius_gev_inverse':log_norm,
        'C5_finite_scalar_bubble_frobenius_gev_inverse':bubble_norm,
        'log_over_formal_LLSS_norm':log_norm/total,
        'source_additivity_residual':residual,
        'saved_muF_heavy_mass_relative_discrepancy':masses_residual,
        'saved_muS_scalar_mixing_relative_discrepancy':mixing_residual,
        'source_contributions':source_rows,
        'source_status':boundary.source_status,
        'matching_details':run,
        'matching_normalization_certified':False,
        'complete_hard_fermion_C5_muF':'not calculated',
        'full_C5_muF':'unknown',
        'reference_cases_sha256':digest(case_path),
        'source_seed_sha256':digest(boundary.path),
    }
    (folder / 'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return summary, source_rows


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-run',default='thesis16_merged_v1')
    p.add_argument('--reference-run',default='corrected_llss16_v4')
    p.add_argument('--models',default=None,help='comma-separated model IDs; default all 16')
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--output',type=Path,default=None)
    args=p.parse_args(argv)
    models=tuple(x.strip() for x in args.models.split(',') if x.strip()) if args.models else seq.MODELS
    if not models or len(set(models))!=len(models) or set(models)-set(seq.MODELS):
        p.error('Invalid models selection')
    source_root=general.OUTPUT/args.source_run
    settings=read(source_root/'settings.json')
    canonical=settings.get('canonical_run_id','full16_general_v1')
    out=args.output or (ROOT/'Reports/output/thesis_qualitative/'
                        'eft1_llss_contributions'/args.reference_run/'fermion_threshold_all16')
    out.mkdir(parents=True,exist_ok=True)
    summaries=[];sources=[];errors=[]
    for model in models:
        try:
            row,source_rows=analyze(model,args,source_root,canonical,out)
            summaries.append(row);sources.extend(source_rows)
            print(f"{model}: LLSS={row['C_LLSS_formal_frobenius_gev_inverse']:.9e}; "
                  f"C5log={row['C5_LLSS_log_frobenius_gev_inverse']:.9e}; "
                  f"ratio={row['log_over_formal_LLSS_norm']:.6e}; "
                  f"sources={len(source_rows)}",flush=True)
        except Exception as e:
            errors.append({'model':model,'error':f'{type(e).__name__}: {e}',
                           'traceback':traceback.format_exc(limit=8)})
            print(f'{model}: FAILED {type(e).__name__}: {e}',flush=True)
    summary_fields=['model','C_LLSS_formal_frobenius_gev_inverse',
                    'C5_LLSS_log_frobenius_gev_inverse',
                    'C5_finite_scalar_bubble_frobenius_gev_inverse',
                    'log_over_formal_LLSS_norm','source_additivity_residual',
                    'saved_muF_heavy_mass_relative_discrepancy',
                    'saved_muS_scalar_mixing_relative_discrepancy']
    write_csv(out/'all16_threshold_summary.csv',summaries and
              [{k:x.get(k) for k in summary_fields} for x in summaries] or [],summary_fields)
    write_csv(out/'all16_source_contributions.csv',sources,
              ['model','source_index','source_norm_gev_inverse',
               'signed_projection_percent','source_norm_over_total'])
    manifest={
        'models_requested':list(models),'models_computed':[r['model'] for r in summaries],
        'errors':errors,'source_run':args.source_run,'reference_run':args.reference_run,
        'fermion_matching':'tree LLSS formal tensor only',
        'hard_C5_at_muF':'not calculated; full C5(muF) unknown',
        'absolute_spinor_normalization_certified':False,
        'mixed_model_caveat':('B_-1 and C_-1 use the corrected campaign\'s '
                              'muS scalar mixing U inside the fermion-source '
                              'tensor; not a strictly local muF Wilson boundary'),
        'log_over_LLSS_norm_warning':('Ratio of norms across different tensor '
                                      'spaces, not an anomalous dimension'),
        'no_fits':True,'no_scalar_matching_or_SMEFT_rerun':True,
    }
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print('Output:',out)
    return 0 if not errors else 2


if __name__=='__main__':
    raise SystemExit(main())
