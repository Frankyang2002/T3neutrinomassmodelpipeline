#!/usr/bin/env python3
"""Batch 7: source-preserving all-16 EFT1 CT / native Laurent-pole audit.

Inputs: completed saved Batches 1, 2, 4 and 5 (unchanged).
Outputs: a NEW Batch-7 directory (no modification of previous batch files).
Use --skip-low-energy-running for fast pole-only diagnostics, or the default
for all three conditional downstream SMEFT predictions at three scales.

The existing project pipeline is imported and reused; no fabricated masses,
UV counterterms, or reconstructed UV spinor normalisation.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.wilson.basis import ALL_MODELS
from Numerical.parallel_t3.wilson.source_audited_all16 import hash_file
from Numerical.parallel_t3.wilson.quartic_channels_batch2 import named_quartics
from Numerical.parallel_t3.wilson.msbar_counterterm_batch6 import (
    SCALE_FACTORS, SCENARIOS, source_scalar_counterterm,
    candidate_adm_counterterms, analytic_ct_from_soft, matching_scale_scenarios,
)
from Numerical.parallel_t3.wilson.native_pole_closure_batch7 import (
    SCHEMA, STATUS, MODELS, relative, physical_total, matrix_record,
    source_by_native_label, source_cancellation_ledger,
    native_pole_evidence, scalar_integral_checks, certify_matching,
)
from Numerical.parallel_t3.threshold.pipeline import evolve_existing_smeft
from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables
from scripts.RunT3CountertermBatch6 import (
    _inputs, read_json, save_json, save_npz,
    MF, MS, BASE2, BASE4, BASE5,
)

DEST = Path('output/parallel_t3/native_pole_closure_all16_batch7_v1')


def _input_data(root: Path, model: str) -> dict:
    """Use Batch6's strictly checked source hashes, tensor and provenance checks."""
    src = _inputs(root, model)
    original = src['MF_quartic']
    t = np.asarray(original['mu_gev'], float)
    ix = np.flatnonzero(np.isclose(t, MF, atol=1e-6, rtol=1e-11))
    if len(ix) != 1:
        raise ValueError('No unique MF scalar-tensor source')
    qfull = np.asarray(original['lambda4'][ix[0]], float)
    names, channels, info = named_quartics(root, model, 'MF', qfull)
    q = channels[:, :4, :4, 4:, 4:]
    if relative(q.sum(axis=0), src['hard']['HHSS_quartic_MF']) > 1e-9:
        raise ValueError('MF RGBeta source quartics do not reconstruct hard triangle')
    return {'source': src, 'quartics': q, 'quartic_names': names, 'quartic_info': info}


def _native_constants(src: dict) -> tuple[np.ndarray, np.ndarray, list[str]]:
    hard, native, pole = src['hard'], src['native'], src['pole']
    masses = np.asarray(hard['heavy_masses_gev'], float)
    factor = 1 + 2 * np.log(MF / masses)
    if masses.shape != (3,) or np.min(abs(factor)) < 0.05:
        raise ValueError('Native Laurent inversion singular or wrong mass spectrum')
    raw = np.asarray(native['C5_native_by_kernel_and_heavy'], complex)
    constants = raw / factor[None, :, None, None]
    if relative(constants.sum(axis=(0, 1)), pole['C5_native_raw_pole_residue']) > 1e-9:
        raise ValueError('Native Laurent pole reconstruction fails Batch5 reference')
    records = src['native_meta']['kernel_records']
    if len(records) != len(constants):
        raise ValueError('Native Matchete kernel record count changed')
    labels = [x['native_quartic'] for x in records]
    return constants, masses, labels


def _soft_samples(pole: dict) -> dict[float, np.ndarray]:
    return {0.5: pole['EFT1_soft_finite_mu_m05'],
            1.0: pole['EFT1_soft_finite_mu_m10'],
            2.0: pole['EFT1_soft_finite_mu_m20']}


def _compare_existing_batch6(root: Path, model: str, ctnorm: float,
                             nativegap: float) -> dict:
    """Cross-check saved 15/16 report without requiring B6 to have succeeded."""
    p = root / 'output/parallel_t3/eft1_msbar_counterterms_all16_batch6_v1' / tag(model) / 'counterterm_and_scale_ledger.json'
    if not p.is_file():
        return {'batch6_output_present': False, 'batch6_comparison_passed': None}
    data = read_json(p)
    old_ct = float(data['original_triangle_B0_CT_relative'])
    old_gap = float(data['native_raw_pole_vs_explicit_EFT1_CT_relative'])
    # Differences are compared absolutely because these are small residual ratios.
    passed = abs(old_ct - ctnorm) <= 1e-9 and abs(old_gap - nativegap) <= 1e-9
    if not passed:
        raise ValueError('Batch7 diagnostics disagree with existing Batch6 physical totals')
    return {'batch6_output_present': True, 'batch6_comparison_passed': True,
            'batch6_CT_B0_relative': old_ct,
            'batch6_native_CT_relative': old_gap}


def _scan(src: dict, factors, constants: np.ndarray,
          masses: np.ndarray, *,
          mu_low: float, vev: float, skip_running: bool) -> tuple[dict, dict, list[dict]]:
    hard, strict, native, pole = src['hard'], src['strict'], src['native'], src['pole']
    scanned = matching_scale_scenarios(
        native_constant_by_kernel_fermion=constants,
        masses_gev=masses,
        c5_analytic_hard_mf=hard['C5_hard_F_triangle_MF'],
        c5_native_formal_mf=native['C5_native_finite_MF'],
        c5_running_delta_mf_to_ms=strict['C5_EFT1_strict_running_delta'],
        c5_scalar_bubble_ms=src['B5_MS']['C5_scalar_bubble_MS'],
        soft_finite=factors, mu_f=MF, mu_s=MS, factors=SCALE_FACTORS)
    references = {
        'analytic_triangle': src['B2_MS']['C5_SMEFT_MS'],
        'native_formal': src['B4_MS']['C5_SMEFT_MS'],
        'native_as_UV_minus_soft_hypothesis': src['B5_MS']['C5_SMEFT_MS'],
    }
    nominal_residuals = {
        k: relative(scanned['scales'][1.0][k]['C5_SMEFT_MS'], v)
        for k, v in references.items()}
    if any(v > 2e-8 for v in nominal_residuals.values()):
        raise ValueError('Batch7 nominal C5 differs from saved Batches 2, 4, 5')

    branch_matrices = {}
    outcome = {}
    all_rows = []
    for factor, branches in scanned['scales'].items():
        outcome[str(factor)] = {}
        for scenario, val in branches.items():
            key = f'{scenario}_scale_{factor:g}'
            c5 = physical_total(val['C5_SMEFT_MS'], key+' at MS')
            branch_matrices[key+'_C5_SMEFT_MS'] = c5
            entry = {'model': None, 'subtraction_scale_factor': factor,
                     'mu_subtraction_gev': MF * factor,
                     'scenario': scenario,
                     'C5_SMEFT_MS_norm_GeV_inverse': float(np.linalg.norm(c5)),
                     'mass_results_conditional_not_physical_MSbar': True,
                     'fixed_MF_source_couplings': True}
            if not skip_running:
                solved = evolve_existing_smeft(
                    c5, src['background_arrays'], src['background_meta'],
                    mu_s=MS, mu_low=mu_low)
                if solved.get('solver_success') is not True:
                    raise RuntimeError('Existing SMEFT Wilson/RG evolution failed')
                low = physical_total(solved['c5_low'], key+' low C5')
                mass = -0.5 * vev * vev * low
                obs, _ = neutrino_observables(mass, src['background_arrays']['ye'])
                entry.update({'masses_eV': obs['masses_eV'],
                    'dm21_sq_eV2': obs['dm21_sq_eV2'],
                    'dm31_sq_eV2': obs['dm31_sq_eV2'],
                    'mixing_defined': obs['mixing_defined'],
                    'sin2_theta12': obs['sin2_theta12'],
                    'sin2_theta13': obs['sin2_theta13'],
                    'sin2_theta23': obs['sin2_theta23'],
                    'SMEFT_solver_nfev': int(solved['nfev'])})
                branch_matrices[key+'_C5_at_low'] = low
                all_rows.append(entry)
            outcome[str(factor)][scenario] = entry
    return {'by_factor': outcome, 'nominal_reproduction': nominal_residuals,
            'matching_scale_only_not_new_threshold': True}, branch_matrices, all_rows


def run_one(root: Path, out: Path, model: str, *, mu_low: float,
            vev: float, skip_running: bool = False) -> dict:
    row = {'model': model, 'status': 'BLOCKED', 'stage': 'load_source_provenance'}
    try:
        data = _input_data(root, model)
        src, q, names = data['source'], data['quartics'], data['quartic_names']
        hard, pole, strict = src['hard'], src['pole'], src['strict']

        row['stage'] = 'derive_and_audit_all_EFT1_scalar_source_counterterms'
        ct = source_scalar_counterterm(hard['C_LLSS_per_fermion_canonical'], q, names)
        # CORE BUG FIX: no source-level flavour symmetry requirement.
        # The original B6 manifest cm(...) incorrectly uses symmetric3(...) on
        # every quartic summand. Per-source antisymmetric contributions are kept.
        ct_sum = physical_total(ct['CT_scalar_LLHH'], 'EFT1 scalar CT total')
        source_ledger = source_cancellation_ledger(
            ct['CT_per_channel'], names, ct_sum, label='EFT1 RGBeta quartics')
        f_names = ['F_1', 'F_2', 'F_3']
        heavy_ledger = source_cancellation_ledger(
            ct['CT_per_heavy'], f_names, ct_sum,
            label='EFT1 three integrated-out fermions')

        row['stage'] = 'independent_analytic_B0_and_ADM_checks'
        samples = _soft_samples(pole)
        slope = analytic_ct_from_soft(samples, SCALE_FACTORS)
        scalar_check = scalar_integral_checks(
            explicit_ct=ct_sum,
            batch5_soft_pole=pole['C5_EFT1_scalar_loop_pole_residue'],
            b0_scale_ct=slope)
        adm = candidate_adm_counterterms(strict, float(src['strict_meta']['log_ms_over_mf']))

        row['stage'] = 'reconstruct_native_source_Laurent_pole'
        constants, masses, native_labels = _native_constants(src)
        groups = source_by_native_label(constants, native_labels)
        native_sum = physical_total(sum(groups.values(), np.zeros((3, 3), complex)),
                                    'native total pole')
        native_sources = source_cancellation_ledger(
            list(groups.values()), list(groups.keys()), native_sum,
            label='native Matchete quartic kernels')
        pole_evidence = native_pole_evidence(
            native=native_sum, eft1_ct=ct_sum,
            candidate_adm=adm['CT_candidate_total_7channel'],
            candidate_scalar_adm=adm['CT_candidate_scalar_pair'])
        old_comparison = _compare_existing_batch6(
            root, model, scalar_check['CT_vs_B0_scale_slope_relative'],
            pole_evidence['native_vs_scalar_CT_relative'])

        row['stage'] = 'reproduce_low_energy_scenarios_without_fitted_pole'
        observables, scans, rows = _scan(
            src, samples, constants, masses,
            mu_low=mu_low, vev=vev, skip_running=skip_running)
        for r in rows:
            r['model'] = model
        source_info = data['quartic_info']
        certification = certify_matching()
        if certification['complete_matching_claim_allowed']:
            raise RuntimeError('Incomplete native UV matching unexpectedly certified')

        row['stage'] = 'write_source_resolved_audit'
        folder = out/tag(model)
        folder.mkdir(parents=True, exist_ok=True)
        npz_path = folder/'source_poles_and_counterterms.npz'
        npz_content = {
            'C5_EFT1_CT_by_RGBeta_quartic_raw': ct['CT_per_channel'],
            'C5_EFT1_CT_by_RGBeta_quartic_and_fermion_raw': ct['CT_per_channel_per_heavy'],
            'C5_EFT1_CT_by_fermion_raw': ct['CT_per_heavy'],
            'C5_EFT1_scalar_total': ct_sum,
            'C5_native_Laurent_raw_by_kernel_and_fermion': constants,
            'C5_native_Laurent_raw_by_named_quartic': np.stack(list(groups.values())),
            'C5_native_raw_pole_total': native_sum,
            'C5_native_minus_EFT1_scalar_CT': native_sum-ct_sum,
            'C5_native_minus_EFT1_candidate_ADM': native_sum-adm['CT_candidate_total_7channel'],
            'C5_EFT1_candidate_ADM_total': adm['CT_candidate_total_7channel'],
            'C5_EFT1_candidate_ADM_scalar_pair': adm['CT_candidate_scalar_pair'],
            'C5_B0_determined_scalar_CT': slope,
        }
        for name, val in adm['CT_per_ADM_channel'].items():
            npz_content['C5_EFT1_candidate_ADM_'+name] = val
        save_npz(npz_path, npz_content)
        scan_path = folder/'matching_scale_scenario_matrices.npz'
        save_npz(scan_path, scans)
        if rows:
            csv_path = folder/'conditional_mass_and_mixing_scan.csv'
            tmp = csv_path.with_suffix('.partial.csv')
            with tmp.open('w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            tmp.replace(csv_path)
        else:
            csv_path = None

        ledger_path = folder/'source_pole_and_matching_ledger.json'
        ledger = {
            'schema': SCHEMA, 'model': model, 'status': STATUS,
            'source_SHA256': src['h'], 'quartic_source_metadata': source_info,
            'MF_GeV': MF, 'MS_GeV': MS, 'mu_low_GeV': mu_low, 'vev_GeV': vev,
            'fermion_threshold_count': 1, 'integrated_out_fermions_at_MF': 3,
            'scalar_threshold_count': 1, 'integrated_out_scalars_at_MS': 2,
            'EFT_order': 'UV -> ALL F out -> EFT1 BOTH S -> BOTH S out -> SMEFT',
            'native_matchete_PL_kernels': len(native_labels),
            'EFT1_source_quartics': source_ledger,
            'EFT1_source_fermions': heavy_ledger,
            'native_source_quartics': native_sources,
            'EFT1_scalar_CT_vs_B0': scalar_check,
            'EFT1_scalar_tensor_purity_residual': float(ct['projected_tensor_purity_residual']),
            'candidate_ADM_seven_channel': {
                name: matrix_record(val, name, physical=True)
                for name, val in adm['CT_per_ADM_channel'].items()},
            'pole_diagnostics': pole_evidence,
            'existing_Batch6_comparison': old_comparison,
            'matching_scale_and_low_energy': observables,
            'physical_matching_certification': certification,
            'source_matrices_never_symmetrised_or_fitted': True,
            'complete_one_loop_UV_to_EFT1_matching': False,
            'complete_native_UV_counterterms': False,
            'low_energy_observables_conditional_not_physical_MSbar': True,
            'data_files_SHA256': {
                'source_poles_and_counterterms.npz': hash_file(npz_path),
                'matching_scale_scenario_matrices.npz': hash_file(scan_path),
                **({'conditional_mass_and_mixing_scan.csv': hash_file(csv_path)} if csv_path else {})},
        }
        save_json(ledger_path, ledger)
        return {
            'model': model, 'status': STATUS, 'stage': 'complete',
            'EFT1_scalar_CT_vs_B0_relative': scalar_check['CT_vs_B0_scale_slope_relative'],
            'native_vs_explicit_scalar_CT_relative': pole_evidence['native_vs_scalar_CT_relative'],
            'native_vs_candidate_seven_channel_ADM_relative': pole_evidence['native_vs_full_EFT1_ADM_relative'],
            'EFT1_scalar_CT_vs_ADM_scalar_pair_relative': pole_evidence['explicit_scalar_CT_vs_candidate_scalar_ADM_relative'],
            'EFT1_source_antisymmetric_cancellation_fraction': source_ledger['antisymmetric_cancellation_fraction'],
            'EFT1_source_quartic_count': len(names),
            'low_energy_running_computed': not skip_running,
            'mu_low_gev': mu_low, 'vev_gev': vev,
            'native_source_quartic_count': len(groups),
            'native_PL_kernel_count': len(native_labels),
            'nominal_analytic_masses_eV': observables['by_factor']['1.0']['analytic_triangle'].get('masses_eV'),
            'batch6_compare': old_comparison,
            'output_path': str(ledger_path.relative_to(root)),
            'output_sha256': hash_file(ledger_path),
        }
    except Exception as exc:
        row['error'] = f'{type(exc).__name__}: {exc}'
        print(model, 'BLOCKED', row['stage'], row['error'], file=sys.stderr)
        return row


def write_summary(out: Path, results: dict, *, selected: list[str],
                  mu_low: float, vev: float, skip_running: bool) -> dict:
    complete = sum(r.get('status') == STATUS for r in results.values())
    summary = {
        'schema': SCHEMA, 'status': 'NUMERICAL_AUDIT_NOT_PHYSICAL_MSbar',
        'completed_count': complete, 'blocked_count': len(MODELS)-complete,
        'results': results, 'executed_in_this_call': selected,
        'EFT_order': 'UV -> ALL F out -> EFT1 with both S -> both S out -> SMEFT',
        'full_one_loop_matching_complete': False,
        'native_UV_pole_cancellation_certified': False,
        'mass_results_physical': False,
        'mass_scan_run': not skip_running,
        'mu_low_gev': mu_low, 'vev_gev': vev,
    }
    save_json(out/'summary.json', summary)
    rows = []
    for model in MODELS:
        r = results.get(model, {'model': model, 'status': 'NOT_RUN'})
        rows.append({k: r.get(k) for k in (
            'model', 'status', 'EFT1_scalar_CT_vs_B0_relative',
            'native_vs_explicit_scalar_CT_relative',
            'native_vs_candidate_seven_channel_ADM_relative',
            'EFT1_scalar_CT_vs_ADM_scalar_pair_relative',
            'EFT1_source_antisymmetric_cancellation_fraction',
            'EFT1_source_quartic_count', 'native_source_quartic_count',
            'native_PL_kernel_count', 'nominal_analytic_masses_eV', 'error')})
    csv_path = out/'all16_source_pole_closure.csv'
    tmp = csv_path.with_suffix('.partial.csv')
    with tmp.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    tmp.replace(csv_path)
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--out', type=Path, default=DEST)
    parser.add_argument('--model', choices=MODELS)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--skip-low-energy-running', action='store_true')
    parser.add_argument('--mu-low', type=float, default=173.0)
    parser.add_argument('--vev', type=float, default=246.22)
    parser.add_argument('--show', action='store_true')
    args = parser.parse_args(argv)
    if args.resume and args.overwrite:
        raise ValueError('Use --resume or --overwrite, never both')
    if not (0 < args.mu_low <= MS and 0 < args.vev < 1000):
        raise ValueError('Invalid low renormalisation scale or Higgs vev')
    root = args.root.resolve()
    out = (root / args.out).resolve() if not args.out.is_absolute() else args.out.resolve()
    if not out.is_relative_to(root/'output'):
        raise ValueError('Batch7 output must remain inside local project output/')
    summary_path = out/'summary.json'
    old = read_json(summary_path) if summary_path.is_file() else None
    if old and old.get('schema') != SCHEMA:
        raise ValueError('Existing output schema differs; refuse to overwrite')
    if args.show:
        if not old:
            print('No Batch7 report:', summary_path)
            return 1
        for model in MODELS:
            r = old.get('results', {}).get(model, {})
            print(f'{model:6} {r.get("status", "NOT_RUN"):68} '
                  f'B0={r.get("EFT1_scalar_CT_vs_B0_relative", "-")} '
                  f'native_gap={r.get("native_vs_explicit_scalar_CT_relative", "-")} '
                  f'masses={r.get("nominal_analytic_masses_eV")}')
        print('completed=', old['completed_count'], 'blocked=', old['blocked_count'])
        return 0
    if args.overwrite and out.is_dir():
        shutil.rmtree(out)
    if out.is_dir() and any(out.iterdir()) and not args.resume and not args.overwrite:
        raise FileExistsError('Batch7 output exists: specify --resume or --overwrite')
    out.mkdir(parents=True, exist_ok=True)
    results = dict(old.get('results', {})) if old and not args.overwrite else {}
    selected = []
    requested = [args.model] if args.model else list(MODELS)
    for model in requested:
        r = results.get(model, {})
        p = out/tag(model)/'source_pole_and_matching_ledger.json'
        good = (r.get('status') == STATUS and p.is_file() and
                r.get('output_sha256') == hash_file(p) and
                (args.skip_low_energy_running or r.get('low_energy_running_computed') is True) and
                r.get('mu_low_gev') == args.mu_low and r.get('vev_gev') == args.vev)
        if args.resume and good:
            continue
        selected.append(model)
    print('Selected:', ', '.join(selected) or 'none (already completed)', flush=True)
    for model in selected:
        r = run_one(root, out, model, mu_low=args.mu_low,
                    vev=args.vev, skip_running=args.skip_low_energy_running)
        results[model] = r
        write_summary(out, results, selected=selected, mu_low=args.mu_low,
                      vev=args.vev, skip_running=args.skip_low_energy_running)
        print(model, r.get('status'),
              'CT/B0', r.get('EFT1_scalar_CT_vs_B0_relative'),
              'native-pole', r.get('native_vs_explicit_scalar_CT_relative'),
              'nominal-masses', r.get('nominal_analytic_masses_eV'), flush=True)
    report = write_summary(out, results, selected=selected, mu_low=args.mu_low,
                           vev=args.vev, skip_running=args.skip_low_energy_running)
    print('Batch7 summary:', summary_path,
          'complete=', report['completed_count'], 'blocked=', report['blocked_count'])
    return 0 if all(results.get(m, {}).get('status') == STATUS for m in requested) else 1


if __name__ == '__main__':
    sys.exit(main())
