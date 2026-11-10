"""Final deterministic report for the completed conditional T3 Batch-7 pipeline.

Read-only source inspection. No change to Batch 1/2/4/5/6/7 numerical files,
no new physics contribution, no fitting of the native pole. Source hashes and
numerical boundary identities must pass for all 16 models before publication.

EFT order: UV -> integrate out all three F -> EFT1 with S1 and S2 ->
integrate out both S -> SMEFT -> neutrino observables.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

SCHEMA = 't3-final-conditional-reproducibility-all16-batch8-v1'
B7_SCHEMA = 't3-native-source-pole-closure-batch7-v1'
B7_STATUS = 'NUMERICAL_SOURCE_RESOLVED_COUNTERTERM_AND_POLE_AUDIT_CONDITIONAL'
MODELS = ('A_-4', 'A_-2', 'A_+0', 'B_-3', 'B_-1', 'B_+1',
          'C_-3', 'C_-1', 'C_+1', 'D_-2', 'D_+0', 'D_+2',
          'E_-4', 'E_-2', 'E_+0', 'E_+2')
SCENARIOS = ('analytic_triangle', 'native_formal', 'native_as_UV_minus_soft_hypothesis')
FACTORS = (0.5, 1.0, 2.0)
MF, MS = 1e5, 1e3
B1 = Path('output/parallel_t3/fhard_triangle_all16_batch1_v1')
B2 = Path('output/parallel_t3/strict_one_loop_all16_batch2_v1')
B4 = Path('output/parallel_t3/native_finite_scenario_all16_batch4_v1')
B5 = Path('output/parallel_t3/hard_soft_matching_all16_batch5_v1')
B7 = Path('output/parallel_t3/native_pole_closure_all16_batch7_v1')
DEST = Path('output/parallel_t3/final_repro_all16_batch8_v1')
SOURCE_CODE = (
    'scripts/RunT3PoleClosureBatch7.py',
    'Numerical/parallel_t3/wilson/native_pole_closure_batch7.py',
    'scripts/RunT3CountertermBatch6.py',
    'Numerical/parallel_t3/wilson/msbar_counterterm_batch6.py',
    'Numerical/parallel_t3/wilson/matching_soft_subtraction_batch5.py',
    'Numerical/parallel_t3/wilson/quartic_channels_batch2.py',
    'Numerical/parallel_t3/threshold/matching.py',
    'Numerical/parallel_t3/threshold/pipeline.py',
    'Numerical/parallel_t3/threshold/rough_numerics.py',
    'Numerical/parallel_t3/background/contracts.py',
    'Numerical/parallel_t3/wilson/final_repro_batch8.py',
    'scripts/RunT3FinalRepro.py',
)


class ReproError(ValueError):
    """A saved output cannot be independently reproduced from its saved inputs."""


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def _load_json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise ReproError(f'Expected JSON object: {path}')
    return obj


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {k: np.asarray(data[k]).copy() for k in data.files}


def tag(model: str) -> str:
    if model not in MODELS:
        raise ReproError(f'Unknown model: {model}')
    return model.replace('-', 'm').replace('+', 'p')


def source_paths(root: Path, model: str) -> dict[str, Path]:
    """Exactly mirrors Batch 6's _source_filenames inspected on 2026-11-10."""
    slug = tag(model)
    return {
        'B1_hard': root / B1 / 'hard_fermion' / slug / 'fermion_hard_triangle.npz',
        'B1_hard_json': root / B1 / 'hard_fermion' / slug / 'fermion_hard_triangle.json',
        'B1_MF_quartic': root / 'output/parallel_t3/rough_inputs/fhard_triangle_all16_batch1_v1/a' / slug / 'background.npz',
        'B2_strict': root / B2 / 'inputs' / slug / 'strict_wilson_channels.npz',
        'B2_strict_json': root / B2 / 'inputs' / slug / 'strict_wilson_channels.json',
        'B2_MS_C5': root / B2 / 'results' / slug / 'C5_smeft_at_MS.npz',
        'B4_native': root / B4 / 'inputs' / slug / 'native_MF_LLHH_finite_source_candidate.npz',
        'B4_native_json': root / B4 / 'inputs' / slug / 'native_MF_LLHH_finite_source_candidate.json',
        'B4_MS_C5': root / B4 / 'results' / slug / 'C5_smeft_at_MS.npz',
        'B5_pole': root / B5 / 'inputs' / slug / 'hard_soft_and_native_Laurent.npz',
        'B5_pole_json': root / B5 / 'inputs' / slug / 'hard_soft_and_native_Laurent.json',
        'B5_MS_C5': root / B5 / 'results' / slug / 'C5_smeft_at_MS.npz',
        'B5_background': root / B5 / 'inputs' / slug / 'background.npz',
        'B5_background_json': root / B5 / 'inputs' / slug / 'background.json',
    }


def relative(a: Any, b: Any) -> float:
    x, y = np.asarray(a), np.asarray(b)
    if x.shape != y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ReproError('Nonfinite or incompatible comparison tensors')
    return float(np.linalg.norm(x - y) / max(float(np.linalg.norm(y)), 1e-100))


def c5(x: Any, label: str) -> np.ndarray:
    a = np.asarray(x, dtype=complex)
    if a.shape != (3, 3) or not np.isfinite(a).all():
        raise ReproError(f'{label}: require finite complex 3x3')
    if relative(a, a.T) > 1e-8:
        raise ReproError(f'{label}: total C5 is not flavour symmetric')
    return a  # no silent projection or symmetrisation


def _require_close(label: str, actual: Any, expected: Any,
                   tolerance: float = 2e-8) -> None:
    difference = relative(actual, expected)
    if difference > tolerance:
        raise ReproError(f'{label} mismatch (relative {difference:.6g})')


def _check_finite_number(val: Any, label: str) -> float:
    if isinstance(val, bool) or val is None:
        raise ReproError(f'{label} missing or not numeric')
    try:
        n = float(val)
    except (ValueError, TypeError) as e:
        raise ReproError(f'{label} not numeric') from e
    if not math.isfinite(n):
        raise ReproError(f'{label} nonfinite')
    return n


def check_source_hashes(root: Path, model: str, ledger: dict) -> dict:
    recorded = ledger.get('source_SHA256')
    paths = source_paths(root, model)
    if not isinstance(recorded, dict) or set(recorded) != set(paths):
        raise ReproError(f'{model}: missing or unexpected Batch-7 source hash keys')
    result = {}
    for key, p in paths.items():
        if not p.is_file():
            raise ReproError(f'{model}: missing source file {p}')
        actual = digest(p)
        if recorded[key] != actual:
            raise ReproError(f'{model}: source SHA256 mismatch: {p}')
        result[key] = {'path': str(p.relative_to(root)).replace('\\', '/'), 'sha256': actual}
    return result


def _check_artifacts(root: Path, folder: Path, row: dict, ledger: dict) -> dict:
    expected = ledger.get('data_files_SHA256')
    needed = {'source_poles_and_counterterms.npz',
              'matching_scale_scenario_matrices.npz',
              'conditional_mass_and_mixing_scan.csv'}
    if not isinstance(expected, dict) or set(expected) != needed:
        raise ReproError('Expected the complete Batch-7 9-scenario low-energy data files')
    hashes = {}
    for name in sorted(needed):
        p = folder / name
        if not p.is_file() or digest(p) != expected[name]:
            raise ReproError(f'Batch-7 data file missing or hash mismatch: {p}')
        hashes[name] = {'path': str(p.relative_to(root)).replace('\\', '/'),
                        'sha256': expected[name]}
    ledger_path = folder / 'source_pole_and_matching_ledger.json'
    if not ledger_path.is_file() or digest(ledger_path) != row.get('output_sha256'):
        raise ReproError(f'Batch-7 per-model ledger hash mismatch: {ledger_path}')
    hashes['source_pole_and_matching_ledger.json'] = {
        'path': str(ledger_path.relative_to(root)).replace('\\', '/'),
        'sha256': digest(ledger_path)}
    return hashes


def _check_observation(model: str, scenario: str, factor: float,
                       record: dict, low: np.ndarray, vev: float) -> dict:
    if record.get('scenario') != scenario or float(record.get('subtraction_scale_factor', 0)) != factor:
        raise ReproError(f'{model}: wrong scenario/scale record')
    if record.get('mass_results_conditional_not_physical_MSbar') is not True:
        raise ReproError(f'{model}: scenario result wrongly labelled physical')
    if record.get('fixed_MF_source_couplings') is not True:
        raise ReproError(f'{model}: non-frozen boundary input changed')
    mass = np.asarray(record.get('masses_eV'), dtype=float)
    if mass.shape != (3,) or not np.isfinite(mass).all() or np.any(mass < 0):
        raise ReproError(f'{model}: missing/nonfinite masses')
    derived = np.sort(np.linalg.svd(-0.5 * vev**2 * low, compute_uv=False))*1e9
    _require_close(f'{model} {scenario} {factor} masses vs low C5', mass, derived, 1e-7)
    dm21 = _check_finite_number(record.get('dm21_sq_eV2'), 'dm21')
    dm31 = _check_finite_number(record.get('dm31_sq_eV2'), 'dm31')
    if not math.isclose(dm21, float(mass[1]**2-mass[0]**2), rel_tol=1e-7, abs_tol=1e-8):
        raise ReproError(f'{model}: delta m21 squared disagrees with saved masses')
    if not math.isclose(dm31, float(mass[2]**2-mass[0]**2), rel_tol=1e-7, abs_tol=1e-8):
        raise ReproError(f'{model}: delta m31 squared disagrees with saved masses')
    mixed = record.get('mixing_defined')
    if not isinstance(mixed, bool):
        raise ReproError(f'{model}: mixing_defined missing')
    angles = {}
    for key in ('sin2_theta12', 'sin2_theta13', 'sin2_theta23'):
        val = record.get(key)
        if mixed:
            n = _check_finite_number(val, key)
            if not 0 <= n <= 1 + 1e-10:
                raise ReproError(f'{model}: invalid {key}')
            angles[key] = n
        else:
            if val is not None:
                # Some code may return partial angles; do not infer missing entries.
                angles[key] = _check_finite_number(val, key)
            else:
                angles[key] = None
    return {'m1_eV': float(mass[0]), 'm2_eV': float(mass[1]),
            'm3_eV': float(mass[2]), 'sum_m_eV': float(mass.sum()),
            'dm21_sq_eV2': dm21, 'dm31_sq_eV2': dm31,
            'mixing_defined': mixed, **angles}


def collect(root: Path) -> tuple[list[dict], list[dict], dict[str, np.ndarray], dict]:
    """Fail closed; no report is written until every saved model passes."""
    root = Path(root).resolve()
    summary_path = root / B7 / 'summary.json'
    summary = _load_json(summary_path)
    if (summary.get('schema') != B7_SCHEMA or summary.get('completed_count') != 16 or
            summary.get('blocked_count') != 0 or summary.get('mass_scan_run') is not True):
        raise ReproError('Batch 7 must report 16 complete, 0 blocked, full low-energy running')
    sr = summary.get('results', {})
    if not isinstance(sr, dict) or set(sr) != set(MODELS):
        raise ReproError('Batch-7 summary has missing, duplicated or unexpected models')
    vev = _check_finite_number(summary.get('vev_gev'), 'vev_gev')
    mu_low = _check_finite_number(summary.get('mu_low_gev'), 'mu_low_gev')
    if not (0 < mu_low <= MS and 0 < vev < 1000):
        raise ReproError('Invalid saved final scale or vev')
    output_rows, overview, all_matrices = [], [], {}
    provenance = {}
    for model in MODELS:
        row = sr[model]
        if (not isinstance(row, dict) or row.get('model') != model or
                row.get('status') != B7_STATUS or row.get('low_energy_running_computed') is not True):
            raise ReproError(f'{model}: incomplete model status or missing low-energy running')
        slug = tag(model)
        folder = root / B7 / slug
        ledger = _load_json(folder / 'source_pole_and_matching_ledger.json')
        if ledger.get('schema') != B7_SCHEMA or ledger.get('model') != model or ledger.get('status') != B7_STATUS:
            raise ReproError(f'{model}: ledger schema, model, or status changed')
        if (ledger.get('EFT_order') != 'UV -> ALL F out -> EFT1 BOTH S -> BOTH S out -> SMEFT'
                or ledger.get('integrated_out_fermions_at_MF') != 3
                or ledger.get('integrated_out_scalars_at_MS') != 2
                or ledger.get('fermion_threshold_count') != 1
                or ledger.get('scalar_threshold_count') != 1):
            raise ReproError(f'{model}: incorrect simultaneous fermion/scalar EFT sequence')
        for k, val in [('MF_GeV', MF), ('MS_GeV', MS), ('mu_low_GeV', mu_low), ('vev_GeV', vev)]:
            if not math.isclose(_check_finite_number(ledger.get(k), k), val, rel_tol=1e-12):
                raise ReproError(f'{model}: inconsistent scale or vev: {k}')
        source_hashes = check_source_hashes(root, model, ledger)
        data_hashes = _check_artifacts(root, folder, row, ledger)
        hard = c5(_load_npz(source_paths(root, model)['B1_hard'])['C5_hard_F_triangle_MF'], 'B1 hard')
        strict = c5(_load_npz(source_paths(root, model)['B2_strict'])['C5_EFT1_strict_running_delta'], 'B2 running')
        bubble = c5(_load_npz(source_paths(root, model)['B5_MS_C5'])['C5_scalar_bubble_MS'], 'B5 scalar')
        refs = {
            'analytic_triangle': _load_npz(source_paths(root, model)['B2_MS_C5'])['C5_SMEFT_MS'],
            'native_formal': _load_npz(source_paths(root, model)['B4_MS_C5'])['C5_SMEFT_MS'],
            'native_as_UV_minus_soft_hypothesis': _load_npz(source_paths(root, model)['B5_MS_C5'])['C5_SMEFT_MS'],
        }
        npz = _load_npz(folder / 'source_poles_and_counterterms.npz')
        ct = c5(npz['C5_EFT1_scalar_total'], model+' explicit scalar CT')
        b0 = c5(npz['C5_B0_determined_scalar_CT'], model+' B0 CT')
        native = c5(npz['C5_native_raw_pole_total'], model+' native pole')
        _require_close(model+' scalar CT vs B0', ct, b0, 1e-8)
        _require_close(model+' native gap metric',
                       np.asarray([relative(native, ct)]),
                       np.asarray([float(row['native_vs_explicit_scalar_CT_relative'])]), 1e-8)
        _require_close(model+' CT/B0 metric',
                       np.asarray([relative(ct, b0)]),
                       np.asarray([float(row['EFT1_scalar_CT_vs_B0_relative'])]), 1e-5)
        scan = _load_npz(folder / 'matching_scale_scenario_matrices.npz')
        scans = ledger.get('matching_scale_and_low_energy', {}).get('by_factor')
        if not isinstance(scans, dict) or set(scans) != {'0.5', '1.0', '2.0'}:
            raise ReproError(f'{model}: incomplete 3-scale output')
        model_rows = []
        for factor in FACTORS:
            group = scans[str(factor)]
            if not isinstance(group, dict) or set(group) != set(SCENARIOS):
                raise ReproError(f'{model}: incomplete scenario inventory at {factor}')
            for scenario in SCENARIOS:
                key = f'{scenario}_scale_{factor:g}'
                upper = c5(scan[key+'_C5_SMEFT_MS'], model+key+' C5 at MS')
                lower = c5(scan[key+'_C5_at_low'], model+key+' C5 at low scale')
                record = group[scenario]
                observables = _check_observation(model, scenario, factor, record, lower, vev)
                _require_close(model+key+' stored C5 norm',
                               np.asarray([float(record['C5_SMEFT_MS_norm_GeV_inverse'])]),
                               np.asarray([float(np.linalg.norm(upper))]), 1e-8)
                if factor == 1.0:
                    _require_close(model+' '+scenario+' nominal vs saved B2/4/5',
                                   upper, refs[scenario], 2e-8)
                    if scenario == 'analytic_triangle':
                        _require_close(model+' analytic final = hard + running + scalar bubble',
                                       upper, hard+strict+bubble, 2e-8)
                full = {
                    'model': model, 'scenario': scenario, 'scale_factor': factor,
                    'mu_subtraction_GeV': factor*MF,
                    'MF_GeV': MF, 'MS_GeV': MS, 'mu_low_GeV': mu_low, 'vev_GeV': vev,
                    **observables,
                    'C5_hard_F_triangle_MF_fro_GeVinv': float(np.linalg.norm(hard)),
                    'C5_EFT1_strict_running_delta_fro_GeVinv': float(np.linalg.norm(strict)),
                    'C5_scalar_bubble_MS_fro_GeVinv': float(np.linalg.norm(bubble)),
                    'C5_EFT1_scalar_CT_fro_GeVinv': float(np.linalg.norm(ct)),
                    'C5_native_raw_pole_fro_GeVinv': float(np.linalg.norm(native)),
                    'C5_SMEFT_MS_fro_GeVinv': float(np.linalg.norm(upper)),
                    'C5_low_fro_GeVinv': float(np.linalg.norm(lower)),
                    'CT_vs_B0_relative': float(row['EFT1_scalar_CT_vs_B0_relative']),
                    'native_vs_EFT1_scalar_CT_relative': float(row['native_vs_explicit_scalar_CT_relative']),
                    'native_vs_candidate_ADM_relative': float(row['native_vs_candidate_seven_channel_ADM_relative']),
                    'native_PL_kernel_count': int(row['native_PL_kernel_count']),
                    'EFT1_quartic_source_count': int(row['EFT1_source_quartic_count']),
                    'prescription': 'CONDITIONAL_NOT_PHYSICAL_MSbar',
                }
                output_rows.append(full)
                model_rows.append(full)
                all_matrices[f'{slug}__{key}__C5_SMEFT_MS'] = upper
                all_matrices[f'{slug}__{key}__C5_low'] = lower
        nominal = next(v for v in model_rows if v['scale_factor'] == 1.0 and v['scenario']=='analytic_triangle')
        _require_close(model+' nominal analytic masses vs Batch-7 summary',
                       np.asarray(row['nominal_analytic_masses_eV'], float),
                       np.asarray([nominal['m1_eV'], nominal['m2_eV'], nominal['m3_eV']]), 1e-8)
        overview.append({k: nominal[k] for k in (
            'model','m1_eV','m2_eV','m3_eV','sum_m_eV','dm21_sq_eV2','dm31_sq_eV2',
            'mixing_defined','sin2_theta12','sin2_theta13','sin2_theta23',
            'CT_vs_B0_relative','native_vs_EFT1_scalar_CT_relative',
            'native_vs_candidate_ADM_relative','native_PL_kernel_count',
            'EFT1_quartic_source_count','prescription')})
        all_matrices[f'{slug}__C5_hard_F_triangle_MF'] = hard
        all_matrices[f'{slug}__C5_EFT1_running_delta'] = strict
        all_matrices[f'{slug}__C5_scalar_bubble_MS'] = bubble
        all_matrices[f'{slug}__C5_EFT1_scalar_CT'] = ct
        all_matrices[f'{slug}__C5_native_raw_pole'] = native
        provenance[model] = {'source_files': source_hashes, 'batch7_data_files': data_hashes,
                             'batch7_summary_row': row,
                             'native_pole_diagnostic_only': True}
    if len(output_rows) != 16*3*3 or len(overview) != 16:
        raise ReproError('Incomplete 144-branch model inventory')
    detail = {'batch7_summary': {
        'path': str((B7/'summary.json')).replace('\\','/'),
        'sha256': digest(summary_path), 'schema': summary['schema'],
        'completed_count': summary['completed_count'],
        'blocked_count': summary['blocked_count']},
        'provenance': provenance, 'mu_low_GeV': mu_low, 'vev_GeV': vev}
    return output_rows, overview, all_matrices, detail


def _atomic_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ReproError('Cannot write empty report table')
    tmp = path.with_name(path.name+'.partial')
    with tmp.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    tmp.replace(path)


def _atomic_json(path: Path, record: dict) -> None:
    tmp = path.with_name(path.name+'.partial')
    tmp.write_text(json.dumps(record, sort_keys=True, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    tmp.replace(path)


def process_document() -> str:
    return '''# Final T3 numerical benchmark process (conditional prescription)\n\n1. Saved UV/fermion matching inputs: completed Batch 1 source-resolved heavy-F triangle.\n2. All three heavy fermions are integrated out at the common MF = 100000 GeV.\n3. EFT1 retains both heavy scalar multiplets; strict one-loop Wilson channels are evolved to MS = 1000 GeV.\n4. Both scalars are integrated out at MS using the existing full-tensor scalar bubble threshold.\n5. The resulting SMEFT C5 evolves to the saved low-energy scale and gives mnu = -v^2 C5/2.\n6. Batch 7 retains named native Laurent poles, explicit EFT1 scalar counterterms, source matrices, and three conditional matching prescriptions over three subtraction scales.\n7. This report validates all 16 model outputs, 14 saved prerequisite source hashes per model, per-model data hashes, nominal threshold identities and masses derived again from stored C5(low). No source is recalculated or altered during collection.\n\n**Status:** Complete and reproducible *conditional numerical benchmark*. No independent bare-UV amplitude, fitted missing counterterm, full physical MSbar matching, or experimental fit is asserted.\n\n**Scenarios:** `analytic_triangle` (the primary conditional source-triangle result); `native_formal` (saved Matchete pole-dropped prescription); `native_as_UV_minus_soft_hypothesis` (comparison hypothesis). Subtraction-scale variation uses MF-frozen couplings and is *not* a new physical threshold.\n\n**Eigenvalues:** m1 <= m2 <= m3 are ascending Takagi singular values; they are not a normal/inverted-ordering identification. Neutrino masses are in eV, C5 in GeV^-1.\n\n**Outputs:** `all16_scenarios.csv` (144 rows), `all16_nominal.csv` (48 rows, nominal subtraction scale), `all16_overview.csv` (16 primary nominal rows), `all16_wilson_matrices.npz` (complex 3x3 Wilson matrices), `manifest.json` (provenance and hashes).\n'''


def write_report(root: Path, out: Path = DEST) -> dict:
    root = Path(root).resolve()
    out = (root / out).resolve() if not Path(out).is_absolute() else Path(out).resolve()
    allowed = (root / 'output/parallel_t3').resolve()
    if out == allowed or not out.is_relative_to(allowed):
        raise ReproError('Final report must be a separate folder inside output/parallel_t3/')
    if any(out.is_relative_to((root / src).resolve()) or
           (root / src).resolve().is_relative_to(out) for src in (B1, B2, B4, B5, B7)):
        raise ReproError('Final reporting must not overwrite or be inside a source batch directory')
    # All validation is performed before writing the first output byte.
    rows, overview, matrices, inputs = collect(root)
    code_hashes = {}
    for name in SOURCE_CODE:
        p = root/name
        if not p.is_file():
            raise ReproError(f'Missing installed source file: {name}')
        code_hashes[name] = digest(p)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        'all16_scenarios.csv': out/'all16_scenarios.csv',
        'all16_nominal.csv': out/'all16_nominal.csv',
        'all16_overview.csv': out/'all16_overview.csv',
        'all16_wilson_matrices.npz': out/'all16_wilson_matrices.npz',
        'PROCESS.md': out/'PROCESS.md',
    }
    _atomic_csv(paths['all16_scenarios.csv'], rows)
    _atomic_csv(paths['all16_nominal.csv'], [r for r in rows if r['scale_factor'] == 1.0])
    _atomic_csv(paths['all16_overview.csv'], overview)
    tmp = paths['all16_wilson_matrices.npz'].with_name('all16_wilson_matrices.partial.npz')
    np.savez_compressed(tmp, **matrices)
    tmp.replace(paths['all16_wilson_matrices.npz'])
    proc_tmp = paths['PROCESS.md'].with_name('PROCESS.md.partial')
    proc_tmp.write_text(process_document(), encoding='utf-8')
    proc_tmp.replace(paths['PROCESS.md'])
    manifest = {
        'schema': SCHEMA, 'status': 'COMPLETE_REPRODUCIBLE_CONDITIONAL_NUMERICAL_BENCHMARK',
        'completed_models': 16, 'blocked_models': 0, 'total_scenario_rows': len(rows),
        'nominal_scenario_rows': 48, 'primary_nominal_rows': len(overview),
        'model_order': list(MODELS), 'scenario_order': list(SCENARIOS),
        'subtraction_scale_factors': list(FACTORS),
        'EFT_order': 'UV -> ALL F out -> EFT1 BOTH S -> BOTH S out -> SMEFT',
        'integration_scales_GeV': {'MF': MF, 'MS': MS, 'mu_low': inputs['mu_low_GeV']},
        'vev_GeV': inputs['vev_GeV'],
        'source_code_SHA256': code_hashes,
        'saved_inputs_and_batch7_outputs': inputs,
        'report_file_SHA256': {name: digest(p) for name, p in paths.items()},
        'python_version': sys.version.split()[0], 'numpy_version': np.__version__,
        'all_file_paths_relative_to_project_root': True,
        'existing_batch_outputs_unchanged': True,
        'physical_MSbar_matching_claimed': False,
        'native_pole_fitted_as_counterterm': False,
        'conditional_results_not_physical_mass_predictions': True,
    }
    _atomic_json(out/'manifest.json', manifest)
    return manifest


def verify_report(root: Path, out: Path = DEST) -> dict:
    root = Path(root).resolve()
    out = (root / out).resolve() if not Path(out).is_absolute() else Path(out).resolve()
    manifest = _load_json(out/'manifest.json')
    if manifest.get('schema') != SCHEMA or manifest.get('completed_models') != 16:
        raise ReproError('Final report schema or model count mismatch')
    for name, expected in manifest['report_file_SHA256'].items():
        if Path(name).name != name or digest(out/name) != expected:
            raise ReproError(f'Final report content changed: {name}')
    hashes = manifest.get('source_code_SHA256', {})
    if set(hashes) != set(SOURCE_CODE):
        raise ReproError('Recorded source code manifest is incomplete')
    for name, recorded_hash in hashes.items():
        p = root/name
        if not p.is_file() or digest(p) != recorded_hash:
            raise ReproError(f'Installed physics/report code differs from recorded run: {name}')
    # Re-read all 14 prerequisite source hashes and three saved outcomes per model.
    # This prevents a stale manifest from appearing verified after recomputation.
    _, _, _, current_inputs = collect(root)
    if current_inputs != manifest.get('saved_inputs_and_batch7_outputs'):
        raise ReproError('Batch-7 source or results changed since final report generation')
    return manifest
