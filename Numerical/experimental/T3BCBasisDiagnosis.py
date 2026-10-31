"""Read-only diagnosis of B/C heavy-flavour basis covariance.

Does not alter the RGBeta exports, UV integrator, or C5 matching.
Requires the existing T3BCComplexUVStress implementation from the repository.
"""
from __future__ import annotations
import json
from pathlib import Path
from typing import Any
import numpy as np

EXPECTED_MODELS = ('B_-1', 'C_-1')


def analyse_report(report: dict[str, Any]) -> dict[str, Any]:
    if report.get('model') not in EXPECTED_MODELS:
        raise ValueError('Unexpected model')
    stress = report.get('native_RGBeta_heavy_basis_stress', {})
    if stress.get('status') != 'NativeRGBetaComplexHeavyBasisCovarianceComputed':
        raise ValueError('Covariance calculation absent or blocked')
    tests = stress['beta_tests']
    candidate = tests['source_CConj_orientation']['beta_relative_residual_by_parameter']
    alternative = tests['same_U_negative_control']['beta_relative_residual_by_parameter']
    keys = sorted(set(candidate) & set(alternative))
    if not {'MF', 'y1', 'y2', 'mS1Sq'}.issubset(keys):
        raise ValueError('Missing decisive beta blocks')
    threshold = 1e-6
    failed = {k: float(candidate[k]) for k in keys if float(candidate[k]) >= threshold}
    alternative_max = max(float(alternative[k]) for k in keys)
    return {
        'model': report['model'],
        'original_status': report['status'],
        'CConj_orientation_failed_blocks_at_1e-6': failed,
        'same_U_covariance_max_relative': alternative_max,
        'same_U_passes_at_1e-6': alternative_max < threshold,
        'CConj_orientation_passes_at_1e-6': not failed,
        'interpretation': (
            'RGBeta export covariant under y1@U,y2@U for this sampled state; '
            'CConj-based external change of basis is incompatible with at least one beta block. '
            'Native chiral vertex mapping remains uncertified.'),
        'changes_to_production_mapping': False,
        'physical_C5_certified': False,
    }


def diagnose_live(model: str, *, seed: int, integrate: bool = False) -> dict[str, Any]:
    from Numerical.experimental.T3BCComplexUVStress import run_complex_rge_covariance
    out = run_complex_rge_covariance(model, seed=seed, integrate=integrate)
    return analyse_report({
        'model': model, 'status': 'LiveCovarianceProbe',
        'native_RGBeta_heavy_basis_stress': out,
    })


def run(root: Path, output: Path, *, live: bool = False) -> dict[str, Any]:
    base = root / 'output/studies/t3_bc_mega_completion_v1'
    rows = {}
    for model in EXPECTED_MODELS:
        if live:
            result = diagnose_live(model, seed=135 if model == 'B_-1' else 136)
        else:
            filename = base / (model.replace('-', 'm') + '_mega.json')
            result = analyse_report(json.loads(filename.read_text(encoding='utf-8-sig')))
        rows[model] = result
    summary = {
        'format': 't3-bc-basis-diagnosis-v1',
        'status': 'RGBetaSameUMapCovariantButNativeSpinorDictionaryUncertified'
                  if all(r['same_U_passes_at_1e-6'] for r in rows.values())
                  else 'RGBetaSameUMapAlsoFails',
        'analysis_mode': 'live_original_rgbeta_rhs' if live else 'saved_native_rgbeta_results',
        'models': rows,
        'production_code_modified': False,
        'physical_matching_certified': False,
        'recommended_decisive_next_gate': 'Native Matchete four-spinor y1/y2 field-redefinition vertex equality; do not revise C5 before its result',
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    return summary
