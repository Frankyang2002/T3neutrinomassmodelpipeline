"""Regression tests for experimental C5 physics summary (not Mathematica tests)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.SummarizeGeneralT3Physics import (  # noqa: E402
    C5_CLASSES, EXPECTED_CLASSES, compile_summary,
)


def test_expected_model_sets():
    assert len(EXPECTED_CLASSES) == 16
    assert C5_CLASSES == {"B_m1", "C_m1"}


def test_missing_artifacts_do_not_pass(tmp_path):
    result = compile_summary(tmp_path / 'empty')
    assert not result['all_zero_extra_limits_formally_proven']
    assert not result['symbolically_consistent']
    assert all(row['zero_extra_proof_status'] == 'NotRun' for row in result['models'])


def test_formal_proof_status_cannot_be_inferred_from_text_comparison(tmp_path):
    folder = tmp_path / 'testroot'
    folder.mkdir()
    data = {'results': [
        {'model': 'B_m1', 'status': 'ZeroExtraLimitNotProven',
         'zero_extra_limit_proven': False},
        {'model': 'A_m4', 'status': 'TextIdentical',
         'zero_extra_limit_proven': True},
    ]}
    (folder / 'c5_algebraic_report.json').write_text(json.dumps(data))
    report = compile_summary(folder)
    by_model = {row['model']: row for row in report['models']}
    assert not by_model['B_m1']['zero_extra_limit_proven']
    assert by_model['A_m4']['zero_extra_limit_proven']
    assert not report['all_zero_extra_limits_formally_proven']


def test_existing_original_matcher_report_fields(tmp_path):
    root = tmp_path / 'trial'
    root.mkdir()
    original = {
        'results': [
            {'model': model,
             'baseline_extended_c5_text_equal_ignoring_whitespace':
                model not in C5_CLASSES,
             'variants': {'extended': {'extra_couplings_observed':
                ['lambdaT3ExtraQ020020I1'] if model in C5_CLASSES else []}}}
            for model in EXPECTED_CLASSES
        ]
    }
    (tmp_path / 'general_matching_symbolic_report.json').write_text(json.dumps(original))
    result = compile_summary(root)
    assert result['symbolically_consistent']
    assert set(result['c5_changed_models']) == C5_CLASSES
    assert not result['all_zero_extra_limits_formally_proven']
