"""Smoke-test the whole-benchmark continuation audit, using saved-format fixtures."""
import json
from pathlib import Path

from Numerical.diagnostics import VerifyFullContinuation as verify


def test_empty_input_reports_incomplete(tmp_path: Path):
    result = verify.audit(tmp_path / 'missing', tmp_path / 'report')
    assert result['status'] == 'INCOMPLETE_OR_FAILED'
    assert result['observed_cases'] == 0
    assert (tmp_path / 'report' / 'continuation_verification.csv').exists()


def test_bad_case_preserves_failure_reason(tmp_path: Path):
    root = tmp_path / 'input'
    path = root / 'smallY_smallL' / 'T3_A_alpha_m2' / 'data' / 'running_diagnostics.json'
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({'status':'Success'}), encoding='utf-8')
    result = verify.audit(root, tmp_path / 'report')
    assert result['observed_cases'] == 1
    assert result['failed_count'] == 1
    assert result['rows'][0]['explanation'].startswith('KeyError:')
