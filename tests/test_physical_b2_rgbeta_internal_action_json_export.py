"""Batch 4b: regression guards for fail-closed native RawJSON transport.

Static Wolfram guards do not execute installed RGBeta; local Python tests
validate malformed reports and the witness boundary. Native Wolfram evaluation
is required on the user's checkout.
"""
from __future__ import annotations

from pathlib import Path
import json
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'Numerical/lo_certification/physical_normalisation'))
import b2_rgbeta_internal_action as B
from b2_gauge_coupling import AuditBlocked

SCRIPT = ROOT / 'scripts/RunPhysicalB2RGBetaInternalActionProbe.wl'


def script() -> str:
    return SCRIPT.read_text(encoding='utf-8')


def test_nonjson_values_preflighted_before_export():
    s = script()
    assert 'B4JSONUnexpected[$B4Probe,"$"]' in s
    assert '"reason"->"NonJSONValue"' in s
    assert 'B4JSONSafe[$B4Probe]' in s
    assert '"status"->"NativeInternalTensorDataNotJSONSafe"' in s
    assert '"json_preflight_issues"->$B4JSONIssues' in s


def test_malformed_native_values_cannot_pass_status_guard():
    bogus = {
        'format': B.FORMAT,
        'status': 'NativeInternalTensorDataNotJSONSafe',
        'json_preflight_issues': [{'path': '$.git_head_at_capture',
                                   'inputform': 'Missing["GitCommandFailed"]'}],
        'physical_source_certified': False,
        'internal_prefactor_is_physical_action_factor': False,
        'numeric_component_axes_verified': False,
        'classes': [],
    }
    with pytest.raises(AuditBlocked, match='wrong native internal-witness format/status'):
        B.validate_native_data(bogus)


def test_rawjson_encoded_in_memory_before_disk_write():
    s = script()
    assert s.index('ExportString[$B4Safe,"RawJSON"]') < s.index('Export[$B4Tmp,$B4Encoded,"Text"')
    assert 'ImportString[$B4Encoded,"RawJSON"]' in s
    assert 'Import[$B4Tmp,"RawJSON"]' in s
    assert 'RenameFile[$B4Tmp,$B4Out]' in s


def test_partial_file_and_existing_output_fail_closed():
    s = script()
    assert 'If[FileExistsQ[$B4Out]' in s
    assert '$B4Tmp=$B4Out<>".partial"' in s
    assert 'If[FileExistsQ[$B4Tmp]' in s
    assert 'FileByteCount[$B4Tmp]<16' in s


def test_zero_byte_native_report_has_explanatory_failure(tmp_path):
    p = tmp_path / 'internal_action_registry.json'
    p.touch()
    with pytest.raises(AuditBlocked, match='empty native output'):
        B.read_report(p)


def test_malformed_json_has_explanatory_failure(tmp_path):
    p = tmp_path / 'internal_action_registry.json'
    p.write_text('{incomplete', encoding='utf-8')
    with pytest.raises(AuditBlocked, match='malformed/incomplete native JSON'):
        B.read_report(p)


def test_valid_json_reader_does_not_claim_native_verification(tmp_path):
    p = tmp_path / 'evidence.json'
    p.write_text(json.dumps({'physical_source_certified': False}), encoding='utf-8')
    assert B.read_report(p) == {'physical_source_certified': False}


def test_no_physical_certificate_and_six_expected_native_operators():
    s = script()
    for name in ('y1', 'y2', 'MF', 'lambdaT3', 'lambdaH1', 'lambdaH2'):
        assert f',{name},' in s
    assert '"physical_source_certified"->False' in s
    assert '"physical_Lagrangian_coefficient_verified"->False' in s
    assert '"numeric_component_axes_verified"->False' in s
