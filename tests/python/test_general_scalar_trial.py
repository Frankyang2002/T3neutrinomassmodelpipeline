"""Static/regression tests for the opt-in general T3 scalar-quartic trial.

These tests require no Mathematica. They do NOT establish a successful
Matchete or RGBeta matching calculation; that requires local Wolfram.
"""

from __future__ import annotations

import json
from pathlib import Path

from check_t3_scalar_quartics import audit, conjugate

PROJECT = Path(__file__).resolve().parents[2]
EXTRA = PROJECT / "Lagrangian/interactions/GeneralScalarOperators.wl"
TOPOLOGY = PROJECT / "Lagrangian/interactions/T3Topology.wl"
MANIFEST = PROJECT / "docs/general_scalar_trial_operator_manifest.json"
AUDIT_WL = PROJECT / "Lagrangian/RunGeneralScalarAudit.wl"

EXPECTED_EXTRA_COMPLEX_QUARTICS = {
    "A_-4": 0, "A_-2": 5, "A_+0": 5,
    "B_-3": 1, "B_-1": 7, "B_+1": 1,
    "C_-3": 1, "C_-1": 7, "C_+1": 1,
    "D_-2": 5, "D_+0": 5, "D_+2": 0,
    "E_-4": 0, "E_-2": 6, "E_+0": 6, "E_+2": 0,
}


def test_all_sixteen_model_operator_counts():
    data = audit()
    assert set(data) == set(EXPECTED_EXTRA_COMPLEX_QUARTICS)
    for key, expected in EXPECTED_EXTRA_COMPLEX_QUARTICS.items():
        entry = data[key]
        missing = entry["quartic_missing_matchete"]
        assert sum(v["invariants"] for v in missing) == expected, key
        assert all(v["complex_coefficient"] for v in missing), key
        assert entry["quartic_missing_matchete_real_parameters"] == 2 * expected
        assert len(entry["dimension_two_missing"]) <= 1


def test_manifest_matches_audit():
    data = audit()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert manifest["model_count"] == 16
    for model, body in manifest["models"].items():
        expected = data[model]
        assert body["extra_quartic_complex_couplings"] == sum(
            q["invariants"] for q in expected["quartic_missing_matchete"]
        )
        assert body["quadratic_complex_couplings_not_enabled"] == sum(
            q["invariants"] for q in expected["dimension_two_missing"]
        )
        for quartic in body["extra_quartic_structures"]:
            assert tuple(quartic["counts"]) != conjugate(tuple(quartic["counts"]))


def test_runtime_opt_in_does_not_change_default():
    ext = EXTRA.read_text(encoding="utf-8")
    topo = TOPOLOGY.read_text(encoding="utf-8")
    assert 'Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] === "1"' in ext
    assert "T3ExtendedBuildExtraQuarticCandidates[model]" in topo
    assert "TrueQ[T3ExperimentalGeneralQuarticsQ[]]" in topo
    assert "BuildT3YukawaCandidates[model, 1]" in topo
    assert "BuildT3YukawaCandidates[model, 2]" in topo
    assert "BuildGeneralScalarPotentialCandidates[model]" in topo
    assert "BuildT3MixingCandidates[model]" in topo
    assert "T3ExtendedDescribeQuadraticOperators" in ext
    assert "T3ExtendedBuildExtraQuarticCandidates" in ext
    assert "T3ExtendedSymmetrizeGroup" in ext
    assert "T3ExtendedConjugateCounts" in ext


def test_symbolic_names_no_underscore_in_generated_symbols():
    ext = EXTRA.read_text(encoding="utf-8")
    assert 'key = StringJoin[ToString /@ counts]' in ext
    assert 'Symbol["lambdaT3Extra"' in ext
    assert 'Symbol["T3ExtraCG"' in ext


def test_wolfram_delimiters_balanced():
    # Lexer ignores WL comments (nested) and quoted strings. This is only a
    # syntax smoke check, not a full parser or Mathematica evaluation.
    for path in (EXTRA, TOPOLOGY, AUDIT_WL,
                 PROJECT / 'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl'):
        text = path.read_text(encoding="utf-8")
        stack = []
        comments = 0
        in_string = False
        i = 0
        pairs = {")": "(", "]": "[", "}": "{", "|>": "<|"}
        while i < len(text):
            chunk = text[i:i + 2]
            if comments:
                if chunk == "(*":
                    comments += 1
                    i += 2
                    continue
                if chunk == "*)":
                    comments -= 1
                    i += 2
                    continue
                i += 1
                continue
            if in_string:
                if text[i] == "\\":
                    i += 2
                    continue
                if text[i] == '"':
                    in_string = False
                i += 1
                continue
            if chunk == "(*":
                comments += 1
                i += 2
                continue
            if text[i] == '"':
                in_string = True
                i += 1
                continue
            if chunk in ("<|", "|>"):
                if chunk == "<|":
                    stack.append(chunk)
                else:
                    assert stack.pop() == pairs[chunk], path
                i += 2
                continue
            if text[i] in "([{":
                stack.append(text[i])
            elif text[i] in ")]}":
                assert stack and stack.pop() == pairs[text[i]], (path, i)
            i += 1
        assert not comments and not in_string and not stack, (path, stack)


def test_validation_script_excludes_matching_and_rg():
    wl = AUDIT_WL.read_text(encoding="utf-8")
    assert "T3ExtendedBuildExtraQuarticCandidates[model]" in wl
    assert "ValidateT3Candidate" in wl
    assert "CheckLagrangian[combinedLag]" in wl
    assert "RunSequentialT3Matching[" not in wl
    assert "RGBetaExtended\" -> False" in wl


def test_audit_frontend_scope_and_error_separation():
    audit_text = AUDIT_WL.read_text(encoding="utf-8")
    assert 'UsingFrontEnd[T3TrialAuditRun[trialModel]]' in audit_text
    assert 'feReady = MatchQ[$FrontEnd, _FrontEndObject]' in audit_text
    assert '"Status" -> "FrontendUnavailable"' in audit_text
    assert '"AcceptedExtraQuarticCount" -> Null' in audit_text
    assert '"Status" -> "BaselineValidationFailed"' in audit_text
    assert 'CheckLagrangian[baselineLag]' in audit_text
    assert 'CheckLagrangian[combinedLag]' in audit_text
    assert 'Total[Lookup[extra, "Expression", {}]]' in audit_text


def test_quartic_builder_fails_loudly_on_missing_invariants():
    text = EXTRA.read_text(encoding="utf-8")
    assert '"ERROR: SU(2) invariant basis failed for "' in text
    assert '"ERROR: CG or coupling registration failed for "' in text
    assert 'If[data === $Failed,' in text
    assert 'Return[$Failed]' in text


def test_python_audit_report_success_failure(tmp_path):
    import sys
    sys.path.insert(0, str(PROJECT / 'scripts'))
    from ReportGeneralScalarTrial import evaluate
    root = tmp_path / 'general'
    path = root / 'B_m1' / 'audit.json'
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({
        'Status': 'FrontendUnavailable',
        'ExtraQuarticCandidateCount': None,
        'AcceptedExtraQuarticCount': None,
        'BaselineLagrangianValid': None,
        'CombinedLagrangianValid': None,
    }))
    report = evaluate(root)
    assert report['model_count_checked'] == 1
    assert report['failed_or_unrun'] == 1
    assert report['models']['B_-1']['status'] == 'FrontendUnavailable'
    path.write_text(json.dumps({
        'Status': 'Success',
        'ExtraQuarticCandidateCount': 7,
        'AcceptedExtraQuarticCount': 7,
        'BaselineLagrangianValid': True,
        'CombinedLagrangianValid': True,
    }))
    report = evaluate(root)
    assert report['passed'] == 1
    path.write_text(json.dumps({
        'Status': 'Success',
        'ExtraQuarticCandidateCount': 6,
        'AcceptedExtraQuarticCount': 6,
        'BaselineLagrangianValid': True,
        'CombinedLagrangianValid': True,
    }))
    report = evaluate(root)
    assert report['passed'] == 0
    assert 'disagrees' in report['models']['B_-1']['diagnostic']


def test_matching_report_detects_new_c5_coupling(tmp_path):
    import sys
    sys.path.insert(0, str(PROJECT / 'scripts'))
    from InspectGeneralScalarMatching import inspect_model
    root = tmp_path / 'B_m1'
    (root / 'baseline').mkdir(parents=True)
    (root / 'extended').mkdir(parents=True)
    (root / 'baseline' / 'c5_coefficient.txt').write_text('lambdaT3[]')
    (root / 'extended' / 'c5_coefficient.txt').write_text(
        'lambdaT3[] + lambdaT3ExtraQ000202I1[]')
    value = inspect_model(root)
    assert value['status'] == 'ExtendedC5DependsOnNewCouplings'
    assert value['variants']['extended']['extra_couplings_observed'] == [
        'lambdaT3ExtraQ000202I1']
    assert value['variants']['baseline']['extra_couplings_observed'] == []


def test_batch_script_restores_environment_and_uses_separate_outputs():
    script = (PROJECT / 'scripts/RunGeneralScalarTrial.ps1').read_text(encoding='utf-8')
    assert "[ValidateSet('B','BC','Problems','All')]" in script
    assert "[ValidateSet('Audit','Matching','Both')]" in script
    assert 'finally {' in script
    assert "SetEnvironmentVariable('T3_EXPERIMENTAL_GENERAL_QUARTICS', $originalFlag" in script
    assert "foreach ($variant in @('baseline','extended'))" in script
    assert 'trial_runs.json' in script


def test_batch_script_windows_powershell_scoped_variable_bug():
    import re
    source = (PROJECT / 'scripts/RunGeneralScalarTrial.ps1').read_text(encoding='utf-8')
    assert 'Skipping matching for ${name}: symbolic validation failed.' in source
    # A bare $name: inside a string is a parse error under Windows PowerShell.
    assert not re.search(r'(?<!\{)\$name:', source)


def test_reporter_explicitly_handles_empty_audit_folder(tmp_path, capsys, monkeypatch):
    import sys
    sys.path.insert(0, str(PROJECT / 'scripts'))
    import ReportGeneralScalarTrial as report_module
    monkeypatch.setattr(sys, 'argv', [
        'ReportGeneralScalarTrial.py', '--root', str(tmp_path / 'nothing')
    ])
    assert report_module.main() == 1
    stdout = capsys.readouterr().out
    assert 'No audit JSON files found' in stdout


def test_verified_wolfram_catalog_matches_static_manifest():
    import sys
    sys.path.insert(0, str(PROJECT / 'scripts'))
    from GenerateGeneralScalarWLManifest import render
    catalog = PROJECT / 'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl'
    assert catalog.read_text(encoding='utf-8') == render()
    assert '"2_2_-3" -> {{{0, 2, 0, 0, 0, 2}, 1}}' in render()
    assert '"2_2_1" -> {{{0, 2, 2, 0, 0, 0}, 1}}' in render()
    assert '"2_2_-1" ->' in render()


def test_generator_requires_complete_su2_operator_basis():
    source = EXTRA.read_text(encoding='utf-8')
    assert 'GeneralScalarVerifiedQuartics.wl' in source
    assert 'T3ExtendedVerifiedEntries[model]' in source
    assert 'T3ExtendedVerifiedStructureCount[model]' in source
    assert 'data["PhysicalCount"] =!= expected' in source
    assert 'Length[candidates] =!= expected' in source
    assert 'Return[$Failed]' in source
    runner = (PROJECT / 'scripts/RunGeneralScalarTrial.ps1').read_text(encoding='utf-8')
    assert "if ($Scope -eq 'Problems')" in runner
    assert "$expectedQuartics = @{" in runner
    assert "'IncompleteOperatorBasis'" in runner
    assert '[int]$candidateCount -eq $expected' in runner


def test_matching_inspector_checks_extended_log_and_text_comparison(tmp_path):
    import sys
    sys.path.insert(0, str(PROJECT / 'scripts'))
    from InspectGeneralScalarMatching import inspect_model
    p = tmp_path / 'B_m3'
    for variant in ('baseline', 'extended'):
        (p / variant).mkdir(parents=True)
        (p / variant / 'c5_coefficient.txt').write_text('lambdaT3[] + x')
    (p / 'extended_matching.log').write_text(
        'Accepted interactions: {Yukawa1_GenericCG, T3ScalarMix_GenericCG}')
    result = inspect_model(p)
    assert result['expected_extra_quartics'] == 1
    assert result['accepted_extra_quartics_in_extended_log'] == 0
    assert result['status'] == 'IncompleteExtendedInteractions'
    assert result['baseline_extended_c5_text_equal_ignoring_whitespace'] is True
    (p / 'extended_matching.log').write_text(
        'Accepted interactions: {Yukawa1_GenericCG, T3ExtraQ_0_2_0_0_0_2_1}')
    result = inspect_model(p)
    assert result['accepted_extra_quartics_in_extended_log'] == 1
    assert result['extended_log_operator_count_matches_manifest'] is True
    assert result['status'] == 'NoNewCouplingsSeenInAvailableC5Artifacts'


def test_verified_catalog_covers_every_supported_scalar_point():
    import sys
    sys.path.insert(0, str(PROJECT / 'scripts'))
    from GenerateGeneralScalarWLManifest import DIMS
    manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))['models']
    set_of_scalar_points = {
        f'{DIMS[cls][0]}_{DIMS[cls][1]}_{int(alpha)}'
        for cls, alpha in (name.split('_') for name in manifest)
    }
    text = (PROJECT / 'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl').read_text(encoding='utf-8')
    assert len(set_of_scalar_points) == 13
    assert all(f'"{point}" ->' in text for point in set_of_scalar_points)
