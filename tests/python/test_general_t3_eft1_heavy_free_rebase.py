"""Offline source/regression checks for the experimental scalar free-Lagrangian rebase.

Real Matchete verification still requires wolframscript on the local project.
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import RunGeneralT3EFT1ScalarThresholdProbe as task
from tests.python.test_general_t3_eft1_scalar_threshold_probe import write_artifacts


ROOT = Path(__file__).resolve().parents[2]
WL = ROOT / 'Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl'


def source() -> str:
    return WL.read_text(encoding='utf-8')


def test_source_captures_free_lagrangians_on_both_sides_of_reclassification():
    src = source()
    before = src.index('lightEigenFree=CheckAbort[FreeLag[NewScalarX,NewScalarY]')
    change = src.index('T3MECheckpoint["scalar-threshold-reclassify-X-Y-started"]')
    after = src.index('heavyEigenFree=CheckAbort[FreeLag[NewScalarX,NewScalarY]')
    assert before < change < after
    assert 'Mass->{Light,MSX}' in src and 'Mass->{Heavy,MSX}' in src
    assert 'Mass->{Light,MSY}' in src and 'Mass->{Heavy,MSY}' in src


def test_rebase_limited_to_explicit_c5_mode():
    src = source()
    start = src.index('If[TrueQ[$T3MELeadingC5Only],', src.index('leadingC5Partition=<|'))
    end = src.index('T3MECheckpoint["scalar-threshold-EFT1-lagrangian-check-started"]')
    body = src[start:end]
    assert 'scalarEFT=leadingC5Partition["candidate"]' in body
    assert 'rebasedSource=Expand[scalarEFT-lightEigenFree+heavyEigenFree]' in body
    assert 'scalarEFT=rebasedSource' in body
    assert '"full_saved_eft_unchanged"->True' in body
    assert 'T3MEFail["ScalarThresholdFreeLagRebaseInvalid"]' in body
    assert 'T3MECheckpoint["scalar-threshold-light-to-heavy-free-rebase-returned"]' in body
    assert 'CheckLagrangian[scalarEFT]' in src[end:]


def test_match_audit_is_written_before_loop_guard():
    src = source()
    assert src.index('matchAuditPath=FileNameJoin') < src.index('treeResult=CheckAbort[RunT3Matching')
    assert src.index('Export[matchAuditPath') < src.index('If[correction === 0')
    assert src.index('If[correction === 0') < src.index('T3MECheckpoint["scalar-threshold-LLHH-extract-started"]')
    assert '"RawEFT","EOMEFT","MatchedEFT"' in src
    assert '"physical_sequential_C5_computed"->False' in src
    assert '"full_eft_validated"->!TrueQ[$T3MELeadingC5Only]' in src


def test_explicit_hbar_and_absent_heavy_fields_remain_mandatory():
    src = source()
    assert 'correction === 0 || !StringContainsQ[correctionText,"hbar"]' in src
    assert 'StringContainsQ[completeText,"Field[NewScalarX"]' in src
    assert 'StringContainsQ[completeText,"Field[NewScalarY"]' in src
    assert 'StringContainsQ[completeText,"Field[NewFermion"]' in src
    assert 'T3MEFail["ScalarThresholdLoopCorrectionInvalid"' in src


def test_failure_diagnostic_reports_stage_summary(monkeypatch, tmp_path):
    root = tmp_path / 'repo'
    (root / 'Lagrangian').mkdir(parents=True)
    (root / 'output').mkdir()
    wolfram = root / 'Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl'
    wolfram.write_text('Matchete source')
    matcher = root / 'Lagrangian/RunMatching.wl'
    matcher.write_text('RunT3Matching')
    params = root / 'output/parameters.json'
    params.write_text('{}')
    cert = root / 'output/cert.json'; cert.write_text('{}')
    tree = root / 'output/source_eft1_tree.inputform.txt'; tree.write_text('0123456789'*4)
    monkeypatch.setattr(task, 'ROOT', root)
    monkeypatch.setattr(task, 'WOLFRAM', wolfram)
    monkeypatch.setattr(task, 'MATCHING_SOURCE', matcher)
    monkeypatch.setattr(task, 'prepare_model', lambda *a, **k: {'payload': {'unitary_U': {}}, 'parameters': params})
    monkeypatch.setattr(task, '_prior_certificate', lambda *a, **k: cert)
    monkeypatch.setattr(task, '_verified_tree_input', lambda *a, **k: (tree,cert,cert))
    monkeypatch.setattr(task.shutil, 'which', lambda _: 'wolframscript')

    def fake_run(cmd, **kwargs):
        target = Path(cmd[-1]); stem = target.stem
        target.write_text(json.dumps({'status':'ScalarThresholdLoopCorrectionInvalid'}))
        target.with_name(stem+'.progress.json').write_text(json.dumps({'last_stage':'scalar-threshold-one-loop-Match-returned'}))
        target.with_name(stem+'_matching_stage_audit.json').write_text(json.dumps({
            'tree_stage': {'RawEFT': {'head':'Times','characters':1,'explicit_hbar':False,
                                       'contains_X':True,'contains_Y':True,'contains_GammaCC':True}},
            'one_loop_stage': {'MatchedEFT': {'head':'Times','characters':1,
                    'explicit_hbar':False,'contains_X':True,'contains_Y':True,'contains_GammaCC':True}},
        }))
        return SimpleNamespace(returncode=5, stdout='zero correction',stderr='')

    monkeypatch.setattr(task.subprocess, 'run', fake_run)
    with pytest.raises(ValueError) as error:
        task.run_model('B_-1', outdir=root/'output'/'stage',c5_relevant_only=True)
    msg = str(error.value)
    assert 'ScalarThresholdLoopCorrectionInvalid' in msg
    assert 'Matchete stage summaries' in msg
    assert '"explicit_hbar":false' in msg
    assert 'matching_stage_audit.json' in msg


def test_old_matching_failure_is_not_misreported_as_success(tmp_path):
    from tests.python.test_general_t3_eft1_scalar_threshold_probe import certificate, make_tree
    stem = task._stem('B_-1','final',c5_relevant_only=True)
    tree = make_tree(tmp_path)
    write_artifacts(tmp_path,stem)
    report=certificate(stem=stem)
    report['scalar_one_loop_correction_hbar_present']=False
    with pytest.raises(ValueError):
        task.validate_report(report,'B_-1',tmp_path,stem,tree,c5_relevant_only=True)
