"""Offline fail-closed tests for the opt-in one-loop LLHH subtheory.

These do not execute Wolfram or certify its one-loop power counting.
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import RunGeneralT3EFT1ScalarThresholdProbe as task
from tests.python.test_general_t3_eft1_scalar_threshold_probe import (
    certificate, make_tree, write_artifacts,
)


def restricted_report(stem: str):
    report = certificate(stem=stem)
    report.update({
        'matching_scope': 'LeadingC5OneLoopSectorOnly',
        'full_eft_validated': False,
        'full_saved_eft_unchanged_on_disk': True,
        'hermiticity_applies_to_matching_subtheory_only': True,
        'excluded_four_heavy_scalar_quartic_terms': 122,
        'full_saved_eft_hermiticity': False,
    })
    return report


def restricted_fixtures(tmp_path: Path):
    stem = task._stem('B_-1', 'final', c5_relevant_only=True)
    tree = make_tree(tmp_path)
    write_artifacts(tmp_path, stem)
    diagnostic = tmp_path / (stem + '_hermiticity_diagnostics.json')
    info = json.loads(diagnostic.read_text())
    info.update({
        'matching_scope': 'LeadingC5OneLoopSectorOnly',
        'full_saved_eft_hermiticity': False,
        'excluded_four_heavy_scalar_quartic_terms': 122,
    })
    diagnostic.write_text(json.dumps(info))
    return stem, tree


def test_distinct_stem_prevents_overwriting_full_eft():
    normal = task._stem('B_-1', 'final')
    restricted = task._stem('B_-1', 'final', c5_relevant_only=True)
    assert restricted == normal + '_c5_relevant'


def test_restricted_scope_validated(tmp_path):
    stem, tree = restricted_fixtures(tmp_path)
    files = task.validate_report(restricted_report(stem), 'B_-1', tmp_path,
                                 stem, tree, c5_relevant_only=True)
    assert len(files) == 4


def test_restricted_result_rejected_as_full_eft(tmp_path):
    stem, tree = restricted_fixtures(tmp_path)
    with pytest.raises(ValueError, match='restricted C5'):
        task.validate_report(restricted_report(stem), 'B_-1', tmp_path, stem, tree)


@pytest.mark.parametrize('key,bad', [
    ('matching_scope', 'FullEFT'),
    ('excluded_four_heavy_scalar_quartic_terms', 121),
    ('full_eft_validated', True),
    ('full_saved_eft_unchanged_on_disk', False),
    ('hermiticity_applies_to_matching_subtheory_only', False),
    ('full_saved_eft_hermiticity', True),
])
def test_restricted_certificate_cannot_overclaim(tmp_path, key, bad):
    stem, tree = restricted_fixtures(tmp_path)
    report = restricted_report(stem)
    report[key] = bad
    with pytest.raises(ValueError):
        task.validate_report(report, 'B_-1', tmp_path, stem, tree,
                             c5_relevant_only=True)


def test_restricted_diagnostic_cannot_overclaim(tmp_path):
    stem, tree = restricted_fixtures(tmp_path)
    path = tmp_path / (stem + '_hermiticity_diagnostics.json')
    info = json.loads(path.read_text())
    info['excluded_four_heavy_scalar_quartic_terms'] = 123
    path.write_text(json.dumps(info))
    with pytest.raises(ValueError, match='provenance'):
        task.validate_report(restricted_report(stem), 'B_-1', tmp_path,
                             stem, tree, c5_relevant_only=True)


def test_source_only_projects_in_explicit_opt_in():
    wl = Path(__file__).resolve().parents[2] / 'Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl'
    src = wl.read_text()
    assert '$T3MELeadingC5Only = (Environment["T3_EFT1_LEADING_C5_ONLY"] === "1")' in src
    assert 'If[TrueQ[$T3MELeadingC5Only],' in src
    assert 'T3MELeadingC5OneLoopProjection[scalarEFT]' in src
    assert 'Length[spectator]=!=122' in src
    assert 'Length[higgsQuartics]=!=52' in src
    assert 'Length[background]=!=21' in src
    assert 'retCount=!=21+52+expCount' in src
    assert 'StringCount[text,"Field[H, Scalar"]===0' in src
    assert 'FreeQ[expr,GammaCC]' in src
    assert 'full_saved_eft_unchanged_on_disk' in src
    assert 'full_eft_validated' in src
    assert 'RunT3Matching[scalarEFT,5,0]' in src
    assert 'RunT3Matching[scalarEFT,5,1]' in src


def test_python_opt_in_passed_only_to_child(monkeypatch, tmp_path):
    root = tmp_path / 'repo'
    (root / 'Lagrangian').mkdir(parents=True)
    (root / 'output').mkdir()
    wl = root / 'Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl'
    wl.write_text('Matchete source')
    match = root / 'Lagrangian/RunMatching.wl'
    match.write_text('RunT3Matching')
    params = root / 'output/parameters.json'; params.write_text('{}')
    mass = root / 'output/mass.json'; mass.write_text('{}')
    tcert = root / 'output/tree.json'; tcert.write_text('{}')
    lcert = root / 'output/oneloop.json'; lcert.write_text('{}')
    tree = root / 'output/source_eft1_tree.inputform.txt'
    tree.write_text('0123456789' * 4)
    for name, val in [('ROOT', root), ('WOLFRAM', wl), ('MATCHING_SOURCE', match)]:
        monkeypatch.setattr(task, name, val)
    monkeypatch.setattr(task, 'prepare_model',
                        lambda *a, **kw: {'payload': {'unitary_U': {}}, 'parameters': params})
    monkeypatch.setattr(task, '_prior_certificate', lambda *a, **kw: mass)
    monkeypatch.setattr(task, '_verified_tree_input',
                        lambda *a, **kw: (tree, tcert, lcert))
    monkeypatch.setattr(task.shutil, 'which', lambda x: 'wolframscript')
    def fake_run(cmd, **kwargs):
        assert kwargs['env']['T3_EFT1_LEADING_C5_ONLY'] == '1'
        target = Path(cmd[-1]); stem = target.stem
        assert stem.endswith('_c5_relevant')
        write_artifacts(target.parent, stem)
        diag = target.parent / (stem + '_hermiticity_diagnostics.json')
        obj = json.loads(diag.read_text())
        obj.update({
            'matching_scope': 'LeadingC5OneLoopSectorOnly',
            'full_saved_eft_hermiticity': False,
            'excluded_four_heavy_scalar_quartic_terms': 122,
        })
        diag.write_text(json.dumps(obj))
        target.write_text(json.dumps(restricted_report(stem)))
        return SimpleNamespace(returncode=0, stdout='limited success', stderr='')
    monkeypatch.setattr(task.subprocess, 'run', fake_run)
    row = task.run_model('B_-1', outdir=root / 'output/restricted',
                         c5_relevant_only=True)
    assert row['matching_scope'] == 'LeadingC5OneLoopSectorOnly'
    assert row['physical_sequential_C5_computed'] is False
