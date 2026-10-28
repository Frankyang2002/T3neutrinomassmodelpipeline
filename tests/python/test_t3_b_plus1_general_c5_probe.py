"""Offline interface tests for Batch 19. Not a substitute for live Matchete."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import ProbeT3BPlus1GeneralC5 as probe


def valid_registry():
    models = {f'X_{i}': {} for i in range(15)}
    models['B_+1'] = {
        'new_rgbeta_symbols': [
            {'name': 'lambdaT3GeneralQ022000I1', 'counts': [0, 2, 2, 0, 0, 0]}],
        'preexisting_rgbeta_symbols': [], 'quadratics_not_implemented': 0}
    return {'models': models}


def fixture_tree(tmp_path):
    for path in probe.SOURCES:
        p = tmp_path / path
        p.parent.mkdir(parents=True, exist_ok=True)
        if path == probe.REGISTRY:
            p.write_text(json.dumps(valid_registry()), encoding='utf-8')
        else:
            p.write_text('source', encoding='utf-8')
    return tmp_path


def write_artifacts(folder, *, extended=False, changes=True, model_class='T3-d2-d2-F1'):
    folder.mkdir(exist_ok=True, parents=True)
    names = ['Yukawa1_GenericCG', 'Yukawa2_GenericCG', 'T3ScalarMix_GenericCG']
    if extended:
        names.append(probe.INTERACTION)
    data = {
        'ModelClass': model_class, 'Alpha': 1,
        'Scalar1SU2': 2, 'Scalar2SU2': 2, 'FermionSU2': 1,
        'Scalar1Hypercharge': '1/2', 'Scalar2Hypercharge': '3/2',
        'FermionHypercharge': '1', 'BuildStatus': 'Success',
        'MatchingStatus': 'Success', 'WeinbergExtractionStatus': 'Success',
        'T3IngredientsPresent': True, 'WeinbergOperatorPresent': True,
        'WeinbergHolomorphicTermCount': 2, 'AcceptedInteractions': names}
    (folder / 'comparison_summary.json').write_text(json.dumps(data), encoding='utf-8')
    text = ('Coupling[y1, {1}, 0]*Coupling[y2, {2}, 0] *'
            'Coupling[lambdaT3,{},0]/Coupling[MF,{},0]')
    if extended and changes:
        text += ' + Coupling[lambdaT3ExtraQ022000I1,{},0]*Coupling[y1,{},0]*Coupling[y2,{},0]'
    (folder / 'c5_coefficient.txt').write_text(text, encoding='utf-8')


def valid_numeric():
    return {
        'format': 't3-b-plus1-c5-wolfram-numeric-audit-v1',
        'status': 'Success', 'sample_count': 3,
        'baseline_agrees_with_general_at_zero_extra': True,
        'nonzero_extra_changes_c5_at_any_sample': False,
        'samples': [{'zero_extra_agrees': True, 'all_numeric_values_finite': True} for _ in range(3)],
    }


def test_registry_valid(tmp_path):
    fixture_tree(tmp_path)
    assert probe.check_registry(tmp_path)['model'] == 'B_+1'


def test_registry_missing_quartic_fails(tmp_path):
    fixture_tree(tmp_path)
    p = tmp_path / probe.REGISTRY
    d = json.loads(p.read_text()); d['models']['B_+1']['new_rgbeta_symbols'] = []
    p.write_text(json.dumps(d))
    with pytest.raises(RuntimeError, match='quartic'):
        probe.check_registry(tmp_path)


def test_registry_mixing_fails(tmp_path):
    fixture_tree(tmp_path)
    p = tmp_path / probe.REGISTRY
    d = json.loads(p.read_text()); d['models']['B_+1']['quadratics_not_implemented'] = 1
    p.write_text(json.dumps(d))
    with pytest.raises(RuntimeError, match='mixing'):
        probe.check_registry(tmp_path)


def test_fingerprints(tmp_path):
    fixture_tree(tmp_path)
    assert len(probe.fingerprints(tmp_path)) == len(probe.SOURCES)
    (tmp_path / probe.WOLFRAM_MODEL).unlink()
    with pytest.raises(RuntimeError, match='missing source'):
        probe.fingerprints(tmp_path)


@pytest.mark.parametrize('bad', ('../bad', 'a/b', '', 'ab c', '.bad'))
def test_bad_run_ids(bad, tmp_path):
    with pytest.raises(RuntimeError, match='run ID'):
        probe.expected_destination(bad, tmp_path)


def test_good_run_id(tmp_path):
    assert probe.expected_destination('batch19_trial_1', tmp_path).name == 'batch19_trial_1'


def test_artifacts(tmp_path):
    base = tmp_path / 'baseline'; general = tmp_path / 'general'
    write_artifacts(base); write_artifacts(general, extended=True)
    assert not probe.checked_artifacts(base, extended=False)['extra_quartic_symbol_occurs_in_ordered_c5']
    assert probe.checked_artifacts(general, extended=True)['extra_quartic_symbol_occurs_in_ordered_c5']


@pytest.mark.parametrize('label', ('T3-d2-d2-F1', 'T3-B'))
def test_b_plus1_model_class_aliases(tmp_path, label):
    folder = tmp_path / 'baseline'
    write_artifacts(folder, model_class=label)
    report = probe.checked_artifacts(folder, extended=False)
    assert report['reported_model_class'] == label


@pytest.mark.parametrize('label', ('T3-d2-d2-F1', 'T3-B'))
def test_b_plus1_class_aliases_still_reject_wrong_representation(tmp_path, label):
    folder = tmp_path / 'baseline'
    write_artifacts(folder, model_class=label)
    path = folder / 'comparison_summary.json'
    doc = json.loads(path.read_text(encoding='utf-8'))
    doc['FermionSU2'] = 3
    path.write_text(json.dumps(doc), encoding='utf-8')
    with pytest.raises(RuntimeError, match=r'not B\(\+1\)'):
        probe.checked_artifacts(folder, extended=False)


def test_extended_interaction_required(tmp_path):
    folder = tmp_path / 'general';write_artifacts(folder)
    with pytest.raises(RuntimeError, match='extra quartic'):
        probe.checked_artifacts(folder, extended=True)


def test_restricted_no_extra(tmp_path):
    folder = tmp_path / 'baseline';write_artifacts(folder, extended=True)
    with pytest.raises(RuntimeError):
        probe.checked_artifacts(folder, extended=False)


def test_incomplete_matching_rejected(tmp_path):
    folder = tmp_path / 'baseline';write_artifacts(folder)
    p=folder/'comparison_summary.json';d=json.loads(p.read_text());d['MatchingStatus']='MatchFailed';p.write_text(json.dumps(d))
    with pytest.raises(RuntimeError, match='not successful'):
        probe.checked_artifacts(folder, extended=False)


def test_numeric_report(tmp_path):
    p=tmp_path/'a.json';p.write_text(json.dumps(valid_numeric()))
    assert not probe.checked_numeric_report(p)['nonzero_extra_changes_c5_at_any_sample']


def test_numeric_mismatch_rejected(tmp_path):
    d=valid_numeric();d['samples'][2]['zero_extra_agrees']=False
    p=tmp_path/'a.json';p.write_text(json.dumps(d))
    with pytest.raises(RuntimeError, match='failed'):
        probe.checked_numeric_report(p)


def test_numeric_false_certification_rejected(tmp_path):
    d=valid_numeric();d['status']='Failed'
    p=tmp_path/'a.json';p.write_text(json.dumps(d))
    with pytest.raises(RuntimeError, match='did not succeed'):
        probe.checked_numeric_report(p)


def test_execute_mocked_success(tmp_path):
    fixture_tree(tmp_path)
    calls=[]
    def fake_run(cmd,env,log,timeout,*,cwd):
        calls.append((cmd,env.get('T3_EXPERIMENTAL_GENERAL_QUARTICS')))
        if len(calls)<=2:
            assert cmd[1] == '-file' and cmd[-5:-1] != ['bad']
            assert env['T3_EXPERIMENTAL_GENERAL_QUARTICS'] == ('0' if len(calls)==1 else '1')
            write_artifacts(Path(cmd[3]), extended=len(calls)==2)
        else:
            Path(cmd[-1]).write_text(json.dumps(valid_numeric()))
    report=probe.execute('trial1',root=tmp_path,runner=fake_run)
    assert report['status']=='OneScaleSourceAndNumericAuditPassed'
    assert len(calls)==3
    p=probe.expected_destination('trial1',tmp_path)
    assert (p/'summary.md').is_file() and (p/'report.json').is_file()
    assert report['three_study_neutrino_refits_completed'] is False
    assert report['physical_matchete_rgbeta_map_certified'] is False
    with pytest.raises(RuntimeError, match='already exists'):
        probe.execute('trial1',root=tmp_path,runner=fake_run)


def test_execute_accepts_local_short_class_label(tmp_path):
    fixture_tree(tmp_path)
    calls = []
    def fake_run(cmd, env, log, timeout, *, cwd):
        calls.append(cmd)
        if len(calls) < 3:
            write_artifacts(Path(cmd[3]), extended=len(calls) == 2,
                            model_class='T3-B')
        else:
            Path(cmd[-1]).write_text(json.dumps(valid_numeric()))
    report = probe.execute('localclass1', root=tmp_path, runner=fake_run)
    assert len(calls) == 3
    assert report['baseline']['reported_model_class'] == 'T3-B'
    assert report['general']['reported_model_class'] == 'T3-B'


def test_execute_preserves_failed_diagnostics(tmp_path):
    fixture_tree(tmp_path)
    def failed_run(cmd,env,log,timeout,*,cwd):
        raise RuntimeError('native matching failed')
    with pytest.raises(RuntimeError, match='native matching failed'):
        probe.execute('failure1',root=tmp_path,runner=failed_run)
    assert probe.expected_destination('failure1',tmp_path).is_dir()


@pytest.mark.parametrize('key,bad', [
    ('ModelClass', 'B'), ('Alpha', -1), ('Scalar1SU2', 3),
    ('FermionSU2', 3), ('Scalar2Hypercharge', '1/2'),
])
def test_wrong_model_rejected(tmp_path, key, bad):
    folder = tmp_path / 'baseline'
    write_artifacts(folder)
    path = folder / 'comparison_summary.json'
    doc = json.loads(path.read_text(encoding='utf-8'))
    doc[key] = bad
    path.write_text(json.dumps(doc), encoding='utf-8')
    with pytest.raises(RuntimeError, match=r'not B\(\+1\)'):
        probe.checked_artifacts(folder, extended=False)


def test_resume_existing_legacy_baseline(tmp_path):
    fixture_tree(tmp_path)
    dest = probe.expected_destination('legacy1', tmp_path)
    dest.mkdir(parents=True)
    write_artifacts(dest / 'baseline')
    (dest / 'baseline.log').write_text('Matchete complete', encoding='utf-8')
    calls = []
    def fake_run(cmd, env, log, timeout, *, cwd):
        calls.append((cmd, env.get('T3_EXPERIMENTAL_GENERAL_QUARTICS')))
        if len(calls) == 1:
            assert env['T3_EXPERIMENTAL_GENERAL_QUARTICS'] == '1'
            write_artifacts(Path(cmd[3]), extended=True)
        else:
            Path(cmd[-1]).write_text(json.dumps(valid_numeric()))
    result = probe.execute('legacy1', root=tmp_path, runner=fake_run, resume=True)
    assert len(calls) == 2
    assert result['status'] == 'OneScaleSourceAndNumericAuditPassed'
    assert result['baseline_source_hash_certified'] is False
    assert result['resumed_interrupted_run'] is True
    assert result['baseline_source_provenance'].startswith('LegacyBaselineReused')


def test_resume_local_t3_b_label_with_saved_fingerprints(tmp_path):
    fixture_tree(tmp_path)
    dest = probe.expected_destination('shortlabel1', tmp_path)
    dest.mkdir(parents=True)
    write_artifacts(dest / 'baseline', model_class='T3-B')
    (dest / 'baseline.log').write_text('Matchete complete', encoding='utf-8')
    (dest / 'source_fingerprints.json').write_text(
        json.dumps(probe.fingerprints(tmp_path)), encoding='utf-8')
    calls = []
    def fake_run(cmd, env, log, timeout, *, cwd):
        calls.append(cmd)
        if len(calls) == 1:
            assert env['T3_EXPERIMENTAL_GENERAL_QUARTICS'] == '1'
            write_artifacts(Path(cmd[3]), extended=True, model_class='T3-B')
        else:
            Path(cmd[-1]).write_text(json.dumps(valid_numeric()))
    result = probe.execute('shortlabel1', root=tmp_path, runner=fake_run, resume=True)
    assert len(calls) == 2
    assert result['baseline_source_hash_certified'] is True
    assert result['baseline']['reported_model_class'] == 'T3-B'
    assert result['general']['reported_model_class'] == 'T3-B'


def test_resume_rejects_changed_source(tmp_path):
    import os
    fixture_tree(tmp_path)
    dest = probe.expected_destination('legacy2', tmp_path)
    dest.mkdir(parents=True)
    write_artifacts(dest / 'baseline')
    (dest / 'baseline.log').write_text('Matchete complete', encoding='utf-8')
    source = tmp_path / probe.WOLFRAM_MODEL
    boundary = (dest / 'baseline' / 'comparison_summary.json').stat().st_mtime_ns
    os.utime(source, ns=(boundary + 10000000000, boundary + 10000000000))
    with pytest.raises(RuntimeError, match='source modified after original baseline'):
        probe.execute('legacy2', root=tmp_path, resume=True)


def test_resume_rejects_completed_run(tmp_path):
    fixture_tree(tmp_path)
    dest = probe.expected_destination('finished1', tmp_path)
    dest.mkdir(parents=True)
    write_artifacts(dest / 'baseline')
    (dest / 'baseline.log').write_text('Matchete complete', encoding='utf-8')
    (dest / 'summary.md').write_text('complete', encoding='utf-8')
    with pytest.raises(RuntimeError, match='already finalised'):
        probe.execute('finished1', root=tmp_path, resume=True)


def setup_complete_native_match(tmp_path, run_id):
    fixture_tree(tmp_path)
    dest = probe.expected_destination(run_id, tmp_path)
    dest.mkdir(parents=True)
    write_artifacts(dest / 'baseline', model_class='T3-B')
    write_artifacts(dest / 'general', extended=True, model_class='T3-B')
    (dest / 'baseline.log').write_text('baseline success', encoding='utf-8')
    (dest / 'general.log').write_text('general success', encoding='utf-8')
    (dest / 'source_fingerprints.json').write_text(
        json.dumps(probe.fingerprints(tmp_path)), encoding='utf-8')
    return dest


def test_resume_malformed_json_replays_only_numeric_audit(tmp_path):
    dest = setup_complete_native_match(tmp_path, 'jsonfail1')
    bad = '{"samples": [1.23*^-15]}'
    (dest / 'native_numeric_audit.json').write_text(bad, encoding='utf-8')
    (dest / 'numeric_audit.log').write_text('first audit log', encoding='utf-8')
    commands = []

    def fake_run(command, env, log, timeout, *, cwd):
        commands.append(command)
        assert Path(command[2]) == tmp_path / probe.WOLFRAM_COMPARE
        Path(command[-1]).write_text(json.dumps(valid_numeric()), encoding='utf-8')

    result = probe.execute('jsonfail1', root=tmp_path, runner=fake_run, resume=True)
    assert len(commands) == 1
    assert result['numeric_audit_recomputed'] is True
    assert result['malformed_prior_numeric_audit_detected'] is True
    assert result['baseline_source_hash_certified'] is True
    assert (dest / 'native_numeric_audit.before_json_fix.json').read_text() == bad
    assert (dest / 'numeric_audit.before_json_fix.log').read_text() == 'first audit log'


def test_resume_changed_audit_script_without_redoing_matchete(tmp_path):
    dest = setup_complete_native_match(tmp_path, 'jsonfix1')
    original = probe.fingerprints(tmp_path)
    (tmp_path / probe.WOLFRAM_COMPARE).write_text('new JSON exporter', encoding='utf-8')
    (dest / 'native_numeric_audit.json').write_text('invalid', encoding='utf-8')
    commands = []

    def fake_run(command, env, log, timeout, *, cwd):
        commands.append(command)
        assert Path(command[2]) == tmp_path / probe.WOLFRAM_COMPARE
        Path(command[-1]).write_text(json.dumps(valid_numeric()), encoding='utf-8')

    result = probe.execute('jsonfix1', root=tmp_path, runner=fake_run, resume=True)
    assert len(commands) == 1
    assert result['status'] == 'OneScaleSourceAndNumericAuditPassed'
    assert result['comparison_script_updated_since_match'] is True
    assert result['comparison_script_sha256_at_matching'] == original[probe.WOLFRAM_COMPARE]
    assert result['baseline_source_hash_certified'] is True
    assert result['full_original_source_manifest_unchanged'] is False
    assert (dest / 'native_numeric_audit.before_json_fix.json').read_text() == 'invalid'


def test_matching_source_change_blocks_resume(tmp_path):
    dest = setup_complete_native_match(tmp_path, 'matchsource1')
    (tmp_path / probe.WOLFRAM_MODEL).write_text('matching has changed', encoding='utf-8')
    (dest / 'native_numeric_audit.json').write_text('invalid', encoding='utf-8')
    with pytest.raises(RuntimeError, match='matching source SHA-256 inputs changed'):
        probe.execute('matchsource1', root=tmp_path, resume=True)


def test_native_writer_explicit_machine_numbers_and_json_roundtrip():
    source = (ROOT / probe.WOLFRAM_COMPARE).read_text(encoding='utf-8')
    assert 'toJSONReal[value_]' in source
    assert 'N[value, MachinePrecision]' in source
    assert '"zero_extra_relative_error" -> errorOut' in source
    assert 'ExportString[report, "RawJSON"]' in source
    assert 'ImportString[jsonText, "RawJSON"]' in source
