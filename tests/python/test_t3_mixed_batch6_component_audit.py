"""Offline Batch 6 v6 residual-evidence tests; not live Wolfram tests."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import types

import pytest

from scripts import RunT3MixedScalarHardBatch6 as native
from scripts import RunT3MixedComponentAuditBatch6 as audit


def fixture_dir(tmp_path: Path, model='B_-1') -> Path:
    stem = audit.MODELS[model]
    folder = tmp_path / stem
    folder.mkdir()
    prefix = stem + '_full_eft1_residual'
    files = {key: prefix + '_' + basename + '.inputform.txt' for key, basename in (
        ('saved_eomeft', 'saved_eomeft_full'),
        ('reconstructed_candidate', 'reconstructed_candidate_full'),
        ('symbolic_residual', 'full_symbolic_residual'),
    )}
    for fn in files.values():
        (folder/fn).write_text('0', encoding='utf-8')
    groups = []
    for i in range(25):
        name = prefix+f'_residual_group_{i+1:02}.inputform.txt'
        (folder/name).write_text('0', encoding='utf-8')
        groups.append(dict(structural_sector=f'lambdaS1__XX{i}',
                           additive_terms=6 if i < 12 else 5, file=name))
    assert sum(z['additive_terms'] for z in groups) == 137
    manifest = {
        'format':'t3-bc-complete-eft1-reconstruction-residual-v1',
        'status':'CompleteSymbolicResidualExportedEquivalenceUnproven',
        'model':model, 'matching_scope':'FullEFT',
        'coupling_group_partition_conserved':True,
        'contains_majorana_residual':False,
        'physical_matching_certified':False,
        'CG_and_index_tensor_equivalence_proven':False,
        'residual_additive_terms':137, 'coupling_groups':groups,'files':files,
    }
    def json_write(name,obj):
        (folder/name).write_text(json.dumps(obj), encoding='utf-8')
    json_write(prefix+'_full_residual_manifest.json',manifest)
    json_write('hard_scalar_eft1_validation.json',{
        'model':model,'initial_check':{'checks':{
            'Hermiticity':'False','CanonicallyNormalized':'True',
            'HeavyMassBasis':'True','GaugeInvariant':'True',
            'AllObjectsDefined':'True'}}, 'full_EFT_validated':False,
    })
    json_write('hard_scalar_reconstruction_audit.json',{
        'algebraic_equality_proven':False,
        'independent_matchete_check_passed':False,
        'component_equivalence_proven':False,
        'residual_manifest_file':prefix+'_full_residual_manifest.json',
        'residual_additive_terms':137,'residual_structural_sectors':25,
    })
    json_write('native_result.progress.json',{'stage_history':[
        'script-initialized', 'batch6-full-component-residual-export-started',
        'batch6-full-component-residual-exported',
        'batch6-native-cross-CG-array-exported',
    ]})
    return folder


def test_v6_runner_exports_full_residual_before_match_and_stays_fail_closed(tmp_path):
    src = native.NATIVE_HELPER
    assert 'T3B6ExportFullResidual[rebased,reconstructed,out]' in src
    assert src.index('T3B6ExportFullResidual[rebased,reconstructed,out]') < \
           src.index('T3B6Validate[reconstructed]')
    assert 'batch6-full-component-residual-exported' in src
    assert 'T3B6ExportNativeCrossCG[out]' in src
    assert 'CrossInv1_coeff_A_inputform' in src
    assert 'CrossInv2_coeff_B_inputform' in src
    assert '"contains_majorana_residual"->False' in src
    assert '"physical_matching_certified"->False' in src
    assert '"coupling_group_partition_conserved"->True' in src
    assert '"component_equivalence_proven"->False' in src
    generated = native.instrument(
        'scriptDir = DirectoryName @ ExpandFileName[$InputFileName];\n'
        'T3MEPerform[] := Module[{},\n'
        ' T3MECheckpoint["fermion-first-lagrangian-checks-returned"];\n'
        ' T3MECheckpoint["fermion-first-production-Match-loader-started"];\n'
        ' T3MECheckpoint["fermion-first-tree-Match-started"];\n'
        ' matchingResult=CheckAbort[RunT3Matching[lag,5,0],$Aborted];\n'
        '];\nrun = CheckAbort[old[]];\n',tmp_path)
    assert generated.index('Return[T3B6Hard[lag]];') < \
           generated.index('RunT3Matching[lag,5,0]')
    assert generated.count('T3B6ExportFullResidual[') == 2  # definition + call


def test_directory_evidence_is_read_only_and_blocks_traversal(tmp_path):
    ev = audit.DirectoryEvidence(tmp_path)
    (tmp_path/'local.inputform.txt').write_text('a', encoding='utf-8')
    assert ev.text('local.inputform.txt') == 'a'
    assert len(ev.sha('local.inputform.txt')) == 64
    for bad in ('../another', 'x/y', '', 'name with spaces', 'x\\y'):
        with pytest.raises(audit.AuditError):
            ev.text(bad)


@pytest.mark.parametrize('model',list(audit.MODELS))
def test_audit_inputs_require_failed_hermiticity_and_full_partition(tmp_path,model):
    folder=fixture_dir(tmp_path,model)
    rec=audit.inspect_inputs(folder,model)
    assert rec['terms']==137 and rec['sectors']==25
    check_file=folder/'hard_scalar_eft1_validation.json'
    check=json.loads(check_file.read_text())
    check['initial_check']['checks']['GaugeInvariant']='False'
    check_file.write_text(json.dumps(check))
    with pytest.raises(audit.AuditError,match='another Matchete invariant'):
        audit.inspect_inputs(folder,model)


def test_audit_rejects_unexpected_sectors(tmp_path):
    folder=fixture_dir(tmp_path)
    stem=audit.MODELS['B_-1']
    path=folder/f'{stem}_full_eft1_residual_full_residual_manifest.json'
    man=json.loads(path.read_text())
    man['coupling_groups'].pop()
    path.write_text(json.dumps(man))
    with pytest.raises(audit.AuditError,match='137-term/25-sector'):
        audit.inspect_inputs(folder,'B_-1')


def test_audit_rejects_post_match_progress(tmp_path):
    folder=fixture_dir(tmp_path)
    path=folder/'native_result.progress.json'
    progress=json.loads(path.read_text())
    progress['stage_history'].append('batch6-native-scalar-one-loop-match-started')
    path.write_text(json.dumps(progress))
    with pytest.raises(audit.AuditError,match='before scalar Match'):
        audit.inspect_inputs(folder,'B_-1')


def test_conditions_only_pass_with_old_exact_native_tensor_solver(monkeypatch,tmp_path):
    folder=fixture_dir(tmp_path)
    fake_comp=types.ModuleType('scripts.T3BCSU2ComponentEquivalence')
    fake_comp.SYM={f'Other{i}':() for i in range(7)}
    def fake_solver(ev,model):
        assert ev.text(f'{model}_full_eft1_residual_full_residual_manifest.json')
        other={k:{'matrix_rank':0,'augmented_rank':0} for k in fake_comp.SYM}
        return ([{'raw_terms':6 if i<12 else 5} for i in range(25)],
            other | {
                'T3MENativeCrossCGInv1':{
                    'status':'UNIQUE','required_coefficients':['1','-1']},
                'T3MENativeCrossCGInv2':{
                    'status':'NONUNIQUE','matrix_rank':1,
                    'required_coefficients':['tau0','tau0']},
            })
    fake_comp.audit_model=fake_solver
    monkeypatch.setitem(sys.modules,'scripts.T3BCSU2ComponentEquivalence',fake_comp)
    fake_verify=types.ModuleType('scripts.ResolveT3BCEFT1HermiticityAndYGate')
    fake_verify.check_actual_native_arrays=lambda native_dir,stem: {
        'native_entries_independently_verified':32}
    monkeypatch.setitem(sys.modules,'scripts.ResolveT3BCEFT1HermiticityAndYGate',fake_verify)
    row=audit.audit_model_run(folder,'B_-1')
    assert row['component_equivalence_under_native_CG'] is True
    assert row['scalar_Match_performed'] is False
    assert row['physical_C5_certified'] is False
    def no_native(*_):
        raise ValueError('Missing exact native arrays')
    fake_verify.check_actual_native_arrays=no_native
    row=audit.audit_model_run(folder,'B_-1')
    assert row['component_equivalence_under_native_CG'] is False
    assert row['status']=='ConditionalComponentIdentityNativeArraysUnavailable'


def test_legacy_absolute_import_works_under_python_m_without_github_edits(tmp_path):
    """Reproduce the real helper's `from T3BC...` in a clean process.

    The old mocked solver test did not exercise ResolveT3BCEFT1HermiticityAndYGate's
    top-level import, which is why the v6 import failure escaped offline checks.
    """
    scripts = tmp_path / 'scripts'
    scripts.mkdir()
    shutil.copy2(Path(audit.__file__), scripts / 'RunT3MixedComponentAuditBatch6.py')
    (scripts / 'T3BCSU2ComponentEquivalence.py').write_text(
        'MARKER = "top_level_legacy_import_resolved"\n', encoding='utf-8')
    (scripts / 'ResolveT3BCEFT1HermiticityAndYGate.py').write_text(
        'from T3BCSU2ComponentEquivalence import MARKER\n', encoding='utf-8')
    command = [sys.executable, '-c',
        'from scripts import RunT3MixedComponentAuditBatch6 as audit; '
        'audit._enable_legacy_absolute_imports(); '
        'from scripts.ResolveT3BCEFT1HermiticityAndYGate import MARKER; '
        'print(MARKER)']
    env = dict(os.environ)
    env['PYTHONPATH'] = str(tmp_path)
    proc = subprocess.run(command, cwd=tmp_path, env=env,
                          capture_output=True, text=True, timeout=20)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == 'top_level_legacy_import_resolved'


def test_component_audit_resolves_legacy_import_before_loading_verifier():
    import inspect
    text = inspect.getsource(audit.audit_model_run)
    assert text.index('_enable_legacy_absolute_imports()') < \
        text.index('from scripts.ResolveT3BCEFT1HermiticityAndYGate import')
