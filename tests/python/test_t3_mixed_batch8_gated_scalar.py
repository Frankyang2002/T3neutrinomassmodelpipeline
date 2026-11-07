"""Offline/fail-closed Batch 8 proof-replay tests. Native Matchete not available here."""
from __future__ import annotations

import copy
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import RunT3MixedScalarGatedBatch8 as v8
from scripts import RunT3MixedScalarHardBatch6 as v6


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(payload),encoding='utf-8')


def fixture(tmp_path: Path, model: str = 'B_-1'):
    source_root = tmp_path/'mixed_hard_06'
    target=source_root/v8.audit.MODELS[model]
    target.mkdir(parents=True)
    prefix=v8.audit.MODELS[model]+'_full_eft1_residual'
    paths={
        'saved_eomeft': target/(prefix+'_saved_eomeft_full.inputform.txt'),
        'reconstructed_candidate':target/(prefix+'_reconstructed_candidate_full.inputform.txt'),
        'symbolic_residual':target/(prefix+'_full_symbolic_residual.inputform.txt'),
    }
    for key,path in paths.items():
        path.write_text(key*5,encoding='utf-8')
    arr=target/(v8.audit.MODELS[model]+'_nativeCG_native_cross_cg_arrays.inputform.txt')
    arr.write_text('{{1,2},{3,4}}',encoding='utf-8')
    files_sha={key:sha256(path.read_bytes()).hexdigest() for key,path in paths.items()}
    row={
        'model':model,'status':v8.PROOF_STATUS,
        'component_equivalence_under_native_CG':True,
        'native_arrays_independently_checked':True,
        'arrays_from_current_native_run':True,
        'scalar_Match_performed':False,
        'matchete_Hermiticity_check_passed':False,
        'physical_C5_certified':False,
        'inputs':{'terms':137,'sectors':25,'files_sha256':files_sha},
        'native_array_evidence':{
            'native_entries_independently_verified':32,
            'native_tensor_sha256':sha256(arr.read_bytes()).hexdigest(),
        },
        'component_conditions':{'exact_component_polynomial_conditions_solved':True},
    }
    report={
        'format':v8.PROOF_FORMAT,'status':v8.PROOF_STATUS,
        'source_hard_run':source_root.name,'scalar_Match_performed':False,
        'fermion_Match_recomputed':False,
        'matchete_Hermiticity_check_passed':False,
        'physical_C5_certified':False,'models':[row],
    }
    current={'source_fingerprints':{'driver':'abcdef','canonical_eom':'fedcba'},
             'canonical_eom':str(tmp_path/'some-eom.inputform.txt')}
    settings={'source_fingerprints':{model:copy.deepcopy(current['source_fingerprints'])}}
    write_json(source_root/'component_equivalence_audit.json',report)
    write_json(source_root/'settings.json',settings)
    return SimpleNamespace(root=source_root,model=model,row=row,report=report,
                           current=current,settings=settings,files=paths,arrays=arr)


def test_certificate_rechecks_component_solver_and_source_hashes(tmp_path,monkeypatch):
    f=fixture(tmp_path)
    monkeypatch.setattr(v8.audit,'audit_model_run',lambda folder,model:copy.deepcopy(f.row))
    result=v8.verify_certificate(f.model,f.root,f.report,f.settings,f.current)
    assert result['component_equivalence_under_native_CG'] is True
    assert result['matchete_Hermiticity_check_passed'] is False
    assert result['physical_C5_certified'] is False
    assert result['source_component_run']=='mixed_hard_06'
    assert Path(result['reconstructed_witness'])==f.files['reconstructed_candidate']


@pytest.mark.parametrize('change,pattern',[ 
    ('source','witness modified'),
    ('arrays','array SHA256'),
    ('wrong_run','incompatible'),
    ('wrong_model','exactly one'),
    ('failed_conditions','not established'),
    ('fingerprint','source inputs changed'),
    ('physical_flag','incompatible'),
])
def test_certificate_fail_closed(tmp_path,monkeypatch,change,pattern):
    f=fixture(tmp_path)
    monkeypatch.setattr(v8.audit,'audit_model_run',lambda folder,model:copy.deepcopy(f.row))
    if change=='source':
        f.files['reconstructed_candidate'].write_text('INVALID',encoding='utf-8')
    elif change=='arrays':
        f.arrays.write_text('INVALID',encoding='utf-8')
    elif change=='wrong_run':
        f.report['source_hard_run']='another_run'
    elif change=='wrong_model':
        f.report['models'][0]['model']='C_-1'
    elif change=='failed_conditions':
        f.report['models'][0]['component_equivalence_under_native_CG']=False
    elif change=='fingerprint':
        f.current['source_fingerprints']['driver']='changed'
    elif change=='physical_flag':
        f.report['physical_C5_certified']=True
    with pytest.raises(v8.GateError,match=pattern):
        v8.verify_certificate(f.model,f.root,f.report,f.settings,f.current)


def test_certificate_cannot_be_self_attested(tmp_path,monkeypatch):
    f=fixture(tmp_path)
    candidate=copy.deepcopy(f.row)
    candidate['native_array_evidence']['native_tensor_sha256']='0'*64
    monkeypatch.setattr(v8.audit,'audit_model_run',lambda folder,model:candidate)
    with pytest.raises(v8.GateError,match='independent current component recheck'):
        v8.verify_certificate(f.model,f.root,f.report,f.settings,f.current)


def test_new_wolfram_script_checks_both_replay_witnesses_and_full_native_gate(tmp_path):
    helper=v8.gated_helper()
    assert helper.count('Match[tree,EFTOrder->5,LoopOrder->1]')==1
    assert helper.count('T3B6ExportFullResidual[rebased,reconstructed,out]')==1
    assert helper.count('T3B6ExportNativeCrossCG[out]')==1
    for expected in (
        'T3_BATCH8_PROOF_SAVED','T3_BATCH8_PROOF_RECONSTRUCTED',
        'T3_BATCH8_PROOF_ARRAYS','sameSaved','sameNative','sameTensor',
        'component_equivalence_proven"->True',
        'batch8-reconstruction-independent-CheckLagrangian-started',
        'CertifiedComponentReconstructionAndMatcheteChecksPassed',
    ):
        assert expected in helper
    gate=helper.index('If[chosen=!="VerifiedComponentEquivalentNativeReconstruction"')
    assert gate<helper.index('Match[tree,EFTOrder->5,LoopOrder->1]')
    assert 'algebraic_equality_proven"->TrueQ[equality]' in helper
    assert '"full_pole_cancellation_certified"->False' in helper
    assert '"physical_C5_certified"->False' in helper
    old=v6.NATIVE_HELPER
    try:
        v6.NATIVE_HELPER=helper
        generated=v6.instrument(
            'scriptDir = DirectoryName @ ExpandFileName[$InputFileName];\n'
            'T3MEPerform[] := Module[{},\n'
            ' T3MECheckpoint["fermion-first-lagrangian-checks-returned"];\n'
            ' T3MECheckpoint["fermion-first-production-Match-loader-started"];\n'
            ' T3MECheckpoint["fermion-first-tree-Match-started"];\n'
            ' matchingResult=CheckAbort[RunT3Matching[lag,5,0],$Aborted];\n'
            '];\nrun = CheckAbort[old[]];\n',tmp_path)
    finally:
        v6.NATIVE_HELPER=old
    assert generated.index('Return[T3B6Hard[lag]];') < generated.index('RunT3Matching[lag,5,0]')
    assert generated.index('batch8-component-certificate-replay-started') < generated.index('Match[tree,EFTOrder->5,LoopOrder->1]')


def test_preflight_only_requires_proof_but_does_not_run_wolfram(tmp_path,monkeypatch):
    root=tmp_path/'mixed_scalar_hard_batch6'
    f=fixture(root)
    monkeypatch.setattr(v8.base,'OUTPUT',root)
    monkeypatch.setattr(v8.base,'preflight',lambda model,*a:{'source_fingerprints':f.current['source_fingerprints'], 'canonical_eom':f.current['canonical_eom']})
    monkeypatch.setattr(v8.audit,'audit_model_run',lambda folder,model:copy.deepcopy(f.row))
    monkeypatch.setattr(v8.base,'run',lambda argv: (_ for _ in ()).throw(AssertionError('native run must not execute')))
    assert v8.main(['--preflight-only','--models',f.model])==0
    f.arrays.write_text('tampered')
    assert v8.main(['--preflight-only','--models',f.model])==2


def test_gated_launcher_passes_proof_environment_and_restores_functions(tmp_path,monkeypatch):
    root=tmp_path/'mixed_scalar_hard_batch6'
    f=fixture(root)
    monkeypatch.setattr(v8.base,'OUTPUT',root)
    monkeypatch.setattr(v8.base,'preflight',lambda model,*a:{'source_fingerprints':f.current['source_fingerprints'], 'canonical_eom':f.current['canonical_eom']})
    monkeypatch.setattr(v8.audit,'audit_model_run',lambda folder,model:copy.deepcopy(f.row))
    old_helper=v8.base.NATIVE_HELPER
    old_env=v8.base.native_subprocess_environment
    def fake_run(argv):
        assert '--run-id' in argv and 'mixed_hard_08' in argv
        assert v8.base.NATIVE_HELPER != old_helper
        env=v8.base.native_subprocess_environment({},f.current['canonical_eom'])
        assert Path(env['T3_BATCH8_PROOF_RECONSTRUCTED'])==f.files['reconstructed_candidate']
        assert env['T3_BATCH8_GATE_VERIFIED']=='1'
        assert 'T3_EXPERIMENTAL_GENERAL_QUADRATICS' not in env
        (root/'mixed_hard_08').mkdir()
        return 2  # A native failure is still propagated, not relabeled PASS.
    monkeypatch.setattr(v8.base,'run',fake_run)
    assert v8.main(['--models',f.model])==2
    assert v8.base.NATIVE_HELPER==old_helper
    assert v8.base.native_subprocess_environment==old_env
    gate=json.loads((root/'mixed_hard_08'/'component_gate.json').read_text())
    assert gate['source_component_run']=='mixed_hard_06'
    assert gate['physical_C5_certified'] is False


def test_generated_wolfram_helper_has_balanced_unquoted_delimiters():
    """Regression: injection must not swallow the outer Module's closing `]`."""
    script=v8.gated_helper()
    stack=[]
    in_string=False
    escaped=False
    comment_depth=0
    i=0
    opening='([{' 
    closing=')]} '
    while i<len(script):
        c=script[i]
        n=script[i+1] if i+1<len(script) else ''
        if comment_depth:
            if c=='(' and n=='*':
                comment_depth+=1; i+=2; continue
            if c=='*' and n==')':
                comment_depth-=1; i+=2; continue
            i+=1; continue
        if in_string:
            if escaped:escaped=False
            elif c=='\\':escaped=True
            elif c=='"':in_string=False
            i+=1; continue
        if c=='(' and n=='*':
            comment_depth+=1; i+=2; continue
        if c=='"':
            in_string=True;i+=1;continue
        if c in opening:
            stack.append(c)
        elif c in ')]}':
            assert stack, 'unexpected closing delimiter'
            assert opening.index(stack.pop())==')]}'.index(c), 'mismatched delimiter'
        i+=1
    assert not stack and not in_string and comment_depth==0
