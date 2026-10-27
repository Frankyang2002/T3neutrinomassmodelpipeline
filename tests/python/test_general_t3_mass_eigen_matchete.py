"""Offline checks; these do not pretend to have run a Wolfram kernel."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
from scripts import RunGeneralT3MassEigenMatchete as runner

QUARTICS={
    'lambdaH','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12',
    'lambdaT3','lambdaH1Adj','lambdaH2Adj','lambda12Adj',
    'lambdaHHdagS2S2','lambdaHHdagS1barS1bar','lambdaS1bar2S2bar2',
    'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar',
    'lambdaHHdagS1barS2barCross','lambdaT3GeneralQ110101I1',
}


def sample(model:str, mu:complex=2500.0+1300.0j)->dict:
    basis=MixedDoubletMassBasis.build(model,1050.0**2,1190.0**2,mu)
    assert basis.eigenvalues_sq[0]>0
    return {
        'format':'t3-mixed-doublet-mass-basis-transport-v1',
        'status':'ExperimentalMassBasisTensorTransportNoPhysicalMatching',
        'model':model,'endpoint':'final',
        'source_quartic_count':17,'source_new_quartics':1,'source_reused_quartics':6,
        'MatcheteRGBetaPhysicalCouplingMapVerified':False,
        'physical_mass_matching_performed':False,
        'quartics':{q:{} for q in QUARTICS},
        'yukawa_scalar_leg_transforms':{'y1':{},'y2':{}},
        'basis':basis.info(),
        'max_source_tensor_transport_relative_error':1e-15,
        'quadratic_hessian_agreement':True,
        'scale_gev':1e3,
        'source_sha256':{},
    }


def fixtures(tmp_path:Path, monkeypatch, model:str, mu=2500.+1300.j):
    root=tmp_path/'project'
    (root/'RGE/general/generated').mkdir(parents=True)
    (root/'docs').mkdir(parents=True)
    source_sha={}
    for rel in ('RGE/general/generated/general_rgbeta_quartic_registry.json',
                'docs/general_scalar_trial_operator_manifest.json'):
        file=root/rel
        file.write_text('fixture registry/manifest',encoding='utf-8')
        source_sha[rel]=hashlib.sha256(file.read_bytes()).hexdigest()
    inputdir=root/'output/general_t3_mass_basis'
    inputdir.mkdir(parents=True)
    body=sample(model,mu=mu)
    body['source_sha256']=source_sha
    (inputdir/runner._model_file(model,'final')).write_text(json.dumps(body),encoding='utf-8')
    monkeypatch.setattr(runner,'ROOT',root)
    return body,inputdir,root/'output/general_t3_mass_eigen_matchete',root


@pytest.mark.parametrize('model',['B_-1','C_-1'])
@pytest.mark.parametrize('mu',[2200.+0j,0.+2100.j,2300.+1600.j,0j])
def test_prepare_valid_mixed_mass_eigenstates(tmp_path,monkeypatch,model,mu):
    _,source,target,_=fixtures(tmp_path,monkeypatch,model,mu)
    obj=runner.prepare_model(model,folder=source,outdir=target)
    p=obj['payload']; assert p['status']==runner.STATUS_PREPARED
    assert p['quartic_source_count']==17 and p['yukawa_source_count']==2
    assert p['heavy_masses_gev'][0]>0 and p['heavy_masses_gev'][1]>p['heavy_masses_gev'][0]
    assert p['threshold_plan']['mass_scales_gev']==list(reversed(p['heavy_masses_gev']))
    assert p['threshold_plan']['sequential_matching_performed'] is False
    assert p['matching_next_actions']
    assert p['physical_matching_performed'] is False
    assert json.loads(obj['parameters'].read_text())==p
    r=runner.run_model(model,folder=source,outdir=target,prepare_only=True)
    assert r['status']==runner.STATUS_PREPARED


@pytest.mark.parametrize('model',['B_-1','C_-1'])
@pytest.mark.parametrize('bad',(lambda v:v.update(status='Unverified'),
                                lambda v:v.update(source_quartic_count=16),
                                lambda v:v.update(physical_mass_matching_performed=True),
                                lambda v:v.update(endpoint='initial'),
                                lambda v:v.update(quartics={})))
def test_rejects_bad_upstream_status_or_inventory(tmp_path,monkeypatch,model,bad):
    body,source,target,_=fixtures(tmp_path,monkeypatch,model)
    bad(body)
    (source/runner._model_file(model,'final')).write_text(json.dumps(body))
    with pytest.raises(ValueError):runner.prepare_model(model,folder=source,outdir=target)


@pytest.mark.parametrize('model',['B_-1','C_-1'])
@pytest.mark.parametrize('entry',('unitary_flavour_rotation_U','mass_squared_flavour_matrix_gev2'))
def test_rejects_nonfinite_matrices(tmp_path,monkeypatch,model,entry):
    body,source,target,_=fixtures(tmp_path,monkeypatch,model)
    body['basis'][entry]['real'][0][0]=float('nan')
    (source/runner._model_file(model,'final')).write_text(json.dumps(body))
    with pytest.raises(ValueError):runner.prepare_model(model,folder=source,outdir=target)


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_rejects_changed_manifest_provenance(tmp_path,monkeypatch,model):
    body,source,target,root=fixtures(tmp_path,monkeypatch,model)
    (root/'docs/general_scalar_trial_operator_manifest.json').write_text('different')
    with pytest.raises(ValueError,match='manifest changed'):
        runner.prepare_model(model,folder=source,outdir=target)


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_rejects_wrong_eigenvectors(tmp_path,monkeypatch,model):
    body,source,target,_=fixtures(tmp_path,monkeypatch,model)
    body['basis']['unitary_flavour_rotation_U']['real'][0][0]+=0.1
    (source/runner._model_file(model,'final')).write_text(json.dumps(body))
    with pytest.raises(ValueError):runner.prepare_model(model,folder=source,outdir=target)


def valid_wolfram(model):
    return {'status':runner.STATUS_PASSED,'model':model,
        'mass_basis_ordinary_check_passed':True,
        'heavy_mass_basis_detailed_check_passed':True,
        'all_old_scalar_fields_removed':True,
        'all_17_quartics_present':True,
        'both_yukawa_interactions_present':True,
        'native_topology_epsilon_tensor':{
            'OriginalTimesEpsilonTensorExact':True,
            'SymmetricMassScalarCG':True,
            'OriginalHHScalarContractions':1,
            'NativeMassScalarContractions':1,
            'MappedRGBetaQuarticNormalizationCertified':False,
            'InvariantCG':'T3MENativeHHMassCG',
            'OriginalTimesEpsilonCGFactor':'1'
        },
        'native_H2_portal_epsilon_tensors':{
            'OldInvariantCount':2, 'NewInvariantCount':2,
            'TransportedRank':2, 'CombinedBasisRank':2,
            'ExactTensorProjectionVerified':True,
            'PhysicalCouplingMapCertified':False,
        },
        'native_cross_scalar_epsilon_tensors':{
            'OldInvariantCount':2,'NewInvariantCount':2,
            'TransportedRank':2,'CombinedBasisRank':2,
            'ExactTensorProjectionVerified':True,
            'PhysicalCouplingMapCertified':False,
        },
        'native_S2_self_epsilon_tensor':{
            'OldInvariantCount':1,'NewInvariantCount':1,
            'TransportedRank':1,'CombinedBasisRank':1,
            'ExactTensorProjectionVerified':True,
            'PhysicalCouplingMapCertified':False,
        },
        'native_yukawa2_epsilon_tensor':{
            'OldInvariantCount':1,'NewInvariantCount':1,
            'ExactTensorProjectionVerified':True,
            'OldTimesEpsilonCGFactor':'1',
            'FermionSU2Dimension': {'B_-1':1,'C_-1':3}[model],
            'NativeCG':'T3MENativeY2CG',
            'PhysicalCouplingMapCertified':False,
        },
        'native_yukawa1_tensor':{
            'OldInvariantCount':1,'NewInvariantCount':1,
            'ExactTensorProjectionVerified':True,'OldToNewCGFactor':'1',
            'FermionSU2Dimension':{'B_-1':1,'C_-1':3}[model],
            'NativeCG':'T3MENativeY1CG',
            'PhysicalCouplingMapCertified':False,
        },
        'native_additional_quartic_tensors':{
            'SourceOperatorCount':7,'RecoupledOperatorCount':7,
            'AllProjectedTensorsInvariant':True,
            'PhysicalCouplingMapCertified':False,
            'Certificates':{
                'T3ExtraQ_'+ '_'.join(map(str,counts))+'_'+str(inv):{
                    'SourceCoupling':'lambdaT3ExtraQ'+''.join(map(str,counts))+'I'+str(inv),
                    'SourceCounts':list(counts),
                    'SourceTensorIndex':inv,
                    'NativeCG':'T3MENativeExtraQ'+''.join(map(str,counts))+'I'+str(inv)+'CG',
                    'S2EpsilonProjectedLegs':counts[4]+counts[5],
                    'NativeInvariantSpanRank':2,
                    'CombinedBasisRank':2,
                    'ExactTensorProjectionVerified':True,
                }
                for counts,inv in (
                    ((0,0,0,1,1,2),1),((0,0,0,2,0,2),1),
                    ((0,0,1,2,0,1),1),((0,2,0,0,2,0),1),
                    ((0,2,0,2,0,0),1),((1,1,0,1,0,1),1),
                    ((1,1,0,1,0,1),2),
                )
            }
        },
        'matchete_rotation_uses_exact_symbolic_spurions':True,
        'rotation_spurion_values_from_numerical_U':{
            f'T3MEU{i}{j}':{'real':float(i+j),'imag':0.0}
            for i in (1,2) for j in (1,2)
        },
        'physical_matching_performed':False,
        'component_basis_map_certified':False,
        'epsilon_bridge_BB':{'phase_relative_to_matchete_component_epsilon_verified':True},
        'epsilon_bridge_UU':{'phase_relative_to_matchete_component_epsilon_verified':True}}


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_live_contract_success_flag_required(model):
    runner.validate_wolfram_output(valid_wolfram(model),model)
    for name in ('status','mass_basis_ordinary_check_passed','heavy_mass_basis_detailed_check_passed','all_old_scalar_fields_removed',
                 'all_17_quartics_present','both_yukawa_interactions_present',
                 'physical_matching_performed','component_basis_map_certified','model'):
        data=valid_wolfram(model)
        if name=='status':data[name]='Pretend'
        elif name=='model':data[name]='Other'
        else:data[name]=not data[name]
        with pytest.raises(ValueError):runner.validate_wolfram_output(data,model)
    for bridge in ('epsilon_bridge_BB','epsilon_bridge_UU'):
        data=valid_wolfram(model)
        data[bridge]={'phase_relative_to_matchete_component_epsilon_verified':False}
        with pytest.raises(ValueError):runner.validate_wolfram_output(data,model)


def test_mixed_model_selection_fails_closed():
    for model in ('A_+0','B_+1','C_+1','E_-2'):
        with pytest.raises(ValueError):runner._model_file(model,'final')


def test_mass_scale_boundary_safety(tmp_path,monkeypatch):
    model='B_-1'
    body,source,target,_=fixtures(tmp_path,monkeypatch,model)
    body['basis']['masses_squared_gev2'][0]=-1
    (source/runner._model_file(model,'final')).write_text(json.dumps(body))
    with pytest.raises(ValueError):runner.prepare_model(model,folder=source,outdir=target)


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_windows_backslash_source_hashes_are_accepted(tmp_path, monkeypatch, model):
    body, source, target, _ = fixtures(tmp_path, monkeypatch, model)
    body['source_sha256'] = {
        key.replace('/', chr(92)): digest
        for key, digest in body['source_sha256'].items()
    }
    (source/runner._model_file(model,'final')).write_text(json.dumps(body))
    result = runner.prepare_model(model, folder=source, outdir=target)
    assert result['payload']['status'] == runner.STATUS_PREPARED


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_windows_backslash_source_hashes_still_detect_tampering(tmp_path,monkeypatch,model):
    body, source, target, root = fixtures(tmp_path, monkeypatch, model)
    body['source_sha256'] = {
        key.replace('/', chr(92)): digest
        for key, digest in body['source_sha256'].items()
    }
    (source/runner._model_file(model,'final')).write_text(json.dumps(body))
    (root/'RGE/general/generated/general_rgbeta_quartic_registry.json').write_text('tampered')
    with pytest.raises(ValueError, match='manifest changed'):
        runner.prepare_model(model, folder=source, outdir=target)


def test_provenance_mixed_separators_normalize_and_duplicates_rejected():
    digest='a'*64
    result=runner._normalised_source_hashes({
        r'RGE\general/generated\registry.json': digest,
        'docs/file.json':'b'*64,
    })
    assert result == {'RGE/general/generated/registry.json':digest,'docs/file.json':'b'*64}
    with pytest.raises(ValueError, match='Duplicate cross-platform'):
        runner._normalised_source_hashes({'docs/file.json':digest,r'docs\file.json':digest})


@pytest.mark.parametrize('bad',[
    {'../docs/file.json':'a'*64},
    {'docs/./file.json':'a'*64},
    {'docs/file.json':'not-a-digest'},
    {'C:/docs/file.json':'a'*64},
])
def test_provenance_bad_paths_or_hashes_rejected(bad):
    with pytest.raises(ValueError):runner._normalised_source_hashes(bad)

# Regression for the live-Wolfram failure which previously exported the
# useless status UncaughtMatcheteError with no execution-stage diagnostic.
@pytest.mark.parametrize('with_progress',[False,True])
def test_wolfram_uncaught_error_reports_current_checkpoint(tmp_path,monkeypatch,with_progress):
    from types import SimpleNamespace
    model='B_-1'
    _,source,target,root=fixtures(tmp_path,monkeypatch,model)
    wolfram=root/'Lagrangian'/'RunGeneralT3MassEigenMatchete.wl'
    wolfram.parent.mkdir(parents=True)
    wolfram.write_text('Mock script')
    monkeypatch.setattr(runner,'WOLFRAM',wolfram)
    monkeypatch.setattr(runner.shutil,'which',lambda name: '/fake/wolframscript')

    name='B_m1_final_matchete'
    target.mkdir(parents=True,exist_ok=True)
    # Stale checkpoints must be discarded rather than misreported.
    (target/(name+'.progress.json')).write_text(json.dumps({'last_stage':'stale'}))

    def pretend_wolfram(command,**kwargs):
        result_file=Path(command[-1])
        result_file.write_text(json.dumps({'status':'UncaughtMatcheteError','model':model}))
        if with_progress:
            (target/(name+'.progress.json')).write_text(json.dumps({
                'last_stage':'rotate-original-interactions-started',
                'stage_history':['baseline-build-started','rotate-original-interactions-started']
            }))
        return SimpleNamespace(returncode=5,
            stdout='Matchete::mockerr: failed during rotation',stderr='Mock kernel failure')
    monkeypatch.setattr(runner.subprocess,'run',pretend_wolfram)
    with pytest.raises(ValueError) as err:
        runner.run_model(model,folder=source,outdir=target)
    msg=str(err.value)
    assert 'UncaughtMatcheteError' in msg
    assert 'Matchete::mockerr' in msg
    assert ('rotate-original-interactions-started' if with_progress else 'no-stage-checkpoint') in msg
    assert 'stale' not in msg


def test_wolfram_diagnostic_source_has_all_checkpoints():
    # Enforces the instrumentation contract without claiming to execute Matchete.
    source=(Path(__file__).resolve().parents[2]/'Lagrangian'/'RunGeneralT3MassEigenMatchete.wl').read_text()
    for checkpoint in (
        'script-initialized','baseline-build-started','extra-quartic-generation-started',
        'bridge-barred-started','bridge-unbarred-started','diagonal-FreeLag-started',
        'rotate-original-interactions-started',
        'ordinary-CheckLagrangian-started','detailed-CheckLagrangian-started',
    ):
        assert checkpoint in source
    assert 'stage_history' in source
    assert '"kernel_result_inputform"' in source


def test_wolfram_unitary_matrix_does_not_shadow_standard_model_u_field():
    """Regression: assigning a matrix to Global`u breaks SM DefineField[u]."""
    import re

    source = (Path(__file__).resolve().parents[2] /
              'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    # The original live failure was at LoadModel["SM"] with
    # DefineField::FieldLabel: the field label was a 2x2 numerical matrix.
    assert re.search(r'(?m)^\$T3MEUnitaryMatrix\s*=\s*Quiet\s*@\s*Check\[', source)
    assert not re.search(r'(?m)^\s*u\s*=', source)
    assert not re.search(r'(?<![A-Za-z0-9_`])u\s*\[\[', source)
    assert 'T3MECheckpoint["baseline-build-started"]' in source
    assert 'built=BuildT3Lagrangian[model]' in source


def test_wolfram_unitary_rename_applies_to_entire_scalar_rotation():
    """All uses of U must refer to the same input matrix after the rename."""
    source = (Path(__file__).resolve().parents[2] /
              'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    assert 'MatrixQ[$T3MEUnitaryMatrix,NumericQ]' in source
    assert 'ConjugateTranspose[$T3MEUnitaryMatrix].$T3MEUnitaryMatrix' in source
    assert 'ConjugateTranspose[$T3MEUnitaryMatrix].aMatrix.$T3MEUnitaryMatrix' in source
    rotate = source.split('T3MERotateInteraction[expr_] :=',1)[1].split('T3MEBuildNativeTopology[] :=',1)[0]
    for a,b in ((1,1),(1,2),(2,1),(2,2)):
        assert f'$T3MEInteractionMixingMatrix[[{a},{b}]]' in rotate
    assert '$T3MEUnitaryMatrix[[' not in rotate
    assert 'NewScalarX' in rotate and 'NewScalarY' in rotate


# Regression against the LIVE B_-1 and C_-1 failure at
# original-interaction-assembly-started: the old With/FreeQ helper
# returned False for present dynamically named scalar coupling heads.
# These are source contract tests, NOT a replacement for live Matchete.
def test_wolfram_detects_internal_matchete_coupling_labels():
    source = (Path(__file__).resolve().parents[2] /
              'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    region = source.split('T3MECouplingPresentQ[expr_, symbol_Symbol] :=', 1)[1].split(
        'T3MECouplingPresenceSelfTest[] :=', 1)[0]
    assert 'Verbatim[label]' in region
    assert '!FreeQ[expr' in region
    assert 'HoldPattern[head_Symbol[___]]' not in region
    assert 'Coupling[present, {}, 0]' in source
    assert 'Coupling[absent, {}, 0]' in source
    assert 'FreeQ[expr, HoldPattern[coupling[___]]]' not in source


def test_wolfram_rotates_internal_field_objects_not_shortcut_heads():
    source = (Path(__file__).resolve().parents[2] /
              'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    rotate = source.split('T3MERotateInteraction[expr_] :=',1)[1].split('T3MEPerform[] :=',1)[0]
    for label in ('NewScalar1', 'NewScalar2'):
        assert f'Field[{label}, Scalar, {{i_}}, rest___]' in rotate
        assert f'Bar[Field[{label}, Scalar, {{i_}}, rest___]]' in rotate
        assert f'HoldPattern[{label}[i_]]' not in rotate
    assert 'FreeQ[expr, Verbatim[NewScalar1] | Verbatim[NewScalar2]' in source
    assert 'T3MEScalarSignatureReport[before]' in source
    assert 'OriginalScalarLabelsAbsentInInternalSource' in source


def test_wolfram_coupling_head_selftest_catches_both_false_positives_and_negatives():
    source = (Path(__file__).resolve().parents[2] /
              'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    assert 'T3MECouplingPresenceSelfTest[] := Module[' in source
    assert 'TrueQ[T3MECouplingPresentQ[test, present]]' in source
    assert 'TrueQ[T3MECouplingPresentQ[internal, present]]' in source
    assert '!TrueQ[T3MECouplingPresentQ[test, absent]]' in source
    assert 'If[!TrueQ[T3MECouplingPresenceSelfTest[]]' in source


def test_wolfram_coupling_check_fails_closed_with_missing_symbol_diagnostic():
    source = (Path(__file__).resolve().parents[2] /
              'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    assert 'missingExtra = Select[extraSymbols,' in source
    assert 'Length[DeleteDuplicates[extraSymbols]] =!= 7 || missingExtra =!= {}' in source
    assert 'T3MEFail["AdditionalQuarticCouplingsAbsentInOriginalSource"' in source
    assert 'T3MECouplingPresentQ[rotated,#]' in source
    assert 'T3MECouplingPresentQ[sourceInt,y1]' in source
    assert 'T3MECouplingPresentQ[sourceInt,y2]' in source

# Regression for the live B_-1/C_-1 Matchete failure:
# DefineCG::invldIndType: j$... is not defined as a representation.
# In the already-internalized expression the old index is Index[_,T3BSMd2].
# Fresh bridge indices must use the same Matchete type, not bare Module locals.
def test_bridge_cg_fresh_indices_are_matchete_typed():
    source = (Path(__file__).resolve().parents[2] /
              'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    rotate = source.split('T3MERotateInteraction[expr_] :=', 1)[1].split('T3MEPerform[] :=', 1)[0]
    for kind in ('BB', 'UU'):
        typed = f'CG[T3MEBridge{kind},{{i,Index[j,T3BSMd2]}}]'
        untyped = f'CG[T3MEBridge{kind},{{i,j}}]'
        assert rotate.count(typed) == 2, f'{kind} must be typed in both X and Y legs'
        assert untyped not in rotate
    assert rotate.count('NewScalarX[j]') == 2
    assert rotate.count('NewScalarY[j]') == 2
    assert 'T3MEBridgeIndicesTypedQ[transformed]' in rotate
    assert 'Return[$Failed]' in rotate


def test_bridge_index_preflight_fails_closed_on_untyped_or_wrong_representation():
    source = (Path(__file__).resolve().parents[2] /
              'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    preflight = source.split('T3MEBridgeIndicesTypedQ[expr_] :=', 1)[1].split(
        'T3MERotateInteraction[expr_] :=', 1)[0]
    assert 'HoldPattern[CG[T3MEBridgeBB | T3MEBridgeUU, {a_, b_}]]' in preflight
    assert '!MatchQ[a, Index[_, T3BSMd2]]' in preflight
    assert '!MatchQ[b, Index[_, T3BSMd2]]' in preflight
    assert 'Return[False]' in preflight
    assert 'T3MEBridgeIndicesTypedQ[transformed]' in source


def test_typed_index_fix_preserves_physics_checks():
    source = (Path(__file__).resolve().parents[2] /
              'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    for token in ('T3MEBridge[True,T3MEBridgeBB]',
                  'T3MEBridge[False,T3MEBridgeUU]',
                  'T3MEOldFieldsAbsentQ[lag]',
                  'all_17_quartics_present',
                  'both_yukawa_interactions_present',
                  'CheckLagrangian[lag]',
                  'CheckLagrangian[lag,DetailedOutput->True]',
                  'NativeHeavyMassBasisCheckFailed',
                  '"physical_matching_performed"->False'):
        assert token in source


# Matchete's exact invariant algebra must never receive floating values as
# tensor coefficients. Original JSON U remains numerically validated first.
def test_mixing_entries_are_registered_as_matchete_symbolic_couplings():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    block=source.split('T3MEDefineMixingSpurions[] :=',1)[1].split(
        'T3MEMixingSpurionNumericalMap[] :=',1)[0]
    for i in (1,2):
        for j in (1,2):
            symbol=f'T3MEU{i}{j}'
            assert f'DefineCoupling[{symbol}, SelfConjugate -> False]' in block
            assert f'{symbol}[]' in block
    assert '$T3MEUnitaryMatrix' not in block
    assert 'T3MECheckpoint["register-mixing-spurions-started"]' in source
    assert 'spurionStatus = T3MEDefineMixingSpurions[]' in source
    assert 'SymbolicMixingSpurionRegistrationFailed' in source


def test_symbolic_couplings_instead_of_machinenumbers_in_interactions():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    rotate=source.split('T3MERotateInteraction[expr_] :=',1)[1].split(
        'T3MEBuildNativeTopology[] :=',1)[0]
    assert rotate.count('Bar[$T3MEInteractionMixingMatrix[[')==4
    assert rotate.count('$T3MEInteractionMixingMatrix[[')==8
    assert '$T3MEUnitaryMatrix' not in rotate
    assert 'Conjugate[$T3MEUnitaryMatrix' not in rotate
    assert 'T3MEMixingSpurionNumericalMap[]' in source
    assert '"rotation_spurion_values_from_numerical_U"' in source
    assert 'missingSpurions = Select[' in source
    assert 'MixingSpurionsMissingAfterRotation' in source


def test_numeric_mass_diagonalization_still_checked_before_symbolic_rotation():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    assert 'ConjugateTranspose[$T3MEUnitaryMatrix].aMatrix.$T3MEUnitaryMatrix' in source
    assert 'unitarityError > 2 10^-10' in source
    assert 'massError > 2 10^-9 matScale' in source
    assert '"component_basis_map_certified"->False' in source
    assert '"physical_matching_performed"->False' in source


def test_failing_full_check_partitions_baseline_without_regenerating_cgs():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    assert 'T3MEIsolateNonInvariantSource[built_Association, free_, extra_List, nativeTopology_, nativeH2_Association, nativeS2_Association, nativeY2_Association, nativeY1_Association, nativeExtra_Association]' in source
    assert 'T3MEPartitionBaselineInteractions[built["LInt"]]' in source
    assert 'T3MEBaselineTermKey[term_]' in source
    assert '"Name" -> ("BaselinePart" <>' in source
    assert 'Length[baseline] < 8' in source
    assert 'Expand[Total[Lookup[parts, "Expression", {}]]-expanded] === 0' in source
    assert 'T3MECheckpoint["failure-isolation-started"]' in source
    assert 'T3MECheckpoint["failure-isolation-operator-" <> ToString[k]]' in source
    assert 'CheckLagrangian[(built["LSM"]+free) // RelabelIndices]' in source
    assert '"DiagonalFreeOnlyCheckFailed"' in source
    assert 'T3MECheckpoint["final-index-preflight-started"]' in source
    assert 'FinalEigenstateTensorIndexPreflightFailed' in source
    assert '"source_signature" -> Lookup[entry,"Signature","ExtraQuartic"]' in source
    assert '"isolated_source_failure"->failureDetail' in source
    assert 'Return[Join[results,<|"status"->"NativeHeavyMassBasisCheckFailed"' in source
    assert 'BuildT3InteractionCandidates[model]' not in source


def test_existing_typed_index_is_not_double_wrapped_by_public_accessor():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    assert 'T3MENewTypedScalar[label_Symbol, index_Index]' in source
    assert 'Field[label, Scalar, {index}, {}]' in source
    rotate=source.split('T3MERotateInteraction[expr_] :=',1)[1].split(
        'T3MEBaselineTermKey[term_]',1)[0]
    assert 'T3MENewTypedScalar[NewScalarX,i]' in rotate
    assert 'T3MENewTypedScalar[NewScalarY,i]' in rotate
    assert 'NewScalarX[i]' not in rotate
    assert 'NewScalarY[i]' not in rotate
    assert 'NewScalarX[j]' in rotate and 'NewScalarY[j]' in rotate
    assert 'T3MENewScalarIndexAudit[transformed]' in rotate
    assert 'Rotated Field eigen-scalar index audit failed' in rotate


def test_eigenfield_index_guard_checks_representation_and_nested_wrapper():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    audit=source.split('T3MENewScalarIndexAudit[expr_] :=',1)[1].split(
        'T3MEDefineMixingSpurions[] :=',1)[0]
    assert 'Field[(NewScalarX | NewScalarY), Scalar, inds_, rest___]' in audit
    assert 'MatchQ[#, {Index[name_, T3BSMd2]} /; Head[name] =!= Index]' in audit
    assert 'fields === {} || bad =!= {}' in audit
    assert 'Return[False]' in audit
    assert 'T3MEBridgeIndicesTypedQ[transformed]' in source
    assert 'T3MEOldFieldsAbsentQ[transformed]' in source


def test_bridge_orientation_and_typed_indices_not_weakend():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    for token in ('T3MEBridge[True,T3MEBridgeBB]',
                  'T3MEBridge[False,T3MEBridgeUU]',
                  'T3MEBridgeIndicesTypedQ[transformed]',
                  'If[!T3MEOldFieldsAbsentQ[transformed]',
                  'If[!T3MEOldFieldsAbsentQ[lag]',
                  'CheckLagrangian[lag]',
                  'CheckLagrangian[lag,DetailedOutput->True]'):
        assert token in source


def test_wolfram_output_rejects_unsupported_spurion_claim_or_incomplete_map():
    baseline=valid_wolfram('B_-1')
    expected={'real':[[2.,3.],[3.,4.]],'imag':[[0.,0.],[0.,0.]]}
    runner.validate_wolfram_output(baseline,'B_-1',expected)
    for broken in (
        {'matchete_rotation_uses_exact_symbolic_spurions':False},
        {'rotation_spurion_values_from_numerical_U':{}},
        {'rotation_spurion_values_from_numerical_U':{
            **baseline['rotation_spurion_values_from_numerical_U'],
            'T3MEU11':{'real':float('nan'),'imag':0.0}
        }},
    ):
        payload=valid_wolfram('B_-1');payload.update(broken)
        with pytest.raises(ValueError): runner.validate_wolfram_output(payload,'B_-1',expected)
    with pytest.raises(ValueError,match='does not match prepared'):
        runner.validate_wolfram_output(baseline,'B_-1',{'real':[[9.,3.],[3.,4.]],'imag':[[0.,0.],[0.,0.]]})


def test_native_topology_must_have_exact_tensor_certificate():
    for key,value in (
        ('OriginalTimesEpsilonTensorExact',False),
        ('SymmetricMassScalarCG',False),
        ('OriginalHHScalarContractions',2),
        ('NativeMassScalarContractions',0),
        ('MappedRGBetaQuarticNormalizationCertified',True),
        ('InvariantCG','IncorrectCG')):
        payload=valid_wolfram('B_-1')
        payload['native_topology_epsilon_tensor'][key]=value
        with pytest.raises(ValueError,match='native invariant'):
            runner.validate_wolfram_output(payload,'B_-1')
    payload=valid_wolfram('C_-1')
    payload.pop('native_topology_epsilon_tensor')
    with pytest.raises(ValueError,match='native invariant'):
        runner.validate_wolfram_output(payload,'C_-1')


def test_topology_is_rebuilt_from_native_hhxx_cg_and_not_composite_bridges():
    source=(Path(__file__).resolve().parents[2] /
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    assert 'T3MEBuildNativeTopology[] := Module[' in source
    region=source.split('T3MEBuildNativeTopology[] := Module[',1)[1].split(
        'T3MESeparateOriginalTopology[expr_]',1)[0]
    assert '{False,False,False,True}' in region
    assert '{False,False,False,False}' in region
    assert 'DefineCG[T3MENativeHHMassCG' in region
    assert 'projected == alphaCG newT' in region
    assert 'Transpose[newT,{1,2,4,3}]' in region
    assert 'u11 u21 NewScalarX[k] NewScalarX[m]' in region
    assert '(u11 u22 + u12 u21) NewScalarX[k] NewScalarY[m]' in region
    assert 'u12 u22 NewScalarY[k] NewScalarY[m]' in region
    assert 'CG[T3MENativeHHMassCG,{i,j,k,m}]' in region
    assert 'T3MECheckpoint["native-topology-basis-started"]' in source
    assert 'nativeTopology["Expression"]' in source
    assert 'T3MESeparateOriginalTopology[sourceInt]' in source
    assert 'T3MESeparateOriginalYukawa1[yukawaParts["RemainingInteractions"]]' in source
    assert 'T3MECouplingPresentQ[entry["Expression"],lambdaT3]' in source


def test_topology_reconstruction_preserves_inventory_and_matching_flags():
    source=(Path(__file__).resolve().parents[2] /
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text(encoding='utf-8')
    assert 'Total[selected]+remainder-expanded' in source
    assert 'native_topology_epsilon_tensor' in source
    assert '"physical_matching_performed"->False' in source
    assert '"component_basis_map_certified"->False' in source
    assert 'T3MEFail["NativeTopologyTensorConstructionFailed"]' in source



@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_native_H2_portal_certificate_is_checked(model):
    runner.validate_wolfram_output(valid_wolfram(model), model)
    for key, bad in (
        ('OldInvariantCount',1),
        ('NewInvariantCount',3),
        ('TransportedRank',1),
        ('CombinedBasisRank',3),
        ('ExactTensorProjectionVerified',False),
        ('PhysicalCouplingMapCertified',True)):
        payload = valid_wolfram(model)
        payload['native_H2_portal_epsilon_tensors'][key]=bad
        with pytest.raises(ValueError,match='H2 portal tensor'):
            runner.validate_wolfram_output(payload,model)
    payload = valid_wolfram(model)
    payload.pop('native_H2_portal_epsilon_tensors')
    with pytest.raises(ValueError,match='H2 portal tensor'):
        runner.validate_wolfram_output(payload,model)


def test_native_H2_portals_replace_both_CG_structures_exactly():
    src = (Path(__file__).resolve().parents[2]/
           'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    assert 'T3MEBuildNativeH2Portals[] := Module[' in src
    region=src.split('T3MEBuildNativeH2Portals[] := Module[',1)[1].split(
        'T3MESeparateOriginalH2Portals[expr_]',1)[0]
    assert '{True,False,True,False}' in region
    assert '{True,False,False,True}' in region
    assert 'eps[[m,k]] eps[[n,l]]' in region
    assert 'DefineCG[cgNames[[inv]],newData["CGReps"],projected[[inv]]]' in region
    assert '{originalRank,newRank,unionRank}=!={2,2,2}' in region
    assert 'coefficients[[a]] Bar[coefficients[[b]]]' in region
    assert 'Bar[H[i]] H[j] CG[cg,{i,j,k,l}]' in region
    assert '"lambdaH2Inv1","lambdaH2Inv2"' in region
    assert 'T3MECheckpoint["native-H2-portals-started"]' in src
    assert 'nativeH2=T3MEBuildNativeH2Portals[]' in src
    assert 'portalsRotated=Total[Values[nativeH2["Expressions"]]]' in src
    assert 'T3MERotateInteraction[extraParts["RemainingInteractions"]]' in src
    assert 'native_H2_portal_epsilon_tensors' in src


def test_original_H2_source_terms_are_not_dropped():
    src = (Path(__file__).resolve().parents[2]/
           'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    region=src.split('T3MESeparateOriginalH2Portals[expr_] := Module[',1)[1].split(
        '(* Diagnose gauge-tensor failures',1)[0]
    assert 'required={lambdaH2Inv1,lambdaH2Inv2}' in region
    assert 'Expand[Total[portalTerms]+rest-expanded]===0' in region
    assert 'T3MECouplingPresentQ[Total[portalTerms],#]' in region
    assert 'AnyTrue[required,T3MECouplingPresentQ[rest,#]&]' in region
    assert 'T3MEFail["OriginalH2PortalSeparationFailed"]' in src
    assert 'T3MEFail["NativeH2PortalTensorConstructionFailed"]' in src


def test_isolation_uses_native_H2_for_both_portals():
    src = (Path(__file__).resolve().parents[2]/
           'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    region=src.split('T3MEIsolateNonInvariantSource[built_Association',1)[1].split(
        'T3MEPerform[] :=',1)[0]
    assert 'nativeH2["Expressions"]["lambdaH2Inv1"]' in region
    assert 'nativeH2["Expressions"]["lambdaH2Inv2"]' in region
    assert 'T3MERotateInteraction[entry["Expression"]]' in region


@pytest.mark.parametrize('model',['B_-1','C_-1'])
@pytest.mark.parametrize('label,n',[
    ('native_cross_scalar_epsilon_tensors',2),
    ('native_S2_self_epsilon_tensor',1),
])
def test_native_S2_scalar_quartic_certificates_fail_closed(model,label,n):
    payload=valid_wolfram(model)
    runner.validate_wolfram_output(payload,model)
    payload.pop(label)
    with pytest.raises(ValueError,match='S2 scalar-quartic tensor'):
        runner.validate_wolfram_output(payload,model)
    for key,value in [
        ('OldInvariantCount',n+1),
        ('NewInvariantCount',n+1),
        ('TransportedRank',n+1),
        ('CombinedBasisRank',n+1),
        ('ExactTensorProjectionVerified',False),
        ('PhysicalCouplingMapCertified',True),
    ]:
        payload=valid_wolfram(model)
        payload[label][key]=value
        with pytest.raises(ValueError,match='S2 scalar-quartic tensor'):
            runner.validate_wolfram_output(payload,model)


def test_native_S2_scalar_quartics_replace_all_three_source_cgs():
    src=(Path(__file__).resolve().parents[2]/
         'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    region=src.split('T3MEBuildNativeS2ScalarQuartics[] := Module[',1)[1].split(
        'T3MESeparateOriginalS2ScalarQuartics[expr_]',1)[0]
    assert '{True,False,True,False}' in region
    assert '{True,False,False,True}' in region
    assert '{False,True,False,True}' in region
    assert 'eps[[m,k]] eps[[n,l]]' in region
    assert 'eps[[m,i]] eps[[n,j]]' in region
    assert 'eps[[o,k]] eps[[p,l]]' in region
    assert 'DefineCG[crossCG[[q]],newCross["CGReps"],' in region
    assert 'DefineCG[selfCG,newSelf["CGReps"],' in region
    assert '{2,2,2,1,1,1}' in region
    assert '"lambda12Inv1","lambda12Inv2","lambdaS2Inv1"' in region
    assert 'nativeS2=T3MEBuildNativeS2ScalarQuartics[]' in src
    assert 'scalarRotated=Total[Values[nativeS2["Expressions"]]]' in src
    assert 'T3MECheckpoint["native-S2-scalar-quartics-started"]' in src
    assert 'T3MERotateInteraction[extraParts["RemainingInteractions"]]' in src
    assert 'native_cross_scalar_epsilon_tensors' in src
    assert 'native_S2_self_epsilon_tensor' in src


def test_native_S2_scalar_source_separation_preserves_all_terms():
    src=(Path(__file__).resolve().parents[2]/
         'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    region=src.split('T3MESeparateOriginalS2ScalarQuartics[expr_] := Module[',1)[1].split(
        '(* Diagnose gauge-tensor failures',1)[0]
    assert 'required={lambda12Inv1,lambda12Inv2,lambdaS2Inv1}' in region
    assert 'Expand[Total[selected]+rest-expanded]===0' in region
    assert 'AnyTrue[required,T3MECouplingPresentQ[rest,#]&]' in region
    assert 'T3MEFail["OriginalS2ScalarQuarticSeparationFailed"]' in src
    assert 'T3MEFail["NativeS2ScalarQuarticTensorConstructionFailed"]' in src


def test_isolation_uses_native_S2_scalar_quartics():
    src=(Path(__file__).resolve().parents[2]/
         'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    region=src.split('T3MEIsolateNonInvariantSource[built_Association',1)[1].split(
        'T3MEPerform[] :=',1)[0]
    for label in ('lambda12Inv1','lambda12Inv2','lambdaS2Inv1'):
        assert 'nativeS2["Expressions"]["'+label+'"]' in region


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_native_yukawa2_certificate_required(model):
    runner.validate_wolfram_output(valid_wolfram(model),model)
    for key,bad in (
        ('OldInvariantCount',0),('NewInvariantCount',2),
        ('ExactTensorProjectionVerified',False),
        ('FermionSU2Dimension',2),('NativeCG','OldY2CG'),
        ('OldTimesEpsilonCGFactor',''),
        ('PhysicalCouplingMapCertified',True)):
        payload=valid_wolfram(model)
        payload['native_yukawa2_epsilon_tensor'][key]=bad
        with pytest.raises(ValueError,match='Yukawa2 tensor'):
            runner.validate_wolfram_output(payload,model)
    payload=valid_wolfram(model)
    payload.pop('native_yukawa2_epsilon_tensor')
    with pytest.raises(ValueError,match='Yukawa2 tensor'):
        runner.validate_wolfram_output(payload,model)


def test_yukawa2_new_cg_complete_replacement():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    section=source.split('T3MEBuildNativeYukawa2[model_Association] := Module[',1)[1].split(
        'T3MESeparateOriginalYukawa2[expr_]',1)[0]
    assert '{2,dF,2},{True,False,True}' in section
    assert '{2,dF,2},{True,False,False}' in section
    assert 'oldT[[i,m]] eps[[m,k]]' in section
    assert 'oldT[[i,j,m]] eps[[m,k]]' in section
    assert 'projected==alphaCG newT' in section
    assert 'DefineCG[T3MENativeY2CG,newData["CGReps"],newT]' in section
    assert 'If[dF===1,NewFermion[r],NewFermion[j,r]]' in section
    assert 'If[dF===1,{i,k},{i,j,k}]' in section
    assert 'y2[p,r] NCM[Bar[l[i,p]],fermion]' in section
    assert 'u21 NewScalarX[k]+u22 NewScalarY[k]' in section
    assert 'T3MECheckpoint["native-Yukawa2-basis-started"]' in source
    assert 'T3MESeparateOriginalYukawa2[' in source
    assert 'T3MESeparateOriginalYukawa1[yukawaParts["RemainingInteractions"]]' in source
    assert 'nativeTopology["Expression"]+nativeY2["Expression"]+' in source
    assert 'native_yukawa2_epsilon_tensor' in source


def test_yukawa2_source_separation_is_lossless_and_isolates():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    section=source.split('T3MESeparateOriginalYukawa2[expr_] := Module[',1)[1].split(
        '(* Diagnose gauge-tensor failures',1)[0]
    assert 'Total[selected]+rest-expanded' in section
    assert 'T3MECouplingPresentQ[rest,y2]' in section
    assert '"OriginalYukawa2SeparationFailed"' in source
    assert '"NativeYukawa2TensorConstructionFailed"' in source
    isolator=source.split('T3MEIsolateNonInvariantSource[built_Association',1)[1].split('T3MEPerform[] :=',1)[0]
    assert 'T3MECouplingPresentQ[entry["Expression"],y2]' in isolator
    assert 'nativeY2["Expression"]' in isolator
    assert 'nativeS2,nativeY2,nativeY1,nativeExtra];' in source


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_native_yukawa1_certificate_required(model):
    runner.validate_wolfram_output(valid_wolfram(model),model)
    for key,bad in (
        ('OldInvariantCount',0),('NewInvariantCount',2),
        ('ExactTensorProjectionVerified',False),
        ('FermionSU2Dimension',2),('NativeCG','T3Y1CG'),
        ('OldToNewCGFactor',''),('PhysicalCouplingMapCertified',True)):
        payload=valid_wolfram(model)
        payload['native_yukawa1_tensor'][key]=bad
        with pytest.raises(ValueError,match='Yukawa1 tensor'):
            runner.validate_wolfram_output(payload,model)
    payload=valid_wolfram(model)
    payload.pop('native_yukawa1_tensor')
    with pytest.raises(ValueError,match='Yukawa1 tensor'):
        runner.validate_wolfram_output(payload,model)


def test_yukawa1_reconstructed_from_original_complete_yukawa():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    region=source.split('T3MEBuildNativeYukawa1[model_Association] := Module[',1)[1].split(
        'T3MESeparateOriginalYukawa1[expr_]',1)[0]
    assert region.count('{2,dF,2},{True,True,False}') == 2
    assert 'oldT==alphaCG newT' in region
    assert 'DefineCG[T3MENativeY1CG,newData["CGReps"],newT]' in region
    assert 'CConj[NewFermion[r]]' in region
    assert 'CConj[NewFermion[j,r]]' in region
    assert 'y1[p,r] NCM[Bar[l[i,p]],fermion]' in region
    assert 'u11 NewScalarX[k]+u12 NewScalarY[k]' in region
    assert 'T3MECheckpoint["native-Yukawa1-basis-started"]' in source
    assert '"OriginalYukawa1SeparationFailed"' in source
    assert '"NativeYukawa1TensorConstructionFailed"' in source
    assert 'native_yukawa1_tensor' in source


def test_yukawa1_full_source_term_coverage():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    region=source.split('T3MESeparateOriginalYukawa1[expr_] := Module[',1)[1].split(
        '(* General direct-CG transformation',1)[0]
    for token in ('Total[selected]+rest-expanded','T3MECouplingPresentQ[rest,y1]',
                  'T3MECouplingPresentQ[Total[selected],y1]'):
        assert token in region
    isolator=source.split('T3MEIsolateNonInvariantSource[built_Association',1)[1].split('T3MEPerform[] :=',1)[0]
    assert 'T3MECouplingPresentQ[entry["Expression"],y1]' in isolator
    assert 'nativeY1["Expression"]' in isolator


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_native_extra_quartic_certificates_strict(model):
    runner.validate_wolfram_output(valid_wolfram(model),model)
    for key,bad in (('SourceOperatorCount',6),('RecoupledOperatorCount',8),
                    ('AllProjectedTensorsInvariant',False),
                    ('PhysicalCouplingMapCertified',True)):
        data=valid_wolfram(model)
        data['native_additional_quartic_tensors'][key]=bad
        with pytest.raises(ValueError,match='additional-quartic'):
            runner.validate_wolfram_output(data,model)
    for key,bad in (('ExactTensorProjectionVerified',False),('CombinedBasisRank',3),
                    ('S2EpsilonProjectedLegs',0),('SourceTensorIndex',99),
                    ('SourceCoupling','wrong')):
        data=valid_wolfram(model)
        name='T3ExtraQ_0_0_0_1_1_2_1'
        data['native_additional_quartic_tensors']['Certificates'][name][key]=bad
        with pytest.raises(ValueError,match='additional-quartic'):
            runner.validate_wolfram_output(data,model)
    data=valid_wolfram(model)
    data['native_additional_quartic_tensors']['Certificates'].pop('T3ExtraQ_1_1_0_1_0_1_2')
    with pytest.raises(ValueError,match='additional-quartic'):
        runner.validate_wolfram_output(data,model)


def test_extra_quartics_use_one_native_cg_per_complete_source():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    section=source.split('T3MEBuildNativeExtraQuartics[model_Association,extra_List] := Module[',1)[1].split(
        'T3MESeparateOriginalExtraQuartics[expr_,symbols_List]',1)[0]
    for token in ('T3ExtendedBosonicTensorBasis[model,names]',
                  'T3ExtendedSpeciesList[counts]',
                  'PhysicalSU2InvariantBasis[',
                  'T3MEExtraProjectedTensorSlot[#1,#2]&',
                  'newData["CGReps"]',
                  'rankUnion=MatrixRank[',
                  'DefineCG[cgName,cgReps,projected]',
                  'PlusHc[coupling[] CG[cgName,indices] product]',
                  '"AllProjectedTensorsInvariant"->True',
                  '"PhysicalCouplingMapCertified"->False'):
        assert token in section
    assert 'T3MECheckpoint["native-extra-quartics-started"]' in source
    assert '"OriginalExtraQuarticSeparationFailed"' in source
    assert '"NativeExtraQuarticTensorConstructionFailed"' in source
    assert 'native_additional_quartic_tensors' in source
    assert 'nativeExtra["Expressions"][entry["Name"]]' in source


def test_extra_quartic_lossless_separation_and_epsilon_rules():
    source=(Path(__file__).resolve().parents[2]/
            'Lagrangian/RunGeneralT3MassEigenMatchete.wl').read_text()
    body=source.split('T3MESeparateOriginalExtraQuartics[expr_,symbols_List] := Module[',1)[1].split(
        '(* Diagnose gauge-tensor failures',1)[0]
    assert 'Total[selected]+rest-expanded' in body
    assert 'AnyTrue[symbols,T3MECouplingPresentQ[part,#]&]' in body
    assert 'Length[DeleteDuplicates[symbols]] =!= 7' in body
    factors=source.split('T3MEExtraMassField[species_String,idx_] :=',1)[1].split(
        'T3MEBuildNativeExtraQuartics[',1)[0]
    for token in ('"S1",u11 NewScalarX[idx]+u12 NewScalarY[idx]',
                  '"S1dag",Bar[u11] Bar[NewScalarX[idx]]+Bar[u12] Bar[NewScalarY[idx]]',
                  '"S2",-Bar[u21] Bar[NewScalarX[idx]]-Bar[u22] Bar[NewScalarY[idx]]',
                  '"S2dag",-u21 NewScalarX[idx]-u22 NewScalarY[idx]'):
        assert token in factors
