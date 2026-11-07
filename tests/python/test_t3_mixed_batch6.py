"""Batch 6 offline algebra and native-wrapper contract regression tests.

Test mock native source only; these are not live Matchete matching tests.
"""
import numpy as np
import pytest

from Numerical.running.T3MixedBatch6Numerics import (
    NativeSourceUnsupported, safe_compiler, evaluate_hard,
)
from scripts.RunT3MixedScalarHardBatch6 import (
    HOOK, instrument, driver_anchors, validate_eom_sources,
)


def z(x):
    return {'real':float(np.real(x)),'imag':float(np.imag(x))}


def cmat(a):
    return [[z(x) for x in row] for row in a]


def fixture():
    m=np.diag([150_000.,180_000.,200_000.]).astype(complex)
    y1=np.diag([.12,.13,.14]).astype(complex)
    y2=np.diag([.22,.23,.24]).astype(complex)
    a={'model':'B_-1','status':'NumericalRenormalisableRunningOnly',
       'scales_gev':{'F':100_000.,'S':1000.},
       'uv_at_fermion_threshold':{'MF':cmat(m),'y1':cmat(y1),'y2':cmat(y2)},
       'eft1_at_fermion_threshold':{
           'lambdaT3':z(.003),
           'lambdaHHdagS2S2':z(.0009),
           'lambdaHHdagS1barS1bar':z(-.0007)}}
    mu2=complex(13000.,5000.)
    A=np.array([[1_800_000.,mu2],[mu2.conjugate(),2_400_000.]])
    w,U=np.linalg.eigh(A)
    params={'unitary_U':{'real':U.real.tolist(),'imag':U.imag.tolist()},
            'heavy_mass_squared_gev2':w.tolist()}
    return a,params


def synthetic_source():
    y1='Bar[Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]'
    y2='Bar[Coupling[y2, {Index[d$$2, Flavor], Index[d$$1, NFlavor]}, 0]]'
    return ('hbar*'+y1+'*'+y2+'*Coupling[lambdaT3, {}, 0]'
            '*(1/ϵ + 1 + Log[μbar2/Coupling[MSX, {}, 0]^2])'
            '/Coupling[MF, {}, 0]')


def test_synthetic_c5_three_heavy_masses_and_pole_scale_identity(tmp_path):
    inp,param=fixture()
    src=tmp_path/'hard_c5_candidate.txt';src.write_text(synthetic_source())
    result=evaluate_hard('B_-1',src,inp,param,conditional_quartic_map=False)
    assert result['checks']['three_nondegenerate_fermions_summed']
    assert result['checks']['native_hard_log_vs_pole_relative_residual']<1e-10
    assert result['checks']['first_threshold_inherited_pole_cancelled'] is False
    assert result['physical_C5_certified'] is False
    pole=np.array(result['hard_pole_residue_c5_gev_inverse']['real'])
    expected=np.array([(.12+.01*i)*(.22+.01*i)*.003 /
                       (16*np.pi**2 * [150000.,180000.,200000.][i]) for i in range(3)])
    assert np.allclose(np.diag(pole),expected,rtol=1e-12)
    assert np.max(np.abs(pole-np.diag(np.diag(pole))))<1e-28
    assert len(result['hard_finite_only_conditional_takagi_masses_ev'])==3


def test_full_source_evaluation_rejects_uncertified_extra_quartic(tmp_path):
    inp,param=fixture()
    s=synthetic_source().replace('Coupling[lambdaT3, {}, 0]',
              'Coupling[lambdaT3ExtraQ020020I1, {}, 0]')
    path=tmp_path/'source.txt';path.write_text(s)
    with pytest.raises(NativeSourceUnsupported,match='map'):
        evaluate_hard('B_-1',path,inp,param,conditional_quartic_map=False)
    out=evaluate_hard('B_-1',path,inp,param,conditional_quartic_map=True)
    assert out['quartic_map_is_conditional'] is True
    assert out['physical_C5_certified'] is False


def test_unsupported_source_fails_closed():
    with pytest.raises(NativeSourceUnsupported):
        safe_compiler(synthetic_source()+'+Coupling[unmappedNewCoupling, {}, 0]')
    with pytest.raises(NativeSourceUnsupported):
        safe_compiler(synthetic_source()+'+NCM[Field[NewScalarX]]')


def test_scalar_eigenbasis_not_unrotated_or_real_split():
    _,param=fixture()
    U=np.array(param['unitary_U']['real'])+1j*np.array(param['unitary_U']['imag'])
    eig=np.array(param['heavy_mass_squared_gev2'])
    A=np.array([[1800000.,13000+5000j],[13000-5000j,2400000.]])
    assert np.linalg.norm(U.conj().T@A@U-np.diag(eig))<1e-6
    assert abs(U[1,0])>0
    assert np.linalg.norm(U.conj().T@U-np.eye(2))<1e-12


def mock_driver(indent=' ',*,include_first_f_match=True):
    return ('scriptDir = DirectoryName @ ExpandFileName[$InputFileName];\n'
        'T3MEPerform[] := Module[{},\n'
        + indent + 'T3MECheckpoint["fermion-first-lagrangian-checks-returned"];\n'
        + indent + HOOK + '\n'
        + (indent+'T3MECheckpoint["fermion-first-tree-Match-started"];\n'
           +indent+'matchingResult=CheckAbort[RunT3Matching[lag,5,0],$Aborted];\n'
           if include_first_f_match else '')
        + '];\nrun = CheckAbort[old[]];\n')


@pytest.mark.parametrize('indent', [' ', '  ', '\t', '   '])
def test_instrument_accepts_driver_indentation_but_not_removed_match(tmp_path, indent):
    original=mock_driver(indent)
    anchors=driver_anchors(original)
    assert len(anchors)==4
    revised=instrument(original,tmp_path)
    assert revised.count('Return[T3B6Hard[lag]];')==1
    assert revised.index('Return[T3B6Hard[lag]];')<revised.index('RunT3Matching[lag,5,0]')
    assert 'T3_BATCH6_EFT1_EOM' in revised
    assert 'OriginalFermionSourceUnmodified' in revised
    assert 'Get[FileNameJoin[{scriptDir,"RunMatching.wl"}]]' in revised
    assert 'old[]' not in revised
    assert original == mock_driver(indent)  # original source unmodified
    with pytest.raises(ValueError,match='anchor hook'):
        instrument(original.replace(HOOK,'removed'),tmp_path)
    with pytest.raises(ValueError,match='matching call'):
        instrument(mock_driver(indent,include_first_f_match=False),tmp_path)


def mock_eom(count, *, mass_mode='unindexed'):
    # Real B/C Matchete archives use Coupling[MF, {}, 0], not an NFlavor index.
    mf = ('Coupling[MF, {}, 0]' if mass_mode == 'unindexed'
          else 'Coupling[MF, {Index[d$$1, NFlavor]}, 0]')
    single=('GammaCC*' + mf + '*'
            'Field[NewScalarX, Scalar]*Field[NewScalarY, Scalar]')
    return '( ' + ' + '.join([single]*count) + ' )' + ' '*6000


@pytest.mark.parametrize('model,eom_count,llss_count',[('B_-1',20,20),('C_-1',8,24)])
def test_canonical_eom_and_expanded_llss_have_distinct_valid_counts(model,eom_count,llss_count):
    eom=mock_eom(eom_count); llss=mock_eom(llss_count)
    report=validate_eom_sources(model,eom,llss,{'llss_term_count':llss_count})
    assert report['canonical_eom_GammaCC_terms']==eom_count
    assert report['canonical_eom_indexed_MF_factors']==0
    assert report['canonical_eom_unindexed_MF_factors']==eom_count
    assert report['expanded_llss_unindexed_MF_factors']==llss_count
    assert report['mf_source_mode']=='SavedUnindexedSingleMassTemplate' 
    assert report['expanded_llss_GammaCC_terms']==llss_count
    assert report['canonical_EOMEFT_used_for_next_threshold']
    with pytest.raises(ValueError,match='LLSS GammaCC count'):
        validate_eom_sources(model,eom,llss.replace('GammaCC','GAMMA',1),{'llss_term_count':llss_count})
    with pytest.raises(ValueError,match='LLSS/MF source inconsistent'):
        validate_eom_sources(model,eom.replace('Coupling[MF,','Coupling[MX,',1),llss,{'llss_term_count':llss_count})

@pytest.mark.parametrize('model,eom_count,llss_count',[('B_-1',20,20),('C_-1',8,24)])
def test_full_preflight_accepts_saved_b_c_layout_without_recomputing_f(
    tmp_path, monkeypatch, model, eom_count, llss_count
):
    """Simulate the repository file layout and strict saved F-stage provenance."""
    import json
    from scripts import RunT3MixedScalarHardBatch6 as b6
    sync = tmp_path / 'sync'
    running = tmp_path / 'running'
    driver = tmp_path / 'Lagrangian' / 'RunGeneralT3FermionFirstEFT1Tree.wl'
    driver.parent.mkdir(parents=True)
    driver.write_text('T3MENativeHHMassCG; DefineField[NewScalarX,Scalar];\n'+mock_driver(' '))
    monkeypatch.setattr(b6,'DRIVER',driver)
    ps=b6.paths(model,sync,running)
    A=np.array([[1_800_000, 13_000],[13_000,2_400_000]],complex)
    eig,U=np.linalg.eigh(A)
    c2=lambda x: {'real':x.real.tolist(),'imag':x.imag.tolist()}
    docs={
      'status':{'status':'CurrentRunningMassBasisFThresholdTreeLLSSExported'},
      'tree_report':{
        'status':'FermionFirstTreeLLSSExtractedNotSequentialEFT1',
        'integrated_heavy_fields':['NewFermion'],
        'scalar_light_fields':['NewScalarX','NewScalarY'],
        'llss_term_count':llss_count,'transition_representation':'EOMEFT'},
      'parameters':{
        'quartic_source_count':17,'yukawa_source_count':2,
        'verified_17_tensor_rotation':True,
        'dimensions_alpha':[2,2,1 if model.startswith('B') else 3,-1],
        'original_mass_squared_A_gev2':c2(A),
        'unitary_U':c2(U), 'heavy_mass_squared_gev2':eig.tolist()},
      'matchete':{
        'all_17_quartics_present':True,'mass_basis_ordinary_check_passed':True,
        'heavy_mass_basis_detailed_check_passed':True,
        'all_old_scalar_fields_removed':True,
        'both_yukawa_interactions_present':True},
      'mass_basis':{'model':model,'basis':{
        'tachyonic':False,'degenerate_basis_nonunique':False,
        'unitary_flavour_rotation_U':c2(U)}},
      'running':{'model':model,'status':'NumericalRenormalisableRunningOnly'},
    }
    for key,path in ps.items():
        if key=='driver': continue
        path.parent.mkdir(parents=True,exist_ok=True)
        if key=='canonical_eom':
            path.write_text(mock_eom(eom_count))
        elif key=='matched_llss':
            path.write_text(mock_eom(llss_count))
        else:
            path.write_text(json.dumps(docs[key]))
    status=b6.preflight(model,sync,running)
    assert status['status']=='Ready'
    assert status['eft1_source_audit']['canonical_eom_GammaCC_terms']==eom_count
    assert status['eft1_source_audit']['expanded_llss_GammaCC_terms']==llss_count
    assert 'matched_llss' in status['source_fingerprints']

@pytest.mark.parametrize('model,eom_count,llss_count',[('B_-1',20,20),('C_-1',8,24)])
def test_legacy_indexed_mf_source_distinguished_not_silently_mixed(model,eom_count,llss_count):
    indexed=mock_eom(eom_count,mass_mode='indexed')
    indexed_llss=mock_eom(llss_count,mass_mode='indexed')
    report=validate_eom_sources(model,indexed,indexed_llss,{'llss_term_count':llss_count})
    assert report['canonical_eom_indexed_MF_factors']==eom_count
    assert report['canonical_eom_unindexed_MF_factors']==0
    assert report['mf_source_mode']=='IndexedMFRequiresSingleMassConversion'
    with pytest.raises(ValueError,match='representations differ'):
        validate_eom_sources(model,indexed,mock_eom(llss_count),{'llss_term_count':llss_count})
    # MF with a different index representation is not accepted.
    bad=mock_eom(eom_count).replace('Coupling[MF, {}, 0]',
            'Coupling[MF, {Index[d$$1, Flavor]}, 0]',1)
    with pytest.raises(ValueError,match='LLSS/MF source inconsistent'):
        validate_eom_sources(model,bad,mock_eom(llss_count),{'llss_term_count':llss_count})


def test_native_helper_has_both_mass_modes_and_reports_inferred_lift():
    from scripts.RunT3MixedScalarHardBatch6 import NATIVE_HELPER
    assert 'nUnindexed=Count[saved,Matchete`Coupling[MF,{},0],Infinity]' in NATIVE_HELPER
    assert 'nIndexed=Count[saved,mfPattern,Infinity]' in NATIVE_HELPER
    assert 'SavedUnindexedSingleMassTemplate' in NATIVE_HELPER
    assert 'IndependentThreeMassNativeMatching' in NATIVE_HELPER
    assert 'nIndexed>0 && nUnindexed>0' in NATIVE_HELPER


@pytest.mark.parametrize('model', ['B_-1', 'C_-1'])
def test_native_runner_does_not_reenable_original_mixed_quadratic(model):
    """Regression: Batch 6 v3 failed before matching at baseline-build-started."""
    from scripts.RunT3MixedScalarHardBatch6 import (
        native_subprocess_environment, BASELINE_OPT_INS, NATIVE_HELPER, TAIL,
    )
    inherited_env={k:'1' for k in BASELINE_OPT_INS}
    inherited_env['UNRELATED_SETTING']='retained'
    canonical=f'output/general_observable/mixed_f_synchronised_v1/{model}/EOMEFT.inputform.txt'
    env=native_subprocess_environment(inherited_env, canonical)
    assert all(k not in env for k in BASELINE_OPT_INS)
    assert env['T3_BATCH6_EFT1_EOM']==canonical
    assert env['UNRELATED_SETTING']=='retained'
    assert all(inherited_env[k]=='1' for k in BASELINE_OPT_INS)
    # Both native steps must still occur: native rotation, then saved EOMEFT.
    assert 'T3MEPerform[]' not in NATIVE_HELPER
    assert 'PreserveIntegratedT3MassCouplings[{"S1","S2"}]' in NATIVE_HELPER
    assert 'Match[tree,EFTOrder->5,LoopOrder->1]' in NATIVE_HELPER
    assert 'If[!FileExistsQ[outPath],' in TAIL


def test_driver_anchor_keeps_original_both_matching_stages():
    from scripts.RunT3MixedScalarHardBatch6 import NATIVE_HELPER
    generated=instrument(mock_driver(' '), __import__('pathlib').Path.cwd())
    assert 'Return[T3B6Hard[lag]];' in generated
    assert generated.index('Return[T3B6Hard[lag]];')<generated.index('RunT3Matching[lag,5,0]')
    assert generated.count('T3B6Hard[lag]') == 1
    assert 'T3B6Hard[nativeUV_]' in generated
    assert 'PreserveIntegratedT3MassCouplings' in NATIVE_HELPER


def test_v5_uses_saved_eom_with_heavy_free_rebase_and_equivalence_gate():
    from scripts.RunT3MixedScalarHardBatch6 import NATIVE_HELPER, instrument
    from pathlib import Path
    generated = instrument(mock_driver(' '), Path.cwd())
    assert 'Return[T3B6Hard[lag]];' in generated
    assert 'nativeUV_' in NATIVE_HELPER
    assert 'lightFree=CheckAbort[FreeLag[NewScalarX,NewScalarY]' in NATIVE_HELPER
    assert 'heavyFree=CheckAbort[FreeLag[NewScalarX,NewScalarY]' in NATIVE_HELPER
    assert 'saved-lightFree+heavyFree' in NATIVE_HELPER
    assert 'T3B6EquivalentQ[rebased,reconstructed]' in NATIVE_HELPER
    assert 'CheckLagrangian[expr,DetailedOutput->True]' in NATIVE_HELPER
    assert 'Hermiticity' in NATIVE_HELPER
    assert 'CanonicallyNormalized' in NATIVE_HELPER
    assert 'GaugeInvariant' in NATIVE_HELPER
    assert 'AllObjectsDefined' in NATIVE_HELPER
    assert 'OriginalFermionSourceUnmodified' in NATIVE_HELPER
    assert '(L+HC[L])/2' in NATIVE_HELPER
    assert 'Match[tree,EFTOrder->5,LoopOrder->1]' in NATIVE_HELPER
    assert generated.index('Return[T3B6Hard[lag]];') < generated.index('RunT3Matching[lag,5,0]')


def test_v5_fail_closed_no_unverified_operator_projection():
    from scripts.RunT3MixedScalarHardBatch6 import NATIVE_HELPER
    # Never paper over a failed Matchete source check.
    assert 'ReconstructionNotAlgebraicallyIdentical' in NATIVE_HELPER
    assert 'EquivalentReconstructionFailedChecks' in NATIVE_HELPER
    assert 'full_EFT_validated' in NATIVE_HELPER
    assert 'full_pole_cancellation_certified"->False' in NATIVE_HELPER
    assert 'If[!TrueQ[first["passed"]] && chosen==="None",' in NATIVE_HELPER
    assert 'hard_scalar_eft1_validation.json' in NATIVE_HELPER
    assert 'hard_scalar_reconstruction_audit.json' in NATIVE_HELPER
