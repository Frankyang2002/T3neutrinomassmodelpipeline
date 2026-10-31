"""Fail-closed B/C mass-eigenstate Yukawa and saved one-scale C5 checks."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from scripts import RunT3BCMassEigenYukawaC5Handoff as gate


def matchete(model):
    return {'native_yukawa1_tensor':{
                'OldInvariantCount':1,'NewInvariantCount':1,
                'ExactTensorProjectionVerified':True,
                'NativeCG':'T3MENativeY1CG',
                'FermionSU2Dimension':1 if model=='B_-1' else 3,
                'OldToNewCGFactor':'1',
                'PhysicalCouplingMapCertified':False},
            'native_yukawa2_epsilon_tensor':{
                'OldInvariantCount':1,'NewInvariantCount':1,
                'ExactTensorProjectionVerified':True,
                'NativeCG':'T3MENativeY2CG',
                'FermionSU2Dimension':1 if model=='B_-1' else 3,
                'OldTimesEpsilonCGFactor':'-1' if model=='B_-1' else '1',
                'PhysicalCouplingMapCertified':False}}


def U_example():
    theta=0.16
    return np.array([[np.cos(theta),np.sin(theta)],
                     [-np.sin(theta),np.cos(theta)]],dtype=complex)


def test_yukawa_phase_b_vs_c():
    u=U_example()
    w_b=gate.yukawa_weights('B_-1',u,matchete('B_-1'))
    w_c=gate.yukawa_weights('C_-1',u,matchete('C_-1'))
    np.testing.assert_allclose(w_b['y1'],u[0],atol=1e-14)
    np.testing.assert_allclose(w_c['y1'],u[0],atol=1e-14)
    np.testing.assert_allclose(w_b['y2'],u[1],atol=1e-14)
    np.testing.assert_allclose(w_c['y2'],-u[1],atol=1e-14)
    assert w_b['y2_source_field_substitution_minus_sign']==-1
    assert w_c['y2_epsilon_CG_factor']==+1


@pytest.mark.parametrize('model',gate.MODELS)
def test_orthogonality_and_norm(model):
    u=U_example()
    w=gate.yukawa_weights(model,u,matchete(model))
    assert abs(np.vdot(w['y1'],w['y2']))<1e-14
    assert np.isclose(np.vdot(w['y1'],w['y1']),1)
    assert np.isclose(np.vdot(w['y2'],w['y2']),1)


def test_reject_wrong_epsilon_sign():
    doc=matchete('B_-1')
    doc['native_yukawa2_epsilon_tensor']['OldTimesEpsilonCGFactor']='1'
    with pytest.raises(ValueError,match='epsilon CG factor'):
        gate.yukawa_weights('B_-1',U_example(),doc)


def test_reject_not_verified_native_tensor():
    doc=matchete('C_-1')
    doc['native_yukawa1_tensor']['ExactTensorProjectionVerified']=False
    with pytest.raises(ValueError,match='projection unverified'):
        gate.yukawa_weights('C_-1',U_example(),doc)


def test_reject_wrong_fermion_representation():
    doc=matchete('B_-1')
    doc['native_yukawa1_tensor']['FermionSU2Dimension']=3
    with pytest.raises(ValueError,match='dimension changed'):
        gate.yukawa_weights('B_-1',U_example(),doc)


def test_reject_nonunitary_u():
    u=U_example();u[0,0]*=1.01
    with pytest.raises(ValueError,match='unitary'):
        gate.yukawa_weights('B_-1',u,matchete('B_-1'))


def endpoint_example(model='B_-1'):
    m0=np.array([[0.05,0.002,0],[0.002,0.01,0],[0,0,0.001]],dtype=complex)
    m1=m0+np.diag([0.005,0.004,0.0001])
    payload=lambda m:gate.pairing_matrix(m)
    return {'status':'ConditionalEFT1ScalarEndpointOneScaleC5NotPhysicalSequential',
            'model':model,'endpoint':'final',
            'physical_matchete_rgbeta_mapping_certified':False,
            'scalar_threshold_matching_performed':False,
            'sequential_EFT1_wilson_running_performed':False,
            'yukawas_evolved_with_candidate_scalar_background':False,
            'physical_neutrino_mass_prediction':False,
            'baseline_one_scale_mnu_ev':payload(m0),
            'candidate_one_scale_mnu_ev':payload(m1),
            'baseline_takagi_masses_ev_descending':np.linalg.svd(m0,compute_uv=False).tolist(),
            'candidate_takagi_masses_ev_descending':np.linalg.svd(m1,compute_uv=False).tolist(),
            'delta_mnu_frobenius_ev':float(np.linalg.norm(m1-m0)),
            'relative_delta_mnu_frobenius':float(np.linalg.norm(m1-m0)/np.linalg.norm(m0))}


def test_saved_candidate_one_scale_matrix_spectra():
    row=gate.check_one_scale(endpoint_example(),'B_-1')
    assert len(row['masses'])==3 and row['relative']>0


def test_reject_inconsistent_saved_spectrum():
    ep=endpoint_example();ep['candidate_takagi_masses_ev_descending'][1]*=1.1
    with pytest.raises(ValueError,match='Takagi mass spectrum mismatch'):
        gate.check_one_scale(ep,'B_-1')


def test_reject_wrong_saved_mnu_response():
    ep=endpoint_example();ep['relative_delta_mnu_frobenius']*=1.2
    with pytest.raises(ValueError,match='response norm mismatch'):
        gate.check_one_scale(ep,'B_-1')


def test_reject_false_physical_certification():
    ep=endpoint_example();ep['physical_neutrino_mass_prediction']=True
    with pytest.raises(ValueError,match='unexpectedly certified'):
        gate.check_one_scale(ep,'B_-1')


def cg_example():
    return {'format':'t3-bc-live-registered-cg-110101-cross-tool-identity-v1',
            'status':'LiveMatcheteRegisteredCGExactWithGatedRGBetaIdentityNoEFTMatching',
            'live_registered_matchete_components_compared':64,
            'signed_direct_and_conjugate_coefficients_compared':128,
            'live_matchete_registered_cg_component_identity_verified':True,
            'cross_tool_signed_component_identity_verified_in_exported_su2_basis':True,
            'physical_coupling_mapping_certified':False,
            'mass_basis_matching_certified':False,
            'sequential_eft_matching_certified':False,
            'direct_live_matchete_action_subtraction_performed':False,
            'models':[{'model':model,
                       'live_registered_components_equal_native_stage22_exact':32,
                       'live_matchete_registered_cg_direct_and_hc_coefficients_checked':64,
                       'exact_registered_component_disagreements':0} for model in gate.MODELS]}


def test_saved_live_quartic_gate():
    gate.check_live_cg(cg_example())


def test_reject_missing_live_quartic_component():
    row=cg_example();row['models'][0]['live_registered_components_equal_native_stage22_exact']=31
    with pytest.raises(ValueError,match='model check failed'):
        gate.check_live_cg(row)


def test_reject_claimed_live_matchete_action_subtraction():
    row=cg_example();row['direct_live_matchete_action_subtraction_performed']=True
    with pytest.raises(ValueError,match='overclaims certification'):
        gate.check_live_cg(row)


def test_complex_matrix_roundtrip():
    z=np.array([[1+2j,3-4j],[2j,-1j]])
    np.testing.assert_allclose(gate.complex_matrix(gate.pairing_matrix(z),'z'),z)


def test_reject_malformed_complex_pair():
    with pytest.raises(ValueError,match='missing complex pair'):
        gate.pair({'real':2},'z')


def test_no_overwrite_existing_run(tmp_path):
    # Even a valid in-memory report cannot reuse an existing output folder.
    # The live CG file is independently checked before any writing occurs.
    directory=tmp_path/'output/audits/t3_bc_matchete_registered_cg/bc_live_registered_cg_v1'
    directory.mkdir(parents=True)
    (directory/'summary.json').write_text(json.dumps(cg_example()))
    out=tmp_path/'existing';out.mkdir()
    with pytest.raises(ValueError):
        gate.run(tmp_path,out,(),recompute=False)


def test_only_existing_candidate_no_matching_flag():
    assert gate.STATUS=='ExistingNativeYukawaMassBasisAndConditionalOneScaleC5Checked'
    # This stage must not claim a new physical prediction.
    assert 'Matching' not in gate.STATUS


def test_full_readonly_handoff_synthetic_files(tmp_path, monkeypatch):
    """Exercise both files, U alignment, frozen flavour, and complete writer."""
    import types
    from hashlib import sha256
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis

    def write(p,data):
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps(data),encoding='utf-8')
    input_paths={}
    for model in gate.MODELS:
        tag=model.replace('-','m')
        base=tmp_path/model
        input_paths[model]={k:base/f'{k}.json' for k in
              ('matchete','parameters','one_loop','endpoint','mass_basis','frozen','frozen_manifest')}
        input_paths[model]['one_loop_expression']=base/'source.txt'
        p=input_paths[model]
        b=MixedDoubletMassBasis.build(model,1.15e6,1.38e6,60000.)
        U=b.U
        # The fixture source is generated independently of the handoff.
        p['mass_basis'].parent.mkdir(parents=True,exist_ok=True)
        p['mass_basis'].write_text('{"fixture":1}')
        p['one_loop_expression'].write_text('simple-matchete-expression')
        p['frozen'].write_text('{"fixture":2}')
        write(p['frozen_manifest'],{'models':[]})
        para={'format':'t3-experimental-mass-eigen-matchete-parameters-v1',
              'model':model,'endpoint':'final','quartic_source_count':17,
              'yukawa_source_count':2,
              'unitary_U':{'real':U.real.tolist(),'imag':U.imag.tolist()},
              'original_mass_squared_A_gev2':{'real':b.A.real.tolist(),
                                                'imag':b.A.imag.tolist()},
              'heavy_masses_gev':np.sqrt(b.eigenvalues_sq).tolist(),
              'source_eft1_scale_gev':1000.,
              'source_sha256':sha256(p['mass_basis'].read_bytes()).hexdigest()}
        mat=dict(matchete(model),
                 status='MatcheteMassEigenHeavyBasisPassedPhysicalMatchingNotPerformed',
                 model=model,mass_basis_ordinary_check_passed=True,
                 heavy_mass_basis_detailed_check_passed=True,
                 both_yukawa_interactions_present=True,
                 all_old_scalar_fields_removed=True,all_17_quartics_present=True,
                 physical_matching_performed=False,
                 rotation_spurion_values_from_numerical_U={
                   f'T3MEU{i+1}{j+1}':{'real':U[i,j].real,'imag':U[i,j].imag}
                   for i in range(2) for j in range(2)})
        one={'model':model,
             'status':'OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending',
             'one_loop_match_stage_completed':True,
             'physical_matching_performed':False,
             'production_sequential_matching_performed':False}
        ep=endpoint_example(model)
        ep.update({'scalar_masses_gev':np.sqrt(b.eigenvalues_sq).tolist(),
                   'one_scale_expression_sha256':sha256(p['one_loop_expression'].read_bytes()).hexdigest(),
                   'frozen_yukawa_source_sha256':sha256(p['frozen'].read_bytes()).hexdigest()})
        write(p['parameters'],para);write(p['matchete'],mat)
        write(p['one_loop'],one);write(p['endpoint'],ep)
    monkeypatch.setattr(gate,'paths',lambda root,model:input_paths[model])
    mod=types.ModuleType('scripts.RunGeneralT3ReferenceC5')
    mod.make_transfer=lambda *args:{'y1':gate.pairing_matrix(np.ones((3,1))),
                                   'y2':gate.pairing_matrix(2*np.ones((3,1))),
                                   'vev_gev':246.22}
    monkeypatch.setitem(sys.modules,'scripts.RunGeneralT3ReferenceC5',mod)
    loc=tmp_path/'output/audits/t3_bc_matchete_registered_cg/bc_live_registered_cg_v1'
    write(loc/'summary.json',cg_example())
    out=tmp_path/'handoff_v1'
    summary=gate.run(tmp_path,out,gate.MODELS)
    assert summary['models_checked']==2
    for model in gate.MODELS:
        output=json.loads((out/(model.replace('-','m')+'_yukawa_c5.json')).read_text())
        assert output['full_physical_neutrino_mass_prediction_certified'] is False
        assert output['Matchete_one_loop_coefficient_numerically_recomputed'] is False
        assert len(output['frozen_flavour_yukawa_matrices_mass_eigenbasis']['X']['y1_CConjF'])==3
        assert output['native_epsilon_projection_factor_y2']==gate.EXPECTED_Y2_CG[model]
    with pytest.raises(ValueError,match='refusing to overwrite'):
        gate.run(tmp_path,out,gate.MODELS)
