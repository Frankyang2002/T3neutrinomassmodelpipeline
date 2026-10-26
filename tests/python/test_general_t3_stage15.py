"""Stage 15: source-convention contract and conditional numeric bridge."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.experimental.GeneralT3C5Pilot import ResponseKernel
from Numerical.experimental.GeneralT3ConditionalUVBridge import conditional_uv_to_c5
from Numerical.experimental.GeneralT3LagrangianConvention import (
    _check_matchete_source,
    _check_rgbeta_source,
    _check_tensor_report,
    derive_lagrangian_convention,
)
from Numerical.experimental.GeneralT3SourceMappedBridge import source_mapped_uv_to_c5


ACTIVE = {'lambdaT3ExtraQ020020I1': 'lambdaHHdagS2S2',
          'lambdaT3ExtraQ020200I1': 'lambdaHHdagS1barS1bar'}


def _z(z):
    return {'real': float(complex(z).real), 'imag': float(complex(z).imag)}


def _fixture(model='B_m1'):
    params = {
        'MF': [1500., 2100., 2800.],
        'MS1': 1350., 'MS2': 1780., 'scale_uv_gev': 20000.,
        'y1_real': (np.eye(3)*1e-5).tolist(),
        'y1_imag': np.zeros((3,3)).tolist(),
        'y2_real': (np.eye(3)*2e-5).tolist(),
        'y2_imag': np.zeros((3,3)).tolist(),
    }
    tensor = {
        'status': 'DoubletGroupTensorRatiosVerified',
        'full_lagrangian_coupling_identification_verified': False,
        'mappings': [{
            'model': m, 'matchete_coupling': mat, 'rgbeta_coupling': rg,
            'rgbeta_over_matchete_tensor_factor': 1.0,
            'full_lagrangian_coupling_identification_verified': False,
        } for m in ('B_m1','C_m1') for mat,rg in ACTIVE.items()],
    }
    data = dict(zip(ACTIVE.values(), (complex(2e-12, 1e-12), complex(-1e-12,2e-12))))
    stage13 = {
        'model': model, 'status': 'CoupledUVOneLoopRGBetaOnly',
        'solver_success': True,
        'full_matchete_rgbeta_coupling_mapping_verified': False,
        'includes_mass_threshold_matching': False,
        'uv_scale_gev': 20000., 'endpoint_gev': 3500.,
        'active_couplings': [{
            'coupling': rg, 'uv_initial_rgbeta': _z(0),
            'uv_endpoint_rgbeta': _z(z), 'change_coupled': _z(z),
        } for rg,z in data.items()],
    }
    convention = {
        'status': 'SourceLagrangianConventionDerived',
        'source_checks_passed': True,
        'full_threshold_matching_verified': False,
        'same_field_basis_assumed': True,
        'results': [{
            'model': m,
            'matchete_coupling': mat,
            'rgbeta_coupling': rg,
            'matchete_in_terms_of_rgbeta': 'negative',
            'relative_sign': -1,
            'requires_complex_conjugation': False,
        } for m in ('B_m1','C_m1') for mat,rg in ACTIVE.items()],
    }
    kernels = tuple(ResponseKernel(
        coupling=mat, conjugated=True,
        yukawa={'lambdaT3ExtraQ020020I1': 'y2',
                'lambdaT3ExtraQ020200I1': 'y1'}[mat],
        expression=sp.Symbol('hbar') / sp.Symbol('MF'),
        raw_expression='synthetic hbar / MF * conjugate(y)**2'
    ) for mat in ACTIVE)
    reference = {'model': model, 'status': 'ExperimentalHardBoundaryC5',
                 'units': 'GeV^-1', 'delta_c5': {'frobenius_norm': 1e-17}}
    return params, tensor, stage13, convention, kernels, reference


def _run_fixture(model='B_m1'):
    params, tensor, s13, convention, kernels, ref = _fixture(model)
    s14 = conditional_uv_to_c5(model, params, s13, tensor, kernels, ref)
    result = source_mapped_uv_to_c5(model, params, s13, tensor, convention,
                                     kernels, ref, stage14=s14)
    return result, (params, tensor, s13, convention, kernels, ref, s14)


def test_matchete_real_source_contract():
    source=(ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl').read_text()
    _check_matchete_source(source)


def test_rgbeta_source_contract_on_real_repository_if_available():
    source=ROOT/'RGE/running/rgbeta/T3RGBetaModel.wl'
    if not source.exists():
        pytest.skip('Production RGBeta file available only in user repository')
    _check_rgbeta_source(source.read_text())


def test_matchete_changed_sign_rejected():
    source=(ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl').read_text()
    with pytest.raises(ValueError, match='Matchete source contract'):
        _check_matchete_source(source.replace('PlusHc[contribution]','PlusHc[-contribution]'))


def test_matchete_changed_coupling_coefficient_rejected():
    source=(ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl').read_text()
    with pytest.raises(ValueError,match='Matchete source contract'):
        _check_matchete_source(source.replace('contribution = coupling[] product',
                                              'contribution = coupling[] product / 2'))


def test_rgbeta_field_orientation_and_vertex_factor_rejected():
    source='''AddQuartic[lambdaHHdagS2S2, {Bar @ H, Bar @ H, S2, S2},
        GroupInvariant -> (delS2[SU2L,a,#1,#2]*delS2[SU2L,a,#3,#4]&),
        SelfConjugate->False];
      AddQuartic[lambdaHHdagS1barS1bar,
        {Bar @ H, Bar @ H, Bar @ S1, Bar @ S1},
        GroupInvariant -> (delS2[SU2L,a,#1,#2]*delS2[SU2L,a,#3,#4]&),
        SelfConjugate->False];'''
    _check_rgbeta_source(source)
    with pytest.raises(ValueError,match='RGBeta source contract'):
        _check_rgbeta_source(source.replace('Bar @ S1, Bar @ S1','S1, S1'))
    with pytest.raises(ValueError,match='RGBeta source contract'):
        _check_rgbeta_source(source.replace('delS2[SU2L,a,#3,#4]&',
                                            '(1/2)*delS2[SU2L,a,#3,#4]&'))


def test_convention_report_derives_negative_without_conjugation(tmp_path):
    msrc=(ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl')
    rgb='''AddQuartic[lambdaHHdagS2S2, {Bar @ H, Bar @ H, S2, S2},
        GroupInvariant -> (delS2[SU2L,a,#1,#2]*delS2[SU2L,a,#3,#4]&),
        SelfConjugate->False];
      AddQuartic[lambdaHHdagS1barS1bar,
        {Bar @ H, Bar @ H, Bar @ S1, Bar @ S1},
        GroupInvariant -> (delS2[SU2L,a,#1,#2]*delS2[SU2L,a,#3,#4]&),
        SelfConjugate->False];'''
    rgb_path=tmp_path/'RGBeta.wl';rgb_path.write_text(rgb)
    tensor=_fixture()[1]
    tensor_path=tmp_path/'tensor.json';tensor_path.write_text(json.dumps(tensor))
    output=derive_lagrangian_convention(msrc,rgb_path,tensor_path)
    assert output['status']=='SourceLagrangianConventionDerived'
    assert len(output['results'])==4
    assert {r['relative_sign'] for r in output['results']}=={-1}
    assert all(r['requires_complex_conjugation'] is False for r in output['results'])
    assert output['full_threshold_matching_verified'] is False


def test_tensor_report_fail_closed():
    tensor=_fixture()[1]
    _check_tensor_report(tensor)
    tensor['mappings'][1]['rgbeta_over_matchete_tensor_factor']=-1
    with pytest.raises(ValueError,match='ratio'):
        _check_tensor_report(tensor)


@pytest.mark.parametrize('model',['B_m1','C_m1'])
def test_mapped_response_reproduces_negative_negative_case(model):
    result, _=_run_fixture(model)
    assert result['status']=='SourceMappedConditionalUVInducedC5'
    assert result['stage14_negative_case_reproduced']
    assert not result['threshold_matching']
    assert not result['physical_neutrino_spectrum']
    singular=result['delta_mnu_ev']['takagi_singular_values_ev']
    assert np.isclose(np.linalg.norm(singular), result['delta_mnu_ev']['frobenius_norm'])
    for entry in result['mapped_matchete_extra_couplings'].values():
        assert abs(complex(entry['real'],entry['imag']))>0


def test_mapping_is_negative_not_conjugate():
    result, (_, _, s13, _, _, _, _)=_run_fixture()
    expected={row['coupling']:complex(row['change_coupled']['real'],row['change_coupled']['imag'])
              for row in s13['active_couplings']}
    for mat,rg in ACTIVE.items():
        actual=result['mapped_matchete_extra_couplings'][mat]
        assert complex(actual['real'],actual['imag']) == -expected[rg]
        assert complex(actual['real'],actual['imag']) != -expected[rg].conjugate()


def test_reject_changed_source_convention_or_bad_stage14():
    params,tensor,s13,convention,kernels,reference=_fixture()
    convention['results'][0]['requires_complex_conjugation']=True
    with pytest.raises(ValueError,match='conjugation'):
        source_mapped_uv_to_c5('B_m1',params,s13,tensor,convention,kernels,reference)
    convention['results'][0]['requires_complex_conjugation']=False
    bad=conditional_uv_to_c5('B_m1',params,s13,tensor,kernels,reference)
    for case in bad['mapping_cases']:
        if all(v=='negative' for v in case['mapping_choices'].values()):
            case['extra_only_mnu_norm_ev']*=3
    with pytest.raises(ValueError,match='Stage 14'):
        source_mapped_uv_to_c5('B_m1',params,s13,tensor,convention,kernels,reference,stage14=bad)


def test_cli_outputs_both_models_with_synthetic_input(tmp_path):
    # Production sources from Stage 15 package; minimal exact RGBeta source in isolated test.
    root=tmp_path/'run'
    paths={p:root/p for p in ('stage13','stage14','matchete','trial/B_m1/general_c5_algebraic',
                             'trial/C_m1/general_c5_algebraic')}
    for d in paths.values():d.mkdir(parents=True)
    pars,tensor,_,_,_,_=_fixture()
    (root/'input.json').write_text(json.dumps(pars))
    (root/'stage13/doublet_tensor_map.json').write_text(json.dumps(tensor))
    rgb='''AddQuartic[lambdaHHdagS2S2,{Bar@H,Bar@H,S2,S2},
        GroupInvariant->(delS2[SU2L,a,#1,#2]*delS2[SU2L,a,#3,#4]&),SelfConjugate->False];
       AddQuartic[lambdaHHdagS1barS1bar,{Bar@H,Bar@H,Bar@S1,Bar@S1},
        GroupInvariant->(delS2[SU2L,a,#1,#2]*delS2[SU2L,a,#3,#4]&),SelfConjugate->False];'''
    (root/'RGBeta.wl').write_text(rgb)
    for model in ('B_m1','C_m1'):
        _,_,s13,_,kernels,ref=_fixture(model)
        (root/f'stage13/{model}_coupled_uv.json').write_text(json.dumps(s13))
        (root/f'matchete/{model}_matchete_hard_c5.json').write_text(json.dumps(ref))
        s14=conditional_uv_to_c5(model,pars,s13,tensor,kernels,ref)
        (root/f'stage14/{model}_conditional_induced_c5.json').write_text(json.dumps(s14))
        atom_data={
            'model':model,
            'response_reconstruction_residual_proven_zero':True,
            'response_at_zero_proven_zero':True,
            'polynomial_in_extra_atoms':True,
            'atoms':[{'coupling_atom':f'Bar[Coupling[{mat},{{}},0]]',
                      'linear_coefficient_inputform':f'hbar*Bar[Coupling[{k.yukawa},{{}},0]]^2/Coupling[MF,{{}},0]',
                      'coefficient_formally_zero':False,
                      'depends_on_other_extra_atoms':False} for mat,k in zip(ACTIVE,kernels)],
        }
        (root/f'trial/{model}/general_c5_algebraic/delta_c5_atoms.json').write_text(json.dumps(atom_data))
    command=[sys.executable,str(ROOT/'scripts/RunGeneralT3Stage15.py'),
             '--input',str(root/'input.json'),
             '--rgbeta-source',str(root/'RGBeta.wl'),
             '--matchete-source',str(ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl'),
             '--trial-root',str(root/'trial'), '--stage13-root',str(root/'stage13'),
             '--stage14-root',str(root/'stage14'), '--matchete-root',str(root/'matchete'),
             '--output',str(root/'output')]
    run=subprocess.run(command,capture_output=True,text=True,cwd=ROOT)
    assert run.returncode==0,run.stdout+'\n'+run.stderr
    assert (root/'output/B_m1_source_mapped_extra_c5.json').exists()
    assert (root/'output/C_m1_source_mapped_extra_c5.json').exists()
    assert json.loads((root/'output/stage15_summary.json').read_text())['status']=='Success'
    assert (root/'output/stage15_summary.csv').exists()


def test_powershell_runner_isolation():
    text=(ROOT/'scripts/RunGeneralT3Stage15.ps1').read_text()
    assert 'RunGeneralT3Stage15.py' in text
    assert 'ConventionOnly' in text
    assert 'RGE/running/' not in text
