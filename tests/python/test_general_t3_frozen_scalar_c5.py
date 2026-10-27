"""Frozen scalar diagonal scenario, explicit mu2 and native Matchete C5 evaluation."""
from __future__ import annotations

import copy
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
from scripts import RunGeneralT3FrozenScalarC5 as aligned
from scripts import RunGeneralT3NumericalC5 as num
from scripts import RunGeneralT3ReferenceC5 as ref


def fixture_files(tmp_path,model='B_-1'):
    slug,rep=ref.REFERENCE_ID[model]
    folder=tmp_path/'Numerical/benchmarks/fitted_reference';folder.mkdir(parents=True)
    z=lambda x:{'real':x,'imag':0.}
    y=lambda x:{'real':[[x]*3 for _ in range(3)],'imag':[[0]*3 for _ in range(3)]}
    mass={'real':[[100000.,0.,0.],[0.,110000.,0.],[0.,0.,120000.]],
          'imag':[[0]*3 for _ in range(3)]}
    frozen={'representation':dict(zip(('d_s1','d_s2','d_f','alpha'),rep),shared_scalar=False),
            'numerical':{'vev_gev':246.22},
            'benchmark_search':{'enabled':False},'sensitivity':{'enabled':False},
            'base_state':{'ordinary':{'y1':y(1.e-4),'y2':y(1.2e-4),'MF':mass,
                'mS1Sq':1300000.,'mS2Sq':1000000.,'lambdaT3':z(.1),
                'lambdaHHdagS2S2':z(0.),'lambdaHHdagS1barS1bar':z(0.)}}}
    p=folder/(slug+'.json');p.write_text(json.dumps(frozen))
    manifest={'models':[{'model':slug,'sha256':sha256(p.read_bytes()).hexdigest(),
             'config_path':'Numerical/benchmarks/fitted_reference/'+slug+'.json'}]+[{'model':f'other{i}'} for i in range(15)]}
    parm={'model':model,'verified_eight_real_mode_hessian':True,
          'original_mass_squared_A_gev2':{'real':[[1400000.,11000.],[11000.,1080000.]],'imag':[[0,0],[0,0]]}}
    pp=tmp_path/'params.json';pp.write_text(json.dumps(parm))
    return manifest,p,pp


def test_frozen_diagonals_zero_mu2_not_old_hybrid(tmp_path):
    manifest,path,params=fixture_files(tmp_path)
    conf,basis=aligned.build_scenario('B_-1',manifest,path,params,0j)
    assert conf['format']==num.ALIGNED_FORMAT
    assert conf['reference_transfer_only'] is False
    assert conf['quadratic_mixing_is_explicit_scenario'] is True
    assert conf['physical_input_certified'] is False
    assert conf['reference_provenance']['quadratic_mixing_source'].startswith('zero_by_explicit')
    assert np.array_equal(basis.eigenvalues_sq,np.array([1000000.,1300000.]))
    assert np.allclose(basis.U,[[0.,1.],[1.,0.]])
    assert np.array_equal(np.diag(basis.A).real,[1300000.,1000000.])
    assert not np.array_equal(np.diag(basis.A).real,[1400000.,1080000.])


def test_nonzero_complex_mu2_is_diagonalized_and_certified(tmp_path):
    manifest,path,params=fixture_files(tmp_path,'C_-1')
    mu2=10000.+2000.j
    config,basis=aligned.build_scenario('C_-1',manifest,path,params,mu2,mixing_explicit_on_cli=True)
    assert np.allclose(basis.A,[[1.3e6,mu2],[mu2.conjugate(),1e6]])
    assert np.allclose(basis.U.conj().T@basis.A@basis.U,np.diag(basis.eigenvalues_sq))
    assert config['reference_provenance']['quadratic_mixing_mu2_gev2']=={'real':10000.,'imag':2000.}
    assert config['reference_provenance']['quadratic_mixing_source'].startswith('user_supplied')


def test_tachyonic_mixing_refused(tmp_path):
    manifest,path,params=fixture_files(tmp_path)
    with pytest.raises(ValueError,match='non-positive|tachyonic'):
        aligned.build_scenario('B_-1',manifest,path,params,2e6)


def test_mismatched_reference_hash_refused(tmp_path):
    manifest,path,params=fixture_files(tmp_path)
    data=json.loads(path.read_text());data['base_state']['ordinary']['mS1Sq']=1200000.
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='hash changed'):
        aligned.build_scenario('B_-1',manifest,path,params,0j)


def test_actual_matchete_expression_with_scalar_aligned_override(tmp_path):
    source=ROOT/'output/general_t3_mass_eigen_matchete/one_loop_probe'
    params=ROOT/'output/general_t3_mass_eigen_matchete'
    if not (source/'B_m1_final_one_loop_probe.json').is_file():
        pytest.skip('actual symbolic one-loop outputs unavailable')
    illustrative=num.read_json(ROOT/'configs/general_t3/c5_illustrative.json')
    aligned_config=copy.deepcopy(illustrative)
    aligned_config.update({'format':num.ALIGNED_FORMAT,'frozen_reference_scalar_diagonals':True,
        'quadratic_mixing_is_explicit_scenario':True,
        'mass_basis_compatible_with_frozen_reference':True,'reference_transfer_only':False,
        'reference_provenance':{
            'frozen_scalar_diagonal_masses_squared_gev2':[1.3e6,1e6],
            'quadratic_mixing_mu2_gev2':{'real':0.,'imag':0.}}})
    b=MixedDoubletMassBasis.build('B_-1',1.3e6,1e6,0j)
    result,_=num.numerical_model('B_-1',aligned_config,source,params,tmp_path,scalar_basis_override=b)
    assert result['status']==num.ALIGNED_STATUS
    assert result['frozen_scalar_aligned_one_scale_scenario'] is True
    assert result['source_mass_parameters_used_as_numeric_scalar_input'] is False
    assert result['scalar_masses_gev']==pytest.approx([1000,np.sqrt(1.3e6)])
    assert result['largest_takagi_mass_eV']>0
    assert result['physical_c5_certified'] is False
    assert result['sequential_EFT1_performed'] is False
    result2,_=num.numerical_model('B_-1',illustrative,source,params,tmp_path)
    assert result2['scalar_masses_gev']!=result['scalar_masses_gev']
    assert result2['frozen_scalar_aligned_one_scale_scenario'] is False


def test_override_requires_aligned_format_and_matching_model(tmp_path):
    source=ROOT/'output/general_t3_mass_eigen_matchete/one_loop_probe'
    params=ROOT/'output/general_t3_mass_eigen_matchete'
    if not (source/'B_m1_final_one_loop_probe.json').is_file():pytest.skip('actual coefficient unavailable')
    template=num.read_json(ROOT/'configs/general_t3/c5_illustrative.json')
    b=MixedDoubletMassBasis.build('B_-1',1.3e6,1e6,0j)
    with pytest.raises(ValueError,match='requires explicit validated'):
        num.numerical_model('B_-1',template,source,params,tmp_path,scalar_basis_override=b)
    template.update({'format':num.ALIGNED_FORMAT,'frozen_reference_scalar_diagonals':True,
                     'quadratic_mixing_is_explicit_scenario':True,
                     'mass_basis_compatible_with_frozen_reference':True,'reference_transfer_only':False,
                     'reference_provenance':{'frozen_scalar_diagonal_masses_squared_gev2':[1.3e6,1e6],
                        'quadratic_mixing_mu2_gev2':{'real':0.,'imag':0.}}})
    c=MixedDoubletMassBasis.build('C_-1',1.3e6,1e6,0j)
    with pytest.raises(ValueError,match='selected model'):
        num.numerical_model('B_-1',template,source,params,tmp_path,scalar_basis_override=c)
    template['reference_provenance']['quadratic_mixing_mu2_gev2']={'real':5000.,'imag':0.}
    with pytest.raises(ValueError,match='quadratic mixing'):
        num.numerical_model('B_-1',template,source,params,tmp_path,scalar_basis_override=b)


def test_new_frozen_quadratic_mixing_value_must_not_be_overridden(tmp_path):
    manifest,path,params=fixture_files(tmp_path)
    frozen=json.loads(path.read_text())
    frozen['base_state']['ordinary']['mu2']={'real':123.,'imag':0.}
    path.write_text(json.dumps(frozen))
    manifest['models'][0]['sha256']=sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError,match='frozen reference now supplies'):
        aligned.build_scenario('B_-1',manifest,path,params,0j)
