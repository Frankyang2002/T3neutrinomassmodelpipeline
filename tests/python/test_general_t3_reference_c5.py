"""Reference-Yukawa transfer: frozen provenance, 3 nondegenerate F masses, no EFT1 claims."""
from __future__ import annotations
import copy
from hashlib import sha256
import json
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts import RunGeneralT3NumericalC5 as num
from scripts import RunGeneralT3ReferenceC5 as ref


def test_manifest_extract_and_mismatch(tmp_path):
    slug=ref.REFERENCE_ID['B_-1'][0]
    folder=tmp_path/'Numerical/benchmarks/fitted_reference';folder.mkdir(parents=True)
    z=lambda value:{'real':value,'imag':0.}
    y=lambda v:{'real':[[v]*3 for _ in range(3)],'imag':[[0]*3 for _ in range(3)]}
    mat={'real':[[100000.,0.,0.],[0.,110000.,0.],[0.,0.,120000.]],
         'imag':[[0]*3 for _ in range(3)]}
    frozen={
        'representation':{'d_s1':2,'d_s2':2,'d_f':1,'alpha':-1,'shared_scalar':False},
        'numerical':{'vev_gev':246.22},
        'benchmark_search':{'enabled':False},'sensitivity':{'enabled':False},
        'base_state':{'ordinary':{'y1':y(.0001),'y2':y(.0002),'MF':mat,
            'mS1Sq':1.3e6,'mS2Sq':1.e6,'lambdaT3':z(0.1),
            'lambdaHHdagS2S2':z(0.),'lambdaHHdagS1barS1bar':z(0.)}}
    }
    fpath=folder/(slug+'.json');fpath.write_text(json.dumps(frozen),encoding='utf-8')
    h=sha256(fpath.read_bytes()).hexdigest()
    entries=[{'model':slug,'sha256':h,'config_path':'Numerical/benchmarks/fitted_reference/'+slug+'.json'}]
    entries.extend({'model':f'other_{i}'} for i in range(15))
    manifest={'models':entries}
    orig={'model':'B_-1','verified_eight_real_mode_hessian':True,
          'original_mass_squared_A_gev2':{'real':[[1.4e6,1e4],[1e4,1.07e6]],'imag':[[0,0],[0,0]]}}
    ppath=tmp_path/'parameters.json';ppath.write_text(json.dumps(orig),encoding='utf-8')
    config=ref.make_transfer('B_-1',manifest,fpath,ppath)
    assert config['fermion_masses_gev']==[100000.,110000.,120000.]
    assert len(config['y1'])==3 and len(config['y1'][0])==3
    assert config['reference_provenance']['mass_squared_difference_gev2']==[100000.,70000.]
    assert config['mass_basis_compatible_with_frozen_reference'] is False
    assert config['reference_transfer_only'] is True
    assert config['quartics']['lambdaT3']['real']==.1
    with pytest.raises(ValueError,match='hash changed'):
        bad=copy.deepcopy(frozen);bad['base_state']['ordinary']['mS2Sq']=9.e5
        fpath.write_text(json.dumps(bad),encoding='utf-8')
        ref.make_transfer('B_-1',manifest,fpath,ppath)
    fpath.write_text(json.dumps(frozen),encoding='utf-8')
    bad=copy.deepcopy(frozen);bad['base_state']['ordinary']['MF']['real'][0][1]=4
    fpath.write_text(json.dumps(bad),encoding='utf-8')
    manifest['models'][0]['sha256']=sha256(fpath.read_bytes()).hexdigest()
    with pytest.raises(ValueError,match='diagonal'):
        ref.make_transfer('B_-1',manifest,fpath,ppath)


def test_exploratory_original_masses_and_labels(tmp_path):
    source=ROOT/'output/general_t3_mass_eigen_matchete/one_loop_probe'
    params=ROOT/'output/general_t3_mass_eigen_matchete'
    if not (source/'B_m1_final_one_loop_probe.json').is_file():
        pytest.skip('Actual Matchete coefficient exports not installed')
    config=num.read_json(ROOT/'configs/general_t3/c5_illustrative.json')
    base,_=num.numerical_model('B_-1',config,source,params,tmp_path)
    assert base['fermion_masses_gev']==[100000]
    assert base['max_abs_mnu_eV']==pytest.approx(61.2038969873076,rel=1e-9)
    assert base['largest_takagi_mass_eV']==pytest.approx(148.42297259166764,rel=1e-9)
    assert base['max_abs_mnu_eV']<base['largest_takagi_mass_eV']
    wrong=copy.deepcopy(config)
    wrong['fermion_masses_gev']=[100000.,110000.]
    with pytest.raises(ValueError,match='count does not match'):
        num.numerical_model('B_-1',wrong,source,params,tmp_path)


def test_multiple_fermion_mass_flavours_are_independently_evaluated(tmp_path):
    source=ROOT/'output/general_t3_mass_eigen_matchete/one_loop_probe'
    params=ROOT/'output/general_t3_mass_eigen_matchete'
    if not (source/'B_m1_final_one_loop_probe.json').is_file():
        pytest.skip('Actual Matchete coefficient exports not installed')
    base=num.read_json(ROOT/'configs/general_t3/c5_illustrative.json')
    config=copy.deepcopy(base)
    for key in ('y1','y2'):
        for row in config[key]:
            val=row[0]
            row[:]=[copy.deepcopy(val) for _ in range(3)]
    config['fermion_masses_gev']=[100000.,110000.,120000.]
    combined,_=num.numerical_model('B_-1',config,source,params,tmp_path)
    assert combined['fermion_masses_gev']==[100000.,110000.,120000.]
    expected=0j
    for mass in (100000.,110000.,120000.):
        part=copy.deepcopy(base);part['fermion_mass_gev']=mass
        result,_=num.numerical_model('B_-1',part,source,params,tmp_path)
        z=result['symmetrized_C5_gev_inv'][0][0]
        expected+=complex(z['real'],z['imag'])
    actual=combined['symmetrized_C5_gev_inv'][0][0]
    assert complex(actual['real'],actual['imag'])==pytest.approx(expected,rel=1e-10)
    assert combined['symmetry_max_abs_C5_gev_inv']==0
    assert combined['physical_c5_certified'] is False
