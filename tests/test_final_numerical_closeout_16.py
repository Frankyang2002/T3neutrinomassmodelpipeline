"""Self-contained closeout tests; no production data or native installations needed."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.lo_patterns_16.final_numerical_closeout import (
    BASE,B5,B6,STUDY,DEST,MODEL_NAMES,TEXTURES,SCENARIOS,SCAN_GRID,
    CloseoutBlocked,close,compute,write,
)


def jsonfile(p:Path,obj:dict) -> str:
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8',newline='\n')
    return hashlib.sha256(p.read_bytes()).hexdigest()


def csvfile(p:Path,rows:list[dict]) -> None:
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='',encoding='utf-8') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


@pytest.fixture
def fake_project(tmp_path):
    root=tmp_path
    catalogue=[]; controls=[]; scenarios=[]; scans=[]; b5models=[]
    modelmanifest={};model6={}
    for index,name in enumerate(MODEL_NAMES):
        family,alpha=name.split('_');alpha=int(alpha)
        ty='Majorana' if (family in ('B','C') and alpha==-1) else 'Dirac'
        masses=(0.001+index*1e-5, 0.008+index*1e-5,0.05+index*1e-5)
        frozen={'model':name,'family':family,'alpha':alpha,
                'fermion_type':ty,'status':'conditional_uv_evaluated',
                'physical_normalisation_certified':False,
                'candidate_topology_charge_channel_count':2 if family=='E' else 1,
                'benchmark_sha256':f'{index+1:064x}',
                'component_input_sha256':f'{index+20:064x}',
                'source_reference_relative_error':2e-11,
                'm1_ev':masses[0],'m2_ev':masses[1],'m3_ev':masses[2]}
        catalogue.append(frozen)
        modelmanifest[name]={'benchmark_sha256':frozen['benchmark_sha256'],
                             'candidate_weyl_input_sha256':frozen['component_input_sha256'],
                             'study_model_sha256':f'{index+100:064x}'}
        model6[name]={'benchmark_sha256':frozen['benchmark_sha256'],
                      'input_sha256':frozen['component_input_sha256'],
                      'study_model_sha256':f'{index+100:064x}'}
        ctr={}
        for j,tex in enumerate(TEXTURES):
            cs=(0.004+j*0.001,0.012+j*0.001,0.055+j*0.001)
            ctr[tex]={'m1_ev':cs[0],'m2_ev':cs[1],'m3_ev':cs[2],
                      'm1_over_m3':cs[0]/cs[2],
                      'm1_geometry_A':0.3,'m1_coherence_C':0.4}
        b5models.append({'model':name,'family':family,'alpha':alpha,
                         'fermion_type':ty,'physical_normalisation_certified':False,
                         'conditional_only':True,
                         'observed':{'m1_ev':masses[0],'m2_ev':masses[1],'m3_ev':masses[2],
                                     'm1_over_m3':masses[0]/masses[2],
                                     'm1_geometry_A':0.2,'m1_coherence_C':0.4,
                                     'm1_mechanism':'mode_geometry_dominant'},
                         'common_mixing_controlled':ctr})
        for mode,ts in SCAN_GRID.items():
            for t in ts:
                fac=0 if t==0 else abs(t) if mode!='heavy_mass' else (1.0 if t==1 else t)
                scans.append({'model':name,'mode':mode,'t':t,
                              'status':'conditional_uv_scanned',
                              'm1_ev':fac*masses[0],
                              'm2_ev':fac*masses[1],
                              'm3_ev':fac*masses[2]})
        for sc in SCENARIOS:
            scenarios.append({'scenario':sc,'model':name,'reference_model':'A_-2',
                              'coefficient_residual':0.0,
                              'max_group_fractional_spread':0.0})
    counts={'catalogue':16,'conditional_evaluated':16,'blocked':0,
            'physical_source_certified':0,'scan_passes':288,'scan_blocked':0}
    rev={'local_head':'0'*40}
    ss={'format':'t3-conditional-pattern-study-all16-v1','counts':counts,
        'provenance':{'revision':rev},'models':catalogue,
        'physical_normalisation_certified':False}
    p=root/STUDY
    studysha=jsonfile(p/'all16_summary.json',ss)
    jsonfile(p/'HANDOFF.json',{'checks':counts,'physical_source_certified':False,
                               'safe_to_integrate_for_physical_mass':False})
    csvfile(p/'all16_summary.csv',[{**x,'physical_normalisation_certified':'False'} for x in catalogue])
    csvfile(p/'scans.csv',scans)
    b5={'format':'t3-chatA-batch05-analytical-yukawa-texture-v1',
        'physical_normalisation_certified':False,'source_conditional_only':True,
        'safe_for_physical_source_certification':False,
        'model_count':16,'reference_texture_models':list(TEXTURES),
        'common_scalar_mixing_abs':0.1,'models':b5models,
        'source_manifest_sha256':{'study':studysha}}
    p5=root/B5
    b5hash=jsonfile(p5/'batch05_synthesis.json',b5)
    manifest={'format':'t3-chatA-batch05-source-manifest-v1',
              'study_sha256':studysha,'physical_normalisation_certified':False,
              'models':modelmanifest}
    m5hash=jsonfile(p5/'SOURCE_MANIFEST.json',manifest)
    sourcehash={'study':studysha,'batch05':b5hash,'batch05_manifest':m5hash}
    b6={'format':'t3-chatA-batch06-heavy-kernel-factorization-v1',
        'source_conditional_only':True,'physical_normalisation_certified':False,
        'safe_for_physical_source_certification':False,'model_count':16,
        'models':list(MODEL_NAMES),'reference_model':'A_-2',
        'scenarios':list(SCENARIOS),'source_sha256':sourcehash,
        'coefficient_comparisons':scenarios,
        'tensor_reconstruction_checks':[
            {'scenario':x['scenario'],'model':x['model'],
             'exact_tensor_reconstruction_error':1e-15} for x in scenarios],
        'model_summary':[{'model':n,'baseline_coefficient_residual':0.0}
                         for n in MODEL_NAMES]}
    p6=root/B6
    jsonfile(p6/'batch06_factorization.json',b6)
    jsonfile(p6/'SOURCE_MANIFEST.json',{
        'format':'t3-chatA-batch06-manifest-v1',
        'physical_normalisation_certified':False,
        'source_sha256':sourcehash,'models':model6})
    csvfile(p6/'batch06_coefficient_tests.csv',scenarios)
    return root


def test_final_closeout_success(fake_project):
    report,out=compute(fake_project)
    assert len(report['results'])==16
    assert report['scan_count']==288
    assert len(report['controls'])==64
    assert len(report['scenarios'])==5
    assert report['zero_checks']==48 and report['baseline_checks']==64
    assert len(out)==4 and 'FINAL_NUMERICAL_CLOSEOUT.md' in out
    assert '0/16' in out['FINAL_NUMERICAL_CLOSEOUT.md']
    assert 'not a physical' in out['FINAL_NUMERICAL_CLOSEOUT.md'].lower()


def test_write_new_dir_and_refuse_overwrite(fake_project):
    report,files=write(fake_project,DEST)
    assert len(files)==4 and (fake_project/DEST/'FINAL_NUMERICAL_CLOSEOUT.md').is_file()
    with pytest.raises(CloseoutBlocked,match='Refusing to overwrite'):
        write(fake_project,DEST)


def test_reject_output_outside_study(fake_project):
    with pytest.raises(CloseoutBlocked,match='Output must be'):
        write(fake_project,Path('outside'))


def test_reject_modified_frozen_study_bytes(fake_project):
    p=fake_project/STUDY/'all16_summary.json'
    obj=json.loads(p.read_text());obj['models'][0]['m1_ev']*=0.9
    jsonfile(p,obj)
    with pytest.raises(CloseoutBlocked,match='SHA256 differs'):
        compute(fake_project)


def test_reject_missing_scans(fake_project):
    p=fake_project/STUDY/'scans.csv'
    rows=list(csv.DictReader(p.open(encoding='utf-8')))
    csvfile(p,rows[:-1])
    with pytest.raises(CloseoutBlocked,match='Scan grid not complete'):
        compute(fake_project)


def test_reject_zero_yukawa_nonzero_mass(fake_project):
    p=fake_project/STUDY/'scans.csv'
    rows=list(csv.DictReader(p.open(encoding='utf-8')))
    rows[1]['m1_ev']='0.001'
    csvfile(p,rows)
    with pytest.raises(CloseoutBlocked,match='Zero Yukawa'):
        compute(fake_project)


def test_reject_b6_coefficient_changed(fake_project):
    p=fake_project/B6/'batch06_coefficient_tests.csv'
    rows=list(csv.DictReader(p.open(encoding='utf-8')))
    rows[0]['coefficient_residual']='0.42'
    csvfile(p,rows)
    with pytest.raises(CloseoutBlocked,match='JSON/CSV inconsistent'):
        compute(fake_project)


def test_reject_physical_certification_invented(fake_project):
    p=fake_project/B6/'SOURCE_MANIFEST.json'
    obj=json.loads(p.read_text());obj['physical_normalisation_certified']=True
    jsonfile(p,obj)
    with pytest.raises(CloseoutBlocked,match='physical-normalisation flag'):
        compute(fake_project)


def test_reject_wrong_model_yukawa_manifest(fake_project):
    p=fake_project/B6/'SOURCE_MANIFEST.json'
    obj=json.loads(p.read_text());obj['models']['C_+1']['input_sha256']='0'*64
    jsonfile(p,obj)
    with pytest.raises(CloseoutBlocked,match='benchmark/input lineage mismatch'):
        compute(fake_project)


def test_accept_pure_crlf_checkout(fake_project):
    p=fake_project/STUDY/'all16_summary.json'
    p.write_bytes(p.read_bytes().replace(b'\n',b'\r\n'))
    report,_=compute(fake_project)
    assert report['digests']['hash_verification_modes'][
        'provisional_01/all16_summary.json']=='CRLF-to-LF-only'


def test_reject_non_line_ending_mutation(fake_project):
    p=fake_project/STUDY/'all16_summary.json'
    p.write_bytes(p.read_bytes().replace(b'\n',b'\r\n')+b' ')
    with pytest.raises(CloseoutBlocked,match='SHA256 differs'):
        compute(fake_project)


def test_reject_controlled_masses_not_sorted(fake_project):
    p=fake_project/B5/'batch05_synthesis.json'
    obj=json.loads(p.read_text());obj['models'][0]['common_mixing_controlled']['C_+1']['m1_ev']=1
    # Rebuild chained lineage to isolate the *physics data check* from hash guards.
    b5hash=jsonfile(p,obj)
    p6=fake_project/B6/'batch06_factorization.json'
    b6=json.loads(p6.read_text());b6['source_sha256']['batch05']=b5hash
    jsonfile(p6,b6)
    man=fake_project/B6/'SOURCE_MANIFEST.json';m=json.loads(man.read_text())
    m['source_sha256']['batch05']=b5hash
    jsonfile(man,m)
    with pytest.raises(CloseoutBlocked,match='controlled masses not ordered'):
        compute(fake_project)


@pytest.mark.parametrize('x,y',[(0,0),(1.0,1.0+1e-9),(0.005,0.005+1e-12)])
def test_close_numeric_contract(x,y):
    assert close(x,y)
