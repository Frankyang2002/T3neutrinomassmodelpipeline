"""Batch 2 synthetic schema and negative controls; NOT 16 physical benchmarks."""
from __future__ import annotations

import csv
import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
from Numerical.lo_patterns_16.result_analysis import analyse,AnalysisError,MODES,T_VALUES,HEAVY_VALUES


def encode(a):
    z=np.asarray(a,dtype=complex)
    return {'real':z.real.tolist(),'imag':z.imag.tolist()}


def put_json(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data),encoding='utf-8')


def h(path):
    return sha256(path.read_bytes()).hexdigest()


def catalogue_fixture(tmp_path):
    root=tmp_path/'project';root.mkdir()
    study=root/'output/general_observable/lo_pattern_study16/provisional_01'
    study.mkdir(parents=True)
    models=model_dictionary()
    names=[m['model'] for m in models]
    family={'A':(1,3,2),'B':(2,2,1),'C':(2,2,3),
            'D':(3,1,2),'E':(3,3,2)}
    summary=[];scans=[];weight_rows={}
    for idx,item in enumerate(models):
        model=item['model'];fam=model[0];alpha=int(model[2:]);a,b,f=family[fam]
        sign='m' if alpha<0 else 'p'
        key=f'T3_dS1_{a}_dS2_{b}_dF_{f}_alpha_{sign}{abs(alpha)}'
        benchmark=root/'Numerical/benchmarks/fitted_reference'/f'{key}.json'
        uv=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'/f'{model}.json'
        put_json(benchmark,{'model':model,'kind':'SYNTHETIC'})
        put_json(uv,{'model':model,'numerically_certified':False})
        factor=1+0.02*idx
        masses=np.array([.001,.009,.045])*factor
        matrix=np.diag(masses/1e9).astype(complex)
        blank=np.zeros((3,3),complex)
        if fam=='E':
            charges=[alpha/2,alpha/2+1]
            qmatrix={f'{charges[0]:g}|{charges[0]:g}':matrix*.4,
                     f'{charges[1]:g}|{charges[1]:g}':matrix*.6}
        else:
            q=alpha/2 if fam=='A' else (alpha+1)/2 if fam in 'BC' else (alpha+2)/2
            charges=[q]
            qmatrix={f'{q:g}|{q:g}':matrix}
        ch=[{'charge_numeric':q,'complex_weight':str(1.),'nonzero':True} for q in charges]
        bsha=h(benchmark);usha=h(uv)
        summary.append({'model':model,'family':fam,'alpha':alpha,
                        'm1_ev':masses[0],'m2_ev':masses[1],'m3_ev':masses[2],
                        'candidate_topology_charge_channel_count':len(charges),
                        'cancellation_factor':1.5,'source_reference_relative_error':1e-12,
                        'benchmark_sha256':bsha,'component_input_sha256':usha,
                        'physical_normalisation_certified':False,'status':'conditional_uv_evaluated'})
        detail={'model':model,'status':'conditional_uv_evaluated',
                'physical_normalisation_certified':False,
                'source':{'benchmark':str(benchmark.relative_to(root)),'benchmark_sha256':bsha,
                          'input':str(uv),'input_sha256':usha},
                'mnu_gev':encode(matrix),'mnu_ev':encode(matrix*1e9),
                'takagi_masses_ev_sorted':masses.tolist(),
                'source_vertex_pair_matrices_gev':{
                    'y1_HC x y1_HC':encode(blank),
                    'y1_HC x y2_HC':encode(matrix),
                    'y2_HC x y2_HC':encode(blank)},
                'source_charge_and_group':{
                    'native_group_audit_channels':ch,
                    'charge_partition_direct_exact':{
                        'charge_channel_matrices':{k:{'matrix_gev':encode(v)} for k,v in qmatrix.items()}}}}
        put_json(study/'models'/f'{model}.json',detail)
        weight_rows[model]={'channels':ch,'relative_physical_normalization_certified':False}
        for mode in MODES:
            for t in (HEAVY_VALUES if mode=='heavy_mass' else T_VALUES):
                scale=(1/t if mode=='heavy_mass' else abs(t))
                mm=masses*scale
                scans.append({'model':model,'mode':mode,'t':t,'status':'conditional_uv_scanned',
                              'reason':'','m1_ev':mm[0],'m2_ev':mm[1],'m3_ev':mm[2],
                              'mass_norm_gev':float(np.linalg.norm(matrix))*scale,
                              'cancellation_factor':1.5,
                              'matrix_scaling_residual':'' if mode=='heavy_mass' else 0.})
    report={'format':'t3-conditional-pattern-study-all16-v1',
            'physical_normalisation_certified':False,
            'counts':{'catalogue':16,'conditional_evaluated':16,'blocked':0,
                      'physical_source_certified':0,'scan_passes':288,'scan_blocked':0},
            'models':summary,
            'provenance':{'revision':{'local_head':'synthetic fixture'},
                          'source_audit':{'status':'source_forms_verified','source_files':{}},
                          'source_audits_verified_for_this_input':True,
                          'source_audit_sha256':{}}}
    put_json(study/'all16_summary.json',report)
    with (study/'all16_summary.csv').open('w',encoding='utf-8',newline='') as fp:
        wr=csv.DictWriter(fp,fieldnames=list(summary[0]));wr.writeheader();wr.writerows(summary)
    with (study/'scans.csv').open('w',encoding='utf-8',newline='') as fp:
        wr=csv.DictWriter(fp,fieldnames=list(scans[0]));wr.writeheader();wr.writerows(scans)
    put_json(study/'channel_weights.json',{'status':'native_group_CG_candidate_only_not_physical_normalization',
                                           'models':weight_rows})
    put_json(study/'common_kernel_group_only.json',
             {'physical_normalisation_certified':False,'rows':[{'model':x} for x in names]})
    return root,study


def test_synthetic_full_grid_analysis_preserves_source_results(tmp_path):
    root,source=catalogue_fixture(tmp_path)
    before={str(x.relative_to(source)):h(x) for x in source.rglob('*') if x.is_file()}
    result=analyse(root,figures=False)
    out=root/'output/general_observable/lo_pattern_study16/analysis_batch02_01'
    assert len(result['models'])==16
    assert result['counts']['source_scan_rows']==288
    assert result['counts']['physical_source_certified']==0
    assert all(r['yukawa_polynomial_max_mass_residual_relative_to_m3']['y1']<1e-12
               for r in result['models'])
    assert all(r['lambda_sign_flip_matrix_residual']==0 for r in result['models'])
    assert all(len(r['channel_overlaps'])==1 for r in result['models'] if r['family']=='E')
    assert (out/'batch02_model_summary.csv').is_file()
    assert (out/'BATCH02_REPORT.md').is_file()
    assert (out/'HANDOFF.json').is_file()
    after={str(x.relative_to(source)):h(x) for x in source.rglob('*') if x.is_file()}
    assert after==before


def test_blocks_overwrite_existing_report(tmp_path):
    root,_=catalogue_fixture(tmp_path)
    analyse(root,figures=False)
    with pytest.raises(FileExistsError):analyse(root,figures=False)


def test_blocks_changed_frozen_component_hash_before_writing(tmp_path):
    root,source=catalogue_fixture(tmp_path)
    p=root/'output/general_observable/lo_direct_one_loop/provisional_inputs/A_-4.json'
    p.write_text(p.read_text()+'\n')
    with pytest.raises(AnalysisError,match='Weyl input SHA256 drift'):
        analyse(root,figures=False)
    assert not (source.parent/'analysis_batch02_01').exists()


def test_blocks_false_physical_certification_promotion(tmp_path):
    root,source=catalogue_fixture(tmp_path)
    p=source/'models/B_-1.json'
    obj=json.loads(p.read_text())
    obj['physical_normalisation_certified']=True
    put_json(p,obj)
    with pytest.raises(AnalysisError,match='false certification'):
        analyse(root,figures=False)


def test_blocks_missing_scan_point(tmp_path):
    root,source=catalogue_fixture(tmp_path)
    p=source/'scans.csv'
    rows=list(csv.DictReader(p.open(newline='',encoding='utf-8')))
    with p.open('w',newline='',encoding='utf-8') as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows[1:])
    with pytest.raises(AnalysisError,match='Incomplete scan grid'):
        analyse(root,figures=False)


def test_blocks_yukawa_polynomial_violation(tmp_path):
    root,source=catalogue_fixture(tmp_path)
    p=source/'scans.csv'
    rows=list(csv.DictReader(p.open(newline='',encoding='utf-8')))
    for row in rows:
        if row['model']=='A_-4' and row['mode']=='y1' and float(row['t'])==2.:
            row['m3_ev']=str(float(row['m3_ev'])*1.1)
            break
    with p.open('w',newline='',encoding='utf-8') as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with pytest.raises(AnalysisError,match='quadratic identity failed'):
        analyse(root,figures=False)


def test_blocks_charge_partition_drift(tmp_path):
    root,source=catalogue_fixture(tmp_path)
    p=source/'models/E_-4.json'
    obj=json.loads(p.read_text())
    partition=obj['source_charge_and_group']['charge_partition_direct_exact']['charge_channel_matrices']
    first=next(iter(partition))
    partition[first]['matrix_gev']=encode(np.zeros((3,3)))
    put_json(p,obj)
    with pytest.raises(AnalysisError,match='charge partition cannot reconstruct'):
        analyse(root,figures=False)


def test_detects_executable_source_guard_drift(tmp_path):
    root,source=catalogue_fixture(tmp_path)
    guard=root/'Numerical/direct_one_loop/core.py'
    guard.parent.mkdir(parents=True,exist_ok=True)
    guard.write_text('trusted source snapshot',encoding='utf-8')
    p=source/'all16_summary.json'
    report=json.loads(p.read_text())
    report['provenance']['source_audit']['source_files']={
        'core':{'path':str(guard.relative_to(root)),'sha256':h(guard)}}
    put_json(p,report)
    guard.write_text('changed source snapshot',encoding='utf-8')
    with pytest.raises(AnalysisError,match='Executable source guard SHA256 drift'):
        analyse(root,figures=False)


def test_wrong_vertex_pair_relative_sign_fails_reconstruction(tmp_path):
    root,source=catalogue_fixture(tmp_path)
    p=source/'models/B_-3.json'
    obj=json.loads(p.read_text())
    z=obj['source_vertex_pair_matrices_gev']['y1_HC x y2_HC']
    z['real']=[[-v for v in row] for row in z['real']]
    z['imag']=[[-v for v in row] for row in z['imag']]
    put_json(p,obj)
    with pytest.raises(AnalysisError,match='vertex-pair sum cannot reconstruct'):
        analyse(root,figures=False)
